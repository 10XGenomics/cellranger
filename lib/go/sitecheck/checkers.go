// Copyright (c) 2024 10x Genomics, Inc. All rights reserved.

package sitecheck

import (
	"bufio"
	"bytes"
	"context"
	"errors"
	"fmt"
	"os"
	"os/exec"
	"path"
	"path/filepath"
	"regexp"
	"runtime"
	"sort"
	"strconv"
	"strings"
	"syscall"
	"time"

	"golang.org/x/sys/unix"
)

// #include <gnu/libc-version.h>
import "C"

func getLibcVersion() string {
	return C.GoString(C.gnu_get_libc_version())
}

// Convert a null-terminated byte array to a string.
func cstring(b []byte) string {
	if before, _, found := bytes.Cut(b, nullByte); !found {
		return ""
	} else {
		return string(before)
	}
}

// Collect the uname data.
func checkUname(_ context.Context, index int, results chan<- siteCheckSection) {
	var buf unix.Utsname
	if err := unix.Uname(&buf); err != nil {
		fmt.Fprintf(os.Stderr, "Could not get OS info: %v\n",
			err)
		return
	}
	results <- siteCheckSection{
		Section: "System Info",
		Cmd:     "uname -a",
		Output: fmt.Sprintf("%s %s %s %s %s",
			cstring(buf.Sysname[:]),
			cstring(buf.Nodename[:]),
			cstring(buf.Release[:]),
			cstring(buf.Version[:]),
			cstring(buf.Machine[:]),
		),
		Index: index,
	}
}

// Collect information about the OS distribution.
func checkDistro(ctx context.Context, index int, results chan<- siteCheckSection) {
	matches, err := filepath.Glob("/etc/*-release")
	if err != nil {
		return
	}
	lines := make(map[string]struct{}, 16)
	for _, f := range matches {
		b, err := os.ReadFile(f)
		if err != nil || ctx.Err() != nil {
			continue
		}
		for _, line := range strings.Split(string(bytes.TrimSpace(b)), "\n") {
			lines[line] = struct{}{}
		}
	}
	result := make([]string, 0, len(lines))
	for line := range lines {
		result = append(result, line)
	}
	sort.Strings(result)

	results <- siteCheckSection{
		Section: "Linux Distro",
		Cmd:     "cat /etc/*-release | sort -u",
		Output:  strings.Join(result, "\n"),
		Index:   index,
	}
}

func catFiles(section string, index int, results chan<- siteCheckSection, names ...string) error {
	if len(names) == 1 {
		b, err := os.ReadFile(names[0])
		if err != nil {
			return err
		}
		results <- siteCheckSection{
			Section: section,
			Cmd:     "cat " + names[0],
			Output:  string(bytes.TrimSpace(b)),
			Index:   index,
		}
		return nil
	}
	var result strings.Builder
	var lastErr error
	for _, name := range names {
		b, err := os.ReadFile(name)
		if err != nil {
			lastErr = err
			continue
		}
		b = bytes.TrimSpace(b)
		if result.Len() > 0 && len(b) > 0 {
			result.WriteByte('\n')
		}
		result.Write(b)
	}
	if result.Len() > 0 {
		results <- siteCheckSection{
			Section: section,
			Cmd:     "cat " + strings.Join(names, " "),
			Output:  result.String(),
			Index:   index,
		}
	}
	return lastErr
}

// Collect the kernel version information.
func checkKernel(_ context.Context, index int, results chan<- siteCheckSection) {
	if err := catFiles("Kernel Build", index, results, "/proc/version"); err != nil {
		fmt.Fprintf(os.Stderr, "Could not get kernel info: %v\n",
			err)
	}
}

// Collect glibc version.
func checkGlibc(_ context.Context, index int, results chan<- siteCheckSection) {
	results <- siteCheckSection{
		Section: "glibc version",
		Cmd:     "ldd --version | head -n 1",
		Output:  "ldd (GNU libc) " + getLibcVersion(),
		Index:   index,
	}
}

// Reusable byte arrays
var (
	colon    = []byte{':'}
	slash    = []byte{'/'}
	space    = []byte{'/'}
	newline  = []byte{'\n'}
	nullByte = []byte{0}
)

func lineValue(line, search []byte) ([]byte, bool) {
	after, found := bytes.CutPrefix(line, search)
	if !found {
		return nil, false
	}
	_, after, found = bytes.Cut(after, colon)
	if !found {
		return nil, false
	}
	return bytes.TrimSpace(after), true
}

// Collect information about the CPU model, features, and count.
func checkCpu(_ context.Context, index int, results chan<- siteCheckSection) {
	f, err := os.Open("/proc/cpuinfo")
	if err != nil {
		fmt.Fprintf(os.Stderr, "Could not get cpu info: %v\n",
			err)
		return
	}
	defer f.Close()
	scanner := bufio.NewScanner(f)
	var procs int
	sockets := make(map[string]struct{}, runtime.NumCPU())
	var foundModel, foundFlags bool
	modelSearch := []byte("model name")
	flagSearch := []byte("flags")
	physIdSearch := []byte("physical id")
	procSearch := []byte("processor")
	for scanner.Scan() {
		if !foundModel {
			if v, found := lineValue(scanner.Bytes(), modelSearch); found {
				foundModel = true
				results <- siteCheckSection{
					Section: "CPU Model",
					Cmd: "grep -m 1 'model name' /proc/cpuinfo | " +
						"cut -d ':' -f 2 | sed 's/^[ \t]*//'",
					Output: string(v),
				}
				continue
			}
		}
		if !foundFlags {
			if v, found := lineValue(scanner.Bytes(), flagSearch); found {
				foundFlags = true
				results <- siteCheckSection{
					Section: "CPU Support",
					Cmd: "grep -m 1 'flags' /proc/cpuinfo | " +
						"cut -d ':' -f 2 | sed 's/^[ \t]*//'",
					Output: string(v),
					Index:  index,
				}
				continue
			}
		}
		if v, found := lineValue(scanner.Bytes(), physIdSearch); found {
			sockets[string(v)] = struct{}{}
		} else if _, found = lineValue(scanner.Bytes(), procSearch); found {
			procs++
		}
	}
	results <- siteCheckSection{
		Section: "CPU Sockets",
		Cmd:     "grep 'physical id' /proc/cpuinfo | sort -u | wc -l",
		Output:  strconv.Itoa(len(sockets)),
		Index:   index,
	}
	results <- siteCheckSection{
		Section: "CPU Cores",
		Cmd:     "grep -c processor /proc/cpuinfo",
		Output:  strconv.Itoa(procs),
		Index:   index,
	}
}

// Get total system memory amount
func checkMemory(_ context.Context, index int, results chan<- siteCheckSection) {
	f, err := os.Open("/proc/meminfo")
	if err != nil {
		fmt.Fprintf(os.Stderr, "Could not get mem info: %v\n",
			err)
		return
	}
	defer f.Close()
	scanner := bufio.NewScanner(f)
	totalSearch := []byte("MemTotal")
	for scanner.Scan() {
		if v, found := lineValue(scanner.Bytes(), totalSearch); found {
			results <- siteCheckSection{
				Section: "Memory Total",
				Cmd:     "grep MemTotal /proc/meminfo | cut -d ':' -f 2 | sed 's/^[ \t]*//'",
				Output:  string(v),
				Index:   index,
			}
			return
		}
	}
	fmt.Fprintln(os.Stderr, "Could not find MemTotal in /proc/meminfo")
}

// Check disk free space and mount options.
func checkDisk(_ context.Context, index int, results chan<- siteCheckSection) {
	f, err := os.Open("/proc/self/mountinfo")
	if err != nil {
		fmt.Fprintf(os.Stderr, "Could not get mount info: %v\n",
			err)
		return
	}
	defer f.Close()
	scanner := bufio.NewScanner(f)
	var spaceBuf, optsBuf strings.Builder
	spaceBuf.WriteString("Size Used Avail")
	spaceStartLen := spaceBuf.Len()
	autofs := []byte("autofs")
	for scanner.Scan() {
		// From `man 5 proc`:
		// /proc/[pid]/mountinfo (since Linux 2.6.26)
		// This file contains information about mount points.  It contains lines of the form:
		//
		// 36 35 98:0 /mnt1 /mnt2 rw,noatime master:1 - ext3 /dev/root rw,errors=continue
		// (1)(2)(3)   (4)   (5)      (6)      (7)   (8) (9)   (10)         (11)
		//
		// The numbers in parentheses are labels for the descriptions below:
		// (1)  mount ID: unique identifier of the mount (may be reused after umount(2)).
		// (2)  parent ID: ID of parent mount (or of self for the top of the mount tree).
		// (3)  major:minor: value of st_dev for files on file system (see stat(2)).
		// (4)  root: root of the mount within the file system.
		// (5)  mount point: mount point relative to the process's root.
		// (6)  mount options: per-mount options.
		// (7)  optional fields: zero or more fields of the form "tag[:value]".
		// (8)  separator: marks the end of the optional fields.
		// (9)  file system type: name of file system in the form "type[.subtype]".
		// (10) mount source: file system-specific information or "none".
		// (11) super options: per-super block options.
		fields := bytes.Fields(scanner.Bytes())
		if len(fields) < 9 {
			continue
		}
		sep := 7
		for i, f := range fields[5:] {
			if len(f) == 1 && f[0] == '-' {
				sep = 5 + i
				break
			}
		}
		kind := fields[sep+1]
		optsBuf.Write(kind)
		optsBuf.WriteString(" (")
		optsBuf.Write(fields[5])
		optsBuf.WriteString(")\n")
		if !bytes.Equal(kind, autofs) {
			getSpace(&spaceBuf, string(fields[4]))
		}
	}
	if space := spaceBuf.String(); len(space) > spaceStartLen {
		results <- siteCheckSection{
			Section: "Disk Space",
			Cmd:     "df -Ph | awk '{print $2, $3, $4}'",
			Output:  space,
			Index:   index,
		}
	}
	if opts := strings.TrimSpace(optsBuf.String()); len(opts) > 0 {
		results <- siteCheckSection{
			Section: "Filesystem Options",
			Cmd:     "mount | cut -d ' ' -f 5,6",
			Output:  opts,
			Index:   index,
		}
	}
}

// Convert block size and number of blocks into human-readable value.
func humanSpace(bsize int64, amount uint64) (string, rune) {
	total := uint64(bsize) * amount
	if total < 1024 {
		return strconv.FormatUint(total, 10), 0
	}
	for i, suffix := range sizeSuffixes[:len(sizeSuffixes)-1] {
		unit := ((total >> (9 + 10*i)) + 1) / 2
		if unit < 1024 {
			return strconv.FormatUint(unit, 10), suffix
		}
	}
	unit := ((total >> (9 + 40)) + 1) / 2
	return strconv.FormatUint(unit, 10), 'P'
}

// Write the total/used/avail space for a mount to the given string builder.
func getSpace(builder *strings.Builder, mount string) {
	if strings.HasPrefix(mount, "/sys/kernel/debug") {
		return
	}
	var buf unix.Statfs_t
	err := unix.Statfs(mount, &buf)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Could not get disk stats for %s: %v\n",
			mount, err)
		return
	}
	if buf.Blocks == 0 {
		return
	}
	builder.WriteByte('\n')
	a, s := humanSpace(buf.Bsize, buf.Blocks)
	builder.WriteString(a)
	if s != 0 {
		builder.WriteRune(s)
	}
	builder.WriteByte(' ')
	a, s = humanSpace(buf.Bsize, buf.Blocks-buf.Bfree)
	builder.WriteString(a)
	if s != 0 {
		builder.WriteRune(s)
	}
	builder.WriteByte(' ')
	a, s = humanSpace(buf.Bsize, buf.Bavail)
	builder.WriteString(a)
	if s != 0 {
		builder.WriteRune(s)
	}
}

func writeRlim(v, scaleN, scaleD uint64, buf *strings.Builder) {
	if v == unix.RLIM_INFINITY {
		buf.WriteString("unlimited")
	} else {
		buf.WriteString(strconv.FormatUint((v*scaleN)/scaleD, 10))
	}
	buf.WriteByte('\n')
}

func writeRlims(resource int, head string, soft, hard *strings.Builder,
	scaleN, scaleD uint64) {
	var rlim unix.Rlimit
	if unix.Getrlimit(resource, &rlim) != nil {
		return
	}
	soft.WriteString(head)
	writeRlim(rlim.Cur, scaleN, scaleD, soft)
	hard.WriteString(head)
	writeRlim(rlim.Max, scaleN, scaleD, hard)
}

// Report rlimit values.
func checkUlimits(_ context.Context, index int, results chan<- siteCheckSection) {
	var soft, hard strings.Builder
	writeRlims(unix.RLIMIT_CORE,
		"core file size          (blocks, -c) ",
		&soft, &hard, 1, 1)
	writeRlims(unix.RLIMIT_DATA,
		"data seg size           (kbytes, -d) ",
		&soft, &hard, 1, 1024)
	writeRlims(unix.RLIMIT_NICE,
		"scheduling priority             (-e) ",
		&soft, &hard, 1, 1)
	writeRlims(unix.RLIMIT_FSIZE,
		"file size               (blocks, -f) ",
		&soft, &hard, 1, 1)
	writeRlims(unix.RLIMIT_SIGPENDING,
		"pending signals                 (-i) ",
		&soft, &hard, 1, 1)
	writeRlims(unix.RLIMIT_MEMLOCK,
		"max locked memory       (kbytes, -l) ",
		&soft, &hard, 1, 1024)
	writeRlims(unix.RLIMIT_RSS,
		"max memory size         (kbytes, -m) ",
		&soft, &hard, uint64(unix.Getpagesize()), 1024)
	writeRlims(unix.RLIMIT_NOFILE,
		"open files                      (-n) ",
		&soft, &hard, 1, 1)
	writeRlims(unix.RLIMIT_MSGQUEUE,
		"POSIX message queues     (bytes, -q) ",
		&soft, &hard, 1, 1)
	writeRlims(unix.RLIMIT_RTPRIO,
		"real-time priority              (-r) ",
		&soft, &hard, 1, 1)
	writeRlims(unix.RLIMIT_STACK,
		"stack size              (kbytes, -s) ",
		&soft, &hard, 1, 1024)
	writeRlims(unix.RLIMIT_CPU,
		"cpu time               (seconds, -t) ",
		&soft, &hard, 1, 1)
	writeRlims(unix.RLIMIT_NPROC,
		"max user processes              (-u) ",
		&soft, &hard, 1, 1)
	writeRlims(unix.RLIMIT_AS,
		"virtual memory          (kbytes, -v) ",
		&soft, &hard, 1, 1024)
	writeRlims(unix.RLIMIT_LOCKS,
		"file locks                      (-x) ",
		&soft, &hard, 1, 1)
	results <- siteCheckSection{
		Section: "User Limits",
		Cmd:     "bash -c 'ulimit -a'",
		Output:  strings.TrimSpace(soft.String()),
		Index:   index,
	}
	results <- siteCheckSection{
		Section: "User Limits (hard)",
		Cmd:     "bash -c 'ulimit -aH'",
		Output:  strings.TrimSpace(hard.String()),
		Index:   index,
	}
}

func checkFileLimit(_ context.Context, index int, results chan<- siteCheckSection) {
	if err := catFiles("Global File Limit", index, results,
		"/proc/sys/fs/file-max", "/proc/sys/fs/file-nr"); err != nil {
		fmt.Fprintf(os.Stderr, "Could not get file max: %v\n",
			err)
	}
}

// Report kernel memory settings.
func checkVmConfig(_ context.Context, index int, results chan<- siteCheckSection) {
	if entries, err := os.ReadDir("/proc/sys/vm"); err != nil {
		fmt.Fprintf(os.Stderr, "Could not list vm proc directory: %v\n",
			err)
	} else {
		var buf strings.Builder
		for _, ent := range entries {
			b, err := os.ReadFile(path.Join("/proc/sys/vm", ent.Name()))
			if err != nil {
				continue
			}
			buf.WriteString("vm.")
			buf.WriteString(ent.Name())
			buf.WriteString(" = ")
			buf.Write(bytes.TrimSpace(b))
			buf.WriteByte('\n')
		}
		results <- siteCheckSection{
			Section: "Memory config",
			Cmd:     "sysctl vm",
			Output:  strings.TrimSpace(buf.String()),
			Index:   index,
		}
	}
	if thp, err := filepath.Glob("/sys/kernel/mm/*transparent_hugepage/enabled"); err != nil {
		fmt.Fprintf(os.Stderr, "Could not list hugepage config files: %v\n",
			err)
	} else if len(thp) > 0 {
		if err := catFiles("THP memory config", index, results, thp...); err != nil {
			fmt.Fprintf(os.Stderr, "Could not read hugepage config file: %v\n",
				err)
		}
	}
}

// Report information about cgroups.
func checkCgroups(_ context.Context, index int, results chan<- siteCheckSection) {
	b, err := os.ReadFile("/proc/self/cgroup")
	if err != nil {
		// no cgroup
		return
	}
	content := string(bytes.TrimSpace(b))
	results <- siteCheckSection{
		Section: "cgroups",
		Cmd:     "cat /proc/self/cgroup",
		Output:  content,
		Index:   index,
	}
	for _, line := range strings.Split(content, "\n") {
		if _, after, found := strings.Cut(line, ":memory:"); found {
			checkMemoryCgroup(after, false, index, results)
			return
		} else if _, after, found := strings.Cut(line, "::"); found {
			checkMemoryCgroup(after, true, index, results)
			return
		}
	}
}

func findMemCgroup() []byte {
	f, err := os.Open("/proc/self/mountinfo")
	if err != nil {
		return nil
	}
	defer f.Close()
	scanner := bufio.NewScanner(f)
	search := regexp.MustCompile(`\bmemory\b`)
	for scanner.Scan() {
		fields := bytes.Fields(scanner.Bytes())
		if len(fields) < 11 {
			continue
		}
		if search.Match(fields[10]) {
			return fields[4]
		}
	}
	return nil
}

func findCgroupV2() []byte {
	f, err := os.Open("/proc/self/mountinfo")
	if err != nil {
		return nil
	}
	defer f.Close()
	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		fields := bytes.Fields(scanner.Bytes())
		if len(fields) < 11 {
			continue
		}
		if string(fields[8]) == "cgroup2" {
			return fields[4]
		}
	}
	return nil
}

// See https://github.com/golang/go/issues/67758 for why this needs to exist.
func makeString(b []byte) string {
	return string(b)
}

func checkMemoryCgroup(cgroup string, v2 bool, index int, results chan<- siteCheckSection) {
	var cgroupRoot []byte
	var (
		softLimit     string
		softSwapLimt  string
		hardLimit     string
		hardSwapLimit string
	)
	if v2 {
		cgroupRoot = findCgroupV2()
		softLimit = "high"
		softSwapLimt = "swap.high"
		hardLimit = "max"
		hardSwapLimit = "swap.max"
	} else {
		cgroupRoot = findMemCgroup()
		softLimit = "soft_limit_in_bytes"
		softSwapLimt = "memsw.soft_limit_in_bytes"
		hardLimit = "limit_in_bytes"
		hardSwapLimit = "memsw.limit_in_bytes"
	}
	if len(cgroupRoot) == 0 {
		return
	}
	cgroupRoot = bytes.TrimSuffix(cgroupRoot, slash)
	cgroupBuf := make([]byte, len(cgroupRoot),
		len(cgroupRoot)+len(cgroup)+1+len("memory.memsw.soft_limit_in_bytes"))
	copy(cgroupBuf, cgroupRoot)
	cgroupBuf = append(cgroupBuf, strings.TrimSuffix(cgroup, "/")...)
	cgroupBuf = append(cgroupBuf, "/memory."...)
	if err := catFiles("cgroup mem stats", index, results,
		string(append(cgroupBuf, "stat"...))); err != nil {
		fmt.Fprintf(os.Stderr, "Could not read cgroup stats: %v\n",
			err)
	}
	if err := catFiles("memory soft limit", index, results,
		makeString(append(cgroupBuf, softLimit...)),
		makeString(append(cgroupBuf, softSwapLimt...))); err != nil &&
		!errors.Is(err, os.ErrNotExist) {
		fmt.Fprintf(os.Stderr, "Could not read cgroup soft limit: %v\n",
			err)
	}
	if err := catFiles("memory hard limit", index, results,
		string(append(cgroupBuf, hardLimit...))); err != nil &&
		!errors.Is(err, os.ErrNotExist) {
		fmt.Fprintf(os.Stderr, "Could not read cgroup hard limit: %v\n",
			err)
	}
	if err := catFiles("memory swap limit", index, results,
		string(append(cgroupBuf, hardSwapLimit...))); err != nil &&
		!errors.Is(err, os.ErrNotExist) {
		fmt.Fprintf(os.Stderr, "Could not read cgroup swap limit: %v\n",
			err)
	}
}

func fExists(p string) bool {
	_, err := os.Stat(p)
	return !errors.Is(err, os.ErrNotExist)
}

func isContainer() bool {
	if fExists("/.dockerenv") || fExists("/.dockerinit") ||
		os.Getenv("container") != "" {
		return true
	}
	b, err := os.ReadFile("/proc/1/cgroup")
	if err != nil {
		return false
	}
	if bytes.Contains(b, []byte("docker")) ||
		bytes.Contains(b, []byte("lxc")) {
		return true
	}
	return false
}

// Report our best guess as to whether we're running in a container.
func checkContainer(_ context.Context, index int, results chan<- siteCheckSection) {
	var output string
	if isContainer() {
		output = "Detected"
	}
	results <- siteCheckSection{
		Section: "Container",
		Cmd: `[ -e /.dockerenv ] || [ -e /.dockerinit ] \
		|| [ ! -z "$container" ] || grep -m 1 -E 'docker|lxc' /proc/1/cgroup \
		> /dev/null && echo 'Detected'`,
		Output: output,
		Index:  index,
	}
}

// Report the current init method (e.g. systemd).
func checkInit(_ context.Context, index int, results chan<- siteCheckSection) {
	f, err := os.Open("/proc/1/sched")
	if err != nil {
		fmt.Fprintf(os.Stderr, "Could not read init process info: %v\n",
			err)
	}
	defer f.Close()
	scanner := bufio.NewScanner(f)
	if scanner.Scan() {
		prefix, _, _ := bytes.Cut(scanner.Bytes(), space)
		results <- siteCheckSection{
			Section: "init process",
			Cmd:     "head -n 1 /proc/1/sched | cut -d ' ' -f 1",
			Output:  string(prefix),
			Index:   index,
		}
	}
}

// Run an executable and return its output.
func execOut(ctx context.Context, stderr bool, e string, args ...string) []byte {
	ctx, cancel := context.WithTimeout(ctx, time.Second*20)
	defer cancel()
	cmd := exec.CommandContext(ctx, e, args...)
	cmd.SysProcAttr = &syscall.SysProcAttr{
		Pdeathsig: syscall.SIGKILL,
	}
	var buf bytes.Buffer
	cmd.Stdout = &buf
	if stderr {
		cmd.Stderr = &buf
	}
	if err := cmd.Run(); err != nil && errors.Is(err, exec.ErrNotFound) {
		return nil
	}
	if ctx.Err() != nil {
		fmt.Fprintln(os.Stderr, e, "timed out")
	}
	return bytes.TrimSpace(buf.Bytes())
}

type hostMemData struct {
	name []byte
	mem  float64
}

const sizeSuffixes = "kMGTP"

var sizeSuffixesBytes = []byte(sizeSuffixes)

func humanParse(b []byte) float64 {
	var v float64
	for i := range sizeSuffixesBytes {
		if n, found := bytes.CutSuffix(b, sizeSuffixesBytes[i:i+1]); found {
			v, _ = strconv.ParseFloat(string(n), 64)
			return v * float64(int64(1)<<((i+1)*10))
		}
	}
	v, _ = strconv.ParseFloat(string(b), 64)
	return v
}

func parseQhost(line []byte) hostMemData {
	f := bytes.Fields(line)
	if len(f) == 0 {
		return hostMemData{}
	}
	if len(f) < 5 {
		return hostMemData{
			name: f[0],
		}
	}
	return hostMemData{
		name: f[0],
		mem:  humanParse(f[4]),
	}
}

// Report SGE parameters.
func checkSGE(ctx context.Context, index int, results chan<- siteCheckSection) {
	p, _ := exec.LookPath("qsub")
	results <- siteCheckSection{
		Section: "SGE Submit",
		Cmd:     "which qsub",
		Output:  p,
		Index:   index,
	}
	if p == "" {
		return
	}
	results <- siteCheckSection{
		Section: "SGE CLUSTER_NAME",
		Cmd:     "echo $SGE_CLUSTER_NAME",
		Output:  os.Getenv("SGE_CLUSTER_NAME"),
		Index:   index,
	}
	results <- siteCheckSection{
		Section: "SGE JOB_NAME",
		Cmd:     "echo $JOB_NAME",
		Output:  os.Getenv("JOB_NAME"),
		Index:   index,
	}

	p, _ = exec.LookPath("qconf")
	results <- siteCheckSection{
		Section: "qconf",
		Cmd:     "which qconf",
		Output:  p,
		Index:   index,
	}
	if p == "" {
		return
	}
	qsc := execOut(ctx, false, p, "-sc")
	sconf := bytes.Split(execOut(ctx, false, p, "-sconf"), newline)
	var qconfBuf strings.Builder
	qconfBuf.Write(qsc)
	for _, line := range sconf {
		if bytes.Contains(line, []byte("shell_start_mode")) ||
			bytes.Contains(line, []byte("login_shells")) ||
			bytes.Contains(line, []byte("max_jobs")) {
			qconfBuf.WriteByte('\n')
			qconfBuf.Write(line)
		}
	}
	results <- siteCheckSection{
		Section: "qconf -sc",
		Cmd:     "qconf -sc && qconf -sconf | grep -E '(shell_start_mode|login_shells|max_jobs)'",
		Output:  qconfBuf.String(),
		Index:   index,
	}

	p, _ = exec.LookPath("qhost")
	results <- siteCheckSection{
		Section: "qhost",
		Cmd:     "which qhost",
		Output:  p,
		Index:   index,
	}
	if p == "" {
		return
	}
	qhostOut := execOut(ctx, false, p, "-l", "mem_total=20G")
	_, qhostOut, _ = bytes.Cut(qhostOut, newline)
	_, qhostOut, _ = bytes.Cut(qhostOut, newline)
	if len(qhostOut) == 0 {
		results <- siteCheckSection{
			Section: "qhost count",
			Cmd:     "qhost -l \"mem_total=20G\" | tail -n +3 | wc -l",
			Output:  "0",
			Index:   index,
		}
		return
	}
	hosts := bytes.Split(qhostOut, newline)
	results <- siteCheckSection{
		Section: "qhost count",
		Cmd:     "qhost -l \"mem_total=20G\" | tail -n +3 | wc -l",
		Output:  strconv.Itoa(len(hosts)),
		Index:   index,
	}
	var biggest hostMemData
	for _, h := range hosts {
		data := parseQhost(h)
		if data.mem > biggest.mem {
			biggest = data
		}
	}
	results <- siteCheckSection{
		Section: "qhost -F",
		Cmd:     "qhost -F -q -h $(qhost | sort -h -k 5 -r | head -n 1 | cut -d \" \" -f 1)",
		Output: string(execOut(ctx, false, p,
			"-F", "-q", "-h", string(biggest.name))),
		Index: index,
	}
}

// Report LSF information.
func checkLSF(_ context.Context, index int, results chan<- siteCheckSection) {
	p, _ := exec.LookPath("bsub")
	results <- siteCheckSection{
		Section: "LSF Submit",
		Cmd:     "which bsub",
		Output:  p,
		Index:   index,
	}
	if p == "" {
		return
	}
	results <- siteCheckSection{
		Section: "LSF LSB_JOBNAME",
		Cmd:     "echo $LSB_JOBNAME",
		Output:  os.Getenv("LSB_JOBNAME"),
		Index:   index,
	}
}

// Report HTCondor information.
func checkCondor(_ context.Context, index int, results chan<- siteCheckSection) {
	p, _ := exec.LookPath("condor_submit")
	results <- siteCheckSection{
		Section: "HTCondor Submit",
		Cmd:     "which condor_submit",
		Output:  p,
		Index:   index,
	}
	results <- siteCheckSection{
		Section: "Batch system",
		Cmd:     "echo $BATCH_SYSTEM",
		Output:  os.Getenv("BATCH_SYSTEM"),
		Index:   index,
	}
}

func listDir(dir string, n int) string {
	f, err := os.Open(dir)
	if err != nil {
		return ""
	}
	defer f.Close()
	names, _ := f.Readdirnames(n)
	return strings.Join(names, "\n")
}

// Report information relevant to BCL2FASTQ.
func checkBCL2FASTQ(ctx context.Context, index int, results chan<- siteCheckSection) {
	p, _ := exec.LookPath("configureBclToFastq.pl")
	results <- siteCheckSection{
		Section: "BCL2FASTQ 1",
		Cmd:     "which configureBclToFastq.pl",
		Output:  p,
		Index:   index,
	}
	if p != "" {
		results <- siteCheckSection{
			Section: "BCL2FASTQ 1 Version",
			Cmd:     "ls $(dirname $(which configureBclToFastq.pl))/../etc",
			Output:  listDir(path.Join(path.Dir(path.Dir(p)), "etc"), 10),
			Index:   index,
		}
		p, _ := exec.LookPath("perl")
		results <- siteCheckSection{
			Section: "Perl",
			Cmd:     "which perl",
			Output:  p,
			Index:   index,
		}
		if p != "" {
			if o := execOut(ctx, false, p, "-v"); len(o) > 0 {
				results <- siteCheckSection{
					Section: "Perl Version",
					Cmd:     "perl -v",
					Output:  string(o),
					Index:   index,
				}
			}
		}
	}
	p, _ = exec.LookPath("bcl2fastq")
	results <- siteCheckSection{
		Section: "BCL2FASTQ 1",
		Cmd:     "which bcl2fastq",
		Output:  p,
		Index:   index,
	}
	if p != "" {
		if o := execOut(ctx, false, p, "--version"); len(o) > 0 {
			results <- siteCheckSection{
				Section: "BCL2FASTQ 2 Version",
				Cmd:     "bcl2fastq --version",
				Output:  string(o),
				Index:   index,
			}
		}
	}
}

// Report java information.
func checkJava(ctx context.Context, index int, results chan<- siteCheckSection) {
	p, _ := exec.LookPath("java")
	results <- siteCheckSection{
		Section: "Java",
		Cmd:     "which java",
		Output:  p,
		Index:   index,
	}
	if p != "" {
		if o := execOut(ctx, true, p, "-version"); len(o) > 0 {
			results <- siteCheckSection{
				Section: "Java Version",
				Cmd:     "java -version 2>&1 | cat",
				Output:  string(o),
				Index:   index,
			}
		}
	}
}

func checkRefdata(_ context.Context, index int, results chan<- siteCheckSection) {
	refdata := os.Getenv("TENX_REFDATA")
	results <- siteCheckSection{
		Section: "10X Refdata",
		Cmd:     "echo $TENX_REFDATA",
		Output:  refdata,
		Index:   index,
	}
	if refdata != "" {
		if err := catFiles("10X Refdata Version", index, results,
			path.Join(refdata, "version")); err != nil {
			fmt.Fprintf(os.Stderr, "Could not read refdata version: %v\n",
				err)
		}
	}
}

func checkSlurm(ctx context.Context, index int, results chan<- siteCheckSection) {
	o := execOut(ctx, false, "sinfo", "-O",
		"nodes,maxcpuspernode,memory,time")
	results <- siteCheckSection{
		Section: "slurm info",
		Cmd:     "sinfo -O nodes,maxcpuspernode,memory,time",
		Output:  string(o),
		Index:   index,
	}
}

// Report mrp path and version.
func checkMrp(ctx context.Context, index int, results chan<- siteCheckSection) {
	p, _ := exec.LookPath("mrp")
	if p == "" {
		results <- siteCheckSection{
			Section: "MRP",
			Cmd:     "mrp --version",
			Output:  "",
			Index:   index,
		}
		return
	}
	o := execOut(ctx, false, "mrp", "--version")
	results <- siteCheckSection{
		Section: "MRP",
		Cmd:     "mrp --version",
		Output:  string(o),
		Index:   index,
	}
	if len(o) == 0 {
		return
	}
	if ctx.Err() != nil {
		return
	}
	dir := path.Join(path.Dir(path.Dir(p)), "jobmanagers")
	entries, _ := os.ReadDir(dir)
	var listBuf strings.Builder
	for _, ent := range entries {
		if !strings.HasSuffix(ent.Name(), ".template") {
			continue
		}
		if listBuf.Len() > 0 {
			listBuf.WriteByte('\n')
		}
		f := path.Join(dir, ent.Name())
		listBuf.WriteString(f)
		_ = catFiles(f, index, results, f)
	}
	results <- siteCheckSection{
		Section: "mrp templates",
		Cmd:     "ls $(dirname $(dirname $(which mrp)))/jobmanagers/*.template",
		Output:  listBuf.String(),
		Index:   index,
	}
}
