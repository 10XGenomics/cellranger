package collector

import (
	"bufio"
	"bytes"
	"context"
	"encoding/json"
	"os"
	"path"
	"runtime/trace"
	"strconv"
	"strings"
	"sync"
	"syscall"

	"golang.org/x/sys/unix"
)

func doStatfs(ctx context.Context, p string, cs []chan<- unix.Statfs_t) {
	defer trace.StartRegion(ctx, "doStatfs").End()
	for _, c := range cs {
		cc := c // don't capture the loop variable
		defer close(cc)
	}
	if ctx.Err() != nil {
		return
	}
	var buf unix.Statfs_t
	if err := unix.Statfs(p, &buf); err == nil {
		for _, c := range cs {
			c <- buf
		}
	}
}

type mountInfo struct {
	mType string
	opts  string
}

func parseMountLine(line []byte, want map[string][]chan<- mountInfo, cg map[string]string) {
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

	fields := bytes.Fields(line)
	if len(fields) < 9 {
		return
	}
	results := want[string(fields[2])]
	if len(results) == 0 && cg == nil {
		return
	}
	var fsType string
	opts := fields[5][0:len(fields[5]):len(fields[5])]
	for i, f := range fields[6:] {
		if len(f) == 1 && f[0] == '-' {
			if len(fields) >= i+7 {
				fsType = string(fields[i+7])
				if len(results) > 0 && len(fields) >= i+9 && len(fields[i+9]) > 0 {
					opts = append(opts, ',')
					opts = append(opts, fields[i+9]...)
				} else if fsType == "cgroup" || fsType == "cgroup2" {
					opts = fields[i+9]
				}
			}
			break
		}
	}
	if len(results) > 0 {
		info := mountInfo{
			mType: fsType,
			opts:  string(opts),
		}
		for _, r := range results {
			r <- info
		}
	}
	if cg != nil {
		if fsType == "cgroup2" {
			cg[""] = string(fields[4])
		} else if fsType == "cgroup" {
			mount := string(fields[4])
			for _, ctrl := range strings.Split(string(bytes.TrimPrefix(opts, []byte("rw,"))), ",") {
				cg[ctrl] = mount
			}
		}
	}
}

func makeMountId(p string) []byte {
	mountId := make([]byte, 0, 21)
	if info, err := os.Stat(p); err != nil || info == nil {
		return nil
	} else if sysInfo, ok := info.Sys().(*syscall.Stat_t); !ok {
		return nil
	} else {
		itoa := func(i uint32) {
			if i == 0 {
				mountId = append(mountId, '0')
			} else {
				mountId = append(mountId, strconv.Itoa(int(i))...)
			}
		}
		itoa(unix.Major(sysInfo.Dev))
		mountId = append(mountId, ':')
		itoa(unix.Minor(sysInfo.Dev))
		return mountId
	}
}

func doMountInfos(ctx context.Context, request map[string][]chan<- mountInfo,
	cgroupMounts chan<- map[string]string) {
	defer trace.StartRegion(ctx, "doMountInfos").End()
	if cgroupMounts != nil {
		defer close(cgroupMounts)
	}
	for _, cs := range request {
		for _, c := range cs {
			cc := c // don't capture the loop variable
			defer close(cc)
		}
	}
	if ctx.Err() != nil {
		return
	}
	devs := make(map[string][]chan<- mountInfo, len(request))
	for p, cs := range request {
		id := makeMountId(p)
		if len(id) == 0 {
			continue
		}
		sid := string(id)
		if ecs := devs[sid]; len(ecs) == 0 {
			devs[sid] = cs
		} else {
			devs[sid] = append(ecs, cs...)
		}
	}
	if ctx.Err() != nil {
		return
	}

	getMountInfos(devs, cgroupMounts)
}

func getMountInfos(devs map[string][]chan<- mountInfo,
	cgroupMounts chan<- map[string]string) {
	f, err := os.Open("/proc/self/mountinfo")
	if err != nil {
		return
	}
	defer f.Close()
	scan := bufio.NewScanner(f)
	var cg map[string]string
	if cgroupMounts != nil {
		cg = make(map[string]string)
	}
	for scan.Scan() {
		parseMountLine(scan.Bytes(), devs, cg)
	}
	if cgroupMounts != nil {
		cgroupMounts <- cg
	}
}

type FilesystemProperty int

const (
	FsKind = FilesystemProperty(iota)
	FsOpts
	FsCapacity
	FsInodeCap
	FsAvail
	FsInodeAvail
)

func (p FilesystemProperty) String() string {
	switch p {
	case FsKind:
		return "kind"
	case FsOpts:
		return "mnt_opts"
	case FsCapacity:
		return "bytes_cap"
	case FsInodeCap:
		return "inode_cap"
	case FsAvail:
		return "bytes_avail"
	case FsInodeAvail:
		return "inode_avail"
	}
	return "unknown"
}

var exeDir = sync.OnceValue[string](func() string {
	exe, err := os.Executable()
	if err != nil || exe == "" {
		return ""
	}
	if !strings.ContainsRune(exe, os.PathSeparator) {
		return exe
	}
	return path.Dir(exe)
})

type statfsInfoExtractor struct {
	simpleExtractor
	kind       FilesystemProperty
	pipestance bool
}

func (m *statfsInfoExtractor) MarshalJSON() ([]byte, error) {
	return json.Marshal(struct {
		Kind       string `json:"filesystem"`
		Pipestance bool   `json:"pipestance,omitempty"`
	}{
		Kind:       m.kind.String(),
		Pipestance: m.pipestance,
	})
}

func (fs *statfsInfoExtractor) Start(ctx context.Context, vc ValueContext,
	fc *FileContext, result chan<- any) {
	var r <-chan unix.Statfs_t
	if fs.pipestance {
		r = fc.statFs(vc.PipestanceDir())
	} else if exe := exeDir(); exe == "" {
		close(result)
		return
	} else {
		r = fc.statFs(path.Dir(exe))
	}
	go fs.get(ctx, r, result)
}

func (fs *statfsInfoExtractor) get(ctx context.Context,
	r <-chan unix.Statfs_t, result chan<- any) {
	defer close(result)
	select {
	case s := <-r:
		switch fs.kind {
		case FsAvail:
			result <- int64(s.Bavail) * s.Bsize
		case FsCapacity:
			result <- int64(s.Blocks) * s.Bsize
		case FsInodeAvail:
			result <- int64(s.Ffree)
		case FsInodeCap:
			result <- int64(s.Files)
		}
	case <-ctx.Done():
	}
}

func (fs *statfsInfoExtractor) RequiresPipeline() bool {
	return fs.pipestance
}

type mountInfoExtractor struct {
	simpleExtractor
	kind       FilesystemProperty
	pipestance bool
}

func (m *mountInfoExtractor) MarshalJSON() ([]byte, error) {
	return json.Marshal(struct {
		Kind       string `json:"filesystem"`
		Pipestance bool   `json:"pipestance,omitempty"`
	}{
		Kind:       m.kind.String(),
		Pipestance: m.pipestance,
	})
}

func (fs *mountInfoExtractor) Start(ctx context.Context, vc ValueContext,
	fc *FileContext, result chan<- any) {
	var r <-chan mountInfo
	if fs.pipestance {
		r = fc.mountInfo(vc.PipestanceDir())
	} else if exe := exeDir(); exe == "" {
		close(result)
		return
	} else {
		r = fc.mountInfo(exe)
	}
	go fs.get(ctx, r, result)
}

func (fs *mountInfoExtractor) get(ctx context.Context, r <-chan mountInfo,
	result chan<- any) {
	defer close(result)
	select {
	case s := <-r:
		switch fs.kind {
		case FsKind:
			result <- s.mType
		case FsOpts:
			result <- s.opts
		}
	case <-ctx.Done():
	}
}

func (fs *mountInfoExtractor) RequiresPipeline() bool {
	return fs.pipestance
}

func MakeFilesystem(pipestance bool, prop FilesystemProperty) ValueExtractor {
	switch prop {
	case FsKind, FsOpts:
		return &mountInfoExtractor{
			kind:       prop,
			pipestance: pipestance,
		}
	case FsAvail, FsInodeAvail, FsCapacity, FsInodeCap:
		return &statfsInfoExtractor{
			kind:       prop,
			pipestance: pipestance,
		}
	}
	panic("Invalid filesystem property")
}
