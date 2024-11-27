package config

import (
	"errors"
	"fmt"
	"path"
	"strings"
	"unicode"
	"unicode/utf8"

	"github.com/10XDev/cellranger/lib/go/telemetry/collector"
	"github.com/10XDev/cellranger/lib/go/telemetry/converter"
	"golang.org/x/sys/unix"
)

func getSpecialSource(name string) (ValueExtractor, error) {
	return collector.MakeSpecial(name)
}

func (fs *FilesystemValue) compile() (ValueExtractor, error) {
	var prop collector.FilesystemProperty
	switch fs.Property {
	case "type":
		prop = collector.FsKind
	case "options":
		prop = collector.FsOpts
	case "available_bytes":
		prop = collector.FsAvail
	case "available_inodes":
		prop = collector.FsInodeAvail
	case "total_bytes":
		prop = collector.FsCapacity
	case "total_inodes":
		prop = collector.FsInodeCap
	default:
		return nil, fmt.Errorf("invalid property %q", fs.Property)
	}
	switch fs.Which {
	case "bins":
		return collector.MakeFilesystem(false, prop), nil
	case "pipestance":
		return collector.MakeFilesystem(true, prop), nil
	default:
		return nil, fmt.Errorf("invalid filesystem kind %q", fs.Which)
	}
}

func (fc *FileContentConfig) lineFilter() (*collector.LineFilter, error) {
	line := fc.Line
	if (fc.LastLine == "") == (fc.Line == "") {
		return nil, errors.New(
			"exactly one of line or last_line must be specified")
	}
	if line == "" {
		line = fc.LastLine
	}
	return collector.MakeLineFilter(line, fc.LastLine != "")
}

func (fc *FileContentConfig) compile() (ValueExtractor, error) {
	conv, err := converter.Make(fc.Converter)
	if err != nil {
		return nil, err
	}
	if fc.Cgroup != "" {
		if len(fc.Path) == 0 {
			return nil, errors.New("cgroup requres path")
		}
		if fc.Flag != "" {
			return nil, errors.New(
				"exactly one of flag or path may be specified")
		}
		if fc.StageRel {
			return nil, errors.New("cgroup cannot be stage relative")
		}
		for _, p := range fc.Path {
			if p == "" || p == "." || p == ".." || strings.ContainsAny(p, `/\`) {
				return nil, fmt.Errorf("invalid cgroup path %q", p)
			}
		}
		if len(fc.JsonPath) != 0 {
			return nil, errors.New("cgroup files are not json")
		}
		if conv != nil {
			return nil, errors.New("no converters with cgroup files")
		}
		filter, err := fc.lineFilter()
		if err != nil {
			return nil, err
		}
		return collector.Cgroup(fc.Cgroup, fc.Path, filter), nil
	}
	var filter collector.FileFilter
	if len(fc.JsonPath) != 0 {
		if fc.Line != "" || fc.LastLine != "" {
			return nil, errors.New("cannot specify both json_path and line")
		}
		filter = collector.JsonPathFilter(fc.JsonPath)
	} else if fc.Line == "" && fc.LastLine == "" {
		return nil, errors.New("must specify exactly one of line, last_line, or json_path")
	} else {
		filter, err = fc.lineFilter()
		if err != nil {
			return nil, err
		}
	}
	if fc.Flag != "" {
		if len(fc.Path) != 0 {
			return nil, errors.New(
				"exactly one of flag or path may be specified")
		}
		if err := validateFlag(fc.Flag); err != nil {
			return nil, err
		}
		if fc.StageRel {
			return nil, errors.New("flag-specified file cannot be stage relative")
		}
		return collector.MakeFileFlag(fc.Flag, conv, filter), nil
	}
	if len(fc.Path) == 0 {
		return nil, errors.New("must specify a source for file_content")
	}
	var isAbs bool
	for i, p := range fc.Path {
		p = path.Clean(p)
		if p == "." {
			return nil, fmt.Errorf("empty path %q",
				fc.Path[i])
		}
		if p == ".." || strings.HasPrefix(p, "../") {
			return nil, fmt.Errorf(
				"path %q is outside of the pipestance directory",
				p)
		}
		if path.IsAbs(p) {
			if !strings.HasPrefix(p, "/proc/self/") &&
				!strings.HasPrefix(p, "/proc/sys/") &&
				!strings.HasPrefix(p, "/sys/kernel/") &&
				p != "/proc/cpuinfo" && p != "/proc/meminfo" {
				return nil, fmt.Errorf("disallowed absolute path %q", p)
			}
			if i > 0 && !isAbs {
				return nil, errors.New("mixed absolute and relative paths")
			}
			if fc.StageRel {
				return nil, errors.New("absolute path cannot be stage relative")
			}
			isAbs = true
		} else if isAbs {
			return nil, errors.New("mixed absolute and relative paths")
		}
		fc.Path[i] = p
	}
	if isAbs {
		return collector.MakeAbsFile(fc.Path, conv, filter), nil
	}
	return collector.MakeRelativeFile(
		fc.Path, fc.StageRel, conv, filter)
}

func (fs *FileSizeConfig) compile(count bool) (ValueExtractor, error) {
	if (fs.Flag == "") == (len(fs.Path) == 0) {
		return nil, errors.New("exactly one of flag or path must be specified")
	}
	if fs.Flag != "" {
		if err := validateFlag(fs.Flag); err != nil {
			return nil, err
		}
		return nil, errors.ErrUnsupported
		// return collector.MakeFileSizeFlag(fs.Flag, count), nil
	}
	for i, p := range fs.Path {
		p = path.Clean(p)
		fs.Path[i] = p
		if p == "." {
			return nil, errors.New("empty path")
		}
		if p == ".." || p[0] == '/' || strings.HasPrefix(p, "../") {
			return nil, fmt.Errorf(
				"%q is outside of the pipestance directory", p)
		}
	}
	if count {
		return nil, fmt.Errorf("%w: file count not yet implemented",
			errors.ErrUnsupported)
	} else {
		return nil, fmt.Errorf("%w: file size not yet implemented",
			errors.ErrUnsupported)
	}
	// return collector.MakeFileSize(fs.Path, count), nil
}

func parseRlimit(res string) (int, error) {
	switch res {
	case "as":
		return unix.RLIMIT_AS, nil
	case "core":
		return unix.RLIMIT_CORE, nil
	case "cpu":
		return unix.RLIMIT_CPU, nil
	case "data":
		return unix.RLIMIT_DATA, nil
	case "fsize":
		return unix.RLIMIT_FSIZE, nil
	case "memlock":
		return unix.RLIMIT_MEMLOCK, nil
	case "msgqueue":
		return unix.RLIMIT_MSGQUEUE, nil
	case "nice":
		return unix.RLIMIT_NICE, nil
	case "nofile":
		return unix.RLIMIT_NOFILE, nil
	case "nproc":
		return unix.RLIMIT_NPROC, nil
	case "rttime":
		return unix.RLIMIT_RTTIME, nil
	case "sigpending":
		return unix.RLIMIT_SIGPENDING, nil
	case "stack":
		return unix.RLIMIT_STACK, nil
	}
	return 0, fmt.Errorf("unrecognized rlimit resource %s",
		res)
}

func (rlim *RlimitValue) compile() (ValueExtractor, error) {
	if res, err := parseRlimit(rlim.Resource); err != nil {
		return nil, err
	} else {
		return collector.MakeRlimit(res, rlim.Hard), nil
	}
}

func compileFlagSource(flag string, mro bool) (ValueExtractor, error) {
	if err := validateFlag(flag); err != nil {
		return nil, err
	}
	return collector.MakeFlag(flag, mro), nil
}

func validateFlag(flag string) error {
	if !strings.HasPrefix(flag, "--") {
		return fmt.Errorf("flag %q does not start with --", flag)
	}
	if len(flag) < 3 {
		return errors.New("empty flag name")
	}
	if strings.ContainsFunc(flag[2:], func(r rune) bool {
		// flag names should only contain ascii letters and digits
		// or - or _ as word separators.
		return r > utf8.RuneSelf ||
			(!unicode.IsLetter(r) && !unicode.IsDigit(r) &&
				r != '-' && r != '_')
	}) {
		return fmt.Errorf("invalid flag name %s", flag)
	}
	return nil
}

func compileRecentCountSource(group string, groups []string) (ValueExtractor, error) {
	for _, g := range groups {
		if group == g {
			return collector.MakeRecentCounter(group), nil
		}
	}
	return nil, fmt.Errorf(
		"recent counter for group %s must be included in a group of that name",
		group)
}
