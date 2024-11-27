package collector

import (
	"bytes"
	"context"
	"io"
	"os"
	"path"
	"path/filepath"
	"strings"

	"github.com/google/shlex"
)

type distroExtractor struct{ noBucketsExtractor }

func (distroExtractor) MarshalJSON() ([]byte, error) {
	return []byte(`{"special":"distro"}`), nil
}

func (d *distroExtractor) Start(ctx context.Context,
	_ ValueContext, _ *FileContext, result chan<- any) {
	go d.get(ctx, result)
}

func openOsReleaseFile() *os.File {
	etc := os.Getenv("UNIXCONFDIR")
	if etc == "" {
		etc = "/etc"
	}
	usrLib := os.Getenv("UNIXUSRLIBDIR")
	if usrLib == "" {
		usrLib = "/usr/lib"
	}
	f, err := os.Open(path.Join(etc, "os-release"))
	if err == nil {
		return f
	}
	f, err = os.Open(path.Join(usrLib, "os-release"))
	if err == nil {
		return f
	}
	if etc != "/etc" {
		f, err = os.Open("/etc/os-release")
		if err == nil {
			return f
		}
	}
	if usrLib != "/usr/lib" {
		f, err = os.Open("os-release")
		if err == nil {
			return f
		}
	}
	f, err = os.Open("/lib/os-release")
	if err == nil {
		return f
	}
	f, err = os.Open("/etc/lsb-release")
	if err == nil {
		return f
	}
	return nil
}

func parseOsRelease(f io.Reader) string {
	lex := shlex.NewLexer(f)
	var id, name string
	for word, err := lex.Next(); err == nil; word, err = lex.Next() {
		if word == "" {
			continue
		}
		key, val, ok := strings.Cut(word, "=")
		if !ok || val == "" {
			continue
		}
		if strings.EqualFold(key, "ID_LIKE") {
			for _, d := range strings.Fields(val) {
				if n := normalizeDistro(d); n != "" {
					return n
				}
			}
		} else if strings.EqualFold(key, "ID") {
			id = val
		} else if strings.EqualFold(key, "NAME") {
			name = val
		} else if strings.EqualFold(key, "DISTRIB_ID") && id == "" {
			id = val
		}
	}
	if id != "" {
		if n := normalizeDistro(id); n != "" {
			return n
		}
		id, _, _ = strings.Cut(id, " ")
		if n := normalizeDistro(id); n != "" {
			return n
		}
	}
	if n := normalizeDistro(name); n != "" {
		return n
	}
	name, _, _ = strings.Cut(name, " ")
	if n := normalizeDistro(name); n != "" {
		return n
	}
	return "unknown"
}

func getDistroFromOsRelease() string {
	f := openOsReleaseFile()
	if f == nil {
		return ""
	}
	defer f.Close()
	return parseOsRelease(f)
}

func normalizeDistro(distro string) string {
	distro = strings.ToLower(distro)
	// Some distros can't be ID'ed from the first word
	switch distro {
	case "linux mint":
		return "debian"
	case "ibm powerkvm":
		return "kvm"
	case "debian", "devuan", "ubuntu", "raspbian":
		return "debian"
	case "redhat", "rhel", "rocky", "almalinux",
		"alt", "altlinux",
		"centos", "fedora",
		"amazon", "amzn",
		"cloudlinux", "pidora",
		"oracle",
		"springdale", "scientific", "red":
		return "rhel"
	case "suse", "openSUSE":
		return "suse"
	case "arch":
		return "arch"
	case "gentoo", "exherbo":
		return "gentoo"
	case "kvm":
		return "kvm"
	case "mandriva", "mageia":
		return "mandriva"
	case "parallels":
		return "parallels"
	case "xenserver":
		return "xenserver"
	case "slackware":
		return "slackware"
	case "guix":
		return "guix"
	case "nixos":
		return "nixos"
	}
	return ""
}

func (d *distroExtractor) get(ctx context.Context, result chan<- any) {
	defer close(result)

	if distro := getDistroFromOsRelease(); distro != "" {
		result <- distro
		return
	}
	if ctx.Err() != nil {
		return
	}
	releaseFiles, _ := filepath.Glob("/etc/*-release")
	for _, fn := range releaseFiles {
		if n := normalizeDistro(strings.TrimSuffix(
			strings.TrimPrefix(fn, "/etc/"), "-release")); n != "" {
			result <- n
			return
		}
	}
}

type containerExtractor struct{ noBucketsExtractor }

func (containerExtractor) MarshalJSON() ([]byte, error) {
	return []byte(`{"special":"container"}`), nil
}

func (d *containerExtractor) Start(ctx context.Context,
	_ ValueContext, _ *FileContext, result chan<- any) {
	go d.get(ctx, result)
}

func (d *containerExtractor) get(ctx context.Context, result chan<- any) {
	defer close(result)
	if ctx.Err() != nil {
		return
	}
	content, _ := os.ReadFile("/proc/1/sched")
	if len(content) > 0 {
		content, _, _ := bytes.Cut(content, []byte{'\n'})
		content, _, _ = bytes.Cut(content, []byte{'('})
		content, _, _ = bytes.Cut(bytes.TrimSpace(content), []byte{' '})
		if string(content) == "init" ||
			string(content) == "systemd" ||
			string(content) == "system" {
			result <- false
		} else {
			result <- true
		}
	}
}
