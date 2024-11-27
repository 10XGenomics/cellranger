package collector

import (
	"bytes"
	"context"
	"fmt"
	"path"
	"regexp"
	"strings"

	"golang.org/x/sys/unix"
)

func MakeSpecial(name string) (ValueExtractor, error) {
	switch name {
	case "glibc":
		return glibcExtractor{}, nil
	case "kernel":
		return kernelExtractor{}, nil
	case "subcommand":
		return new(subcommandExtractor), nil
	case "runtime":
		return new(runtimeExtractor), nil
	case "distro":
		return new(distroExtractor), nil
	case "container":
		return new(containerExtractor), nil
	case "success":
		return new(successExtractor), nil
	case "failed_stage":
		return new(failedStageExtractor), nil
	}
	return nil, fmt.Errorf("unknown special collector type %s", name)
}

// Unit struct to be embedded in othe extractor types to avoid re-implementing
// common functions all over the place.
type simpleExtractor struct{}

func (simpleExtractor) RequiredGlobKeys() []string {
	return nil
}

// Because functions from embedded types don't have access to the object
// they're embedded in, this can't just return that object by default,
// so we just return nil and let the public `WithGlob` function deal with it,
// to avoid needing to reimplement this everywhere.
func (simpleExtractor) withGlob(map[string]string) ValueExtractor {
	return nil
}

func (simpleExtractor) RequiresPipeline() bool {
	return false
}

func (simpleExtractor) RequiresSuccess() bool {
	return false
}

func (simpleExtractor) RequiresFailure() bool {
	return false
}

func (simpleExtractor) BucketsOptional() bool {
	return false
}

type noBucketsExtractor struct{ simpleExtractor }

func (noBucketsExtractor) BucketsOptional() bool {
	return true
}

type kernelExtractor struct{ simpleExtractor }

func (kernelExtractor) MarshalJSON() ([]byte, error) {
	return []byte(`{"special":"kernel"}`), nil
}

var semverRegex = regexp.MustCompile(`\d+\.\d+`)

func (kernelExtractor) Start(_ context.Context, _ ValueContext, _ *FileContext, result chan<- any) {
	go func(result chan<- any) {
		defer close(result)
		var buf unix.Utsname
		if unix.Uname(&buf) != nil {
			return
		}
		beforeNull, _, _ := bytes.Cut(buf.Release[:], []byte{0})
		if v := semverRegex.Find(beforeNull); len(v) > 0 {
			result <- string(v)
		}
	}(result)
}

type subcommandExtractor struct{ noBucketsExtractor }

func (subcommandExtractor) MarshalJSON() ([]byte, error) {
	return []byte(`{"special":"subcommand"}`), nil
}

func (subcommandExtractor) Start(_ context.Context,
	vc ValueContext, _ *FileContext, result chan<- any) {
	defer close(result)
	if cl := vc.CommandLine(); len(cl) > 0 {
		result <- cl[0]
	}
}

type runtimeExtractor struct{ simpleExtractor }

func (runtimeExtractor) MarshalJSON() ([]byte, error) {
	return []byte(`{"special":"runtime"}`), nil
}

func (runtimeExtractor) Start(_ context.Context,
	vc ValueContext, _ *FileContext, result chan<- any) {
	defer close(result)
	if rt := vc.Runtime(); rt != 0 {
		result <- rt.Seconds()
	}
}

func (runtimeExtractor) RequiresPipeline() bool {
	return true
}

type successExtractor struct{ noBucketsExtractor }

func (successExtractor) MarshalJSON() ([]byte, error) {
	return []byte(`{"special":"success"}`), nil
}

func (successExtractor) RequiresPipeline() bool {
	return true
}

// Start implements ValueExtractor.
func (successExtractor) Start(_ context.Context,
	vc ValueContext, _ *FileContext, result chan<- any) {
	defer close(result)
	result <- vc.Success()
}

type failedStageExtractor struct{ noBucketsExtractor }

func (failedStageExtractor) MarshalJSON() ([]byte, error) {
	return []byte(`{"special":"failed_stage"}`), nil
}

// Start implements ValueExtractor.
func (f *failedStageExtractor) Start(ctx context.Context,
	vc ValueContext, fc *FileContext, result chan<- any) {
	r := fc.errorFile()
	go func(result chan<- any, r <-chan errorFile) {
		defer close(result)
		select {
		case ef := <-r:
			result <- stagePathToFqid(ef.path)
		case <-ctx.Done():
		}
	}(result, r)
}

func stagePathToFqid(s string) string {
	return strings.ReplaceAll(path.Dir(path.Dir(s)), "/", ".")
}

func (failedStageExtractor) RequiresPipeline() bool {
	return true
}

func (failedStageExtractor) RequiresFailure() bool {
	return true
}
