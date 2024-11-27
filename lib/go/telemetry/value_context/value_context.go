// Package value_context parses the command line for the telemetry tool.
//
// The results are used to provide context for collection of telemetry
// information.
package value_context

import (
	"flag"
	"os"
	"strings"
	"sync"
	"time"

	"github.com/10XDev/cellranger/lib/go/telemetry/collector"
	"github.com/martian-lang/martian/martian/util"
)

type CollectFlags struct {
	args     []string
	mroflags string
	walltime float64
	exitCode int
	Verbose  bool
}

func ParseFlags(args []string) CollectFlags {
	flagset := flag.NewFlagSet("telemetry collect",
		flag.ExitOnError)
	var f CollectFlags
	flagset.Float64Var(&f.walltime, "walltime", 0,
		"Wall time for the pipestance run, if any")
	flagset.StringVar(&f.mroflags, "mrp_flags", "",
		"Flags passed to mrp")
	flagset.IntVar(&f.exitCode, "mrp_exit_code", -1,
		"Exit code from mrp.")
	flagset.BoolVar(&f.Verbose, "verbose", false,
		"Log error messages")
	if err := flagset.Parse(args); err != nil {
		panic(err)
	}
	f.args = flagset.Args()
	return f
}

type valueContext struct {
	eventCounts func() map[string]int
	args        []string
	mroflags    []string
	pipestance  string
	runtime     time.Duration
	success     bool
	mrpRan      bool
}

func (f CollectFlags) MakeValueContext(
	eventCounts func() map[string]int,
) (collector.ValueContext, error) {
	vc := new(valueContext)
	vc.args = f.args

	vc.mroflags = strings.Fields(f.mroflags)
	// Look for the pipestance path.  It'll be either the `-psdir` flag or it'll
	// be the last argument.
	var foundPsdir bool
	for _, f := range vc.mroflags {
		if after, found := strings.CutPrefix(f, "--psdir="); found {
			foundPsdir = true
			vc.pipestance = after
			break
		} else if !strings.HasPrefix(f, "-") && util.ValidateID(f) == nil {
			vc.pipestance = f
		}
	}
	if !foundPsdir {
		// `--psdir` may be set in `MROFLAGS`, which would override the
		// pipestance ID.
		for _, f := range strings.Fields(os.Getenv("MROFLAGS")) {
			if after, found := strings.CutPrefix(f, "--psdir="); found {
				vc.pipestance = after
				break
			}
		}
	}

	vc.runtime = time.Duration(float64(time.Second) * f.walltime)
	vc.success = f.exitCode == 0
	vc.mrpRan = f.exitCode >= 0
	vc.eventCounts = sync.OnceValue(eventCounts)
	return vc, nil
}

// CommandLine implements collector.ValueContext.
func (v *valueContext) CommandLine() []string {
	return v.args
}

// MroFlags implements collector.ValueContext.
func (v *valueContext) MroFlags() []string {
	return v.mroflags
}

// PipestanceDir implements collector.ValueContext.
func (v *valueContext) PipestanceDir() string {
	if v.mrpRan {
		return v.pipestance
	} else {
		// There's nothing useful to be done with the pipestance directory
		// in this case.
		return ""
	}
}

// Runtime implements collector.ValueContext.
func (v *valueContext) Runtime() time.Duration {
	return v.runtime
}

// Success implements collector.ValueContext.
func (v *valueContext) Success() bool {
	return v.success
}

func (v *valueContext) MrpHasRun() bool {
	return v.mrpRan
}

// RecentCount implements collector.ValueContext.
func (v *valueContext) RecentCount(group string) (int, bool) {
	r, ok := v.eventCounts()[group]
	return r, ok
}
