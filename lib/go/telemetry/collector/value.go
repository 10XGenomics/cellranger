package collector

import (
	"context"
	"encoding/json"
	"time"
)

// This allows e.g. command line flags to be passed in to
// the value extractor.
type ValueContext interface {
	CommandLine() []string
	MroFlags() []string
	PipestanceDir() string
	Runtime() time.Duration
	MrpHasRun() bool
	Success() bool
	RecentCount(group string) (int, bool)
}

// A source of metric data.
//
// This gets a little bit complicated because of performance requirements.
// We need to be able to collect multiple values in parallel.
// It is also going to be fairly common for the source of multiple
// metric values to involve shared work,
// e.g. different keys from the same json file,
// or different files in the same cgroup controller.
// We also need to support timeouts.
// There is also context that needs to be passed down from
// the top level, e.g. command line flags, pipestance location,
// etc.
type ValueExtractor interface {
	json.Marshaler

	// Begin the search for value data.
	//
	// Implementations MUST close the result channel,
	// either synchronously or from a goroutine they start.
	//
	// Implementations SHOULD supply at most one value as
	// a result on the channel.
	// They MAY supply that result synchronously.
	// It therefore follows that callers SHOULD supply a channel
	// with a capacity of one.
	// If they do not, then there MUST be some other goroutine
	// reading from the channel.
	//
	// Implementations SHOULD start a goroutine to produce
	// results and return them on the given channel if producing
	// the result may block, e.g. if making a syscall.
	//
	// Implementations MUST NOT call any methods on the
	// `FileContext` after it returns (e.g. from a
	// goroutine that it starts).
	//
	// Implementations MUST NOT synchronously block on
	// any result channels returned by the `FileContext`.
	Start(ctx context.Context, vc ValueContext,
		fc *FileContext, result chan<- any)

	// If true, this collector is only valid when a pipestance has completed.
	RequiresPipeline() bool
	// If true, this collector is only valid when a pipestance
	// has completed successfully.
	RequiresSuccess() bool
	// If true, this collector is only valid when a pipestance
	// has completed with a failure.
	RequiresFailure() bool
	// If non-empty, this collector expects the metric group(s)
	// that it participates in to have the given glob keys defined.
	RequiredGlobKeys() []string

	withGlob(map[string]string) ValueExtractor

	// Returns true if buckets are optional for this metric.
	BucketsOptional() bool
}

func WithGlob(ve ValueExtractor, globs map[string]string) ValueExtractor {
	if r := ve.withGlob(globs); r != nil {
		return r
	}
	return ve
}
