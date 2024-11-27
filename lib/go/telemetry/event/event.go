// Package event defines a telemetry event type and code for collecting them.
package event

import (
	"context"
	"fmt"
	"io"
	"io/fs"
	"maps"
	"os"
	"runtime/trace"
	"sort"
	"sync"

	"github.com/10XDev/cellranger/lib/go/telemetry/collector"
	"github.com/10XDev/cellranger/lib/go/telemetry/config"
)

type Event struct {
	Properties    config.SortedAnyMap `json:"properties"`
	GroupHash     []byte              `json:"group_hash"`
	Name          string              `json:"name"`
	Product       string              `json:"product"`
	Version       string              `json:"version"`
	Internal      string              `json:"internal_user,omitempty"`
	ConfigVersion int                 `json:"config_version"`
	Replicas      int                 `json:"replicas,omitempty"`
}

type StringOrBytesWriter interface {
	io.Writer
	io.StringWriter
}

// Print outputs a human-readable presentation of the data.
func (ev *Event) Print(w StringOrBytesWriter) error {
	if ev.Replicas != 0 {
		if _, err := fmt.Fprintf(w, "(one of %d)\n", ev.Replicas); err != nil {
			return err
		}
	}
	if len(ev.Properties) == 0 {
		_, err := w.WriteString("(no data)\n")
		return err
	}
	keys := make([]string, 0, len(ev.Properties))
	for k := range ev.Properties {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	for _, key := range keys {
		if _, err := w.WriteString(key); err != nil {
			return err
		}
		if _, err := w.WriteString(":\t"); err != nil {
			return err
		}
		if v := ev.Properties[key]; v == nil {
			if _, err := w.WriteString("null\n"); err != nil {
				return err
			}
		} else if _, err := fmt.Fprintln(w, v); err != nil {
			return err
		}
	}
	return nil
}

func Collect(ctx context.Context, c *config.CompiledConfig,
	vc collector.ValueContext) <-chan Event {
	ctx, task := trace.NewTask(ctx, "Collect")
	defer trace.StartRegion(ctx, "beginCollect").End()
	results := make(chan Event, len(c.Groups))
	var doneWait sync.WaitGroup
	var fc collector.FileContext
	internalTag := os.Getenv("_TENX_TELEMETRY_INTERNAL")
	for gn := range c.Groups {
		collectGroup(ctx, &doneWait, c, gn, vc, &fc,
			internalTag, results)
	}
	fc.Close(ctx, vc.PipestanceDir())
	go func(wg *sync.WaitGroup, r chan<- Event) {
		wg.Wait()
		task.End()
		close(r)
	}(&doneWait, results)
	return results
}

// groupActive returns true if the given group is active given the context.
func groupActive(g *config.MetricGroup, vc collector.ValueContext) bool {
	switch g.When {
	case config.WhenCommand:
		if vc.MrpHasRun() {
			return false
		}
	case config.WhenSuccess:
		if vc.PipestanceDir() == "" || !vc.MrpHasRun() || !vc.Success() {
			return false
		}
	case config.WhenFailure:
		if vc.PipestanceDir() == "" || !vc.MrpHasRun() || vc.Success() {
			return false
		}
	}
	if len(g.Subcommands) != 0 {
		subcommand := vc.CommandLine()[0]
		for _, sc := range g.Subcommands {
			if subcommand == sc {
				return true
			}
		}
		return false
	}
	return true
}

func makeGlobs(ps string, globs map[string]string) []map[string]string {
	psfs := os.DirFS(ps)
	result := make([]map[string]string, 0, 2*len(globs))
	for key, val := range globs {
		matches, err := fs.Glob(psfs, val)
		if err != nil {
			// This should have been caught in config validation, unless the
			// pipestance directory itself is an invalid glob.
			continue
		}
		if len(matches) == 0 {
			// If there's no match, we have a few choices.
			// We could just not run this group at all, but that's not great
			// because there may be metrics which don't depend on the missing
			// key.
			// We could use an empty string, except that much of the time
			// the glob is going to be used in expressions like `${key}/foo`,
			// so using `.` is safer.
			if len(result) == 0 {
				result = append(result, map[string]string{key: "."})
			} else {
				for _, m := range result {
					m[key] = "."
				}
			}
			continue
		}
		if len(result) == 0 {
			for _, match := range matches {
				result = append(result, map[string]string{
					key: match,
				})
			}
		} else {
			for _, m := range result {
				m[key] = matches[0]
			}
			existing := result
			for _, match := range matches[1:] {
				for _, m := range existing {
					mc := maps.Clone(m)
					mc[key] = match
					result = append(result, mc)
				}
			}
		}
	}
	return result
}

func collectGroup(ctx context.Context,
	doneWait *sync.WaitGroup,
	c *config.CompiledConfig, groupName string,
	vc collector.ValueContext,
	fc *collector.FileContext,
	internalTag string,
	results chan<- Event) {
	g := c.Groups[groupName]
	if g == nil {
		panic("invalid group name")
	}
	if !groupActive(g, vc) {
		return
	}
	hash := g.Hash()
	ev := Event{
		GroupHash:     hash[:],
		Name:          groupName,
		ConfigVersion: c.ConfigVersion,
		Product:       c.Product,
		Version:       c.Version,
		Internal:      internalTag,
	}
	if len(g.Globs) == 0 {
		collectWithGlob(ctx, doneWait, g, nil, vc, fc, ev, results)
		return
	}
	globs := makeGlobs(vc.PipestanceDir(), g.Globs)
	ev.Replicas = len(globs)
	for _, glob := range globs {
		collectWithGlob(ctx, doneWait, g, glob, vc, fc, ev, results)
	}
}

func collectWithGlob(ctx context.Context,
	doneWait *sync.WaitGroup,
	g *config.MetricGroup,
	glob map[string]string,
	vc collector.ValueContext,
	fc *collector.FileContext,
	ev Event, results chan<- Event) {
	ev.Properties = make(config.SortedAnyMap, len(g.Metrics))
	metrics := g.Metrics
	if len(glob) > 0 {
		metrics = make(config.SortedMetricMap, len(g.Metrics))
		for n, m := range g.Metrics {
			metrics[n] = m.WithGlob(glob)
		}
	}
	for _, m := range metrics {
		m.Start(ctx, vc, fc)
	}
	doneWait.Add(1)
	go func(ev Event,
		doneWait *sync.WaitGroup, results chan<- Event) {
		defer doneWait.Done()
		for n, m := range metrics {
			ev.Properties[n] = m.Get()
		}
		results <- ev
	}(ev, doneWait, results)
}
