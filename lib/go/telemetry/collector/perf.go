package collector

import (
	"context"
	"encoding/json"
	"fmt"
	"sort"

	"github.com/martian-lang/martian/martian/core"
	"github.com/martian-lang/martian/martian/syntax"
)

// A collector that parses the _perf file and emits top N stage metrics.
//
// Reports the top 5 stages in terms of wall time (longest fork), RSS
// (largest fork), and CPU time (total across all forks).
type perfStatsCollector struct{ noBucketsExtractor }

func (perfStatsCollector) MarshalJSON() ([]byte, error) {
	return []byte(`{"special":"perf_stats"}`), nil
}

func (perfStatsCollector) RequiresPipeline() bool { return true }
func (perfStatsCollector) RequiresSuccess() bool  { return true }

func (perfStatsCollector) IsMultiValue() bool { return true }

func (perfStatsCollector) Start(ctx context.Context, _ ValueContext,
	fc *FileContext, result chan<- any) {
	if fc.files == nil {
		fc.files = make(multiFileResultMap)
	}
	bytesCh := fc.files.getSink([]string{"_perf"}).subscribeRaw()
	go func() {
		defer close(result)
		select {
		case data, ok := <-bytesCh:
			if !ok {
				return
			}
			if stats, err := computePerfStats(data); err == nil {
				result <- stats
			}
		case <-ctx.Done():
		}
	}()
}

type stagePerfSummary struct {
	name     string
	wallTime float64
	cpuTime  float64
	rss      int
}

const TOP_N = 5

// Gather perf stats, sort by top N, and report as a collection.
//
// Pairs of (stage name, value) are reported as two entries with related keys:
// top_foo_1_stage: STAGE_NAME, top_foo_1_val: 12345.
func computePerfStats(data []byte) (map[string]any, error) {
	var nodes []*core.NodePerfInfo
	if err := json.Unmarshal(data, &nodes); err != nil {
		return nil, err
	}
	stages := make([]stagePerfSummary, 0, len(nodes))
	for _, node := range nodes {
		if node == nil || node.Type != syntax.KindStage {
			continue
		}
		var wallTime float64
		var cpuTime float64
		var rss int
		for _, fork := range node.Forks {
			if fork == nil || fork.ForkStats == nil {
				continue
			}
			fs := fork.ForkStats
			wallTime = max(wallTime, fs.WallTime)
			cpuTime += fs.UserTime + fs.SystemTime
			rss = max(rss, fs.MaxRss)
		}
		stages = append(stages, stagePerfSummary{
			name:     node.Name,
			wallTime: wallTime,
			cpuTime:  cpuTime,
			rss:      rss,
		})
	}
	// TOP_N vals, two entries for each (stage name and val), 3 categories.
	result := make(map[string]any, TOP_N*2*3)
	emitTop(result, stages, "wall_time",
		func(s stagePerfSummary) float64 { return s.wallTime },
		func(s stagePerfSummary) any { return s.wallTime })
	emitTop(result, stages, "cpu_time",
		func(s stagePerfSummary) float64 { return s.cpuTime },
		func(s stagePerfSummary) any { return s.cpuTime })
	emitTop(result, stages, "rss",
		func(s stagePerfSummary) float64 { return float64(s.rss) },
		func(s stagePerfSummary) any { return s.rss })
	return result, nil
}

// Given a stage perf summary and getter functions, emit top N.
func emitTop(result map[string]any, stages []stagePerfSummary,
	metric string,
	key func(stagePerfSummary) float64,
	val func(stagePerfSummary) any) {
	sorted := make([]stagePerfSummary, len(stages))
	copy(sorted, stages)
	sort.Slice(sorted, func(i, j int) bool {
		return key(sorted[i]) > key(sorted[j])
	})
	n := min(TOP_N, len(sorted))
	for i := range n {
		prefix := fmt.Sprintf("top_%s_%d", metric, i+1)
		result[prefix+"_stage"] = sorted[i].name
		result[prefix+"_val"] = val(sorted[i])
	}
}
