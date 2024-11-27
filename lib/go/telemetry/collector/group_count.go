package collector

import "context"

// MakeRecentCounter returns a ValueExtractor to count recent events
// from the given group.
func MakeRecentCounter(group string) ValueExtractor {
	return &recentCountExtractor{groupName: group}
}

type recentCountExtractor struct {
	simpleExtractor
	groupName string
}

func (recentCountExtractor) MarshalJSON() ([]byte, error) {
	return []byte(`{"special":"recent_count"}`), nil
}

func (f *recentCountExtractor) Start(ctx context.Context,
	vc ValueContext, _ *FileContext, result chan<- any) {
	go func(result chan<- any, name string) {
		defer close(result)
		c, ok := vc.RecentCount(name)
		if ok {
			result <- c
		}
	}(result, f.groupName)
}
