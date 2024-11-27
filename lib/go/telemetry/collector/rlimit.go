package collector

import (
	"context"
	"encoding/json"

	"golang.org/x/sys/unix"
)

type rlimitCollector struct {
	simpleExtractor
	resource int
	hard     bool
}

func (r *rlimitCollector) MarshalJSON() ([]byte, error) {
	return json.Marshal(struct {
		Res  int  `json:"rlimit"`
		Hard bool `json:"hard,omitempty"`
	}{
		Res:  r.resource,
		Hard: r.hard,
	})
}

// Get implements ValueExtractor.
func (r *rlimitCollector) Start(_ context.Context,
	_ ValueContext, _ *FileContext, result chan<- any) {
	go func(result chan<- any) {
		defer close(result)
		var rlim unix.Rlimit
		if err := unix.Getrlimit(r.resource, &rlim); err != nil {
			return
		}
		if r.hard {
			result <- rlim.Max
		} else {
			result <- rlim.Cur
		}
	}(result)
}

func MakeRlimit(resource int, hard bool) ValueExtractor {
	return &rlimitCollector{
		resource: resource,
		hard:     hard,
	}
}
