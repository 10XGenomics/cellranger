package collector

import (
	"bytes"
	"context"
	"encoding/json"
	"os"
	"path"
	"runtime/trace"
)

func doCgroups(ctx context.Context,
	requests map[string][]cgroupRequest, mounts <-chan map[string]string) {
	defer trace.StartRegion(ctx, "doCgroups").End()
	groupFile, err := os.ReadFile("/proc/self/cgroup")
	if err != nil {
		for _, cs := range requests {
			for _, c := range cs {
				close(c.result)
			}
		}
		return
	}
	var cgMounts map[string]string
	select {
	case cgMounts = <-mounts:
	case <-ctx.Done():
		for _, cs := range requests {
			for _, c := range cs {
				close(c.result)
			}
		}
		return
	}
	if len(cgMounts) == 0 {
		for _, cs := range requests {
			for _, c := range cs {
				close(c.result)
			}
		}
		return
	}
	var v2 string
	dirs := make(map[string]string, len(requests))
	for _, line := range bytes.Split(
		bytes.TrimSpace(groupFile),
		[]byte{'\n'}) {
		_, after, ok := bytes.Cut(line, []byte{':'})
		if !ok || len(after) == 0 {
			continue
		}
		if after[0] == ':' {
			if p := cgMounts[""]; p == "" {
				continue
			} else {
				v2 = p + string(after[1:])
			}
		}
		ctrl, pb, ok := bytes.Cut(after, []byte{':'})
		if !ok {
			continue
		}
		p := string(pb)
		for _, c := range bytes.Split(ctrl, []byte{','}) {
			if _, ok := requests[string(c)]; !ok {
				continue
			}
			if m := cgMounts[string(c)]; m != "" {
				dirs[string(c)] = m + p
			}
		}
	}
	for ctrl, rs := range requests {
		p := dirs[ctrl]
		if p == "" {
			if v2 != "" {
				go doCgReqs(ctx, v2, rs)
			} else {
				for _, r := range rs {
					close(r.result)
				}
			}
		} else {
			go doCgReqs(ctx, p, rs)
		}
	}
}

func doCgReqs(ctx context.Context, dir string, reqs []cgroupRequest) {
	for _, req := range reqs {
		doCgReq(ctx, dir, req)
	}
}

func doCgReq(ctx context.Context, dir string, req cgroupRequest) {
	defer close(req.result)
	if ctx.Err() != nil {
		return
	}
	for _, fn := range req.names {
		if c, err := os.ReadFile(path.Join(dir, fn)); err == nil {
			req.result <- c
			return
		}
	}
}

type cgroupExtractor struct {
	simpleExtractor
	filter *LineFilter
	paths  []string
	ctrl   string
}

func (c *cgroupExtractor) Start(ctx context.Context, _ ValueContext,
	fc *FileContext, result chan<- any) {
	cr := fc.cgroup(c.ctrl, c.paths)
	go c.get(ctx, cr, result)
}

func (c *cgroupExtractor) get(ctx context.Context, cr <-chan []byte, result chan<- any) {
	defer close(result)
	select {
	case content := <-cr:
		if r, ok := c.filter.compile().matchContent(string(content)); ok {
			result <- r
		}
	case <-ctx.Done():
	}
}

func (c *cgroupExtractor) MarshalJSON() ([]byte, error) {
	b, err := json.Marshal(struct {
		Filter FileFilter `json:"filter"`
		Paths  []string   `json:"paths"`
		Ctrl   string     `json:"controller"`
	}{
		Filter: c.filter,
		Paths:  c.paths,
		Ctrl:   c.ctrl,
	})
	if err != nil {
		return nil, err
	}
	return append(append([]byte(`{"cgroup":`), b...), '}'), nil
}

func Cgroup(ctrl string, paths []string, filter *LineFilter) ValueExtractor {
	return &cgroupExtractor{
		filter: filter,
		paths:  paths,
		ctrl:   ctrl,
	}
}
