package collector

import (
	"fmt"
	"regexp"
	"strings"
)

type FileFilter interface {
	// We've only got two implementations of this
	// interface, and we want to keep them private to
	// this package, because to be efficient, they
	// have to be special-cased regardless.
	isJson() bool
}

type JsonPathFilter []string

func (JsonPathFilter) isJson() bool {
	return true
}

type LineFilter struct {
	Filter string `json:"filter"`
	Last   bool   `json:"last,omitempty"`
}

func (LineFilter) isJson() bool {
	return false
}

func (lf LineFilter) compile() compiledLineFilter {
	// We store the filter string rather than the compiled regexp
	// because we want to be able to use `LineFilter` as a map key.
	return compiledLineFilter{
		filter: regexp.MustCompile(lf.Filter),
		last:   lf.Last,
	}
}

type compiledLineFilter struct {
	filter *regexp.Regexp
	last   bool
}

func (lf compiledLineFilter) matchContent(content string) (string, bool) {
	lines := strings.Split(content, "\n")
	if lf.last {
		for i := len(lines) - 1; i >= 0; i-- {
			sm := lf.filter.FindStringSubmatch(lines[i])
			if len(sm) == 1 {
				return sm[0], true
			} else if len(sm) >= 2 {
				return sm[1], true
			}
		}
	} else {
		for _, line := range lines {
			sm := lf.filter.FindStringSubmatch(line)
			if len(sm) == 1 {
				return sm[0], true
			} else if len(sm) >= 2 {
				return sm[1], true
			}
		}
	}
	return "", false
}

func (lf compiledLineFilter) matchLine(line []byte) ([]byte, bool) {
	sm := lf.filter.FindSubmatch(line)
	if len(sm) == 1 {
		return sm[0], true
	} else if len(sm) >= 2 {
		return sm[1], true
	}
	return nil, false
}

func MakeLineFilter(filter string, last bool) (*LineFilter, error) {
	_, err := regexp.Compile(filter)
	if err != nil {
		return nil, fmt.Errorf(
			"invalid regex for line filter: %w",
			err)
	}
	return &LineFilter{
		Filter: filter,
		Last:   last,
	}, nil
}
