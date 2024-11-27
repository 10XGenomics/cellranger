package config

import (
	"bytes"
	"encoding/json"
	"sort"
)

// In a few places, for convenience we allow
// either a list or a string, where we want the string
// to be interpreted as a singleton list.
type StringOrList []string

func (sl *StringOrList) UnmarshalJSON(b []byte) error {
	if len(b) == 0 {
		*sl = nil
		return nil
	} else if b[0] == '"' {
		*sl = []string{""}
		return json.Unmarshal(b, &(*sl)[0])
	} else {
		var ll []string
		err := json.Unmarshal(b, &ll)
		*sl = ll
		return err
	}
}

// marshalSorted marshals a string-keyed map to json,
// with the keys in sorted order.
func marshalSorted[T any](sm map[string]T) ([]byte, error) {
	if sm == nil {
		return []byte("null"), nil
	}
	if len(sm) == 0 {
		return []byte("{}"), nil
	}
	keys := make([]string, 0, len(sm))
	for k := range sm {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	var buf bytes.Buffer
	buf.WriteByte('{')
	for i, k := range keys {
		if i != 0 {
			buf.WriteByte(',')
		}
		if b, err := json.Marshal(k); err != nil {
			return nil, err
		} else {
			buf.Write(b)
		}
		buf.WriteByte(':')
		if b, err := json.Marshal(sm[k]); err != nil {
			return nil, err
		} else {
			buf.Write(b)
		}
	}
	buf.WriteByte('}')
	return buf.Bytes(), nil
}

// Just a map[string]string but serializes
// with keys in sorted order for reproducibility.
type SortedStringMap map[string]string

func (sm SortedStringMap) MarshalJSON() ([]byte, error) {
	return marshalSorted(sm)
}

// Just a map[string]*Metric but serializes
// with keys in sorted order for reproducibility.
type SortedMetricMap map[string]*Metric

func (sm SortedMetricMap) MarshalJSON() ([]byte, error) {
	return marshalSorted(sm)
}

// Add adds a metric to the map.  It returns true
// if the metric was added, or false if a metric
// with that name was already present.
func (sm SortedMetricMap) Add(m *Metric) bool {
	if om := sm[m.Name]; om != nil {
		return false
	}
	sm[m.Name] = m
	return true
}

type SortedAnyMap map[string]any

func (sm SortedAnyMap) MarshalJSON() ([]byte, error) {
	return marshalSorted(sm)
}
