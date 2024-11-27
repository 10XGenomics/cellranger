package config

import (
	"bytes"
	"cmp"
	"encoding/json"
	"errors"
	"fmt"
	"math"
	"regexp"
	"slices"
	"sort"
	"strconv"
	"strings"
)

func nTrue(b ...bool) int {
	n := 0
	for _, v := range b {
		if v {
			n++
		}
	}
	return n
}

func (pc *BucketSpec) compile(name string) (Bucketizer, error) {
	// Ensure we have exactly one option set.
	nOpt := nTrue(
		len(pc.Thresholds) > 0,
		pc.Range != nil,
		pc.Exponential != nil,
		len(pc.Semver) > 0,
		pc.Bool,
		pc.Float,
		len(pc.Match) > 0,
	)

	if nOpt > 1 {
		return nil, fmt.Errorf(
			"%w %s: more than one option set: %+v",
			BucketConfigurationError,
			name,
			*pc,
		)
	}

	if len(pc.Thresholds) > 0 {
		return thresholdBuckets(name, pc.Thresholds)
	}
	if pc.Range != nil {
		return pc.compileRange(name)
	}
	if pc.Exponential != nil {
		return pc.complileExponential(name)
	}
	if len(pc.Semver) > 0 {
		return pc.compileSemver(name)
	}
	if len(pc.Match) > 0 {
		return pc.compileMatch(name)
	}
	if pc.Bool {
		return boolBuckets{}, nil
	}
	if pc.Float {
		return rawFloat{}, nil
	}
	return nil, fmt.Errorf("%w %s: must specify bucket parameters",
		BucketConfigurationError, name)
}

func (pc *BucketSpec) compileRange(name string) (Bucketizer, error) {
	if pc.Range.Buckets == 0 {
		pc.Range.Buckets = 10
	} else if pc.Range.Buckets < 1 {
		return nil, fmt.Errorf("%w %s: invalid bucket count %d",
			BucketConfigurationError, name, pc.Range.Buckets)
	} else if pc.Range.Buckets > 20 {
		return nil, fmt.Errorf("%w %s: too many buckets %d",
			BucketConfigurationError, name, pc.Range.Buckets)
	}
	bucketSize := (pc.Range.Max - pc.Range.Min) / float64(pc.Range.Buckets)
	if !(bucketSize > 0) { // catch NaN as well here
		return nil, fmt.Errorf("%w %s: max <= min",
			BucketConfigurationError, name)
	}
	thresholds := make([]float64, pc.Range.Buckets+1)
	thresholds[0] = pc.Range.Min
	for i, v := range thresholds[:len(thresholds)-1] {
		thresholds[i+1] = v + bucketSize
	}
	// Set explicitly because of the potential for round-off error
	thresholds[len(thresholds)-1] = pc.Range.Max
	return thresholdBuckets(name, thresholds)
}

func (pc *BucketSpec) complileExponential(name string) (Bucketizer, error) {
	if !(pc.Exponential.Min > math.SmallestNonzeroFloat32) {
		return nil, fmt.Errorf(
			"%w %s: invalid exponential min %g",
			BucketConfigurationError, name, pc.Exponential.Min)
	}
	if !(pc.Exponential.Multiple >= 2) {
		return nil, fmt.Errorf(
			"%w %s: invalid exponential multiple %g",
			BucketConfigurationError, name, pc.Exponential.Multiple)
	}
	if pc.Exponential.Max == 0 {
		return &unboundedExponentialBuckets{
			Min:      pc.Exponential.Min,
			Multiple: pc.Exponential.Multiple,
		}, nil
	}
	if !(pc.Exponential.Max > pc.Exponential.Min*pc.Exponential.Multiple) {
		return nil, fmt.Errorf(
			"%w %s: max %g too small compared to min %g",
			BucketConfigurationError, name,
			pc.Exponential.Max, pc.Exponential.Min)
	}
	// estimate required bucket count
	bucketCount := int(math.Ceil(math.Log(pc.Exponential.Max/pc.Exponential.Min) /
		math.Log(pc.Exponential.Multiple)))
	thresholds := make([]float64, 1, bucketCount+1)
	thresholds[0] = pc.Exponential.Min
	for thresholds[len(thresholds)-1] < pc.Exponential.Max {
		thresholds = append(thresholds,
			thresholds[len(thresholds)-1]*pc.Exponential.Multiple)
	}
	// Round-off error, or in case the max threshold is not aligned to
	// a multiple of min.
	thresholds[len(thresholds)-1] = pc.Exponential.Max
	if len(thresholds) > 2 &&
		(pc.Exponential.Max-thresholds[len(thresholds)-2] < pc.Exponential.Min) {
		thresholds = thresholds[:len(thresholds)-1]
		thresholds[len(thresholds)-1] = pc.Exponential.Max
	}
	return thresholdBuckets(name, thresholds)
}

// thresholdBuckets includes common checks and implementation for floatBuckets.
// When a range or bounded exponential is encountered, it is compiled to a
// standard float bucket to cannonicalize the representation.
func thresholdBuckets(name string, thresholds []float64) (Bucketizer, error) {
	if len(thresholds) < 2 {
		return nil, fmt.Errorf("%w %s: at least two buckets required",
			BucketConfigurationError, name)
	}
	sort.Float64s(thresholds)
	for i, v := range thresholds {
		if math.IsNaN(v) {
			return nil, fmt.Errorf("%w %s: threshold value is not a number",
				BucketConfigurationError, name)
		}
		if math.IsInf(v, 0) {
			return nil, fmt.Errorf("%w %s: threshold value is not finite",
				BucketConfigurationError, name)
		}
		if i > 0 {
			diff := v - thresholds[i-1]
			if diff < math.SmallestNonzeroFloat32 {
				return nil, fmt.Errorf(
					"%w %s: threshold values %g and %g are too close together",
					BucketConfigurationError, name,
					thresholds[i-1], v)
			}
		}
	}
	return &floatBuckets{thresholds: thresholds}, nil
}

type floatBuckets struct {
	thresholds []float64
}

func (f *floatBuckets) bucket(v float64) any {
	for i, t := range f.thresholds {
		if t > v {
			if i == 0 {
				return fmt.Sprintf("<%g", t)
			}
			return fmt.Sprintf("%g-%g", f.thresholds[i-1], t)
		}
	}
	if f.thresholds[len(f.thresholds)-1] == v {
		return fmt.Sprintf("%g",
			f.thresholds[len(f.thresholds)-1])
	}

	return fmt.Sprintf(">%g",
		f.thresholds[len(f.thresholds)-1])
}

func coerceFloat(value any) any {
	if value == nil {
		return nil
	}
	switch v := value.(type) {
	case float64:
		if math.IsNaN(v) {
			return "NaN"
		} else if math.IsInf(v, -1) {
			return "-Inf"
		} else if math.IsInf(v, 1) {
			return "Inf"
		}
		return v
	case float32:
		return coerceFloat(float64(v))
	case int:
		return coerceFloat(float64(v))
	case uint:
		return coerceFloat(float64(v))
	case int8:
		return coerceFloat(float64(v))
	case uint8:
		return coerceFloat(float64(v))
	case int16:
		return coerceFloat(float64(v))
	case uint16:
		return coerceFloat(float64(v))
	case int32:
		return coerceFloat(float64(v))
	case uint32:
		return coerceFloat(float64(v))
	case int64:
		return coerceFloat(float64(v))
	case uint64:
		return coerceFloat(float64(v))
	case string:
		v = strings.TrimSpace(v)
		if n, err := strconv.ParseFloat(v, 64); err != nil {
			return "NaN"
		} else {
			return coerceFloat(n)
		}
	case json.Number:
		if n, err := v.Float64(); err != nil {
			return "NaN"
		} else {
			return coerceFloat(n)
		}
	case json.RawMessage:
		v = bytes.TrimSpace(v)
		if len(v) == 0 {
			return nil
		}
		if v[0] == '"' {
			var s string
			if json.Unmarshal(v, &s) == nil {
				return coerceFloat(s)
			}
		}
		if string(v) == "true" || string(v) == "false" {
			return "NaN"
		}
		if string(v) == "null" {
			return nil
		}
		return coerceFloat(json.Number(v))
	case []byte:
		return coerceFloat(string(bytes.TrimSpace(v)))
	}
	return "NaN"
}

func (f *floatBuckets) Bucket(value any) any {
	value = coerceFloat(value)
	if value == nil {
		return nil
	}
	if v, ok := value.(float64); ok {
		return f.bucket(v)
	}
	return value
}

func (f *floatBuckets) MarshalJSON() ([]byte, error) {
	return json.Marshal(struct {
		Thresholds []float64 `json:"thresholds"`
	}{f.thresholds})
}

func (floatBuckets) private() {}

type unboundedExponentialBuckets struct {
	Min      float64 `json:"min"`
	Multiple float64 `json:"multiple"`
}

func (u *unboundedExponentialBuckets) bucket(v float64) any {
	if v < u.Min {
		return fmt.Sprintf("<%g", u.Min)
	}
	bound := u.Min
	for bound*u.Multiple <= v {
		bound *= u.Multiple
	}
	return fmt.Sprintf(">=%g", bound)
}

// Bucket implements Bucketizer.
func (u *unboundedExponentialBuckets) Bucket(value any) any {
	value = coerceFloat(value)
	if v, ok := value.(float64); ok {
		return u.bucket(v)
	}
	return value
}

func (u *unboundedExponentialBuckets) MarshalJSON() ([]byte, error) {
	return json.Marshal(struct {
		Exp unboundedExponentialBuckets `json:"exponential"`
	}{*u})
}

func (unboundedExponentialBuckets) private() {}

func (pc *BucketSpec) compileSemver(name string) (Bucketizer, error) {
	th := make([][2]int, len(pc.Semver))
	for i, ver := range pc.Semver {
		maj, min, ok := strings.Cut(strings.TrimSpace(ver), ".")
		if !ok {
			return nil, fmt.Errorf("%w %s: version %q did not contain a .",
				BucketConfigurationError, name, ver)
		}
		mav, err := strconv.Atoi(maj)
		if err != nil {
			return nil, fmt.Errorf("%w %s: invalid major version %q: %w",
				BucketConfigurationError, name, maj, err)
		}
		miv, err := strconv.Atoi(min)
		if err != nil {
			return nil, fmt.Errorf("%w %s: invalid minor version %q: %w",
				BucketConfigurationError, name, min, err)
		}
		th[i] = [2]int{mav, miv}
	}
	slices.SortFunc(th, func(a, b [2]int) int {
		if x := cmp.Compare(a[0], b[0]); x != 0 {
			return -x
		}
		return -cmp.Compare(a[1], b[1])
	})
	for i, v := range pc.Semver[1:] {
		if v == pc.Semver[i] {
			return nil, fmt.Errorf("%w %s: duplicate version %d.%d",
				BucketConfigurationError, name, v[0], v[1])
		}
	}
	return &semverBuckets{
		thresholds: th,
	}, nil
}

type semverBuckets struct {
	thresholds [][2]int
}

var semverRe = regexp.MustCompile(`^(\d+)\.(\d+)(?:$|[a-z._-])`)

func coerceString(value any) (string, bool) {
	switch v := value.(type) {
	case string:
		return v, true
	case json.RawMessage:
		var s string
		if json.Unmarshal(v, &s) == nil {
			return s, true
		}
	case []byte:
		return string(v), true
	}
	return "", false
}

func (s *semverBuckets) Bucket(value any) any {
	if value == nil {
		return nil
	}
	v, ok := coerceString(value)
	if !ok {
		return nil
	}
	m := semverRe.FindStringSubmatch(strings.TrimSpace(v))
	if len(m) < 3 {
		return "invalid"
	}
	maj, err := strconv.Atoi(m[1])
	if err != nil {
		return "invalid"
	}
	min, err := strconv.Atoi(m[2])
	if err != nil {
		return "invalid"
	}

	for _, t := range s.thresholds {
		if t[0] < maj || t[0] == maj && t[1] <= min {
			return fmt.Sprintf(">=%d.%d", t[0], t[1])
		}
	}
	return fmt.Sprintf("<%d.%d",
		s.thresholds[len(s.thresholds)-1][0],
		s.thresholds[len(s.thresholds)-1][1])
}

func (s *semverBuckets) MarshalJSON() ([]byte, error) {
	th := make([]string, len(s.thresholds))
	for i, t := range s.thresholds {
		th[i] = fmt.Sprintf("%d.%d", t[0], t[1])
	}
	return json.Marshal(struct {
		Sv []string `json:"semver_thresholds"`
	}{th})
}

func (semverBuckets) private() {}

func (pc *BucketSpec) compileMatch(name string) (Bucketizer, error) {
	var errs []error
	buckets := make([]compiledReBucket, len(pc.Match))
	for i, exp := range pc.Match {
		if re, err := regexp.Compile(exp.Exp); err != nil {
			errs = append(errs, fmt.Errorf("%q: %w", exp.Exp, err))
		} else {
			if exp.Key == "" {
				buckets[i] = compiledReBucket{
					RegexpBucket: RegexpBucket{
						Key: exp.Exp,
						Exp: exp.Exp,
					},
					re: re,
				}
			} else {
				buckets[i] = compiledReBucket{
					RegexpBucket: exp,
					re:           re,
				}
			}
		}
	}
	if len(errs) > 0 {
		return nil, fmt.Errorf("%w %s: %w", BucketConfigurationError,
			name, errors.Join(errs...))
	}
	return &regexBuckets{
		buckets: buckets,
	}, nil
}

type compiledReBucket struct {
	RegexpBucket
	re *regexp.Regexp
}

type regexBuckets struct {
	buckets []compiledReBucket
}

// Bucket implements Bucketizer.
func (r *regexBuckets) Bucket(value any) any {
	if value == nil {
		return nil
	}
	v, ok := coerceString(value)
	if !ok {
		return nil
	}
	for _, b := range r.buckets {
		if b.re.MatchString(v) {
			return b.Key
		}
	}
	return "unknown"
}

// MarshalJSON implements Bucketizer.
func (r *regexBuckets) MarshalJSON() ([]byte, error) {
	return json.Marshal(struct {
		Match []compiledReBucket `json:"match"`
	}{r.buckets})
}

func (regexBuckets) private() {}

type boolBuckets struct{}

func (boolBuckets) Bucket(value any) any {
	if value == nil {
		return nil
	}
	switch v := value.(type) {
	case bool:
		return v
	case string:
		if v == "" || strings.EqualFold(v, "false") {
			return false
		}
		return true
	case float64:
		if v == 0 {
			return false
		}
		return true
	case float32:
		if v == 0 {
			return false
		}
		return true
	case int:
		if v == 0 {
			return false
		}
		return true
	case uint:
		if v == 0 {
			return false
		}
		return true
	case int8:
		if v == 0 {
			return false
		}
		return true
	case uint8:
		if v == 0 {
			return false
		}
		return true
	case int16:
		if v == 0 {
			return false
		}
		return true
	case uint16:
		if v == 0 {
			return false
		}
		return true
	case int32:
		if v == 0 {
			return false
		}
		return true
	case uint32:
		if v == 0 {
			return false
		}
		return true
	case int64:
		if v == 0 {
			return false
		}
		return true
	case uint64:
		if v == 0 {
			return false
		}
		return true
	case json.RawMessage:
		v = bytes.TrimSpace(v)
		if string(v) == "null" {
			return nil
		}
		var uv any
		if json.Unmarshal(v, &uv) == nil {
			return boolBuckets{}.Bucket(uv)
		}
	case []byte:
		return boolBuckets{}.Bucket(string(v))
	}
	return true
}

func (boolBuckets) MarshalJSON() ([]byte, error) {
	return []byte(`{"boolean": true}`), nil
}

func (boolBuckets) private() {}

type rawFloat struct {
}

func (f rawFloat) Bucket(value any) any {
	value = coerceFloat(value)
	if value == nil {
		return nil
	}
	if v, ok := value.(float64); ok {
		return v
	}
	return value
}

func (f rawFloat) MarshalJSON() ([]byte, error) {
	return []byte(`{"float": true}`), nil
}

func (rawFloat) private() {}
