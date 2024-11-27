package config

import (
	"context"
	"crypto/sha256"
	"encoding/json"
	"errors"
	"fmt"
	"path/filepath"
	"strings"
	"sync"
	"time"

	"github.com/10XDev/cellranger/lib/go/telemetry/collector"
)

type CompiledConfig struct {
	ConfigVersion int       `json:"config_version"`
	ConfigTime    time.Time `json:"config_time"`
	Product       string    `json:"product"`
	Version       string    `json:"version"`
	Groups        map[string]*MetricGroup
}

// A compiled, validated metric group.
type MetricGroup struct {
	MetricGroupConfig
	Metrics SortedMetricMap
}

func (g *MetricGroup) Hash() [sha256.Size]byte {
	hashContent, err := json.Marshal(g)
	if err != nil {
		panic(err)
	}
	return sha256.Sum256(hashContent)
}

// A compiled metric.
type Metric struct {
	Name   string
	Source ValueExtractor
	Bucket Bucketizer
	get    func() any
}

// A source of metric data.
type ValueExtractor = collector.ValueExtractor

// Converts a raw metric to a bucketed value.
type Bucketizer interface {
	json.Marshaler

	// Converts a raw metric source to a bucketed value.
	Bucket(any) any

	// Cannot be implemented outside of this package.
	private()
}

func (m *Metric) WithGlob(globs map[string]string) *Metric {
	if ve := collector.WithGlob(m.Source, globs); ve == m.Source {
		return m
	} else {
		return &Metric{
			Name:   m.Name,
			Source: ve,
			Bucket: m.Bucket,
		}
	}
}

func (m *Metric) Start(ctx context.Context, vc collector.ValueContext, fc *collector.FileContext) {
	result := make(chan any, 1)
	m.Source.Start(ctx, vc, fc, result)
	m.get = sync.OnceValue(func() any {
		select {
		case r := <-result:
			if r == nil {
				return nil
			}
			if m.Bucket != nil {
				return m.Bucket.Bucket(r)
			}
			return r
		case <-ctx.Done():
			return nil
		}
	})
}

func (m *Metric) Get() any {
	if m.get == nil {
		panic("Get called without a call to Start")
	}
	return m.get()
}

var (
	GroupConfigurationError  = errors.New("invalid group")
	MetricConfigurationError = errors.New("invalid metric")
	BucketConfigurationError = errors.New("invalid bucketization")
)

const (
	WhenCommand = "command"
	WhenSuccess = "success"
	WhenFailure = "failure"
)

// Compile validates a config and converts it to a canonical
// form to be used in actually collecting metrics.
func (c *Config) Compile() (CompiledConfig, error) {
	var errs []error
	if c.ConfigVersion < 0 {
		errs = append(errs, fmt.Errorf("invalid version %d",
			c.ConfigVersion))
	}
	if c.Product == "" {
		errs = append(errs, errors.New("product not specified"))
	}
	if c.Version == "" {
		errs = append(errs, errors.New("version not specified"))
	}
	if len(c.Groups) == 0 {
		errs = append(errs, errors.New("no groups specified"))
	}
	if len(c.Metrics) == 0 {
		errs = append(errs, errors.New("no metrics specified"))
	}
	cc := CompiledConfig{
		ConfigVersion: c.ConfigVersion,
		ConfigTime:    c.ConfigTime,
		Product:       c.Product,
		Version:       c.Version,
		Groups:        make(map[string]*MetricGroup, len(c.Groups)),
	}
	for n, mg := range c.Groups {
		g, err := mg.compile(n)
		if err != nil {
			errs = append(errs, err)
			continue
		}
		if n == "" {
			errs = append(errs, fmt.Errorf("%w: blank group name",
				GroupConfigurationError))
		}
		cc.Groups[n] = g
	}
	var presets map[string]Bucketizer
	if len(c.Presets) > 0 {
		presets = make(map[string]Bucketizer, len(c.Presets))
		for n, pc := range c.Presets {
			if p, err := pc.compile(n); err != nil {
				errs = append(errs, err)
			} else {
				presets[n] = p
			}
		}
	}
	for _, mc := range c.Metrics {
		m, err := mc.compile(presets)
		if err != nil {
			errs = append(errs, err)
			continue
		}
		for _, gn := range mc.Groups {
			g := cc.Groups[gn]
			if g == nil {
				errs = append(errs, fmt.Errorf(
					"%w %s: unknown group %s",
					MetricConfigurationError, mc.Name, gn))
				continue
			}
			if err := g.addMetric(gn, m); err != nil {
				errs = append(errs, err)
			}
		}
	}
	for name, mg := range cc.Groups {
		if len(mg.Metrics) == 0 {
			errs = append(errs, fmt.Errorf("%w: no metrics defined for group %s",
				GroupConfigurationError, name))
		}
	}
	return cc, errors.Join(errs...)
}

func (mg *MetricGroupConfig) compile(name string) (*MetricGroup, error) {
	if mg == nil {
		return nil, fmt.Errorf("%w %s: null group config",
			GroupConfigurationError, name)
	}
	var errs []error
	g := MetricGroup{
		MetricGroupConfig: *mg,
		Metrics:           make(SortedMetricMap),
	}
	switch mg.When {
	case WhenCommand:
		if len(mg.Globs) != 0 {
			errs = append(errs,
				fmt.Errorf(
					"%w %s: globs are not permitted for groups which run on command",
					GroupConfigurationError, name))
		}
		g.Globs = nil
	case WhenSuccess:
		if len(mg.Subcommands) == 0 {
			errs = append(errs, fmt.Errorf(
				"%w %s: success groups must specify a subcommand",
				GroupConfigurationError, name))
		}
	case WhenFailure:
	default:
		errs = append(errs, fmt.Errorf(
			"%w %s: invalid value %q for when",
			GroupConfigurationError, name, mg.When))
		return nil, errors.Join(errs...)
	}
	for n, glob := range g.Globs {
		if n == "" {
			errs = append(errs, fmt.Errorf(
				"%w %s: empty glob key",
				GroupConfigurationError, name))
		} else if glob == "" {
			errs = append(errs, fmt.Errorf(
				"%w %s: empty glob for key %s",
				GroupConfigurationError, name, n,
			))
		} else if _, err := filepath.Match(glob, "test"); err != nil {
			errs = append(errs, fmt.Errorf(
				"%w %s: glob %s: %w",
				GroupConfigurationError, name, n, err,
			))
		} else if strings.Contains(glob, "..") || strings.HasPrefix(glob, "/") {
			errs = append(errs, fmt.Errorf(
				"%w %s: glob %s: not relative to pipestance directory",
				GroupConfigurationError, name, n,
			))
		}
	}
	return &g, errors.Join(errs...)
}

func (mg *MetricGroup) addMetric(gn string, m *Metric) error {
	var errs []error
	if !mg.Metrics.Add(m) {
		return fmt.Errorf(
			"%w: multiple metrics named %q in group %s",
			MetricConfigurationError,
			m.Name, gn)
	}
	if m.Source.RequiresPipeline() {
		if m.Source.RequiresSuccess() && mg.When != WhenSuccess {
			errs = append(errs, fmt.Errorf(
				"%w: metric %s cannot be added to group %s because"+
					" its source is only valid on success, but %s "+
					"triggers on %s",
				MetricConfigurationError,
				m.Name, gn, gn, mg.When))
		} else if m.Source.RequiresFailure() && mg.When != WhenFailure {
			errs = append(errs, fmt.Errorf(
				"%w: metric %s cannot be added to group %s because"+
					" its source is only valid on failure, but %s "+
					"triggers on %s",
				MetricConfigurationError,
				m.Name, gn, gn, mg.When))
		} else if mg.When != WhenSuccess && mg.When != WhenFailure {
			errs = append(errs, fmt.Errorf(
				"%w: metric %s cannot be added to group %s because"+
					" its source is only valid after a pipeline completes, "+
					"but %s triggers on %s",
				MetricConfigurationError,
				m.Name, gn, gn, mg.When))
		}
		for _, glob := range m.Source.RequiredGlobKeys() {
			if v := mg.Globs[glob]; v == "" {
				errs = append(errs, fmt.Errorf(
					"%w: metric %s cannot be added to group %s because"+
						" its source expects a glob key %s",
					MetricConfigurationError,
					m.Name, gn, glob))
			}
		}
	}
	return errors.Join(errs...)
}

func (mc *MetricConfig) checkOneSource() error {
	sources := make([]string, 0, 8)
	if mc.Special != "" {
		sources = append(sources, "special")
	}
	if mc.Filesystem != nil {
		sources = append(sources, "filesystem")
	}
	if mc.FileContent != nil {
		sources = append(sources, "file_content")
	}
	if mc.FileSize != nil {
		sources = append(sources, "file_size")
	}
	if mc.FileCount != nil {
		sources = append(sources, "file_count")
	}
	if mc.Rlimit != nil {
		sources = append(sources, "rlimit")
	}
	if mc.Flag != "" {
		sources = append(sources, "flag")
	}
	if mc.MroFlag != "" {
		sources = append(sources, "mro_flag")
	}
	if mc.RecentCount != "" {
		sources = append(sources, "recent_count")
	}
	if len(sources) == 0 {
		return errors.New("no sources specified")
	} else if len(sources) != 1 {
		return fmt.Errorf("%d sources specified: %s",
			len(sources), strings.Join(sources, ","))
	}
	return nil
}

func (mc *MetricConfig) compileSource() (ValueExtractor, error) {
	if err := mc.checkOneSource(); err != nil {
		return nil, err
	}
	switch {
	case mc.Special != "":
		return getSpecialSource(mc.Special)
	case mc.Filesystem != nil:
		return mc.Filesystem.compile()
	case mc.FileContent != nil:
		return mc.FileContent.compile()
	case mc.FileSize != nil:
		return mc.FileSize.compile(false)
	case mc.FileCount != nil:
		return mc.FileCount.compile(true)
	case mc.Rlimit != nil:
		return mc.Rlimit.compile()
	case mc.Flag != "":
		return compileFlagSource(mc.Flag, false)
	case mc.MroFlag != "":
		return compileFlagSource(mc.MroFlag, true)
	case mc.RecentCount != "":
		return compileRecentCountSource(mc.RecentCount, mc.Groups)
	}
	return nil, errors.New("no sources specified")
}

func (mc *MetricConfig) compile(presets map[string]Bucketizer) (*Metric, error) {
	if mc.Name == "" {
		return nil, fmt.Errorf("%w: empty name",
			MetricConfigurationError)
	}
	v, err := mc.compileSource()
	if err != nil {
		return nil, fmt.Errorf("%w %s: %w",
			MetricConfigurationError, mc.Name, err)
	}
	m := Metric{
		Name:   mc.Name,
		Source: v,
	}
	if mc.Preset != "" {
		if mc.Values != nil {
			return nil, fmt.Errorf(
				"%w %s: both preset and values were specified",
				MetricConfigurationError, mc.Name)
		}
		if b := presets[mc.Preset]; b == nil {
			return nil, fmt.Errorf("%w %s: unknown preset %s",
				MetricConfigurationError, mc.Name, mc.Preset)
		} else {
			m.Bucket = b
		}
	} else if mc.Values == nil {
		if !v.BucketsOptional() {
			return nil, fmt.Errorf("%w %s: either values or preset must be specified",
				MetricConfigurationError, mc.Name)
		}
	} else if b, err := mc.Values.compile(mc.Name); err != nil {
		return nil, fmt.Errorf("%w %s: %w",
			MetricConfigurationError, mc.Name, err)
	} else {
		m.Bucket = b
	}
	return &m, nil
}
