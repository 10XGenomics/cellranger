package config

import (
	"time"
)

type Config struct {
	ConfigVersion int                           `json:"config_version"`
	ConfigTime    time.Time                     `json:"config_version_time"`
	Product       string                        `json:"product"`
	Version       string                        `json:"version"`
	Groups        map[string]*MetricGroupConfig `json:"groups"`
	Metrics       []*MetricConfig               `json:"metrics"`
	Presets       map[string]*BucketSpec        `json:"presets"`
}

type MetricGroupConfig struct {
	When        string          `json:"when"`
	Subcommands StringOrList    `json:"subcommands,omitempty"`
	Globs       SortedStringMap `json:"globs,omitempty"`
}

type MetricConfig struct {
	Name   string       `json:"name"`
	Groups StringOrList `json:"groups"`

	Values *BucketSpec `json:"values,omitempty"`
	Preset string      `json:"preset,omitempty"`

	// Data source
	Special     string             `json:"special,omitempty"`
	Filesystem  *FilesystemValue   `json:"filesystem,omitempty"`
	FileContent *FileContentConfig `json:"file_content,omitempty"`
	FileSize    *FileSizeConfig    `json:"file_size,omitempty"`
	FileCount   *FileSizeConfig    `json:"file_count,omitempty"`
	Rlimit      *RlimitValue       `json:"rlimit,omitempty"`
	Flag        string             `json:"flag,omitempty"`
	MroFlag     string             `json:"mro_flag,omitempty"`
	RecentCount string             `json:"recent_count,omitempty"`
}

type FilesystemValue struct {
	Which    string `json:"which"`
	Property string `json:"property"`
}

type FileContentConfig struct {
	FileSizeConfig
	Cgroup    string       `json:"cgroup,omitempty"`
	Converter string       `json:"converter,omitempty"`
	Line      string       `json:"line,omitempty"`
	LastLine  string       `json:"last_line,omitempty"`
	JsonPath  StringOrList `json:"json_path,omitempty"`
	StageRel  bool         `json:"stage_relative,omitempty"`
}

type FileSizeConfig struct {
	Path StringOrList `json:"path,omitempty"`
	Flag string       `json:"flag,omitempty"`
}

type CgroupSource struct {
	Controller string       `json:"controller"`
	File       StringOrList `json:"file"`
}

type RlimitValue struct {
	Hard     bool   `json:"hard"`
	Resource string `json:"resource"`
}

type RegexpBucket struct {
	Exp string `json:"exp"`
	Key string `json:"key,omitempty"`
}

type RangeBucketSpec struct {
	Min     float64 `json:"min"`
	Max     float64 `json:"max"`
	Buckets int     `json:"buckets"`
}

type ExpBucketSpec struct {
	Multiple float64 `json:"multiple"`
	Min      float64 `json:"min"`
	Max      float64 `json:"max,omitempty"`
}

type BucketSpec struct {
	Thresholds  []float64        `json:"thresholds,omitempty"`
	Range       *RangeBucketSpec `json:"range,omitempty"`
	Exponential *ExpBucketSpec   `json:"exponential,omitempty"`
	Semver      []string         `json:"semver_thresholds,omitempty"`
	Bool        bool             `json:"boolean,omitempty"`
	Float       bool             `json:"float,omitempty"`
	Match       []RegexpBucket   `json:"match,omitempty"`
}
