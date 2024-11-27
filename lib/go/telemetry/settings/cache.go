package settings

import (
	"errors"
	"fmt"
	"io/fs"
	"os"
	"os/user"
	"path"
	"path/filepath"
	"time"
)

// CacheDir returns the path to the directory containing previously
// collected data.
func CacheDir() (string, error) {
	c, err := os.UserCacheDir()
	if err != nil {
		// os.UserCacheDir can fail if $HOME isn't set, but there's another
		// way to get the user's home directory if you're willing to take a
		// dependency on libc, which we are, here.
		if u, err := user.Current(); err == nil && u.HomeDir != "" {
			return path.Join(u.HomeDir, ".cache", "tenx", "telemetry"), nil
		}
		return "", err
	}
	return path.Join(c, "tenx", "telemetry"), nil
}

func ProductCacheDir() (string, error) {
	c, err := CacheDir()
	if err != nil {
		return "", fmt.Errorf(
			"could not determine location of cache directory: %w", err)
	}
	p, err := getProduct()
	if err != nil {
		return "", err
	}
	v, err := getVersion()
	if err != nil {
		return "", err
	}
	return filepath.Join(c, p, v), nil
}

func ListCache() ([]string, error) {
	d, err := ProductCacheDir()
	if err != nil {
		return nil, err
	}
	listing, err := os.ReadDir(d)
	if err != nil {
		if errors.Is(err, fs.ErrNotExist) {
			return nil, nil
		}
		return nil, fmt.Errorf("could not list cache directory: %w", err)
	}
	r := make([]string, 0, len(listing))
	for _, ent := range listing {
		if !ent.IsDir() {
			r = append(r, filepath.Join(d, ent.Name()))
		}
	}
	return r, nil
}

// CountRecentCacheFiles returns the number of cache files that have been saved
// for each group in the past 30 days.
func CountRecentCacheFiles() map[string]int {
	listing, err := ListCache()
	if err != nil {
		return nil
	}
	if len(listing) == 0 {
		return nil
	}
	result := make(map[string]int)
	for _, fn := range listing {
		group, tm := parseEventFileName(fn)
		if group != "" && !tm.IsZero() && time.Since(tm) <= maxEventTTL {
			result[group] += 1
		}
	}
	return result
}
