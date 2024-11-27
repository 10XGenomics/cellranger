package collector

import (
	"context"
	"io/fs"

	"golang.org/x/sys/unix"
)

type fileResult struct {
	info     fs.FileInfo
	fullPath string
	relPath  string
}

type errorFile struct {
	path    string
	content string
}

type cgroupRequest struct {
	result chan<- []byte
	names  []string
}

// In order to facilitate having multiple metrics collected
// from the same file (e.g. different lines or json keys),
// a shared instance of this gets passed in to all of the
// value extractor's Start functions.
//
// Essentially this works a bit like a sync.OnceValue,
// except with several keyed maps, and without the concurrency
// protections.
//
// Once all active metrics' Start functions have been called,
// the FileContext's `Close()` method should be called.
// It is undefined behavior for any other methods to be called
// after that point.
// The `Close()` method will launch goroutines to start
// collecting the requested data and distribute it to the waiting
// value extractors.
type FileContext struct {
	statfs  map[string][]chan<- unix.Statfs_t
	mounts  map[string][]chan<- mountInfo
	cgroups map[string][]cgroupRequest

	files multiFileResultMap
	// We don't want to walk the directory structure more than once.
	stageRelFiles multiFileResultMap
	errorFiles    []chan<- errorFile
}

func (fc *FileContext) errorFile() <-chan errorFile {
	c := make(chan errorFile, 1)
	fc.errorFiles = append(fc.errorFiles, c)
	return c
}

func (fc *FileContext) stageRelative(p []string) *fileSink {
	if fc.stageRelFiles == nil {
		fc.stageRelFiles = make(multiFileResultMap)
	}
	return fc.stageRelFiles.getSink(p)
}

func (fc *FileContext) statFs(p string) <-chan unix.Statfs_t {
	if fc.statfs == nil {
		fc.statfs = make(map[string][]chan<- unix.Statfs_t, 2)
	}
	c := make(chan unix.Statfs_t, 1)
	fc.statfs[p] = append(fc.statfs[p], c)
	return c
}

func (fc *FileContext) mountInfo(p string) <-chan mountInfo {
	if fc.statfs == nil {
		fc.mounts = make(map[string][]chan<- mountInfo, 2)
	}
	c := make(chan mountInfo, 1)
	fc.mounts[p] = append(fc.mounts[p], c)
	return c
}

func (fc *FileContext) cgroup(ctrl string, paths []string) <-chan []byte {
	if fc.cgroups == nil {
		fc.cgroups = make(map[string][]cgroupRequest)
	}
	c := make(chan []byte, 1)
	fc.cgroups[ctrl] = append(fc.cgroups[ctrl], cgroupRequest{
		result: c,
		names:  paths,
	})
	return c
}

func (fc *FileContext) Close(ctx context.Context, pipestance string) {
	for p, cs := range fc.statfs {
		go doStatfs(ctx, p, cs)
	}
	fc.statfs = nil // Allow gc to clean up
	var cgm chan map[string]string
	if len(fc.cgroups) > 0 {
		cgm = make(chan map[string]string, 1)
	}
	if len(fc.mounts) > 0 || len(fc.cgroups) > 0 {
		go doMountInfos(ctx, fc.mounts, cgm)
	}
	fc.mounts = nil
	if len(fc.cgroups) > 0 {
		go doCgroups(ctx, fc.cgroups, cgm)
	}
	fc.cgroups = nil

	if len(fc.errorFiles) > 0 || len(fc.stageRelFiles) > 0 {
		go findRelativeFiles(ctx, pipestance,
			fc.errorFiles, fc.stageRelFiles)
	}
	fc.errorFiles = nil
	fc.stageRelFiles = nil
	if len(fc.files) > 0 {
		go findFiles(ctx, pipestance, fc.files)
	}
	fc.files = nil
}
