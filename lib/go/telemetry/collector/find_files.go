package collector

import (
	"archive/zip"
	"bufio"
	"bytes"
	"context"
	"encoding/json"
	"errors"
	"fmt"
	"io"
	"io/fs"
	"os"
	"path"
	"path/filepath"
	"regexp"
	"runtime/trace"
	"sort"
	"strings"
	"sync"

	"github.com/10XDev/cellranger/lib/go/telemetry/converter"
)

var chunkDirNameRe = regexp.MustCompile(
	`^(?:(?:chnk[0-9]+)|split|join)(?:-u[0-9a-z]+)?$`,
)

func isChunkDir(name string) bool {
	return chunkDirNameRe.MatchString(name)
}

type jsonSink struct {
	// Listeners for sub-keys.
	children map[string]*jsonSink
	// Listeners for the value at this level of hierarchy
	listeners []chan<- any
}

func (j *jsonSink) subscribe(p []string, r chan<- any) {
	if len(p) == 0 {
		j.listeners = append(j.listeners, r)
		return
	}
	if j.children == nil {
		child := new(jsonSink)
		j.children = map[string]*jsonSink{
			p[0]: child,
		}
		child.subscribe(p[1:], r)
	} else if child := j.children[p[0]]; child == nil {
		child = new(jsonSink)
		j.children[p[0]] = child
		child.subscribe(p[1:], r)
	} else {
		child.subscribe(p[1:], r)
	}
}

func (j *jsonSink) close() {
	for _, c := range j.listeners {
		close(c)
	}
	for _, c := range j.children {
		c.close()
	}
}

func readJson(ctx context.Context, p string, sinks map[converter.Converter]*jsonSink) {
	defer trace.StartRegion(ctx, "readJson").End()
	content, err := os.ReadFile(p)
	if err != nil {
		for _, sink := range sinks {
			sink.close()
		}
		return
	}
	for conv, sink := range sinks {
		if ctx.Err() != nil {
			sink.close()
			continue
		}
		if conv != nil {
			converted, err := conv.Convert(content)
			if err == nil {
				sink.process(converted)
			}
		} else {
			sink.process(content)
		}
	}
}

func unmarshalPrimative(b json.RawMessage) (any, bool) {
	if len(b) == 0 {
		return nil, false
	}
	switch b[0] {
	case '[', '{':
		return nil, false
	case '"':
		var s string
		if json.Unmarshal(b, &s) == nil {
			return s, true
		}
	case 'n':
		if string(b) == "null" {
			return nil, true
		}
	case 't', 'f':
		var v bool
		if json.Unmarshal(b, &v) == nil {
			return v, true
		}
	default:
		var v json.Number
		if json.Unmarshal(b, &v) == nil {
			return v, true
		}
	}
	return nil, false
}

func (j *jsonSink) process(content json.RawMessage) {
	content = bytes.TrimSpace(content)
	if content[0] == '[' {
		// Take the first array element.
		// This is mostly to simplify handling of `_perf`, from which
		// we may wish to extract the storage high-water-mark.
		var arr []json.RawMessage
		if json.Unmarshal(content, &arr) == nil && len(arr) > 0 {
			j.process(arr[0])
		} else {
			j.close()
		}
		return
	}
	defer j.close()
	if len(j.listeners) > 0 {
		if v, ok := unmarshalPrimative(content); ok {
			for _, r := range j.listeners {
				r <- v
				close(r)
			}
		} else {
			for _, r := range j.listeners {
				close(r)
			}
		}
		j.listeners = nil
	}
	if len(j.children) > 0 {
		var m map[string]json.RawMessage
		if json.Unmarshal(content, &m) == nil {
			for k, s := range j.children {
				if v, ok := m[k]; ok {
					s.process(v)
				} else {
					s.close()
				}
			}
			j.children = nil
		}
	}
}

// A fileSink represents subscribers waiting for results from a file.
type fileSink struct {
	// Listeners waiting for specific line results.
	lineListeners map[converter.Converter]lineFilterSet
	// listeners for json keys
	jsonListeners  map[converter.Converter]*jsonSink
	key            string
	alternatePaths []string
}

func (sink *fileSink) subscribe(filter FileFilter, conv converter.Converter, result chan<- any) {
	switch f := filter.(type) {
	case *LineFilter:
		if sink.lineListeners == nil {
			sink.lineListeners = map[converter.Converter]lineFilterSet{
				conv: {
					*f: []chan<- any{result},
				},
			}
		} else {
			m := sink.lineListeners[conv]
			if m == nil {
				sink.lineListeners[conv] = lineFilterSet{
					*f: []chan<- any{result},
				}
			} else {
				m[*f] = append(m[*f], result)
			}
		}
	case JsonPathFilter:
		var root *jsonSink
		if sink.jsonListeners == nil {
			root = new(jsonSink)
			sink.jsonListeners = map[converter.Converter]*jsonSink{
				conv: root,
			}
		} else {
			root = sink.jsonListeners[conv]
			if root == nil {
				root = new(jsonSink)
				sink.jsonListeners[conv] = root
			}
		}
		root.subscribe(f, result)
	default:
		validateFilterType(filter)
	}
}

func (sink *fileSink) notify(ctx context.Context, res fileResult) {
	if sink.jsonListeners != nil {
		go readJson(ctx, res.fullPath, sink.jsonListeners)
		sink.jsonListeners = nil
	}
	if len(sink.lineListeners) > 0 {
		go readLines(ctx, res.fullPath, sink.lineListeners)
		sink.lineListeners = nil
	}
}

func readLines(ctx context.Context, p string,
	results map[converter.Converter]lineFilterSet) {
	for conv, filters := range results {
		if conv == nil {
			readRawLines(ctx, p, filters)
		} else {
			convertLines(ctx, p, conv, filters)
		}
	}
}

type lineFilterSet map[LineFilter][]chan<- any

func (lf lineFilterSet) close() {
	for _, cs := range lf {
		for _, c := range cs {
			close(c)
		}
	}
}

func readRawLines(ctx context.Context, p string, filters lineFilterSet) {
	defer trace.StartRegion(ctx, "readRawLines").End()
	f, err := os.Open(p)
	if err != nil || ctx.Err() != nil {
		filters.close()
		return
	}
	defer f.Close()

	first := make(map[*compiledLineFilter][]chan<- any)
	type lastSeen struct {
		saw string
		r   []chan<- any
	}
	last := make(map[*compiledLineFilter]*lastSeen)

	for lf, cs := range filters {
		cf := lf.compile()
		if cf.last {
			last[&cf] = &lastSeen{r: cs}
		} else {
			first[&cf] = cs
		}
	}

	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		if ctx.Err() != nil {
			for _, cs := range first {
				for _, c := range cs {
					close(c)
				}
			}
			for _, seen := range last {
				for _, c := range seen.r {
					close(c)
				}
			}
			return
		}
		for lf, cs := range first {
			if m, ok := lf.matchLine(scanner.Bytes()); ok {
				s := string(m)
				for _, c := range cs {
					c <- s
					close(c)
				}
				delete(first, lf)
			}
		}
		if len(first) == 0 && len(last) == 0 {
			return
		}
		for lf, seen := range last {
			if m, ok := lf.matchLine(scanner.Bytes()); ok {
				seen.saw = string(m)
			}
		}
	}
	for _, cs := range first {
		for _, c := range cs {
			close(c)
		}
	}
	for _, seen := range last {
		for _, c := range seen.r {
			if seen.saw != "" && scanner.Err() == nil {
				c <- seen.saw
			}
			close(c)
		}
	}
}

func convertLines(ctx context.Context, p string,
	conv converter.Converter, filters lineFilterSet) {
	defer trace.StartRegion(ctx, "convertLines").End()
	b, err := os.ReadFile(p)
	if err != nil || ctx.Err() != nil {
		filters.close()
		return
	}
	converted, err := conv.Convert(b)
	if err != nil || ctx.Err() != nil {
		filters.close()
		return
	}
	cStr := string(converted)
	for filter, cs := range filters {
		v, ok := filter.compile().matchContent(cStr)
		for _, c := range cs {
			if ok {
				c <- v
			}
			close(c)
		}
	}
}

func (sink *fileSink) close() {
	for _, ls := range sink.lineListeners {
		for _, cs := range ls {
			for _, c := range cs {
				close(c)
			}
		}
	}
	sink.lineListeners = nil
	for _, j := range sink.jsonListeners {
		j.close()
	}
	sink.jsonListeners = nil
}

// We can have multiple requests for a single file, and we can have multiple
// files for a single request.
// This type encapsulates the complexity of making sure we
// 1. Don't notify a given requestor more than once.
// 2. Don't keep searching for files we no longer need.
type multiFileResultMap map[string]map[string]*fileSink

func (m multiFileResultMap) close() {
	for _, sinks := range m {
		for _, sink := range sinks {
			sink.close()
		}
	}
}

func (m multiFileResultMap) getSink(paths []string) *fileSink {
	if len(paths) == 0 {
		panic("at least one path required")
	}
	var key string
	if len(paths) > 1 {
		sort.Strings(paths)
		key = strings.Join(paths, "\x00")
	}
	cs := m[paths[0]]
	if sink := cs[key]; sink != nil {
		return sink
	}
	sink := &fileSink{
		key: key,
	}
	if len(paths) > 1 {
		// Re-split rather than simply using `paths` so we don't
		// retain a reference to `paths`.
		sink.alternatePaths = strings.SplitN(key, "\x00", len(paths))
	}
	if cs == nil {
		m[paths[0]] = map[string]*fileSink{key: sink}
	} else {
		cs[key] = sink
	}
	for _, p := range paths[1:] {
		cs := m[p]
		if cs == nil {
			m[p] = map[string]*fileSink{key: sink}
		} else {
			cs[key] = sink
		}
	}
	return sink
}

func (m multiFileResultMap) notify(ctx context.Context, fn string,
	sinks map[string]*fileSink, result fileResult) {
	// assert: sinks == m[fn].  This is perhaps a premature optimization.
	for _, sink := range sinks {
		sink.notify(ctx, result)
		for _, alt := range sink.alternatePaths {
			if alt != fn {
				o := m[alt]
				delete(o, sink.key)
				if len(o) == 0 {
					delete(m, alt)
				}
			}
		}
	}
}

func findFiles(ctx context.Context, pipestance string, requests multiFileResultMap) {
	defer trace.StartRegion(ctx, "findFiles").End()
	defer requests.close()
	psfs := os.DirFS(pipestance)
	sfs := psfs.(fs.StatFS)
	for p, sinks := range requests {
		if !path.IsAbs(p) {
			paths, _ := fs.Glob(psfs, p)
			if len(paths) > 0 {
				fp := paths[0]
				info, _ := sfs.Stat(fp)
				requests.notify(ctx, p, sinks, fileResult{
					info:     info,
					fullPath: path.Join(pipestance, fp),
					relPath:  fp,
				})
			}
		} else if info, err := os.Stat(p); err == nil {
			requests.notify(ctx, p, sinks, fileResult{
				info:     info,
				fullPath: p,
				relPath:  p,
			})
		}
	}
}

func findRelativeFiles(ctx context.Context, pipestance string,
	errFiles []chan<- errorFile, requests multiFileResultMap) {
	defer trace.StartRegion(ctx, "findRelativeFiles").End()
	for _, c := range errFiles {
		cc := c // don't capture the loop variable
		defer close(cc)
	}
	defer requests.close()
	_ = filepath.WalkDir(pipestance, func(p string, d fs.DirEntry, err error) error {
		defer trace.StartRegion(ctx, "findRelativeFiles_walk").End()
		if ctx.Err() != nil {
			return fs.SkipAll
		}
		rp := strings.TrimPrefix(strings.TrimPrefix(p, pipestance), "/")
		if !d.IsDir() {
			return nil
		}
		// Don't walk top-level outs or extras directories.
		// They don't contain stage-relative files.
		if rp == "outs" || rp == "extras" {
			return fs.SkipDir
		}
		if isChunkDir(d.Name()) {
			if len(errFiles) > 0 {
				// Must do this synchronously because we need to take just the
				// first result.
				if content := readErrorFile(p); len(content) > 0 {
					efr := errorFile{
						path:    rp,
						content: content,
					}
					for _, c := range errFiles {
						c <- efr
					}
					errFiles = nil
				}
			}
			for fn, sinks := range requests {
				fp := path.Join(p, fn)
				if info, err := os.Stat(fp); err == nil {
					requests.notify(ctx, fn, sinks, fileResult{
						info:     info,
						fullPath: fp,
						relPath:  strings.TrimPrefix(strings.TrimPrefix(fp, pipestance), "/"),
					})
					break
				}
			}
			if len(errFiles) == 0 && len(requests) == 0 {
				return fs.SkipAll
			}
			return fs.SkipDir
		}
		return nil
	})
	if len(errFiles) == 0 {
		return
	}
	// We might have metadata files in a `_metadata.zip`
	// It's a bit tricky to support this in the general case
	// but also largely unnecessary except in the case of error
	// files.
	zr, err := zip.OpenReader(path.Join(pipestance, "_metadata.zip"))
	if err != nil {
		return
	}
	defer zr.Close()
	for _, zf := range zr.File {
		if zf.FileInfo().IsDir() {
			continue
		}
		bn := path.Base(zf.Name)
		if bn == "_assert" || bn == "_error" {
			if content := readErrorFromZip(zf); content != "" {
				efr := errorFile{
					path:    zf.Name,
					content: content,
				}
				for _, c := range errFiles {
					c <- efr
				}
				return
			}
		}
	}
}

func readErrorFile(p string) string {
	content, err := os.ReadFile(path.Join(p, "_assert"))
	if err == nil {
		return string(bytes.TrimSpace(content))
	}
	content, err = os.ReadFile(path.Join(p, "_errors"))
	if err != nil {
		return ""
	}
	content = bytes.TrimSpace(content)
	if string(content) == "Caught signal terminated" {
		// Usually, this is what you get when `mrp` shuts down
		// due to some other stage encountering an error.
		// It's not helpful to report it as "the" error.
		return ""
	}
	return string(content)
}

func readErrorFromZip(zf *zip.File) string {
	zfr, err := zf.Open()
	if err != nil {
		return ""
	}
	content, err := io.ReadAll(zfr)
	zfr.Close()
	if err != nil {
		return ""
	}
	content = bytes.TrimSpace(content)
	if string(content) == "Caught signal terminated" &&
		path.Base(zf.Name) == "_errors" {
		return ""
	}
	return string(content)
}

type fileContentExtractor struct {
	filter FileFilter
	conv   converter.Converter
	paths  []string
}

func (f *fileContentExtractor) MarshalJSON() ([]byte, error) {
	b, err := json.Marshal(struct {
		Filter FileFilter          `json:"filter"`
		Conv   converter.Converter `json:"converter,omitempty"`
		Paths  []string            `json:"paths"`
	}{
		Filter: f.filter,
		Conv:   f.conv,
		Paths:  f.paths,
	})
	if err != nil {
		return nil, err
	}
	return append(append([]byte(`{"file_content":`), b...), '}'), nil
}

// Start implements ValueExtractor.
func (f *fileContentExtractor) Start(ctx context.Context, _ ValueContext,
	fc *FileContext, result chan<- any) {
	if fc.files == nil {
		fc.files = make(multiFileResultMap)
	}
	fc.files.getSink(f.paths).subscribe(f.filter, f.conv, result)
}

func (f *fileContentExtractor) RequiresPipeline() bool {
	return !path.IsAbs(f.paths[0])
}

func (fileContentExtractor) RequiresFailure() bool {
	return false
}

func (f *fileContentExtractor) RequiresSuccess() bool {
	return f.RequiresPipeline()
}

var globKeyRegex = regexp.MustCompile(`\$\{([^}]+)\}`)

func appendGlobKeys(keys []string, p string) []string {
	for _, m := range globKeyRegex.FindAllStringSubmatch(p, -1) {
		keys = append(keys, m[1])
	}
	return keys
}

func (f *fileContentExtractor) RequiredGlobKeys() []string {
	if !f.RequiresPipeline() {
		return nil
	}
	var keys []string
	for _, p := range f.paths {
		keys = appendGlobKeys(keys, p)
	}
	return keys
}

func (f *fileContentExtractor) withGlob(globs map[string]string) ValueExtractor {
	keys := f.RequiredGlobKeys()
	if len(keys) == 0 {
		return f
	}
	replacements := make([]string, 0, 2*len(keys))
	for _, key := range keys {
		replace, ok := globs[key]
		if !ok {
			panic("this should have failed configuration validation")
		}
		replacements = append(replacements, "${"+key+"}", replace)
	}
	replacer := strings.NewReplacer(replacements...)
	newPaths := make([]string, len(f.paths))
	for i, p := range f.paths {
		newPaths[i] = replacer.Replace(p)
	}
	return &fileContentExtractor{
		filter: f.filter,
		conv:   f.conv,
		paths:  newPaths,
	}
}

func (fileContentExtractor) BucketsOptional() bool { return false }

func MakeAbsFile(ps []string, conv converter.Converter,
	filter FileFilter) ValueExtractor {
	return &fileContentExtractor{
		filter: filter,
		conv:   conv,
		paths:  ps,
	}
}

type stageRelativeExtractor struct {
	simpleExtractor
	filter FileFilter
	conv   converter.Converter
	paths  []string
}

func (f *stageRelativeExtractor) MarshalJSON() ([]byte, error) {
	b, err := json.Marshal(struct {
		Filter FileFilter          `json:"filter"`
		Conv   converter.Converter `json:"converter,omitempty"`
		Paths  []string            `json:"paths"`
		Rel    bool                `json:"stage_relative"`
	}{
		Filter: f.filter,
		Conv:   f.conv,
		Paths:  f.paths,
		Rel:    true,
	})
	if err != nil {
		return nil, err
	}
	return append(append([]byte(`{"file_content":`), b...), '}'), nil
}

func (f *stageRelativeExtractor) Start(ctx context.Context, _ ValueContext,
	fc *FileContext, result chan<- any) {
	fc.stageRelative(f.paths).subscribe(f.filter, f.conv, result)
}

func (stageRelativeExtractor) RequiresPipeline() bool { return true }

type errorContentExtractor struct {
	simpleExtractor
	filter *LineFilter
}

func (ec errorContentExtractor) MarshalJSON() ([]byte, error) {
	b, err := json.Marshal(ec.filter)
	if err != nil {
		return nil, err
	}
	return append(append([]byte(`{"error_content":`), b...), '}'), nil
}

func (e errorContentExtractor) RequiresFailure() bool {
	return true
}

func (e errorContentExtractor) RequiresPipeline() bool {
	return true
}

// Start implements ValueExtractor.
func (e errorContentExtractor) Start(ctx context.Context, _ ValueContext,
	fc *FileContext, result chan<- any) {
	r := fc.errorFile()
	go func(result chan<- any, r <-chan errorFile) {
		defer close(result)
		select {
		case ef := <-r:
			if m, ok := e.filter.compile().matchContent(ef.content); ok {
				result <- m
			}
		case <-ctx.Done():
		}
	}(result, r)
}

func isErrors(paths []string) bool {
	if len(paths) == 1 {
		if paths[0] == "_errors" || paths[0] == "_assert" {
			return true
		}
	} else if len(paths) == 2 {
		if paths[0] == "_errors" && paths[1] == "_assert" ||
			paths[0] == "_assert" && paths[1] == "_errors" {
			return true
		}
	}
	return false
}

func validateFilterType(filter FileFilter) {
	switch filter.(type) {
	case *LineFilter:
	case JsonPathFilter:
	default:
		panic(fmt.Sprintf("invalid filter type %T", filter))
	}
}

func MakeRelativeFile(paths []string, stageRel bool,
	conv converter.Converter, filter FileFilter) (ValueExtractor, error) {
	validateFilterType(filter)
	if stageRel {
		if isErrors(paths) {
			// Special case this for a couple of reasons.
			// One is that we want to error out in config validation
			// if someone looks for this in non-failure situations.
			// Another is that we have special-case code for looking
			// at the errors files because we want to ignore the
			// "Caught signal terminate" errors.
			lf, ok := filter.(*LineFilter)
			if !ok {
				return nil, errors.New("only line filters are permitted for error files")
			}
			if conv != nil {
				return nil, errors.New("converters are not supported for error files")
			}
			return errorContentExtractor{
				filter: lf,
			}, nil
		}
		return &stageRelativeExtractor{
			paths:  paths,
			filter: filter,
			conv:   conv,
		}, nil
	}

	return &fileContentExtractor{
		paths:  paths,
		filter: filter,
		conv:   conv,
	}, nil
}

type fileFlagExtractor struct {
	simpleExtractor
	filter FileFilter
	conv   converter.Converter
	flag   string
}

func (f *fileFlagExtractor) MarshalJSON() ([]byte, error) {
	b, err := json.Marshal(struct {
		Filter FileFilter          `json:"filter"`
		Conv   converter.Converter `json:"converter,omitempty"`
		Flag   string              `json:"flag"`
	}{
		Filter: f.filter,
		Conv:   f.conv,
		Flag:   f.flag,
	})
	if err != nil {
		return nil, err
	}
	return append(append([]byte(`{"file_content":`), b...), '}'), nil
}

var getWd = sync.OnceValue(func() string {
	dir, _ := os.Getwd()
	return dir
})

func (f *fileFlagExtractor) Start(ctx context.Context, vc ValueContext,
	fc *FileContext, result chan<- any) {
	for _, arg := range vc.CommandLine() {
		flag, value, ok := strings.Cut(arg, "=")
		if !ok {
			continue
		}
		if flag == f.flag {
			if value == "" {
				close(result)
				return
			}
			if path.IsAbs(value) {
				value = path.Clean(value)
			} else if cwd := getWd(); cwd != "" {
				value = path.Join(cwd, value)
			}
			if fc.files == nil {
				fc.files = make(multiFileResultMap)
			}
			fc.files.getSink([]string{value}).subscribe(
				f.filter, f.conv, result)
			return
		}
	}
	close(result)
}

func MakeFileFlag(flag string,
	conv converter.Converter, filter FileFilter) ValueExtractor {
	return &fileFlagExtractor{
		filter: filter,
		conv:   conv,
		flag:   flag,
	}
}
