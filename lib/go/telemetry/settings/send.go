package settings

import (
	"bytes"
	"context"
	"crypto/sha256"
	"encoding/base64"
	"encoding/json"
	"errors"
	"fmt"
	"io"
	"net/http"
	"os"
	"path"
	"regexp"
	"runtime/trace"
	"strings"
	"sync"
	"time"

	"github.com/10XDev/cellranger/lib/go/telemetry/event"
)

type eventPayload struct {
	content []byte
	product string
	version string
	name    string
}

func (status *TelemetryStatus) GatherEvents(ctx context.Context, psDir string,
	events <-chan event.Event) error {
	ctx, task := trace.NewTask(ctx, "GatherEvents")
	defer task.End()
	errs, errc := func() ([]error, <-chan error) {
		defer trace.StartRegion(ctx, "GatherEvents.forwarding").End()
		var wg sync.WaitGroup
		var errs []error
		errc := make(chan error, 4)
		// Separate save and send channels so they can run in parallel.
		var save chan eventPayload
		cache, _ := CacheDir()
		if cache != "" || psDir != "" {
			if cache != "" {
				if err := os.MkdirAll(cache, 0755); err != nil {
					errs = append(errs,
						fmt.Errorf("creating save location: %w", err))
					cache = ""
				}
			}
			if psDir != "" {
				if err := os.MkdirAll(psDir, 0777); err != nil {
					errs = append(errs,
						fmt.Errorf("creating save location: %w", err))
					psDir = ""
				}
			}
			if cache != "" || psDir != "" {
				save = make(chan eventPayload, 4)
				defer close(save)
				wg.Add(1)
				go saveEvents(ctx, cache, psDir, save, errc, &wg)
			}
		}
		for event := range events {
			b, err := json.Marshal(event)
			if err != nil {
				errs = append(errs, err)
				continue
			}
			p := eventPayload{
				content: b,
				product: event.Product,
				version: event.Version,
				name:    event.Name,
			}
			if !status.DisableUpload {
				wg.Add(1)
				// Upload in parallel, because network requests to the same URL
				// can be combined efficiently.
				go uploadEvent(ctx, p, errc, &wg)
			}
			if cache != "" || psDir != "" {
				// But save serially because they cannot.
				select {
				case save <- p:
				case <-ctx.Done():
				}
			}
		}
		go func() {
			wg.Wait()
			close(errc)
		}()
		return errs, errc
	}()
	for err := range errc {
		errs = append(errs, err)
	}
	return errors.Join(errs...)
}

func uploadEvent(ctx context.Context, event eventPayload,
	errc chan<- error, wg *sync.WaitGroup) {
	defer wg.Done()
	defer trace.StartRegion(ctx, "uploadEvent").End()
	if ctx.Err() != nil {
		return
	}
	ctx, cancel := context.WithTimeout(ctx, time.Second*5)
	defer cancel()
	req, err := http.NewRequestWithContext(ctx, http.MethodPost,
		telemetryUploadUrl, bytes.NewReader(event.content))
	if err != nil {
		select {
		case errc <- fmt.Errorf(
			"creating request telementry upload request: %w", err):
		default:
			// Enough errors to overwhelm reporting anyway.
		}
		return
	}
	req.Header.Set("Content-Type", "application/json")
	// This isn't a "key", really.
	// There's no way to protect an actually secret key in this binary, since it
	// is being distributed to users.
	// It's set here simply to prevent accidental or casual abuse of the upload
	// API endpoint (in addition to server-side DDoS protections etc.)
	pvHash := sha256.Sum256(append([]byte(event.product), event.version...))
	req.Header.Set("X-Api-Key",
		base64.RawStdEncoding.EncodeToString(pvHash[:]))
	resp, err := http.DefaultClient.Do(req)
	if err != nil {
		select {
		case errc <- fmt.Errorf("sending telementry: %w", err):
		default:
			// Enough errors to overwhelm reporting anyway.
		}
		return
	}
	b, _ := io.ReadAll(resp.Body)
	resp.Body.Close()
	b = bytes.TrimSpace(b)
	if resp.StatusCode < http.StatusOK || resp.StatusCode >= http.StatusMultipleChoices {
		if len(b) > 0 {
			errc <- fmt.Errorf("sending telementry: http %d: %s\n%s",
				resp.StatusCode, resp.Status, string(b))
		} else {
			errc <- fmt.Errorf("sending telementry: http %d: %s",
				resp.StatusCode, resp.Status)
		}
	}
}

func saveEvents(ctx context.Context, cache, psDir string, events <-chan eventPayload,
	errc chan<- error, wg *sync.WaitGroup) {
	defer trace.StartRegion(ctx, "saveEvents").End()
	defer wg.Done()
	var cleanupDirs map[string]struct{}
	var extrasDir string
	if psDir != "" {
		extrasDir = path.Join(psDir, "extras", "telemetry")
		for _, d := range [...]string{path.Dir(extrasDir), extrasDir} {
			if err := os.Mkdir(d, 0777); err != nil && !errors.Is(err, os.ErrExist) {
				select {
				case errc <- fmt.Errorf(
					"creating save location in pipestance directory file: %w",
					err):
				default:
					// Enough errors to overwhelm reporting anyway.
				}
			}
		}
	}
	for event := range events {
		if ctx.Err() != nil {
			return
		}
		var fn string
		if cache != "" {
			fn = saveEventCache(cache, event, errc)
			if fn != "" {
				if cleanupDirs == nil {
					cleanupDirs = make(map[string]struct{}, 1)
				}
				cleanupDirs[path.Dir(fn)] = struct{}{}
			}
		}
		if psDir != "" {
			saveEventInPipestance(extrasDir, path.Base(fn), event, errc)
		}
	}
	for dir := range cleanupDirs {
		cleanOldEvents(dir, errc)
	}
}

func saveEventCache(cache string, event eventPayload, errc chan<- error) string {
	dir := path.Join(cache, event.product, event.version)
	if err := os.MkdirAll(dir, 0755); err != nil {
		select {
		case errc <- fmt.Errorf("creating telementry directory: %w", err):
		default:
			// Enough errors to overwhelm reporting anyway.
		}
	}
	f, err := os.CreateTemp(dir,
		event.name+"_"+
			// Include the date so that we can reliably remove old files.
			// Convert `:` to `_` so that the file paths remain valid
			// on Windows, because sometimes people browse linux directories
			// from windows via WSL or SMB.
			strings.Replace(time.Now().UTC().Format(time.RFC3339),
				":", "_", 2)+
			"*.json")
	if err != nil {
		select {
		case errc <- fmt.Errorf("creating telementry file: %w", err):
		default:
			// Enough errors to overwhelm reporting anyway.
		}
		return ""
	}
	defer f.Close()
	if _, err := f.Write(event.content); err != nil {
		select {
		case errc <- fmt.Errorf("writing telementry file: %w", err):
		default:
			// Enough errors to overwhelm reporting anyway.
		}
	}
	if err := f.Close(); err != nil {
		select {
		case errc <- fmt.Errorf("closing telementry file: %w", err):
		default:
			// Enough errors to overwhelm reporting anyway.
		}
	}
	return f.Name()
}

var eventFileDateRe = regexp.MustCompile(`([^/]*)_(\d{4}-\d{2}-\d{2}T\d{2}_\d{2}_\d{2}Z).*\.json`)

func parseEventFileName(fn string) (string, time.Time) {
	match := eventFileDateRe.FindStringSubmatch(fn)
	if len(match) != 3 {
		return "", time.Time{}
	}
	dateString := strings.Replace(match[2], "_", ":", 2)
	t, err := time.Parse(time.RFC3339, dateString)
	if err != nil {
		return match[1], time.Time{}
	}
	return match[1], t
}

const maxEventTTL = 30 * 24 * time.Hour

func cleanOldEvents(dir string, errc chan<- error) {
	entries, err := os.ReadDir(dir)
	if err != nil {
		select {
		case errc <- fmt.Errorf(
			"checking telementry output directory: %w", err):
		default:
		}
		return
	}
	for _, entry := range entries {
		if entry.IsDir() {
			continue
		}
		_, t := parseEventFileName(entry.Name())
		if !t.IsZero() && time.Since(t) > maxEventTTL {
			if err := os.Remove(path.Join(dir, entry.Name())); err != nil {
				select {
				case errc <- fmt.Errorf(
					"removing old telementry event: %w", err):
				default:
				}
			}
		}
	}
}

func saveEventInPipestance(dir, fn string, event eventPayload, errc chan<- error) string {
	if err := os.MkdirAll(dir, 0755); err != nil {
		select {
		case errc <- fmt.Errorf("creating telementry directory: %w", err):
		default:
			// Enough errors to overwhelm reporting anyway.
		}
	}
	var f *os.File
	var err error
	if fn == "" && fn != "." {
		f, err = os.CreateTemp(dir,
			event.name+"_"+
				// Include the date so that we can reliably remove old files.
				// Convert `:` to `_` so that the file paths remain valid
				// on Windows, because sometimes people browse linux directories
				// from windows via WSL or SMB.
				strings.Replace(time.Now().UTC().Format(time.RFC3339),
					":", "_", -1)+
				"*.json")
		if err == nil {
			// Ignore errors; not much we can do about it.
			_ = os.Chmod(f.Name(), 0644)
		}
	} else {
		f, err = os.Create(path.Join(dir, fn))
	}
	if err != nil {
		select {
		case errc <- fmt.Errorf("creating telementry file: %w", err):
		default:
			// Enough errors to overwhelm reporting anyway.
		}
		return ""
	}
	defer f.Close()
	if _, err := f.Write(event.content); err != nil {
		select {
		case errc <- fmt.Errorf("writing telementry file: %w", err):
		default:
			// Enough errors to overwhelm reporting anyway.
		}
	}
	if err := f.Close(); err != nil {
		select {
		case errc <- fmt.Errorf("closing telementry file: %w", err):
		default:
			// Enough errors to overwhelm reporting anyway.
		}
	}
	return f.Name()
}
