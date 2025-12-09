// Command telemetry collects and transmits usage statistics for the product.
package main

import (
	"bytes"
	"context"
	"encoding/json"
	"errors"
	"flag"
	"fmt"
	"io/fs"
	"os"
	"path"
	"regexp"
	"runtime/pprof"
	"runtime/trace"
	"sort"
	"strings"
	"time"

	"github.com/10XDev/cellranger/lib/go/telemetry/event"
	"github.com/10XDev/cellranger/lib/go/telemetry/settings"
	"github.com/10XDev/cellranger/lib/go/telemetry/value_context"
)

func getUsage() string {
	return `Usage: telemetry [help] (collect|check|disable|enable|list|show)

collect: Collect telemetry data, if enabled.
check:   Show whether telemetry is currently enabled and
         configuration information.
disable: Disable telemetry collection for this user.
enable:  Enable telemetry collection for this user.
list:    List files containing saved telemetry data for this product.
show:    Display contents of saved telemetry data for this product.

For more information about what data is collected and how it's used, visit
` + settings.PrivacyStatementUrl()
}

func main() {
	f := flag.NewFlagSet("telemetry", flag.ExitOnError)
	f.Usage = func() {
		fmt.Fprintln(os.Stderr, getUsage())
	}
	var cpuProfile, traceFile string
	f.StringVar(&cpuProfile, "profile", "",
		"Save CPU profiling information for this process to the given file.")
	f.StringVar(&traceFile, "trace", "",
		"Save an execution trace to the given file.")

	if err := f.Parse(os.Args[1:]); err != nil {
		panic(err)
	}
	if f.NArg() < 1 {
		f.Usage()
		os.Exit(1)
	}

	os.Exit(func(f *flag.FlagSet) int {
		if cpuProfile != "" {
			f, err := os.Create(cpuProfile)
			if err != nil {
				fmt.Fprintln(os.Stderr,
					"Error opening profile destination:", err)
				os.Exit(7)
			}
			if err := pprof.StartCPUProfile(f); err != nil {
				f.Close()
				fmt.Fprintln(os.Stderr, "Error recording CPU profile:", err)
				os.Exit(7)
			}
			defer func(f *os.File) {
				pprof.StopCPUProfile()
				if err := f.Close(); err != nil {
					fmt.Fprintln(os.Stderr, "Error recording CPU profile:", err)
				}
			}(f)
		}
		if traceFile != "" {
			f, err := os.Create(traceFile)
			if err != nil {
				fmt.Fprintln(os.Stderr,
					"Error opening trace destination:", err)
				os.Exit(7)
			}
			if err := trace.Start(f); err != nil {
				f.Close()
				fmt.Fprintln(os.Stderr, "Error recording trace:", err)
				os.Exit(7)
			}
			defer func(f *os.File) {
				trace.Stop()
				if err := f.Close(); err != nil {
					fmt.Fprintln(os.Stderr, "Error recording trace:", err)
				}
			}(f)
		}

		switch f.Arg(0) {
		case "help":
			return doHelp(f.Args()[1:])
		case "collect":
			return doCollect(f.Args()[1:])
		case "disable":
			return doDisable(f.Args()[1:])
		case "enable":
			return doEnable(f.Args()[1:])
		case "check":
			return doCheck()
		case "list":
			return doList(f.Args()[1:])
		case "show":
			return doShow(f.Args()[1:])
		default:
			fmt.Fprintln(os.Stderr, getUsage())
			return 1
		}
	}(f))
}

func doHelp(args []string) int {
	if len(args) == 0 {
		fmt.Fprintln(os.Stderr, getUsage())
		return 0
	}
	switch args[0] {
	case "collect":
		value_context.ParseFlags([]string{"--help"})
	case "disable":
		fmt.Fprintln(os.Stderr, disableUsage())
	case "enable":
		fmt.Fprintln(os.Stderr, enableUsage())
	case "check":
		fmt.Fprintln(os.Stderr, `Show the whether telemetry is enabled.`)
	case "list":
		fmt.Fprintln(os.Stderr, `List saved telemetry files.`)
	case "show":
		return doShow([]string{"-h"})
	default:
		fmt.Fprintln(os.Stderr, getUsage())
		return 1
	}
	return 0
}

func disableUsage() string {
	conf, _ := settings.ConfigDir()
	cache, _ := settings.CacheDir()
	if conf == "" {
		conf = "~/.config/tenx/telemetry"
	}
	if cache == "" {
		cache = "~/.cache/tenx/telemetry"
	}
	return `Usage: telemetry disable [upload|update]

This will write a file to ` + conf + `
to disable part or all of telemetry collection.

With no arguments, all telemetry collection will be disabled.

Disabling 'update' will prevent the tool from checking for an updated
configuration file.

Disabling 'upload' will prevent collected telemetry data from being sent to
10X Genomics, but will still collect it locally into
` + cache + `

Disabling both upload and update will prevent any communication between this
tool and 10XGenomics, but allow one to view the data that would have been sent.

Note that telemetry may also be disabled by setting the environment variables
` + settings.DisableEnvVar + `, ` + settings.DisableEnvVar + `_UPDATE,
or ` + settings.DisableEnvVar + `_UPLOAD to non-empty values.
Config updates or upload can be disabled for all users by creating the files
/etc/tenx/telemetry/disable_update or /etc/tenx/telemetry/disable_upgrade.`
}

func enableUsage() string {
	conf, _ := settings.ConfigDir()
	if conf == "" {
		conf = "~/.config/tenx/telemetry"
	}
	return `Usage: telemetry enable [upload|update]

This will remove the file in ` + conf + `
that would disable part or all of telemetry collection.
Telemetry may remain disabled based on environment variables or
system-wide settings in /etc/tenx/telemetry/.

With no arguments, all telemetry collection will be enabled.

Enabling 'update' will allow the tool to check for an updated
configuration file.

Allowing 'upload' will permit collected telemetry data from being sent
to 10X Genomics.

For more information about what data is collected and how it's used, visit
` + settings.PrivacyStatementUrl()
}

func doCollect(args []string) int {
	ctx, cancel := context.WithTimeout(context.Background(),
		15*time.Second)
	defer cancel()
	f := value_context.ParseFlags(args)
	status := settings.GetTelemetryStatus()
	if status.DisableAll {
		return 0
	}
	config, err := status.GetConfig(ctx, f.Verbose)
	if err != nil {
		if f.Verbose {
			fmt.Fprintln(os.Stderr, err)
		}
		return 1
	}
	c, err := config.Compile()
	if err != nil {
		if f.Verbose {
			fmt.Fprintln(os.Stderr, err)
		}
		return 1
	}
	vc, err := f.MakeValueContext(settings.CountRecentCacheFiles)
	if err != nil {
		if f.Verbose {
			fmt.Fprintln(os.Stderr, err)
		}
		return 1
	}
	if err := status.GatherEvents(
		ctx, vc.PipestanceDir(),
		event.Collect(ctx, &c, vc)); err != nil {
		if f.Verbose {
			fmt.Fprintln(os.Stderr, err)
		}
		return 1
	}
	return 0
}

func doDisable(args []string) int {
	var all, upload, update bool
	if len(args) == 0 {
		all = true
	}
	for _, arg := range args {
		if arg == "-h" || arg == "--help" {
			fmt.Println(disableUsage())
			return 0
		}
		if strings.EqualFold(arg, "all") {
			all = true
		} else if strings.EqualFold(arg, "upload") {
			upload = true
		} else if strings.EqualFold(arg, "update") {
			update = true
		} else {
			fmt.Println(disableUsage())
			return 1
		}
	}
	if err := settings.Disable(all, upload, update); err != nil {
		fmt.Fprintln(os.Stderr, "Failed to disable telemetry:\n", err)
		return 3
	}
	return 0
}

func doEnable(args []string) int {
	var all, upload, update bool
	if len(args) == 0 {
		all = true
	}
	for _, arg := range args {
		if arg == "-h" || arg == "--help" {
			fmt.Println(enableUsage())
			return 0
		}
		if strings.EqualFold(arg, "all") {
			all = true
		} else if strings.EqualFold(arg, "upload") {
			upload = true
		} else if strings.EqualFold(arg, "update") {
			update = true
		} else {
			fmt.Println(disableUsage())
			return 1
		}
	}
	if err := settings.Enable(all, upload, update); err != nil {
		fmt.Fprintln(os.Stderr, "Failed to disable telemetry:\n", err)
		return 3
	}
	return 0
}

func doCheck() int {
	status := settings.GetTelemetryStatus()
	if status.DisableAll {
		fmt.Println("Telemetry collection is disabled.")
		return 0
	} else {
		if status.DisableConfigUpdate {
			fmt.Println("Telemetry configuration updates are disabled.")
		}
		if status.DisableUpload {
			fmt.Println("Telemetry uploads are disabled.\n" +
				"Data is still collected locally.")
		}
	}
	if conf, err := settings.ConfigDir(); err != nil {
		fmt.Fprintln(os.Stderr,
			"Could not locate configuration directory:", err)
	} else {
		fmt.Println("Configuration data is stored in", conf)
	}
	if cache, err := settings.ProductCacheDir(); err != nil {
		if cache, err := settings.CacheDir(); cache == "" {
			fmt.Fprintln(os.Stderr,
				"Could not locate cache directory:", err)
		} else {
			fmt.Println("Collected telemetry data is being saved to", cache)
		}
	} else {
		fmt.Println("Collected telemetry data is being saved to", cache)
	}
	ctx, cancel := context.WithTimeout(context.Background(), time.Second*5)
	config, err := status.GetConfig(ctx, true)
	cancel()
	if err != nil {
		fmt.Fprintln(os.Stderr, "Error retrieving configuration:", err)
		return 1
	}
	if config.ConfigTime.IsZero() {
		fmt.Printf("Configuration for %s v%s version %d\n",
			config.Product, config.Version, config.ConfigVersion)
	} else {
		fmt.Printf("Configuration for %s v%s version %d, dated %s\n",
			config.Product, config.Version, config.ConfigVersion,
			config.ConfigTime.Format(time.RFC1123))
	}
	return 0
}

func doList(args []string) int {
	if len(args) > 0 && (args[0] == "-h" || args[0] == "--help") {
		return doHelp([]string{"list"})
	}
	listing, err := settings.ListCache()
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 3
	}
	for _, f := range listing {
		fmt.Println(f)
	}
	return 0
}

func doShow(args []string) int {
	flagset := flag.NewFlagSet("telemetry show", flag.ExitOnError)
	var limit int
	var eventType string
	var asJson bool
	flagset.IntVar(&limit, "limit", 0,
		"Limit to at most N entries.")
	flagset.StringVar(&eventType, "type", "",
		"Limit to events of the given type.")
	flagset.BoolVar(&asJson, "json", false,
		"Output as json, one entry per line.")
	if err := flagset.Parse(args); err != nil {
		panic(err)
	}
	listing, err := settings.ListCache()
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 3
	}
	eventRe := regexp.MustCompile(`^(\w+)_\d{4}-\d{2}-\d{2}T\w+.json$`)
	if eventType != "" {
		filteredListing := make([]string, 0, len(listing))
		for _, f := range listing {
			m := eventRe.FindStringSubmatch(path.Base(f))
			if len(m) > 1 && m[1] == eventType {
				filteredListing = append(filteredListing, f)
			}
		}
		if limit > 0 && limit < len(filteredListing) {
			filteredListing = filteredListing[len(filteredListing)-limit:]
		}
		showEvents(filteredListing, asJson)
	} else {
		groups := make(map[string][]string)
		for _, f := range listing {
			m := eventRe.FindStringSubmatch(path.Base(f))
			if len(m) > 1 {
				groups[m[1]] = append(groups[m[1]], f)
			}
		}
		groupNames := make([]string, 0, len(groups))
		for g := range groups {
			groupNames = append(groupNames, g)
		}
		sort.Strings(groupNames)
		for i, g := range groupNames {
			if !asJson {
				if i != 0 {
					fmt.Println()
				}
				fmt.Printf("%s:\n", g)
			}
			list := groups[g]
			if limit > 0 && limit < len(list) {
				list = list[len(list)-limit:]
			}
			showEvents(list, asJson)
		}
	}
	return 0
}

func showEvents(list []string, asJson bool) {
	for i, f := range list {
		b, err := os.ReadFile(f)
		if err != nil {
			if !errors.Is(err, fs.ErrNotExist) {
				fmt.Fprintln(os.Stderr, "Could not read file:", err)
			}
			continue
		}
		if asJson {
			var buf bytes.Buffer
			if err := json.Indent(&buf, b, "", "  "); err != nil {
				fmt.Fprintln(os.Stderr, "Invalid json found in",
					f, ":", err)
				continue
			}
			os.Stdout.Write(buf.Bytes())
			os.Stdout.Write([]byte{'\n'})
		} else {
			var ev event.Event
			if err := json.Unmarshal(b, &ev); err != nil {
				fmt.Fprintln(os.Stderr, "Invalid json found in",
					f, ":", err)
				continue
			}
			if i != 0 {
				fmt.Println()
			}
			if err := ev.Print(os.Stdout); err != nil {
				panic(err)
			}
		}
	}
}
