// Command telemetry_converter is used to convert from yaml to json.
//
// It will also validate the config file.
//
// The command does not accept any command-line flags.
// It reads from stdin and writes to stdout.
package main

import (
	"bytes"
	"encoding/json"
	"flag"
	"fmt"
	"io"
	"os"
	"strconv"
	"time"

	"github.com/10XDev/cellranger/lib/go/telemetry/config"
	"gopkg.in/yaml.v3"
)

func main() {
	var hashes bool
	flag.BoolVar(&hashes, "hash", false,
		"Output group hashes instead of json.")
	flag.Parse()
	in, err := io.ReadAll(os.Stdin)
	if err != nil {
		fmt.Fprintln(os.Stderr, "error reading from stdin:", err)
		os.Exit(1)
	}
	if len(in) == 0 {
		fmt.Fprintln(os.Stderr, "empty input")
		os.Exit(1)
	}
	var jsonify any
	if err := yaml.Unmarshal(in, &jsonify); err != nil {
		fmt.Fprintln(os.Stderr, "invalid yaml:", err)
		os.Exit(1)
	}
	b, err := json.Marshal(jsonify)
	if err != nil {
		fmt.Fprintln(os.Stderr, "could not marshal yaml to json:", err)
		os.Exit(1)
	}
	dec := json.NewDecoder(bytes.NewReader(b))
	dec.DisallowUnknownFields()
	var conf config.Config
	if err := dec.Decode(&conf); err != nil {
		fmt.Fprintln(os.Stderr, "malformed config:", err)
		os.Exit(3)
	}
	if t, ok := findBuildTime(); ok {
		conf.ConfigTime = t
	}
	if c, err := conf.Compile(); err != nil {
		fmt.Fprintln(os.Stderr, "invalid config:", err)
		os.Exit(4)
	} else if hashes {
		groups := make(config.SortedAnyMap, len(c.Groups))
		for gn, gc := range c.Groups {
			h := gc.Hash()
			groups[gn] = h[:]
		}
		b, err := json.MarshalIndent(groups, "", "  ")
		if err != nil {
			fmt.Fprintln(os.Stderr, "Error marshaling result:", err)
			os.Exit(5)
		}
		os.Stdout.Write(b)
		return
	}
	if b, err := json.Marshal(&conf); err != nil {
		fmt.Fprintln(os.Stderr, "could not marshal config to json:", err)
		os.Exit(1)
	} else {
		os.Stdout.Write(b)
	}
}

func findBuildTime() (time.Time, bool) {
	b, err := os.ReadFile("bazel-out/volatile-status.txt")
	if err != nil {
		return time.Time{}, false
	}
	i := bytes.Index(b, []byte("BUILD_TIMESTAMP "))
	if i < 0 {
		return time.Time{}, false
	}
	ts, _, _ := bytes.Cut(b[i+len("BUILD_TIMESTAMP "):], []byte("\n"))
	ti, err := strconv.ParseInt(string(ts), 10, 64)
	if err != nil {
		return time.Time{}, false
	}
	return time.Unix(ti, 0), true
}
