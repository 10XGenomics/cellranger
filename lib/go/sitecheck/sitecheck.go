// Copyright (c) 2024 10x Genomics, Inc. All rights reserved.

// Package sitecheck implements a function for generating sitecheck telemetry file.
package sitecheck

import (
	"context"
	"fmt"
	"io"
	"os"
	"sort"
	"sync"
	"time"
)

type siteCheckSection struct {
	// Title of the section, e.g. "System Information"
	Section string
	// Ideally, a shell script which would generate equivalent output.
	Cmd    string
	Output string
	// The order in which this should show up in the output, ideally.
	Index int
	// This is the last section with this index.
	Last bool
}

func SiteCheck(ctx context.Context, dest io.Writer) error {
	ctx, cancel := context.WithTimeout(ctx, 2*time.Minute)
	defer cancel()
	results := collectData(ctx)

	product := os.Getenv("TENX_PRODUCT")
	version := os.Getenv("TENX_VERSION")
	subcommand := os.Getenv("TENX_SUBCMD")
	copyright := os.Getenv("TENX_COPYRIGHT")
	if copyright == "" {
		copyright = "Copyright (c) 2021 10x Genomics, Inc.  All rights reserved."
	}

	if version != "" {
		fmt.Fprintf(dest, "%s %s (%s)\n",
			product, subcommand, version)
	} else if subcommand != "" {
		fmt.Fprintf(dest, "%s %s (unknown version)\n",
			product, subcommand)
	} else if product != "" {
		fmt.Fprintf(dest, "%s (unknown version)\n", product)
	} else {
		fmt.Fprintln(dest, "10X Genomics sitecheck")
	}
	fmt.Fprintln(dest, copyright)
	fmt.Fprintln(dest,
		"-------------------------------------------------------------------------------")
	fmt.Fprintln(dest, time.Now().Format(time.UnixDate))
	fmt.Fprintln(dest)

	return orderedPrint(ctx, dest, results)
}

// Print sections, ordered by index.
//
// The contract here is is that each checker will have a unique index >= 0,
// and will emit zero or more sections with that index, followed by a section
// with no content and `Last` set to true.
//
// This prints the sections as soon as it can without violating the ordering
// constraint.  If the context times out, it will print any sections which
// arrived out of order.
func orderedPrint(ctx context.Context, dest io.Writer, results <-chan siteCheckSection) error {
	index := 0
	// In order to maintain a consistent ordering for sections,
	// when we receive one that shouldn't be reported yet, we add it to a
	// buffer.  Otherwise we print it immediately - better to get it out there
	// ASAP in case the program is terminated.
	// If we time out, we dump out any buffered sections regardless of whether
	// we're still missing some.
	buffer := make([]siteCheckSection, 0, 24)
	for {
		select {
		case section, ok := <-results:
			if !ok {
				clearBuffer(buffer, -1, dest)
				return nil
			}
			if section.Index > index {
				// Not ready to report this one yet, just stick it in the buffer.
				buffer = append(buffer, section)
			} else if section.Last {
				// Sentinel indicating that we're done with this index.
				index++
				// See if we can remove any buffered items.
				buffer, index = clearBuffer(buffer, index, dest)
			} else {
				section.Print(dest)
			}

		case <-ctx.Done():
			fmt.Fprintln(os.Stderr, "timed out")
			clearBuffer(buffer, -1, dest)
			return ctx.Err()
		}
	}
}

func clearBuffer(buffer []siteCheckSection, index int, dest io.Writer) ([]siteCheckSection, int) {
	if len(buffer) == 0 {
		return buffer, index
	}
	sort.SliceStable(buffer, func(i, j int) bool {
		return buffer[i].Index < buffer[j].Index
	})
	if index < 0 {
		// Just dump everything.
		for _, section := range buffer {
			if !section.Last {
				section.Print(dest)
			}
		}
	} else {
		// Dump sections which aren't out of order.
		for i, section := range buffer {
			if section.Index > index {
				return append(buffer[:0], buffer[i:]...), index
			}
			if section.Last {
				index++
			} else {
				section.Print(dest)
			}
		}
	}
	return buffer[:0], index
}

func (section *siteCheckSection) Print(dest io.Writer) {
	fmt.Fprintln(dest, "=====================================================================")
	fmt.Fprintln(dest, section.Section)
	fmt.Fprintln(dest, section.Cmd)
	fmt.Fprintln(dest, "---------------------------------------------------------------------")
	if section.Output != "" {
		fmt.Fprintln(dest, section.Output)
	}
	fmt.Fprintln(dest, "=====================================================================")
	fmt.Fprintln(dest)
}

type checkFunc func(context.Context, int, chan<- siteCheckSection)

// Run a check.
func (c checkFunc) Do(ctx context.Context,
	wg *sync.WaitGroup, i int,
	results chan<- siteCheckSection) {
	defer wg.Done()
	if ctx.Err() != nil {
		return
	}
	c(ctx, i, results)
	results <- siteCheckSection{Index: i, Last: true}
}

func collectData(ctx context.Context) <-chan siteCheckSection {
	results := make(chan siteCheckSection, 24)
	var wg sync.WaitGroup
	for i, c := range [...]checkFunc{
		checkUname,
		checkDistro,
		checkKernel,
		checkGlibc,
		checkCpu,
		checkMemory,
		checkDisk,
		checkUlimits,
		checkFileLimit,
		checkVmConfig,
		checkCgroups,
		checkContainer,
		checkInit,
		checkSGE,
		checkLSF,
		checkCondor,
		checkBCL2FASTQ,
		checkJava,
		checkRefdata,
		checkSlurm,
		checkMrp,
	} {
		wg.Add(1)
		go c.Do(ctx, &wg, i, results)
	}
	go func(wg *sync.WaitGroup, results chan<- siteCheckSection) {
		wg.Wait()
		close(results)
	}(&wg, results)
	return results
}
