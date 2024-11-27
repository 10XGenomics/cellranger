// Copyright (c) 2024 10x Genomics, Inc. All rights reserved.

// Command sitecheck generates a telemetry report about a system's capabilities.
//
// Use care when making changes here, to avoid collecting privacy-sensitive
// information.
//
// Also, avoid making changes to the formatting of the output, as we have tools
// meant to parse it.  Those tools can adapt, of course, but we're rather not
// have to.
package main

import (
	"context"
	"flag"
	"fmt"
	"os"

	"github.com/10XDev/cellranger/lib/go/sitecheck"
)

func main() {
	showVersion := flag.Bool("version", false, "Show version and exit.")
	flag.Usage = func() {
		product := os.Getenv("TENX_PRODUCT")
		fmt.Fprintf(flag.CommandLine.Output(),
			"Usage: %s sitecheck [options]\n",
			product)
		flag.PrintDefaults()
	}
	flag.Parse()
	if *showVersion {
		version := os.Getenv("TENX_VERSION")
		fmt.Println(version)
		os.Exit(0)
	}
	if err := sitecheck.SiteCheck(context.Background(), os.Stdout); err != nil {
		os.Exit(1)
	}
}
