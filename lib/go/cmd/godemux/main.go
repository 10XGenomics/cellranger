// Copyright (c) 2014 10X Genomics, Inc. All rights reserved.

// Tenkit sample index demultiplexer.
package main

import (
	"context"
	"encoding/json"
	"errors"
	"io"
	"log"
	"os"
	"runtime/pprof"
	"runtime/trace"
	"strconv"

	"github.com/10XDev/cellranger/lib/go/tenkit/fastq"

	docopt "github.com/martian-lang/docopt.go"
)

type FileGroupJson struct {
	Files map[string]string `json:"files"`
	Lane  int               `json:"lane"`
}

// dual-index: CommonSampleIndices may need to be a map[string][]string
// mapping i7 to legal i5 (nil if any matches are OK).
type InputJson struct {
	FileGroups      []*FileGroupJson `json:"file_groups"`
	I7SampleIndices []string         `json:"primary_sample_indices"`
	I5SampleIndices []string         `json:"dual_sample_indices"`
	//CommonSampleIndices map[string][]string `json:"common_sample_indices"`
}

// Either a map of allowed i5 indicies or nil to indicate a wildcard.
type i7i5Lookup map[string]map[string]struct{}

func (m i7i5Lookup) Contains(i7, i5 []byte) bool {
	if mm, ok := m[string(i7)]; !ok {
		return false
	} else if mm == nil {
		return true
	} else {
		_, ok := mm[string(i5)]
		return ok
	}
}

func makeCommonLookup(i7, i5 []string, isDualIndexed bool) i7i5Lookup {
	lookup := make(i7i5Lookup, len(i7))
	if isDualIndexed {
		for i, i7index := range i7 {
			i5index := i5[i]
			if m := lookup[i7index]; m == nil {
				lookup[i7index] = map[string]struct{}{
					i5index: {},
				}
			} else {
				m[i5index] = struct{}{}
			}
		}
	} else {
		for _, i7index := range i7 {
			lookup[i7index] = nil
		}
	}
	return lookup
}

var siEOF = errors.New("End of SI fastq reader")

func processFastqRow(ctx context.Context,
	siFastqReader, read1FastqReader, read2FastqReader, otherIndexFastqReader *fastq.FastqReader,
	demultFastqWriter *fastq.DemultFastqWriter,
	summaryCounts map[string]int,
	i7i5CommonBcs i7i5Lookup) error {
	defer trace.StartRegion(ctx, "processFastqRow").End()
	var si fastq.FastqRecord
	var read1 fastq.FastqRecord
	var read2 fastq.FastqRecord
	var otherIndexRead fastq.FastqRecord
	err := siFastqReader.ReadRecord(&si)
	if err != nil {
		if err == io.EOF {
			return siEOF
		}
		return err
	}

	if err := read1FastqReader.ReadRecord(&read1); err != nil {
		return err
	}

	if read2FastqReader != nil {
		if err = read2FastqReader.ReadRecord(&read2); err != nil {
			return err
		}
	}

	if otherIndexFastqReader != nil {
		if err = otherIndexFastqReader.ReadRecord(&otherIndexRead); err != nil {
			return err
		}
	}

	// dual-index: will also need to read otherIndexRead.Seq and see if
	// it's in the permitted i5s list
	if i7i5CommonBcs.Contains(si.Seq, otherIndexRead.Seq) {
		// i5 is "*" (single index)
		// i5 exists and it is common
		// this include single indexes, being "*"
		summaryCounts[string(si.Seq)] += 1
		return demultFastqWriter.WriteRecords(&si, &read1, &read2, &otherIndexRead, si.Seq)
	} else {
		xSeq := fastq.INVALID_SAMPLE_INDEX
		summaryCounts[string(xSeq)] += 1
		// prepare the record so it is overridden to X by demultFastqWriter
		return demultFastqWriter.WriteRecords(&si, &read1, &read2, &otherIndexRead, xSeq)
	}
}

func processFastqChunk(ctx context.Context,
	siFastqReader, read1FastqReader, read2FastqReader, otherIndexFastqReader *fastq.FastqReader,
	demultFastqWriter *fastq.DemultFastqWriter,
	summaryCounts map[string]int,
	i7i5CommonBcs i7i5Lookup) error {
	defer trace.StartRegion(ctx, "processFastqChunk").End()
	for {
		err := processFastqRow(ctx, siFastqReader,
			read1FastqReader, read2FastqReader, otherIndexFastqReader,
			demultFastqWriter, summaryCounts, i7i5CommonBcs)
		if err != nil {
			if err == siEOF {
				return nil
			}
			return err
		}
	}
}

func parseConfig(inputJsonFile string,
	isDualIndexed bool) (InputJson, map[string]int, i7i5Lookup) {
	defer trace.StartRegion(context.Background(), "parseConfig").End()
	// Parse input JSON
	bytes, err := os.ReadFile(inputJsonFile)
	if err != nil {
		log.Fatal("ERROR: Unable to read input JSON file: ", err.Error())
	}

	var inputJson InputJson
	if err := json.Unmarshal(bytes, &inputJson); err != nil {
		log.Fatal("ERROR: Unable to parse input JSON file: ", err.Error())
	}

	// Initialize summary counts
	summaryCounts := make(map[string]int, 1+len(inputJson.I7SampleIndices))
	summaryCounts[string(fastq.INVALID_SAMPLE_INDEX)] = 0
	for _, i7SampleIndex := range inputJson.I7SampleIndices {
		summaryCounts[i7SampleIndex] = 0
	}
	// Initialize i7 and i5 seeking index-hopping maps
	i7i5CommonBcs := makeCommonLookup(
		inputJson.I7SampleIndices, inputJson.I5SampleIndices, isDualIndexed)

	return inputJson, summaryCounts, i7i5CommonBcs
}

func main() {
	doc := `Tenkit Sample Index Demultiplexer.

Usage:
    godemux <input_json> <demult_fastq_path> <output_json> [options]
    godemux -h | --help | --version
Options:
    --demult-read=VAL  The index read type to demultiplex on.
                         Defaults to I1.
    --chunk=NUM        Chunk number to use for naming FASTQs.
                         Defaults to 0.
    --trace            Turn on performance tracing.
    --profile          Turn on CPU profiling.
    --rci2read         Use the reverse complement of I2 read.
    -h --help          Show this message.
    --version          Show version.`
	opts, _ := docopt.Parse(doc, nil, true, "0.1", false)

	if opts["--trace"].(bool) {
		// Set up profiling.
		if traceFile, err := os.Create("trace.out"); err != nil {
			log.Printf("failed to create trace output file: %v", err)
		} else {
			defer func(traceFile *os.File) {
				if err := traceFile.Close(); err != nil {
					log.Printf("failed to close trace file: %v", err)
				}
			}(traceFile)
			if err := trace.Start(traceFile); err != nil {
				log.Printf("failed to start trace: %v", err)
			} else {
				defer trace.Stop()
			}
		}
	}
	if opts["--profile"].(bool) {
		if profile, err := os.Create("cpu.pprof"); err != nil {
			log.Fatal("could not create CPU profile: ", err)
		} else {
			defer func(profile *os.File) {
				if err := profile.Close(); err != nil {
					log.Printf("failed to close profile file: %v", err)
				}
			}(profile)
			if err := pprof.StartCPUProfile(profile); err != nil {
				log.Println("could not start CPU profile: ", err)
			} else {
				defer pprof.StopCPUProfile()
			}
		}
	}

	// Parse command line arguments
	inputJsonFile := opts["<input_json>"].(string)
	demultFastqPath := opts["<demult_fastq_path>"].(string)
	outputJsonFile := opts["<output_json>"].(string)

	demultReadType := fastq.INDEX1_TYPE
	if value := opts["--demult-read"]; value != nil {
		demultReadType = value.(string)
	}
	if demultReadType != fastq.INDEX1_TYPE &&
		demultReadType != fastq.INDEX2_TYPE &&
		demultReadType != fastq.DUAL_INDEX_TYPE {
		log.Fatal("ERROR: Demultiplex read type must be a recognized index read type")
	}
	// Do we have a dual index library?
	isDualIndexed := (demultReadType == fastq.DUAL_INDEX_TYPE)
	if isDualIndexed {
		demultReadType = fastq.INDEX1_TYPE
	}

	chunk := 0
	if value := opts["--chunk"]; value != nil {
		if value, err := strconv.Atoi(value.(string)); err == nil {
			chunk = value
		}
	}

	rcI2Read := opts["--rci2read"].(bool)

	inputJson, summaryCounts, i7i5CommonBcs := parseConfig(inputJsonFile, isDualIndexed)

	for _, fileGroup := range inputJson.FileGroups {
		if _, ok := fileGroup.Files[demultReadType]; !ok {
			log.Fatal("ERROR: File group does not contain demult read type: ", demultReadType)
		}
		processFileGroup(fileGroup, chunk, &inputJson,
			demultReadType, demultFastqPath,
			isDualIndexed, rcI2Read,
			summaryCounts, i7i5CommonBcs)
	}

	// Write summary counts
	summaryCountsData, err := json.Marshal(summaryCounts)
	if err != nil {
		log.Fatal("ERROR: Failed to serialize summary counts: ", err.Error())
	}
	if err := os.WriteFile(outputJsonFile, summaryCountsData, 0644); err != nil {
		log.Fatal("ERROR: Failed to write summary counts file: ", err.Error())
	}
}

func processFileGroup(fileGroup *FileGroupJson, chunk int,
	inputJson *InputJson,
	demultReadType, demultFastqPath string,
	isDualIndexed, rcI2Read bool,
	summaryCounts map[string]int,
	i7i5CommonBcs i7i5Lookup) {
	ctx, task := trace.NewTask(context.Background(), "processFileGroup")
	defer task.End()
	// Get fastq readers for demultiplexing read + other reads
	var siFastqReader *fastq.FastqReader
	var read1FastqReader *fastq.FastqReader
	var read2FastqReader *fastq.FastqReader
	var otherIndexFastqReader *fastq.FastqReader
	var otherIndexReadType string

	for readType, filename := range fileGroup.Files {
		rc := (readType == fastq.INDEX2_TYPE) && rcI2Read
		fastqReader, err := fastq.NewFastqReader(filename, rc)
		if err != nil {
			log.Fatal("ERROR: Failed to create fastq reader: ", err.Error())
		}
		defer fastqReader.Close()
		if readType == demultReadType ||
			(isDualIndexed && readType == fastq.INDEX1_TYPE) {
			siFastqReader = fastqReader
		} else if readType == fastq.READ1_TYPE {
			read1FastqReader = fastqReader
		} else if readType == fastq.READ2_TYPE {
			read2FastqReader = fastqReader
		} else if readType == fastq.INDEX1_TYPE || readType == fastq.INDEX2_TYPE {
			otherIndexFastqReader = fastqReader
			otherIndexReadType = readType
		} else {
			log.Fatal("Invalid read type: ", readType)
		}
	}

	// Initialize demultiplexing fastq writer
	demultFastqWriter, err := fastq.NewDemultFastqWriter(demultFastqPath,
		demultReadType, otherIndexReadType,
		fileGroup.Lane, chunk, inputJson.I7SampleIndices)
	if err != nil {
		log.Fatal("ERROR: Failed to create demultiplexing fastq writer: ", err.Error())
	}

	// Demultiplex fastqs
	if err := processFastqChunk(ctx,
		siFastqReader, read1FastqReader, read2FastqReader, otherIndexFastqReader,
		demultFastqWriter, summaryCounts, i7i5CommonBcs); err != nil {
		log.Fatal("ERROR: Failed to process fastq chunk: ", err.Error())
	}

	// Close writer
	if err := demultFastqWriter.Close(); err != nil {
		log.Fatal("ERROR: could not flush output writer: ", err.Error())
	}
}
