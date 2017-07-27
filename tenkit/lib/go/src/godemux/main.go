//
// Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
//
// Tenkit sample index demultiplexer.
//
package main

import (
	"encoding/json"
	"github.com/10XGenomics/docopt"
	"io"
	"io/ioutil"
	"log"
	"strconv"
	"tenkit/fastq"
)

type FileGroupJson struct {
	Lane  int               `json:"lane"`
	Files map[string]string `json:"files"`
}

type InputJson struct {
	FileGroups          []*FileGroupJson `json:"file_groups"`
	CommonSampleIndices []string         `json:"common_sample_indices"`
}

func processFastqChunk(siFastqReader *fastq.FastqReader, read1FastqReader *fastq.FastqReader, read2FastqReader *fastq.FastqReader,
	otherIndexFastqReader *fastq.FastqReader, demultFastqWriter *fastq.DemultFastqWriter, summaryCounts map[string]int) error {

	var si fastq.FastqRecord
	var read1 fastq.FastqRecord
	var read2 fastq.FastqRecord
	var otherIndexRead fastq.FastqRecord
	for {
		err := siFastqReader.ReadRecord(&si)
		if err == io.EOF {
			break
		} else if err != nil {
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

		seq := string(si.Seq)
		if _, ok := summaryCounts[seq]; ok {
			summaryCounts[seq] += 1
		} else {
			summaryCounts[fastq.INVALID_SAMPLE_INDEX] += 1
		}

		demultFastqWriter.WriteRecords(&si, &read1, &read2, &otherIndexRead)
	}

	return nil
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
    --rci2read         Use the reverse complement of I2 read.
    -h --help          Show this message.
    --version          Show version.`
	opts, _ := docopt.Parse(doc, nil, true, "0.1", false)

	// Parse command line arguments
	inputJsonFile := opts["<input_json>"].(string)
	demultFastqPath := opts["<demult_fastq_path>"].(string)
	outputJsonFile := opts["<output_json>"].(string)

	demultReadType := fastq.INDEX1_TYPE
	if value := opts["--demult-read"]; value != nil {
		demultReadType = value.(string)
	}
	if demultReadType != fastq.INDEX1_TYPE && demultReadType != fastq.INDEX2_TYPE {
		log.Fatal("ERROR: Demultiplex read type must be an index read type")
	}

	chunk := 0
	if value := opts["--chunk"]; value != nil {
		if value, err := strconv.Atoi(value.(string)); err == nil {
			chunk = value
		}
	}

	rcI2Read := opts["--rci2read"].(bool)

	// Parse input JSON
	bytes, err := ioutil.ReadFile(inputJsonFile)
	if err != nil {
		log.Fatal("ERROR: Unable to read input JSON file: ", err.Error())
	}

	var inputJson *InputJson
	if err := json.Unmarshal(bytes, &inputJson); err != nil {
		log.Fatal("ERROR: Unable to parse input JSON file: ", err.Error())
	}

	// Initialize summary counts
	summaryCounts := map[string]int{
		fastq.INVALID_SAMPLE_INDEX: 0,
	}
	for _, commonSampleIndex := range inputJson.CommonSampleIndices {
		summaryCounts[commonSampleIndex] = 0
	}

	for _, fileGroup := range inputJson.FileGroups {
		if _, ok := fileGroup.Files[demultReadType]; !ok {
			log.Fatal("ERROR: File group does not contain demult read type: ", demultReadType)
		}

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
			if readType == demultReadType {
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
		demultFastqWriter, err := fastq.NewDemultFastqWriter(demultFastqPath, demultReadType, otherIndexReadType,
			fileGroup.Lane, chunk, inputJson.CommonSampleIndices)
		if err != nil {
			log.Fatal("ERROR: Failed to create demultiplexing fastq writer: ", err.Error())
		}

		// Demultiplex fastqs
		if err := processFastqChunk(siFastqReader, read1FastqReader, read2FastqReader, otherIndexFastqReader,
			demultFastqWriter, summaryCounts); err != nil {
			log.Fatal("ERROR: Failed to process fastq chunk: ", err.Error())
		}

		// Close writer
		demultFastqWriter.Close()

		// Close all readers
		siFastqReader.Close()
		read1FastqReader.Close()
		if read2FastqReader != nil {
			read2FastqReader.Close()
		}
		if otherIndexFastqReader != nil {
			otherIndexFastqReader.Close()
		}
	}

	// Write summary counts
	summaryCountsData, err := json.Marshal(summaryCounts)
	if err != nil {
		log.Fatal("ERROR: Failed to serialize summary counts: ", err.Error())
	}
	if err := ioutil.WriteFile(outputJsonFile, summaryCountsData, 0755); err != nil {
		log.Fatal("ERROR: Failed to write summary counts file: ", err.Error())
	}
}
