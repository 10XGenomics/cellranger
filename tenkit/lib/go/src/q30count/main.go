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

type QualityInfo struct {
	Q30BaseCount   int64 `json:"q30_base_count"`
	TotalBaseCount int64 `json:"total_base_count"`
}

const READ_TO_END = -1

func countQ30Bases(qualityInfo *QualityInfo, fastqReader *fastq.FastqReader, startIndex int, length int) error {
	var read fastq.FastqRecord
	var q30BaseCount int64
	var totalBaseCount int64
	for {
		err := fastqReader.ReadRecord(&read)
		if err == io.EOF {
			break
		} else if err != nil {
			return err
		}

		endIndex := startIndex+length
		if length == READ_TO_END {
			endIndex = len(read.Seq)
		}

		seqQual := read.Qual[startIndex : endIndex]
		for _, qual := range seqQual {
			qualInt := uint32(qual) - 33
			if qualInt >= 30 {
				q30BaseCount += 1
			}
			totalBaseCount += 1
		}
	}
	qualityInfo.Q30BaseCount = q30BaseCount
	qualityInfo.TotalBaseCount = totalBaseCount
	return nil
}


func main() {
	doc := `Tenkit Q30% Counter.

Usage:
    q30count <fastq_path> <output_json> [options]
    q30count -h | --help | --version
Options:
    --read-start-index=NUM  The starting index of the read to consider.
    --read-length=NUM       The length of the barcode in the read.  If not supplied,
                              the Q30 count will be measured from the start index
                              to the end of the read.
    -h --help             Show this message.
    --version             Show version.
`
	opts, _ := docopt.Parse(doc, nil, true, "0.1", false)

	inputFastq := opts["<fastq_path>"].(string)
	outputJsonFile := opts["<output_json>"].(string)

	readStartIndex := 0
	readLength := READ_TO_END
	if value := opts["--read-start-index"]; value != nil {
		if value, err := strconv.Atoi(value.(string)); err != nil {
			log.Fatal("ERROR: --read-start-index argument must be an integer.")
		} else {
			readStartIndex = value
		}
	}
	if value := opts["--read-length"]; value != nil {
		log.Fatal("ERROR: Must supply --bc-length argument")
		if value, err := strconv.Atoi(value.(string)); err != nil {
			log.Fatal("ERROR: --read-length argument must be an integer.")
		} else {
			readLength = value
		}
	}

	// RC doesn't really matter here since we're only counting Q30s
	fastqReader, err := fastq.NewFastqReader(inputFastq, false)
	qualityInfo := &QualityInfo{}
	err = countQ30Bases(qualityInfo, fastqReader, readStartIndex, readLength)
	if err != nil {
		log.Fatal("ERROR: Could not process barcodes from FASTQ file ", err.Error())
	}
	fastqReader.Close()

	qualityInfoData, err := json.Marshal(*qualityInfo)
	if err != nil {
		log.Fatal("ERROR: Failed to serialize read QC information: ", err.Error())
	}
	if err := ioutil.WriteFile(outputJsonFile, qualityInfoData, 0755); err != nil {
		log.Fatal("ERROR: Failed to write read QC file: ", err.Error())
	}
}
