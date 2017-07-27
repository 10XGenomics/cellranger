package main

import (
	"encoding/json"
	"github.com/10XGenomics/docopt"
	"io"
	"io/ioutil"
	"log"
	"strconv"
	"tenkit/barcode"
	"tenkit/fastq"
)

type BarcodeInfo struct {
	BarcodeQ30BaseCount        int64          `json:"q30_base_count"`
	ExactBarcodeMatchCount     int            `json:"exact_match_count"`
	CorrectedBarcodeMatchCount int            `json:"corrected_match_count"`
	TotalBaseCount             int64          `json:"total_base_count"`
	TotalReadCount             int            `json:"total_read_count"`
	ValidatedBarcodeCounts     map[string]int `json:"validated_sequence_counts"`
	MeanBaseQualityScore       float64        `json:"mean_base_qscore"`
}

func countBarcodes(barcodeInfo *BarcodeInfo, fastqReader *fastq.FastqReader, barcodeCounter *barcode.BarcodeCounter, startIndex int, length int) error {
	var read fastq.FastqRecord
	for {
		err := fastqReader.ReadRecord(&read)
		if err == io.EOF {
			break
		} else if err != nil {
			return err
		}

		seq := string(read.Seq)
		bclen := length
		if len(seq) < length {
			bclen = len(seq)
		}
		barcode := seq[startIndex : startIndex+bclen]
		barcodeQual := read.Qual[startIndex : startIndex+bclen]
		barcodeCounter.CountBarcode(barcode, barcodeQual)
	}
	barcodeInfo.BarcodeQ30BaseCount = barcodeCounter.GetQ30BaseCount()
	barcodeInfo.TotalBaseCount = barcodeCounter.GetTotalBaseCount()
	barcodeInfo.TotalReadCount = barcodeCounter.GetTotalCount()
	barcodeInfo.ExactBarcodeMatchCount = barcodeCounter.GetTotalCount() - barcodeCounter.GetMismatchCount()
	barcodeInfo.MeanBaseQualityScore = barcodeCounter.GetMeanBaseQualityScore()
	return nil
}

func countCorrectedBarcodes(barcodeInfo *BarcodeInfo, fastqReader *fastq.FastqReader, validator *barcode.BarcodeValidator, startIndex int, length int) error {
	var read fastq.FastqRecord
	correctedBarcodeCounts := make(map[string]int)
	validBarcodeCount := 0
	for {
		err := fastqReader.ReadRecord(&read)
		if err == io.EOF {
			break
		} else if err != nil {
			return err
		}
		seq := string(read.Seq)
		bclen := length
		if len(seq) < length {
			bclen = len(seq)
		}
		barcode := seq[startIndex : startIndex+bclen]
		correctedBarcode, valid := validator.ValidateBarcode(barcode, read.Qual[startIndex:startIndex+bclen])
		correctedBarcodeString := string(correctedBarcode)
		if valid {
			if _, ok := correctedBarcodeCounts[correctedBarcodeString]; !ok {
				correctedBarcodeCounts[correctedBarcodeString] = 1
			} else {
				correctedBarcodeCounts[correctedBarcodeString] += 1
			}
			validBarcodeCount += 1
		}
	}
	barcodeInfo.CorrectedBarcodeMatchCount = validBarcodeCount
	barcodeInfo.ValidatedBarcodeCounts = correctedBarcodeCounts
	return nil
}

func main() {
	doc := `Tenkit Barcode QC Analyzer.

Usage:
    barcodeqc <fastq_path> <output_json> [options]
    barcodeqc -h | --help | --version
Options:
    --bc-start-index=NUM  The starting index of the barcode in the read.
    --bc-length=NUM       The length of the barcode in the read.
    --rc                  Process the reverse complement of the read.
    --whitelist=PATH      The path to the barcode whitelist.
    -h --help             Show this message.
    --version             Show version.
`
	opts, _ := docopt.Parse(doc, nil, true, "0.1", false)

	inputFastq := opts["<fastq_path>"].(string)
	outputJsonFile := opts["<output_json>"].(string)

	whitelist := opts["--whitelist"].(string)
	rcRead := opts["--rc"].(bool)
	bcStartIndex := 0
	bcLength := 16
	if value := opts["--bc-start-index"]; value != nil {
		if value, err := strconv.Atoi(value.(string)); err != nil {
			log.Fatal("ERROR: --bc-start-index argument must be an integer.")
		} else {
			bcStartIndex = value
		}
	}
	if value := opts["--bc-length"]; value == nil {
		log.Fatal("ERROR: Must supply --bc-length argument")
	} else {
		if value, err := strconv.Atoi(value.(string)); err != nil {
			log.Fatal("ERROR: --bc-length argument must be an integer.")
		} else {
			bcLength = value
		}
	}

	// parse whitelist
	_, err := ioutil.ReadFile(whitelist)
	if err != nil {
		log.Fatal("ERROR: Unable to read barcode whitelist: ", err.Error())
	}

	barcodeCounter := barcode.NewBarcodeCounter(whitelist)
	fastqReader, err := fastq.NewFastqReader(inputFastq, rcRead)
	if err != nil {
		log.Fatal("ERROR: Unable to open FASTQ file: ", err.Error())
	}
	barcodeInfo := &BarcodeInfo{}
	err = countBarcodes(barcodeInfo, fastqReader, barcodeCounter, bcStartIndex, bcLength)
	if err != nil {
		log.Fatal("ERROR: Could not process barcodes from FASTQ file ", err.Error())
	}
	fastqReader.Close()

	// accuracy vs performance tradeoff.  Need to read through FASTQ barcodes again
	// in order to rescue mismatches based on the BarcodeValidator distribution, and
	// agree with the output counts.  If not needed, change this by writing
	// the barcodes in countBarcodes instead.
	fastqReader, err = fastq.NewFastqReader(inputFastq, rcRead)
	if err != nil {
		log.Fatal("ERROR: Unable to open FASTQ file: ", err.Error())
	}
	validator := barcodeCounter.GetBarcodeValidator(1, 0.975)
	err = countCorrectedBarcodes(barcodeInfo, fastqReader, validator, bcStartIndex, bcLength)
	if err != nil {
		log.Fatal("ERROR: Unable to determine corrected Q30 count")
	}

	barcodeInfoData, err := json.Marshal(*barcodeInfo)
	if err != nil {
		log.Fatal("ERROR: Failed to serialize barcode QC information: ", err.Error())
	}
	if err := ioutil.WriteFile(outputJsonFile, barcodeInfoData, 0755); err != nil {
		log.Fatal("ERROR: Failed to write barcode QC file: ", err.Error())
	}
}
