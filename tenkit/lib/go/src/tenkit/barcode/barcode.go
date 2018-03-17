package barcode

// Copyright (c) 2016 10x Genomics, Inc. All rights reserved.

import (
	"bufio"
	"encoding/json"
	"fmt"
	"math"
	"os"
	"sort"
	"strings"
)

type BarcodeCounter struct {
	whitelistExists bool
	whitelistArray  []string
	barcodeCounts   map[string]int
	mismatchCount   int
	totalCount      int
	totalBaseCount  int64
	q30BaseCount    int64
	baseQualSum     int64
}

const BLANK_WHITELIST_FILE = "none"

func NewBarcodeCounter(whitelistFile string) *BarcodeCounter {
	if whitelistFile == BLANK_WHITELIST_FILE {
		return &BarcodeCounter{
			barcodeCounts:   make(map[string]int),
			whitelistArray:  []string{},
			whitelistExists: false,
			mismatchCount:   0, totalCount: 0}
	}
	file, _ := os.Open(whitelistFile)
	defer file.Close()
	reader := bufio.NewReader(file)
	barcodeCounts := make(map[string]int)
	whitelistArray := []string{}
	for {
		line, err := reader.ReadString('\n')
		if err != nil {
			break
		}
		bc := strings.TrimSpace(line)
		if !strings.HasPrefix(line, "#") {
			whitelistArray = append(whitelistArray, bc)
			barcodeCounts[bc] = 0
		}
	}
	return &BarcodeCounter{
		barcodeCounts:   barcodeCounts,
		whitelistArray:  whitelistArray,
		whitelistExists: true,
		mismatchCount:   0,
		totalCount:      0}
}

func (c *BarcodeCounter) CountBarcode(barcode string, barcodeQual []byte) {
	if c.whitelistExists {
		if num, ok := c.barcodeCounts[barcode]; !ok {
			c.mismatchCount += 1
		} else {
			c.barcodeCounts[barcode] = num + 1
		}
	} else if !c.whitelistExists && len(barcode) > 1 {
		if num, ok := c.barcodeCounts[barcode]; !ok {
			c.barcodeCounts[barcode] = 1
		} else {
			c.barcodeCounts[barcode] = num + 1
		}
	} else {
		c.mismatchCount += 1
	}
	c.totalCount += 1

	for _, qual := range barcodeQual {
		qualInt := uint32(qual) - 33
		c.baseQualSum += int64(qualInt)
		if qualInt >= 30 {
			c.q30BaseCount += 1
		}
		c.totalBaseCount += 1
	}
}

func (c *BarcodeCounter) GetSortedBarcodeList() []string {
	barcodeArray := []string{}
	for barcode, _ := range c.barcodeCounts {
		barcodeArray = append(barcodeArray, barcode)
	}
	sort.Strings(barcodeArray)
	return barcodeArray
}

func (c *BarcodeCounter) GetSortedBarcodeCount() []int {
	barcodeList := c.GetSortedBarcodeList()
	sortedBarcodeCounts := make([]int, len(barcodeList))
	for idx, barcode := range barcodeList {
		sortedBarcodeCounts[idx] = c.barcodeCounts[barcode]
	}
	return sortedBarcodeCounts
}

func (c *BarcodeCounter) GetBarcodeWhitelist() []string {
	return c.whitelistArray
}

func (c *BarcodeCounter) GetWhitelistBarcodeCount() []int {
	if !c.whitelistExists {
		return []int{}
	}
	whitelistBarcodeCounts := make([]int, len(c.whitelistArray))
	for idx, barcode := range c.whitelistArray {
		whitelistBarcodeCounts[idx] = c.barcodeCounts[barcode]
	}
	return whitelistBarcodeCounts
}

func (c *BarcodeCounter) GetMismatchCount() int {
	return c.mismatchCount
}

func (c *BarcodeCounter) GetTotalCount() int {
	return c.totalCount
}

func (c *BarcodeCounter) GetTotalBaseCount() int64 {
	return c.totalBaseCount
}

func (c *BarcodeCounter) GetQ30BaseCount() int64 {
	return c.q30BaseCount
}

func (c *BarcodeCounter) GetMeanBaseQualityScore() float64 {
	if c.totalBaseCount <= 0 {
		return 0
	} else {
		return float64(c.baseQualSum) / float64(c.totalBaseCount)
	}
}

/**
 * Generate a BarcodeValidator object using the existing barcode count and
 * whitelist as the basis for barcode correction.
 */
func (c *BarcodeCounter) GetBarcodeValidator(max_expected_barcode_errors float64, bcConfidenceThreshold float64) *BarcodeValidator {
	// TODO: may just want to use the detected barcodes in the future for added
	// flexibility but for now make sure there's a vetted whitelist to ensure
	// behavioral parity
	if !c.whitelistExists {
		return &BarcodeValidator{whitelist: make(map[string]bool), max_expected_errors: max_expected_barcode_errors, whitelist_exists: false}
	}
	// use same initialization algorithm as NewBarcodeValidator
	barcodeSum := 0.0
	barcodeRates := make(map[string]float64)
	barcodeSum += float64(c.totalCount + 1)
	sortedBarcodes := c.GetBarcodeWhitelist()
	sortedCounts := c.GetWhitelistBarcodeCount()
	for index, count := range sortedCounts {
		barcodeRates[sortedBarcodes[index]] = float64(count+1) / barcodeSum
	}
	whitelistMap := make(map[string]bool)
	for _, barcode := range sortedBarcodes {
		whitelistMap[barcode] = true
	}
	return &BarcodeValidator{whitelist: whitelistMap, max_expected_errors: max_expected_barcode_errors,
		barcodeRates: barcodeRates, bcConfidenceThreshold: bcConfidenceThreshold, whitelist_exists: true}
}

type BarcodeValidator struct {
	whitelist_exists      bool
	whitelist             map[string]bool
	max_expected_errors   float64
	barcodeRates          map[string]float64
	bcConfidenceThreshold float64
}

func loadWhitelist(whitelist_filename string) (map[string]bool, []string) {
	file, _ := os.Open(whitelist_filename)
	defer file.Close()
	reader := bufio.NewReader(file)
	whitelist_map := map[string]bool{}
	whitelistArray := []string{}
	for {
		line, err := reader.ReadString('\n') // will contain "\n" but so should the barcodes we are asked to validate
		if err != nil {
			break
		}

		if !strings.HasPrefix(line, "#") {
			whitelist_map[line] = true
			whitelistArray = append(whitelistArray, line)
		}
	}
	sort.Strings(whitelistArray)
	return whitelist_map, whitelistArray
}

func getCountsForGem(barcodeCounts string, gemGroup string) []int {
	//load barcode counts
	bcCountFile, _ := os.Open(barcodeCounts)
	defer bcCountFile.Close()
	decoder := json.NewDecoder(bcCountFile)
	var barcodeCountsObj map[string]BarcodeCounts
	decoder.Decode(&barcodeCountsObj)
	return barcodeCountsObj[gemGroup].Bc_counts
}

func NewBarcodeValidator(whitelist_filename string, max_expected_barcode_errors float64, barcodeCounts string, bcConfidenceThreshold float64, gemGroup string) *BarcodeValidator {
	if whitelist_filename == "none" || barcodeCounts == "none" {
		return &BarcodeValidator{
			whitelist:           map[string]bool{},
			max_expected_errors: max_expected_barcode_errors,
			whitelist_exists:    false,
		}
	}
	whitelist_map, whitelistArray := loadWhitelist(whitelist_filename)
	bc_counts := getCountsForGem(barcodeCounts, gemGroup)
	if len(bc_counts) > len(whitelistArray) {
		panic(fmt.Sprintf(
			"Expected no more than %d barcode counts for gem group %s in %s, "+
				"but found %d.  There should not be counts for barcodes not in %s",
			len(whitelistArray), gemGroup, barcodeCounts,
			len(bc_counts),
			whitelist_filename))
	}
	barcodeRates := make(map[string]float64, len(bc_counts))
	barcodeSum := 0.0
	for _, count := range bc_counts {
		barcodeSum += float64(count + 1)
	}
	for index, count := range bc_counts {
		barcodeRates[whitelistArray[index]] = float64(count+1) / barcodeSum
	}
	return &BarcodeValidator{
		whitelist:             whitelist_map,
		max_expected_errors:   max_expected_barcode_errors,
		whitelist_exists:      true,
		barcodeRates:          barcodeRates,
		bcConfidenceThreshold: bcConfidenceThreshold,
	}
}

type BarcodeCounts struct {
	Bad_bc_count  int
	Bc_error_rate float64
	Bc_counts     []int
}

func (v *BarcodeValidator) correctBarcode(barcode string, qual []byte) ([]byte, bool) {
	a := make([]byte, len(barcode))
	for i, val := range []byte(barcode) {
		a[i] = val
	}
	whitelistCandidates := []string{}
	likelihoods := []float64{}
	totalLikelihood := 0.0
	vals := []byte{[]byte("A")[0], []byte("C")[0], []byte("G")[0], []byte("T")[0]}
	for pos, qv := range qual {
		existing := a[pos]
		for _, val := range vals {
			if val == existing {
				continue
			}
			a[pos] = val
			barcodeRate, isBarcode := v.barcodeRates[string(a)]
			if isBarcode {
				probEdit := math.Max(0.0005, probability(qv))
				whitelistCandidates = append(whitelistCandidates, string(a))
				likelihood := barcodeRate * probEdit
				totalLikelihood += likelihood
				likelihoods = append(likelihoods, likelihood)
			}
		}
		a[pos] = existing
	}

	maxProbIndex := -1
	maxProb := -math.MaxFloat64
	for index := range likelihoods {
		likelihoods[index] /= totalLikelihood
		if likelihoods[index] > maxProb {
			maxProb = likelihoods[index]
			maxProbIndex = index
		}
	}
	if maxProbIndex != -1 && maxProb > v.bcConfidenceThreshold {
		return []byte(whitelistCandidates[maxProbIndex]), true
	}
	return []byte(barcode), false

}

func (v *BarcodeValidator) NoWhitelist() bool {
	return !v.whitelist_exists
}

func (v *BarcodeValidator) ValidateBarcode(barcode string, qual []byte) ([]byte, bool) {
	if !v.whitelist_exists && len(barcode) > 1 {
		return []byte(barcode), true
	}
	expected_errors := 0.0
	for i := 0; i < len(qual); i++ {
		expected_errors += probability(qual[i])
	}
	_, in_whitelist := v.whitelist[barcode]
	if !in_whitelist {
		newBarcode, corrected := v.correctBarcode(barcode, qual)
		if corrected && expected_errors < v.max_expected_errors {
			return newBarcode, true
		}
	}
	if in_whitelist && expected_errors < v.max_expected_errors {
		return []byte(barcode), true
	}
	return []byte(barcode), false
}

func probability(qual byte) float64 {
	return math.Pow(10, -(float64(uint32(qual)-33.0))/10.0) //33 is the illumina qual offset
}
