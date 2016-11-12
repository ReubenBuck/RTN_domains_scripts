package main

// the basic idea is to merge the ungapped alignmnets over the small gaps
// this is done based on what we consider a reasonably large gap
// for starters we will just choose a minimum size that seems reasonable
// afterwards we will chose a minnimum size based on what seems reasonable
// given the sizes of the ungapped alignmnets

import (
	"bufio"
	"fmt"
	"log"
	"os"
	"strconv"
	"strings"
)

type alignRange struct {
	refAlignLen int
	queAlignLen int
	refGapLen   int
	queGapLen   int
}

type alignRangeMerged struct {
	refAlignLen     int
	queAlignLen     int
	refGapLen       int
	queGapLen       int
	refMergedGapLen int
	queMergedGapLen int
}

type alignRanges []alignRange

type chainInfo struct {
	score     int
	refName   string
	refSize   int
	refStrand string
	refStart  int
	refEnd    int
	queName   string
	queSize   int
	queStrand string
	queStart  int
	queEnd    int
	id        int
}

type alignChain struct {
	header     chainInfo
	alignments alignRanges
}

type alignChains []alignChain

type rangeObject struct {
	refName   string
	refSize   int
	refStrand string
	refStart  int
	refEnd    int
	queName   string
	queSize   int
	queStrand string
	queStart  int
	queEnd    int
	id        int
	refGap    int
	queGap    int
}

type rangeObjects []rangeObject

var aligns alignRanges
var chains alignChains

func main() {

	f, err := os.Open("/Users/labadmin/Desktop/RTN_domains/data/chainAlignments/test/testAlignGo.chain")
	if err != nil {
		log.Fatalf("reading file: %v", err)
	}
	defer f.Close()

	var chain alignChain
	var align alignRange

	sc := bufio.NewScanner(f)
	sc.Split(bufio.ScanLines)
	for sc.Scan() {

		line := sc.Text()
		if line == "" {
			continue
		}
		words := strings.Split(line, " ")
		if words[0] == "chain" {
			chain.header.score, err = strconv.Atoi(words[1])
			chain.header.refName = words[2]
			chain.header.refSize, err = strconv.Atoi(words[3])
			chain.header.refStrand = words[4]
			chain.header.refStart, err = strconv.Atoi(words[5])
			chain.header.refEnd, err = strconv.Atoi(words[6])
			chain.header.queName = words[7]
			chain.header.queSize, err = strconv.Atoi(words[8])
			chain.header.queStrand = words[9]
			chain.header.queStart, err = strconv.Atoi(words[10])
			chain.header.queEnd, err = strconv.Atoi(words[11])
			chain.header.id, err = strconv.Atoi(words[12])
		} else {
			wordsAli := strings.Split(words[0], "\t")
			align.refAlignLen, err = strconv.Atoi(wordsAli[0])
			align.queAlignLen, err = strconv.Atoi(wordsAli[0])
			if len(wordsAli) == 3 {
				align.refGapLen, err = strconv.Atoi(wordsAli[1])
				align.queGapLen, err = strconv.Atoi(wordsAli[2])
				aligns = append(aligns, align)
			} else if len(wordsAli) == 1 {
				align.refGapLen = 0
				align.queGapLen = 0
				aligns = append(aligns, align)
				chain.alignments = aligns
				chains = append(chains, chain)
				aligns = nil
			}

		}
	}
	// for i := 0; i < len(chains); i++ {
	// 	fmt.Println(chains[i].header)
	// 	fmt.Println(chains[i].alignments)
	// }

	var rangeGenome rangeObject
	var rangeGenomes rangeObjects
	minAlign := 10

	for i := 0; i < len(chains); i++ {
		var newAlign alignRangeMerged
		var refGap int
		var queGap int
		fmt.Println(chains[i].header)
		firstTime := true
		for j := 0; j < len(chains[i].alignments); j++ {

			newAlign.refAlignLen = newAlign.refAlignLen + chains[i].alignments[j].refAlignLen
			newAlign.queAlignLen = newAlign.queAlignLen + chains[i].alignments[j].queAlignLen

			if chains[i].alignments[j].queGapLen > minAlign || chains[i].alignments[j].refGapLen > minAlign {
				newAlign.refGapLen = chains[i].alignments[j].refGapLen
				newAlign.queGapLen = chains[i].alignments[j].queGapLen

				// there will probably need to be an if conditional here
				rangeGenome.refName = chains[i].header.refName
				rangeGenome.queName = chains[i].header.queName
				rangeGenome.refStrand = chains[i].header.refStrand
				rangeGenome.queStrand = chains[i].header.queStrand
				rangeGenome.refSize = chains[i].header.refSize
				rangeGenome.queSize = chains[i].header.queSize
				rangeGenome.id = chains[i].header.id
				rangeGenome.refGap = newAlign.refMergedGapLen
				rangeGenome.queGap = newAlign.queMergedGapLen

				fmt.Println(newAlign)

				if firstTime {
					rangeGenome.refStart = chains[i].header.refStart
					rangeGenome.refEnd = rangeGenome.refStart + newAlign.refAlignLen

					rangeGenome.queStart = chains[i].header.queStart
					rangeGenome.queEnd = rangeGenome.queStart + newAlign.queAlignLen

					rangeGenomes = append(rangeGenomes, rangeGenome)
					refGap = newAlign.refGapLen
					queGap = newAlign.queGapLen

					firstTime = false

				} else {
					rangeGenome.refStart = refGap + rangeGenomes[len(rangeGenomes)-1].refEnd
					rangeGenome.refEnd = newAlign.refAlignLen + rangeGenome.refStart

					rangeGenome.queStart = queGap + rangeGenomes[len(rangeGenomes)-1].queEnd
					rangeGenome.queEnd = newAlign.queAlignLen + rangeGenome.queStart

					rangeGenomes = append(rangeGenomes, rangeGenome)

					refGap = newAlign.refGapLen
					queGap = newAlign.queGapLen
				}

				newAlign.refAlignLen = 0
				newAlign.queAlignLen = 0
				newAlign.refGapLen = 0
				newAlign.queGapLen = 0
				newAlign.refMergedGapLen = 0
				newAlign.queMergedGapLen = 0
				continue

				// Here we can push out ranged data
			}

			newAlign.refAlignLen = newAlign.refAlignLen + chains[i].alignments[j].refGapLen
			newAlign.queAlignLen = newAlign.queAlignLen + chains[i].alignments[j].queGapLen
			newAlign.refMergedGapLen = newAlign.refMergedGapLen + chains[i].alignments[j].refGapLen
			newAlign.queMergedGapLen = newAlign.queMergedGapLen + chains[i].alignments[j].queGapLen

			// just get this summing working

		}

		rangeGenome.refName = chains[i].header.refName
		rangeGenome.queName = chains[i].header.queName
		rangeGenome.refStrand = chains[i].header.refStrand
		rangeGenome.queStrand = chains[i].header.queStrand
		rangeGenome.refSize = chains[i].header.refSize
		rangeGenome.queSize = chains[i].header.queSize
		rangeGenome.id = chains[i].header.id
		rangeGenome.refGap = newAlign.refMergedGapLen
		rangeGenome.queGap = newAlign.queMergedGapLen

		if firstTime {
			rangeGenome.refStart = chains[i].header.refStart
			rangeGenome.refEnd = rangeGenome.refStart + newAlign.refAlignLen

			rangeGenome.queStart = chains[i].header.queStart
			rangeGenome.queEnd = rangeGenome.queStart + newAlign.queAlignLen

			rangeGenomes = append(rangeGenomes, rangeGenome)
			refGap = newAlign.refGapLen
			queGap = newAlign.queGapLen

			firstTime = false

		} else {
			rangeGenome.refStart = refGap + rangeGenomes[len(rangeGenomes)-1].refEnd
			rangeGenome.refEnd = newAlign.refAlignLen + rangeGenome.refStart

			rangeGenome.queStart = queGap + rangeGenomes[len(rangeGenomes)-1].queEnd
			rangeGenome.queEnd = newAlign.queAlignLen + rangeGenome.queStart

			rangeGenomes = append(rangeGenomes, rangeGenome)

			refGap = newAlign.refGapLen
			queGap = newAlign.queGapLen
		}

		fmt.Println(newAlign)
		// here we also need to push out ranged data
	}

	for i := 0; i < len(rangeGenomes); i++ {
		fmt.Println(rangeGenomes[i])

	}

	// will read out as two seperate ranges on the same line

	// data is read in and intervals are merging
	// do we want to read them out as their merged version or would it be better to do it as two seperate ranges
	// it probably will also be worth while keeping a tally on the number of gaps that have been merged over

	// probably worth printing it out as something that is easy to read into R as well
	// ie a consitent tabular format

	// running it at zero will just reformat the data, cool
	// what is the best thing to produce from this?
	// maybe start by producing a dataset that contains everything that may be needed.
	// we can use the bufio package to just scan the 100 best chains

}
