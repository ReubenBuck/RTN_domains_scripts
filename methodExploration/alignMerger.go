package main

// the basic idea is to merge the ungapped alignmnets over the small gaps
// this is done based on what we consider a reasonably large gap
// for starters we will just choose a minimum size that seems reasonable
// afterwards we will chose a minnimum size based on what seems reasonable
// given the sizes of the ungapped alignmnets
// we could add some flage to determine how much of the file is read in
// realisticly we probably only care about the top 1000 chains.

// now we need to put in our command line variables

//os.Args[1] chainFile
//os.Args[2] chainLimit
//os.Args[3] alignmnetGap

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
var chainLimit int

func main() {

	f, err := os.Open(os.Args[1])
	if err != nil {
		log.Fatalf("reading file: %v", err)
	}
	defer f.Close()

	var chain alignChain
	var align alignRange

	sc := bufio.NewScanner(f)
	sc.Split(bufio.ScanLines)
	// could insert a counter here to determne how many chains to read in
	chainLimit, err = strconv.Atoi(os.Args[2])
	if err != nil {
		log.Fatalf("chain limit did not convert: %v", err)
	}
	chainCount := 0
	for sc.Scan() {

		line := sc.Text()
		if line == "" {
			continue
		}
		if strings.Contains(line, "#") {
			continue
		}
		words := strings.Split(line, " ")
		if words[0] == "chain" {
			if chainCount == chainLimit {
				break
			}
			chainCount++
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
		if err != nil {
			log.Fatalf("chain header did not convert: %v", err)
		}
	}

	var rangeGenome rangeObject
	var rangeGenomes rangeObjects
	minAlign, err := strconv.Atoi(os.Args[3])
	if err != nil {
		log.Fatalf("min Align convert fail: %v", err)
	}

	for i := 0; i < len(chains); i++ {
		var newAlign alignRangeMerged
		var refGap int
		var queGap int
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

		// here we also need to push out ranged data
	}

	for i := 0; i < len(rangeGenomes); i++ {
		fmt.Println(
			rangeGenomes[i].refName,
			rangeGenomes[i].refSize,
			rangeGenomes[i].refStart,
			rangeGenomes[i].refEnd,
			rangeGenomes[i].refStrand,
			rangeGenomes[i].refGap,
			rangeGenomes[i].queName,
			rangeGenomes[i].queSize,
			rangeGenomes[i].queStart,
			rangeGenomes[i].queEnd,
			rangeGenomes[i].queStrand,
			rangeGenomes[i].queGap,
			rangeGenomes[i].id,
		)
	}
}
