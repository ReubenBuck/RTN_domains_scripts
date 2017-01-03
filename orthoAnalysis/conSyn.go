// simple idea for this one
// read in three files and pull out orthologs with conserved adjacency.

// maybe its safe to assume these gene lists are sorted.
// if so it has made our job a hell of a lot easier.

package main

import (
	"bufio"
	"bytes"
	"flag"
	"fmt"
	"log"
	"os"
	"strconv"
	"strings"
)

type gene struct {
	orthoName string
	ensName   string
	chr       string
	start     int
	end       int
	strand    string
}

type genes []gene

//type chromosomes map[string]genes
//type genome []chromosomes

const (
	ensNameField = 1
	chrField     = 2
	strandField  = 3
	startField   = 4
	endField     = 5

	orthoOneField = 0
	orthoTwoField = 1
)

func main() {

	// start reading in pairs

	orthoPtr := flag.String("orthoPairs", "", "a string")
	refBasePtr := flag.String("refBase", "ref", "basename for reference genome output")
	queBasePtr := flag.String("queBase", "que", "basename for query genome output")
	outDirPtr := flag.String("outDir", "./", "set output directory, defults to current")
	// something to check if outdir exists

	flag.Parse()

	if len(flag.Args()) != 2 {
		fmt.Println("incorrect number of positional args as input gene models")
		os.Exit(1)
	}

	outRef, err := os.Create(*outDirPtr + "/" + *refBasePtr + ".out")
	if err != nil {
		log.Fatalf("creating file: %v", err)
	}
	defer outRef.Close()

	outQue, err := os.Create(*outDirPtr + "/" + *queBasePtr + ".out")
	if err != nil {
		log.Fatalf("creating file: %v", err)
	}
	defer outQue.Close()

	orthoFile, err := os.Open(*orthoPtr)
	if err != nil {
		log.Fatalf("reading file: %v", err)
	}
	defer orthoFile.Close()

	// put pairs in a map and assign each pair a unique orthoID
	orthoPairs := make(map[string]string)
	scOrtho := bufio.NewScanner(orthoFile)
	orthoID := 0
	for scOrtho.Scan() {
		f := bytes.Fields(scOrtho.Bytes())
		orthoPairs[string(f[orthoOneField])] = "ortho" + strconv.Itoa(orthoID)
		orthoPairs[string(f[orthoTwoField])] = "ortho" + strconv.Itoa(orthoID)
		orthoID++
	}

	// read through both gene files and assign each gene its ortho ID
	// could also include a sort
	// just figure out a way to sort on chr

	var genomes []genes
	for i := 0; i < 2; i++ {

		var genome genes

		geneFile, err := os.Open(flag.Args()[i])
		if err != nil {
			log.Fatalf("reading file: %v", err)
		}
		defer geneFile.Close()

		scGene := bufio.NewScanner(geneFile)
		for scGene.Scan() {
			f := bytes.Fields(scGene.Bytes())

			if orthoPairs[string(f[ensNameField])] != "" {

				g := gene{
					orthoName: orthoPairs[string(f[ensNameField])],
					ensName:   string(f[ensNameField]),
					chr:       string(f[chrField]),
					strand:    string(f[strandField]),
				}
				g.start, err = strconv.Atoi(string(f[startField]))
				g.end, err = strconv.Atoi(string(f[endField]))
				if err != nil {
					log.Fatalf("trouble parsing bytes to ints", err)
				}
				genome = append(genome, g)
			}

		}
		genomes = append(genomes, genome)
	}

	// build a map with the ortho IDs in order to index the data
	genomeIndex := make(map[string]int)
	for i := 0; i < len(genomes[1]); i++ {
		genomeIndex[genomes[1][i].orthoName] = i
	}

	// set up each writer
	wRef := bufio.NewWriter(outRef)
	wQue := bufio.NewWriter(outQue)

	for i := 1; i < len(genomes[0]); i++ {
		genomeI := genomeIndex[genomes[0][i].orthoName]
		//	fmt.Println(genomes[0][i], genomes[1][genomeI], i, genomeI)

		// ensure that indexing is inbounds
		var upper string
		var lower string
		if len(genomes[1]) > genomeI+1 {
			upper = genomes[1][genomeI+1].orthoName
		} else {
			upper = genomes[1][genomeI-1].orthoName
		}

		if genomeI-1 >= 0 {
			lower = genomes[1][genomeI-1].orthoName
		} else {
			lower = genomes[1][genomeI+1].orthoName
		}

		refOutStr := []string{
			genomes[0][i-1].orthoName + "_" + genomes[0][i].orthoName,
			genomes[0][i-1].chr,
			"",
			"",
			genomes[0][i-1].strand,
			genomes[0][i-1].orthoName,
			genomes[0][i-1].ensName,
			genomes[0][i].chr,
			"",
			"",
			genomes[0][i].strand,
			genomes[0][i].orthoName,
			genomes[0][i].ensName + "\n",
		}
		refOutStr[2] = strconv.Itoa(genomes[0][i-1].start)
		refOutStr[3] = strconv.Itoa(genomes[0][i-1].end)
		refOutStr[8] = strconv.Itoa(genomes[0][i].start)
		refOutStr[9] = strconv.Itoa(genomes[0][i].end)

		if genomes[0][i-1].orthoName == lower &&
			genomes[0][i].chr == genomes[0][i-1].chr &&
			genomes[1][genomeI].chr == genomes[1][genomeI-1].chr {
			// create the string then send it to the writer
			wRef.WriteString(strings.Join(refOutStr, "\t"))

			queOutStr := []string{
				genomes[1][genomeI-1].orthoName + "_" + genomes[1][genomeI].orthoName,
				genomes[1][genomeI-1].chr,
				"",
				"",
				genomes[1][genomeI-1].strand,
				genomes[1][genomeI-1].orthoName,
				genomes[1][genomeI-1].ensName,
				genomes[1][genomeI].chr,
				"",
				"",
				genomes[1][genomeI].strand,
				genomes[1][genomeI].orthoName,
				genomes[1][genomeI].ensName + "\n",
			}
			queOutStr[2] = strconv.Itoa(genomes[1][genomeI-1].start)
			queOutStr[3] = strconv.Itoa(genomes[1][genomeI-1].end)
			queOutStr[8] = strconv.Itoa(genomes[1][genomeI].start)
			queOutStr[9] = strconv.Itoa(genomes[1][genomeI].end)

			wQue.WriteString(strings.Join(queOutStr, "\t"))

		} else if genomes[0][i-1].orthoName == upper &&
			genomes[0][i].chr == genomes[0][i-1].chr &&
			genomes[1][genomeI].chr == genomes[1][genomeI+1].chr {
			wRef.WriteString(strings.Join(refOutStr, "\t"))

			queOutStr := []string{
				genomes[1][genomeI+1].orthoName + "_" + genomes[1][genomeI].orthoName,
				genomes[1][genomeI+1].chr,
				"",
				"",
				genomes[1][genomeI+1].strand,
				genomes[1][genomeI+1].orthoName,
				genomes[1][genomeI+1].ensName,
				genomes[1][genomeI].chr,
				"",
				"",
				genomes[1][genomeI].strand,
				genomes[1][genomeI].orthoName,
				genomes[1][genomeI].ensName + "\n",
			}
			queOutStr[2] = strconv.Itoa(genomes[1][genomeI+1].start)
			queOutStr[3] = strconv.Itoa(genomes[1][genomeI+1].end)
			queOutStr[8] = strconv.Itoa(genomes[1][genomeI].start)
			queOutStr[9] = strconv.Itoa(genomes[1][genomeI].end)

			wQue.WriteString(strings.Join(queOutStr, "\t"))
		}
		//fmt.Println(i)

	}

	wRef.Flush()
	wQue.Flush()
	// make sure chromosome matches
	// After two files are read in and sorted
	// go through and identify conserved synteny
	// chr has ro match previous entry within species and ortho names have to match across species.

}
