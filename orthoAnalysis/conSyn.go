// simple idea for this one
// read in three files and pull out orthologs with conserved adjacency.

package main

import (
	"bufio"
	"bytes"
	"fmt"
	"log"
	"os"
	"strconv"
)

type gene struct {
	orthoName string
	ensName   string
	chr       string
	start     int
	end       int
	strand    string
}

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
	orthoFile, err := os.Open(os.Args[1])
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
		orthoPairs[string(f[orthoOneField])] = "ortho_" + strconv.Itoa(orthoID)
		orthoPairs[string(f[orthoTwoField])] = "ortho_" + strconv.Itoa(orthoID)
		orthoID++
	}

	fmt.Println(orthoPairs["ENSMMUT00000078431.1"])
	fmt.Println(orthoPairs["ENSMMUT00000026535.3"])

	// read through both gene files and assign each gene its ortho ID
	// could also include a sort
	// just figure out a way to sort on chr
	for i := 2; i < 4; i++ {
		geneFile := os.Open(os.Args[i])
		if err != nil {
			log.Fatalf("reading file: %v", err)
		}
		defer geneFile.Close()

		scGene := bufio.NewScanner(os.Stdin)
		for scGene {
			f := bytes.Fields(sc.Bytes())

		}
	}

	// After two files are read in and sorted
	// go through and identify conserved synteny
	// chr has ro match previous entry within species and ortho names have to match across species.

}
