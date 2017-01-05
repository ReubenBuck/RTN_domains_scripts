package main

// program is run per species and finds the interaction frequency of our conserved syntenic pairs

// two inuts:
// 1. pairs
// 2. paired reads hicpipe input

import (
	"bufio"
	"flag"
	"fmt"
)

func main() {

	orthoPtr := flag.String("orthoPairsCoordinates", "", "output from conSyn.go")
	interactionPtr := flag.String("pairedReads", "", "paired reads used for hiclib input")
	outPtr := flag.String("out", "./interaction.out", "output file")
	flag.Parse()

	outFile, err := os.Create(*outPtr)
	if err != nil {
		log.Fatalf("creating file: %v", err)
	}
	defer outFile.Close()

	// basicaly I'm going to have to build an interval tree

	// not exactly sure how to build a tree

}
