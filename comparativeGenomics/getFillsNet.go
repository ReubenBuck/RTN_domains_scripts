// extract fills from net files, min ref size is 10
// Reuben
// need to see what is happening with the stacking
package main

import (
	"bufio"
	"fmt"
	"log"
	"os"
	"strconv"
	"strings"
)

type fill struct {
	RefChr   string
	RefStart int
	RefEnd   int
	QueChr   string
	QueStart int
	QueEnd   int
	Strand   string
	//Type     string
}

type gap struct {
	RefChr   string
	RefStart int
	RefEnd   int
}

type net struct {
	Chr string
	Len int
}

func main() {

	var newFill fill
	var newGap gap
	var chrStack []string
	var netElem net
	var fillStartStack []int
	var fillEndStack []int
	var quePosStack []int
	var queEndStack []int
	var strandStack []string
	var fillPrint fill

	fillFile, err := os.Open(os.Args[1])
	if err != nil {
		log.Fatalf("reading file: %v", err)
	}
	defer fillFile.Close()

	scFill := bufio.NewScanner(fillFile)
	scFill.Split(bufio.ScanLines)

	for scFill.Scan() {

		var refLen int
		line := scFill.Text()
		line = strings.TrimSpace(line)
		words := strings.Split(line, " ")
		switch words[0] {
		case "net":

			netElem.Chr = words[1]
			netElem.Len, err = strconv.Atoi(words[2])

		case "gap":

			// get information for next gap
			refLen, err = strconv.Atoi(words[2])
			if refLen <= 9 {
				continue
			}
			newGap.RefStart, err = strconv.Atoi(words[1])
			newGap.RefStart = newGap.RefStart + 1
			newGap.RefEnd = newGap.RefStart + refLen - 1

			//fmt.Println(newGap, "newFILL")

			// if the gap is totaly left of our stack, we should print it
			// after that we find we have the fill thats interupted

			if len(fillStartStack) > 0 && len(fillEndStack) > 0 {
				for i := len(fillStartStack) - 1; i > -1; i-- {
					if fillEndStack[i] < newGap.RefStart {

						// send fill info to print
						fillPrint.RefStart = fillStartStack[len(fillStartStack)-1]
						fillPrint.RefEnd = fillEndStack[len(fillEndStack)-1]
						fillPrint.RefChr = netElem.Chr
						fillPrint.QueChr = chrStack[len(chrStack)-1]
						fillPrint.Strand = strandStack[len(strandStack)-1]
						fillPrint.QueStart = quePosStack[len(quePosStack)-1]
						fillPrint.QueEnd = queEndStack[len(queEndStack)-1]

						//	if fillPrint.RefEnd-fillPrint.RefStart > 8 {
						fmt.Println(strings.Trim(fmt.Sprintf("%v", fillPrint), "{}"))
						//	}

						// remove from stack
						fillStartStack = fillStartStack[0 : len(fillStartStack)-1]
						fillEndStack = fillEndStack[0 : len(fillEndStack)-1]
						chrStack = chrStack[0 : len(chrStack)-1]
						quePosStack = quePosStack[0 : len(quePosStack)-1]
						queEndStack = queEndStack[0 : len(queEndStack)-1]
						strandStack = strandStack[0 : len(strandStack)-1]
					}
				}
			}

			if len(fillEndStack) > 0 && newGap.RefStart < fillEndStack[len(fillEndStack)-1] {

				fillPrint.RefStart = fillStartStack[len(fillStartStack)-1]
				fillPrint.RefEnd = newGap.RefStart - 1 //to get it in position
				fillPrint.RefChr = netElem.Chr
				fillPrint.QueChr = chrStack[len(chrStack)-1]
				fillPrint.Strand = strandStack[len(strandStack)-1]
				fillPrint.QueStart = quePosStack[len(quePosStack)-1]
				fillPrint.QueEnd = queEndStack[len(queEndStack)-1]

				//	if fillPrint.RefEnd-fillPrint.RefStart > 8 {
				fmt.Println(strings.Trim(fmt.Sprintf("%v", fillPrint), "{}"))
				//	}

				fillStartStack = fillStartStack[0 : len(fillStartStack)-1]
				fillStartStack = append(fillStartStack, newGap.RefEnd+1) // to get it in position
			}

			//now we need to think what a gap interruption looks like
			// mape our gap interuptions need to work recursivly
			// this may be the problem

		case "fill":

			// get new fill information
			// feed back on position
			newFill.RefStart, err = strconv.Atoi(words[1])
			newFill.RefStart = newFill.RefStart + 1
			refLen, err = strconv.Atoi(words[2])
			newFill.RefEnd = newFill.RefStart + refLen - 1
			newFill.QueChr = words[3]
			newFill.Strand = words[4]
			newFill.QueStart, err = strconv.Atoi(words[5])
			newFill.QueStart = newFill.QueStart + 1
			refLen, err = strconv.Atoi(words[6])
			newFill.QueEnd = newFill.QueStart + refLen - 1

			//fmt.Println(newFill, "newGAP")

			// here we enter fill this means we got to the end of our last fill non interrupted
			// here we need to think about our position a bit more
			// If we can get to a start point that is beyond our end point we need to print
			if len(fillStartStack) > 0 && len(fillEndStack) > 0 {
				for i := len(fillStartStack) - 1; i > -1; i-- {
					if fillEndStack[i] < newFill.RefStart {

						// send fill info to print
						fillPrint.RefStart = fillStartStack[len(fillStartStack)-1]
						fillPrint.RefEnd = fillEndStack[len(fillEndStack)-1]
						fillPrint.RefChr = netElem.Chr
						fillPrint.QueChr = chrStack[len(chrStack)-1]
						fillPrint.Strand = strandStack[len(strandStack)-1]
						fillPrint.QueStart = quePosStack[len(quePosStack)-1]
						fillPrint.QueEnd = queEndStack[len(queEndStack)-1]

						//				if fillPrint.RefEnd-fillPrint.RefStart > 8 {
						fmt.Println(strings.Trim(fmt.Sprintf("%v", fillPrint), "{}"))
						//				}

						// remove from stack
						fillStartStack = fillStartStack[0 : len(fillStartStack)-1]
						fillEndStack = fillEndStack[0 : len(fillEndStack)-1]
						chrStack = chrStack[0 : len(chrStack)-1]
						quePosStack = quePosStack[0 : len(quePosStack)-1]
						queEndStack = queEndStack[0 : len(queEndStack)-1]
						strandStack = strandStack[0 : len(strandStack)-1]
					}
				}
			}

			// add to stack
			fillStartStack = append(fillStartStack, newFill.RefStart)
			fillEndStack = append(fillEndStack, newFill.RefEnd)
			chrStack = append(chrStack, newFill.QueChr)
			quePosStack = append(quePosStack, newFill.QueStart)
			queEndStack = append(queEndStack, newFill.QueEnd)
			strandStack = append(strandStack, newFill.Strand)

		}
		if err != nil {
			log.Fatalf("value assignment error", err)
		}
		//fmt.Println(fillStartStack, fillEndStack, "fillStacks")

	}

	for i := len(fillStartStack) - 1; i > -1; i-- {
		// if fillEndStack[i]-fillStartStack[i] < 10 {
		// 	fillStartStack = fillStartStack[0 : len(fillStartStack)-1]
		// 	fillEndStack = fillEndStack[0 : len(fillEndStack)-1]

		//	} else

		// send fill info to print
		fillPrint.RefStart = fillStartStack[len(fillStartStack)-1]
		fillPrint.RefEnd = fillEndStack[len(fillEndStack)-1]
		fillPrint.RefChr = netElem.Chr
		fillPrint.QueChr = chrStack[len(chrStack)-1]
		fillPrint.Strand = strandStack[len(strandStack)-1]
		fillPrint.QueStart = quePosStack[len(quePosStack)-1]
		fillPrint.QueEnd = queEndStack[len(queEndStack)-1]

		//	if fillPrint.RefEnd-fillPrint.RefStart > 8 {
		fmt.Println(strings.Trim(fmt.Sprintf("%v", fillPrint), "{}"))
		//	}

		// remove from stack
		fillStartStack = fillStartStack[0 : len(fillStartStack)-1]
		fillEndStack = fillEndStack[0 : len(fillEndStack)-1]
		chrStack = chrStack[0 : len(chrStack)-1]
		quePosStack = quePosStack[0 : len(quePosStack)-1]
		queEndStack = queEndStack[0 : len(queEndStack)-1]
		strandStack = strandStack[0 : len(strandStack)-1]

	}

}
