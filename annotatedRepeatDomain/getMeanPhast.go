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

// given a set of bed intervals extract mean phastcon scores

// start by ordering bed intervals
// then go through required phastcon files
