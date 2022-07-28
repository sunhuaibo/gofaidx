package gofaidx

import (
	"testing"
)

// type Sequence struct {
// 	Chrom string
// 	Start int
// 	End   int
// 	Seq   string
// 	Len   int
// }

var seq = "CGATGATGACGGTCTGATGATGCTGGGCGT"
var seqLen = len(seq)
var sequence = Sequence{Chrom: "chr1", Start: 1, End: 100, Seq: seq, Len: seqLen}

func TestCalculateGC(t *testing.T) {

	gc := sequence.CalculateGC()
	t.Log(gc)
}

func TestCalculateComplexity(t *testing.T) {
	compl := sequence.CalculateComplexity()
	t.Log(compl)

}

func TestCountCpG(t *testing.T) {

	cpg := sequence.CountCpG()

	t.Log(cpg)

}
