/*
	Author: Huaibo Sun
	Data:   2022-07-28
	E-mail: huaibo_sun@foxmail.com
	license: MIT license
*/

package gofaidx

import (
	"errors"
	"io"
	"os"
	"strings"

	"github.com/biogo/hts/fai"
)

type Faidx struct {
	*fai.File

	idx fai.Index
	rdr io.Reader
}

type Sequence struct {
	Chrom string
	Start int
	End   int
	Seq   string
	Len   int
}

// NewFaidx returns a new Faidx, reading from the given fasta.
func NewFaidx(fasta string) (*Faidx, error) {
	// read fai
	fid, err := os.Open(fasta + ".fai")
	if err != nil {
		return &Faidx{}, err
	}
	defer fid.Close()

	idx, err := fai.ReadFrom(fid)
	if err != nil {
		return &Faidx{}, err
	}

	// read fasta
	fna, err := os.Open(fasta)
	if err != nil {
		return &Faidx{}, err
	}

	file := fai.NewFile(fna, idx)

	return &Faidx{File: file, idx: idx, rdr: fna}, nil
}

// Fetch takes a region and return a *Sequence.
// 0-base
func (f *Faidx) Fetch(chrom string, start, end int) (*Sequence, error) {

	fetchError := errors.New("error: unknow sequence name or start less than zero or end less than start")

	if start < 0 && end < start {
		return &Sequence{}, fetchError
	}

	_, ok := f.idx[chrom]
	if !ok {
		return &Sequence{}, fetchError
	}
	seq, err := f.SeqRange(chrom, start, end)

	if err != nil {
		return &Sequence{}, err
	}

	buf := make([]byte, end-start)
	_, err = seq.Read(buf)

	if err != nil {
		return &Sequence{}, err
	}

	seqFetch := string(buf)
	seqLength := len(seqFetch)
	sequence := &Sequence{Chrom: chrom, Start: start, End: end, Seq: seqFetch, Len: seqLength}
	return sequence, nil
}

// At takes a single point and returns the single base.
// 0-base
func (f *Faidx) At(chrom string, pos int) (string, error) {
	seq, err := f.Fetch(chrom, pos, pos+1)
	if err != nil {
		return "", err
	}

	return seq.Seq, nil
}

// Close the associated Reader.
func (f *Faidx) Close() {
	f.rdr.(io.Closer).Close()
}

// Calculate GC content
func (s *Sequence) CalculateGC() float32 {
	if s.Len == 0 {
		return 0
	}
	posAs, posTs, posGs, posCs := 0, 0, 0, 0
	for _, base := range s.Seq {
		switch string(base) {
		case "A", "a":
			posAs++
		case "T", "t":
			posTs++
		case "G", "g":
			posGs++
		case "C", "c":
			posCs++
		}
	}

	gcNum := posGs + posCs
	gc := float32(gcNum) / float32(s.Len)
	return gc
}

// The complexity is defined as the percentage of base that is different from its before base (base[i] != base[i-1]). For example:
// a 21-bp sequence, with 3 bases that is different from its next base
// seq = "AAAAAATTTTTCCCCCCCCCC"
// complexity = 3 / (21-1) = 15%
// if len(seq) <= 1, return 0
func (s *Sequence) CalculateComplexity() float32 {
	num := 0

	if len(s.Seq) <= 1 {
		return 0
	}

	for i := 1; i < len(s.Seq); i++ {
		if s.Seq[i] != s.Seq[i-1] {
			num++
		}
	}

	complexity := float32(num) / float32(s.Len-1)
	return complexity
}

// count the number of CpGs from given Sequence
// CGC or GCG counts as 1 CpG
func (s *Sequence) CountCpG() int {
	cpgCnt := 0
	lastIsCpG := false

	for i := 0; i < len(s.Seq)-1; i++ {
		dinucleotides := strings.ToUpper(s.Seq[i : i+2])
		if (dinucleotides == "CG" || dinucleotides == "GC") && (!lastIsCpG) {
			cpgCnt++
			lastIsCpG = true
		} else {
			lastIsCpG = false
		}
	}

	return cpgCnt
}
