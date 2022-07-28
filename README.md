# gofaidx

get a subsequence from fasta

```go
cd /path/to/your/

samtools faidx xx.fasta

idx, err := NewFaidx("/path/to/your/xx.fasta")

subseq := idx.Fetch("chr1", 100, 200)

gc := subseq.CalculateGC()

```
