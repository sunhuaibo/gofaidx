# gofaidx

get a subsequence from fasta

```go
idx, err := NewFaidx("fasta")

subseq := idx.Fetch("chr1", 100, 200)

gc := subseq.CalculateGC()

```
