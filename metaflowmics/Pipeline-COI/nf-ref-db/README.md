# Download from BOLD (http queries)
Download all metazoan classes in parallel:

```bash
$ head ~/db/BOLD/taxa_to_dl
Diptera
Lepidoptera
Chordata
Coleoptera
```

Manually downloaded:
```
Annelida
Arachnida
Chordata
Coleoptera
Collembola
Diptera
Hemiptera
Hymenoptera
Lepidoptera
Malacostraca
Mollusca
```

# Preprocessing

  --> 7,342,377 entries
- Remove entries with "SUPPRESSED"
  --> 7,200,558
- Remove entries with no annotation at order level
  --> 7,079,725
- Remove entries with ambiguous nucleotides
  --> 6,533,927
- Remove entries shorter than 600bp
  --> 3,761,906
- If two entries have the same sequence, keep the one with the most annotated taxonomic levels
  --> 2,114,394 sequences
- Discard sequences with inner stop codons in the invertebrate mitochondrial code
  --> 2,073,540 sequences
- Cluster with CD-HIT at 98% similarity
  --> 191,049 sequences 
