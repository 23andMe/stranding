stranding [![Build Status](https://travis-ci.org/23andMe/stranding.svg?branch=master)](https://travis-ci.org/23andMe/stranding)
=================
Determines genome stranding for sequences mapped to a human reference assembly 

Motivation
====
When handling DNA sequences from third-party data vendors it is common
to need to determine whether a given sequence lies on the forward or
reverse strand of a human genome reference assembly. This is often crucial to
determining which allele is the reference allele and which is the alternate
allele. 

This software package can:

* Determine whether a DNA sequence is on the forward or reverse strand 
  of a reference assembly. 
* Check for strand-flip errors in public datasets such as HapMap
* Validate genomic coordinates of SNV nucleotide probes that are mapped to a
  reference assembly 
* Determine the reference and alternate alleles of probe sequences on
  Illumina Infinium assays

All use-cases require a DNA sequence and an educated guess at where that
sequence is mapped in either the GRCh37/hg19 or GRCh38/hg38 reference genome
assembly.  


Algorithm 
=====

We use the provided coordinates to fetch the surrounding reference sequence,
computes the reverse complement, and scores pairwise local alignments between
these sequences and the provided query sequence. 

* If the query sequence aligns to the forward reference sequence then we
  determine that the query sequence lies on the forward strand of the reference
  sequence. 
* If the query sequence aligns to the reverse complement of the reference
  sequence then we determine that the query sequence lies on the reverse strand
  of the reference sequence and a strand flip is required to get to forward
  genome stranding. 
* If the query sequence aligns to both the forward and reverse complement of
  the reference sequence then an exception is raised. 
* If the query sequence aligns to neither the forward nor reverse complement of
  the reference sequence then an `Unstrandable` exception is raised. 
 

Installation 
====

```
pip install stranding 
```

[SeqSeek](https://github.com/23andMe/seqseek) is installed as a dependency to pull 
the reference sequences but the sequence files themselves are not automatically 
downloaded by default. Download them with either:

```
download_build_37
```

or 

```
download_build_38
```

depending on which assembly you require. 


Interface 
====
Consider a probe sequence for the famous rs429358 marker in the APOE gene in Illumina 
TOP stranding:

```
apoe_probe = 'TGCACCAGGCGGCCGC[A/G]CACGTCCTCCATGTCCGCG'
apoe_5p = 'TACTGCACCAGGCGGCCGC'
apoe_3p = 'CACGTCCTCCATGTCCGCG'
```

We can visually determine that this probe sequence aligns to the negative strand of the 
reference by examining the tracks on the 
[dbSNP site](https://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=429358).

The interbase (0-based) coordinate for the position being assayed is `45411940`
on chromosome 19 of the build 37 reference assembly. 

```
pos = 45411940
chr = 19
```

We use the SeqSeek constants `BUILD37` and `BUILD38` to designate the assembly.  

```
from seqseek import BUILD37, BUILD38
```

We can determine which strand of the reference assembly these flanking sequences 
lie on as follows:

```
from stranding import GenomeStranding 
GenomeStranding().strand_flanks(apoe_5p, apoe_3p, BUILD37, chr, pos)
>>> -1 
```

-1 indicates that the sequences align to the reverse/minus strand of the reference assembly.

+1 indicates that the sequences align to the forward/plus strand of the reference assembly.  


Alignment tolerance
====

By default `GenomeStranding(...).strand_flanks(...)` allows the query sequences to 
differ slightly from the reference sequence because it is common to observe
slight variations in third-party sequence data. This is especially true when working 
with locus-specific databases that account for specific haplotypes. 

```
GenomeStranding().strand_flanks('AT' + apoe_5p, apoe_3p + 'CG', BUILD37, chr, pos)
>>> -1 
```

The same result is produced with these altered sequences because a suitable alignment 
is still found at the specified locus. 


You can restrict this to find exact alignments only by passing `tolerance=1.0` to the 
constructor:

```
GenomeStranding(tolerance=1.0).strand_flanks('AT' + apoe_5p, apoe_3p + 'CG', BUILD37, chr, pos)
```


Alignment offsets
====

It is also common to find that sequences from third-party data sources don't align 
exactly where they claim. 

Some reasons for this include:

* assaying a nearby locus to detect an insertion or deletion
* off-by-one errors or discrepancies due to converting to/from interbase and one-based 
  coordinate systems

By default this package requires that the query sequence aligns to the exact location
specified within the allowed tolerance. 

```
GenomeStranding().strand_flanks(apoe_5p, apoe_3p, BUILD37, chr, pos - 10, window=10)
```

Large window sizes will take longer to compute than small windows. 
