# MethylKit
_Locating differentially methylated regions with methylKit_


methylKit is an R package for analysis and annotation of DNA methylation information obtained by high-throughput bisulfite sequencing.
The package is designed to deal with sequencing data from RRBS and its variants. But, it can potentially handle whole-genome bisulfite sequencing data if proper input format is provided.

## Input Requirements
### BED File
This workflow currently accepts 3 kinds of input files:
* Bismark Coverage Files
* Bismark Cytosine Reports
* bedMethyl Files

#### Bismark Coverage Files

This file type is an output of the ```bismark2bedGraph``` function of [Bismark](https://felixkrueger.github.io/Bismark/). It contains coverage in the file format below (using 1-based genomic genomic coordinates):

```
<chromosome> <start position> <end position> <methylation percentage> <count methylated> <count unmethylated>
```

Typically, these files end in the suffix ```.bismark.cov```. They are not stranded and more information on how they are used by Bismark can be found on page 41 of this [guide](https://www.bioconductor.org/packages/release/bioc/manuals/methylKit/man/methylKit.pdf).

#### Bismark Cytosine Reports

This file is a [genome-wide cytosine report](https://felixkrueger.github.io/Bismark/bismark/methylation_extraction/) output from Bismark. Starting from the coverage output, the Bismark methylation extractor can optionally also output a genome-wide cytosine methylation report.
The module ```coverage2cytosine``` (part of the Bismark package) may also be run individually. It is also sorted by chromosomal coordinates but also contains the sequence context and is in the following format:

```
<chromosome> <position> <strand> <count methylated> <count unmethylated> <C-context> <trinucleotide context>

```

#### bedMethyl Format
These are BED files in the [bedMethyl format](https://www.encodeproject.org/data-standards/wgbs/) from ENCODE.
The bedMethyl file is a bed9+2 file containing the number of reads and the percent methylation. Each column represents the following:

* Reference chromosome or scaffold
* Start position in chromosome
* End position in chromosome
* Name of item
* Score from 0-1000. Capped number of reads
* Strandedness, plus (+), minus (-), or unknown (.)
* Start of where display should be thick (start codon)
* End of where display should be thick (stop codon)
* Color value (RGB)
* Coverage, or number of reads
* Percentage of reads that show methylation at this position in the genome

### Treatment Information

For each BED file, there is a boolean flag used to indicate which samples should be compared against other samples. Group samples by this toggle.

### Sample Name

For each sample, give it a name (this will be used for graphing purposes).

### Track Name

The Track Name will be used in IGV to label the track of DMR Regions.

## Parameters

The following parameters are used in the workflow to process the BED files and assign DMRs:

* Minimum Base Coverage (default = 3): Filters a methylRawList and discards bases that have coverage below this threshold (e.g. if set to 10, CpGs with coverage less than 10X will be discarded).
* Tiling Parameter (default = 200): Used in TileMethylCounts(). This value tiles the genome with windows of this specified length and step-size and summarizes the methylation information on those tiles.
* Tiling Coverage Cutoff (default = 10): Filter based on the number of bases (cytosines) per tiled region.
* Percent Methylation Difference Cutoff (default = 25): Used in getMethylDiff(). This value defines what the percent methylation difference must be greater than between samples for a region in order to be counted as a DMR.
* q-Value Cutoff (default = 0.05): Used in getMethylDiff(). This value defines what the q-value cutoff is in order for a region to be counted as a DMR.

## Output Files

This workflow outputs three folders:
* processed_beds: Contains BED files that are processed from methylBED format with a new column representing frequency of methylation as a float percent. These are used as input to MethylKit.
* methylkit_results: Contains a list of DMRs and multiple plots on coverage.
* IGV_tracks: Contains the IGV track of DMRs between samples.

## Test Data

The test data includes Chromosome 1 from the following WGBS samples from ENCODE:

**WGBS from skeletal muscle myoblast (ENCLB587BLQ) - Hg38**

[GEO Link](https://0-www-ncbi-nlm-nih-gov.brum.beds.ac.uk/geo/query/acc.cgi?acc=GSM2138752), [ENCODE Link](https://www.encodeproject.org/experiments/ENCSR328TBS/)


**WGBS from body of pancreas (ENCLB777VZI) - Hg38**

[GEO Link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3633686), [ENCODE Link](https://www.encodeproject.org/experiments/ENCSR344YUA/)

## Resources

This workflow was based on the following [guide](https://bioconductor.org/packages/release/bioc/vignettes/methylKit/inst/doc/methylKit.html) by MethylKit

[MethylKit GitHub](https://github.com/al2na/methylKit) |
[MethylKit Vignette](https://www.bioconductor.org/packages/release/bioc/manuals/methylKit/man/methylKit.pdf) |
[Bismark Methylation Extraction](https://felixkrueger.github.io/Bismark/bismark/methylation_extraction/)
