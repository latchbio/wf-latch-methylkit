"""
Minimal template workflow to show how to run R programs in Latch

For a more comprehensive template, see the assemble_and_sort workflow
For examples on how to use R in Latch, see https://docs.latch.bio/examples/workflows_examples.html
"""

from dataclasses import dataclass
from enum import Enum
from typing import Annotated, Iterable, List, Optional, Tuple, Union

from latch.resources.launch_plan import LaunchPlan
from latch.resources.tasks import small_task
from latch.resources.workflow import workflow
from latch.types import LatchDir, LatchFile
from latch.types.directory import LatchOutputDir
from latch.types.file import LatchFile
from latch.types.metadata import (
    Fork,
    ForkBranch,
    LatchAuthor,
    LatchMetadata,
    LatchParameter,
    Params,
    Section,
    Spoiler,
    Text,
)

from wf.task import create_track, format_bed_files, methyl_task


@dataclass
class Sample:
    sample_name: str
    bed_file: LatchFile
    treatment: bool = False


@dataclass
class ProcessedBED:
    sample_name: str
    bed_file: LatchFile
    treatment: bool = False


class FileFormat(Enum):
    bismark_cov = "Bismark Coverage File"
    bismark_cytosine = "Bismark Cytosine Report"
    bedmethyl = "bedMethyl Format"


class Genome(Enum):
    hg38 = "hg38"
    hg19 = "hg19"


flow = [
    Section(
        "Sample Information",
        Text(
            "Add BED files per for each sample. Use the treatment flag to indicate which samples should be compared against each other."
        ),
        Params(
            "samples",
            "bed_format",
            "genome",
            "track_name",
        ),
    ),
    Section(
        "Parameter Settings",
        Text(
            "Adjust the parameters below to adjust how Differentially Methylated Regions (DMRs) are classified"
        ),
        Params(
            "base_cov_val",
            "tiling_val",
            "tile_coverage",
            "difference_val",
            "q_val",
        ),
    ),
    Section(
        "Output",
        Params(
            "output_directory",
        ),
    ),
]


"""Minimal metadata object - fill in fields with your own values"""
metadata = LatchMetadata(
    display_name="MethylKit DMR Annotator",
    author=LatchAuthor(
        name="LatchBio",
    ),
    repository="https://github.com/latchbio/wf-latch-methylkit",
    parameters={
        "samples": LatchParameter(
            display_name="Input Samples",
            batch_table_column=True,  # Show this parameter in batched mode.
            samplesheet=True,
        ),
        "genome": LatchParameter(
            display_name="Select Genome",
            batch_table_column=True,  # Show this parameter in batched mode.
        ),
        "bed_format": LatchParameter(
            display_name="Select input format",
            batch_table_column=True,  # Show this parameter in batched mode.
        ),
        "track_name": LatchParameter(
            display_name="Track Name",
            batch_table_column=True,  # Show this parameter in batched mode.
        ),
        "output_directory": LatchParameter(
            display_name="Output Directory",
            batch_table_column=True,  # Show this parameter in batched mode.
        ),
        "base_cov_val": LatchParameter(
            display_name="Minimum Base Coverage",
            batch_table_column=True,  # Show this parameter in batched mode.
            description="Filters a methylRawList and discards bases that have coverage below this \
            threshold (e.g. if set to 10, CpGs with coverage less than 10X will be discarded)",
        ),
        "tiling_val": LatchParameter(
            display_name="Tiling Parameter",
            batch_table_column=True,  # Show this parameter in batched mode.
            description="Used in TileMethylCounts(). This value tiles the genome with \
            windows of this specified length and step-size and summarizes the methylation information on those tiles.",
        ),
        "tile_coverage": LatchParameter(
            display_name="Tiling Coverage Cutoff",
            batch_table_column=True,  # Show this parameter in batched mode.
            description="Filter based on the number of bases (cytosines) per tiled region.",
        ),
        "difference_val": LatchParameter(
            display_name="Percent Methylation Difference Cutoff",
            batch_table_column=True,  # Show this parameter in batched mode.
            description="Used in getMethylDiff(). This value defines what the percent  \
            methylation difference must be greater than between samples for a region in order to be counted as a DMR.",
        ),
        "q_val": LatchParameter(
            display_name="q-Value Cutoff",
            batch_table_column=True,  # Show this parameter in batched mode.
            description="Used in getMethylDiff(). This value defines what the q-value  \
            cutoff is in order for a region to be counted as a DMR.",
        ),
    },
    tags=[],
    flow=flow,
)


@workflow(metadata)
def dmr_methylkit(
    samples: List[Sample],
    track_name: str,
    output_directory: LatchOutputDir,
    base_cov_val: int = 3,
    tiling_val: int = 200,
    tile_coverage: int = 10,
    difference_val: int = 25,
    q_val: float = 0.05,
    bed_format: FileFormat = FileFormat.bismark_cov,
    genome: Genome = Genome.hg19,
) -> LatchDir:
    """Locating differentially methylated regions with methylKit


    # MethylKit
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

    """

    processed_bed_files = format_bed_files(
        samples=samples,
        output_directory=output_directory,
        track_name=track_name,
        bed_format=bed_format,
    )

    methyl_results = methyl_task(
        samples=processed_bed_files,
        output_directory=output_directory,
        base_cov_val=base_cov_val,
        tiling_val=tiling_val,
        tile_coverage=tile_coverage,
        difference_val=difference_val,
        q_val=q_val,
        track_name=track_name,
        bed_format=bed_format,
        genome=genome,
    )

    return create_track(
        DMR_results=methyl_results,
        track_name=track_name,
        output_directory=output_directory,
    )


LaunchPlan(
    dmr_methylkit,
    "WGBS Test Data from ENCODE - chr1",
    {
        "samples": [
            Sample(
                sample_name="ENCLB587BLQ",
                bed_file=LatchFile(
                    "s3://latch-public/test-data/22353/chr1_entries.bed"
                ),
                treatment=True,
            ),
            Sample(
                sample_name="ENCSR344YUA",
                bed_file=LatchFile(
                    "s3://latch-public/test-data/22353/chr1_entries_2.bed"
                ),
                treatment=False,
            ),
        ],
        "base_cov_val": 3,
        "tiling_val": 200,
        "tile_coverage": 10,
        "difference_val": 25,
        "q_val": 0.05,
    },
)
