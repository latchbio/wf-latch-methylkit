import json
import os
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Annotated, Iterable, List, Optional, Tuple, Union

import pandas as pd
from latch.resources.launch_plan import LaunchPlan
from latch.resources.tasks import medium_task, small_task
from latch.resources.workflow import workflow
from latch.types import LatchDir, LatchFile
from latch.types.directory import LatchOutputDir
from latch.types.file import LatchFile
from latch.types.metadata import LatchAuthor, LatchMetadata, LatchParameter


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


@medium_task
def format_bed_files(
    samples: List[Sample], output_directory: LatchOutputDir, track_name: str
) -> List[ProcessedBED]:

    output_dir = f"{track_name}/processed_beds"
    output_dirpath = Path(output_dir).resolve()
    output_dirpath.mkdir(parents=True, exist_ok=True)
    print("BED File Output Directory: ", output_dirpath)

    bed_results = []

    for s in samples:
        sample_name = s.sample_name
        table = pd.read_csv(s.bed_file.local_path, sep="\t", header=None)
        table[11] = table[10] * 0.01 * table[9]
        table[12] = table[10] * 0.01
        table[11] = table[11].round()
        table[11] = pd.to_numeric(table[11], downcast="integer")
        bed_table = table[[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12]]
        bed_table = bed_table.fillna(0)
        # table[11] = table[10] * 0.01 * table[9]
        # table[11] = table[11].round()
        # table[11] = pd.to_numeric(table[11], downcast="integer")
        # bed_table = table[[0, 1, 2, 5, 9, 11]]
        # bed_table.columns = ["chr", "start", "end", "strand", "coverage", "numC"]
        # bed_table = bed_table.fillna(0)
        s_path = str(output_dirpath) + f"/{sample_name}.bed"
        bed_table.to_csv(
            s_path,
            index=False,
            header=False,
            sep="\t",
        )
        print(s_path)
        output_location = f"{output_directory.remote_directory}/{track_name}/processed_beds/{s.sample_name}.bed"
        print(output_location)

        bed_results.append(
            ProcessedBED(
                sample_name=s.sample_name,
                bed_file=LatchFile(
                    s_path,
                    output_location,
                ),
                treatment=s.treatment,
            )
        )

    return bed_results


@medium_task
def methyl_task(
    samples: List[ProcessedBED],
    output_directory: LatchOutputDir,
    track_name: str,
    base_cov_val: int = 3,
    tiling_val: int = 200,
    tile_coverage: int = 10,
    difference_val: int = 25,
    q_val: float = 0.01,
) -> LatchDir:

    file_names = []
    file_paths = []
    treatments = []

    for s in samples:
        file_names.append(s.sample_name)
        file_paths.append(s.bed_file.local_path)
        if s.treatment:
            treatments.append(1)
        else:
            treatments.append(0)

    delimiter = ","
    file_names = delimiter.join(file_names)
    file_paths = delimiter.join(file_paths)
    treatements_str = ",".join(map(str, treatments))
    print("TREATMENTS", treatments, treatements_str)

    output_dir = f"{track_name}/methylkit"
    output_dirpath = Path(output_dir).resolve()
    output_dirpath.mkdir(parents=True, exist_ok=True)
    print("MethylKit Output Directory: ", output_dirpath)

    subprocess.run(
        " ".join(
            [
                "Rscript",
                "methylkit_task.R",
                str(file_names),
                str(file_paths),
                str(output_dirpath),
                str(base_cov_val),
                str(tiling_val),
                str(tile_coverage),
                str(difference_val),
                str(q_val),
                str(treatements_str),
            ]
        ),
        check=True,
        shell=True,
    )

    output_location = (
        f"{output_directory.remote_directory}/{track_name}/methylkit_results"
    )
    return LatchDir(str(output_dirpath), output_location)


def interpolate_color(normalized_val):
    # Red (255,0,0) to Blue (0,0,255)
    red = int((1 - normalized_val) * 255)
    blue = int(normalized_val * 255)
    return f"{red},0,{blue}"


@small_task
def create_track(
    DMR_results: LatchDir,
    output_directory: LatchOutputDir,
    track_name: str,
) -> LatchDir:

    output_dir = f"{track_name}/IGV_tracks"
    output_dirpath = Path(output_dir).resolve()
    output_dirpath.mkdir(parents=True, exist_ok=True)
    print("IGV Track Output Directory: ", output_dirpath)

    file_path = ""
    for file in DMR_results.iterdir():
        current = file.local_path
        if "DMR_regions.csv" in current:
            file_path = current

    table = pd.read_csv(file_path)

    table.columns = ["chr", "start", "end", "strand", "pvalue", "qvalue", "meth.diff"]
    df = table.copy()
    min_pval = df["pvalue"].min()
    max_pval = df["pvalue"].max()

    # Avoid division by zero in case all p-values are the same
    if min_pval == max_pval:
        df["normalized"] = 1
    else:
        df["normalized"] = (df["pvalue"] - min_pval) / (max_pval - min_pval)

    # Function to interpolate between red and blue based on the normalized p-value
    df["color"] = df["normalized"].apply(interpolate_color)
    df["score"] = 0
    df["strand"] = "."
    df["1_o"] = df["start"]
    df["2_o"] = df["end"]
    df["pvalue"] = df["pvalue"].apply(str)
    df["qvalue"] = df["qvalue"].apply(str)
    df["meth.diff"] = df["meth.diff"].apply(str)
    df["name"] = df["qvalue"]
    df = df[["chr", "start", "end", "name", "score", "strand", "1_o", "2_o", "color"]]
    df_path = str(output_dirpath) + f"/IGV_track.bed"
    df.to_csv(df_path, index=False, header=False, sep="\t")
    str_track_name = f'track name="{track_name}" description="." itemRgb="On"'
    with open(df_path, "r") as file:
        original_content = file.read()

    new_content = f"{str_track_name}\n{original_content}"

    with open(df_path, "w") as file:
        file.write(new_content)

    output_location = f"{output_directory.remote_directory}/{track_name}/IGV_tracks"
    return LatchDir(str(output_dirpath), output_location)
