#!/usr/bin/env python3

"""
  Adapted from nf-core/rnaseq and nf-core/circrna

  Strandedness commented out should you wish to include this in future releases.
"""

import os
import sys
import errno
import argparse
import csv

def parse_args(args=None):
    Description = "Reformat nf-core/rnaseq samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py --FILE_IN <FILE_IN> --FILE_OUT <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input samplesheet file.")
    parser.add_argument("FILE_OUT", help="Output file.")
    return parser.parse_args(args)


def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception


def sniff_format(handle):
    """
    Detect the tabular format.

    Args:
        handle (text file): A handle to a `text file`_ object. The read position is
        expected to be at the beginning (index 0).

    Returns:
        csv.Dialect: The detected tabular format.

    .. _text file:
        https://docs.python.org/3/glossary.html#term-text-file

    """
    peek = read_head(handle)
    handle.seek(0)
    sniffer = csv.Sniffer()
    dialect = sniffer.sniff(peek)
    return dialect


def check_samplesheet(file_in, file_out):
    """
    This function checks that the samplesheet follows one of the following structures:
    1. For FASTQ format:
       sample,fastq_1,fastq_2
       SAMPLE_PE,SAMPLE_PE_RUN1_1.fastq.gz,SAMPLE_PE_RUN1_2.fastq.gz
       SAMPLE_PE,SAMPLE_PE_RUN2_1.fastq.gz,SAMPLE_PE_RUN2_2.fastq.gz
       SAMPLE_SE,SAMPLE_SE_RUN1_1.fastq.gz,
    2. For BAM format:
       sample,bam
       SAMPLE_PE,SAMPLE_PE_RUN1_1.bam
       SAMPLE_PE,SAMPLE_PE_RUN2_1.bam
       SAMPLE_SE,SAMPLE_SE_RUN1_1.bam
    """
    sample_mapping_dict = {}
    with open(file_in, "r", encoding="utf-8-sig") as fin:
        HEADER_BAM = ["sample", "bam"]
        HEADER_FASTQ = ["sample", "fastq_1", "fastq_2"]
        header = [x.strip('"') for x in fin.readline().strip().split(",")]
        
        if header == HEADER_BAM:
            format_type = "BAM"
        elif header == HEADER_FASTQ:
            format_type = "FASTQ"
        else:
            print("ERROR: Invalid samplesheet header. It should be either 'sample,bam' or 'sample,fastq_1,fastq_2'.")
            sys.exit(1)

        for line in fin:
            if line.strip():
                lspl = [x.strip().strip('"') for x in line.strip().split(",")]

                # Check valid number of columns per row
                if format_type == "FASTQ":
                    MIN_COLS = 3  # Minimum columns for FASTQ format
                else:
                    MIN_COLS = 2  # Minimum columns for BAM format
                    
                num_cols = len([x for x in lspl if x])
                if num_cols < MIN_COLS:
                    print_error(f"Invalid number of columns (minimum = {MIN_COLS})!", "Line", line)

                if format_type == "FASTQ":
                    sample, fastq_1, fastq_2 = lspl[: len(header)]
                else:
                    sample, bam = lspl[: len(header)]

                # Check sample name entries
                if sample.find(" ") != -1:
                    print(f"WARNING: Spaces have been replaced by underscores for sample: {sample}")
                    sample = sample.replace(" ", "_")
                if not sample:
                    print_error("Sample entry has not been specified!", "Line", line)

                # Check BAM file extension
                if format_type == "BAM":
                    if bam:
                        if bam.find(" ") != -1:
                            print_error("BAM file contains spaces!", "Line", line)
                        if not bam.endswith(".bam"):
                            print_error("BAM file does not have extension '.bam'!", "Line", line)

                sample_info = []  # [format_type, bam/fastq_1, fastq_2]
                if format_type == "BAM":
                    if sample and bam:
                        sample_info = [format_type, bam]
                    else:
                        print_error("Invalid combination of columns provided!", "Line", line)
                else:
                    if sample and fastq_1:
                        sample_info = [format_type, fastq_1]
                        if fastq_2:
                            sample_info.append(fastq_2)
                    else:
                        print_error("Invalid combination of columns provided!", "Line", line)

                # Create sample mapping dictionary
                if sample not in sample_mapping_dict:
                    sample_mapping_dict[sample] = [sample_info]
                else:
                    if sample_info in sample_mapping_dict[sample]:
                        print_error("Samplesheet contains duplicate rows!", "Line", line)
                    else:
                        sample_mapping_dict[sample].append(sample_info)

    # Write validated samplesheet with appropriate columns
    if len(sample_mapping_dict) > 0:
        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:
            if format_type == "FASTQ":
                fout.write(",".join(["sample", "format_type", "fastq_1", "fastq_2"]) + "\n")
            else:
                fout.write(",".join(["sample", "format_type", "bam"]) + "\n")
            for sample in sorted(sample_mapping_dict.keys()):
                if not all(x[0] == sample_mapping_dict[sample][0][0] for x in sample_mapping_dict[sample]):
                    print_error(f"Multiple runs of a sample must be of the same datatype!", "Sample", sample)

                for idx, val in enumerate(sample_mapping_dict[sample]):
                    if format_type == "FASTQ":
                        fout.write(",".join([f"{sample}_T{idx+1}"] + val) + "\n")
                    else:
                        fout.write(",".join([f"{sample}_T{idx+1}", val[0], val[1]]) + "\n")
    else:
        print_error(f"No entries to process!", f"Samplesheet: {file_in}")


def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN, args.FILE_OUT)


if __name__ == "__main__":
    sys.exit(main())
