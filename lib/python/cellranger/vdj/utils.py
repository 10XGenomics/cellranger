# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.

# Do not add new things to this module.
# Instead, either find or create a module with a name that better describes
# the functionality implemented by the methods or classes you want to add.

import json
import os
import sys

import numpy as np
import pandas as pd
from six import ensure_binary, ensure_str

import cellranger.vdj.constants as vdj_constants
import tenkit.cache as tk_cache


def get_genes_in_pair(gene_pair_str):
    return gene_pair_str.split("_")


def load_contig_summary_table(filename):
    df = pd.read_csv(ensure_str(filename))
    # Index on contig id
    df.set_index("contig_id", drop=False, inplace=True, verify_integrity=True)
    return df


def get_barcode_from_contig_name(contig_name):
    return ensure_binary(contig_name).split(b"_")[0]


def is_contig_filtered(contig_df, contig_id):
    """Return true if a contig passed filter.

    contig_df (pd.DataFrame): Contains contig ids and counts.
    """
    if "pass_filter" in contig_df.columns:
        return contig_df["pass_filter"][contig_id]
    else:
        # This is expected for the new contig summary coming from the ASSEMBLE stage
        sys.stderr.write(
            "Checked if contig passes filter but contig summary did not have the column\n"
        )
        return True


def save_contig_summary_table(data_frame, filename):
    data_frame.to_csv(filename, index=False)


def load_cell_barcodes_json(filename):
    with open(filename) as f:
        return json.load(f)


def write_csv_optional(df, filename):
    """Write a pandas dataframe to CSV if it is not None."""
    if df is not None:
        write_csv(df, filename)


def write_csv(df, filename):
    """Write a pandas dataframe to CSV in a standard way."""
    # Verify that the data do not contain commas
    for col in df.select_dtypes([object]):
        if df[col].str.contains(",").any():
            raise ValueError(f"Failed write to {filename}: Column {col} contains commas")

    df.to_csv(filename, header=True, index=False, sep=",")


def write_csv_row(row, f):
    """Write a standard CSV row to an open file."""
    # Verify that the data do not contain commas
    row = [str(x) for x in row]
    for i, v in enumerate(row):
        if "," in v:
            raise ValueError("Failed write to csv file: Column %d contains commas" % i)
    f.write("{}\n".format(",".join(row)))


def format_clonotype_id(clonotype_index, inferred):
    """Takes a 0-based clonotype index and formats it into a clonotype id string."""
    if clonotype_index is None:
        return None
    prefix = "inferred_clonotype" if inferred else "clonotype"
    return "%s%d" % (prefix, 1 + clonotype_index)


def get_mem_gb_from_annotations_json(filename):
    """Estimate mem request for loading an entire annotations json into memory."""
    return np.ceil(
        vdj_constants.MEM_GB_PER_ANNOTATIONS_JSON_GB * float(os.path.getsize(filename)) / 1e9
    )


def get_json_obj_iter(f):
    """Generator that streams items from a list of dicts [{}, {},...]."""
    brace_count = 0
    x = ""
    in_str = False
    backslash_count = 0  # Track consecutive backslashes encountered preceding this char
    bufsize = 4096

    while True:
        buf = f.read(bufsize)
        if len(buf) == 0:
            return
        for c in buf:
            # Handle escaped quotes, possibly preceded by escaped backslashes
            # "foo\" => stay in string, 'foo"...'
            # "foo\\" => leave string, 'foo\'
            # "foo\\\" => stay in string, 'foo\"...'
            if c == '"' and (backslash_count % 2) == 0:
                in_str = not in_str

            if c == "\\":
                backslash_count += 1
            else:
                backslash_count = 0

            if not in_str:
                if c.isspace():
                    continue
                if brace_count == 0 and c == "{":
                    brace_count = 1
                elif brace_count == 1 and c == "}":
                    brace_count = 0
                elif c in ("{", "}"):
                    brace_count += 1 if c == "{" else -1
            if brace_count != 0 or len(x) > 0:
                x += c

            if brace_count == 0 and len(x) > 0:
                yield json.loads(x)
                x = ""


class JsonDictListWriter:
    """Streams a list of dicts as json to a file."""

    def __init__(self, fh):
        """Fh - file handle."""
        self.wrote_any = False  # Track whether we need to prepend with a comma
        fh.write("[\n")

    def write(self, x, fh):
        """Write a dict."""
        if self.wrote_any:
            fh.write("\n,")
        json.dump(x, fh)
        self.wrote_any = True

    def finish(self, fh):
        """Finish writing all dicts."""
        fh.write("\n]")


class CachedJsonDictListWriters:
    """Stream multiple json DictLists.

    Using a FileHandleCache to limit the number of open file handles.
    This is useful when you're writing to many files from the same process.
    """

    def __init__(self, filenames: list[str]):
        """Filenames: list of filenames to write to."""
        self.filenames = filenames
        self.cache = tk_cache.FileHandleCache(mode="w")
        self.writers = [JsonDictListWriter(self.cache.get(fn)) for fn in filenames]

    def __enter__(self):
        return self

    def __exit__(self, e_type, e_val, e_tb):
        for writer, filename in zip(self.writers, self.filenames):
            writer.finish(self.cache.get(filename))

    def write(self, d: dict, file_idx: int):
        """Write content.

        Args:
            d: data to write.
            file_idx: index of filename/writer in original given list
        """
        self.writers[file_idx].write(d, self.cache.get(self.filenames[file_idx]))
