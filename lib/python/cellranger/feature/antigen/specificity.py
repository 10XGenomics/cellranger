#!/usr/bin/env python
#
# Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
#
#
"""Generate antigen specificity scores and antigen assignments."""
from __future__ import annotations

import csv
from ast import literal_eval
from collections import Counter, OrderedDict, defaultdict, namedtuple
from enum import Enum
from typing import AnyStr

import numpy as np
from scipy.stats import beta

import cellranger.rna.library as rna_library
import tenkit.stats as tk_stats
from cellranger.report import PercentMetric  # pylint: disable=no-name-in-module

SIGNAL_PRIOR = 1
NOISE_PRIOR = 3

FEATURE_SEPARATOR = "|"
UNASSIGNED = "Unassigned"
BLANK = "Blank"
FUNCTIONAL_NAME = "functional_name"
TARGETING_ANTIGEN = "targeting_antigen"
MHC_ALLELE = "mhc_allele"
BEAM_AB = "beam_ab"
BEAM_T = "beam_t"
NO_ALLELE = "no_allele"
SAMPLE_ID = "sample_id"

CANONICAL_VDJ_GENE_PAIRS = ["TRA_TRB", "IGH_IGK", "IGH_IGL"]

ANTIGEN_SPECIFICITY_CSV_HEADER = [
    "barcode",
    "antigen",
    "antigen_umi",
    "control",
    "control_umi",
    "antigen_specificity_score",
    "mhc_allele",
    "raw_clonotype_id",
    "exact_subclonotype_id",
]

AGGR_ANTIGEN_SPECIFICITY_CSV_HEADER = [
    "barcode",
    "sample_id",
    "antigen",
    "antigen_umi",
    "control",
    "control_umi",
    "antigen_specificity_score",
    "mhc_allele",
    "raw_clonotype_id",
    "exact_subclonotype_id",
]

AssignmentPerCell = namedtuple(
    "AssignmentPerCell",
    [
        "barcode",
        "num_antigens",
        "assignment",
        "functional_name",
        "assignment_and_functional_name",
        "clonotype_id",
        "exact_subclonotype_id",
    ],
)

Concordance = namedtuple(
    "Concordance",
    [
        "clonotype_key",
        "size",
        "canonical_pair",
        "assigned_antigen",
        "num_bcs_with_assigned_antigen",
        "concordance",
    ],
)


def get_collapsed_functional_names(func_names):
    """Collapse functional_names (such as HEL|HEL|GP120|GP120|GP120) into a more concise format (HEL|GP120)."""
    return sorted(list(set(func_names)))


def get_functional_map(fr_obj):
    """Build a map from feature id to functional map associated with that feature id (if available) using feature ref object."""
    functional_map = {}
    if FUNCTIONAL_NAME not in fr_obj.all_tag_keys:
        return functional_map
    for feature_def in fr_obj.feature_defs:
        functional_map[feature_def.id.decode()] = feature_def.tags.get(FUNCTIONAL_NAME, "")
    return functional_map


def build_clonotype_maps(filtered_contig_annotations):
    """Build a map of vdj cell barcode to clonotype_id / exact_subclonotype_id."""
    bc_clonotype_map = {}
    bc_exact_subclonotype_map = {}
    with open(filtered_contig_annotations) as csvfile:
        for row in csv.DictReader(csvfile):
            barcode = row["barcode"]
            raw_clonotype_id = row["raw_clonotype_id"]
            exact_subclonotype_id = row["exact_subclonotype_id"]
            if raw_clonotype_id == "":
                continue
            elif barcode in bc_clonotype_map:
                assert bc_clonotype_map[barcode] == raw_clonotype_id
                assert bc_exact_subclonotype_map[barcode] == exact_subclonotype_id
            else:
                bc_clonotype_map[barcode] = raw_clonotype_id
                bc_exact_subclonotype_map[barcode] = exact_subclonotype_id
    return bc_clonotype_map, bc_exact_subclonotype_map


class ClonalGroupLevel(Enum):
    """Specifies grouping level for clonotypes."""

    # pylint:disable=invalid-name
    CLONOTYPE = "clonotypes"
    EXACT_SUBCLONOTYPE = "exact_subclonotypes"


class CellsCategory(Enum):
    """Specifies category of cells."""

    # pylint:disable=invalid-name
    ALL_CELLS = "cells"
    GEX_ONLY_CELLS = "gex_only_cells"


# pylint: disable=too-few-public-methods
class AssignmentMetadata:
    """Metadata for feature assignments."""

    def __init__(
        self,
        cells_name: CellsCategory,
        num_categ_cells,
        num_all_cells,
        num_cells_blanks,
        num_cells_unassigned,
        num_cells_singlets,
        num_cells_multiplets,
    ):
        self.cells_name = cells_name
        self.num_categ_cells = num_categ_cells
        self.num_all_cells = num_all_cells
        self.num_cells_blanks = num_cells_blanks
        self.frac_cells_without_antigen_umis = tk_stats.robust_divide(
            self.num_cells_blanks, self.num_all_cells
        )
        self.num_cells_unassigned = num_cells_unassigned
        self.frac_cells_unassigned = tk_stats.robust_divide(
            self.num_cells_unassigned, self.num_all_cells
        )
        self.num_cells_singlets = num_cells_singlets
        self.frac_singlets = tk_stats.robust_divide(self.num_cells_singlets, self.num_all_cells)
        self.num_cells_multiplets = num_cells_multiplets
        self.frac_multiplets = tk_stats.robust_divide(self.num_cells_multiplets, self.num_all_cells)


class AntigenAssigner:
    """Assigns antigens to barcodes."""

    def __init__(
        self,
        barcode_data: BarcodeData,
    ):
        self.barcode_data = barcode_data
        self.barcode_data.sort_barcode_map()
        self._assignment_per_cell = None
        self.include_functional_names = len(barcode_data.functional_map.keys()) != 0
        self.aggr_pipeline = barcode_data.gem_well_map is not None

    @property
    def assignment_per_cell(self):
        """Returns a dict indexed by bc containing metadata on features assigned."""
        if self._assignment_per_cell is None:
            result: dict[str, AssignmentPerCell] = OrderedDict()
            for bc_as in self.barcode_data.barcode_map.values():
                assignment = bc_as.get_assigned_antigen(self.barcode_data.antigen_to_control)
                if assignment in [UNASSIGNED, BLANK]:
                    result[bc_as.barcode] = AssignmentPerCell(
                        bc_as.barcode,
                        0,
                        assignment,
                        "None",
                        "None",
                        bc_as.clonotype_id,
                        bc_as.exact_subclonotype_id,
                    )
                else:
                    antigens = assignment.split(FEATURE_SEPARATOR)
                    num_antigens = len(antigens)
                    if self.include_functional_names:
                        func_names = FEATURE_SEPARATOR.join(
                            get_collapsed_functional_names(
                                [self.barcode_data.functional_map.get(a, "None") for a in antigens]
                            )
                        )
                        assignment_and_func_names = FEATURE_SEPARATOR.join(
                            [
                                "-".join(tup)
                                for tup in zip(
                                    [
                                        self.barcode_data.functional_map.get(a, "None")
                                        for a in antigens
                                    ],
                                    antigens,
                                )
                            ]
                        )
                    else:
                        func_names = ""
                        assignment_and_func_names = ""
                    result[bc_as.barcode] = AssignmentPerCell(
                        bc_as.barcode,
                        num_antigens,
                        assignment,
                        func_names,
                        assignment_and_func_names,
                        bc_as.clonotype_id,
                        bc_as.exact_subclonotype_id,
                    )
            self._assignment_per_cell = result
        return self._assignment_per_cell

    def get_assignment_metadata(self, cells_name: CellsCategory):
        """Returns metadata for feature assignments for cells or gex_only cells."""
        # Number of all cells
        num_all_cells = len(self.barcode_data.barcode_map)

        if cells_name == CellsCategory.ALL_CELLS:
            data = self.barcode_data.barcode_map
        elif cells_name == CellsCategory.GEX_ONLY_CELLS:
            data = {
                k: v for k, v in self.barcode_data.barcode_map.items() if v.clonotype_id == "None"
            }
        # Number of cells in category (== num_all_cells for CellsCategory.ALL_CELLS)
        num_categ_cells = len(data)
        num_cells_blanks = len(
            [v for k, v in self.assignment_per_cell.items() if v.assignment == BLANK and k in data]
        )
        num_cells_unassigned = len(
            [
                v
                for k, v in self.assignment_per_cell.items()
                if v.assignment == UNASSIGNED and k in data
            ]
        )
        num_cells_singlets = len(
            [v for k, v in self.assignment_per_cell.items() if v.num_antigens == 1 and k in data]
        )
        num_cells_multiplets = (
            num_categ_cells - num_cells_blanks - num_cells_unassigned - num_cells_singlets
        )
        return AssignmentMetadata(
            cells_name,
            num_categ_cells,
            num_all_cells,
            num_cells_blanks,
            num_cells_unassigned,
            num_cells_singlets,
            num_cells_multiplets,
        )

    def get_antigen_assignment_metrics(
        self,
        cells_name: CellsCategory,
    ) -> dict[str, float]:
        """Compute summary metrics on antigen assignments."""
        report_prefix = rna_library.get_library_type_metric_prefix(rna_library.ANTIGEN_LIBRARY_TYPE)
        assignment_metadata = self.get_assignment_metadata(cells_name)

        name_fr_ce_w_no_features = report_prefix + f"frac_{cells_name.value}_blank_antigen"

        frac_cells_with_no_features = assignment_metadata.frac_cells_without_antigen_umis

        name_fr_ce_no_features = report_prefix + f"frac_{cells_name.value}_unassigned_antigen"

        frac_cells_unassigned_features = assignment_metadata.frac_cells_unassigned

        name_fr_ce_w_sg_features = report_prefix + f"frac_{cells_name.value}_with_single_antigen"
        frac_cells_with_single_features = assignment_metadata.frac_singlets

        name_fr_ce_w_mult_features = (
            report_prefix + f"frac_{cells_name.value}_with_multiple_antigen"
        )
        frac_cells_multiple_features = assignment_metadata.frac_multiplets

        res = {
            name_fr_ce_w_no_features: frac_cells_with_no_features,
            name_fr_ce_no_features: frac_cells_unassigned_features,
            name_fr_ce_w_sg_features: frac_cells_with_single_features,
            name_fr_ce_w_mult_features: frac_cells_multiple_features,
            f"{report_prefix}{cells_name.value}_num_blanks": assignment_metadata.num_cells_blanks,
            f"{report_prefix}{cells_name.value}_num_unassigned": assignment_metadata.num_cells_unassigned,
            f"{report_prefix}{cells_name.value}_num_singlets": assignment_metadata.num_cells_singlets,
            f"{report_prefix}{cells_name.value}_num_multiplets": assignment_metadata.num_cells_multiplets,
        }
        return res

    @classmethod
    def load_from_file(cls, filename: AnyStr):
        """Load the data from an antigen_specificity_scores.csv file.

        Args:
            filename: the file to load from

        Returns:
            a new AntigenAssigner object
        """
        barcode_map = defaultdict(BarcodeAS)
        gem_well_map = defaultdict(str)
        functional_map = defaultdict(str)
        has_functional_map = False
        aggr_mode = False
        if filename is None:
            return cls(BarcodeData(None, None))
        with open(filename) as csvfile:
            csvreader = csv.reader(csvfile)
            header = next(csvreader)
            # Check if we are in an aggr pipeline
            if SAMPLE_ID in header:
                aggr_mode = True
                csv_row = namedtuple("CsvRow", AGGR_ANTIGEN_SPECIFICITY_CSV_HEADER)
            else:
                csv_row = namedtuple("CsvRow", ANTIGEN_SPECIFICITY_CSV_HEADER)
            if FUNCTIONAL_NAME in header:
                has_functional_map = True

            for row in map(csv_row._make, csvreader):
                if row.barcode in barcode_map:
                    barcode_map[row.barcode].update_barcode(
                        control={row.control: int(row.control_umi)},
                        antigens={row.antigen: int(row.antigen_umi)},
                        scores={row.antigen: float(row.antigen_specificity_score)},
                    )
                else:
                    exact_subclonotype_id = (
                        "None" if row.raw_clonotype_id == "None" else row.exact_subclonotype_id
                    )
                    sample_id = row.sample_id if aggr_mode else None
                    if aggr_mode:
                        gem_idx = row.barcode.split("-")[1]
                        if gem_idx not in gem_well_map:
                            gem_well_map[gem_idx] = sample_id
                        else:
                            assert gem_well_map[gem_idx] == sample_id
                    if has_functional_map:
                        if row.antigen not in functional_map:
                            functional_map[row.antigen] = row.functional_name
                        else:
                            assert functional_map[row.antigen] == row.functional_name

                    barcode_map[row.barcode] = BarcodeAS(
                        barcode=row.barcode,
                        clonotype_id=row.raw_clonotype_id,
                        exact_subclonotype_id=exact_subclonotype_id,
                        control={row.control: int(row.control_umi)},
                        antigens={row.antigen: int(row.antigen_umi)},
                        specificity_scores={row.antigen: float(row.antigen_specificity_score)},
                        allele=row.mhc_allele,
                    )
            return cls(
                BarcodeData.from_barcode_map(
                    barcode_map=barcode_map,
                    count_gem_well_map=gem_well_map,
                    functional_map=functional_map,
                )
            )

    def write_antigen_specificity_csv(self, path: AnyStr):
        """Create antigen specificity CSV."""
        header = ANTIGEN_SPECIFICITY_CSV_HEADER
        if self.barcode_data.is_aggr():
            header = AGGR_ANTIGEN_SPECIFICITY_CSV_HEADER
        with open(path, "w") as csv_handle:
            csv_writer = csv.writer(csv_handle, delimiter=",")
            csv_writer.writerow(header)
            for bc_as in self.barcode_data.barcode_map.values():
                antigens_to_umi = bc_as.antigens(self.barcode_data.antigen_to_control)
                bc = [bc_as.barcode] * len(antigens_to_umi)
                antigen = antigens_to_umi.keys()
                ag_umi = antigens_to_umi.values()
                ctrl = self.barcode_data.antigen_to_control.values()
                ctrl_umi = [bc_as.controls[c] for c in ctrl]
                allele = [self.barcode_data.seen_antigens[a] for a in antigens_to_umi]
                score = [
                    round(score, 3)
                    for score in bc_as.specificity_scores(
                        self.barcode_data.antigen_to_control
                    ).values()
                ]
                c_id = [bc_as.clonotype_id] * len(antigens_to_umi)
                esc_id = [bc_as.exact_subclonotype_id] * len(antigens_to_umi)
                sample_id = (
                    [self.barcode_data.gem_well_map[bc_as.gem_id()]] * len(antigens_to_umi)
                    if self.barcode_data.is_aggr()
                    else [None] * len(antigens_to_umi)
                )
                items = [bc]
                if self.barcode_data.is_aggr():
                    # Aggr run
                    items.append(sample_id)
                items.extend([antigen, ag_umi, ctrl, ctrl_umi, score, allele, c_id, esc_id])
                csv_writer.writerows(zip(*items))

    def write_antigen_assignment_csv(self, path: AnyStr):
        """Create antigen assignemnt CSV."""
        with open(path, "w") as csv_handle:
            csv_writer = csv.writer(csv_handle, delimiter=",")

            csv_writer.writerow(AssignmentPerCell._fields)
            csv_writer.writerows(self.assignment_per_cell.values())

    def generate_cells_by_clonotype(self, grouped_by: ClonalGroupLevel, clonotypes_csv: AnyStr):
        """Returns CellsPerClonotype grouped by clonotype_id."""
        per_clonotype_dict = defaultdict(list[BarcodeAS])
        for barcode_as in self.barcode_data.barcode_map.values():
            if grouped_by == ClonalGroupLevel.CLONOTYPE:
                key = barcode_as.clonotype_id
            elif grouped_by == ClonalGroupLevel.EXACT_SUBCLONOTYPE:
                key = (
                    "None"
                    if barcode_as.clonotype_id == "None"
                    else "_".join([barcode_as.clonotype_id, barcode_as.exact_subclonotype_id])
                )
            per_clonotype_dict[key].append(barcode_as)
        return CellsPerClonotype.from_dictionary(per_clonotype_dict, grouped_by, clonotypes_csv)


class BarcodeData:
    """Data structure representing a group of BarcodeAS objects."""

    def __init__(
        self,
        count_gem_well_map,
        functional_map,
    ):
        self.barcode_map: dict[str, BarcodeAS] = defaultdict(BarcodeAS)
        self.functional_map: dict[str, str] = functional_map  # antigen:functinal_name
        if count_gem_well_map is None:
            self.gem_well_map = None
        else:
            self.gem_well_map = defaultdict(str)
            for gem_idx, gem_info in count_gem_well_map.items():
                if isinstance(gem_info, list):
                    self.gem_well_map[gem_idx] = gem_info[0]
                else:
                    self.gem_well_map[gem_idx] = gem_info

        self.control_by_allele: dict[str, str] = defaultdict(str)  # allele:control
        self.seen_antigens: dict[str, str] = defaultdict(str)  # antigen:allele

    def is_aggr(self):
        return self.gem_well_map is not None

    def add_barcode(self, bc_as: BarcodeAS):
        """Add a barcode to barcode map."""
        self.barcode_map[bc_as.barcode] = bc_as
        self.update_controls(bc_as.allele, bc_as.controls)
        self.update_antigens(bc_as.allele, bc_as.antigens(self.antigen_to_control))

    def contains_barcode(self, bc):
        """Check if barcode present in barcode map."""
        return bc in self.barcode_map

    def update_barcode(self, bc, control, antigens, allele):
        """Updates a BarcodeAS instance in barcode map with additional antigen data."""
        assert bc in self.barcode_map
        self.barcode_map[bc].update_barcode(control, antigens)
        self.update_controls(allele, control)
        self.update_antigens(allele, antigens)

    def update_controls(self, allele, control):
        """Update control_by_allele."""
        if allele in self.control_by_allele:
            for cnt in control:
                assert (
                    cnt == self.control_by_allele[allele]
                ), f"{cnt} != {self.control_by_allele[allele]}"
        else:
            self.control_by_allele[allele] = next(iter(control.keys()))

    def update_antigens(self, allele, antigens):
        """Update seen_antigens."""
        for antigen in antigens:
            if antigen in self.seen_antigens:
                assert self.seen_antigens[antigen] == allele
            else:
                self.seen_antigens[antigen] = allele

    @property
    def antigen_to_control(self):
        ags = OrderedDict(sorted(self.seen_antigens.items()))
        return {antigen: self.control_by_allele[allele] for antigen, allele in ags.items()}

    def sort_barcode_map(self):
        """Sort barcodes in the barcode_map of BarcodeData."""
        self.barcode_map = OrderedDict(sorted(self.barcode_map.items()))

    @classmethod
    def from_barcode_map(
        cls,
        barcode_map: dict[str, BarcodeAS],
        count_gem_well_map,
        functional_map,
    ):
        """Create a barcode_data object from a map of barcode:BarcodeAS."""
        barcode_data = cls(count_gem_well_map, functional_map)
        for bc, bc_as in barcode_map.items():
            assert bc_as.barcode == bc
            barcode_data.add_barcode(bc_as)
        return barcode_data


class BarcodeAS:
    """Data structure representing antigen counts and specificity scores for a single barcode."""

    __slots__ = (
        "barcode",
        "clonotype_id",
        "exact_subclonotype_id",
        "controls",
        "allele",
        "_antigens",
        "_specificity_scores",
        "_assignments",
        # "sample_id",
    )

    def __init__(
        self,
        barcode: bytes,
        clonotype_id: str,
        exact_subclonotype_id: str,
        control: dict[str, int],  # key:value = control_id:umi
        antigens: dict[str, int],  # key:value = feature_id:umi
        allele: str,
        specificity_scores: dict[str, float] | None = None,  # key:value = feature_id:score
    ):
        self.barcode = barcode
        self.clonotype_id = clonotype_id
        self.exact_subclonotype_id = exact_subclonotype_id

        # sanity checks
        assert len(control) == 1
        assert set(control.keys()).isdisjoint(set(antigens.keys()))
        if specificity_scores:
            assert antigens.keys() == specificity_scores.keys()
        self.controls = control
        self.allele = allele
        self._antigens = antigens
        self._specificity_scores = specificity_scores
        self._assignments = None

    def gem_id(self):
        return self.barcode.split("-")[1]

    def antigens(self, antigen_to_control):
        """Generates an orderd dict of antigens:umi_counts for all observed antigens."""
        ags = OrderedDict((a, 0) for a in antigen_to_control)
        if list(self._antigens.keys()) != list(ags.keys()):
            for antigen, umi_count in self._antigens.items():
                ags[antigen] = umi_count
            self._antigens = ags
        return self._antigens

    def specificity_scores(self, antigen_to_control):
        """Generates an orderd dict of antigens:scores for all observed antigens."""
        if self._specificity_scores is None:
            self._specificity_scores = self.calculate_antigen_specificity(antigen_to_control)
        elif list(self._specificity_scores.keys()) != list(
            self.antigens(antigen_to_control).keys()
        ):
            scores = OrderedDict((a, 0.0) for a in self.antigens(antigen_to_control))
            for antigen, score in self._specificity_scores.items():
                scores[antigen] = score
            self._specificity_scores = scores
        return self._specificity_scores

    def assignments(self, antigen_to_control, threshold=75):
        """Assigns antigens with specificity score above a threshold to this barcode."""
        if self._assignments is None:
            self._assignments = {
                k: v >= threshold for k, v in self.specificity_scores(antigen_to_control).items()
            }
        return self._assignments

    def calculate_antigen_specificity(self, antigen_to_control):
        """Calculate specificity scores for each antigen from the beta-distribution."""
        noise = [self.controls[c] for c in antigen_to_control.values()]
        signal = self.antigens(antigen_to_control).values()
        scores = [
            (1 - beta.cdf(0.925, S + SIGNAL_PRIOR, N + NOISE_PRIOR)) * 100
            for S, N in zip(signal, noise)
        ]
        return dict(zip(list(self.antigens(antigen_to_control).keys()), scores))

    def get_assigned_antigen(self, antigen_to_control):
        """Returns antigen(s) assigned to this barcode as a string."""
        if not any(self.assignments(antigen_to_control).values()):
            if sum(self.antigens(antigen_to_control).values()) == 0:
                return BLANK
            else:
                return UNASSIGNED
        else:
            return FEATURE_SEPARATOR.join(
                [k for k, v in self.assignments(antigen_to_control).items() if v]
            )

    def update_barcode(self, control, antigens, scores=None):
        """Updates a BarcodeAS instance with additional antigen data."""
        assert len(control) == 1
        assert all(a not in self._antigens.keys() for a in antigens)
        assert set(control.keys()).isdisjoint(set(antigens.keys()))
        if scores:
            assert antigens.keys() == scores.keys()

        self._antigens.update(antigens)
        self.controls.update(control)
        if scores:
            self._specificity_scores.update(scores)


class CellsPerClonotype:
    """Stores a list of BarcodeAS objects per Clonotype or Exact subclonotype.

    Basically a dictionary like:

            "Clonotype1": [BarcodeAS(object),BarcodeAS(object),...],
            "Clonotype2": [BarcodeAS(object),BarcodeAS(object),...],
    """

    def __init__(self, grouped_by: ClonalGroupLevel, clonotypes_csv: AnyStr):
        """Initialize like a dictionary."""
        self.grouped_by = grouped_by
        self.clonotypes_csv = clonotypes_csv
        self._data: dict[str, list[BarcodeAS]] = {}
        self._num_chains_per_clonotype = None
        self._concordance_per_clonotype = None

    def __getitem__(self, item: str):
        """Returns a list for a given assignment."""
        return self._data[item]

    def __setitem__(self, key: str, value: list[BarcodeAS]):
        assert isinstance(value, list)
        self._data[key] = value

    def __delitem__(self, key: str):
        del self._data[key]

    def __iter__(self):
        return iter(self._data)

    def __len__(self):
        return len(self._data)

    def clonotype_size(self, key):
        return len(self._data[key])

    def clonotype_concordance(self, key, antigen_to_control):
        """Returns clonotype concordance for a given clonotype."""
        assignments = []
        for bc in self._data[key]:
            assignments.append(bc.get_assigned_antigen(antigen_to_control))
        assignments = [UNASSIGNED if a == BLANK else a for a in assignments]
        max_antigen = Counter(sorted(assignments)).most_common(1)[0]
        size = self.clonotype_size(key)
        concordance = tk_stats.robust_divide(max_antigen[1], size)

        clonotype_id = (
            key.rsplit("_", 1)[0] if self.grouped_by == ClonalGroupLevel.EXACT_SUBCLONOTYPE else key
        )
        is_canonical_pair = (
            self.num_chains_per_clonotype[clonotype_id] if clonotype_id != "None" else False
        )
        return Concordance(
            key, size, is_canonical_pair, max_antigen[0], max_antigen[1], concordance
        )

    @property
    def num_chains_per_clonotype(self):
        """Returns a dict indexed by clonotype containing metadata on num chains."""
        if self._num_chains_per_clonotype is None:
            result: dict[str, bool] = {}
            with open(self.clonotypes_csv) as infile:
                csvreader = csv.reader(infile)
                seen_header = False
                for row in csvreader:
                    if not seen_header:
                        assert row[0] == "clonotype_id"
                        assert row[3] == "cdr3s_aa"
                        seen_header = True
                        continue
                    chains = "_".join(sorted(genes.split(":")[0] for genes in row[3].split(";")))
                    result[row[0]] = chains in CANONICAL_VDJ_GENE_PAIRS
            self._num_chains_per_clonotype = result
        return self._num_chains_per_clonotype

    def concordance_per_clonotype(self, antigen_to_control):
        """Returns a dict indexed by clonotype containing metadata on antigen assigned and concordance."""
        if self._concordance_per_clonotype is None:
            result: dict[str, Concordance] = {}
            for clonotype in self._data:
                result.update(
                    {clonotype: self.clonotype_concordance(clonotype, antigen_to_control)}
                )
            order_by = {
                k: [
                    literal_eval(i) if i != "None" else 0
                    for i in k.split("clonotype")[-1].split("_")
                ]
                for k in result
            }
            result = OrderedDict(sorted(result.items(), key=lambda kv: order_by[kv[0]]))
            self._concordance_per_clonotype = result
        return self._concordance_per_clonotype

    def get_clonotype_concordance_metrics(self, antigen_to_control):
        """Compute summary metrics on clonotype (or exact_subclonotype) concordance."""
        report_prefix = rna_library.get_library_type_metric_prefix(rna_library.ANTIGEN_LIBRARY_TYPE)

        name_median_concordance_gt9 = (
            report_prefix + f"median_concordance_of_{self.grouped_by.value}_size_gt9"
        )
        name_min_concordance_gt9 = (
            report_prefix + f"lowest_concordance_of_{self.grouped_by.value}_size_gt9"
        )
        concordance_gt9 = [
            v.concordance
            for k, v in self.concordance_per_clonotype(antigen_to_control).items()
            if v.size > 9 and k != "None"
        ]

        name_median_concordance_gt9_canonical_pair = (
            report_prefix + f"median_concordance_of_{self.grouped_by.value}_size_gt9_canonical_pair"
        )
        name_min_concordance_gt9_canonical_pair = (
            report_prefix + f"lowest_concordance_of_{self.grouped_by.value}_size_gt9_canonical_pair"
        )
        concordance_gt9_canonical_pair = [
            v.concordance
            for k, v in self.concordance_per_clonotype(antigen_to_control).items()
            if v.size > 9 and k != "None" and v.canonical_pair
        ]

        name_aggregate_concordance = (
            report_prefix + f"aggregate_concordance_of_{self.grouped_by.value}"
        )

        name_aggregate_concordance_canonical_pair = (
            report_prefix + f"aggregate_concordance_of_{self.grouped_by.value}_canonical_pair"
        )

        aggregate_concordance = PercentMetric()
        aggregate_concordance_canonical_pair = PercentMetric()
        for clonotype, value in self.concordance_per_clonotype(antigen_to_control).items():
            if clonotype != "None":
                aggregate_concordance.add_value(value.num_bcs_with_assigned_antigen, value.size)
                if value.canonical_pair:
                    aggregate_concordance_canonical_pair.add_value(
                        value.num_bcs_with_assigned_antigen, value.size
                    )

        return {
            name_median_concordance_gt9: np.median(concordance_gt9) if concordance_gt9 else np.nan,
            name_min_concordance_gt9: np.min(concordance_gt9) if concordance_gt9 else np.nan,
            name_median_concordance_gt9_canonical_pair: (
                np.median(concordance_gt9_canonical_pair)
                if concordance_gt9_canonical_pair
                else np.nan
            ),
            name_min_concordance_gt9_canonical_pair: (
                np.min(concordance_gt9_canonical_pair) if concordance_gt9_canonical_pair else np.nan
            ),
            name_aggregate_concordance: aggregate_concordance.report(),
            name_aggregate_concordance_canonical_pair: aggregate_concordance_canonical_pair.report(),
        }

    @classmethod
    def from_dictionary(
        cls, data: dict[str, list[BarcodeAS]], grouped_by: ClonalGroupLevel, clonotypes_csv: AnyStr
    ):
        """Loads from a dictionary."""
        return_value = cls(grouped_by, clonotypes_csv)
        for clonotype, bc_list in data.items():
            return_value[clonotype] = bc_list
        return return_value

    def write_clonotype_concordance_csv(self, path: AnyStr, antigen_to_control):
        """Create clonotype concordance CSV."""
        with open(path, "w") as csv_handle:
            csv_writer = csv.writer(csv_handle, delimiter=",")
            csv_writer.writerow(Concordance._fields)
            for clonotype in self.concordance_per_clonotype(antigen_to_control).values():
                csv_writer.writerow(clonotype)
