#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

"""Types for loading, saving, and using feature reference data."""

from __future__ import annotations

import csv
import os
from collections.abc import Collection, Generator, ItemsView, Iterable, KeysView
from dataclasses import dataclass
from typing import TextIO

import h5py
from six import ensure_binary, ensure_str

import cellranger.h5_constants as h5_constants
import cellranger.hdf5 as cr_h5
from cellranger.feature.antigen.specificity import MHC_ALLELE, TARGETING_ANTIGEN
from cellranger.rna.library import ANTIGEN_LIBRARY_TYPE
from cellranger.targeted.targeted_constants import EXCLUDED_PROBE_ID_PREFIXES

FEATURE_TYPE = "feature_type"
# Required HDF5 datasets
REQUIRED_DATASETS = ["id", "name", FEATURE_TYPE]


@dataclass(frozen=True)
class FeatureDef:
    """A more sophisticated object for FeatureDef that allows defaults."""

    index: int
    id: bytes  # pylint: disable=invalid-name
    name: str | None
    feature_type: str
    tags: dict[str, str]

    def get_attr_list(self, headers: Iterable[str]) -> list[str]:
        """Get list of attributes corresponding to the feature."""
        ret_list: list[str] = []
        for header in headers:
            if header == REQUIRED_DATASETS[0]:
                ret_list.append(self.id.decode("utf-8"))
            elif header == REQUIRED_DATASETS[1]:
                ret_list.append(self.name)
            elif header == REQUIRED_DATASETS[2]:
                ret_list.append(self.feature_type)
            else:
                ret_list.append(self.tags.get(header, ""))
        return ret_list


# These feature tag keys are reserved for internal use.
GENOME_FEATURE_TAG = "genome"
DEFAULT_FEATURE_TAGS = [GENOME_FEATURE_TAG]
RESERVED_TAGS = DEFAULT_FEATURE_TAGS

# Defined here lib/rust/cr_types/src/reference/feature_reference.rs
HASHTAG_TAG = "hashtag"


class FeatureDefException(Exception):
    """Exception type for trying to create a `FeatureReference` with non-distinct IDs."""

    def __init__(self, msg):
        super().__init__(msg)
        self.msg = msg

    def __str__(self):
        return self.msg


class TargetFeatureSets:
    """A set of target features.

    Each target set should be a list of integers with a dictionary key representing the name.
    This name will also exist in the library_info.
    """

    def __init__(
        self, target_features: None | TargetFeatureSets | dict[str, Iterable[int]]
    ) -> None:
        self.target_feature_sets: dict[str, set[int]] = {}
        if target_features is not None:
            for name, _set in target_features.items():
                assert all(isinstance(x, int) for x in _set)
                _set = set(_set)
                self.target_feature_sets[name] = _set

    def __eq__(self, other):
        if not isinstance(other, TargetFeatureSets):
            return False
        return self.target_feature_sets == other.target_feature_sets

    def __hash__(self):
        return hash(self.target_feature_sets)

    @property
    def names(self) -> KeysView[str]:
        return self.target_feature_sets.keys()

    def items(self) -> ItemsView[str, set[int]]:
        return self.target_feature_sets.items()

    def iteritems(self) -> ItemsView[str, set[int]]:
        """Returns a view of the items.

        .. deprecated::4.0
            Use `items` instead.  Python 3 doesn't do `iteritems`.
        """
        return self.target_feature_sets.items()

    def get(self, key: str) -> set[int] | None:
        return self.target_feature_sets.get(key)

    def update(self, other: dict[str, set[int]]) -> None:
        self.target_feature_sets.update(other)

    def count_target_feature_indices(self) -> int:
        """Returns the number of distinct feature indices.

        This is more efficient than len(get_target_feature_indices()
        because it doesn't bother to sort the indices.
        """
        return len(set.union(*self.target_feature_sets.values()))

    def get_target_feature_indices(self) -> list[int]:
        return sorted(set.union(*self.target_feature_sets.values()))

    def subset_to_reduced_features(
        self, old_to_new_index_translation: dict[int, int]
    ) -> TargetFeatureSets:
        """Recreates a list of target features by subsetting it.

        If the Feature Reference holding this object is ever subset, or the index values are changed,
        we need to update the target_feature_sets to hold these new values instead.  This function
        recreates a list of target features by subsetting it to only those present in a translation
        table that maps old index positions to new ones.

        Args:
            old_to_new_index_translation: a dictionary that tells what old index position is being converted
            to a new index position

        Returns:
            A new TargetFeatureSets
        """
        assert isinstance(old_to_new_index_translation, dict)
        new_target_features = {}
        for name, indices in self.target_feature_sets.items():
            new_vals = [
                old_to_new_index_translation[old_index]
                for old_index in indices
                if old_index in old_to_new_index_translation
            ]
            new_target_features[name] = new_vals
        return TargetFeatureSets(new_target_features)


class FeatureReference:  # pylint: disable=too-many-public-methods
    """Store a list of features (genes, antibodies, etc)."""

    def __init__(
        self,
        feature_defs: list[FeatureDef],
        all_tag_keys: list[str],
        target_features: None | TargetFeatureSets | dict[str, Iterable[int]] = None,
    ):
        """Create a FeatureReference.

        Args:
            feature_defs (list of FeatureDef): All feature definitions.
            all_tag_keys (list of str): All optional tag keys.
            target_features (dictionary of list of int): Optional target set(s). Each target set
            should be a list of integers with a dictionary key representing the name. This name
            will also exist in the library_info.
        """
        self.feature_defs = feature_defs
        self.all_tag_keys = all_tag_keys

        if target_features is not None:
            self.target_features = TargetFeatureSets(target_features)
        else:
            self.target_features = None

        # Assert uniqueness of feature IDs
        id_map = {}
        for fdef in self.feature_defs:
            if fdef.id in id_map:
                this_fd_str = f"ID: {fdef.id}; name: {fdef.name}; type: {fdef.feature_type}"
                seen_fd_str = f"ID: {id_map[fdef.id].id}; name: {id_map[fdef.id].name}; type: {id_map[fdef.id].feature_type}"
                raise FeatureDefException(
                    "Found two feature definitions with the same ID: "
                    f"({this_fd_str}) and ({seen_fd_str}). All feature IDs must be distinct."
                )
            id_map[fdef.id] = fdef

        self.id_map: dict[bytes, FeatureDef] = id_map

    def get_count_of_feature_type(self, feature_type: str) -> int:
        """Count how many features in the matrix are of a given type.

        (e.g. "Gene Expression")
        """
        total = 0
        for feature in self.feature_defs:
            if feature.feature_type == feature_type:
                total += 1
        return total

    def __eq__(self, other):
        return (
            self.feature_defs == other.feature_defs
            and self.all_tag_keys == other.all_tag_keys
            and (self.target_features is None) == (other.target_features is None)
            and self.target_features == other.target_features
        )

    def __hash__(self):
        if self.target_features is None:
            return hash((self.feature_defs, self.all_tag_keys))
        return hash((self.feature_defs, self.all_tag_keys, self.target_features))

    def __ne__(self, other):
        return not self == other

    def get_feature_ids_excluding_deprecated_probes(self) -> list[bytes]:
        """Return the list of feature IDs excluding deprecated probes."""
        return [f.id for f in self.feature_defs if not f.id.startswith(EXCLUDED_PROBE_ID_PREFIXES)]

    def has_deprecated_probes(self) -> bool:
        """Return true if there are deprecated probes in features."""
        return any(f.id.startswith(EXCLUDED_PROBE_ID_PREFIXES) for f in self.feature_defs)

    def get_feature_types_excluding_deprecated_probes(self) -> list[str]:
        """Return the list of feature types excluding deprecated probes."""
        return [
            f.feature_type
            for f in self.feature_defs
            if not f.id.startswith(EXCLUDED_PROBE_ID_PREFIXES)
        ]

    def get_antigen_control(self) -> tuple | None:
        """Extract antigen control feature from feature reference."""
        control = None
        if TARGETING_ANTIGEN in self.all_tag_keys:
            control = [
                (f.index, f.id, f.name, f.tags.get(MHC_ALLELE, None))
                for f in self.select_features_by_type(ANTIGEN_LIBRARY_TYPE).feature_defs
                if f.tags[TARGETING_ANTIGEN] == "False"
            ]
            assert len(control) <= 1
            control = control[0] if len(control) == 1 else None
        return control

    def get_antigen_capture_content(self) -> frozenset[list[tuple]]:
        """Extract the antigen capture content from feature reference.

        The conent extracted in this function is what we expect to match between
        inputs to aggr.
        """
        lst = []
        for f in self.select_features_by_type(ANTIGEN_LIBRARY_TYPE).feature_defs:
            lst.append(
                (
                    f.index,
                    f.id,
                    f.name,
                    f.tags.get(MHC_ALLELE, None),
                    f.tags.get(TARGETING_ANTIGEN, None),
                )
            )

        return frozenset(lst)

    def equals_antigen_capture_content(
        self, other: FeatureReference, is_pd: bool
    ) -> tuple[bool, str | None]:
        """Checks if two feature references have the same antigen capture features.

        returns appropriate error message in case of a mismatch
        """
        antigen_content_1 = self.get_antigen_capture_content()
        antigen_content_2 = other.get_antigen_capture_content()
        if len(antigen_content_1.difference(antigen_content_2)) != 0:
            sorted_antigen_content_1 = sorted(list(antigen_content_1), key=lambda x: x[0])
            sorted_antigen_content_2 = sorted(list(antigen_content_2), key=lambda x: x[0])
            # number of antigen capture features should match, checked in equals_ignore_target_set
            assert len(sorted_antigen_content_1) == len(sorted_antigen_content_2)
            for i, _ in enumerate(sorted_antigen_content_1):
                # index, id, and name must match since they are checked in equals_ignore_target_set
                assert sorted_antigen_content_1[i][0] == sorted_antigen_content_2[i][0]
                assert sorted_antigen_content_1[i][1] == sorted_antigen_content_2[i][1]
                assert sorted_antigen_content_1[i][2] == sorted_antigen_content_2[i][2]
                # check matching mhc_allele
                if sorted_antigen_content_1[i][3] != sorted_antigen_content_2[i][3]:
                    return (
                        False,
                        "The datasets you are trying to aggregate have incompatible "
                        "MHC alleles for the same control feature id. Please re-run the original "
                        "multi pipelines with uniform [antigen-specificity] sections.",
                    )
                # check matching control
                if (sorted_antigen_content_1[i][4] != sorted_antigen_content_2[i][4]) and not is_pd:
                    return (
                        False,
                        "The datasets you are trying to aggregate have incompatible "
                        "control feature ids. Please re-run the original multi pipelines with "
                        "uniform [antigen-specificity] sections.",
                    )

        return True, None

    def equals_ignore_target_set(self, other: FeatureReference) -> bool:
        """Checks if two feature references are equal, ignoring missing targets.

        Only takes into account equality of the target set if it is present in
        both feature references. Useful for targeted feature reference
        compatibility checks in preflights. Excludes deprecated probes.
        """
        ignore_target_set = (not self.has_target_features()) or (not other.has_target_features())
        compatible_target_sets = ignore_target_set or (
            self.get_target_feature_indices() == other.get_target_feature_indices()
        )
        # exclude hashtag
        self_all_tags = [tag for tag in self.all_tag_keys if tag != HASHTAG_TAG]
        other_all_tags = [tag for tag in other.all_tag_keys if tag != HASHTAG_TAG]

        return (
            self.get_feature_types_excluding_deprecated_probes()
            == other.get_feature_types_excluding_deprecated_probes()
            and self.get_feature_ids_excluding_deprecated_probes()
            == other.get_feature_ids_excluding_deprecated_probes()
            and self_all_tags == other_all_tags
            and compatible_target_sets
        )

    def has_compatible_target_set(self, other: FeatureReference) -> bool:
        """Checks if two feature references have compatible target sets.

        Useful for targeted feature reference compatibility checks in preflights.

        Returns:
            True if either target set does not exist or if they are identical.
        """
        ignore_target_set = (not self.has_target_features()) or (not other.has_target_features())
        compatible_target_sets = ignore_target_set or (
            self.get_target_feature_indices() == other.get_target_feature_indices()
        )
        return compatible_target_sets

    def add_tags(
        self,
        new_tags: list[str],
        new_labels: Collection[Collection[str]] | None = None,
    ):
        """Add new tags and corresponding labels to self.

        If new labels are None, empty strings are supplied by default.

        Args:
            new_tags (list[str]):  a list of new tags
            new_labels (Optional[Collection[Collection[str]]], optional): per feature
                list of label values corresponding to the new tags. Defaults to None.

        Raises:
            ValueError: Tag already present
            ValueError: Number of labels is not equal to number of features in feature ref
        """
        assert len(new_tags) > 0
        for tag in new_tags:
            if tag in self.all_tag_keys:
                raise ValueError("tag {} is already present in feature_ref")

        if new_labels is None:
            for tag in new_tags:
                self.add_tag(tag, {})
        else:
            if len(self.feature_defs) != len(new_labels):
                raise ValueError(
                    "number of labels does not match number of features in feature_ref"
                )
            for labels in new_labels:
                assert len(labels) == len(new_tags)
            for ind, tag in enumerate(new_tags):
                label_list = (x[ind] for x in new_labels)
                feature_id_to_label_dict = {
                    f.id: label for (f, label) in zip(self.feature_defs, label_list)
                }
                self.add_tag(tag, feature_id_to_label_dict)

    def update_tag(
        self,
        tag_to_update: str,
        feature_id_to_label_dict: dict[str, str] | dict[bytes, str],
    ):
        """Updates a tag in the feature reference.

        The features that need to
        be updated with their new tag values are passed in.

        Args:
            tag_to_update (str): a tag to update. Should be in feature_ref
            feature_id_to_label_dict (Union[dict[str, str], dict[bytes, str]]): a
                dict of {feature_id: new_value_of_tag} to update. The features that need to
                be updated mapped to their new tag values.

        Raises:
            ValueError: If tag passed in to be updated does not exist in the feature ref
        """
        if tag_to_update not in self.all_tag_keys:
            raise ValueError(
                "Tag to update not in Feature ref. "
                + f"Tag to update {tag_to_update} . "
                + f"Tags in matrix {self.all_tag_keys}"
            )

        for f_def in self.feature_defs:
            if f_def.id in feature_id_to_label_dict:
                f_def.tags[tag_to_update] = feature_id_to_label_dict.get(f_def.id, "")

    def add_tag(
        self,
        tag_to_add: str,
        feature_id_to_label_dict: dict[str, str] | dict[bytes, str],
        default_tag_value: str | None = "",
    ):
        """Adds a tag to a Feature Reference.

        Args:
            tag_to_add (str): name of the tag to add
            feature_id_to_label_dict (Union[dict[str, str], dict[bytes, str]]):dict of
                {feature_id: new_value_of_tag} to update. The features that need to
                be updated mapped to their new tag values
            default_tag_value (Optional[str], optional): Value of the tag for features
            not in feature_id_to_label_dict. Defaults to "".

        Raises:
            ValueError: If tag passed in to be added already exists in the feature ref
        """
        if tag_to_add in self.all_tag_keys:
            raise ValueError(
                "Tag to be added already in Feature ref. "
                + f"Tag to update {tag_to_add} . "
                + f"Tags in matrix {self.all_tag_keys}"
            )

        for f_def in self.feature_defs:
            f_def.tags[tag_to_add] = feature_id_to_label_dict.get(f_def.id, default_tag_value)
        self.all_tag_keys.append(tag_to_add)

    @staticmethod
    def join(feature_ref1: FeatureReference, feature_ref2: FeatureReference) -> FeatureReference:
        """Concatenate two feature references, requires unique ids and identical tags."""
        assert feature_ref1.all_tag_keys == feature_ref2.all_tag_keys
        feature_defs1 = feature_ref1.feature_defs
        feature_defs2 = feature_ref2.feature_defs

        if feature_ref1.target_features is None:
            combined_target_features = feature_ref2.target_features
        elif feature_ref2.target_features is None:
            combined_target_features = feature_ref1.target_features
        else:
            combined_target_features = feature_ref1.target_features
            # if feature_ref2 has the same keys, they will be over-written
            combined_target_features.update(feature_ref2.target_features)
        return FeatureReference(
            feature_defs=feature_defs1 + feature_defs2,
            all_tag_keys=feature_ref1.all_tag_keys,
            target_features=combined_target_features,
        )

    @classmethod
    def empty(cls) -> FeatureReference:
        return cls(feature_defs=[], all_tag_keys=[], target_features=None)

    def get_num_features(self) -> int:
        return len(self.feature_defs)

    def get_feature_ids_by_type(self, feature_type: str) -> list[bytes]:
        """Return a list of feature ids of a particular feature type (e.g. "Gene Expression")."""
        return [f.id for f in self.feature_defs if f.feature_type == feature_type]

    def get_indices_for_type(self, feature_type: str) -> list[int]:
        return [
            feature.index for feature in self.feature_defs if feature.feature_type == feature_type
        ]

    def get_genomes(self, feature_type: str | None = None) -> list[str]:
        """Get sorted list of genomes.

        Empty string is for reverse compatibility.

        A specific feature type can optionally be specified.
        """
        genomes = {
            f.tags.get(GENOME_FEATURE_TAG, "")
            for f in self.feature_defs
            if (feature_type is None or f.feature_type == feature_type)
        }
        return sorted(genomes)

    def has_target_features(self) -> bool:
        return self.target_features is not None

    def count_target_feature_indices(self):
        """Returns the number of feature indices.

        If there are no target features, returns zero.

        This is more efficient than len(get_target_feature_indices()) because
        it doesn't bother sorting the indices.
        """
        return (
            self.target_features.count_target_feature_indices()
            if self.target_features is not None
            else 0
        )

    def get_target_feature_indices(self) -> list[int] | None:
        """Gets the indices of on-target features within the FeatureReference.

        Returns None if there is no target set.
        """
        if not self.has_target_features():
            return None
        else:
            return self.target_features.get_target_feature_indices()

    def get_target_feature_ids(self) -> list[bytes] | None:
        """Gets the feature ids of on-target features.

        Returns None if there is no target set.
        """
        if not self.has_target_features():
            return None
        else:
            return [self.feature_defs[i].id for i in self.get_target_feature_indices()]

    def select_features_by_type(self, feature_type: str) -> FeatureReference:
        """Returns all features that match a given type.

        Args:
            feature_type:
        """
        indices = [i for i, fd in enumerate(self.feature_defs) if fd.feature_type == feature_type]
        return self.select_features(indices)

    def get_feature_indices_by_genome(self, genome: str) -> list[int]:
        """Returns the indices of features within the FeatureReference with genome.

        Args:
            genome: Genome name
        """
        return [
            i
            for i, fd in enumerate(self.feature_defs)
            if fd.tags.get(GENOME_FEATURE_TAG, "") == genome
        ]

    def select_features_by_genome(self, genome: str) -> FeatureReference:
        """Returns all features that match a given genome.

        Args:
            genome: Genome name
        """
        return self.select_features(self.get_feature_indices_by_genome(genome))

    def has_feature_type(self, feature_type: str) -> bool:
        """Determines if a feature type is present in the FeatureRef.

        Args:
            feature_type:

        Returns:
            True | False
        """
        for feat in reversed(self.feature_defs):
            if feat.feature_type == feature_type:
                return True
        return False

    def get_feature_names(self) -> list[str]:
        """Get a list of feature names.

        Returns:
            The names of all of the featture_defs.
        """
        return [x.name for x in self.feature_defs]

    def get_feature_ids(self) -> list[str] | list[bytes]:
        """Get a list of feature IDs.

        Returns:
            The IDs of all of the feature_defs.
        """
        return [x.id for x in self.feature_defs]

    def select_features(self, indices: Iterable[int]) -> FeatureReference:
        """Create a new FeatureReference that only contains features with the given indices.

        Any target sets present are updated, but the target_features attribute in the output
        is set to None if all target features were removed.
        """
        old_defs = [self.feature_defs[i] for i in indices]
        new_defs: list[FeatureDef] = []
        translation_map = {}
        for i, f_def in enumerate(old_defs):
            translation_map[f_def.index] = i
            new_defs.append(
                FeatureDef(
                    index=i,
                    id=f_def.id,
                    name=f_def.name,
                    feature_type=f_def.feature_type,
                    tags=f_def.tags,
                )
            )

        if self.has_target_features():
            new_target_features = self.target_features.subset_to_reduced_features(translation_map)
            # If we have removed all targeted features, we toss the empty target sets
            if len(new_target_features.get_target_feature_indices()) == 0:
                new_target_features = None
        else:
            new_target_features = None
        return FeatureReference(
            feature_defs=new_defs,
            all_tag_keys=self.all_tag_keys,
            target_features=new_target_features,
        )

    def to_hdf5(self, group: h5py.Group) -> None:
        """Write to an HDF5 group."""
        # Write required datasets
        for col in REQUIRED_DATASETS:
            data = [getattr(f, col) for f in self.feature_defs]
            cr_h5.create_hdf5_string_dataset(group, col, data, compression=True)

        # Write tag datasets
        for col in self.all_tag_keys:
            # Serialize missing data as empty unicode string
            data = [f.tags.get(col, "") for f in self.feature_defs]
            cr_h5.create_hdf5_string_dataset(group, col, data, compression=True)

        # Write target_features as a new sub-group
        if self.target_features is not None:
            target_set_group = group.create_group(h5_constants.H5_TARGET_SET_ATTR)
            for key, val in self.target_features.items():
                cr_h5.create_hdf5_string_dataset(
                    target_set_group, key, [str(x) for x in sorted(val)], compression=True
                )

        # Record names of all tag columns
        cr_h5.create_hdf5_string_dataset(group, "_all_tag_keys", self.all_tag_keys)

    def to_csv(self, handle: TextIO) -> None:
        """Write to a csv file (only columns specified in CSV_HEADER).

        Args:
            handle: the file-like object the feature ref will be written to
        """
        csv_headers = REQUIRED_DATASETS + self.all_tag_keys
        writer = csv.writer(handle)
        # Write header:
        writer.writerow(csv_headers)
        for feature in self.feature_defs:
            writer.writerow(feature.get_attr_list(csv_headers))

    @classmethod
    def from_hdf5(cls, group: h5py.Dataset) -> FeatureReference:
        """Load from an HDF5 group.

        Args:
            group (h5py.Dataset): Group to load from.

        Returns:
            feature_ref (FeatureReference): New object.
        """

        # FIXME: ordering may not be guaranteed in python3
        def _load_str(node: h5py.Dataset) -> list[str]:
            if node.shape is None:
                return []
            if node.dtype.char == cr_h5.STR_DTYPE_CHAR:
                memoize = node.name in ("/features/feature_type", "/features/genome")
                return cr_h5.read_hdf5_string_dataset(node, memoize)
            else:
                return node[:]

        def _format_path(path: str | bytes) -> str:
            """Strip off leading slash and convert to str (could be unicode)."""
            path = path[1:]
            return path.decode() if isinstance(path, bytes) else path

        def _h5py_dataset_iterator(
            group: h5py.Group | h5py.Dataset, prefix: str = ""
        ) -> Generator[tuple[str, list[str]], None, None]:
            for key in group:
                item = group[key]
                path = f"{prefix}/{key}"
                if isinstance(item, h5py.Dataset):
                    yield _format_path(path), _load_str(item)
                elif isinstance(item, h5py.Group):
                    yield from _h5py_dataset_iterator(item, path)

        data: dict[str, list[str]] = dict(_h5py_dataset_iterator(group))

        # Load Tag Keys
        all_tag_keys = data["_all_tag_keys"]

        # Load Target Sets, if they exist
        target_features = {
            os.path.basename(key): [int(x) for x in val]
            for key, val in data.items()
            if h5_constants.H5_TARGET_SET_ATTR in ensure_binary(key)
        }
        if len(target_features) == 0:
            target_features = None

        # Build FeatureDefs
        feature_defs: list[FeatureDef] = []
        num_features = len(data[REQUIRED_DATASETS[0]])

        for i in range(num_features):
            tags = {ensure_str(k): ensure_str(data[k][i]) for k in all_tag_keys if data[k][i]}
            feature_defs.append(
                FeatureDef(
                    id=ensure_binary(data["id"][i]),
                    index=i,
                    name=data["name"][i],
                    feature_type=data[FEATURE_TYPE][i],
                    tags=tags,
                )
            )

        return cls(feature_defs, all_tag_keys=all_tag_keys, target_features=target_features)
