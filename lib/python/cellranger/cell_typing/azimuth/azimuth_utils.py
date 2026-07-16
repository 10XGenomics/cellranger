"""Utility functions used with azimuth. Code taken from panhumanpy and edited."""

import numpy as np
from scipy.sparse import csr_matrix
from sklearn.calibration import LabelEncoder

from cellranger.cell_typing.azimuth.azimuth_ref_paths import (
    AZIMUTH_INFERENCE_ENCODERS_PATHS_SORTED_LIST,
)

#######################################################################
################## functions ##########################################
#######################################################################
#######################################################################


def check_normalization(matrix, normalization_override, norm_check_batch_size):
    """Check if a matrix is normalized based on examining a subset of values.

    This function determines if a matrix has been normalized by checking
     if any values are non-integer, which would indicate normalization
     has occurred.

    Parameters
    ----------
    matrix : array-like
        The matrix to check for normalization.
    normalization_override : bool
        If True, bypass the normalization check and return True.
    norm_check_batch_size : int
        Maximum number of rows to check for efficiency.

    Returns:
    -------
    bool
        True if the matrix is normalized (contains non-integer values)
        or if normalization_override is True, False otherwise.
    """
    if matrix.shape[0] > norm_check_batch_size:
        mat = matrix[:norm_check_batch_size, :]
    else:
        mat = matrix

    mat = mat.toarray()
    mat_floor = np.floor(mat)

    integer_counts_detect = not np.any((mat_floor - mat) != 0.0)
    print("Integer counts detected by panhumanpy.ANNotate : " f"{integer_counts_detect}")

    if normalization_override:
        return True
    else:
        return not integer_counts_detect


def normalize(matrix, normalization_override, norm_check_batch_size):
    """Normalizes a gene expression count matrix using log1p transformation.

    This function first checks whether the provided matrix is already
    normalized by inspecting a batch of cells (specified by
    norm_check_batch_size). If the matrix is not normalized and
    normalization_override is False, the function scales each cell to
    10,000 total counts and applies a log1p transformation to the data.
    The input matrix is assumed to be a sparse matrix (e.g., CSR format)
    that supports the operations .sum(), .multiply(), and .tocsr().

    A caveat of the current formulation of this function is that it
    merely checks if a certain batch of cells (as specified) has
    integer counts data or not, and if not it assumes that the data is
    log1p normalized, which need not be the case.

    Parameters
    ----------
    matrix : scipy.sparse matrix
        Gene expression count matrix with shape (num_cells, num_genes).
    normalization_override : bool
        If True, assumes that the matrix is already normalized and
        bypasses the scaling and log1p transformation.
    norm_check_batch_size : int, optional
        The number of cells used to determine whether the matrix is
        normalized.

    Returns:
    -------
    scipy.sparse.csr_matrix
        The normalized gene expression matrix.
    bool
        Whether the provided matrix was normalized initially.
    """
    check_norm = check_normalization(matrix, normalization_override, norm_check_batch_size)

    if not check_norm:
        total_counts = matrix.sum(axis=1)
        total_counts = np.array(total_counts).reshape(-1, 1)

        if 0 in total_counts.flatten():
            print("Cells with 0 counts across the entire feature panel found.")
        total_counts[total_counts == 0] = 1

        # Scale each nonzero by its row factor directly on the CSR data array,
        # then log1p. Numerically identical to matrix.multiply(10000/total)
        # followed by log1p (same data[k] * (10000/total[row(k)])), but avoids
        # building an intermediate COO and converting it back to CSR.
        matrix = matrix.tocsr()
        row_scale = np.repeat((10000 / total_counts).ravel(), np.diff(matrix.indptr))
        new_data = np.log1p(matrix.data * row_scale)
        # Reuse the original index arrays (the source matrix is discarded here);
        # only the data array is new.
        matrix = csr_matrix((new_data, matrix.indices, matrix.indptr), shape=matrix.shape)

    print(f"Log-normalization performed by panhumanpy.ANNotate: {not check_norm}")
    if not check_norm:
        print("To override log-normalization, set normalization_override=False")

    return matrix, check_norm


def reorder_subset_data_matrix(
    data_matrix, query_features, feature_panel_template, common_features=None
):
    """Reorders, subsets, and extends a sparse data matrix to match a template feature order.

    This function takes a data matrix with features defined by
    query_features and transforms it to match the feature order
    specified in feature_panel_template. If features in
    feature_panel_template are missing from query_features, the
    function adds zero columns for these features. Extra query features
    not present in the feature_panel_template are dropped.

    Parameters
    ----------
    data_matrix : scipy.sparse.csr_matrix
        The input data matrix to be reordered.
    query_features : list
        List of feature names corresponding to the columns in
        data_matrix.
    feature_panel_template : list
        List of feature names defining the target order for output.
    common_features : set, optional
        Set of features present in both query_features and
        feature_panel_template. If None, will be computed as the
        intersection of both lists.

    Returns:
    -------
    scipy.sparse.csr_matrix
        A reordered, subsetted, and possibly extended data matrix with
        the same number of rows (cells) as the input, but with columns
        matching the order and features of feature_panel_template.

    Notes:
    -----
    - The function assumes that data_matrix is a scipy sparse matrix in
        CSR format.
    - Zero columns are added for features that are in
        feature_panel_template but not in query_features.
    - Columns are reordered to match the exact order in
        feature_panel_template.
    """
    if not common_features:
        common_features = set(query_features).intersection(set(feature_panel_template))

    # First-occurrence column index for each query feature (matches the
    # semantics of the original `list.index`, which returns the first match).
    first_idx = {}
    for i, name in enumerate(query_features):
        if name not in first_idx:
            first_idx[name] = i

    # For each template column, the source query column (if the feature is
    # present) and the destination template position. Template features not in
    # the query become all-zero columns (absent from the CSR).
    present_dst = []
    present_src = []
    for k, name in enumerate(feature_panel_template):
        if name in common_features:
            present_dst.append(k)
            present_src.append(first_idx[name])

    # Select the present query columns (preserving cell rows) and relabel their
    # column indices to the template positions. This is dense-identical to the
    # original "append zero columns then fancy-index" but avoids the O(n^2)
    # list.index scan and the zero-column hstack/reindex.
    sub = data_matrix[:, present_src].tocsr()
    n_cols = len(feature_panel_template)
    if sub.nnz:
        dst = np.asarray(present_dst, dtype=sub.indices.dtype)
        new_indices = dst[sub.indices]
    else:
        new_indices = sub.indices
    reordered_data_matrix = csr_matrix(
        (sub.data, new_indices, sub.indptr),
        shape=(data_matrix.shape[0], n_cols),
    )

    return reordered_data_matrix


def if_full_consistent_hierarchy(cell_label, max_depth):
    """Checks if a label hierarchy is consistent.

    Returns a boolean indicating whether the hierarchical
    predictions returned by the model for a given cell form an
    internally consistent hierarchy or not.
    """
    res = True
    for i in range(max_depth - 1):
        res = cell_label[i].split("|")[-1] == (cell_label[i + 1].split("|")[-2])

        if res is False:
            break

    return res


def _found_parts(array_label, depth, max_depth):
    """Per-level component (or None) used to assemble a combined label.

    Mirrors the inner scan of :func:`comb_label`: for each hierarchy level
    ``i`` in ``range(depth)``, take the first non-empty ``split("|")[i]``
    among ``array_label[i : i + (max_depth - i)]``, returning the found
    component or ``None`` if no non-empty component exists for that level.

    Each row of ``array_label`` is split exactly once and cached, so the
    scan costs O(max_depth) splits instead of the O(max_depth**2) the
    original per-level ``comb_label`` calls incurred.
    """
    split_cache = [s.split("|") for s in array_label]
    parts = []
    for i in range(depth):
        comp = None
        for j in range(max_depth - i):
            row = split_cache[i + j]
            if i < len(row) and row[i] != "":
                comp = row[i]
                break
        parts.append(comp)
    return parts


def comb_label(array_label, depth, max_depth):
    """Combine hierarchical label components from a list of hierarchically split label strings.

    This function takes a list of label strings, passed through the
    argument 'array_label', each expected to contain hierarchical
    components separated by the '|' character in a cumulative manner,
    and produces a single concatenated label. Each list corresponds to a
    hierarchical label for a cell split into components in a cumulative
    manner.

    The i-th element (0-indexed) of the label is anticipated to have i
    '|' characters and none of the strings are anticipated to begin with
     '|'. For each hierarchy level (from 0 to depth-1), the function
    looks for the first non-empty component between two the '|' characters
     approproite for that level among a subset of the labels in the
     `array_label` argument. The search for each level starts at index i
     and considers up to (max_depth - i) subsequent elements. Once a
     non-empty component is found, it is appended to the output string,
     followed by a '|' delimiter. This process is iterated up to the
     specified depth. 'depth' is also anticipated to be less than or
     or equal to 'max_depth'. After processing all levels, the trailing
     delimiter is removed.

    Example 1: an internally consistent tree.
        array_label = ["A", "A|B", "A|B|C", "A|B|C|D","A|B|C|D|",
                        "A|B|C|D||"]
        depth = 5
        max_depth = 6

        Level 0: Search for non-empty at positions [0, 1, 2, 3, 4, 5]
                (buffer = 6)
                - For j = 0: from array_label[0] ("A"), split yields
                    ["A"]. Index 0 is "A" → out = "A|"

        Level 1: Search for non-empty at positions [1, 2, 3, 4, 5]
                (buffer = 5)
                - For j = 0: from array_label[1] ("A|B"), split yields
                    ["A", "B"]. Index 1 is "B" → out = "A|B|"

        Level 2: Search for non-empty at positions [2, 3, 4, 5]
                (buffer = 4)
                - For j = 0: from array_label[2] ("A|B|C"), split yields
                 ["A", "B", "C"]. Index 2 is "C" → out = "A|B|C|"

        Level 3: Search for non-empty at positions [3, 4, 5]
                (buffer = 3)
                - For j = 0: from array_label[3] ("A|B|C|D"), split
                yields ["A", "B", "C", "D"].
                Index 3 is "D" → out = "A|B|C|D|"

        Level 4: Search for non-empty at positions [4, 5] (buffer = 2)
                - For j = 0: from array_label[4] ("A|B|C|D|"), split
                    yields ["A", "B", "C", "D", ""]. Index 4 is empty.
                - For j = 1: from array_label[5] ("A|B|C|D||"), split
                    yields ["A", "B", "C", "D", "", ""]. Index 4 is
                    still empty.
                → No non-empty component found; out remains "A|B|C|D|"

        After processing all levels, the trailing delimiter is removed.

        Result: "A|B|C|D"

    Example 2: an internally inconsistent tree (that can still be output)
                by the model.
        array_label = ["A", "A|B", "A|B|C", "A|B|C1|","A|B|C|D|",
                        "A|B|C1|||"]
        depth = 6
        max_depth = 6

        Level 0: Search for non-empty at positions [0, 1, 2, 3, 4, 5]
                (buffer = 6)
                - For j = 0: from array_label[0] ("A"), splitting yields
                 ["A"]. Index 0 is "A" (non-empty) → out becomes "A|"

        Level 1: Search for non-empty at positions [1, 2, 3, 4, 5]
                (buffer = 5)
                - For j = 0: from array_label[1] ("A|B"), splitting
                    yields ["A", "B"]. Index 1 is "B" (non-empty) → out
                    becomes "A|B|"

        Level 2: Search for non-empty at positions [2, 3, 4, 5]
                (buffer = 4)
                - For j = 0: from array_label[2] ("A|B|C"), splitting
                    yields ["A", "B", "C"].
                    Index 2 is "C" (non-empty) → out becomes "A|B|C|"

        Level 3: Search for non-empty at positions [3, 4, 5]
                (buffer = 3)
                - For j = 0: from array_label[3] ("A|B|C1|"), splitting
                    yields ["A", "B", "C1", ""].
                    Index 3 is "" (empty).
                - For j = 1: from array_label[4] ("A|B|C|D|"), splitting
                     yields ["A", "B", "C", "D", ""].
                    Index 3 is "D" (non-empty) → out becomes "A|B|C|D|"

        Level 4: Search for non-empty at positions [4, 5] (buffer = 2)
                - For j = 0: from array_label[4] ("A|B|C|D|"), splitting
                    yields ["A", "B", "C", "D", ""].
                    Index 4 is "" (empty).
                - For j = 1: from array_label[5] ("A|B|C1|||"),
                    splitting yields ["A", "B", "C1", "", "", ""].
                    Index 4 is "" (empty).
                → No non-empty component found; out remains "A|B|C|D|"

        Level 5: Search for non-empty at position [5] (buffer = 1)
                - For j = 0: from array_label[5] ("A|B|C1|||"),
                    splitting yields ["A", "B", "C1", "", "", ""].
                    Index 5 is "" (empty).
                 → No non-empty component found.

        After processing all levels, the trailing delimiter is removed.

        Result: "A|B|C|D"

    Example 3: a more complicated internally inconsistent tree that can
        still be a model output.

        array_label = ["A", "A|B", "A|B||", "A|B|C1|", "A|B|C1|D1|",
                        "A|B|C1|||"]
        depth = 6
        max_depth = 6

        Level 0 : Search at positions [0, 1, 2, 3, 4, 5] (buffer = 6)
                - j = 0: array_label[0] is "A"
                    → "A".split('|') gives ["A"]
                    → Element at index 0 is "A" (non-empty)
                    → out becomes "A|"

        Level 1 : Search at positions [1, 2, 3, 4, 5] (buffer = 5)
                - j = 0: array_label[1] is "A|B"
                    → "A|B".split('|') gives ["A", "B"]
                    → Element at index 1 is "B" (non-empty)
                    → out becomes "A|B|"

        Level 2 : Search at positions [2, 3, 4, 5] (buffer = 4)
                - j = 0: array_label[2] is "A|B||"
                    → "A|B||".split('|') gives ["A", "B", "", ""]
                    → Element at index 2 is "" (empty)
                - j = 1: array_label[3] is "A|B|C1|"
                    → "A|B|C1|".split('|') gives ["A", "B", "C1", ""]
                    → Element at index 2 is "C1" (non-empty)
                    → out becomes "A|B|C1|"

        Level 3 : Search at positions [3, 4, 5] (buffer = 3)
                - j = 0: array_label[3] is "A|B|C1|"
                    → "A|B|C1|".split('|') gives ["A", "B", "C1", ""]
                    → Element at index 3 is "" (empty)
                - j = 1: array_label[4] is "A|B|C1|D1|"
                    → "A|B|C1|D1|".split('|') gives ["A", "B", "C1",
                                                            "D1", ""]
                    → Element at index 3 is "D1" (non-empty)
                    → out becomes "A|B|C1|D1|"

        Level 4 : Search at positions [4, 5] (buffer = 2)
                - j = 0: array_label[4] is "A|B|C1|D1|"
                    → "A|B|C1|D1|".split('|') gives ["A", "B", "C1",
                                                            "D1", ""]
                    → Element at index 4 is "" (empty)
                - j = 1: array_label[5] is "A|B|C1|||"
                    → "A|B|C1|||".split('|') gives ["A", "B", "C1", "",
                                                                "", ""]
                    → Element at index 4 is "" (empty)
                    → No non-empty component found; out remains
                                                    "A|B|C1|D1|"

        Level 5 : Search at position [5] (buffer = 1)
                - j = 0: array_label[5] is "A|B|C1|||"
                    → "A|B|C1|||".split('|') gives ["A", "B", "C1", "",
                                                                "", ""]
                    → Element at index 5 is "" (empty)
                    → No non-empty component found.

        After processing all levels, the trailing '|' is removed.

        Final Result: "A|B|C1|D1"

    Parameters
    ----------
    array_label : list of str
        List of hierarchical label strings. Each string should contain
        components separated by '|'.
    depth : int
        Number of hierarchical levels to combine.
    max_depth : int
        Maximum number of label elements to consider when building the
        combined label.

    Returns:
    -------
    str
        A concatenated label string built from the appropriate components
        with '|' delimiters between hierarchy levels.

    Notes:
    -----
    - Some internally consistent hierarchies cases may be merged
      to give a sensible output, which hides the inconsistency.
      See Example 2 above.
    - This function assumes that each element in array_label
      contains at least 'depth' number of '|' delimiters.
    - If no non-empty component is found for a level, that level
      will be blank in the output label.
    """
    parts = _found_parts(array_label, depth, max_depth)
    return "|".join(p for p in parts if p is not None)


def abs_labels(hierarchical_labels_array, max_depth):
    """Compute absolute hierarchical labels for each cell up to each level.

    This function processes an array of hierarchical labels (one per
    cell) and computes an "absolute" label for each hierarchy level from
      1 up to max_depth. For each level, the function combines label
      components using the comb_label function. If the combined label
      for a cell does not contain enough components (i.e. the number of
      '|' delimiters is less than the expected level), the label is
      replaced with "NA".

    Parameters
    ----------
    hierarchical_labels_array : list
        A list of hierarchical labels, where each element represents the
        label of a cell. Each cell label should be in a format that is
        compatible with the comb_label function (e.g., a list or string
        with components separated by '|').
    max_depth : int
        The maximum number of hierarchical levels available in the labels.

    Returns:
    -------
    list
        A list of length max_depth. The i-th element is a list of
        absolute labels (strings) for each cell at hierarchical
        level i+1 (1-indexed).

    Notes:
    -----
    - The function relies on the comb_label function to combine label
      components for each hierarchical level.
    - Labels that don't have enough hierarchical depth are marked as "NA".
    - The returned nested list is organized by level first, then by cell.
    """
    # Each level's absolute label is a prefix of the same per-cell component
    # list, so compute that list once per cell (one `_found_parts` pass) and
    # build all `max_depth` cumulative prefixes from it. This replaces the
    # original num_cells x max_depth `comb_label` calls (each re-scanning and
    # re-splitting the row) and is output-identical: `comb_label(cell, i+1,
    # max_depth)` == join of the non-None components in levels 0..i, and the
    # "NA" rule (len(split("|")) < i+1) fires iff any of those levels has no
    # component (since components never contain "|").
    abs_labels_upto_level = [[] for _ in range(max_depth)]
    for cell_label in hierarchical_labels_array:
        parts = _found_parts(cell_label, max_depth, max_depth)
        prefix = []
        n_components = 0
        for i in range(max_depth):
            if parts[i] is not None:
                prefix.append(parts[i])
                n_components += 1
            combined_label = "|".join(prefix)
            split_len = 1 if combined_label == "" else n_components
            abs_labels_upto_level[i].append("NA" if split_len < i + 1 else combined_label)

    return abs_labels_upto_level


def split_labels_w_final_level(hierarchical_labels_array, max_depth):
    """Split hierarchical labels and extract final level annotations.

    This function processes an array of hierarchical labels for cells,
    where each cell's label is expressed as a series of hierarchical
    components separated by the '|' character. It performs the following
    steps:
    1. Uses the `abs_labels` function to generate absolute labels for
       each hierarchical level (from 1 to max_depth) for each cell.
    2. Constructs a 2D array of these labels (cells x levels) and
       determines the deepest level (final level) for each cell by
       counting non-'NA' entries.
    3. Appends the computed final level (1-indexed) as a new row to the
       list of absolute labels.
    4. For each hierarchical level, extracts only the final component
       (i.e., the substring after the last '|' delimiter).
    5. Extracts the final level label for each cell using the computed
       final level indices and appends this as an additional element.

    Parameters
    ----------
    hierarchical_labels_array : list
        A list (or iterable) where each element corresponds to a cell's
        hierarchical label. Each cell's label is itself a list (or
        similar iterable) of label strings for each level.
    max_depth : int
        The maximum number of hierarchical levels present in the
        predicted labels.

    Returns:
    -------
    tuple
        A tuple containing:
        - abs_labels_list : list
            A list where:
            * The first max_depth number of elements are lists of the
              final label components for each level (i.e., the part
              after the last '|' delimiter in the abs label up to the
              corresponding level).
            * The (max_depth+1)-th element is a list of the computed
              final levels (1-indexed) for each cell.
            * The (max_depth+2)-th element is a list of the final level
              labels for each cell.
        - final_levels_arr : numpy.ndarray
            A 1D numpy array of shape (num_cells,) containing the
            zero-indexed final level for each cell, computed as the
            number of non-'NA' levels minus one.

    Notes:
    -----
    - The function relies on the `abs_labels` function to generate
      absolute labels for each hierarchical level.
    - Final level determination is based on counting non-'NA' entries
      in the absolute labels.
    - The returned abs_labels_list contains processed label components,
      not the original absolute labels.
    """
    abs_labels_list = abs_labels(hierarchical_labels_array, max_depth)

    # Count non-"NA" levels and take the leaf of each abs label in a single
    # Python pass, instead of materialising two large (num_cells x levels)
    # string numpy arrays just to do a `!= "NA"` sum and a fancy-index. Output
    # is identical: level 0 is never "NA" (split_len >= 1 always), so the
    # non-NA count is >= 1 and final_levels_arr in [0, max_depth-1] indexes the
    # leaf levels exactly as the original .T fancy-index did.
    num_cells = len(abs_labels_list[0]) if max_depth else 0
    non_na_counts = [0] * num_cells
    for i in range(max_depth):
        col = abs_labels_list[i]
        leaf_col = [None] * num_cells
        for c in range(num_cells):
            v = col[c]
            if v != "NA":
                non_na_counts[c] += 1
            leaf_col[c] = v.split("|")[-1]
        abs_labels_list[i] = leaf_col

    final_levels_arr = np.array(non_na_counts, dtype=int) - 1
    final_levels_list = list(final_levels_arr + 1)
    final_level_labels = [abs_labels_list[final_levels_arr[c]][c] for c in range(num_cells)]

    abs_labels_list.append(final_levels_list)
    abs_labels_list.append(final_level_labels)

    return abs_labels_list, final_levels_arr


def categorize_refinement_type(value):
    """Categorize refinement type based on the input value.

    Parameters
    ----------
    value : any
        The value to be categorized.

    Returns:
    -------
    str
        The categorization result:
        - "Further": if value is a list
        - "No-match": if value is the string "False"
        - "Match": if value is any other string
        - "Unknown": for any other type
    """
    if isinstance(value, list):
        return "Further"
    elif value == "False":
        return "No-match"
    elif isinstance(value, str):
        return "Match"
    return "Unknown"


def compute_entropy(softmax):
    """Compute entropy from softmax outputs.

    Args:
        softmax: A softmax array

    Returns:
        entropy_vals: Array of entropy values
    """
    # Clip probabilities to avoid log(0)
    probs_clipped = np.clip(softmax, 1e-8, 1.0)

    # Compute entropy: H = -sum(p * log(p))
    entropy_vals = -np.sum(probs_clipped * np.log(probs_clipped), axis=1)

    return entropy_vals


def softmax_to_logits(softmax_probs, temperature=1.0):
    """Convert softmax probabilities back to logits.

    Args:
        softmax_probs: numpy array of softmax probabilities
        temperature: temperature parameter (default=1.0)

    Returns:
        logits: numpy array of logits
    """
    # Clip probabilities to avoid log(0)
    probs_clipped = np.clip(softmax_probs, 1e-8, 1.0)

    # Convert to logits using inverse softmax (log probabilities)
    logits = np.log(probs_clipped) * temperature

    # Normalize by subtracting the mean (softmax is invariant to constant shifts)
    # This helps with numerical stability and follows standard practice
    logits = logits - np.mean(logits, axis=1, keepdims=True)

    return logits


def get_total_number_of_classes() -> int:
    """Get the total number of Azimuth classification models available."""
    num_classes = 0
    for encoder_file in AZIMUTH_INFERENCE_ENCODERS_PATHS_SORTED_LIST:
        list_of_classes = np.loadtxt(encoder_file, dtype=str, delimiter=",")
        label_encoder = LabelEncoder().fit(list_of_classes)
        num_classes += len(label_encoder.classes_)
    return num_classes
