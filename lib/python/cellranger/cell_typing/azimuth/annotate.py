"""Code taken from Azimuth panhuman.

Azimuth Cell Annotation: Neural network-based hierarchical cell type
annotation for single-cell RNA-seq data.

This module provides tools for hierarchical cell type annotation and
interpretation based on single-cell RNA-seq data using the Azimuth
neural network model trained on annotated panhuman scRNA-seq data.

Key Components
-------------
AzimuthNN : Class
    Low-level class providing fine-grained control over the annotation
    process. Suitable for advanced users who need detailed control or
    are processing data in batches to optimize memory usage.

annotate_core : Function
    Core function for script-based automated annotation. Designed for
    batch processing and integration into analysis pipelines.

Usage Examples
-------------
Example usage:
    >>> from cellranger.cell_typing.azimuth.annotate import AzimuthNN
    >>> matrix_path = "/path/to/feature/barcode/matrix.h5"
    >>> azimuth = AzimuthNN()
    >>> azimuth.query_feature_barcode_matrix_h5(matrix_path=matrix_path)
    >>> azimuth.process_query() # get the inputs in format that the model can be run
    >>> _ = azimuth.run_inference_model() # run the model
    >>> azimuth.calibrate_predictions()
    >>> _ = azimuth.process_outputs() # postprocess model outputs
    >>> # (Optional) refinement of labels for consistent granularity and ease of interpretation
    >>> _ = azimuth.azimuth_refine()
    >>> azimuth.write_results_to_csv("/path/to/annotations.csv")  # save results

For more detailed documentation on specific classes and functions:
    >>> help(AzimuthNN)
    >>> help(annotate_core)
"""

import numpy as np

import cellranger.matrix as cr_matrix
import cellranger.rna.library as rna_library
from cellranger.cell_typing.azimuth.annotate_tools import (
    DEFAULT_CALIBRATION_METHOD,
    AutoloadInferenceTools,
    CalibrationSingleClassifier,
    Inference,
    InferenceInputData,
    OutputLabels,
    PostprocessingAzimuthLabels,
)
from cellranger.cell_typing.azimuth.azimuth_models.latest.model_metadata import (
    model_meta as model_metadata,
)
from cellranger.cell_typing.common_cell_typing import (
    AZIMUTH_SUMMARY_SCORE_KEY,
    BARCODE_KEY,
    BROAD_CELL_TYPES_KEY,
    COARSE_CELL_TYPES_KEY,
    FINE_CELL_TYPES_KEY,
    FULL_CELL_TYPES_KEY,
    SUMMARY_SCORE_KEY,
)

BARCODE_AZIMUTH_KEY = "barcode"
FULL_CELL_TYPES_AZIMUTH_KEY = "full_hierarchical_labels"
FINAL_CELL_TYPES_AZIMUTH_KEY = "final_level_labels"
FINE_CELL_TYPES_AZIMUTH_KEY = "azimuth_fine"
MEDIUM_CELL_TYPES_AZIMUTH_KEY = "azimuth_medium"
BROAD_CELL_TYPES_AZIMUTH_KEY = "azimuth_broad"
SOFTMAX_PROB_AZIMUTH_KEY = "final_level_softmax_prob"
FALSE_CELL_TYPE_AZIMUTH_KEY = "false"

AZIMUTH_DATA_NAME_TO_CR_DATA_NAME_DICT = {
    BARCODE_AZIMUTH_KEY: BARCODE_KEY,
    BROAD_CELL_TYPES_AZIMUTH_KEY: BROAD_CELL_TYPES_KEY,
    MEDIUM_CELL_TYPES_AZIMUTH_KEY: COARSE_CELL_TYPES_KEY,
    FINE_CELL_TYPES_AZIMUTH_KEY: FINE_CELL_TYPES_KEY,
    FULL_CELL_TYPES_AZIMUTH_KEY: FULL_CELL_TYPES_KEY,
    SOFTMAX_PROB_AZIMUTH_KEY: AZIMUTH_SUMMARY_SCORE_KEY,
}
IN_ORDER_COLUMN_NAMES = [
    BARCODE_KEY,
    BROAD_CELL_TYPES_KEY,
    COARSE_CELL_TYPES_KEY,
    FINE_CELL_TYPES_KEY,
    FULL_CELL_TYPES_KEY,
    SUMMARY_SCORE_KEY,
]


########################################################################
###### Base class for low level interactive usage ######################


class AzimuthNN(AutoloadInferenceTools):  # pylint: disable=too-many-instance-attributes
    """Base class for low-level interactive usage of the Azimuth neural network annotation pipeline.

    This class provides a comprehensive framework for single-cell
    RNA-seq annotation using neural network models. It handles the
    complete workflow from data loading and preprocessing to inference,
    post-processing, and result visualization. This includes functionality
    for refining annotations at different levels of granularity.

    Parameters
    ----------
    annotation_pipeline : str, default='supervised'
        The type of annotation pipeline to use.
    eval_batch_size : int, default=8192
        Batch size for inference.

    Attributes:
    ----------
    X_query : scipy.sparse.csr_matrix or None
        Expression matrix in CSR format.
    query_features : list or None
        List of feature names.
    cells_meta : pandas.DataFrame or None
        Cell metadata.
    num_cells : int or None
        Number of cells in the query.
    processed_outputs : dict or None
        Processed inference results.

    Raises:
    ------
    TypeError
        If input parameters are not of the correct type.
    RuntimeError
        If model metadata fails to load.

    Notes:
    -----
    This class is designed for programmatic use and provides
    fine-grained control over each step of the annotation pipeline.
    Consider using a higher-level interface for convenience if a
    standard workflow is sufficient for your needs.
    """

    def __init__(self, annotation_pipeline="supervised", eval_batch_size=8192):

        if not isinstance(annotation_pipeline, str):
            raise TypeError("annotation_pipeline must be a string")

        if not isinstance(eval_batch_size, int):
            raise TypeError("eval_batch_size must be an integer")

        self._annotation_pipeline = annotation_pipeline
        self._eval_batch_size = eval_batch_size

        super().__init__(annotation_pipeline)

        if not hasattr(self, "model_meta"):
            raise RuntimeError("Failed to load model metadata")

        self.inference_model_name = self.model_meta["inference_model_name"]
        self.max_depth = self.model_meta["max_depth"]
        self.inference_model_embedding_layer = self.model_meta["inference_model_embedding_layer"]
        self.feature_panel_size = self.model_meta["feature_panel_size"]
        self.calibration_method = self.model_meta["calibration"]
        self.azimuth_version = self.model_meta["version"]

        self.X_query = None  # pylint: disable=invalid-name
        self.query_features = None
        self.cells_meta = None
        self.num_cells = None
        self._inference_input_matrix = None

        self._inference_outputs_unprocessed = None
        self.processed_outputs = None
        self._azimuth_refined_labels = {}

    @staticmethod
    def get_azimuth_version() -> str:
        return model_metadata.get("version", "unknown")

    def query_feature_barcode_matrix_h5(self, matrix_path):
        """Load query data from feature-barcode matrix H5."""
        cmatrix = cr_matrix.CountMatrix.load_h5_file(matrix_path)
        cmatrix = cmatrix.select_features_by_type(rna_library.GENE_EXPRESSION_LIBRARY_TYPE)
        cmatrix = cmatrix.select_nonzero_axis(1)[0]  # Removing barcodes with zero UMIs
        self.X_query = (
            cmatrix.m.T
        )  # needs to be barcode x feature matrix. We store feature x barcode
        self.query_features = cmatrix.feature_ref.get_feature_names()
        self.cells_meta = [x.decode() for x in cmatrix.bcs]
        self.num_cells = self.X_query.shape[0]

    def process_query(self, normalization_override=False, norm_check_batch_size=100):
        """Process the query data to prepare it for inference.

        This method prepares the expression data for the inference model
        according to the specified annotation pipeline. The processing steps
        vary depending on the pipeline type, potentially including
        normalization, feature selection, dimensionality reduction, or
        other transformations.

        Parameters
        ----------
        normalization_override : bool, default=False
            If True, bypasses normalization entirely regardless of
            whether the values are integers or not.
        norm_check_batch_size : int, default=1000
            Batch size for checking normalization status.

        Raises:
        ------
        TypeError
            If parameters are not of the correct type.

        Notes:
        -----
        This method must be called after loading query data and before
        running inference. The specific processing
        steps depend on the annotation_pipeline specified during
        initialization.

        Currently, only the 'supervised' annotation pipeline is implemented,
        which normalizes the expression data and aligns it with a reference
        feature panel.
        """
        if not isinstance(normalization_override, bool):
            raise TypeError("normalization override must be a bool")

        if not isinstance(norm_check_batch_size, int):
            raise TypeError("norm_check_batch_size must be an integer")

        query_processing_class = InferenceInputData(
            self.X_query,
            self.query_features,
            self.inference_feature_panel,
            normalization_override=normalization_override,
            norm_check_batch_size=norm_check_batch_size,
        )

        self._inference_input_matrix = query_processing_class.inference_input(
            annotation_pipeline=self._annotation_pipeline
        )

    def run_inference_model(self):
        """Run the inference model on the processed query data.

        This method executes the neural network inference to generate
        cell type predictions.

        Returns:
        -------
        dict
            Dictionary of raw inference outputs including hierarchical
            label predictions and probabilities.

        Raises:
        ------
        AssertionError
            If input matrix has not been initialized by calling
            process_query().

        Notes:
        -----
        The raw outputs should typically be processed using the
        process_outputs() method before further use downstream.
        """
        assert (
            self._inference_input_matrix is not None
        ), "Input matrix not initialized. Call process_query() first."

        inference_class = Inference(
            self._inference_input_matrix,
            self.inference_model,
            self.inference_encoders,
            self._eval_batch_size,
            self.max_depth,
        )

        self._inference_outputs_unprocessed = inference_class.run_inference()

        return self._inference_outputs_unprocessed

    def calibrate_predictions(self):
        """Apply calibration to softmax outputs using trained calibrators.

        This method calibrates the softmax probability outputs from each
        hierarchical level using the corresponding trained calibration models,
        if available. Calibration improves the reliability and trustworthiness
        of prediction confidence scores by adjusting for overconfidence or
        underconfidence in the original model outputs.

        Returns:
        -------
        dict
            Updated inference outputs dictionary with calibrated results.
            Contains the same keys as the original inference outputs but with
            calibrated values:

            - 'softmax_vals_all': List of calibrated softmax probability arrays,
            one per hierarchical level
            - 'probability_of_preds': Updated maximum probability values from
            the calibrated softmax distributions
            - Other keys remain unchanged from the original inference outputs

        Notes:
        -----
        - Only applies calibration if calibration method is specified in model
        metadata and calibrator models are available
        - If no calibration is configured (calibration method is None), the
        method returns the original inference outputs unchanged
        - Each hierarchical level is calibrated independently using its own
        trained calibrator model
        - Memory management is applied during processing to handle large datasets
        efficiently by cleaning up intermediate results after each level

        Raises:
        ------
        AssertionError
            If calibration method is not None but the number of available
            calibrator models doesn't match the expected number of hierarchical
            levels (max_depth).

        Examples:
        --------
        The method is typically called as part of the inference pipeline:

        >>> azimuth = AzimuthNN_base()
        >>> # ... load data and run inference ...
        >>> raw_outputs = azimuth.run_inference_model()
        >>> calibrated_outputs = azimuth.calibrate_predictions()
        >>> # Calibrated outputs now have adjusted confidence scores

        The calibration process transforms prediction confidence scores:

        - Before calibration: Model might be over(/under)-confident
        (high probabilities for uncertain predictions or vice-versa)
        - After calibration: Probabilities better reflect true prediction
        confidence and uncertainty
        """
        if self.calibration_method is not None:
            assert len(self.calibrators) == self.max_depth, (
                "Calibration method is not None, expected number of "
                f"calibrator models: {self.max_depth}, number of "
                f"calibrator models found: {len(self.calibrators)}"
            )
            assert self.calibration_method == DEFAULT_CALIBRATION_METHOD, (
                f"Only {DEFAULT_CALIBRATION_METHOD} supported for calibration. "
                f"{self.calibration_method} requested."
            )

            softmax_all = self._inference_outputs_unprocessed["softmax_vals_all"]
            calibrated_levels_cache = []
            max_probs_levels_cache = []

            for level in range(self.max_depth):
                sm_array = softmax_all[level]
                calibrator_model = self.calibrators[level]

                calibration_obj = CalibrationSingleClassifier(
                    softmax=sm_array, eval_batch_size=self._eval_batch_size
                )

                calibrated_sm = calibration_obj.calibrate(calibrator_model)

                calibrated_levels_cache.append(calibrated_sm)
                max_probs_levels_cache.append(np.max(calibrated_sm, axis=-1))

            max_probs = np.column_stack(max_probs_levels_cache)
            self._inference_outputs_unprocessed["softmax_vals_all"] = calibrated_levels_cache
            self._inference_outputs_unprocessed["probability_of_preds"] = max_probs

        return self._inference_outputs_unprocessed

    def process_outputs(self, mode="minimal"):
        """Process raw inference outputs into usable predictions.

        This method organizes the raw inference outputs into a structured
        dictionary of predictions at various hierarchical levels.

        Parameters
        ----------
        mode : str, default='minimal'
            Processing mode: 'minimal' provides essential outputs,
            'detailed' includes additional information for all levels.

        Returns:
        -------
        dict
            Dictionary of processed outputs including hierarchical labels,
            level-specific labels, and confidence scores.

        Raises:
        ------
        AssertionError
            If mode is not 'minimal' or 'detailed'.

        Notes:
        -----
        This method should be called after run_inference_model().
        """
        assert mode in ["minimal", "detailed"], (
            "mode for output processing should be either " "'minimal' or 'detailed'"
        )

        labels_pred = self._inference_outputs_unprocessed["hierarchical_label_preds"]
        labels_prob = self._inference_outputs_unprocessed["probability_of_preds"]

        output_processing_class = OutputLabels(
            labels_pred, labels_prob, self.max_depth, self.num_cells
        )

        combined_labels = output_processing_class.combined_labels
        level_zero_labels = output_processing_class.level_zero_labels
        final_level_labels = output_processing_class.final_level_labels
        final_level_softmax_prob = output_processing_class.final_level_softmax_prob
        full_consistent_hierarchy = output_processing_class.full_consistent_hierarchy

        self.processed_outputs = {
            "full_hierarchical_labels": combined_labels,
            "level_zero_labels": level_zero_labels,
            "final_level_labels": final_level_labels,
            "final_level_softmax_prob": final_level_softmax_prob,
            "full_consistent_hierarchy": full_consistent_hierarchy,
        }

        if mode == "detailed":
            for i in range(1, self.max_depth):
                self.processed_outputs[f"level_{i}_labels"] = (
                    output_processing_class.all_level_labels()[i]
                )

        return self.processed_outputs

    def refine_labels(self, refine_level):
        """Refine hierarchical labels to a consistent level of granularity.

        This method applies post-processing rules to standardize
        annotations at the specified level of granularity (broad,
        medium, or fine).

        Parameters:
        ----------
        refine_level : str
            Level of refinement: 'broad', 'medium', or 'fine'.

        Returns:
        -------
        list
            List of refined labels at the specified level.

        Raises:
        ------
        AssertionError
            If refine_level is not valid or inference hasn't been run.

        Notes:
        -----
        For 'broad' level, this returns the top-level annotations.
        For 'medium' and 'fine' levels, specialized refinement is
        applied.
        """
        assert refine_level in [
            "broad",
            "medium",
            "fine",
        ], "refine_level should be 'broad', 'medium', or 'fine'."

        assert (
            self._inference_outputs_unprocessed is not None
        ), "Labels can be refined only after inference model has been run."

        print(
            "Interpreting label predictions for consistent granularity "
            f"at {refine_level} level.\n"
        )

        # All refine levels reuse the base outputs (combined labels + the
        # consistency flag) computed once by process_outputs().
        if self.processed_outputs is None:
            self.process_outputs()

        if refine_level in ["medium", "fine"]:

            softmax_probs = self._inference_outputs_unprocessed["softmax_vals_all"]

            refine_class = PostprocessingAzimuthLabels(
                self.processed_outputs[FULL_CELL_TYPES_AZIMUTH_KEY],
                self.processed_outputs["full_consistent_hierarchy"],
                self.num_cells,
                softmax_probs,
                self.inference_encoders,
                refine_level,
            )

            results = refine_class.refine_labels()

        else:
            results = self.processed_outputs["level_zero_labels"]

        (self._azimuth_refined_labels[f"azimuth_{refine_level}"]) = results

        return results

    def azimuth_refine(self):
        """Refine cell type annotations at multiple granularity levels.

        This method applies a hierarchical refinement of cell type labels,
        progressing from broad to fine classifications. The results are
        stored in the annotations attribute and the cell metadata is
        updated.

        Returns:
        -------
        None
            Updates the annotations attribute in-place and updates cell
            metadata.
        """
        _ = self.refine_labels(refine_level="broad")
        _ = self.refine_labels(refine_level="medium")
        _ = self.refine_labels(refine_level="fine")

    def write_results_to_csv(self, output_csv_path, fill_unannotated=False):
        """Write inference results to a CSV.

        Args:
            output_csv_path (str): Path to save the CSV.
            fill_unannotated (bool): If True, replaces missing values with
                "Coarse Not Annotated" (or "Fine..."). If False, falls back
                to the parent label (Broad -> Medium -> Fine).
        """
        if self.processed_outputs is None:
            raise RuntimeError(
                "Processed outputs not found. Please call process_outputs() before writing results."
            )

        # If any value is False, None replace with value above it OR specific string
        labels = self._azimuth_refined_labels
        azimuth_broad = np.array(labels[BROAD_CELL_TYPES_AZIMUTH_KEY], dtype=object)
        azimuth_medium = np.array(labels[MEDIUM_CELL_TYPES_AZIMUTH_KEY], dtype=object)
        azimuth_fine = np.array(labels[FINE_CELL_TYPES_AZIMUTH_KEY], dtype=object)

        medium_series = azimuth_medium.copy()
        fine_series = azimuth_fine.copy()
        broad_series = azimuth_broad.copy()

        # 1. Handle Coarse/Medium Column
        mask_medium = np.isin(
            medium_series, np.array([None, False, "False", "false"], dtype=object)
        )
        if fill_unannotated:
            medium_series[mask_medium] = "Coarse Not Confidently Annotated"
        else:
            medium_series[mask_medium] = broad_series[mask_medium]

        # 2. Handle Fine Column
        mask_fine = np.isin(fine_series, np.array([None, False, "False", "false"], dtype=object))
        if fill_unannotated:
            # Using "Fine Not Annotated" to distinguish it from the Coarse column
            fine_series[mask_fine] = "Fine Not Confidently Annotated"
        else:
            # Falls back to the (now potentially updated) medium_series
            fine_series[mask_fine] = medium_series[mask_fine]

        ## Update the dict
        labels[MEDIUM_CELL_TYPES_AZIMUTH_KEY] = medium_series.tolist()
        labels[FINE_CELL_TYPES_AZIMUTH_KEY] = fine_series.tolist()

        broad_col = labels[BROAD_CELL_TYPES_AZIMUTH_KEY]
        coarse_col = labels[MEDIUM_CELL_TYPES_AZIMUTH_KEY]
        fine_col = labels[FINE_CELL_TYPES_AZIMUTH_KEY]
        full_labels = self.processed_outputs[FULL_CELL_TYPES_AZIMUTH_KEY]

        # Get the softmax probabilities for the final level predictions
        probs = np.array(self.processed_outputs[f"{SOFTMAX_PROB_AZIMUTH_KEY}"])

        # Combine all columns into a single structured matrix for writing
        if self.cells_meta is None:
            raise ValueError(
                "cells_meta (barcode_col) is None. Please ensure query_feature_barcode_matrix_h5() has been called and completed successfully."
            )
        barcode_col = np.array(self.cells_meta, dtype=object)
        combined_data = np.column_stack(
            (barcode_col, broad_col, coarse_col, fine_col, full_labels, probs.astype(str))
        )

        # Write CSV
        header = f"{BARCODE_KEY},{BROAD_CELL_TYPES_KEY},{COARSE_CELL_TYPES_KEY},{FINE_CELL_TYPES_KEY},{FULL_CELL_TYPES_AZIMUTH_KEY},{AZIMUTH_SUMMARY_SCORE_KEY}"
        np.savetxt(
            output_csv_path, combined_data, delimiter=",", header=header, fmt="%s", comments=""
        )
