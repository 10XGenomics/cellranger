"""Code taken from Azimuth.

Module with inference tools using Azimuth Neural Network trained on
annotated panhuman scRNA-seq data.
"""

# pylint: disable = too-many-lines
import csv
import gc

import numpy as np
import onnxruntime as ort
from scipy.sparse import csr_matrix
from sklearn.preprocessing import LabelEncoder

from cellranger.cell_typing.azimuth.azimuth_models.latest.model_metadata import (
    model_meta as model_metadata,
)
from cellranger.cell_typing.azimuth.azimuth_ref_paths import (
    AZIMUTH_CALIBRATION_MODELS_PATHS_SORTED_LIST,
    AZIMUTH_INFERENCE_ENCODERS_PATHS_SORTED_LIST,
    AZIMUTH_INFERENCE_FEATURE_PANEL_PATH,
    AZIMUTH_INFERENCE_MODEL_PATH,
    AZIMUTH_POSTPROCESSING_BASE_DIR,
    AZIMUTH_POSTPROCESSING_NAME_TEMPLATE,
)
from cellranger.cell_typing.azimuth.azimuth_utils import (
    categorize_refinement_type,
    comb_label,
    compute_entropy,
    if_full_consistent_hierarchy,
    normalize,
    reorder_subset_data_matrix,
    softmax_to_logits,
    split_labels_w_final_level,
)

VALID_REFINE_LEVELS = ["broad", "medium", "fine"]
DEFAULT_CALIBRATION_METHOD = "temperature_scaling_entropy_informed"

########################################################################
################### class objects ######################################
########################################################################
########################################################################


class MemoryContext:
    """A context manager for memory-intensive operations.

    This context manager provides a simple interface to encapsulate
    operations that are memory-intensive. Upon exiting the context, it
    triggers garbage collection to help free up memory.

    Attributes:
    ----------
    description : str
        A description of the memory-intensive operation.
    """

    def __init__(self, description="Memory-intensive operation"):
        """Initialize a MemoryContext instance.

        Parameters
        ----------
        description : str, optional
            A description of the memory-intensive operation.
            Default is "Memory-intensive operation".
        """
        self.description = description

    def __enter__(self):
        """Enter the memory context.

        Returns:
        -------
        MemoryContext
            The current MemoryContext instance.
        """
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Exit the memory context and perform garbage collection.

        Parameters
        ----------
        exc_type : type
            The exception type, if an exception occurred, otherwise None.
        exc_val : Exception
            The exception instance, if an exception occurred, otherwise
            None.
        exc_tb : traceback
            The traceback object, if an exception occurred, otherwise
            None.

        Returns:
        -------
        None
        """
        gc.collect()


class Inference:
    """Run a supervised inference workflow using a confidence-calibrated model.

    Public Interface
    ----------------

    Methods:
      run_inference()
          Execute the complete inference workflow on the input matrix and
          decode predictions to human-readable labels. Wrapped around
          the method run_on_X()
      run_on_X()
          Run inference on the entire input matrix in evaluation batches.

    Private
    --------------------------

    Attributes:
      _X : numpy.ndarray or scipy.sparse matrix
          The input expression matrix of shape (n_samples, n_features).
      _model : keras.Model
          The trained neural network model used for making predictions.
      _label_encoders : list
          List of label encoder objects (one per hierarchical level)
          used for decoding model predictions.
      _eval_batch_size : int
          Batch size used during inference.
      _eval_steps : int
          Number of complete evaluation batches (minus one), computed as
          X.shape[0] // eval_batch_size.
      _max_depth : int
          Maximum number of hierarchical levels to process.

    Private Methods
      run_on_minibatch(minibatch)
          Run model inference on a single minibatch.
      process_minibatch(
                        minibatch,
                        softmax_mb_levels_cache,
                        class_preds_mb_cache,
                        max_probs_mb_cache
                    )
          Process a minibatch to cache softmax outputs, class
          predictions, and maximum probabilities for each hierarchical
          level.
    """

    def __init__(self, X, model, label_encoders, eval_batch_size, max_depth):
        """Initialize the Inference object with input data, model, and evaluation parameters.

        Parameters
        ----------
        X : numpy.ndarray or scipy.sparse matrix
            The input expression matrix of shape (n_samples, n_features).
            Can be sparse for memory efficiency.
        model : keras.Model
            The trained neural network model to be used for inference.
            Should output a list of tensors, one per hierarchical level.
        label_encoders : list
            A list of label encoder objects, one per hierarchical level,
            used for decoding numerical predictions to string labels.
            The length of this list should match max_depth.
        eval_batch_size : int
            The batch size to use during evaluation. Controls memory usage
            during inference.
        max_depth : int
            The maximum number of hierarchical levels to process in the
            classification taxonomy.
        """
        self._X = X  # pylint: disable=invalid-name
        self._model = model
        self._label_encoders = label_encoders
        self._eval_batch_size = eval_batch_size
        self._eval_steps = X.shape[0] // self._eval_batch_size
        self._max_depth = max_depth

    def run_on_minibatch(self, minibatch):
        """[Developer-Facing] Run model inference on a single minibatch.

        Parameters
        ----------
        minibatch : numpy.ndarray or a tensor
            A dense array representing a subset of the input matrix.

        Returns:
        -------
        list of numpy.ndarray/tensor
            The model's predictions (softmax probabilities) for
            the minibatch.
        """
        y_mb = self._model.run(None, {"input_layer": minibatch.astype(np.float32)})
        return y_mb

    def process_minibatch(
        self, minibatch, softmax_mb_levels_cache, class_preds_mb_cache, max_probs_mb_cache
    ):
        """Process a minibatch by running inference and caching predictions.

        This method runs inference on the provided minibatch and, for
        each hierarchical level, caches the softmax outputs, computes
        the predicted class indices, and determines the maximum
        probability values.

        Parameters
        ----------
        minibatch : numpy.ndarray or tensor
            A dense array corresponding to a batch of input data.
        softmax_mb_levels_cache : list of lists
            Cache for storing softmax outputs for each hierarchical level.
        class_preds_mb_cache : list
            Cache for storing predicted class indices for each minibatch.
        max_probs_mb_cache : list
            Cache for storing maximum probability values for each
            minibatch.

        Returns:
        -------
        tuple
            Updated caches:
              - softmax_mb_levels_cache (list of lists)
              - class_preds_mb_cache (list)
              - max_probs_mb_cache (list)

        Notes:
        -----
        This method uses TensorFlow operations (tf) for processing
        predictions.
        """
        class_levels_cache = []
        prob_levels_cache = []

        # This function must now return a list of NumPy arrays
        y_mb = self.run_on_minibatch(minibatch)

        for i in range(self._max_depth):
            # No change here, just appending the NumPy array to a list
            softmax_mb_levels_cache[i].append(y_mb[i])

            # Find the index of the highest probability prediction
            # TF: tf.argmax(y_mb[i], axis=-1) -> NP: np.argmax(...)
            class_indices = np.argmax(y_mb[i], axis=-1)

            # Add a new dimension to make it a column vector for concatenation
            # TF: tf.expand_dims(tf.cast(..., tf.int32), axis=1) -> NP: np.expand_dims(...)
            class_outs = np.expand_dims(class_indices.astype(np.int32), axis=1)
            class_levels_cache.append(class_outs)

            # Find the highest probability value itself
            # TF: tf.reduce_max(..., keepdims=True) -> NP: np.max(...)
            max_probs = np.max(y_mb[i], axis=-1, keepdims=True)
            prob_levels_cache.append(max_probs)

        # Concatenate the results from all levels into single arrays
        # TF: tf.concat(...) -> NP: np.concatenate(...)
        class_preds_mb = np.concatenate(class_levels_cache, axis=1)
        max_probs_mb = np.concatenate(prob_levels_cache, axis=1)

        # Append the results for the entire minibatch to the main caches
        class_preds_mb_cache.append(class_preds_mb)
        max_probs_mb_cache.append(max_probs_mb)

        return (softmax_mb_levels_cache, class_preds_mb_cache, max_probs_mb_cache)

    def run_on_X(self):  # pylint: disable=invalid-name
        """Run inference on the entire input matrix in evaluation batches.

        The input matrix is divided into batches based on
        _eval_batch_size. For each batch, the model is run to obtain
        predictions which are then aggregated. Consistency checks
        are performed to ensure that the number of predictions matches
        the input dimensions.

        Returns:
        -------
        tuple
            A tuple containing:
              - class_preds : numpy.ndarray
                  Array of predicted class indices with shape
                  (n_samples, _max_depth).
              - max_probs : numpy.ndarray
                  Array of maximum probability values at each level
                  with shape (n_samples, _max_depth).
              - softmax_vals_all : list of numpy.ndarray
                  List of softmax outputs for each hierarchical level.

        Raises:
        ------
        ValueError
            If the input matrix is empty.
        RuntimeError
            If no predictions were generated.
        AssertionError
            If the dimensions of the predictions do not match the input.
        """
        if self._X.shape[0] == 0:
            raise ValueError("Input matrix is empty")

        softmax_mb_levels_cache = [[] for _ in range(self._max_depth)]
        class_preds_mb_cache = []
        max_probs_mb_cache = []

        print(f"Splitting query data into {self._eval_steps+1} " "evaluation batches.\n")

        print("Running model:")

        for i in range(self._eval_steps):
            start_idx = i * self._eval_batch_size
            end_idx = (i + 1) * self._eval_batch_size
            if start_idx >= self._X.shape[0]:
                break

            minibatch = self._X[start_idx:end_idx].toarray()
            if minibatch.size == 0:
                continue

            (softmax_mb_levels_cache, class_preds_mb_cache, max_probs_mb_cache) = (
                self.process_minibatch(
                    minibatch, softmax_mb_levels_cache, class_preds_mb_cache, max_probs_mb_cache
                )
            )
            del minibatch

        final_start = self._eval_steps * self._eval_batch_size
        if final_start < self._X.shape[0]:
            final_batch = self._X[final_start:].toarray()
            if final_batch.size > 0:
                (softmax_mb_levels_cache, class_preds_mb_cache, max_probs_mb_cache) = (
                    self.process_minibatch(
                        final_batch,
                        softmax_mb_levels_cache,
                        class_preds_mb_cache,
                        max_probs_mb_cache,
                    )
                )
                del final_batch

        if not class_preds_mb_cache:
            raise RuntimeError("No predictions were generated")

        softmax_vals_all = [np.concatenate(level, axis=0) for level in softmax_mb_levels_cache]
        class_preds = np.concatenate(class_preds_mb_cache, axis=0)
        max_probs = np.concatenate(max_probs_mb_cache, axis=0)

        assert self._max_depth == class_preds.shape[1], (
            "The array of all predictions should have dimensions consistent " "with max_depth."
        )
        assert (
            self._X.shape[0] == class_preds.shape[0]
        ), "The number of predictions != the number of cells in X_query"
        for level in softmax_vals_all:
            assert (
                level.shape[0] == self._X.shape[0]
            ), "Number of cells in softmax_outputs != numbers of cells in query "

        return (np.array(class_preds), np.array(max_probs), softmax_vals_all)

    def run_inference(self):
        """Run the complete inference workflow and decode predictions to string labels.

        This public method executes the full inference process on the
        input matrix (in batches), gathers predictions from the model,
        and decodes the class predictions using the provided
        label encoders. The final output is a dictionary containing
        both the raw numerical predictions and human-readable labels.

        Returns:
        -------
        dict
            Dictionary with keys:
              'hierarchical_label_preds' : numpy.ndarray
                  Array of decoded string labels with shape
                  (n_samples, _max_depth).
              'class_preds' : numpy.ndarray
                  Array of predicted class indices with shape
                  (n_samples, _max_depth).
              'probability_of_preds' : numpy.ndarray
                  Array of maximum probability values with shape
                  (n_samples, _max_depth).
              'softmax_vals_all' : list of numpy.ndarray
                  List of softmax outputs for each hierarchical level.

        Notes:
        -----
        This method uses a MemoryContext context manager to handle memory
        efficient processing, which should be imported from an appropriate
        utility module.
        """
        string_labels_out = []

        with MemoryContext():
            class_preds, max_probs, softmax_vals_all = self.run_on_X()
        for i in range(self._max_depth):
            string_labels_out.append(self._label_encoders[i].inverse_transform(class_preds[:, i]))

        return {
            "hierarchical_label_preds": np.array(string_labels_out).T,
            "class_preds": class_preds,
            "probability_of_preds": max_probs,
            "softmax_vals_all": softmax_vals_all,
        }


class CalibrationSingleClassifier:
    """Calibrate softmax outputs. Can accommodate different calibration methods.

    This class provides functionality to apply calibration techniques to
    softmax probability outputs from neural network classifiers. It supports
    dynamic method selection for different calibration approaches.
    Currently, it hosts a entropy-informed temperature scaling method for
    calibration.

    The calibration process operates on batches of softmax outputs to
    manage memory usage efficiently, making it suitable for large datasets.

    Parameters
    ----------
    softmax : numpy.ndarray or scipy.sparse matrix
        Softmax probability outputs from a neural network model, with
        shape (n_samples, n_classes).
    eval_batch_size : int
        Batch size for processing softmax outputs. Controls memory usage
        during calibration.

    Attributes:
    ----------
    _softmax : numpy.ndarray or scipy.sparse matrix
        Original softmax probability outputs.
    _calibrated_softmax : numpy.ndarray
        Calibrated softmax outputs after applying calibration method.
        Initially set to original softmax values.
    _eval_batch_size : int
        Batch size used for processing.
    _eval_steps : int
        Number of complete evaluation batches, computed as
        softmax.shape[0] // eval_batch_size.

    Examples:
    --------
    >>> import numpy as np
    >>> from tensorflow.keras.models import load_model
    >>>
    >>> # Get softmax outputs from a Keras model
    >>> main_model = load_model('classifier_model.keras')
    >>> softmax_outputs = main_model.predict(X_data)
    >>> # Shape: (n_samples, n_classes)
    >>>
    >>> # Load a trained temperature scaling model
    >>> calibrator = load_model('temperature_scaling_model.keras')
    >>>
    >>> # Initialize calibration object
    >>> calibration = CalibrationSingleClassifier(
    ...     softmax=softmax_outputs,
    ...     eval_batch_size=256
    ... )
    >>>
    >>> # Apply entropy-informed temperature scaling
    >>> calibrated_probs = calibration.calibrate(
    ...     'temperature_scaling_entropy_informed',
    ...     calibrator
    ... )

    Notes:
    -----
    The class uses a dynamic method calling approach via the `calibrate()`
    method, allowing for flexible selection of calibration techniques.
    New calibration methods can be added by implementing them as class
    methods and calling them through `calibrate()`.

    The entropy-informed temperature scaling method computes entropy from
    the softmax outputs and uses it along with converted logits as input
    to the calibration model.
    """

    def __init__(self, softmax, eval_batch_size):
        """Initialize the CalibrationSingleClassifier.

        Parameters
        ----------
        softmax : numpy.ndarray or scipy.sparse matrix
            Softmax probability outputs with shape (n_samples, n_classes).
        eval_batch_size : int
            Batch size for processing softmax outputs.
        """
        self._softmax = softmax
        self._calibrated_softmax = softmax
        self._eval_batch_size = eval_batch_size
        self._eval_steps = softmax.shape[0] // self._eval_batch_size

    def calibrate(self, calibrator):
        """Call calibration using temperature_scaling_entropy_informed."""
        return self.temperature_scaling_entropy_informed(calibrator=calibrator)

    def temperature_scaling_entropy_informed_minibatch(self, softmax_minibatch, calibrator):
        """Apply entropy-informed temperature scaling to a single minibatch.

        This method processes a single batch of softmax outputs by:
        1. Computing entropy from the softmax probabilities
        2. Converting softmax probabilities back to logits
        3. Applying the calibration model using both logits and entropy

        Parameters
        ----------
        softmax_minibatch : numpy.ndarray
            A batch of softmax probability outputs with shape
            (batch_size, n_classes).
        calibrator : callable
            Calibration model or function that takes (logits, entropy)
            as input and returns calibrated probabilities.

        Returns:
        -------
        numpy.ndarray
            Calibrated softmax probabilities for the minibatch with
            the same shape as input.

        Notes:
        -----
        This method expects the calibrator to accept a tuple of
        (logits, entropy) as input. The entropy is computed using
        the Shannon entropy formula: H = -sum(p * log(p)).
        """
        ent = compute_entropy(softmax_minibatch)
        ent = np.array(ent, dtype=np.float32).reshape(-1, 1)

        logits = softmax_to_logits(softmax_minibatch)
        # out_minibatch = calibrator((logits, ent))
        input_name_logits = calibrator.get_inputs()[0].name
        input_name_entropy = calibrator.get_inputs()[1].name
        output_name = calibrator.get_outputs()[0].name
        out_minibatch = calibrator.run(
            [output_name], {input_name_logits: logits, input_name_entropy: ent}
        )[0]

        return out_minibatch

    def temperature_scaling_entropy_informed(self, calibrator):
        """Apply entropy-informed temperature scaling to all softmax outputs.

        This method processes the complete set of softmax outputs in
        batches, applying entropy-informed temperature scaling to each
        batch and concatenating the results. It handles the final
        incomplete batch automatically.

        The calibration process:
        1. Divides softmax outputs into batches based on eval_batch_size
        2. For each batch, computes entropy and converts to logits
        3. Applies the calibration model to get temperature-scaled outputs
        4. Concatenates all calibrated batches into final result

        Parameters
        ----------
        calibrator : callable
            Calibration model or function that takes (logits, entropy)
            as input and returns calibrated probabilities.

        Returns:
        -------
        numpy.ndarray
            Calibrated softmax probabilities for all samples with shape
            (n_samples, n_classes). Also updates the internal
            _calibrated_softmax attribute.

        Notes:
        -----
        This method modifies the internal state by updating
        _calibrated_softmax with the calibrated results. It handles
        both sparse matrices and numpy arrays by checking for the
        .toarray() method.

        The method prints progress information including the number of
        batches being processed.
        """
        print(f"Splitting softmax outputs into {self._eval_steps+1} " "batches.\n")

        calibrated_softmax_cache = []

        for i in range(self._eval_steps):
            start_idx = i * self._eval_batch_size
            end_idx = (i + 1) * self._eval_batch_size
            if start_idx >= self._softmax.shape[0]:
                break

            # Handle both sparse matrices and numpy arrays
            minibatch_source = self._softmax[start_idx:end_idx]
            if hasattr(minibatch_source, "toarray"):
                minibatch = minibatch_source.toarray()
            else:
                minibatch = minibatch_source

            if minibatch.size == 0:
                continue

            calibrated_softmax_cache.append(
                self.temperature_scaling_entropy_informed_minibatch(minibatch, calibrator)
            )

        final_start = self._eval_steps * self._eval_batch_size
        if final_start < self._softmax.shape[0]:
            final_batch_source = self._softmax[final_start:]
            if hasattr(final_batch_source, "toarray"):
                final_batch = final_batch_source.toarray()
            else:
                final_batch = final_batch_source

            if final_batch.size > 0:
                calibrated_softmax_cache.append(
                    self.temperature_scaling_entropy_informed_minibatch(final_batch, calibrator)
                )

        self._calibrated_softmax = np.concatenate(calibrated_softmax_cache, axis=0)

        return self._calibrated_softmax


class InferenceTools:
    """Load and manage the necessary components required for model inference.

    This class handles loading of model artifacts from predefined locations,
    including the trained model, label encoders, feature panel, and model
    metadata. It supports different annotation pipelines and performs
    validation of loaded components.

    Public Interface
    ----------------

    Methods:
      load_inference_model()
          Load the trained Keras model for inference.
      load_model_meta()
          Load and validate model metadata containing configuration
          parameters.
      load_inference_encoders()
          Load label encoders used for decoding model predictions.
      load_inference_feature_panel()
          Load the feature panel for subsetting query data (supervised
          pipeline only).

    Private Attributes
    -----------------
      _inference_model_filename : str
          Filename of the Keras model file to load.
      _model_meta : dict
          Dictionary containing model metadata and configuration
          parameters.
      _inference_encoders_file_basename : str
          Filename of the pickled label encoders.
      _inference_feature_panel_filename : str
          Filename of the text file containing feature names.
      _annotation_pipeline : str
          Type of annotation pipeline to use (e.g., 'supervised').
    """

    def __init__(
        self,
        annotation_pipeline,
        inference_model_filename="inference_model.onnx",
        model_meta=None,
        inference_encoders_file_basename="inference_encoders_level",
        inference_feature_panel_filename="inference_feature_panel.txt",
    ):
        """Initialize the InferenceTools with file paths and configuration.

        Parameters
        ----------
        annotation_pipeline : str
            Type of annotation pipeline to use (e.g., 'supervised',
            'self-supervised', currently the only pipeline implemented
            is 'supervised'). Determines which components are required
            for inference.
        inference_model_filename : str, optional
            Filename of the Keras model file to load, by default
            'inference_model.keras'.
        model_meta : dict
            Dictionary containing model metadata and configuration
            parameters. Must include specific keys depending on the
            annotation pipeline.
        inference_encoders_file_basename : str, optional
            Filename of the pickled label encoders, by default
            'inference_encoders.pkl'.
        inference_feature_panel_filename : str, optional
            Filename of the text file containing feature names,
            by default 'inference_feature_panel.txt'.
        """
        self._inference_model_filename = inference_model_filename
        if model_meta is not None:
            self._model_meta = model_meta
        else:
            self._model_meta = model_metadata
        self._inference_encoders_file_basename = inference_encoders_file_basename
        self._inference_feature_panel_filename = inference_feature_panel_filename

        self._annotation_pipeline = annotation_pipeline

    def load_inference_model(self):
        """Load the trained ONNX model for inference.

        Loads the model from a predefined directory structure using the
        filename specified during initialization. The model is expected
        to be in the ONNX format.

        Returns:
        -------
        onnxruntime.InferenceSession
            The loaded ONNX model session used for inference.

        Notes:
        -----
        The model must be saved in the 'inference_model' directory
        accessible via the 'files' import system. Model name must
        correspond to the name provided to the object at initialization.
        """
        session = ort.InferenceSession(AZIMUTH_INFERENCE_MODEL_PATH)
        return session

    def load_model_meta(self):
        """Load and validate model metadata containing configuration parameters.

        Validates that the model metadata dictionary contains all
        required keys based on the annotation pipeline type. Checks that
         values have the expected types and enforces constraints on the
         metadata structure.

        The model metadata dictionary is defined in the __init__ file
        for the package directory and is imported from there onto this
        module.

        Returns:
        -------
        dict
            Dictionary containing validated model metadata and
            configuration parameters.

        Raises:
        ------
        AssertionError
            If any required key is missing or if values have incorrect
            types.
        """
        meta_dict = self._model_meta

        for meta_key in meta_dict.keys():
            assert isinstance(
                meta_key, str
            ), "keys in the model metadata dict must be provided as strings."

        assert "inference_model_name" in meta_dict.keys(), (
            "model metadata dict must have a " "key 'inference_model_name'"
        )

        assert "max_depth" in meta_dict.keys(), "model metadata dict must have a key 'max_depth'."

        assert isinstance(meta_dict["max_depth"], int), "max_depth must be an integer."

        assert "inference_model_embedding_layer" in meta_dict.keys(), (
            "model metadata dict must have a key " "'inference_model_embedding_layer'"
        )

        assert isinstance(
            meta_dict["inference_model_embedding_layer"], str
        ), "inference_model_embedding_layer name must be a string."

        if self._annotation_pipeline == "supervised":

            assert (
                "feature_panel_size" in meta_dict.keys()
            ), "model metadata dict must have a key 'feature_panel_size'."

            assert isinstance(
                meta_dict["feature_panel_size"], int
            ), "feature_panel_size must be an integer."
        # enforce existence of certain key:vals here

        return meta_dict

    def load_inference_encoders(self):
        """Load label encoders used for decoding model predictions.

        Loads pickled label encoder objects from a predefined directory
        structure using the filename specified during initialization.
        These encoders are used to convert numeric predictions back to
        string labels.

        Returns:
        -------
        list
            List of label encoders, typically one per hierarchical level
            in the classification taxonomy.

        Notes:
        -----
        The encoders must be saved as a pickle file in the
        'inference_encoders' directory accessible via the 'files' import
         system, and the name of the pickled file must correspond to the
         name provided to this object at initialization.
        """
        encoders = []
        for encoder_file in AZIMUTH_INFERENCE_ENCODERS_PATHS_SORTED_LIST:
            list_of_classes = np.loadtxt(encoder_file, dtype=str, delimiter=",")
            label_encoder = LabelEncoder().fit(list_of_classes)
            encoders.append(label_encoder)

        return encoders

    def load_inference_feature_panel(self):
        """Load the feature panel for subsetting query data.

        For supervised annotation pipelines, loads a list of feature names
        (typically gene names) from a text file. These features are used
        to subset the input data before inference. Validates that the
        feature list matches expected size from metadata and contains no
         duplicates.

        Returns:
        -------
        list or None
            List of feature names (strings) if using a supervised pipeline,
            or None for other pipeline types.

        Raises:
        ------
        AssertionError
            If duplicates are found in the feature list or if the number
            of features doesn't match the value in model metadata.

        Notes:
        -----
        The feature panel must be saved as a text file in the
        'inference_feature_panel' directory accessible via the 'files'
        import system, with the filename matching the name passed to this
        object at initialization.

        For non-supervised pipelines, this method returns None.
        """
        if self._annotation_pipeline == "supervised":
            with open(AZIMUTH_INFERENCE_FEATURE_PANEL_PATH) as f:
                feat_list = f.read().splitlines()

            feat_list = [
                gene.decode("utf-8") if isinstance(gene, bytes) else gene for gene in feat_list
            ]
            feat_list = sorted(feat_list)

            assert len(feat_list) == len(
                set(feat_list)
            ), "There are possible duplicates in the inference feature panel."

            if "feature_panel_size" in self.load_model_meta().keys():
                assert len(feat_list) == (self.load_model_meta()["feature_panel_size"]), (
                    "Number of features in inference feature panel provided \n"
                    " does not match with the number provided in inference \n"
                    "model metadata."
                )

            return feat_list
        else:
            return None

    def load_calibrators(self):
        """Load the trained Onnx models for calibration.

        Loads the models from a predefined directory structure using the
        version and the calibration filenames pre-defined for that version.
        The models are expected to be in .onnx format. If available, custom
        objects required for model loading are imported from the version's
        calibration.custom_objects module.

        Returns:
        -------
        list of onnx.Model
            List of loaded onnx models with length = max_depth, if
            calibrators are available.
        """
        calibrators = [
            ort.InferenceSession(str(path)) for path in AZIMUTH_CALIBRATION_MODELS_PATHS_SORTED_LIST
        ]
        return calibrators


class AutoloadInferenceTools(InferenceTools):
    """Enhanced version of InferenceTools that automatically loads all components on initialization.

    This subclass extends InferenceTools by automatically calling all
    loading methods during initialization and storing their results as
    instance attributes. This provides a more convenient interface where
     components can be accessed directly as attributes rather than
     calling loading methods each time they are needed.

    For example, instead of calling `tools.load_inference_model()` each
    time, the model can be accessed directly via `tools.inference_model`.

    Attributes:
    ----------
    inference_model : keras.Model
        The loaded model for inference (from load_inference_model).
    model_meta : dict
        Dictionary of model metadata (from load_model_meta).
    inference_encoders : list
        List of label encoders (from load_inference_encoders).
    inference_feature_panel : list or None
        List of features for subsetting data, or None for non-supervised
        pipelines (from load_inference_feature_panel).

    Notes:
    -----
    All methods that start with 'load_' will be automatically called
    during initialization, and their return values will be stored as
    attributes. The attribute names are derived by removing the 'load_'
    prefix from the method name.
    """

    def __init__(self, annotation_pipeline, *args, **kwargs):
        super().__init__(annotation_pipeline, *args, **kwargs)

        self.inference_encoders = self.load_inference_encoders()
        self.inference_feature_panel = self.load_inference_feature_panel()
        self.inference_model = self.load_inference_model()
        self.calibrators = self.load_calibrators()
        self.model_meta = self.load_model_meta()


class InferenceInputData:
    """Processes a counts data matrix into the standardized format required for the inference model.

    This class handles normalization and feature alignment of single-cell
    expression data, preparing it for downstream inference tasks. It
    validates input data types, performs normalization, and can align
    query features with a reference feature panel template.

    Parameters
    ----------
    X_query : scipy.sparse.csr_matrix
        Expression matrix in CSR format with cells as rows and features
        (genes) as columns.
    query_features : list of str
        List of feature names corresponding to the columns in X_query.
    feature_panel_template : list of str or None
        Optional reference feature list used for standardizing the
        feature space. If provided, the query data will be reordered and
         subset to match this template.
    normalization_override : bool, default=False
        If True, bypasses normalization entirely regardless of whether
        the values are integers or not. When False, normalization will
        be applied based on the data characteristics.
    norm_check_batch_size : int, default=1000
        Batch size used for checking normalization status by checking
        whether the values in this batch are integers or not.

    Attributes:
    ----------
    X_query : scipy.sparse.csr_matrix
        Normalized expression matrix.
    query_features : list of str
        List of feature names from the query data.
    feature_panel_template : list of str or None
        Reference feature list for standardization.
    num_cells : int
        Number of cells in the query data.
    template_size : int
        Number of features in the reference template (if provided).
    common_features : set
        Set of features common to both query and template.
    overlap_percent : float
        Percentage of template features present in the query data.

    Raises:
    ------
    AssertionError
        If inputs don't match expected types or if there are duplicates
        in the feature panel template.
    UserWarning
        If feature overlap between query and template is below 20%.
    """

    def __init__(
        self,
        X_query,  # pylint: disable=invalid-name
        query_features,
        feature_panel_template,
        normalization_override=False,
        norm_check_batch_size=1000,
    ):

        assert isinstance(
            X_query, csr_matrix
        ), "X_query should be in the scipy.sparse.csr_matrix format."
        self.X_query, _ = normalize(  # pylint: disable=invalid-name
            X_query,
            normalization_override=normalization_override,
            norm_check_batch_size=norm_check_batch_size,
        )

        self.num_cells = X_query.shape[0]

        assert isinstance(query_features, list), "query_features must be a list."
        assert all(
            isinstance(item, str) for item in query_features
        ), "Each element of the list of query_features must be a string."
        self.query_features = query_features

        if feature_panel_template:
            assert isinstance(
                feature_panel_template, list
            ), "feature panel template must be a list."
            assert all(isinstance(item, str) for item in feature_panel_template), (
                "Each element of the list of feature panel " "template must be a string."
            )
        self.feature_panel_template = feature_panel_template
        if self.feature_panel_template:
            self.template_size = len(set(self.feature_panel_template))
            assert self.template_size == len(
                self.feature_panel_template
            ), "There are possible duplicates in the feature_panel_template."

            self.common_features = set(self.query_features).intersection(
                set(self.feature_panel_template)
            )

            self.overlap_percent = (len(self.common_features) * 100) / self.template_size

            if self.overlap_percent < 20:
                print(
                    f"Low feature overlap ({self.overlap_percent:.1f}%) may "
                    "lead to suboptimal results.",
                )

    def reorder_subset_on_feature_template(self):
        """Reorders and subsets the query data to match the feature panel template.

        Aligns the feature dimensions of the query data with the
        reference feature panel. Features present in the template but
        missing in the query will be filled with zeros. Features present
         in the query but not in the template will be excluded.

        Returns:
        -------
        scipy.sparse.csr_matrix
            Expression matrix reordered padded and subset to match the
            feature panel template. If no template is provided, returns
            the original expression matrix.

        Warns:
        -----
        UserWarning
            If no feature panel template was provided.

        Notes:
        -----
        This method prints a summary of the query data dimensions and
        feature overlap with the reference panel when successful.
        """
        if self.feature_panel_template:
            reordered_X_query = reorder_subset_data_matrix(  # pylint: disable=invalid-name
                self.X_query,
                self.query_features,
                self.feature_panel_template,
                common_features=self.common_features,
            )

            assert reordered_X_query.shape[1] == self.template_size, (
                "Feature dimension mismatch after reordering " "and subsetting the query."
            )

            print("Query object:")
            print(f"    Total number of cells: {self.num_cells}")
            print("    Total number of features: " f"{len(set(self.query_features))}")
            print(
                "    Overlap w. feature reference: "
                f"{len(self.common_features)} "
                f"(~{int(self.overlap_percent)}%)\n"
            )

            return reordered_X_query
        else:
            print(
                "A reference feature panel has not been provided, returning "
                "query data on all query features."
            )
            return self.X_query

    def inference_input(self, annotation_pipeline="supervised"):
        """Produces the final input matrix for passing to the inference model.

        Based on the specified annotation pipeline, this method produces
        the appropriately formatted input data matrix for the inference
        model. Currently only 'supervised' pipeline is supported, which
        reorders pads and subsets the query data according to the
        feature panel template.

        Parameters
        ----------
        annotation_pipeline : str, default='supervised'
            The type of annotation pipeline to use. Currently only
            'supervised' is supported.

        Returns:
        -------
        scipy.sparse.csr_matrix
            Expression matrix ready for inference.

        Raises:
        ------
        ValueError
            If an unsupported annotation pipeline is specified.
        """
        if annotation_pipeline == "supervised":
            input_matrix = self.reorder_subset_on_feature_template()
        else:
            raise ValueError("annotation_pipeline not recognized.")

        return input_matrix


class OutputLabels:
    """Processes and organizes hierarchical labels produced by inference models.

    This class handles the post-processing of prediction outputs from
    hierarchical classification models. It extracts and organizes
    predictions at different levels of the hierarchy, extracts
    probabilities, and provides methods to access level-specific
    information.

    Parameters
    ----------
    labels_pred : numpy.ndarray
        Predicted hierarchical labels, typically with shape
        (num_cells, max_depth).
    labels_prob : numpy.ndarray
        Probability scores associated with each prediction.
    max_depth : int
        Maximum depth of the hierarchical classification.
    num_cells : int
        Number of cells/samples in the query data.

    Attributes:
    ----------
    combined_labels : list
        Combined hierarchical labels for each cell.
    level_zero_labels : numpy.ndarray
        Labels at the top level (level 0) of the hierarchy.
    final_level_labels : numpy.ndarray
        Labels at the final level determined for each cell.
    level_zero_softmax_prob : numpy.ndarray
        Softmax probabilities for the top level predictions.
    final_level_softmax_prob : numpy.ndarray
        Softmax probabilities for the final level predictions.
    full_consistent_hierarchy : list
        Flags indicating whether each cell has a fully consistent
        hierarchy.

    Raises:
    ------
    AssertionError
        If the number of predicted labels doesn't match the number of
        cells.
    """

    def __init__(self, labels_pred, labels_prob, max_depth, num_cells):

        assert (
            labels_pred.shape[0] == num_cells
        ), "Number of labels predicted != number of cells in query."

        self._labels_predicted = labels_pred
        self._labels_prob = labels_prob
        self._max_depth = max_depth
        self._num_cells = num_cells

        # The combined label, the abs/leaf decomposition and the consistency
        # flag are all pure functions of a cell's predicted-label row. Many
        # cells share identical rows (typically <20% are unique), so compute
        # these on the unique rows only and expand back by inverse index. This
        # is output-identical to the original per-cell computation; the lookup
        # preserves per-cell ordering.
        inverse = [0] * num_cells
        uniq_rows = []
        _uniq_index = {}
        for c in range(num_cells):
            key = tuple(self._labels_predicted[c])
            j = _uniq_index.get(key)
            if j is None:
                j = len(uniq_rows)
                _uniq_index[key] = j
                uniq_rows.append(self._labels_predicted[c])
            inverse[c] = j
        uniq_rows_arr = np.array(uniq_rows)

        uniq_combined = [comb_label(r, self._max_depth, self._max_depth) for r in uniq_rows_arr]
        self.combined_labels = [uniq_combined[inverse[c]] for c in range(num_cells)]

        uniq_level_specific, uniq_final_levels = split_labels_w_final_level(
            uniq_rows_arr, self._max_depth
        )
        self._level_specific_labels_and_final_level = [
            [col[inverse[c]] for c in range(num_cells)] for col in uniq_level_specific
        ]
        self._final_levels_array = uniq_final_levels[np.asarray(inverse)]

        self.level_zero_labels = self._level_specific_labels_and_final_level[0]
        self.final_level_labels = self._level_specific_labels_and_final_level[self._max_depth + 1]

        self.level_zero_softmax_prob = self._labels_prob[:, 0]
        self.final_level_softmax_prob = self._labels_prob[
            np.arange(self._num_cells), self._final_levels_array
        ]

        uniq_consistency = [if_full_consistent_hierarchy(r, self._max_depth) for r in uniq_rows_arr]
        self.full_consistent_hierarchy = [uniq_consistency[inverse[c]] for c in range(num_cells)]

    def level_specific_softmax_prob(self, level):
        """Get softmax probabilities for a specific hierarchical level.

        Parameters
        ----------
        level : int
            Level to get probabilities for (1-indexed).

        Returns:
        -------
        numpy.ndarray
            Softmax probabilities for the specified level.

        Raises:
        ------
        TypeError
            If level is not an integer.
        ValueError
            If level is out of bounds.
        """
        if not isinstance(level, int):
            raise TypeError("level must be an integer")
        if level < 1 or level > self._max_depth:
            raise ValueError(f"level must be between 1 and {self._max_depth}")

        softmax_probs = self._labels_prob[:, level - 1]
        return softmax_probs

    def level_specific_labels(self, level):
        """Get predicted labels for a specific hierarchical level.

        Parameters
        ----------
        level : int
            Level to get labels for (1-indexed).

        Returns:
        -------
        numpy.ndarray
            Predicted labels for the specified level.

        Raises:
        ------
        IndexError
            If level is out of bounds.
        """
        labels = self._level_specific_labels_and_final_level[level - 1]
        return labels

    def all_level_labels(self):
        """Get predicted labels for all hierarchical levels.

        Returns:
        -------
        list of numpy.ndarray
            List containing the predicted labels for each level in the
            hierarchy, organized by level.

        Notes:
        -----
        This method provides a convenient way to access all hierarchical
         levels at once, which can be useful for comprehensive analysis
         or visualization of the hierarchical classification results.
        """
        all_level_labels = [self.level_specific_labels(level) for level in range(self._max_depth)]

        return all_level_labels


class PostprocessingAzimuthLabels:  # pylint: disable=too-few-public-methods
    """Specialized class for refining and postprocessing Azimuth labels.

    This class refines cell type annotations at different levels of granularity
     (broad, medium, fine). It applies post-processing rules to improve
     annotation quality based on hierarchical consistency, reference
     annotations, and confidence scores.

    Parameters
    ----------
    combined_labels : list of str
        Per-cell combined hierarchical label (the `full_hierarchical_labels`
        produced by `OutputLabels`/`process_outputs`). Passed in so the
        expensive base computation is done once and reused across refine
        levels rather than recomputed per level.
    full_consistent_hierarchy : list of bool
        Per-cell flag indicating whether the predicted hierarchy is
        internally consistent (also from `OutputLabels`).
    num_cells : int
        Number of cells/samples in the query data.
    softmax_probs : list of numpy.ndarray
        Softmax probability arrays for each level of the hierarchy.
    encoders : list
        List of label encoders used during model training/inference.
    refine_level : str
        Level of refinement to apply ('broad', 'medium', or 'fine').

    Attributes:
    ----------
    softmax_probs : list of numpy.ndarray
        Softmax probability arrays for each level.
    refine_level : str
        Selected level of refinement.
    encoders : list
        Label encoders for each hierarchical level.
    refined_annotations_df_path : str
        Path to reference annotation mapping data.
    refined_annotations_file_name : str
        Filename of the reference annotation file.

    Raises:
    ------
    ValueError
        If refine_level is not one of the valid options.
    RuntimeError
        If the postprocessing module cannot be found.
    FileNotFoundError
        If the annotation file is not found.
    ValueError
        If the annotation file is missing required columns.

    Notes:
    -----
    This class uses several private methods (prefixed with double
    underscore) for internal processing. Users should primarily interact
     with the public `refine_labels()` method.
    """

    def __init__(
        self,
        combined_labels,
        full_consistent_hierarchy,
        num_cells,
        softmax_probs,
        encoders,
        refine_level,
    ):

        if refine_level not in VALID_REFINE_LEVELS:
            raise ValueError(f"refine_level must be one of {VALID_REFINE_LEVELS}")

        # Reuse the already-computed base outputs instead of re-running the
        # expensive OutputLabels.__init__ (combined-label assembly +
        # split_labels_w_final_level), none of which the refine path needs
        # beyond these two per-cell arrays.
        self.combined_labels = combined_labels
        self.full_consistent_hierarchy = full_consistent_hierarchy
        self._num_cells = num_cells

        self.softmax_probs = softmax_probs
        self.refine_level = refine_level
        self.encoders = encoders

        # Memoizes __map_to_annotate_refine by its (consistency-flagged) label
        # string. The mapping is a pure function of that string + annot_dict,
        # and distinct hierarchy paths are far fewer than cells, so this
        # collapses the per-cell full-dict prefix scan to once per unique path.
        self._map_cache: dict = {}

        self.refined_annotations_file_name = AZIMUTH_POSTPROCESSING_NAME_TEMPLATE.format(
            refine_level
        )
        self.refined_annotations_df_path = (
            AZIMUTH_POSTPROCESSING_BASE_DIR / self.refined_annotations_file_name
        )

        try:
            with open(self.refined_annotations_df_path) as fl:
                csv_reader = csv.DictReader(fl)
                columns_in_csv = csv_reader.fieldnames
        except FileNotFoundError as e:
            raise FileNotFoundError(
                f"Annotation file {self.refined_annotations_file_name} not "
                f"found in postprocessing directory"
            ) from e

        required_columns = ["Orig_Label", f"Annotate_{refine_level}"]
        missing_cols = [col for col in required_columns if col not in columns_in_csv]
        if missing_cols:
            raise ValueError(f"Annotation file missing required columns: {missing_cols}")
        self.annot_dict = self.__annotations_dict()

    def __combined_labels_w_consistency(self):
        """Private method: Creates labels with consistency flag.

        Returns:
        -------
        list of str
            Combined labels with appended consistency flag.
        """
        labels_w_consistency_flag = [
            f"{self.combined_labels[i]}_{self.full_consistent_hierarchy[i]}"
            for i in range(self._num_cells)
        ]

        return labels_w_consistency_flag

    def __encoders_dict(self):
        """Private method: Creates a dictionary mapping of label encoders.

        Returns:
        -------
        dict
            Dictionary mapping level numbers to label encodings.
        """
        enc_dicts = {}
        for i, encoder in enumerate(self.encoders):
            level_name = int(i)
            enc_dicts[level_name] = {label: index for index, label in enumerate(encoder.classes_)}

        return enc_dicts

    def __annotations_dict(self):
        """Private method: Creates a dictionary of annotation mappings.

        Returns:
        -------
        dict
            Dictionary mapping original labels to refined annotations.
        """
        with open(self.refined_annotations_df_path) as fl:
            csv_reader = csv.DictReader(fl)
            annot_dict = dict(
                (
                    row["Orig_Label"],
                    row[f"Annotate_{self.refine_level}"],
                )
                for row in csv_reader
            )
        with open(self.refined_annotations_df_path) as fl:
            csv_reader = csv.DictReader(fl)
            for row in csv_reader:
                annot_dict[row[f"Annotate_{self.refine_level}"]] = row[
                    f"Annotate_{self.refine_level}"
                ]

        return annot_dict

    def __softmax_arrays_dict(self):
        """Private method: Creates a dictionary of softmax probability arrays.

        Returns:
        -------
        dict
            Dictionary mapping array identifiers to softmax arrays.
        """
        softmax = {f"arr_{i}": np.array(row) for i, row in enumerate(self.softmax_probs)}

        return softmax

    def __map_to_annotate_refine(self, label):
        """Private method: Maps a label to its refined annotation (memoized).

        Parameters
        ----------
        label : str
            Original label with consistency flag.

        Returns:
        -------
        str or list
            Refined label or list of potential refined labels.
        """
        cached = self._map_cache.get(label)
        if cached is not None:
            return cached
        result = self.__compute_map_to_annotate_refine(label)
        self._map_cache[label] = result
        return result

    def __compute_map_to_annotate_refine(self, label):
        """Uncached body of :meth:`__map_to_annotate_refine`."""
        label_split = label.split("_")
        flag = label_split[-1]
        label = "_".join(label_split[:-1])

        if flag == "True":
            refined_label = self.annot_dict.get(label)
            if refined_label is None:
                matching_refined_labels = [
                    self.annot_dict[orig_label]
                    for orig_label in self.annot_dict
                    if orig_label.startswith(label)
                ]
                matching_refined_labels = list(set(matching_refined_labels))
                if len(matching_refined_labels) == 1:
                    return matching_refined_labels[0]
                return matching_refined_labels
            if refined_label == label:
                return refined_label
            elif refined_label in label or all(
                p in label.split("|") for p in refined_label.split("|")
            ):
                return refined_label
        return "False"

    def __softmax_classes_dict(self):
        """Private method: Creates a dictionary of softmax probabilities by class.

        Returns:
        -------
        dict
            Nested dictionary mapping hierarchy levels and classes to
            probabilities.
        """
        precomputed_softmax = {}

        # Build the per-level softmax arrays once; the original rebuilt them
        # on every iteration of the loop below (max_depth re-materialisations).
        softmax = self.__softmax_arrays_dict()
        levels_classes_dict = self.__encoders_dict()
        for level, level_dict in levels_classes_dict.items():
            arr_key = f"arr_{level}"
            if arr_key in softmax:
                softmax_array = softmax[arr_key]
                precomputed_softmax[level] = {
                    key: softmax_array[:, idx]
                    for key, idx in level_dict.items()
                    if idx < softmax_array.shape[1]
                }

        return precomputed_softmax

    def refine_labels(self):
        """Refines raw predicted labels using reference annotations and confidence scores.

        This method applies a multi-step refinement process to improve
        annotation quality:
        1. Identifies labels with consistent hierarchies
        2. Maps consistent labels to reference annotations
        3. For ambiguous cases, selects the annotation with highest
            confidence
        4. Returns the most specific (leaf) annotation level

        Returns:
        -------
        list of str
            Refined cell type labels at the specified refinement level.

        Notes:
        -----
        The refinement process prioritizes hierarchical consistency and
        reference annotations while using confidence scores to resolve
        ambiguities.
        """
        labels_w_consistency_flag = self.__combined_labels_w_consistency()
        labels_options = list(map(self.__map_to_annotate_refine, labels_w_consistency_flag))

        refined_labels = []

        # Derive refine types from the already-computed labels_options instead
        # of re-running __combined_labels_w_consistency + the mapping inside
        # __refinement_type (which recomputed both from scratch).
        refine_types = list(map(categorize_refinement_type, labels_options))
        precomputed_softmax = self.__softmax_classes_dict()

        for i in range(self._num_cells):
            if refine_types[i] == "Further":
                annotations = labels_options[i]
                max_value = -np.inf
                best_prediction = None

                for pred in annotations:
                    levels = len(pred.split("|")) - 1

                    if pred in precomputed_softmax[levels]:
                        value = precomputed_softmax[levels][pred][i]
                        if value > max_value:
                            max_value = value
                            best_prediction = pred

            else:
                best_prediction = labels_options[i]

            refined_label = (
                best_prediction.split("|")[-1]
                if isinstance(best_prediction, str)
                else best_prediction
            )

            refined_labels.append(refined_label)

        return refined_labels
