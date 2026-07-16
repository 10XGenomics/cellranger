"""Paths associated with Azimuth."""

from pathlib import Path

AZIMUTH_CALIBRATION_MODELS_BASE_DIR = (
    Path(__file__).parent / "azimuth_models" / "latest" / "calibration"
)
AZIMUTH_CALIBRATION_MODELS_PREFIX = "temp_scaler_L_"
AZIMUTH_CALIBRATION_MODELS_SUFFIX = ".onnx"
AZIMUTH_CALIBRATION_MODELS_PATHS_SORTED_LIST = sorted(
    [
        y
        for y in AZIMUTH_CALIBRATION_MODELS_BASE_DIR.iterdir()
        if not y.is_dir()
        and y.name.startswith(AZIMUTH_CALIBRATION_MODELS_PREFIX)
        and y.name.endswith(AZIMUTH_CALIBRATION_MODELS_SUFFIX)
    ],
    key=lambda x: x.name,
)

AZIMUTH_INFERENCE_ENCODERS_BASE_DIR = (
    Path(__file__).parent / "azimuth_models" / "latest" / "inference_encoders"
)
AZIMUTH_INFERENCE_ENCODERS_PREFIX = "inference_encoders_level"
AZIMUTH_INFERENCE_ENCODERS_SUFFIX = ".txt"
AZIMUTH_INFERENCE_ENCODERS_PATHS_SORTED_LIST = sorted(
    [
        y
        for y in AZIMUTH_INFERENCE_ENCODERS_BASE_DIR.iterdir()
        if not y.is_dir()
        and y.name.startswith(AZIMUTH_INFERENCE_ENCODERS_PREFIX)
        and y.name.endswith(AZIMUTH_INFERENCE_ENCODERS_SUFFIX)
    ],
    key=lambda x: x.name,
)

AZIMUTH_INFERENCE_FEATURE_PANEL_BASE_DIR = (
    Path(__file__).parent / "azimuth_models" / "latest" / "inference_feature_panel"
)
AZIMUTH_INFERENCE_FEATURE_PANEL_NAME = "inference_feature_panel.txt"
AZIMUTH_INFERENCE_FEATURE_PANEL_PATH = (
    AZIMUTH_INFERENCE_FEATURE_PANEL_BASE_DIR / AZIMUTH_INFERENCE_FEATURE_PANEL_NAME
)

AZIMUTH_INFERENCE_MODEL_BASE_DIR = (
    Path(__file__).parent / "azimuth_models" / "latest" / "inference_model"
)
AZIMUTH_INFERENCE_MODEL_NAME = "inference_model.onnx"
AZIMUTH_INFERENCE_MODEL_PATH = AZIMUTH_INFERENCE_MODEL_BASE_DIR / AZIMUTH_INFERENCE_MODEL_NAME

AZIMUTH_POSTPROCESSING_BASE_DIR = (
    Path(__file__).parent / "azimuth_models" / "latest" / "postprocessing"
)
AZIMUTH_POSTPROCESSING_NAME_TEMPLATE = "panhuman_annotate_{}.csv"
