#!/usr/bin/env bash
# Build script for OpenCV.
#
# Run as
# build_opencv.sh \
#     version \
#     path/to/cmake \
#     path/to/ninja \
#     path/to/patchelf \
#     path/to/conda/bin/python \
#     path/to/opencv/CMakeLists.txt \
#     path/to/dest/cv2.so

VERSION="$1"
shift
CMAKE=$(realpath -sL "$1")
shift
NINJA=$(realpath -sL "$1")
shift
PATCHELF=$(realpath -sL "$1")
shift
PYTHON=$(realpath -sL "$1")
shift
SRC=$(realpath -sL "$1")
SRC=$(dirname "$SRC")
shift
DEST=$(realpath -smL "$1")
set -e

CONDA=$(dirname "$PYTHON")
CONDA=$(dirname "$CONDA")
NINJA_DIR=$(dirname "$NINJA")
export PATH="$CONDA/bin:$PATH:$NINJA_DIR"
if [ -n "$CC" ]; then
    CC_DIR=$(dirname "$CC")
    CC_DIR=$(realpath -sL "$CC_DIR")
    export PATH="$CC_DIR:$PATH"
    CXX=$(which c++)
    export CXX
fi
if [ -n "$AR" ]; then
    AR_DIR=$(dirname "$AR")
    export RANLIB="$AR_DIR/ranlib"
fi
export TBBROOT="$CONDA"

# Jan 1 2024.  We shouldn't have the build date embedded anywhere at all, but
# adding this just in case to improve build reproducibility.
export SOURCE_DATE_EPOCH=1704096000

DEST_DIR=$(dirname "$DEST")
mkdir -p "$DEST_DIR"
DIR=$(dirname "$DEST_DIR")
cd "$DIR"
mkdir -p build
rm -rf build/*

cd build

REAL_DIR=$(realpath "$DIR")
REAL_SRC=$(realpath "$SRC")
REAL_CONDA=$(realpath "$CONDA")

declare -a COMMON_OPTS=(
    "-mtune=skylake"
    "-march=core-avx-i"
    "-gno-record-gcc-switches"
    "-fcolor-diagnostics"
    "-g0"
    "-gline-tables-only"
    "-Wno-builtin-macro-redefined"
    "-DOPENCV_INSTALL_PREFIX=install"
    "-DOPENCV_BUILD_DIR=build"
    "-ffile-prefix-map=.=opencv"
    "-ffile-prefix-map=$DIR=opencv"
    "-ffile-prefix-map=${REAL_DIR}=opencv"
    "-ffile-prefix-map=$SRC=opencv"
    "-ffile-prefix-map=${REAL_SRC}=opencv"
    "-ffile-prefix-map=$CONDA=opencv"
    "-ffile-prefix-map=${REAL_CONDA}=opencv"
)

declare -a _C_FLAGS=(
    "${COMMON_OPTS[@]}"
    '-D__DATE__="redacted"'
    '-D__TIMESTAMP__="redacted"'
    '-D__TIME__="redacted"'
)

declare -a _SO_LDFLAGS=(
    "${COMMON_OPTS[@]}"
    "-L$CONDA/lib"
    "-Wl,--icf=safe,--gc-sections"
    "-Wl,-rpath='\$ORIGIN/../..'"
    "-Wl,--compress-debug-sections=zlib"
)

declare -a _CMAKE_FLAGS=(
    # General build setup
    -DCMAKE_BUILD_TYPE="Release"
    -DOPENCV_VCSVERSION="$VERSION"
    -DOPENCV_DOWNLOAD_PATH="$DIR/cache"
    -DCMAKE_INSTALL_PREFIX=/tmp
    -DCMAKE_INCLUDE_PATH="$CONDA/include"
    -DCMAKE_LIBRARY_PATH="$CONDA/lib"
    -DCMAKE_BUILD_RPATH="\$ORIGIN/../.."
    -DOPENCV_TIMESTAMP="redacted"
    -DCMAKE_COLOR_DIAGNOSTICS=ON
    -DDOPENCV_GENERATE_PKGCONFIG=OFF
    -DWITH_PROTOBUF=OFF
    # Dependencies
    -DPython_ROOT_DIR="$CONDA"
    -DPython_VERSION_MAJOR=3
    -DPYTHON_EXECUTABLE="$PYTHON"
    -DPYTHON2_EXECUTABLE="$PYTHON"
    -DPYTHON3_EXECUTABLE="$PYTHON"
    -DPYTHON3_LIBRARY="$CONDA/lib/libpython3.12.so"
    -DMKL_ROOT_DIR="$CONDA"
    -DWITH_TBB=ON
    -DMKL_WITH_TBB=ON
    -DOPENCV_LAPACK_DISABLE_MKL=OFF
    -DLAPACK_INCLUDE_DIR="$CONDA/include"
    -DMKL_INCLUDE_DIRS="$CONDA/include"
    -DMKL_USE_SINGLE_DYNAMIC_LIBRARY=ON
    -DEIGEN_INCLUDE_PATH="$CONDA/include/eigen3"
    -DATLAS_FOUND=OFF
    -DAtlas_LAPACK_LIBRARY=/dev/null
    -DAtlas_CBLAS_INCLUDE_DIR=/dev/null
    -DAtlas_INCLUDE_SEARCH_PATHS=/dev/null
    -DAtlas_LIB_SEARCH_PATHS=/dev/null
    -DOpenBLAS_FOUND=OFF
    -DOpenBLAS_INCLUDE_DIR=/dev/null
    -DOpen_BLAS_INCLUDE_SEARCH_PATHS=/dev/null
    -DOpen_BLAS_LIB_SEARCH_PATHS=/dev/null
    -DUSE_IPPICV=OFF
    -DBUILD_IPP=OFF
    -DWITH_IPP=OFF
    -DBUILD_IPP_IW=OFF
    -DOPENCV_PYTHON2_INSTALL_PATH=python
    -DOPENCV_PYTHON3_INSTALL_PATH=python
    -DBUILD_OPENEXR=OFF
    -DWITH_OPENEXR=OFF
    -DBUILD_JASPER=OFF
    -DWITH_JASPER=OFF
    -DBUILD_WEBP=OFF
    -DWITH_WEBP=OFF
    -DBUILD_QUIRC=OFF
    -DWITH_1394=OFF
    -DWITH_QUIRC=OFF
    -DWITH_WIN32UI=OFF
    -DWITH_QT=OFF
    -DWITH_GSTREAMER=OFF
    -DWITH_VTK=OFF
    -DWITH_OPENJPEG=OFF
    -DWITH_DSHOW=OFF
    -DWITH_MSMF=OFF
    -DWITH_DIRECTX=OFF
    -DWITH_MATLAB=OFF
    -DWITH_HDF5=OFF
    -DLAPACK_IMPL=MKL
    -DWITH_V4L=OFF
    -DWITH_FFMPEG=OFF
    -DBUILD_JPEG=OFF
    -DWITH_JPEG=ON
    -DWITH_GTK=OFF
    -DWITH_EIGEN=ON
    -DWITH_LAPACK=ON
    # Components to build
    -DOPENCV_MODULE_opencv_core_WRAPPERS=python
    -DOPENCV_MODULE_opencv_features2d_WRAPPERS=python
    -DOPENCV_MODULE_opencv_imgcodecs_WRAPPERS=python
    -DOPENCV_MODULE_opencv_imgproc_WRAPPERS=python
    -DOPENCV_PYTHON_PIP_METADATA_INSTALL=OFF
    -DBUILD_opencv_apps=OFF
    -DBUILD_opencv_calib3d=ON
    -DBUILD_opencv_core=ON
    -DBUILD_opencv_dnn=OFF
    -DBUILD_opencv_features2d=ON
    # TODO(azarchs): disable flann.  It shouldn't be required for features2d, but
    # it is, for the python bindings.
    -DBUILD_opencv_flann=ON
    -DBUILD_opencv_gapi=OFF
    -DBUILD_opencv_highgui=OFF
    -DBUILD_opencv_imgcodecs=ON
    -DBUILD_opencv_imgproc=ON
    -DBUILD_opencv_java_bindings_generator=OFF
    -DBUILD_opencv_java=OFF
    -DBUILD_opencv_js=OFF
    -DBUILD_opencv_js_bindings_generator=OFF
    -DBUILD_opencv_ml=OFF
    -DBUILD_opencv_objdetect=OFF
    -DBUILD_opencv_objc=OFF
    -DBUILD_opencv_objc_bindings_generator=OFF
    -DBUILD_opencv_photo=OFF
    -DBUILD_opencv_python_tests=OFF
    -DBUILD_opencv_python3=ON
    -DBUILD_opencv_stitching=OFF
    -DBUILD_opencv_ts=OFF
    -DBUILD_opencv_video=OFF
    -DBUILD_opencv_videoio=OFF
    -DBUILD_opencv_world=OFF
    -DBUILD_DOCS=OFF
    -DBUILD_EXAMPLES=OFF
    -DBUILD_PERF_TESTS=OFF
    -DBUILD_SHARED_LIBS=OFF
    -DBUILD_TESTS=OFF
    -DENABLE_FLAKE8=OFF
    -DENABLE_PYLINT=OFF
    -DINSTALL_C_EXAMPLES=OFF
    -DINSTALL_CREATE_DISTRIB=ON
    # Features to enable
    -DAPPLE_FRAMEWORK=OFF
    -DBUILD_CUDA_STUBS=OFF
    -DBUILD_JAVA=OFF
    -DCV_TRACE=OFF
    -DHIGHGUI_ENABLE_PLUGINS=OFF
    -DOPENCV_SKIP_PYTHON_LOADER=ON
    -DWITH_ADE=OFF
    -DWITH_CUDA=OFF
    -DWITH_IMGCODEC_HDR=OFF
    -DWITH_IMGCODEC_SUNRASTER=OFF
    -DWITH_OPENCL=OFF
    -DWITH_PROTOBUF=OFF
    -DWITH_VA=OFF
    -DWITH_VA_INTEL=OFF
    # Compilation options
    -DCMAKE_BUILD_RPATH_USE_ORIGIN=ON
    -DCMAKE_C_FLAGS="${_C_FLAGS[*]}"
    -DCMAKE_CXX_FLAGS="${_C_FLAGS[*]}"
    -DCMAKE_EXE_LINKER_FLAGS="${COMMON_OPTS[*]}"
    -DCMAKE_MODULE_LINKER_FLAGS="${_SO_LDFLAGS[*]}"
    -DCMAKE_SHARED_LINKER_FLAGS="${_SO_LDFLAGS[*]}"
    -DCPU_BASELINE_REQUIRE="SSE4_2;AVX"
    # CPU dispatch for avx2 changes results.
    -DCPU_DISPATCH="FP16"
    -DCV_DISABLE_OPTIMIZATION=OFF
    -DCMAKE_CXX_STANDARD=17
    -DCMAKE_CXX_EXTENSIONS=OFF
    -DENABLE_CXX11=ON
    -DENABLE_PRECOMPILED_HEADERS=OFF
    -DENABLE_THIN_LTO=ON
)

# Suppress warning about failure to generate python typing stubs.
export PYTHONWARNINGS="ignore::UserWarning"

"$CMAKE" --log-level=WARNING -Wno-dev -GNinja "$SRC" "${_CMAKE_FLAGS[@]}"

# Remove absolute paths from cmake output that is embedded in the binary.
ESCAPED_BUILD_DIR=$(echo "$DIR" | sed -e 's/[]"\/$*.^[]/\\&/g')
ESCAPED_CONDA_DIR=$(echo "$CONDA" | sed -e 's/[]"\/$*.^[]/\\&/g')
ESCAPED_SRC_DIR=$(echo "$SRC" | sed -e 's/[]"\/$*.^[]/\\&/g')
ESCAPED_NINJA_DIR=$(dirname "$NINJA" | sed -e 's/[]"\/$*.^[]/\\&/g')
sed -i "s/${ESCAPED_BUILD_DIR}/opencv/g
s/${ESCAPED_CONDA_DIR}/anaconda/g
s/${ESCAPED_SRC_DIR}/anaconda/g
s/${ESCAPED_NINJA_DIR}/toolchain\\/bin/g
s/\\bTimestamp:\\s*[0-9TZ:-]*/Timestamp:                   REDACTED/
s/\\bHost:\\s*Linux \S*/Host:                        Linux REDACTED/" modules/core/version_string.inc
"$NINJA" --quiet opencv_python3
SOLIB=$(find lib -name '*.so' | head -n 1)
"$PATCHELF" --set-rpath "\$ORIGIN/../.." --output "$DEST" "$SOLIB"
