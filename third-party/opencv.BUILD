load("@conda_package_python//:vars.bzl", "PYTHON_PREFIX")
load(
    "@rules_license//rules:license.bzl",
    "license",
)
load(
    "@tenx_bazel_rules//rules/conda:conda_manifest.bzl",
    "conda_deps",
    "conda_manifest",
)

package(default_applicable_licenses = ["license"])

license(
    name = "license",
    package_name = "opencv",
    additional_info = {
        "homepage": "https://opencv.org/",
        "version": "4.5.4",
        "manifest": "third-party/deps.bzl",
        "pURL": "pkg:github/opencv/opencv@4.5.4",
    },
    copyright_notice = "Copyright (C) 2000-2015, Intel Corporation, " +
                       "2009-2011, Willow Garage Inc., " +
                       "2015, Itseez Inc., " +
                       "2015, OpenCV Foundation",
    license_kinds = [
        "@rules_license//licenses/spdx:Apache-2.0",
    ],
    license_text = "src/LICENSE",
)

alias(
    name = "conda_package_py_opencv_python",
    actual = ":cv2",
    visibility = ["@anaconda//:__pkg__"],
)

filegroup(
    name = "conda_package_py_opencv_hdrs",
    visibility = ["@anaconda//:__pkg__"],
)

filegroup(
    name = "conda_package_py_opencv_libs",
    visibility = ["@anaconda//:__pkg__"],
)

filegroup(
    name = "conda_package_py_opencv_solibs",
    visibility = ["@anaconda//:__pkg__"],
)

filegroup(
    name = "conda_package_py_opencv_data",
    srcs = [
        ":cv2",
    ],
    visibility = ["@anaconda//:__pkg__"],
)

exports_files(
    ["BUILD.bazel"],
    visibility = ["//visibility:public"],
)

filegroup(
    name = "srcs",
    srcs = glob(
        ["src/**/*"],
        exclude = [
            "src/doc/**",
        ],
    ) + ["src/doc/CMakeLists.txt"],
)

genrule(
    name = "cv2",
    srcs = [
        ":srcs",
        "version.txt",
        "src/CMakeLists.txt",
        "@anaconda//:eigen",
        "@anaconda//:jpeg",
        "@anaconda//:libffi",
        "@anaconda//:libgcc-ng",
        "@anaconda//:libpng",
        "@anaconda//:libstdcxx-ng",
        "@anaconda//:libtiff",
        "@anaconda//:mkl-include",
        "@anaconda//:mkl-service",
        "@anaconda//:mkl",
        "@anaconda//:numpy-base",
        "@anaconda//:python_interpreter_exe",
        "@anaconda//:python",
        "@anaconda//:six",
        "@anaconda//:tbb-devel",
        "@anaconda//:tbb",
        "@anaconda//:xz",
        "@anaconda//:zlib",
        "@anaconda//:zstd",
    ],
    outs = [PYTHON_PREFIX + "/site-packages/cv2.so"],
    cmd_bash = r'''
        export AR=$(AR)
        export CC=$(CC)
        export NM=$(NM)
        export OBJCOPY=$(OBJCOPY)
        export STRIP=$(STRIP)
        "$(execpath @cellranger//third-party:build_opencv)" \
            $$(cat "$(execpath version.txt)") \
            "$(execpath @tenx_toolchain//:cmake)" \
            "$(execpath @tenx_toolchain//:bin/ninja)" \
            "$(execpath @tenx_toolchain//:bin/patchelf)" \
            "$(execpath @anaconda//:python_interpreter_exe)" \
            "$(execpath src/CMakeLists.txt)" \
            "$@"''',
    exec_properties = {
        "Pool": "6core",
    },
    message = "Building OpenCV from source",
    toolchains = ["@bazel_tools//tools/cpp:current_cc_toolchain"],
    tools = [
        "@cellranger//third-party:build_opencv",
        "@tenx_toolchain//:bin/ninja",
        "@tenx_toolchain//:bin/patchelf",
        "@tenx_toolchain//:cmake",
        "@tenx_toolchain//:cmake_modules",
    ],
)

conda_manifest(
    name = "conda_metadata",
    info_files = ["info/index.json"],
    manifest = "info/files",
    visibility = ["//visibility:public"],
)

conda_deps(
    name = "conda_deps",
    visibility = ["@anaconda//:__pkg__"],
    deps = [
        "@anaconda//:jpeg",
        "@anaconda//:libgcc-ng",
        "@anaconda//:libpng",
        "@anaconda//:libstdcxx-ng",
        "@anaconda//:libtiff",
        "@anaconda//:mkl",
        "@anaconda//:numpy-base",
        "@anaconda//:python",
        "@anaconda//:tbb",
        "@anaconda//:xz",
        "@anaconda//:zlib",
        "@anaconda//:zstd",
    ],
)
