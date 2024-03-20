load("@conda_package_python//:vars.bzl", "PYTHON_PREFIX")
load(
    "@tenx_bazel_rules//rules/conda:conda_manifest.bzl",
    "conda_deps",
    "conda_manifest",
)

filegroup(
    name = "conda_package_tenxpy_python",
    srcs = [PYTHON_PREFIX + "/site-packages/tenxpy/" + f for f in [
        "__init__.py",
        "constants.py",
        "lena2.py",
        "utils.py",
    ]],
    visibility = ["@anaconda//:__pkg__"],
)

filegroup(
    name = "conda_package_tenxpy_hdrs",
    visibility = ["@anaconda//:__pkg__"],
)

filegroup(
    name = "conda_package_tenxpy_libs",
    visibility = ["@anaconda//:__pkg__"],
)

filegroup(
    name = "conda_package_tenxpy_solibs",
    visibility = ["@anaconda//:__pkg__"],
)

filegroup(
    name = "conda_package_tenxpy_data",
    visibility = ["@anaconda//:__pkg__"],
)

exports_files(
    ["BUILD.bazel"],
    visibility = ["//visibility:public"],
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
        "@anaconda//:h5py",
        "@anaconda//:pandas",
        "@anaconda//:requests",
    ],
)
