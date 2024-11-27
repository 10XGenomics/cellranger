load("@conda_package_python//:vars.bzl", "PYTHON_PREFIX")
load(
    "@tenx_bazel_rules//rules/conda:conda_manifest.bzl",
    "conda_deps",
    "conda_files",
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

conda_files(
    name = "files",
    py_srcs = [PYTHON_PREFIX + "/site-packages/tenxpy/" + f for f in [
        "__init__.py",
        "constants.py",
        "lena2.py",
        "utils.py",
    ]],
    visibility = ["@anaconda//:__pkg__"],
)

conda_manifest(
    name = "conda_metadata",
    info_files = ["info/index.json"],
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
