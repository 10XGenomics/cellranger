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

site_packages = PYTHON_PREFIX + "/site-packages"

license(
    name = "license",
    package_name = "cell-annotation-service-client",
    additional_info = {
        "homepage": "https://github.com/10XGenomics/cell-annotation-service-client",
        "version": "1.4.0.dev",
        "manifest": "third-party/deps.bzl",
        "pURL": "pkg:github/10XGenomics/cell-annotation-service-client@1.4.0.dev",
    },
    license_kinds = [
        "@rules_license//licenses/spdx:BSD-3-Clause",
    ],
    license_text = site_packages + "/LICENSE.txt",
)

filegroup(
    name = "conda_package_cellarium_hdrs",
    visibility = ["@anaconda//:__pkg__"],
)

filegroup(
    name = "conda_package_cellarium_libs",
    visibility = ["@anaconda//:__pkg__"],
)

filegroup(
    name = "conda_package_cellarium_solibs",
    visibility = ["@anaconda//:__pkg__"],
)

filegroup(
    name = "conda_package_cellarium_python",
    srcs = [site_packages + f for f in [
        "/cellarium/cas/__init__.py",
        "/cellarium/cas/_io.py",
        "/cellarium/cas/assets/cellarium_cas_tx_pca_002_grch38_2020_a.json",
        "/cellarium/cas/client.py",
        "/cellarium/cas/data_preparation.py",
        "/cellarium/cas/endpoints.py",
        "/cellarium/cas/exceptions.py",
        "/cellarium/cas/service.py",
    ]],
    visibility = ["@anaconda//:__pkg__"],
)

filegroup(
    name = "conda_package_cellarium_data",
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
        "@anaconda//:aiohttp",
        "@anaconda//:anndata",
        "@anaconda//:certifi",
        "@anaconda//:nest-asyncio",
        "@anaconda//:tqdm",
    ],
)
