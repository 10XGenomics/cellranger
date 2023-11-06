load(
    "@rules_license//rules:license.bzl",
    "license",
)
load(
    "@tenx_bazel_rules//rules/conda:conda_manifest.bzl",
    "conda_deps",
    "conda_manifest",
)
load("@conda_package_python//:vars.bzl", "PYTHON_PREFIX")

package(default_applicable_licenses = ["license"])

site_packages = PYTHON_PREFIX + "/site-packages"

license(
    name = "license",
    package_name = "cell-annotation-service-client",
    additional_info = {
        "homepage": "https://github.com/10XGenomics/cell-annotation-service-client",
        "version": "0.1.0",
        "manifest": "third-party/deps.bzl",
        "pURL": "pkg:github/10XGenomics/cell-annotation-service-client@v0.1.0",
    },
    license_kinds = [
        "@rules_license//licenses/spdx:BSD-3-Clause",
    ],
    license_text = site_packages + "/LICENSE.txt",
)

filegroup(
    name = "conda_package_cas_hdrs",
    visibility = ["@anaconda//:__pkg__"],
)

filegroup(
    name = "conda_package_cas_libs",
    visibility = ["@anaconda//:__pkg__"],
)

filegroup(
    name = "conda_package_cas_solibs",
    visibility = ["@anaconda//:__pkg__"],
)

filegroup(
    name = "conda_package_cas_python",
    srcs = [site_packages + "/cas" + f for f in [
        "/cas_cli/service.py",
        "/cas_cli/__init__.py",
        "/cas_cli/exceptions.py",
        "/cas_helper/__init__.py",
        "/cas_helper/_cli_helper.py",
    ]],
    visibility = ["@anaconda//:__pkg__"],
)

filegroup(
    name = "conda_package_cas_data",
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
