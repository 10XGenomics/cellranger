load("@conda_package_python//:vars.bzl", "PYTHON_PREFIX")
load(
    "@rules_license//rules:license.bzl",
    "license",
)
load(
    "@tenx_bazel_rules//rules/conda:conda_manifest.bzl",
    "conda_deps",
    "conda_files",
    "conda_manifest",
)

package(default_applicable_licenses = ["license"])

site_packages = PYTHON_PREFIX + "/site-packages"

license(
    name = "license",
    package_name = "fbpca",
    additional_info = {
        "homepage": "https://github.com/facebookarchive/fbpca",
        "version": "151c2364cd0ec1f",
        "manifest": "third-party/deps.bzl",
        "pURL": "pkg:github/facebookarchive/fbpca@151c2364cd0ec1fd7da6dec6d3ed2de9c6bfec5d",
    },
    copyright_notice = "Copyright (c) 2014, Facebook, Inc.",
    license_kinds = [
        "@rules_license//licenses/spdx:BSD-3-Clause",
    ],
    license_text = site_packages + "/LICENSE",
)

conda_files(
    name = "files",
    py_srcs = [site_packages + "/fbpca.py"],
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
        "@anaconda//:numpy",
        "@anaconda//:scipy",
    ],
)
