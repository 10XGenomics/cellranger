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
    package_name = "hjson",
    additional_info = {
        "homepage": "https://github.com/hjson/hjson-py",
        "version": "3.0.2",
        "manifest": "third-party/deps.bzl",
        "pURL": "pkg:github/hjson/hjson-py@3.0.2",
    },
    copyright_notice = "Copyright (c) 2006 Bob Ippolito, " +
                       "Copyright (c) 2015 Christian Zangl",
    license_kinds = [
        "@rules_license//licenses/spdx:MIT",
    ],
    license_text = site_packages + "/LICENSE.txt",
)

conda_files(
    name = "files",
    py_srcs = [site_packages + "/hjson/" + f for f in [
        "__init__.py",
        "compat.py",
        "decoder.py",
        "encoder.py",
        "encoderH.py",
        "ordered_dict.py",
        "scanner.py",
        "tool.py",
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
)
