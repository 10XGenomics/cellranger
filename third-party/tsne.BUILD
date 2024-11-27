load("@bazel_skylib//rules:copy_file.bzl", "copy_file")
load("@conda_package_python//:vars.bzl", "PYTHON_PREFIX")
load("@rules_cc//cc:defs.bzl", "cc_library")
load(
    "@rules_license//rules:license.bzl",
    "license",
)
load("@tenx_bazel_rules//rules:cython_library.bzl", "cython_library")
load(
    "@tenx_bazel_rules//rules/conda:conda_manifest.bzl",
    "conda_deps",
    "conda_files",
    "conda_manifest",
)

package(
    default_applicable_licenses = ["license"],
    features = ["thin_lto"],
)

site_packages = PYTHON_PREFIX + "/site-packages"

license(
    name = "license",
    package_name = "tsne",
    additional_info = {
        "homepage": "https://github.com/10XDev/tsne",
        "version": "45360b6",
        "manifest": "third-party/deps.bzl",
        "pURL": "pkg:github/10XDev/tsne@45360b6",
    },
    copyright_notice = "Copyright (c) 2013, Laurens van der Maaten (Delft University of Technology)",
    license_kinds = [
        "@rules_license//licenses/spdx:BSD-4-Clause",
    ],
    license_text = site_packages + "/LICENSE",
)

conda_files(
    name = "files",
    link_safe_runfiles = [
        ":bh_sne_3d_cp",
        ":bh_sne_cp",
    ],
    py_srcs = [
        site_packages + "/tsne/__init__.py",
        site_packages + "/tsne/_version.py",
    ],
    visibility = ["@anaconda//:__pkg__"],
)

copy_file(
    name = "bh_sne_cp",
    src = ":bh_sne.so",
    out = site_packages + "/bh_sne.so",
    allow_symlink = True,
)

copy_file(
    name = "bh_sne_3d_cp",
    src = ":bh_sne_3d.so",
    out = site_packages + "/bh_sne_3d.so",
    allow_symlink = True,
)

cython_library(
    name = "bh_sne",
    srcs = [site_packages + "/tsne/bh_sne.pyx"],
    c_deps = [
        ":tsne_2d",
        "@anaconda//:numpy-base",
    ],
    copts = [
        "-fwrapv",
        "-fvisibility=hidden",
        "-fvisibility-inlines-hidden",
        "-DNPY_NO_DEPRECATED_API=7",
    ],
    cpp = True,
    cython_flags = [
        "--line-directives",
        "-3",
    ],
    imports = [site_packages],
    srcs_version = "PY3",
    deps = [
        "@anaconda//:numpy",
    ],
)

cython_library(
    name = "bh_sne_3d",
    srcs = [site_packages + "/tsne/bh_sne_3d.pyx"],
    c_deps = [
        ":tsne_3d",
        "@anaconda//:numpy-base",
    ],
    copts = [
        "-fwrapv",
        "-fvisibility=hidden",
        "-fvisibility-inlines-hidden",
        "-DNPY_NO_DEPRECATED_API=7",
    ],
    cpp = True,
    cython_flags = [
        "--line-directives",
        "-3",
    ],
    imports = [site_packages],
    srcs_version = "PY3",
    deps = [
        "@anaconda//:numpy",
    ],
)

bh_sne_srcs = [
    site_packages + "/tsne/bh_sne_src/sptree.cpp",
    site_packages + "/tsne/bh_sne_src/tsne.cpp",
]

bh_sne_hdrs = [
    site_packages + "/tsne/bh_sne_src/sptree.h",
    site_packages + "/tsne/bh_sne_src/tsne.h",
    site_packages + "/tsne/bh_sne_src/vptree.h",
]

cc_library(
    name = "tsne_2d",
    srcs = bh_sne_srcs,
    hdrs = bh_sne_hdrs,
    copts = [
        "-fvisibility=hidden",
        "-Wno-unused-variable",
    ],
    linkstatic = True,
    strip_include_prefix = site_packages + "/tsne/bh_sne_src/",
)

cc_library(
    name = "tsne_3d",
    srcs = bh_sne_srcs,
    hdrs = bh_sne_hdrs,
    copts = [
        "-fvisibility=hidden",
        "-Wno-unused-variable",
    ],
    defines = ["TSNE3D"],
    linkstatic = True,
    strip_include_prefix = site_packages + "/tsne/bh_sne_src/",
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
