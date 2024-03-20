load("@conda_package_python//:vars.bzl", "PYTHON_PREFIX")
load("@rules_cc//cc:defs.bzl", "cc_library")
load(
    "@rules_license//rules:license.bzl",
    "license",
)
load(
    "@tenx_bazel_rules//rules/conda:conda_manifest.bzl",
    "conda_deps",
    "conda_manifest",
)

package(
    default_applicable_licenses = ["license"],
    features = ["thin_lto"],
)

site_packages = PYTHON_PREFIX + "/site-packages"

license(
    name = "license",
    package_name = "MOODS",
    additional_info = {
        "homepage": "https://github.com/jhkorhonen/MOODS",
        "version": "1.9.3",
        "manifest": "third-party/deps.bzl",
        "pURL": "pkg:github/jhkorhonen/MOODS@1.9.3",
    },
    copyright_notice = "Copyright (C) 2009 by Petri J Martinmäki, Janne Korhonen, Pasi Rastas",
    license_kinds = [
        "@tenx_bazel_rules//licensing/known:Biopython",
    ],
    license_text = site_packages + "/COPYING.BIOPYTHON",
)

filegroup(
    name = "conda_package_moods_python",
    srcs = [site_packages + "/MOODS/" + f for f in [
        "__init__.py",
        "misc.py",
        "parsers.py",
        "scan.py",
        "tools.py",
    ]],
    visibility = ["@anaconda//:__pkg__"],
)

filegroup(
    name = "conda_package_moods_hdrs",
    visibility = ["@anaconda//:__pkg__"],
)

filegroup(
    name = "conda_package_moods_libs",
    visibility = ["@anaconda//:__pkg__"],
)

filegroup(
    name = "conda_package_moods_solibs",
    visibility = ["@anaconda//:__pkg__"],
)

exports_files(
    ["BUILD.bazel"],
    visibility = ["//visibility:public"],
)

core = site_packages + "/core"

copts = [
    "-fvisibility=hidden",
    "-fvisibility-inlines-hidden",
    "-std=c++11",
    "-Wno-bitwise-op-parentheses",
    "-Wno-delete-abstract-non-virtual-dtor",
    "-Wno-deprecated-declarations",
    "-Wno-unused-variable",
]

TOOLS = site_packages + "/MOODS/_tools.so"

cc_binary(
    name = TOOLS,
    srcs = [
        core + "/tools_wrap.cxx",
    ],
    copts = copts,
    linkopts = [
        "-Wl,--export-dynamic-symbol=PyInit__tools",
        "-Wl,--undefined=PyInit__tools",
        "-Wl,-rpath=$$ORIGIN/../../..",
    ],
    linkshared = True,
    deps = [
        ":common",
        "@anaconda//:python",
    ],
)

SCAN = site_packages + "/MOODS/_scan.so"

cc_binary(
    name = SCAN,
    srcs = [core + "/scan_wrap.cxx"],
    copts = copts,
    linkopts = [
        "-Wl,--export-dynamic-symbol=PyInit__scan",
        "-Wl,--undefined=PyInit__scan",
        "-Wl,-rpath=$$ORIGIN/../../..",
    ],
    linkshared = True,
    deps = [
        ":common",
        ":scan",
        "@anaconda//:python",
    ],
)

cc_library(
    name = "scan",
    srcs = [
        core + "/moods_scan.cpp",
        core + "/motif_0.cpp",
        core + "/motif_h.cpp",
        core + "/scanner.cpp",
    ],
    hdrs = [
        core + "/motif.h",
        core + "/scanner.h",
        core + "/moods_scan.h",
    ],
    copts = copts,
    linkstatic = True,
    deps = [":common"],
)

PARSERS = site_packages + "/MOODS/_parsers.so"

cc_binary(
    name = PARSERS,
    srcs = [core + "/parsers_wrap.cxx"],
    copts = copts,
    linkopts = [
        "-Wl,--export-dynamic-symbol=PyInit__parsers",
        "-Wl,--undefined=PyInit__parsers",
        "-Wl,-rpath=$$ORIGIN/../../..",
    ],
    linkshared = True,
    deps = [
        ":common",
        ":parsers",
        "@anaconda//:python",
    ],
)

cc_library(
    name = "parsers",
    srcs = [
        core + "/moods_parsers.cpp",
    ],
    hdrs = [
        core + "/moods_parsers.h",
    ],
    copts = copts,
    linkstatic = True,
    deps = [":common"],
)

cc_library(
    name = "common",
    srcs = [
        core + "/moods_tools.cpp",
        core + "/match_types.cpp",
        core + "/moods_misc.cpp",
    ],
    hdrs = [
        core + "/moods_tools.h",
        core + "/match_types.h",
        core + "/moods_misc.h",
        core + "/moods.h",
    ],
    copts = copts,
    linkstatic = True,
)

filegroup(
    name = "conda_package_moods_data",
    srcs = [
        PARSERS,
        SCAN,
        TOOLS,
    ],
    visibility = ["@anaconda//:__pkg__"],
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
)
