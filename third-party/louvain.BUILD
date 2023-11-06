load("@rules_cc//cc:defs.bzl", "cc_binary", "cc_library")
load(
    "@rules_license//rules:license.bzl",
    "license",
)

package(
    default_applicable_licenses = ["license"],
    features = ["thin_lto"],
)

license(
    name = "license",
    package_name = "louvain",
    additional_info = {
        "homepage": "https://github.com/10XGenomics/louvain",
        "version": "fb6919c",
        "manifest": "third-party/deps.bzl",
        "pURL": "pkg:github/10XGenomics/louvain@fb6919c",
    },
    copyright_notice = "E. Lefebvre, adapted by J.-L. Guillaume and R. Campigotto",
    license_kinds = [
        "@rules_license//licenses/spdx:GPL-3.0-or-later",
    ],
    license_text = "README.txt",
)

cc_binary(
    name = "louvain",
    srcs = [
        "src/main_louvain.cpp",
    ],
    copts = ["-Wno-unused-but-set-variable"],
    visibility = ["//visibility:public"],
    deps = [":louvain_lib"],
)

cc_binary(
    name = "convert",
    srcs = [
        "src/main_convert.cpp",
    ],
    visibility = ["//visibility:public"],
    deps = [
        ":graph",
    ],
)

cc_binary(
    name = "hierarchy",
    srcs = [
        "src/main_hierarchy.cpp",
    ],
    visibility = ["//visibility:public"],
)

cc_binary(
    name = "matrix",
    srcs = ["src/main_matrix.cpp"],
    visibility = ["//visibility:public"],
)

cc_library(
    name = "graph",
    srcs = ["src/graph.cpp"],
    hdrs = ["src/graph.h"],
    copts = ["-Wno-unused-but-set-variable"],
    linkstatic = 1,
    visibility = ["//visibility:public"],
)

cc_library(
    name = "louvain_lib",
    srcs = [
        "src/balmod.cpp",
        "src/condora.cpp",
        "src/devind.cpp",
        "src/devuni.cpp",
        "src/dp.cpp",
        "src/goldberg.cpp",
        "src/graph_binary.cpp",
        "src/louvain.cpp",
        "src/modularity.cpp",
        "src/owzad.cpp",
        "src/quality.cpp",
        "src/shimalik.cpp",
        "src/zahn.cpp",
    ],
    hdrs = [
        "src/balmod.h",
        "src/condora.h",
        "src/devind.h",
        "src/devuni.h",
        "src/dp.h",
        "src/goldberg.h",
        "src/graph_binary.h",
        "src/louvain.h",
        "src/modularity.h",
        "src/owzad.h",
        "src/quality.h",
        "src/shimalik.h",
        "src/zahn.h",
    ],
    copts = ["-Wno-unused-but-set-variable"],
    linkstatic = 1,
)
