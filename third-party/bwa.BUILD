load("@rules_cc//cc:defs.bzl", "cc_binary", "cc_library")
load(
    "@rules_license//rules:license.bzl",
    "license",
)

package(default_applicable_licenses = ["license"])

license(
    name = "license",
    package_name = "bwa",
    additional_info = {
        "homepage": "https://github.com/lh3/bwa",
        "version": "13b5637fe6b",
        "manifest": "third-party/deps.bzl",
        "pURL": "pkg:github/lh3/bwa@13b5637fe6b",
    },
    copyright_notice = "Copyright (c) 2018- Dana-Farber Cancer Institute, " +
                       "2009-2018 Broad Institute, Inc., " +
                       "2008-2009 Genome Research Ltd. (GRL)",
    license_kinds = [
        "@rules_license//licenses/spdx:GPL-3.0",
    ],
    license_text = "COPYING",
)

cc_binary(
    name = "bwa",
    srcs = ["main.c"],
    visibility = ["//visibility:public"],
    deps = [
        ":binlib",
        ":libbwa",
    ],
)

cc_library(
    name = "libbwa",
    linkstatic = 1,
    visibility = ["//visibility:public"],
    deps = [
        ":QSufSort",
        ":bntseq",
        ":bwa_c",
        ":bwamem",
        ":bwt",
        ":bwt_gen",
        ":bwtindex",
        ":is",
        ":kstring",
        ":ksw",
        ":kthread",
        ":malloc_wrap",
        ":rle",
        ":rope",
        ":utils",
    ],
)

cc_library(
    name = "binlib",
    linkstatic = 1,
    deps = [
        ":bamlite",
        ":bwape",
        ":bwase",
        ":bwashm",
        ":bwt_lite",
        ":bwtaln",
        ":bwtsw2",
        ":fastmap",
        ":kopen",
        ":maxk",
        ":pemerge",
    ],
)

DFLAGS = [
    "-DHAVE_PTHREAD",
    "-DUSE_MALLOC_WRAPPERS",
    "-Wno-unused-function",
    "-Wno-unused-but-set-variable",
    "-Wno-single-bit-bitfield-constant-conversion",
]

cc_library(
    name = "utils",
    srcs = ["utils.c"],
    hdrs = ["utils.h"],
    copts = DFLAGS,
    deps = [
        ":kseq",
        ":ksort",
        "@zlib",
    ],
)

cc_library(
    name = "kthread",
    srcs = ["kthread.c"],
    copts = DFLAGS,
    linkopts = ["-lpthread"],
)

cc_library(
    name = "kstring",
    srcs = ["kstring.c"],
    hdrs = ["kstring.h"],
    copts = DFLAGS,
    deps = [":malloc_wrap"],
)

cc_library(
    name = "ksw",
    srcs = ["ksw.c"],
    hdrs = ["ksw.h"],
    copts = DFLAGS,
    deps = [":malloc_wrap"],
)

cc_library(
    name = "bwt",
    srcs = ["bwt.c"],
    hdrs = ["bwt.h"],
    copts = DFLAGS,
    deps = [
        ":kvec",
        ":malloc_wrap",
        ":utils",
    ],
)

cc_library(
    name = "bntseq",
    srcs = ["bntseq.c"],
    hdrs = ["bntseq.h"],
    copts = DFLAGS,
    deps = [
        ":khash",
        ":kseq",
        ":malloc_wrap",
        ":utils",
    ],
)

cc_library(
    name = "bwa_c",
    srcs = ["bwa.c"],
    hdrs = ["bwa.h"],
    copts = DFLAGS,
    deps = [
        ":bntseq",
        ":bwt",
        ":kstring",
        ":ksw",
        ":kvec",
        ":utils",
    ],
)

cc_library(
    name = "bwamem",
    srcs = [
        "bwamem.c",
        "bwamem_extra.c",
        "bwamem_pair.c",
        "kbtree.h",
    ],
    hdrs = ["bwamem.h"],
    copts = DFLAGS,
    linkopts = ["-lm"],
    deps = [
        ":bntseq",
        ":bwa_c",
        ":bwt",
        ":ksort",
        ":kstring",
        ":ksw",
        ":kthread",
        ":kvec",
        ":malloc_wrap",
        ":utils",
    ],
)

cc_library(
    name = "malloc_wrap",
    srcs = ["malloc_wrap.c"],
    hdrs = ["malloc_wrap.h"],
    copts = DFLAGS,
)

cc_library(
    name = "QSufSort",
    srcs = ["QSufSort.c"],
    hdrs = ["QSufSort.h"],
    copts = DFLAGS,
)

cc_library(
    name = "bwt_gen",
    srcs = ["bwt_gen.c"],
    copts = DFLAGS,
    deps = [
        ":QSufSort",
        ":malloc_wrap",
    ],
)

cc_library(
    name = "rope",
    srcs = ["rope.c"],
    hdrs = ["rope.h"],
    copts = DFLAGS,
    deps = [":rle"],
)

cc_library(
    name = "rle",
    srcs = ["rle.c"],
    hdrs = ["rle.h"],
    copts = DFLAGS,
)

cc_library(
    name = "is",
    srcs = ["is.c"],
    copts = DFLAGS,
    deps = [":malloc_wrap"],
)

cc_library(
    name = "bwtindex",
    srcs = ["bwtindex.c"],
    copts = DFLAGS,
    deps = [
        ":bntseq",
        ":bwa_c",
        ":bwt",
        ":is",
        ":malloc_wrap",
        ":rle",
        ":rope",
        ":utils",
    ],
)

cc_library(
    name = "bwashm",
    srcs = ["bwashm.c"],
    copts = DFLAGS,
    linkopts = ["-lrt"],
    deps = [":bwa_c"],
)

cc_library(
    name = "bwase",
    srcs = ["bwase.c"],
    hdrs = ["bwase.h"],
    copts = DFLAGS,
    deps = [
        ":bntseq",
        ":bwa_c",
        ":bwt",
        ":bwtaln",
        ":kstring",
        ":ksw",
        ":utils",
    ],
)

cc_library(
    name = "bwtaln",
    srcs = [
        "bwaseqio.c",
        "bwtaln.c",
        "bwtgap.c",
    ],
    hdrs = [
        "bwtaln.h",
        "bwtgap.h",
    ],
    copts = DFLAGS,
    deps = [
        ":bamlite",
        ":bwa_c",
        ":bwt",
        ":kseq",
        ":malloc_wrap",
        ":utils",
    ],
)

cc_library(
    name = "bamlite",
    srcs = ["bamlite.c"],
    hdrs = ["bamlite.h"],
    copts = DFLAGS,
    deps = [":malloc_wrap"],
)

cc_library(
    name = "bwape",
    srcs = ["bwape.c"],
    copts = DFLAGS,
    deps = [
        ":bntseq",
        ":bwa_c",
        ":bwase",
        ":bwtaln",
        ":ksw",
        ":kvec",
        ":malloc_wrap",
        ":utils",
    ],
)

cc_library(
    name = "kopen",
    srcs = ["kopen.c"],
    copts = DFLAGS,
    deps = [":malloc_wrap"],
)

cc_library(
    name = "pemerge",
    srcs = ["pemerge.c"],
    copts = DFLAGS,
    deps = [
        ":bwa_c",
        ":kseq",
        ":kstring",
        ":ksw",
        ":malloc_wrap",
        ":utils",
    ],
)

cc_library(
    name = "maxk",
    srcs = ["maxk.c"],
    copts = DFLAGS,
    deps = [
        ":bwa_c",
        ":bwamem",
        ":kseq",
    ],
)

cc_library(
    name = "bwtsw2",
    srcs = [
        "bwtsw2_aux.c",
        "bwtsw2_chain.c",
        "bwtsw2_core.c",
        "bwtsw2_main.c",
        "bwtsw2_pair.c",
    ],
    hdrs = [
        "bwtsw2.h",
    ],
    copts = DFLAGS,
    deps = [
        ":bntseq",
        ":bwa_c",
        ":bwt",
        ":bwt_lite",
        ":bwtindex",
        ":kseq",
        ":ksort",
        ":kstring",
        ":ksw",
        ":malloc_wrap",
        ":utils",
    ],
)

cc_library(
    name = "bwt_lite",
    srcs = ["bwt_lite.c"],
    hdrs = ["bwt_lite.h"],
    copts = DFLAGS,
    deps = [
        ":malloc_wrap",
    ],
)

cc_library(
    name = "fastmap",
    srcs = ["fastmap.c"],
    copts = DFLAGS,
    deps = [
        ":bntseq",
        ":bwa_c",
        ":bwamem",
        ":kopen",
        ":kseq",
        ":kvec",
        ":utils",
    ],
)

cc_library(
    name = "kseq",
    hdrs = ["kseq.h"],
    deps = [":malloc_wrap"],
)

cc_library(
    name = "ksort",
    hdrs = ["ksort.h"],
    deps = [":malloc_wrap"],
)

cc_library(
    name = "khash",
    hdrs = ["khash.h"],
    deps = [":malloc_wrap"],
)

cc_library(
    name = "kvec",
    hdrs = ["kvec.h"],
    deps = [":malloc_wrap"],
)
