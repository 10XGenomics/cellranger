"""Repository rule invocations for c dependencies."""

load(
    "@bazel_tools//tools/build_defs/repo:http.bzl",
    "http_archive",
)
load(
    "@tenx_bazel_rules//rules:git_http_archive.bzl",
    "github_http_archive",
)

def load_c_deps():
    """Loads c dependencies."""

    bwa_commit = "13b5637fe6bd678b5756247a857c7ed9460b8e77"
    http_archive(
        name = "bwa",
        build_file = "@cellranger//third-party:bwa.BUILD",
        urls = [
            "https://github.com/lh3/bwa/archive/{}.tar.gz".format(bwa_commit),
        ],
        sha256 = "fa8700e82167ecba76300845a8a0480097cb252b447915e1ec0bc85164915f43",
        patch_args = ["-p1"],
        patches = [
            "@cellranger//third-party:bwa-unrestrict-mergedBwt.patch",
        ],
        strip_prefix = "bwa-" + bwa_commit,
    )

    # this points to 2.7.2a-tenx branch of the STAR repo
    star_commit = "230afea1013d462e19fb31f56791f76c355d871c"
    http_archive(
        name = "STAR",
        build_file = "@cellranger//third-party:STAR.BUILD",
        urls = [
            "https://github.com/10XGenomics/STAR/archive/{}.tar.gz".format(star_commit),
        ],
        patch_cmds = [
            "mv source/bam_cat.c source/bam_cat.cpp",
            "echo \"#define HTS_VERSION \\\"0.0.1\\\"\" > source/htslib/version.h",
        ],
        sha256 = "cdc36c7d76932099b5937f1c1455401924c08569b8e26487ddb04d913513927f",
        strip_prefix = "STAR-" + star_commit,
    )

    http_archive(
        name = "com_github_samtools",
        build_file = "@cellranger//third-party:samtools.BUILD",
        urls = [
            "https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2",
        ],
        sha256 = "2fa0a25f78594cf23d07c9d32d5060a14f1c5ee14d7b0af7a8a71abc9fdf1d07",
        strip_prefix = "samtools-1.16.1",
    )

    louvain_commit = "fb6919ca85e0284f124f43cead25bafdb0e988df"
    http_archive(
        name = "louvain",
        build_file = "@cellranger//third-party:louvain.BUILD",
        urls = [
            "https://github.com/10XGenomics/louvain/archive/{}.tar.gz".format(louvain_commit),
        ],
        sha256 = "83c66c7acb83e2b6b1e890bcb0cff916af500ab16900ab67d4dd48fb438bccbf",
        strip_prefix = "louvain-" + louvain_commit,
    )

    # tracking the "bazel" branch in 10xdev/graphviz
    github_http_archive(
        name = "com_github_10xdev_graphviz",
        commit = "11a86ea0ca32a79d65102f90eae122dba203ddd6",
        url = "https://github.com/10XDev/graphviz",
        sha256 = "aa7842ce0d400277531ba89016218b6426c2245deee6b2668d72c25aa055c2ef",
    )
