"""Repository rule invocations for third-party go dependencies."""

# gazelle:repo bazel_gazelle
load("@bazel_gazelle//:deps.bzl", "go_repository")
load("@tenx_cloud_cli//:deps.bzl", "cloud_cli_dependencies")
load("@tenx_loupekit//:deps.bzl", "loupekit_dependencies")

def load_go_deps():
    """Loads workspace dependencies for go builds."""
    _maybe_go_repository(
        name = "com_github_youtube_vitess",
        commit = "89cc312b4da3004d3b5382cf2b37d70e901b1c36",
        importpath = "github.com/youtube/vitess",
        patch_args = ["-p1"],
        patches = ["@cellranger//third-party:com_github_youtube_vitess.patch"],
        build_file_proto_mode = "disable_global",
        build_directives = [
            "gazelle:exclude data",
            "gazelle:exclude doc",
            "gazelle:exclude docker",
            "gazelle:exclude docs",
            "gazelle:exclude examples",
            "gazelle:exclude go/acl",
            "gazelle:exclude go/cache",
            "gazelle:exclude go/cacheservice",
            "gazelle:exclude go/cmd",
            "gazelle:exclude go/event",
            "gazelle:exclude go/ewma",
            "gazelle:exclude go/exit",
            "gazelle:exclude go/fileutil",
            "gazelle:exclude go/flagutil",
            "gazelle:exclude go/hack",
            "gazelle:exclude go/history",
            "gazelle:exclude go/ioutil2",
            "gazelle:exclude go/memcache",
            "gazelle:exclude go/mysql",
            "gazelle:exclude go/mysqlconn",
            "gazelle:exclude go/netutil",
            "gazelle:exclude go/pools",
            "gazelle:exclude go/proc",
            "gazelle:exclude go/race",
            "gazelle:exclude go/ratelimiter",
            "gazelle:exclude go/sqldb",
            "gazelle:exclude go/sqltypes",
            "gazelle:exclude go/stats",
            "gazelle:exclude go/streamlog",
            "gazelle:exclude go/sync2",
            "gazelle:exclude go/tb",
            "gazelle:exclude go/testfiles",
            "gazelle:exclude go/timer",
            "gazelle:exclude go/trace",
            "gazelle:exclude go/vt",
            "gazelle:exclude go/zk",
            "gazelle:exclude helm",
            "gazelle:exclude java",
            "gazelle:exclude misc",
            "gazelle:exclude php",
            "gazelle:exclude proto",
            "gazelle:exclude py",
            "gazelle:exclude test",
            "gazelle:exclude third_party",
            "gazelle:exclude tools",
            "gazelle:exclude travis",
            "gazelle:exclude vendor",
            "gazelle:exclude vitess.io",
            "gazelle:exclude web",
        ],
    )

    loupekit_dependencies()

    cloud_cli_dependencies()

    _maybe_go_repository(
        name = "com_github_stretchr_testify",
        commit = "8019298d9fa5a04fc2ad10ae03349df3483096a6",
        importpath = "github.com/stretchr/testify",
    )

    _maybe_go_repository(
        name = "in_gopkg_stretchr_testify_v1",
        commit = "f35b8ab0b5a2cef36673838d662e249dd9c94686",
        importpath = "gopkg.in/stretchr/testify.v1",
    )

    _maybe_go_repository(
        name = "in_gopkg_yaml_v3",
        sum = "h1:fxVm/GzAzEWqLHuvctI91KS9hhNmmWOoWu0XTYJS7CA=",
        version = "v3.0.1",
        importpath = "gopkg.in/yaml.v3",
    )

def _maybe_go_repository(name, **kwargs):
    if name not in native.existing_rules():
        go_repository(name = name, **kwargs)
