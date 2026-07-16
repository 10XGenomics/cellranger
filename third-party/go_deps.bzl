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

    # Dependencies of google_jsonnet_go (rules_jsonnet 0.6.0).
    # jsonnet_go_dependencies() is not called from WORKSPACE because it would
    # invoke go_register_toolchains(); these repos are inlined here instead.
    _maybe_go_repository(
        name = "com_github_davecgh_go_spew",
        importpath = "github.com/davecgh/go-spew",
        sum = "h1:vj9j/u1bqnvCEfJOwUhtlOARqs3+rkHYY13jYWTU97c=",
        version = "v1.1.1",
    )

    _maybe_go_repository(
        name = "com_github_fatih_color",
        importpath = "github.com/fatih/color",
        sum = "h1:s36xzo75JdqLaaWoiEHk767eHiwo0598uUxyfiPkDsg=",
        version = "v1.10.0",
    )

    _maybe_go_repository(
        name = "com_github_kr_pretty",
        importpath = "github.com/kr/pretty",
        sum = "h1:L/CwN0zerZDmRFUapSPitk6f+Q3+0za1rQkzVuMiMFI=",
        version = "v0.1.0",
    )

    _maybe_go_repository(
        name = "com_github_kr_pty",
        importpath = "github.com/kr/pty",
        sum = "h1:VkoXIwSboBpnk99O/KFauAEILuNHv5DVFKZMBN/gUgw=",
        version = "v1.1.1",
    )

    _maybe_go_repository(
        name = "com_github_kr_text",
        importpath = "github.com/kr/text",
        sum = "h1:45sCR5RtlFHMR4UwH9sdQ5TC8v0qDQCHnXt+kaKSTVE=",
        version = "v0.1.0",
    )

    _maybe_go_repository(
        name = "com_github_mattn_go_colorable",
        importpath = "github.com/mattn/go-colorable",
        sum = "h1:c1ghPdyEDarC70ftn0y+A/Ee++9zz8ljHG1b13eJ0s8=",
        version = "v0.1.8",
    )

    _maybe_go_repository(
        name = "com_github_mattn_go_isatty",
        importpath = "github.com/mattn/go-isatty",
        sum = "h1:wuysRhFDzyxgEmMf5xjvJ2M9dZoWAXNNr5LSBS7uHXY=",
        version = "v0.0.12",
    )

    _maybe_go_repository(
        name = "com_github_pmezard_go_difflib",
        importpath = "github.com/pmezard/go-difflib",
        sum = "h1:4DBwDE0NGyQoBHbLQYPwSUPoCMWR5BEzIk/f1lZbAQM=",
        version = "v1.0.0",
    )

    _maybe_go_repository(
        name = "com_github_sergi_go_diff",
        importpath = "github.com/sergi/go-diff",
        sum = "h1:we8PVUC3FE2uYfodKH/nBHMSetSfHDR6scGdBi+erh0=",
        version = "v1.1.0",
    )

    _maybe_go_repository(
        name = "com_github_stretchr_objx",
        importpath = "github.com/stretchr/objx",
        sum = "h1:4G4v2dO3VZwixGIRoQ5Lfboy6nUhCyYzaqnIAPPhYs4=",
        version = "v0.1.0",
    )

    _maybe_go_repository(
        name = "in_gopkg_check_v1",
        importpath = "gopkg.in/check.v1",
        sum = "h1:YR8cESwS4TdDjEe65xsg0ogRM/Nc3DYOhEAlW+xobZo=",
        version = "v1.0.0-20190902080502-41f04d3bba15",
    )

    _maybe_go_repository(
        name = "in_gopkg_yaml_v2",
        importpath = "gopkg.in/yaml.v2",
        sum = "h1:VUgggvou5XRW9mHwD/yXxIYSMtY0zoKQf/v226p2nyo=",
        version = "v2.2.7",
    )

    _maybe_go_repository(
        name = "io_k8s_sigs_yaml",
        importpath = "sigs.k8s.io/yaml",
        sum = "h1:4A07+ZFc2wgJwo8YNlQpr1rVlgUDlxXHhPJciaPY5gs=",
        version = "v1.1.0",
        patches = ["@cellranger//third-party:io_k8s_sigs_yaml.patch"],
        patch_args = ["-p1"],
    )

    _maybe_go_repository(
        name = "org_golang_x_sys",
        importpath = "golang.org/x/sys",
        sum = "h1:/WDfKMnPU+m5M4xB+6x4kaepxRw6jWvR5iDRdvjHgy8=",
        version = "v0.0.0-20200223170610-d5e6a3e2c0ae",
    )

def _maybe_go_repository(name, **kwargs):
    if name not in native.existing_rules():
        go_repository(name = name, **kwargs)
