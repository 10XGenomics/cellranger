"""Repository rule invocations for cargo dependencies."""

load(
    "@tenx_bazel_rules//rules:cargo_repository.bzl",
    "cargo_repository",
)

def load_rust_deps():
    cargo_repository(
        name = "bamtofastq_cargo_dependencies",
        lockfile = "@bamtofastq//:Cargo.lock",
        cargo_config = "@bamtofastq//:.cargo/config.toml",
        srcs = [
            "@bamtofastq//:Cargo.toml",
        ],
    )
    cargo_repository(
        name = "cr_rust_cargo_dependencies",
        lockfile = "//lib/rust:Cargo.lock",
        env = {"PYO3_NO_PYTHON": "1"},
        cargo_config = "//lib/rust:.cargo/config",
        srcs = [
            "//lib/rust:toml",
        ],
    )
    cargo_repository(
        name = "cr_ana_cargo_dependencies",
        lockfile = "//lib/rust/cr_ana:Cargo.lock",
        cargo_config = "//lib/rust/cr_ana:.cargo/config",
        srcs = [
            "//lib/rust/cr_ana:Cargo.toml",
        ],
    )
    cargo_repository(
        name = "pseudoaligner_cargo_dependencies",
        lockfile = "@com_github_10XGenomics_rust_pseudoaligner//:Cargo.lock",
        srcs = ["@com_github_10XGenomics_rust_pseudoaligner//:Cargo.toml"],
    )
