"""Repository rule invocations for third-party dependencies."""

load(":c_deps.bzl", "load_c_deps")
load(":go_deps.bzl", "load_go_deps")
load(":py_deps.bzl", "load_py_deps")
load(":refdata.bzl", "load_refdata")
load(":rust_deps.bzl", "load_rust_deps")

def load_deps():
    """Loads third-party dependencies for the pipeline."""
    load_c_deps()
    load_go_deps()
    load_rust_deps()
    load_py_deps()
    load_refdata()
