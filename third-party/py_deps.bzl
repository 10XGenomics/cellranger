"""Repository rule invocations for third-party python dependencies."""

load(
    "@tenx_bazel_rules//rules:new_conda_package_repository.bzl",
    "new_conda_package_http_repository",
)

def load_py_deps():
    """Loads third-party python dependencies for the pipeline."""

    site_packages = "lib/python3.10/site-packages"

    opencv_version = "4.5.4"

    new_conda_package_http_repository(
        name = "conda_package_py_opencv",
        package_name = "py-opencv",
        build_file = "//third-party:opencv.BUILD",
        urls = [
            "https://github.com/opencv/opencv/archive/refs/tags/{}.tar.gz".format(opencv_version),
        ],
        strip_prefix = "opencv-" + opencv_version,
        add_prefix = "src",
        patch_cmds = [
            "echo '{}' > version.txt".format(opencv_version),
            # Remove unused code.  Better to do it here for more consistency
            # between remote and local build.
            "rm -rf src/apps",
            "rm -rf src/modules/dnn",
            "rm -rf src/modules/gapi",
            "rm -rf src/modules/highgui",
            "rm -rf src/modules/java",
            "rm -rf src/modules/js",
            "rm -rf src/modules/ml",
            "rm -rf src/modules/photo",
            "rm -rf src/modules/stitching",
            "rm -rf src/modules/ts",
            "rm -rf src/modules/video",
            "rm -rf src/modules/videoio",
            "rm -rf src/modules/world",
            "rm -rf src/platforms/android",
            "rm -rf src/platforms/ios",
            "rm -rf src/platforms/js",
            "rm -rf src/platforms/maven",
            "rm -rf src/platforms/win*",
            "rm -rf src/samples",
        ],
        patch_args = ["-p1"],
        patches = ["//third-party:open-cv-warnings.patch"],
        sha256 = "c20bb83dd790fc69df9f105477e24267706715a9d3c705ca1e7f613c7b3bad3d",
        exported_files = [
            site_packages + "/cv2.so",
        ],
        license = "Apache-2.0",
        version = opencv_version,
    )

    tenxpy_commit = "c99de5e11ed7c564752c2c44c9bb8e4890dac1bf"

    new_conda_package_http_repository(
        name = "conda_package_tenxpy",
        build_file = "//third-party:tenxpy.BUILD",
        urls = [
            "https://github.com/10XDev/tenxpy/archive/{}.tar.gz".format(
                tenxpy_commit,
            ),
        ],
        strip_prefix = "tenxpy-" + tenxpy_commit,
        add_prefix = site_packages,
        exported_files = [site_packages + "/tenxpy/" + f for f in [
            "__init__.py",
            "constants.py",
            "lena2.py",
            "utils.py",
        ]],
        sha256 = "b1aa4089b075d02b089e71ac9da8b4a69bb3ea611d82810bf142af2d95f0490a",
        auth_patterns = {
            "github.com": "token <password>",
        },
        license = "10X Genomics",
    )

    tsne_commit = "45360b641183b6d899bdef201455d89de59a03a2"
    new_conda_package_http_repository(
        name = "conda_package_tsne",
        build_file = "//third-party:tsne.BUILD",
        urls = [
            "https://github.com/10XDev/tsne/archive/{}.tar.gz".format(
                tsne_commit,
            ),
        ],
        strip_prefix = "tsne-" + tsne_commit,
        add_prefix = site_packages,
        exported_files = [
            site_packages + "/bh_sne.so",
            site_packages + "/bh_sne_3d.so",
            site_packages + "/tsne/__init__.py",
            site_packages + "/tsne/_version.py",
        ],
        sha256 = "7c7c29c19a703c4ec30d3d5886656f35ae4352e070c4106115336e317279a4c0",
        license = "BSD-4-Clause",
        version = "0.1.5-10x",
    )

    cas_version = "0.1.3"
    cas_prefix = site_packages + "/cas/"
    new_conda_package_http_repository(
        name = "conda_package_cas",
        build_file = "//third-party:cell-annotation-service-client.BUILD",
        urls = [
            "https://github.com/10XGenomics/cell-annotation-service-client/archive/v{}.tar.gz".format(cas_version),
        ],
        strip_prefix = "cell-annotation-service-client-" + cas_version + "/src/",
        add_prefix = cas_prefix,
        exported_files = [cas_prefix + f for f in [
            "cas_cli/service.py",
            "cas_cli/__init__.py",
            "cas_cli/exceptions.py",
            "cas_helper/__init__.py",
            "cas_helper/_cli_helper.py",
        ]],
        sha256 = "35b9bffa792310f4c9066b2f5ce7a2d7d5bd8081178d58ca92e1f8880f578b26",
        license = "BSD-3",
        version = cas_version,
    )

    hjson_version = "3.0.2"

    new_conda_package_http_repository(
        name = "conda_package_hjson",
        build_file = "//third-party:hjson.BUILD",
        urls = [
            "https://github.com/hjson/hjson-py/archive/v{}.tar.gz".format(hjson_version),
        ],
        sha256 = "1c16084568a6328ba404703aad6db8ddf7e3b6afec9e141b375a40a61a209c30",
        strip_prefix = "hjson-py-" + hjson_version,
        add_prefix = site_packages,
        exported_files = [site_packages + "/hjson/" + f for f in [
            "__init__.py",
            "compat.py",
            "decoder.py",
            "encoder.py",
            "encoderH.py",
            "ordered_dict.py",
            "scanner.py",
            "tool.py",
        ]],
        license = "MIT and AFL-2.1",
        version = hjson_version,
    )

    fbpca_commit = "151c2364cd0ec1fd7da6dec6d3ed2de9c6bfec5d"

    new_conda_package_http_repository(
        name = "conda_package_fbpca",
        build_file = "//third-party:fbpca.BUILD",
        urls = [
            "https://github.com/facebookarchive/fbpca/archive/{}.tar.gz".format(fbpca_commit),
        ],
        sha256 = "49ea11626b7ef26ed3b2ff27e6efca04cd16e3e499442202a73d52f2e4a9f345",
        strip_prefix = "fbpca-" + fbpca_commit,
        add_prefix = site_packages,
        exported_files = [site_packages + "/fbpca.py"],
        version = "1.0",
        license = "BSD-3-Clause",
    )

    new_conda_package_http_repository(
        name = "conda_package_moods",
        build_file = "//third-party:moods.BUILD",
        urls = [
            "https://github.com/jhkorhonen/MOODS/releases/download/v1.9.3/MOODS-python-1.9.3.tar.gz",
        ],
        sha256 = "79d9ffe8acb7d32182dd190bfd55ad9e3170d1f69ab53ee7e243d2c1449f50d4",
        strip_prefix = "MOODS-python-1.9.3",
        add_prefix = site_packages,
        exported_files = [
            site_packages + "/MOODS/__init__.py",
            site_packages + "/MOODS/_parsers.so",
            site_packages + "/MOODS/_scan.so",
            site_packages + "/MOODS/_tools.so",
            site_packages + "/MOODS/misc.py",
            site_packages + "/MOODS/parsers.py",
            site_packages + "/MOODS/scan.py",
            site_packages + "/MOODS/tools.py",
        ],
        license = "Biopython License Agreement and GPL-3.0",
        version = "1.9.3",
    )
