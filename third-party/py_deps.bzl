"""Repository rule invocations for third-party python dependencies."""

load(
    "@tenx_bazel_rules//rules:conda_package_repository.bzl",
    "conda_package_repository",
)
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
        build_file = "@cellranger//third-party:opencv.BUILD",
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
        patches = ["@cellranger//third-party:open-cv-warnings.patch", "@cellranger//third-party:open-cv-grabcut.patch"],
        sha256 = "c20bb83dd790fc69df9f105477e24267706715a9d3c705ca1e7f613c7b3bad3d",
    )

    # This package over-constrains libtiff, which makes other updates impossible.
    # So, just pull it in separately without letting it take part in the solve.
    conda_package_repository(
        name = "conda_package_simpleitk",
        base_urls = ["https://conda.anaconda.org/simpleitk/linux-64"],
        dist_name = "simpleitk-2.2.1-py310h3fd9d12_0",
        sha256 = "bbe4fce613459cf15d5c176b30ecc0d8b7c17d594822c36af9b3d735c7283aee",
    )

    tenxpy_commit = "c99de5e11ed7c564752c2c44c9bb8e4890dac1bf"

    new_conda_package_http_repository(
        name = "conda_package_tenxpy",
        build_file = "@cellranger//third-party:tenxpy.BUILD",
        urls = [
            "https://github.com/10XDev/tenxpy/archive/{}.tar.gz".format(
                tenxpy_commit,
            ),
        ],
        strip_prefix = "tenxpy-" + tenxpy_commit,
        add_prefix = site_packages,
        sha256 = "b1aa4089b075d02b089e71ac9da8b4a69bb3ea611d82810bf142af2d95f0490a",
        auth_patterns = {
            "github.com": "token <password>",
        },
    )

    tsne_commit = "45360b641183b6d899bdef201455d89de59a03a2"
    new_conda_package_http_repository(
        name = "conda_package_tsne",
        build_file = "@cellranger//third-party:tsne.BUILD",
        urls = [
            "https://github.com/10XDev/tsne/archive/{}.tar.gz".format(
                tsne_commit,
            ),
        ],
        strip_prefix = "tsne-" + tsne_commit,
        add_prefix = site_packages,
        sha256 = "7c7c29c19a703c4ec30d3d5886656f35ae4352e070c4106115336e317279a4c0",
    )

    hjson_version = "3.0.2"

    new_conda_package_http_repository(
        name = "conda_package_hjson",
        build_file = "@cellranger//third-party:hjson.BUILD",
        urls = [
            "https://github.com/hjson/hjson-py/archive/v{}.tar.gz".format(hjson_version),
        ],
        sha256 = "1c16084568a6328ba404703aad6db8ddf7e3b6afec9e141b375a40a61a209c30",
        strip_prefix = "hjson-py-" + hjson_version,
        add_prefix = site_packages,
    )

    ipytablewidgets_version = "0.3.1"

    new_conda_package_http_repository(
        name = "conda_package_ipytablewidgets",
        build_file = "@cellranger//third-party:ipytablewidgets.BUILD",
        urls = [
            "https://files.pythonhosted.org/packages/a8/24/ba9f1126efb9399d22f27b459b5bdcb30f2e58f6f9f9645ca043896f1421/ipytablewidgets-{}-py2.py3-none-any.whl".format(ipytablewidgets_version),
        ],
        sha256 = "6624cbeb73c73c68c5900a3a7719bf6aab560455e76247d373f3fa50336342a8",
        add_prefix = site_packages,
        patch_args = ["-p1"],
        patches = ["@cellranger//third-party:ipytablewidgets_dependency_fix.patch"],
        type = "zip",
    )

    faicons_version = "0.2.1"

    new_conda_package_http_repository(
        name = "conda_package_faicons",
        build_file = "@cellranger//third-party:faicons.BUILD",
        urls = [
            "https://files.pythonhosted.org/packages/b1/90/1a812cbaab68925d2a331ca83dbaad6ab997cbf6f49c9e98ed1088aa8b76/faicons-{}-py3-none-any.whl".format(faicons_version),
        ],
        sha256 = "8a32638798dfaae1f84c1d215fef2fa4936e51a07ea48b667422b6756943ab4c",
        add_prefix = site_packages,
        type = "zip",
    )

    pyarrow_version = "15.0.2"

    new_conda_package_http_repository(
        name = "conda_package_pyarrow",
        build_file = "@cellranger//third-party:pyarrow.BUILD",
        urls = [
            "https://files.pythonhosted.org/packages/01/e0/13aada7b0af1039554e675bd8c878acb3d86bab690e5a6b05fc8547a9cf2/pyarrow-{}-cp310-cp310-manylinux_2_17_x86_64.manylinux2014_x86_64.whl".format(pyarrow_version),
        ],
        sha256 = "f639c059035011db8c0497e541a8a45d98a58dbe34dc8fadd0ef128f2cee46e5",
        add_prefix = site_packages,
        type = "zip",
    )

    shiny_semantic_version = "0.1.3"

    new_conda_package_http_repository(
        name = "conda_package_shiny_semantic",
        build_file = "@cellranger//third-party:shiny-semantic.BUILD",
        urls = [
            "https://files.pythonhosted.org/packages/bd/1c/89aa57d6b1ca5421d4a6f8207095129aab116ba0ed9bed59a8b9f5b379df/shiny_semantic-{}-py3-none-any.whl".format(shiny_semantic_version),
        ],
        sha256 = "404e7d3496a1f2a15eb94b4fd25a543f3b502466223b0e7c2220aef3eb4d2e0a",
        add_prefix = site_packages,
        type = "zip",
    )

    vega_version = "4.0.0"

    new_conda_package_http_repository(
        name = "conda_package_vega",
        build_file = "@cellranger//third-party:vega.BUILD",
        urls = [
            "https://github.com/vega/ipyvega/archive/refs/tags/v{}.tar.gz".format(vega_version),
        ],
        sha256 = "a40002cc618bbd00fbe6f688024abb6a2805e93aed31df1f5f915d37c0c11688",
        strip_prefix = "ipyvega-" + vega_version,
        add_prefix = site_packages,
    )

    fbpca_commit = "151c2364cd0ec1fd7da6dec6d3ed2de9c6bfec5d"

    new_conda_package_http_repository(
        name = "conda_package_fbpca",
        build_file = "@cellranger//third-party:fbpca.BUILD",
        urls = [
            "https://github.com/facebookarchive/fbpca/archive/{}.tar.gz".format(fbpca_commit),
        ],
        sha256 = "49ea11626b7ef26ed3b2ff27e6efca04cd16e3e499442202a73d52f2e4a9f345",
        strip_prefix = "fbpca-" + fbpca_commit,
        add_prefix = site_packages,
    )

    new_conda_package_http_repository(
        name = "conda_package_moods",
        build_file = "@cellranger//third-party:moods.BUILD",
        urls = [
            "https://github.com/jhkorhonen/MOODS/releases/download/v1.9.3/MOODS-python-1.9.3.tar.gz",
        ],
        sha256 = "79d9ffe8acb7d32182dd190bfd55ad9e3170d1f69ab53ee7e243d2c1449f50d4",
        strip_prefix = "MOODS-python-1.9.3",
        add_prefix = site_packages,
    )

    csbdeep_version = "0.8.0"
    new_conda_package_http_repository(
        name = "conda_package_csbdeep",
        build_file = "@cellranger//third-party:csbdeep.BUILD",
        urls = [
            "https://files.pythonhosted.org/packages/fa/ae/16c5541707c15bfb3ec438dd738846f274899312094b323677e9f47fa3f5/csbdeep-{}-py2.py3-none-any.whl".format(csbdeep_version),
        ],
        sha256 = "366cfd039cc440331095f68786563afb54ecb1f35c6e0b629d045507a213f265",
        add_prefix = site_packages,
        type = "zip",
    )

    stardist_version = "0.9.1"
    new_conda_package_http_repository(
        name = "conda_package_stardist",
        build_file = "@cellranger//third-party:stardist.BUILD",
        urls = [
            "https://files.pythonhosted.org/packages/75/5e/4c0d2e48d064dec8f81be2b88404ce1f3440d0b86ee68dc43ab6924008da/stardist-{}-cp310-cp310-manylinux_2_17_x86_64.manylinux2014_x86_64.whl".format(stardist_version),
        ],
        sha256 = "f0f4542970e1c37472bb185bb36c3e4acb9e4efd65181cbbb89de917b7511e42",
        add_prefix = site_packages,
        type = "zip",
    )

    keras_version = "2.13.1-py3"
    new_conda_package_http_repository(
        name = "conda_package_keras",
        build_file = "@cellranger//third-party:keras.BUILD",
        urls = [
            "https://files.pythonhosted.org/packages/2e/f3/19da7511b45e80216cbbd9467137b2d28919c58ba1ccb971435cb631e470/keras-{}-none-any.whl".format(keras_version),
        ],
        sha256 = "5ce5f706f779fa7330e63632f327b75ce38144a120376b2ae1917c00fa6136af",
        add_prefix = site_packages,
        type = "zip",
    )
