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
        exported_files = [
            site_packages + "/cv2.so",
        ],
        license = "Apache-2.0",
        version = opencv_version,
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
        build_file = "@cellranger//third-party:tsne.BUILD",
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

    cas_commit = "9aa178dae7a0a03a8ac05a2fd3deb75e67388ea3"
    cas_version = "1.4.0.dev"
    cas_prefix = site_packages + "/cellarium/cas/"
    new_conda_package_http_repository(
        name = "conda_package_cellarium",
        build_file = "@cellranger//third-party:cell-annotation-service-client.BUILD",
        urls = [
            "https://github.com/cellarium-ai/cellarium-cas/archive/{}.tar.gz".format(cas_commit),
        ],
        strip_prefix = "cellarium-cas-" + cas_commit,
        add_prefix = site_packages,
        exported_files = [cas_prefix + f for f in [
            "service.py",
            "__init__.py",
            "_io.py",
            "assets/cellarium_cas_tx_pca_002_grch38_2020_a.json",
            "client.py",
            "data_preparation.py",
            "endpoints.py",
            "exceptions.py",
        ]],
        sha256 = "3c9a44141ba01c5c0f8a054cc9d8e8bd6d41acd849dcdb3e35cb974424f45c86",
        license = "BSD-3-Clause",
        version = cas_version,
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

    shinywidgets_version = "0.2.4"

    new_conda_package_http_repository(
        name = "conda_package_shinywidgets",
        build_file = "@cellranger//third-party:shinywidgets.BUILD",
        urls = [
            "https://github.com/posit-dev/py-shinywidgets/archive/refs/tags/v{}.tar.gz".format(shinywidgets_version),
        ],
        sha256 = "b167449404e46827d313973ce7d08ef4a98997fc07726a5e406e0cb791826f16",
        strip_prefix = "py-shinywidgets-" + shinywidgets_version,
        add_prefix = site_packages,
        exported_files = [site_packages + "/shinywidgets/" + f for f in [
            "__init__.py",
            "_as_widget.py",
            "_comm.py",
            "_dependencies.py",
            "_serialization.py",
            "_shinywidgets.py",
            "static/1e59d2330b4c6deb84b340635ed36249.ttf",
            "static/20fd1704ea223900efa9fd4e869efb08.woff2",
            "static/8b43027f47b20503057dfbbaa9401fef.eot",
            "static/c1e38fd9e0e74ba58f7a2b77ef29fdd3.svg",
            "static/f691f37e57f04c152e2315ab7dbad881.woff",
            "static/libembed-amd.js",
            "static/node_modules_codemirror_mode_sync_recursive_js_.output.js",
            "static/output.js",
            "static/vendors-node_modules_codemirror_mode_apl_apl_js-node_modules_codemirror_mode_asciiarmor_ascii-26282f.output.js",
            "static/shinywidgets.css",
        ]],
        license = "MIT",
        version = shinywidgets_version,
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
        exported_files = [site_packages + "/ipytablewidgets/" + f for f in [
            "__init__.py",
            "_frontend.py",
            "_version.py",
            "compressors.py",
            "numpy_adapter.py",
            "pandas_adapter.py",
            "progressivis_adapter.py",
            "serializers.py",
            "source_adapter.py",
            "traitlets.py",
            "widgets.py",
            "labextension/package.json",
            "labextension/static/149.40f4249249bc14416a21.js",
            "labextension/static/261.bc115825a135ccbef0e9.js",
            "labextension/static/405.6e1d012d86b0ed03bddf.js",
            "labextension/static/424.e42e48076aba64887fd2.js",
            "labextension/static/568.46074a0ad0e0cb2bad6c.js",
            "labextension/static/861.78ea6fdd74d109f94c90.js",
            "labextension/static/861.78ea6fdd74d109f94c90.js.LICENSE.txt",
            "labextension/static/remoteEntry.830facb2b8d2798f80f1.js",
            "labextension/static/style.js",
            "labextension/static/third-party-licenses.json",
            "static/extension.js",
            "static/index.js",
            "static/index.js.LICENSE.txt",
        ]],
        license = "BSD-3-Clause",
        version = ipytablewidgets_version,
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
        exported_files = [site_packages + "/faicons/" + f for f in [
            "__init__.py",
            "_core.py",
            "icons.json",
            "py.typed",
        ]],
        license = "MIT",
        version = faicons_version,
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
        exported_files = [site_packages + "/shiny_semantic/" + f for f in [
            "__init__.py",
            "_utils.py",
            "page.py",
            "collections/__init__.py",
            "elements/__init__.py",
            "elements/button.py",
            "elements/container.py",
            "elements/divider.py",
            "elements/emoji.py",
            "elements/flag.py",
            "elements/header.py",
            "elements/icon.py",
            "elements/input.py",
            "elements/segment.py",
            "modules/__init__.py",
            "modules/checkbox.py",
            "modules/dropdown.py",
            "modules/modal.py",
            "modules/slider.py",
            "types/__init__.py",
            "views/__init__.py",
            "views/statistic.py",
            "www/shiny-semantic-bindings.js",
            "www/bindings/semanticButton.js",
            "www/bindings/semanticCheckbox.js",
            "www/bindings/semanticCheckboxGroup.js",
            "www/bindings/semanticDropdown.js",
            "www/bindings/semanticModal.js",
            "www/bindings/semanticSlider.js",
            "www/semantic/fomantic.min.css",
            "www/semantic/fomantic.min.js",
            "www/semantic/themes/default/assets/fonts/Lato-Bold.woff2",
            "www/semantic/themes/default/assets/fonts/Lato-BoldItalic.woff2",
            "www/semantic/themes/default/assets/fonts/Lato-Italic.woff2",
            "www/semantic/themes/default/assets/fonts/Lato-Regular.woff2",
            "www/semantic/themes/default/assets/fonts/LatoLatin-Bold.woff2",
            "www/semantic/themes/default/assets/fonts/LatoLatin-BoldItalic.woff2",
            "www/semantic/themes/default/assets/fonts/LatoLatin-Italic.woff2",
            "www/semantic/themes/default/assets/fonts/LatoLatin-Regular.woff2",
            "www/semantic/themes/default/assets/fonts/brand-icons.woff2",
            "www/semantic/themes/default/assets/fonts/icons.woff2",
            "www/semantic/themes/default/assets/fonts/outline-icons.woff2",
        ]],
        license = "LGPL-3.0",
        version = shiny_semantic_version,
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
        exported_files = [site_packages + "/vega/" + f for f in [
            "__init__.py",
            "_frontend.py",
            "altair.py",
            "base.py",
            "utils.py",
            "vega.py",
            "vegalite.py",
            "widget.py",
            "static/extension.js",
            "static/extension.js.map",
            "static/index.js",
            "static/index.js.LICENSE.txt",
            "static/index.js.map",
            "static/labplugin.js",
            "static/labplugin.js.LICENSE.txt",
            "static/labplugin.js.map",
            "static/vega.js",
        ]],
        license = "BSD-3-Clause",
        version = vega_version,
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
        exported_files = [site_packages + "/fbpca.py"],
        version = "1.0",
        license = "BSD-3-Clause",
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
