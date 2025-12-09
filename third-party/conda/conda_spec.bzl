"""This file contains the workspace rules for the anaconda distribution.

To use, add
```
    load(":conda_spec.bzl", "anaconda_workspace")
    anaconda_workspace()
```
to your `WORKSPACE` file.
"""

load("@tenx_bazel_rules//rules:anaconda_repository.bzl", "anaconda_repository")
load(
    "@tenx_bazel_rules//rules:conda_package_repository.bzl",
    "conda_package_repository",
)

def anaconda_workspace(name = "anaconda"):
    """Create remote repositories to download each conda package.

    Also create the repository rule to generate the complete distribution.  In
    general, other targets should depend on targets in the `@anaconda`
    repository, rather than individual package repositories.

    Args:
        name (string): The name of the top level distribution repo.
    """
    conda_package_repository(
        name = "conda_package__python_abi3_support",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "_python_abi3_support-1.0-hd8ed1ab_2",
        sha256 = "a3967b937b9abf0f2a99f3173fa4630293979bd1644709d89580e7c62a544661",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_altair",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "altair-5.5.0-pyhd8ed1ab_1",
        sha256 = "74e60a5c0af8fa6f15a0e7860ad5f7b7c43c03a29b4ebba1433d24fc28029ebb",
        archive_type = "conda",
        extra_deps = ["pyarrow"],
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_attrs",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "attrs-25.4.0-pyh71513ae_0",
        sha256 = "f6c3c19fa599a1a856a88db166c318b148cac3ee4851a9905ed8a04eeec79f45",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_biopython",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "biopython-1.85-py312h4c3975b_2",
        sha256 = "cf6c4b111d1926524e1304e82bb78404ba08f82e83c63f5c16e8c3aca6dc7b57",
        exclude = ["lib/python*/site-packages/Bio/Entrez/DTDs"],
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_blosc",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "blosc-1.21.6-he440d0b_1",
        sha256 = "e7af5d1183b06a206192ff440e08db1c4e8b2ca1f8376ee45fb2f3a85d4ee45d",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_brotli_python",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "brotli-python-1.1.0-py312h1289d80_4",
        sha256 = "52a9ac412512b418ecdb364ba21c0f3dc96f0abbdb356b3cfbb980020b663d9b",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_bzip2",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "bzip2-1.0.8-hda65f42_8",
        sha256 = "c30daba32ddebbb7ded490f0e371eae90f51e72db620554089103b4a6934b0d5",
        exclude = ["man"],
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_c_ares",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "c-ares-1.34.5-hb9d3cd8_0",
        sha256 = "f8003bef369f57396593ccd03d08a8e21966157269426f71e943f96e4b579aeb",
        archive_type = "conda",
        exclude = ["share/man"],
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_c_blosc2",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "c-blosc2-2.21.3-h4cfbee9_0",
        sha256 = "4a222cff1b3507b289352ab94d110974dad3dace11e2d0eea405ba3147764eba",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_ca_certificates",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "ca-certificates-2025.10.5-hbd8a1cb_0",
        sha256 = "3b5ad78b8bb61b6cdc0978a6a99f8dfb2cc789a451378d054698441005ecbdb6",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_cached_property",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "cached_property-1.5.2-pyha770c72_1",
        sha256 = "6dbf7a5070cc43d90a1e4c2ec0c541c69d8e30a0e25f50ce9f6e4a432e42c5d7",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_certifi",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "certifi-2025.10.5-pyhd8ed1ab_0",
        sha256 = "955bac31be82592093f6bc006e09822cd13daf52b28643c9a6abd38cd5f4a306",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_cffi",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "cffi-2.0.0-py312h35888ee_0",
        sha256 = "f9e906b2cb9ae800b5818259472c3f781b14eb1952e867ac5c1f548e92bf02d9",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_charset_normalizer",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "charset-normalizer-3.4.3-pyhd8ed1ab_0",
        sha256 = "838d5a011f0e7422be6427becba3de743c78f3874ad2743c341accbba9bb2624",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_cpython",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "cpython-3.12.11-py312hd8ed1ab_0",
        sha256 = "7e7bc8e73a2f3736444a8564cbece7216464c00f0bc38e604b0c792ff60d621a",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_cython",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "cython-3.1.4-py312h7c45ced_2",
        sha256 = "7f240b106d117ecb6da873e3dd827b2efc2509b8403a0296d4d030f8daf77d6e",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_docopt",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "docopt-0.6.2-pyhd8ed1ab_2",
        sha256 = "7581a21e9bbe279d73d8ea32333f07ab286d2880edcee76a52480e2e4e53470d",
        exclude = ["lib/python3.*/site-packages/docopt-0.6.2.dist-info/METADATA"],
        conda_repo = name,
        archive_type = "conda",
    )
    conda_package_repository(
        name = "conda_package_h2",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "h2-4.3.0-pyhcf101f3_0",
        sha256 = "84c64443368f84b600bfecc529a1194a3b14c3656ee2e832d15a20e0329b6da3",
        conda_repo = name,
        archive_type = "conda",
    )
    conda_package_repository(
        name = "conda_package_h5py",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "h5py-3.14.0-nompi_py312ha4f8f14_101",
        sha256 = "6736b00b257aecef97e5e607ff275780cacdec48ff85963fe53abeb9ee4fb53f",
        exclude = ["lib/python3.*/site-packages/h5py/tests"],
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_hdf5",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "hdf5-1.14.6-nompi_h6e4c0c1_103",
        sha256 = "4f173af9e2299de7eee1af3d79e851bca28ee71e7426b377e841648b51d48614",
        exclude = ["share"],
        archive_type = "conda",
        conda_repo = name,
        exclude_deps = ["libgfortran"],
    )
    conda_package_repository(
        name = "conda_package_hpack",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "hpack-4.1.0-pyhd8ed1ab_0",
        sha256 = "6ad78a180576c706aabeb5b4c8ceb97c0cb25f1e112d76495bff23e3779948ba",
        conda_repo = name,
        archive_type = "conda",
    )
    conda_package_repository(
        name = "conda_package_htslib",
        base_urls = ["https://conda.anaconda.org/bioconda/linux-64"],
        dist_name = "htslib-1.22.1-h566b1c6_0",
        sha256 = "858e634fea447555dcf9725c6533ec7e0f345cb74b53d0d6c99f0ff9575924a1",
        exclude = ["share/man"],
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_hyperframe",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "hyperframe-6.1.0-pyhd8ed1ab_0",
        sha256 = "77af6f5fe8b62ca07d09ac60127a30d9069fdc3c68d6b256754d0ffb1f7779f8",
        conda_repo = name,
        archive_type = "conda",
    )
    conda_package_repository(
        name = "conda_package_icu",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "icu-75.1-he02047a_0",
        sha256 = "71e750d509f5fa3421087ba88ef9a7b9be11c53174af3aa4d06aff4c18b38e8e",
        archive_type = "conda",
        conda_repo = name,
        exclude = ["bin/gen*", "bin/makeconv", "bin/pkgdata", "lib/icu/*.inc"],
    )
    conda_package_repository(
        name = "conda_package_idna",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "idna-3.10-pyhd8ed1ab_1",
        sha256 = "d7a472c9fd479e2e8dcb83fb8d433fce971ea369d704ece380e876f9c3494e87",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_importlib_metadata",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "importlib-metadata-8.7.0-pyhe01879c_1",
        sha256 = "c18ab120a0613ada4391b15981d86ff777b5690ca461ea7e9e49531e8f374745",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_jinja2",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "jinja2-3.1.6-pyhd8ed1ab_0",
        sha256 = "f1ac18b11637ddadc05642e8185a851c7fab5998c6f5470d716812fae943b2af",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_joblib",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "joblib-1.5.2-pyhd8ed1ab_0",
        sha256 = "6fc414c5ae7289739c2ba75ff569b79f72e38991d61eb67426a8a4b92f90462c",
        exclude = ["site-packages/joblib/test"],
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_jsonschema",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "jsonschema-4.25.1-pyhe01879c_0",
        sha256 = "ac377ef7762e49cb9c4f985f1281eeff471e9adc3402526eea78e6ac6589cf1d",
        archive_type = "conda",
        exclude = ["site-packages/jsonschema/tests"],
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_jsonschema_specifications",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "jsonschema-specifications-2025.9.1-pyhcf101f3_0",
        sha256 = "0a4f3b132f0faca10c89fdf3b60e15abb62ded6fa80aebfc007d05965192aa04",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_keyutils",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "keyutils-1.6.3-hb9d3cd8_0",
        sha256 = "0960d06048a7185d3542d850986d807c6e37ca2e644342dd0c72feefcf26c2a4",
        exclude = ["share/man"],
        conda_repo = name,
        archive_type = "conda",
    )
    conda_package_repository(
        name = "conda_package_krb5",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "krb5-1.21.3-h659f571_0",
        sha256 = "99df692f7a8a5c27cd14b5fb1374ee55e756631b9c3d659ed3ee60830249b238",
        exclude = ["share/man", "share/examples"],
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_lcms2",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "lcms2-2.17-h717163a_0",
        sha256 = "d6a61830a354da022eae93fa896d0991385a875c6bba53c82263a289deda9db8",
        exclude = ["share/man"],
        conda_repo = name,
        archive_type = "conda",
    )
    conda_package_repository(
        name = "conda_package_lerc",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "lerc-4.0.0-h0aef613_1",
        sha256 = "412381a43d5ff9bbed82cd52a0bbca5b90623f62e41007c9c42d3870c60945ff",
        exclude = ["share/man"],
        conda_repo = name,
        archive_type = "conda",
    )
    conda_package_repository(
        name = "conda_package_libaec",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libaec-1.1.4-h3f801dc_0",
        sha256 = "410ab78fe89bc869d435de04c9ffa189598ac15bb0fe1ea8ace8fb1b860a2aa3",
        archive_type = "conda",
        exclude = ["share/man"],
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libblas",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libblas-3.9.0-20_linux64_mkl",
        sha256 = "9e5f27fca79223a5d38ccdf4c468e798c3684ba01bdb6b4b44e61f2103a298eb",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libcblas",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libcblas-3.9.0-20_linux64_mkl",
        sha256 = "841b4d44e20e5207f4a74ca98176629ead5ba590384ed6b0fe3c8600248c9fef",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libcurl",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libcurl-8.14.1-h332b0f4_0",
        sha256 = "b6c5cf340a4f80d70d64b3a29a7d9885a5918d16a5cb952022820e6d3e79dc8b",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libdeflate",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libdeflate-1.24-h86f0d12_0",
        sha256 = "8420748ea1cc5f18ecc5068b4f24c7a023cc9b20971c99c824ba10641fb95ddf",
        conda_repo = name,
        archive_type = "conda",
    )
    conda_package_repository(
        name = "conda_package_libedit",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libedit-3.1.20250104-pl5321h7949ede_0",
        sha256 = "d789471216e7aba3c184cd054ed61ce3f6dac6f87a50ec69291b9297f8c18724",
        exclude = ["share/man"],
        conda_repo = name,
        archive_type = "conda",
    )
    conda_package_repository(
        name = "conda_package_libev",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libev-4.33-hd590300_2",
        sha256 = "1cd6048169fa0395af74ed5d8f1716e22c19a81a8a36f934c110ca3ad4dd27b4",
        archive_type = "conda",
        exclude = ["share/man"],
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libexpat",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libexpat-2.7.1-hecca717_0",
        sha256 = "da2080da8f0288b95dd86765c801c6e166c4619b910b11f9a8446fb852438dc2",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libffi",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libffi-3.4.6-h2dba641_1",
        sha256 = "764432d32db45466e87f10621db5b74363a9f847d2b8b1f9743746cd160f06ab",
        licenses = ["@rules_license//licenses/spdx:MIT"],
        exclude = ["share/info", "share/man"],
        conda_repo = name,
        archive_type = "conda",
    )
    conda_package_repository(
        name = "conda_package_libfreetype",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libfreetype-2.14.1-ha770c72_0",
        sha256 = "4641d37faeb97cf8a121efafd6afd040904d4bca8c46798122f417c31d5dfbec",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libfreetype6",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libfreetype6-2.14.1-h73754d4_0",
        sha256 = "4a7af818a3179fafb6c91111752954e29d3a2a950259c14a2fc7ba40a8b03652",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libgcc",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libgcc-15.2.0-h767d61c_7",
        sha256 = "08f9b87578ab981c7713e4e6a7d935e40766e10691732bba376d4964562bcb45",
        archive_type = "conda",
        exclude_deps = ["_openmp_mutex"],
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libgcc_ng",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libgcc-ng-15.2.0-h69a702a_7",
        sha256 = "2045066dd8e6e58aaf5ae2b722fb6dfdbb57c862b5f34ac7bfb58c40ef39b6ad",
        exclude = [
            "lib/libgomp.so*",
            "x86_64-conda_cos6-linux-gnu/sysroot/lib/libgomp.so*",
            "share/info",
        ],
        exclude_deps = ["_openmp_mutex"],
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libgfortran5",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libgfortran5-15.2.0-hcd61629_7",
        sha256 = "e93ceda56498d98c9f94fedec3e2d00f717cbedfc97c49be0e5a5828802f2d34",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libhwloc",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libhwloc-2.12.1-default_h7f8ec31_1002",
        sha256 = "f7fbc792dbcd04bf27219c765c10c239937b34c6c1a1f77a5827724753e02da1",
        archive_type = "conda",
        exclude = ["share/doc", "share/man"],
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libiconv",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libiconv-1.18-h3b78370_2",
        sha256 = "c467851a7312765447155e071752d7bf9bf44d610a5687e32706f480aad2833f",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libjpeg_turbo",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libjpeg-turbo-3.1.0-hb9d3cd8_0",
        sha256 = "98b399287e27768bf79d48faba8a99a2289748c65cd342ca21033fab1860d4a4",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_liblapack",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "liblapack-3.9.0-20_linux64_mkl",
        sha256 = "21b4324dd65815f6b5a83c15f0b9a201434d0aa55eeecc37efce7ee70bbbf482",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_liblzma",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "liblzma-5.8.1-hb9d3cd8_2",
        sha256 = "f2591c0069447bbe28d4d696b7fcb0c5bd0b4ac582769b89addbcf26fb3430d8",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libnghttp2",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libnghttp2-1.67.0-had1ee68_0",
        sha256 = "a4a7dab8db4dc81c736e9a9b42bdfd97b087816e029e221380511960ac46c690",
        archive_type = "conda",
        exclude = ["share/doc", "share/man"],
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libnsl",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libnsl-2.0.1-hb9d3cd8_1",
        sha256 = "927fe72b054277cde6cb82597d0fcf6baf127dcbce2e0a9d8925a68f1265eef5",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libpng",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libpng-1.6.50-h421ea60_1",
        sha256 = "e75a2723000ce3a4b9fd9b9b9ce77553556c93e475a4657db6ed01abc02ea347",
        exclude = ["share/man"],
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libsqlite",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libsqlite-3.50.4-h0c1763c_0",
        sha256 = "6d9c32fc369af5a84875725f7ddfbfc2ace795c28f246dc70055a79f9b2003da",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libssh2",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libssh2-1.11.1-hcf80075_0",
        sha256 = "fa39bfd69228a13e553bd24601332b7cfeb30ca11a3ca50bb028108fe90a7661",
        exclude = ["share/doc", "share/man"],
        conda_repo = name,
        archive_type = "conda",
    )
    conda_package_repository(
        name = "conda_package_libstdcxx",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libstdcxx-15.2.0-h8f9b012_7",
        sha256 = "1b981647d9775e1cdeb2fab0a4dd9cd75a6b0de2963f6c3953dbd712f78334b3",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libstdcxx_ng",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libstdcxx-ng-15.2.0-h4852527_7",
        sha256 = "024fd46ac3ea8032a5ec3ea7b91c4c235701a8bf0e6520fe5e6539992a6bd05f",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libtiff",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libtiff-4.7.1-h8261f1e_0",
        sha256 = "ddda0d7ee67e71e904a452010c73e32da416806f5cb9145fb62c322f97e717fb",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libuuid",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libuuid-2.41.2-he9a06e4_0",
        sha256 = "e5ec6d2ad7eef538ddcb9ea62ad4346fde70a4736342c4ad87bd713641eb9808",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libwebp_base",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libwebp-base-1.6.0-hd42ef1d_0",
        sha256 = "3aed21ab28eddffdaf7f804f49be7a7d701e8f0e46c856d801270b470820a37b",
        exclude = ["share/man"],
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libxcb",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libxcb-1.17.0-h8a09558_0",
        sha256 = "666c0c431b23c6cec6e492840b176dde533d48b7e6fb8883f5071223433776aa",
        conda_repo = name,
        archive_type = "conda",
    )
    conda_package_repository(
        name = "conda_package_libxcrypt",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libxcrypt-4.4.36-hd590300_1",
        sha256 = "6ae68e0b86423ef188196fff6207ed0c8195dd84273cb5623b85aa08033a410c",
        archive_type = "conda",
        exclude = ["share/man"],
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libxml2_16",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libxml2-16-2.15.0-ha9997c6_1",
        sha256 = "5420ea77505a8d5ca7b5351ddb2da7e8a178052fccf8fca00189af7877608e89",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libxml2",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libxml2-2.15.0-h26afc86_1",
        sha256 = "4310577d7eea817d35a1c05e1e54575b06ce085d73e6dd59aa38523adf50168f",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libzlib",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libzlib-1.3.1-hb9d3cd8_2",
        sha256 = "d4bfe88d7cb447768e31650f06257995601f89076080e76df55e3112d4e47dc4",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_llvm_openmp",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "llvm-openmp-21.1.2-h4922eb0_3",
        sha256 = "2b8d157370cb9202d4970a2353a02517ccf72e81f2d95920570aef934d0508fd",
        archive_type = "conda",
        exclude = ["lib/libomptarget.rtl.*"],
        exclude_deps = ["zstd", "libzlib"],
        patch_cmds = ["ln -s libomp.so lib/libgomp.so.1"],
        patch_args = ["-p1"],
        patches = ["@cellranger//third-party/conda:intel-openmp-symlinks.patch"],
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_lz4",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "lz4-4.4.4-py312h5d89b6d_1",
        sha256 = "672bd94e67feff49461b7eb7a3ca08100681ebf76456e1f98fa0f08b17a04d2e",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_lz4_c",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "lz4-c-1.10.0-h5888daf_1",
        sha256 = "47326f811392a5fd3055f0f773036c392d26fdb32e4d8e7a8197eed951489346",
        exclude = ["share/man"],
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_markupsafe",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "markupsafe-3.0.3-py312h8a5da7c_0",
        sha256 = "f77f9f1a4da45cbc8792d16b41b6f169f649651a68afdc10b2da9da12b9aa42b",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_mkl",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "mkl-2023.2.0-ha770c72_50498",
        sha256 = "cafca5b96c1431757909fa11942cb57056870d3fede2eb6bbc138f5ed86679a4",
        exclude = [
            "lib/libmkl_*avx512_mic.so.1",
            "lib/libmkl_blacs_*mpi_*.so*",
            "lib/libmkl_blacs_sgimpt*.so*",
            "lib/libmkl_pgi_thread.so*",
            "lib/libmkl_scalapack*.so*",
            "lib/libmkl_*_ilp64.so*",
        ],
        patch_args = ["-p1"],
        patches = ["@cellranger//third-party/conda:mkl-libblas-symlinks.patch"],
        patch_cmds = [
            "ln -s libmkl_rt.so.2 lib/libblas.so",
            "ln -s libmkl_rt.so.2 lib/libblas.so.3",
            "ln -s libmkl_rt.so.2 lib/libcblas.so",
            "ln -s libmkl_rt.so.2 lib/libcblas.so.3",
            "ln -s libmkl_rt.so.2 lib/liblapack.so",
            "ln -s libmkl_rt.so.2 lib/liblapack.so.3",
        ],
        exclude_deps = ["_openmp_mutex"],
        extra_deps = ["llvm-openmp"],
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_narwhals",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "narwhals-2.7.0-pyhcf101f3_0",
        sha256 = "b377b79c37a9c42b51b1d701adae968f853f82f7bcce5252e4ccaf56a03f943c",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_ncurses",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "ncurses-6.5-h2d0b736_3",
        sha256 = "3fde293232fa3fca98635e1167de6b7c7fda83caf24b9d6c91ec9eefb4f4d586",
        licenses = ["@rules_license//licenses/spdx:MIT-open-group"],
        exclude = [
            "share/terminfo/[A-Z]",
            "share/terminfo/d/darwin*",
            "share/terminfo/n/n7900",
            "share/terminfo/n/ncr*",
        ],
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_numexpr",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "numexpr-2.10.1-mkl_py312h791fadb_2",
        sha256 = "e36e16a29da93bf1f18d58a49f0319b3ebed91bc994b2e8ee6d4cb59218eea61",
        exclude = ["lib/python*/site-packages/numexpr/tests"],
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_numpy",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "numpy-2.3.3-py312h33ff503_0",
        sha256 = "8443315a60500ea8e3d7ecd9756cee07a60b8c3497e0fc98884963c3108f8bef",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_openjpeg",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "openjpeg-2.5.4-h55fea9a_0",
        sha256 = "3900f9f2dbbf4129cf3ad6acf4e4b6f7101390b53843591c53b00f034343bc4d",
        conda_repo = name,
        archive_type = "conda",
    )
    conda_package_repository(
        name = "conda_package_openssl",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "openssl-3.5.4-h26f9b46_0",
        sha256 = "e807f3bad09bdf4075dbb4168619e14b0c0360bacb2e12ef18641a834c8c5549",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_packaging",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "packaging-25.0-pyh29332c3_1",
        sha256 = "289861ed0c13a15d7bbb408796af4de72c2fe67e2bcb0de98f4c3fce259d7991",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_pandas",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "pandas-2.2.3-py312hf9745cd_1",
        sha256 = "ad275a83bfebfa8a8fee9b0569aaf6f513ada6a246b2f5d5b85903d8ca61887e",
        exclude = ["lib/python3.*/site-packages/pandas/tests"],
        extra_deps = ["pytables"],
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_patsy",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "patsy-1.0.1-pyhd8ed1ab_1",
        sha256 = "ab52916f056b435757d46d4ce0a93fd73af47df9c11fd72b74cc4b7e1caca563",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_pillow",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "pillow-11.3.0-py312h7b42cdd_3",
        sha256 = "ad4a22899819a2bb86550d1fc3833a44e073aac80ea61529676b5e73220fcc2b",
        conda_repo = name,
        archive_type = "conda",
    )
    conda_package_repository(
        name = "conda_package_plotly",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "plotly-6.3.1-pyhd8ed1ab_0",
        sha256 = "077f4a0ebf4d0408a9cd04c6e41b4d0f6bbc2daa5ac62124df0dc6c1c057ed91",
        license_file = "info/licenses/LICENSE0.txt",
        archive_type = "conda",
        exclude = [
            "lib/python3.*/site-packages/plotly/graph_objs/*/__pycache__",
            "lib/python3.*/site-packages/plotly/graph_objs/*/*/__pycache__",
            "lib/python3.*/site-packages/plotly/graph_objs/*/*/*/__pycache__",
            "lib/python3.*/site-packages/plotly/graph_objs/*/*/*/*/__pycache__",
            "lib/python3.*/site-packages/plotly/validators/*/__pycache__",
            "lib/python3.*/site-packages/plotly/validators/*/*/__pycache__",
            "lib/python3.*/site-packages/plotly/validators/*/*/*/__pycache__",
            "lib/python3.*/site-packages/plotly/validators/*/*/*/*/__pycache__",
            "lib/python3.*/python-scripts",
        ],
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_polars_lts_cpu",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "polars-lts-cpu-1.33.1-py310hee84464_1",
        sha256 = "f9c174305a882bf07d76f565d1b786ea6ed000f25fdd044e0a5b14bcdf8b1414",
        archive_type = "conda",
        extra_deps = ["pyarrow"],
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_pthread_stubs",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "pthread-stubs-0.4-hb9d3cd8_1002",
        sha256 = "9c88f8c64590e9567c6c80823f0328e58d3b1efb0e1c539c0315ceca764e0973",
        conda_repo = name,
        archive_type = "conda",
    )
    conda_package_repository(
        name = "conda_package_py_cpuinfo",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "py-cpuinfo-9.0.0-pyhd8ed1ab_1",
        sha256 = "6d8f03c13d085a569fde931892cded813474acbef2e03381a1a87f420c7da035",
        conda_repo = name,
        archive_type = "conda",
    )
    conda_package_repository(
        name = "conda_package_pycparser",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "pycparser-2.22-pyh29332c3_1",
        sha256 = "79db7928d13fab2d892592223d7570f5061c192f27b9febd1a418427b719acc6",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_pyfaidx",
        base_urls = ["https://conda.anaconda.org/bioconda/noarch"],
        dist_name = "pyfaidx-0.9.0.3-pyhdfd78af_0",
        sha256 = "6bc1f9da77cf3d1bffaf1a4bc04f56d3d6d8f0566f1638a5e6a4a62160594bdd",
        conda_repo = name,
        archive_type = "conda",
    )
    conda_package_repository(
        name = "conda_package_pysam",
        base_urls = ["https://conda.anaconda.org/bioconda/linux-64"],
        dist_name = "pysam-0.23.3-py312h47d5410_1",
        sha256 = "435029af0a0e9c211974f16f607f6c326d5b3309bc1d80b42bff9b97e06441c6",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_pysocks",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "pysocks-1.7.1-pyha55dd90_7",
        sha256 = "ba3b032fa52709ce0d9fd388f63d330a026754587a2f461117cac9ab73d8d0d8",
        conda_repo = name,
        archive_type = "conda",
    )
    conda_package_repository(
        name = "conda_package_pytables",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "pytables-3.10.2-py312hefc0c3f_9",
        sha256 = "a0865975778fcc7da8f60ed469c072fd919a04b11e3c92d3460d7e733a273720",
        patch_args = ["-p1"],
        patches = ["@cellranger//third-party/conda:pytables-3.10.2-perms.patch"],
        exclude = [
            "lib/python3.*/site-packages/tables/tests/*.h5",
            "lib/python3.*/site-packages/tables/tests/*.mat",
            "lib/python3.*/site-packages/tables/tests/*/*",
        ],
        exclude_deps = ["lzo"],
        extra_deps = ["setuptools"],
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_python",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "python-3.12.11-h9e4cc4f_0_cpython",
        sha256 = "6cca004806ceceea9585d4d655059e951152fc774a471593d4f5138e6a54c81d",
        exclude = [
            "lib/python*/lib2to3",
            "lib/python*/distutils/tests",
            "lib/python*/ensurepip",
            "lib/python*/idlelib/idle_test",
            "lib/python*/test",
            "doc",
            "share/man",
        ],
        exclude_deps = ["ld_impl_linux-64", "readline", "tk"],
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_python_dateutil",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "python-dateutil-2.9.0.post0-pyhe01879c_2",
        sha256 = "d6a17ece93bbd5139e02d2bd7dbfa80bee1a4261dced63f65f679121686bf664",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_python_gil",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "python-gil-3.12.11-hd8ed1ab_0",
        sha256 = "b8afeaefe409d61fa4b68513b25a66bb17f3ca430d67cfea51083c7bfbe098ef",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_python_tzdata",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "python-tzdata-2025.2-pyhd8ed1ab_0",
        sha256 = "e8392a8044d56ad017c08fec2b0eb10ae3d1235ac967d0aab8bd7b41c4a5eaf0",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_python_abi",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "python_abi-3.12-8_cp312",
        sha256 = "80677180dd3c22deb7426ca89d6203f1c7f1f256f2d5a94dc210f6e758229809",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_pytz",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "pytz-2024.1-pyhd8ed1ab_0",
        sha256 = "1a7d6b233f7e6e3bbcbad054c8fd51e690a67b129a899a056a5e45dd9f00cb41",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_pyvcf3",
        base_urls = ["https://conda.anaconda.org/bioconda/linux-64"],
        dist_name = "pyvcf3-1.0.4-py312h0fa9677_0",
        sha256 = "998f6d9663fc6549ab4f9d9a8c814ad83f8dc3a7fd33b0080f3a3b15c519e3f4",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_referencing",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "referencing-0.36.2-pyh29332c3_0",
        sha256 = "e20909f474a6cece176dfc0dc1addac265deb5fa92ea90e975fbca48085b20c3",
        exclude = ["site-packages/referencing/tests"],
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_requests",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "requests-2.32.5-pyhd8ed1ab_0",
        sha256 = "8dc54e94721e9ab545d7234aa5192b74102263d3e704e6d0c8aa7008f2da2a7b",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_rpds_py",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "rpds-py-0.27.1-py312h868fb18_1",
        sha256 = "76efba673e02d4d47bc2de6e48a8787ed98bae4933233dee5ce810fa3de6ef2b",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_scikit_learn",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "scikit-learn-1.7.2-py312h4f0b9e3_0",
        sha256 = "27e2f65075556804a7524dc6bcc34829602f986621728100c0ef07b404168aa8",
        exclude = [
            "lib/python*/site-packages/sklearn/*/*/tests",
            "lib/python*/site-packages/sklearn/*/tests",
            "lib/python*/site-packages/sklearn/tests",
        ],
        archive_type = "conda",
        exclude_deps = ["_openmp_mutex"],
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_scipy",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "scipy-1.15.2-py312ha707e6e_0",
        sha256 = "b9faaa024b77a3678a988c5a490f02c4029c0d5903998b585100e05bc7d4ff36",
        exclude = [
            "lib/python*/site-packages/scipy/linalg/src/id_dist/doc",
            "lib/python*/site-packages/scipy/*/*/*/*/tests",
            "lib/python*/site-packages/scipy/*/*/*/tests",
            "lib/python*/site-packages/scipy/*/*/tests",
            "lib/python*/site-packages/scipy/*/tests",
        ],
        archive_type = "conda",
        exclude_deps = ["libgfortran"],
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_setuptools",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "setuptools-80.9.0-pyhff2d567_0",
        sha256 = "972560fcf9657058e3e1f97186cc94389144b46dbdf58c807ce62e83f977e863",
        exclude = [
            "lib/python*/site-packages/setuptools/command/launcher manifest.xml",
            "lib/python*/site-packages/setuptools/script (dev).tmpl",
            "lib/python*/site-packages/setuptools/*.exe",
            "lib/python*/site-packages/pkg_resources/tests",
            "site-packages/setuptools/command/launcher manifest.xml",
            "site-packages/setuptools/script (dev).tmpl",
            "site-packages/setuptools/*.exe",
            "site-packages/pkg_resources/tests",
        ],
        patch_args = ["-p1"],
        patches = ["@cellranger//third-party/conda:setuptools-71.0.4-rename.patch"],
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_six",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "six-1.17.0-pyhe01879c_1",
        sha256 = "458227f759d5e3fcec5d9b7acce54e10c9e1f4f4b7ec978f3bfd54ce4ee9853d",
        conda_repo = name,
        archive_type = "conda",
    )
    conda_package_repository(
        name = "conda_package_snappy",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "snappy-1.2.2-h03e3b7b_0",
        sha256 = "8b8acbde6814d1643da509e11afeb6bb30eb1e3004cf04a7c9ae43e9b097f063",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_statsmodels",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "statsmodels-0.14.5-py312h4f23490_1",
        sha256 = "5ea409c808240d4aafd4b90f2b881cab0ae411187caed90008a843748723f115",
        exclude = [
            "lib/python*/site-packages/statsmodels/*/tests",
            "lib/python*/site-packages/statsmodels/*/*/tests",
            "lib/python*/site-packages/statsmodels/*/*/*/tests",
            "lib/python*/site-packages/statsmodels/tests",
            "lib/python*/site-packages/statsmodels/sandbox/distributions/examples",
        ],
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_tbb",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "tbb-2021.13.0-hb60516a_3",
        sha256 = "cf9101d1327de410a844f29463c486c47dfde506d0c0656d2716c03135666c3f",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_threadpoolctl",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "threadpoolctl-3.6.0-pyhecae5ae_0",
        sha256 = "6016672e0e72c4cf23c0cf7b1986283bd86a9c17e8d319212d78d8e9ae42fdfd",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_tk",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "tk-8.6.13-noxft_hd72426e_102",
        sha256 = "a84ff687119e6d8752346d1d408d5cf360dee0badd487a472aa8ddedfdc219e1",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_typing_extensions",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "typing_extensions-4.15.0-pyhcf101f3_0",
        sha256 = "032271135bca55aeb156cee361c81350c6f3fb203f57d024d7e5a1fc9ef18731",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_tzdata",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "tzdata-2025b-h78e105d_0",
        sha256 = "5aaa366385d716557e365f0a4e9c3fca43ba196872abbbe3d56bb610d131e192",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_urllib3",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "urllib3-2.5.0-pyhd8ed1ab_0",
        sha256 = "4fb9789154bd666ca74e428d973df81087a697dbb987775bc3198d2215f240f8",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_xorg_libxau",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "xorg-libxau-1.0.12-hb9d3cd8_0",
        sha256 = "ed10c9283974d311855ae08a16dfd7e56241fac632aec3b92e3cfe73cff31038",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_xorg_libxdmcp",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "xorg-libxdmcp-1.1.5-hb9d3cd8_0",
        sha256 = "6b250f3e59db07c2514057944a3ea2044d6a8cdde8a47b6497c254520fade1ee",
        conda_repo = name,
        archive_type = "conda",
    )
    conda_package_repository(
        name = "conda_package_zipp",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "zipp-3.23.0-pyhd8ed1ab_0",
        sha256 = "7560d21e1b021fd40b65bfb72f67945a3fcb83d78ad7ccf37b8b3165ec3b68ad",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_zlib_ng",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "zlib-ng-2.2.5-hde8ca8f_0",
        sha256 = "3a8e7798deafd0722b6b5da50c36b7f361a80b30165d600f7760d569a162ff95",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_zstandard",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "zstandard-0.25.0-py312h5253ce2_0",
        sha256 = "1a3beda8068b55639edb92da8e0dc2d487e2a11aba627f709aab1d3cd5dd271c",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_zstd",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "zstd-1.5.7-hb8e6e7a_2",
        sha256 = "a4166e3d8ff4e35932510aaff7aa90772f84b4d07e9f6f83c614cba7ceefe0eb",
        exclude = ["share/man"],
        licenses = ["@rules_license//licenses/spdx:BSD-3-Clause"],
        archive_type = "conda",
        conda_repo = name,
    )
    anaconda_repository(
        name = name,
        conda_packages = [
            "_python_abi3_support",
            "altair",
            "attrs",
            "biopython",
            "blosc",
            "bzip2",
            "cached_property",
            "certifi",
            "cffi",
            "cpython",
            "cython",
            "docopt",
            "h2",
            "h5py",
            "hdf5",
            "hpack",
            "htslib",
            "hyperframe",
            "icu",
            "idna",
            "jinja2",
            "joblib",
            "jsonschema",
            "keyutils",
            "krb5",
            "lcms2",
            "lerc",
            "libaec",
            "libblas",
            "libcblas",
            "libcurl",
            "libdeflate",
            "libedit",
            "libev",
            "libexpat",
            "libffi",
            "libfreetype",
            "libfreetype6",
            "libgcc",
            "libgfortran5",
            "libhwloc",
            "libiconv",
            "liblapack",
            "liblzma",
            "libnghttp2",
            "libnsl",
            "libpng",
            "libsqlite",
            "libssh2",
            "libstdcxx",
            "libtiff",
            "libuuid",
            "libxcb",
            "libxcrypt",
            "libxml2",
            "libzlib",
            "lz4",
            "markupsafe",
            "mkl",
            "narwhals",
            "ncurses",
            "numexpr",
            "numpy",
            "openjpeg",
            "openssl",
            "packaging",
            "pandas",
            "patsy",
            "pillow",
            "plotly",
            "pycparser",
            "pyfaidx",
            "pysam",
            "pysocks",
            "pytables",
            "python",
            "python_abi",
            "pytz",
            "pyvcf3",
            "referencing",
            "requests",
            "scipy",
            "setuptools",
            "six",
            "snappy",
            "statsmodels",
            "tbb",
            "threadpoolctl",
            "tk",
            "typing_extensions",
            "tzdata",
            "urllib3",
            "zipp",
            "zstandard",
            "zstd",
            "pyarrow",
        ],
        executable_packages = [
            "black",
            "coverage",
            "cython",
            "fastcat",
            "ipython",
            "jupyterlab",
            "notebook",
            "python",
            "ruff",
            "shellcheck",
        ],
        aliases = {
            "bh_sne_3d": "tsne",
            "bh_sne": "tsne",
            "cached-property": "cached_property",
            "Bio": "biopython",
            "importlib_metadata": "importlib-metadata",
            "lazy_loader": "lazy-loader",
            "matplotlib": "matplotlib-base",
            "melting": "melt",
            "MOODS": "moods",
            "numpy-base": "numpy",
            "numpy_base": "numpy",
            "parasail": "parasail-python",
            "pigeon": "pbpigeon",
            "skera": "pbskera",
            "PIL": "pillow",
            "polars": "polars-lts-cpu",
            "primer3": "primer3-py",
            "progressbar": "progressbar2",
            "prompt_toolkit": "prompt-toolkit",
            "pyBigWig": "pybigwig",
            "SimpleITK": "simpleitk",
            "skimage": "scikit-image",
            "sklearn": "scikit-learn",
            "tables": "pytables",
            "typing-extensions": "typing_extensions",
            "umap": "umap-learn",
        },
        py_version = 3,
    )
    native.register_toolchains("@{}//:python_toolchain".format(name))
