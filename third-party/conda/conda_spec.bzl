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
        name = "conda_package__libgcc_mutex",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "_libgcc_mutex-0.1-conda_forge",
        sha256 = "fe51de6107f9edc7aa4f786a70f4a883943bc9d39b3bb7307c04c41410990726",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_altair",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "altair-5.4.1-pyhd8ed1ab_3",
        sha256 = "fc20cce4f91eeaf7c4f801f635f6fb5f2d29e10ac373ffa3292a4407f9ef73b0",
        archive_type = "conda",
        extra_deps = ["pyarrow"],
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_attrs",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "attrs-24.3.0-pyh71513ae_0",
        sha256 = "750186af694a7130eaf7119fbb56db0d2326d8995ad5b8eae23c622b85fea29a",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_biopython",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "biopython-1.84-py310hc51659f_0",
        sha256 = "cc686ead6074272cc5f1f137eee1fd746b60cc4fbcef070a4a4f2efa0935040b",
        exclude = ["lib/python*/site-packages/Bio/Entrez/DTDs"],
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_blosc",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "blosc-1.21.5-h0f2a231_0",
        sha256 = "e2b15b017775d1bda8edbb1bc48e545e45364edefa4d926732fc5488cc600731",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_brotli_python",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "brotli-python-1.0.9-py310hd8f1fbe_9",
        sha256 = "a4984cb906910850ae979387f0ac4e2623e0a3c8139bc67eb1fff0403bf8388b",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_bzip2",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "bzip2-1.0.8-h4bc722e_7",
        sha256 = "5ced96500d945fb286c9c838e54fa759aa04a7129c59800f0846b4335cee770d",
        exclude = ["man"],
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_c_ares",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "c-ares-1.34.4-hb9d3cd8_0",
        sha256 = "d4f28d87b6339b94f74762c0076e29c8ef8ddfff51a564a92da2843573c18320",
        archive_type = "conda",
        exclude = ["share/man"],
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_c_blosc2",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "c-blosc2-2.12.0-hb4ffafa_0",
        sha256 = "68ae377f7baeb616e5a24facadebd8bf7a9cd48a297124be9d814ba92ff5e40f",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_ca_certificates",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "ca-certificates-2024.12.14-hbcca054_0",
        sha256 = "1afd7274cbc9a334d6d0bc62fa760acc7afdaceb0b91a8df370ec01fd75dc7dd",
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
        dist_name = "certifi-2024.12.14-pyhd8ed1ab_0",
        sha256 = "048c16a9cbcb1fbad02083414d3bc7c1d0eea4b39aee6aa6bf8d1d5089ca8bad",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_cffi",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "cffi-1.17.1-py310h8deb56e_0",
        sha256 = "1b389293670268ab80c3b8735bc61bc71366862953e000efbb82204d00e41b6c",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_charset_normalizer",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "charset-normalizer-3.4.1-pyhd8ed1ab_0",
        sha256 = "4e0ee91b97e5de3e74567bdacea27f0139709fceca4db8adffbe24deffccb09b",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_cython",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "cython-3.0.11-py310h5b1441d_3",
        sha256 = "ab2fc6b4b8c203da4cade3c25334342dc9f3ce846a2ae81591c3d09bf7b7ed4d",
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
        name = "conda_package_freetype",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "freetype-2.12.1-h267a509_2",
        sha256 = "b2e3c449ec9d907dd4656cb0dc93e140f447175b125a3824b31368b06c666bb6",
        licenses = ["@rules_license//licenses/spdx:FTL"],
        license_file = "info/licenses/docs/FTL.TXT",
        exclude = ["share/man"],
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_h2",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "h2-4.1.0-pyhd8ed1ab_1",
        sha256 = "843ddad410c370672a8250470697027618f104153612439076d4d7b91eeb7b5c",
        conda_repo = name,
        archive_type = "conda",
    )
    conda_package_repository(
        name = "conda_package_h5py",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "h5py-3.9.0-nompi_py310h367e799_100",
        sha256 = "6a8a5b595032eff60df322492b865640dba7392348afe318756855f98de94792",
        exclude = ["lib/python3.*/site-packages/h5py/tests"],
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_hdf5",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "hdf5-1.14.0-nompi_h5231ba7_103",
        sha256 = "cfec313dcb1c51a0f6faee9c6af7dcb9cbe86847697f0be9846eaa0058500d29",
        exclude = ["share"],
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_hpack",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "hpack-4.0.0-pyhd8ed1ab_1",
        sha256 = "ec89b7e5b8aa2f0219f666084446e1fb7b54545861e9caa892acb24d125761b5",
        conda_repo = name,
        archive_type = "conda",
    )
    conda_package_repository(
        name = "conda_package_htslib",
        base_urls = ["https://conda.anaconda.org/bioconda/linux-64"],
        dist_name = "htslib-1.17-h6bc39ce_1",
        sha256 = "3236e98be603dabec298c15e493445e2b6e9059648e951b595f032e8cb994cb6",
        exclude = ["share/man"],
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_hyperframe",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "hyperframe-6.0.1-pyhd8ed1ab_1",
        sha256 = "e91c6ef09d076e1d9a02819cd00fa7ee18ecf30cdd667605c853980216584d1b",
        conda_repo = name,
        archive_type = "conda",
    )
    conda_package_repository(
        name = "conda_package_icu",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "icu-73.2-h59595ed_0",
        sha256 = "e12fd90ef6601da2875ebc432452590bc82a893041473bc1c13ef29001a73ea8",
        archive_type = "conda",
        conda_repo = name,
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
        dist_name = "importlib-metadata-8.5.0-pyha770c72_1",
        sha256 = "13766b88fc5b23581530d3a0287c0c58ad82f60401afefab283bf158d2be55a9",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_importlib_resources",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "importlib_resources-6.4.5-pyhd8ed1ab_1",
        sha256 = "461199e429a3db01f0a673f8beaac5e0be75b88895952fb9183f2ab01c5c3c24",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_jinja2",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "jinja2-3.1.5-pyhd8ed1ab_0",
        sha256 = "98977694b9ecaa3218662f843425f39501f81973c450f995eec68f1803ed71c3",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_joblib",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "joblib-1.4.2-pyhd8ed1ab_1",
        sha256 = "51cc2dc491668af0c4d9299b0ab750f16ccf413ec5e2391b924108c1fbacae9b",
        exclude = ["site-packages/joblib/test"],
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_jpeg",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "jpeg-9e-h0b41bf4_3",
        sha256 = "8f73194d09c9ea4a7e2b3562766b8d72125cc147b62c7cf83393e3a3bbfd581b",
        archive_type = "conda",
        licenses = ["@rules_license//licenses/spdx:IJG"],
        exclude = ["share/man"],
        license_file = "info/licenses/LICENSE0.txt",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_jsonschema",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "jsonschema-4.23.0-pyhd8ed1ab_1",
        sha256 = "be992a99e589146f229c58fe5083e0b60551d774511c494f91fe011931bd7893",
        archive_type = "conda",
        exclude = ["site-packages/jsonschema/tests"],
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_jsonschema_specifications",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "jsonschema-specifications-2024.10.1-pyhd8ed1ab_1",
        sha256 = "37127133837444cf0e6d1a95ff5a505f8214ed4e89e8e9343284840e674c6891",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_keyutils",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "keyutils-1.6.1-h166bdaf_0",
        sha256 = "150c05a6e538610ca7c43beb3a40d65c90537497a4f6a5f4d15ec0451b6f5ebb",
        exclude = ["share/man"],
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_krb5",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "krb5-1.20.1-hf9c8cef_0",
        sha256 = "3274ef26e40df1c23bd34adc075e40cc6c335540688c36ac261d64561da56278",
        exclude = ["share/man", "share/examples"],
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_lcms2",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "lcms2-2.14-h6ed2654_0",
        sha256 = "cbadb4150850941bf0518ba948effbbdd89b2c28dfdfed54eae196037e015b43",
        exclude = ["share/man"],
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_lerc",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "lerc-4.0.0-h27087fc_0",
        sha256 = "cb55f36dcd898203927133280ae1dc643368af041a48bcf7c026acb7c47b0c12",
        exclude = ["share/man"],
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libaec",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libaec-1.1.3-h59595ed_0",
        sha256 = "2ef420a655528bca9d269086cf33b7e90d2f54ad941b437fb1ed5eca87cee017",
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
        dist_name = "libcurl-7.87.0-h6312ad2_0",
        sha256 = "4e95c12244a50c8f8e9173e0bd37d6067fd753437ab636afb44ce28382a325eb",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libdeflate",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libdeflate-1.14-h166bdaf_0",
        sha256 = "6f7cbc9347964e7f9697bde98a8fb68e0ed926888b3116474b1224eaa92209dc",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libedit",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libedit-3.1.20191231-he28a2e2_2",
        sha256 = "a57d37c236d8f7c886e01656f4949d9dcca131d2a0728609c6f7fa338b65f1cf",
        exclude = ["share/man"],
        conda_repo = name,
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
        name = "conda_package_libffi",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libffi-3.4.2-h7f98852_5",
        sha256 = "ab6e9856c21709b7b517e940ae7028ae0737546122f83c2aa5d692860c3b149e",
        licenses = ["@rules_license//licenses/spdx:MIT"],
        exclude = ["share/info", "share/man"],
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libgcc",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libgcc-14.2.0-h77fa898_1",
        sha256 = "53eb8a79365e58849e7b1a068d31f4f9e718dc938d6f2c03e960345739a03569",
        archive_type = "conda",
        exclude_deps = ["_openmp_mutex"],
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libgcc_ng",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libgcc-ng-14.2.0-h69a702a_1",
        sha256 = "3a76969c80e9af8b6e7a55090088bc41da4cffcde9e2c71b17f44d37b7cb87f7",
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
        name = "conda_package_libgfortran_ng",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libgfortran-ng-14.2.0-h69a702a_1",
        sha256 = "423f1e2403f0c665748e42d335e421e53fd03c08d457cfb6f360d329d9459851",
        archive_type = "conda",
        exclude_deps = ["libgfortran"],
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libgfortran5",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libgfortran5-14.2.0-hd5240d6_1",
        sha256 = "d149a37ca73611e425041f33b9d8dbed6e52ec506fe8cc1fc0ee054bddeb6d5d",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libhwloc",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libhwloc-2.11.2-default_he43201b_1000",
        sha256 = "75be8732e6f94ff2faa129f44ec4970275e1d977559b0c2fb75b7baa5347e16b",
        archive_type = "conda",
        exclude = ["share/doc", "share/man"],
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libiconv",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libiconv-1.17-hd590300_2",
        sha256 = "8ac2f6a9f186e76539439e50505d98581472fedb347a20e7d1f36429849f05c9",
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
        name = "conda_package_libnghttp2",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libnghttp2-1.51.0-hdcd2b5c_0",
        sha256 = "3f76e99eacfc4ce3deac3b78d5508449efb8b72dea3e31d9f2c6db7f5cf00e75",
        archive_type = "conda",
        exclude = ["share/doc", "share/man"],
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libnsl",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libnsl-2.0.1-hd590300_0",
        sha256 = "26d77a3bb4dceeedc2a41bd688564fe71bf2d149fdcf117049970bc02ff1add6",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libpng",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libpng-1.6.43-h2797004_0",
        sha256 = "502f6ff148ac2777cc55ae4ade01a8fc3543b4ffab25c4e0eaa15f94e90dd997",
        exclude = ["share/man"],
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libsqlite",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libsqlite-3.46.0-hde9e2c9_0",
        sha256 = "daee3f68786231dad457d0dfde3f7f1f9a7f2018adabdbb864226775101341a8",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libssh2",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libssh2-1.10.0-haa6b8db_3",
        sha256 = "3c2ed83502bedf4ec8c5b972accb6ff1b6c018f72fb711cdb65cb8540d5ab89e",
        exclude = ["share/doc", "share/man"],
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libstdcxx",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libstdcxx-14.2.0-hc0a3c3a_1",
        sha256 = "4661af0eb9bdcbb5fb33e5d0023b001ad4be828fccdcc56500059d56f9869462",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libstdcxx_ng",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libstdcxx-ng-14.2.0-h4852527_1",
        sha256 = "25bb30b827d4f6d6f0522cc0579e431695503822f144043b93c50237017fffd8",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libtiff",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libtiff-4.4.0-h82bc61c_5",
        sha256 = "f81d38e7458c6ba2fcf93bef4ed2c12c2977e89ca9a7f936ce53a3338a88352f",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libuuid",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libuuid-2.38.1-h0b41bf4_0",
        sha256 = "787eb542f055a2b3de553614b25f09eefb0a0931b0c87dbcce6efdfd92f04f18",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libwebp_base",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libwebp-base-1.5.0-h851e524_0",
        sha256 = "c45283fd3e90df5f0bd3dbcd31f59cdd2b001d424cf30a07223655413b158eaf",
        exclude = ["share/man"],
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libxcb",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libxcb-1.13-h7f98852_1004",
        sha256 = "8d5d24cbeda9282dd707edd3156e5fde2e3f3fe86c802fa7ce08c8f1e803bfd9",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libxml2",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libxml2-2.12.7-hc051c1a_1",
        sha256 = "576ea9134176636283ff052897bf7a91ffd8ac35b2c505dfde2890ec52849698",
        archive_type = "conda",
        exclude = ["share/man", "share/doc"],
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libzlib",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libzlib-1.2.13-h4ab18f5_6",
        sha256 = "8ced4afed6322172182af503f21725d072a589a6eb918f8a58135c1e00d35980",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_llvm_openmp",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "llvm-openmp-19.1.6-h024ca30_0",
        sha256 = "9e385c2a8169d951cf153221fb7fbb3dc8f1e5ac77371edee7329f8721dbe1ae",
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
        dist_name = "lz4-4.3.3-py310hb259640_1",
        sha256 = "432b919e052fd0e0687ae7395f147629ba708457e17e491db1d7a286a83d6f41",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_lz4_c",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "lz4-c-1.9.4-hcb278e6_0",
        sha256 = "1b4c105a887f9b2041219d57036f72c4739ab9e9fe5a1486f094e58c76b31f5f",
        exclude = ["share/man"],
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_markupsafe",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "markupsafe-3.0.2-py310h89163eb_1",
        sha256 = "0bed20ec27dcbcaf04f02b2345358e1161fb338f8423a4ada1cf0f4d46918741",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_mkl",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "mkl-2023.2.0-h84fe81f_50496",
        sha256 = "046073737bf73153b0c39e343b197cdf0b7867d336962369407465a17ea5979a",
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
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_narwhals",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "narwhals-1.20.1-pyhd8ed1ab_0",
        sha256 = "0f3262d8d1d336dfa1a3b5df42e82aaf616be2790724719eff8cf8fc31aa3df6",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_ncurses",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "ncurses-6.5-he02047a_1",
        sha256 = "6a1d5d8634c1a07913f1c525db6455918cbc589d745fac46d9d6e30340c8731a",
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
        dist_name = "numexpr-2.10.1-mkl_py310h78a85ae_2",
        sha256 = "a251988d4d1a589c0d042ea2db68f905fa83e9ac643213a4bcf935f31a68f6f7",
        exclude = ["lib/python*/site-packages/numexpr/tests"],
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_numpy",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "numpy-1.26.4-py310hb13e2d6_0",
        sha256 = "028fe2ea8e915a0a032b75165f11747770326f3d767e642880540c60a3256425",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_openjpeg",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "openjpeg-2.5.0-h7d73246_1",
        sha256 = "a715cba5649f12a1dca53dfd72fc49577152041f033d7595cf4b6a655a5b93b6",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_openssl",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "openssl-1.1.1w-hd590300_0",
        sha256 = "4fe19885c77f0758084feb54954bd1977dfeeab7134fba0a1d9c0cfff821d6bd",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_packaging",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "packaging-24.2-pyhd8ed1ab_2",
        sha256 = "da157b19bcd398b9804c5c52fc000fcb8ab0525bdb9c70f95beaa0bb42f85af1",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_pandas",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "pandas-2.1.4-py310hcc13569_0",
        sha256 = "d0743541397140a25a89ab0686933005a4c104d95c23ff1c322f903a50b18099",
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
        dist_name = "pillow-9.2.0-py310h454ad03_3",
        sha256 = "202cc5b4c60e32096b67791f822699bf91670584ac3db7e86ebb1b6a4c584218",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_pkgutil_resolve_name",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "pkgutil-resolve-name-1.3.10-pyhd8ed1ab_2",
        sha256 = "adb2dde5b4f7da70ae81309cce6188ed3286ff280355cf1931b45d91164d2ad8",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_plotly",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "plotly-5.24.1-pyhd8ed1ab_1",
        sha256 = "d1bbf2d80105bfc8a7ed9817888f4a1686ed393d6435572921add09cc9347c1c",
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
        ],
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_polars_lts_cpu",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "polars-lts-cpu-1.17.1-py310hf9ba7b5_0",
        sha256 = "5a6700df00cce6d2d610e771ca1b3522532da7288044d8f26d52a9b8435d6eff",
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
        dist_name = "pycparser-2.22-pyhd8ed1ab_0",
        sha256 = "406001ebf017688b1a1554b49127ca3a4ac4626ec0fd51dc75ffa4415b720b64",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_pyfaidx",
        base_urls = ["https://conda.anaconda.org/bioconda/noarch"],
        dist_name = "pyfaidx-0.8.1.3-pyhdfd78af_0",
        sha256 = "00e960f12330c3bc44fccac0417f2602143d145a5ad0e6ac7599b002ec7f1cf4",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_pysam",
        base_urls = ["https://conda.anaconda.org/bioconda/linux-64"],
        dist_name = "pysam-0.21.0-py310hff46b53_0",
        sha256 = "c231557375a2fb4db78ecdd0d63a660c9316244885e9cc8e5eee75d0f5cb97b3",
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
        dist_name = "pytables-3.8.0-py310hde6a235_1",
        sha256 = "2b444e6e7fe42ac9a025661a4bb00b1de4d8e08c260d9e242cb790cb8cd32ba7",
        patch_args = ["-p1"],
        patches = ["@cellranger//third-party/conda:pytables-3.8.0-slurmdir.patch"],
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
        dist_name = "python-3.10.8-h257c98d_0_cpython",
        sha256 = "090a5d0ed7acf75479664f4751f69833c9aa183f73366fe83681f4e85f72670b",
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
        dist_name = "python-dateutil-2.9.0.post0-pyhff2d567_1",
        sha256 = "a50052536f1ef8516ed11a844f9413661829aa083304dc624c5925298d078d79",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_python_tzdata",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "python-tzdata-2024.2-pyhd8ed1ab_1",
        sha256 = "57c9a02ec25926fb48edca59b9ede107823e5d5c473b94a0e05cc0b9a193a642",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_python_abi",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "python_abi-3.10-5_cp310",
        sha256 = "074d2f0b31f0333b7e553042b17ea54714b74263f8adda9a68a4bd8c7e219971",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_pytz",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "pytz-2024.2-pyhd8ed1ab_1",
        sha256 = "0a7c706b2eb13f7da5692d9ddf1567209964875710b471de6f2743b33d1ba960",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_pyvcf3",
        base_urls = ["https://conda.anaconda.org/bioconda/noarch"],
        dist_name = "pyvcf3-1.0.3-pyhdfd78af_0",
        sha256 = "5283bca618cea395f5aabdc40a9a136c0203c247e16935a2b250e84670e1dccb",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_referencing",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "referencing-0.35.1-pyhd8ed1ab_1",
        sha256 = "f972eecb4dc8e06257af37642f92b0f2df04a7fe4c950f2e1045505e5e93985f",
        exclude = ["site-packages/referencing/tests"],
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_requests",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "requests-2.32.3-pyhd8ed1ab_1",
        sha256 = "d701ca1136197aa121bbbe0e8c18db6b5c94acbd041c2b43c70e5ae104e1d8ad",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_rpds_py",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "rpds-py-0.22.3-py310h505e2c1_0",
        sha256 = "e13019600e75707126118cf3f02187e7dd96f475a82e8fa06e59091f76159274",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_scikit_learn",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "scikit-learn-1.5.2-py310h27f47ee_1",
        sha256 = "777580d5ba89c5382fa63807a7981ae2261784258e84f5a9e747f5bd3d3428f3",
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
        dist_name = "scipy-1.10.1-py310ha4c1d20_3",
        sha256 = "c7beb091db82a1be2fa9dafb878695b1e8bd6d7efe7764afa457cabfea2a93d3",
        exclude = [
            "lib/python*/site-packages/scipy/linalg/src/id_dist/doc",
            "lib/python*/site-packages/scipy/*/*/*/*/tests",
            "lib/python*/site-packages/scipy/*/*/*/tests",
            "lib/python*/site-packages/scipy/*/*/tests",
            "lib/python*/site-packages/scipy/*/tests",
        ],
        archive_type = "conda",
        exclude_deps = ["fftw", "pooch"],
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_setuptools",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "setuptools-75.6.0-pyhff2d567_1",
        sha256 = "abb12e1dd515b13660aacb5d0fd43835bc2186cab472df25b7716cd65e095111",
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
        dist_name = "six-1.17.0-pyhd8ed1ab_0",
        sha256 = "41db0180680cc67c3fa76544ffd48d6a5679d96f4b71d7498a759e94edc9a2db",
        conda_repo = name,
        archive_type = "conda",
    )
    conda_package_repository(
        name = "conda_package_snappy",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "snappy-1.1.10-hdb0a2a9_1",
        sha256 = "082eadbc355016e948f1acc2f16e721ae362ecdaa204cbd60136ada19bd43f3a",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_statsmodels",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "statsmodels-0.14.4-py310hf462985_0",
        sha256 = "a060f9b7e9bff3cae075a00e278089893c20cc0663ced09f9c4d92522ce76a21",
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
        dist_name = "tbb-2021.13.0-hceb3a55_1",
        sha256 = "65463732129899770d54b1fbf30e1bb82fdebda9d7553caf08d23db4590cd691",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_tenacity",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "tenacity-9.0.0-pyhd8ed1ab_1",
        sha256 = "dcf2155fb959773fb102066bfab8e7d79aff67054d142716979274a43fc85735",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_threadpoolctl",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "threadpoolctl-3.5.0-pyhc1e730c_0",
        sha256 = "45e402941f6bed094022c5726a2ca494e6224b85180d2367fb6ddd9aea68079d",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_tk",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "tk-8.6.13-noxft_h4845f30_101",
        sha256 = "e0569c9caa68bf476bead1bed3d79650bb080b532c64a4af7d8ca286c08dea4e",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_typing_extensions",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "typing_extensions-4.12.2-pyha770c72_1",
        sha256 = "337be7af5af8b2817f115b3b68870208b30c31d3439bec07bfb2d8f4823e3568",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_tzdata",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "tzdata-2024b-hc8b5060_0",
        sha256 = "4fde5c3008bf5d2db82f2b50204464314cc3c91c1d953652f7bd01d9e52aefdf",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_urllib3",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "urllib3-2.3.0-pyhd8ed1ab_0",
        sha256 = "114919ffa80c328127dab9c8e7a38f9d563c617691fb81fccb11c1e86763727e",
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
        name = "conda_package_xz",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "xz-5.2.6-h166bdaf_0",
        sha256 = "03a6d28ded42af8a347345f82f3eebdd6807a08526d47899a42d62d319609162",
        licenses = ["@rules_license//licenses/spdx:LGPL-2.1"],
        exclude = [
            "bin/*cmp",
            "bin/*diff",
            "bin/*grep",
            "bin/*less",
            "bin/*more",
            "share",
        ],
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_zipp",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "zipp-3.21.0-pyhd8ed1ab_1",
        sha256 = "567c04f124525c97a096b65769834b7acb047db24b15a56888a322bf3966c3e1",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_zlib",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "zlib-1.2.13-h4ab18f5_6",
        sha256 = "534824ea44939f3e59ca8ebb95e3ece6f50f9d2a0e69999fbc692311252ed6ac",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_zlib_ng",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "zlib-ng-2.0.7-h0b41bf4_0",
        sha256 = "6b3a22b7cc219e8d83f16c1ceba67aa51e0b7e3bcc4a647b97a0a510559b0477",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_zstandard",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "zstandard-0.23.0-py310ha39cb0e_1",
        sha256 = "fcd784735205d6c5f19dcb339f92d2eede9bc42a01ec2c384381ee1b6089d4f6",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_zstd",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "zstd-1.5.6-ha6fb4c9_0",
        sha256 = "c558b9cc01d9c1444031bd1ce4b9cff86f9085765f17627a6cd85fc623c8a02b",
        exclude = ["share/man"],
        licenses = ["@rules_license//licenses/spdx:BSD-3-Clause"],
        archive_type = "conda",
        conda_repo = name,
    )
    anaconda_repository(
        name = name,
        conda_packages = [
            "_libgcc_mutex",
            "altair",
            "attrs",
            "biopython",
            "blosc",
            "bzip2",
            "cached_property",
            "certifi",
            "cffi",
            "cython",
            "docopt",
            "freetype",
            "h2",
            "h5py",
            "hdf5",
            "hpack",
            "htslib",
            "hyperframe",
            "icu",
            "idna",
            "importlib_resources",
            "jinja2",
            "joblib",
            "jpeg",
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
            "libffi",
            "libgcc",
            "libgfortran5",
            "libhwloc",
            "libiconv",
            "liblapack",
            "libnghttp2",
            "libnsl",
            "libpng",
            "libsqlite",
            "libssh2",
            "libstdcxx",
            "libtiff",
            "libuuid",
            "libxcb",
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
            "tenacity",
            "threadpoolctl",
            "tk",
            "typing_extensions",
            "tzdata",
            "urllib3",
            "xz",
            "zipp",
            "zlib",
            "zstandard",
            "zstd",
            "pyarrow",
        ],
        executable_packages = [
            "coverage",
            "cython",
            "ipython",
            "jupyterlab",
            "notebook",
            "python",
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
            "PIL": "pillow",
            "pyBigWig": "pybigwig",
            "polars": "polars-lts-cpu",
            "primer3": "primer3-py",
            "progressbar": "progressbar2",
            "prompt_toolkit": "prompt-toolkit",
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
