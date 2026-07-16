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
        dist_name = "altair-6.1.0-pyhd8ed1ab_0",
        sha256 = "fdb791d8310f9969d8d0254bce45620f0d1c3d7a3101c1ef5efb2ffdb1471f8c",
        archive_type = "conda",
        extra_deps = ["pyarrow"],
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_attrs",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "attrs-26.1.0-pyhcf101f3_0",
        sha256 = "1b6124230bb4e571b1b9401537ecff575b7b109cc3a21ee019f65e083b8399ab",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_aws_c_auth",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "aws-c-auth-0.10.1-ha62d5e7_3",
        sha256 = "ccbf2cc4bea4aab6e071d67ecc2743197759f6df855787e7a5f57f7973f913a2",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_aws_c_cal",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "aws-c-cal-0.9.13-h2c9d079_1",
        sha256 = "f21d648349a318f4ae457ea5403d542ba6c0e0343b8642038523dd612b2a5064",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_aws_c_common",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "aws-c-common-0.12.6-hb03c661_0",
        sha256 = "926a5b9de0a586e88669d81de717c8dd3218c51ce55658e8a16af7e7fe87c833",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_aws_c_compression",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "aws-c-compression-0.3.2-h8b1a151_0",
        sha256 = "1838bdc077b77168416801f4715335b65e9223f83641a2c28644f8acd8f9db0e",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_aws_c_http",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "aws-c-http-0.10.13-h4bacb7b_0",
        sha256 = "38cfc8894db6729770ac18f900296c3f7c20f349a5586a8d8e1a62571fce61d5",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_aws_c_io",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "aws-c-io-0.26.3-hb18f61d_2",
        sha256 = "eee7f7aa2c5b9e0a31edba7b81482036fbe751c40bc6697fd057fbd2c656406b",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_aws_c_s3",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "aws-c-s3-0.12.2-he6ee468_1",
        sha256 = "4cecb4d595b7cf558087c37b8131cae5204b2c64d75f6b951dc3731d3f872bb8",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_aws_c_sdkutils",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "aws-c-sdkutils-0.2.4-h8b1a151_4",
        sha256 = "9d62c5029f6f8219368a8665f0a549da572dc777f52413b7d75609cacdbc02cc",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_aws_checksums",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "aws-checksums-0.2.10-h8b1a151_0",
        sha256 = "09472dd5fa4473cffd44741ee4c1112f2c76d7168d1343de53c2ad283dc1efa6",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_backports_zstd",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "backports.zstd-1.5.0-py312h90b7ffd_0",
        sha256 = "a2b08a4e5e549b5f67c38edffd175437e2208547a7e67b5fa5373b67ef419e50",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_biopython",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "biopython-1.87-py312h4c3975b_0",
        sha256 = "8439485bab2f8ae5af4069b624d76186939c5218d3fdf33c49fc67a48605ccc5",
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
        dist_name = "brotli-python-1.2.0-py312hdb49522_1",
        sha256 = "49df13a1bb5e388ca0e4e87022260f9501ed4192656d23dc9d9a1b4bf3787918",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_bzip2",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "bzip2-1.0.8-hda65f42_9",
        sha256 = "0b75d45f0bba3e95dc693336fa51f40ea28c980131fec438afb7ce6118ed05f6",
        exclude = ["man"],
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_c_ares",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "c-ares-1.34.6-hb03c661_0",
        sha256 = "cc9accf72fa028d31c2a038460787751127317dcfa991f8d1f1babf216bb454e",
        archive_type = "conda",
        exclude = ["share/man"],
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_c_blosc2",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "c-blosc2-3.1.2-hc31b594_0",
        sha256 = "ee9b4248b9cad8ded8c5bf03a2624bf8b8670f0779d6cca775c06be7a54dea0c",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_ca_certificates",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "ca-certificates-2026.5.20-hbd8a1cb_0",
        sha256 = "9812a303a1395e1dafbd92e5bc8a1ff6013bcbba0a09c7f03a8d23e43560aa9b",
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
        dist_name = "certifi-2026.5.20-pyhd8ed1ab_0",
        sha256 = "645655a3510e38e625da136595f3f16f2130c3263630cc3bc8f60f619ddbe490",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_charset_normalizer",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "charset-normalizer-3.4.7-pyhd8ed1ab_0",
        sha256 = "3f9483d62ce24ecd063f8a5a714448445dc8d9e201147c46699fc0033e824457",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_cpython",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "cpython-3.12.13-py312hd8ed1ab_0",
        sha256 = "d3e9bbd7340199527f28bbacf947702368f31de60c433a16446767d3c6aaf6fe",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_cython",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "cython-3.2.5-py312h68e6be4_0",
        sha256 = "d89f74179b0758dc2c314b0a3a04a7738ed1f3600529060d78eb26987c22d885",
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
        name = "conda_package_eigen",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "eigen-5.0.1-hc65338a_0",
        sha256 = "f122c5bb618532eb40124f34dc3d467b9142c4a573c206e3e6a51df671345d6a",
        archive_type = "conda",
        conda_repo = name,
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
        dist_name = "h5py-3.16.0-nompi_py312ha829cd9_102",
        sha256 = "578ab0db104435ac003061e103789299720fc49b8b3e8401dabd5130a1ef8663",
        exclude = ["lib/python3.*/site-packages/h5py/tests"],
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_hdf5",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "hdf5-2.1.0-nompi_h87a9417_105",
        sha256 = "beb8a2fb18924ca7b5b82cfb50f008f882f577daef2c00ed88022abea35fec76",
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
        dist_name = "htslib-1.23.1-h633afcb_0",
        sha256 = "d0977efb1885c9f00ca76ee13e7c0f9a8936bf271eec25ffa78606cd816a13c9",
        exclude = ["share/man"],
        conda_repo = name,
        archive_type = "conda",
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
        dist_name = "icu-78.3-h33c6efd_0",
        sha256 = "fbf86c4a59c2ed05bbffb2ba25c7ed94f6185ec30ecb691615d42342baa1a16a",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_idna",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "idna-3.17-pyhcf101f3_0",
        sha256 = "f9fe1f9e539c544405ccb7ba632d4ba79edf243c05554d76ace073158a80b691",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_imageio",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "imageio-2.37.0-pyhfb79c49_0",
        sha256 = "8ef69fa00c68fad34a3b7b260ea774fda9bd9274fd706d3baffb9519fd0063fe",
        exclude = ["site-packages/imageio/resources/images"],
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_importlib_metadata",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "importlib-metadata-9.0.0-pyhcf101f3_0",
        sha256 = "43e2a5497cad1598ff88a3e69f69bc88b7b8f141fa63c60eab5db296317318b8",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_jinja2",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "jinja2-3.1.6-pyhcf101f3_1",
        sha256 = "fc9ca7348a4f25fed2079f2153ecdcf5f9cf2a0bc36c4172420ca09e1849df7b",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_joblib",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "joblib-1.5.3-pyhd8ed1ab_0",
        sha256 = "301539229d7be6420c084490b8145583291123f0ce6b92f56be5948a2c83a379",
        exclude = ["site-packages/joblib/test"],
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_jsonschema",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "jsonschema-4.26.0-pyhcf101f3_0",
        sha256 = "db973a37d75db8e19b5f44bbbdaead0c68dde745407f281e2a7fe4db74ec51d7",
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
        dist_name = "krb5-1.22.2-ha1258a1_0",
        sha256 = "3e307628ca3527448dd1cb14ad7bb9d04d1d28c7d4c5f97ba196ae984571dd25",
        exclude = ["share/man", "share/examples"],
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_lazy_loader",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "lazy-loader-0.5-pyhd8ed1ab_0",
        sha256 = "1a88069ac61d2756ccaf26a6c206ab4d56610fb054bd2fffb5df4cd0744ab78e",
        archive_type = "conda",
        exclude = ["site-packages/lazy_loader/tests/fake_pkg"],
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_lcms2",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "lcms2-2.19.1-h0c24ade_1",
        sha256 = "112b5b9462572d970f4abd2912f76a25ee7db158b1e7260163d91dd8a630db84",
        exclude = ["share/man"],
        conda_repo = name,
        archive_type = "conda",
    )
    conda_package_repository(
        name = "conda_package_lerc",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "lerc-4.1.0-hdb68285_0",
        sha256 = "f84cb54782f7e9cea95e810ea8fef186e0652d0fa73d3009914fa2c1262594e1",
        exclude = ["share/man"],
        conda_repo = name,
        archive_type = "conda",
    )
    conda_package_repository(
        name = "conda_package_libaec",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libaec-1.1.5-h088129d_0",
        sha256 = "822e4ae421a7e9c04e841323526321185f6659222325e1a9aedec811c686e688",
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
        dist_name = "libcurl-8.20.0-hcf29cc6_0",
        sha256 = "75963a5dd913311f59a35dbd307592f4fa754c4808aff9c33edb430c415e38eb",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libdeflate",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libdeflate-1.25-h17f619e_0",
        sha256 = "aa8e8c4be9a2e81610ddf574e05b64ee131fab5e0e3693210c9d6d2fba32c680",
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
        dist_name = "libexpat-2.8.1-hecca717_0",
        sha256 = "363018b25fdb5534c79783d912bd4b685a3547f4fc5996357ad548899b0ee8e7",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libffi",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libffi-3.5.2-h3435931_0",
        sha256 = "31f19b6a88ce40ebc0d5a992c131f57d919f73c0b92cd1617a5bec83f6e961e6",
        licenses = ["@rules_license//licenses/spdx:MIT"],
        exclude = ["share/info", "share/man"],
        conda_repo = name,
        archive_type = "conda",
    )
    conda_package_repository(
        name = "conda_package_libfreetype",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libfreetype-2.14.3-ha770c72_0",
        sha256 = "38f014a7129e644636e46064ecd6b1945e729c2140e21d75bb476af39e692db2",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libfreetype6",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libfreetype6-2.14.3-h73754d4_0",
        sha256 = "16f020f96da79db1863fcdd8f2b8f4f7d52f177dd4c58601e38e9182e91adf1d",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libgcc",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libgcc-15.2.0-he0feb66_19",
        sha256 = "8e0a3b5e41272e5678499b5dfc4cddb673f9e935de01eb0767ce857001229f46",
        archive_type = "conda",
        exclude_deps = ["_openmp_mutex"],
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libgcc_ng",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libgcc-ng-15.2.0-h69a702a_19",
        sha256 = "9dcf54adfaa5e861123c2da4f2f0451a685464ea7e5a41ad91cf67b31d658d98",
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
        dist_name = "libgfortran5-15.2.0-h68bc16d_19",
        sha256 = "057978bb69fea29ed715a9b98adf71015c31baecc4aeb2bfc20d4fd5d83579d4",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libhwloc",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libhwloc-2.13.0-default_he001693_1000",
        sha256 = "5041d295813dfb84652557839825880aae296222ab725972285c5abe3b6e4288",
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
        dist_name = "libjpeg-turbo-3.1.4.1-hb03c661_0",
        sha256 = "10056646c28115b174de81a44e23e3a0a3b95b5347d2e6c45cc6d49d35294256",
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
        dist_name = "liblzma-5.8.3-hb03c661_0",
        sha256 = "ec30e52a3c1bf7d0425380a189d209a52baa03f22fb66dd3eb587acaa765bd6d",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_liblzma_devel",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "liblzma-devel-5.8.3-hb03c661_0",
        sha256 = "7858f6a173206bc8a5bdc8e75690483bb66c0dcc3809ac1cb43c561a4723623a",
        archive_type = "conda",
        conda_repo = name,
        exclude = ["share"],
    )
    conda_package_repository(
        name = "conda_package_libnghttp2",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libnghttp2-1.68.1-h877daf1_0",
        sha256 = "663444d77a42f2265f54fb8b48c5450bfff4388d9c0f8253dd7855f0d993153f",
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
        dist_name = "libpng-1.6.58-h421ea60_0",
        sha256 = "377cfe037f3eeb3b1bf3ad333f724a64d32f315ee1958581fc671891d63d3f89",
        exclude = ["share/man"],
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libsqlite",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libsqlite-3.53.1-h0c1763c_0",
        sha256 = "54cdcd3214313b62c2a8ee277e6f42150d9b748264c1b70d958bf735e420ef8d",
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
        dist_name = "libstdcxx-15.2.0-h934c35e_19",
        sha256 = "dff1058c76ec6b8759e41cefa2508162d00e4a5e6721aa68ec3fd10094e702dc",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libstdcxx_ng",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libstdcxx-ng-15.2.0-hdf11a46_19",
        sha256 = "0672b6b6e1791c92e8eccad58081a99d614fcf82bca5841f9dfa3c3e658f83b9",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libtiff",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libtiff-4.7.1-h9d88235_1",
        sha256 = "e5f8c38625aa6d567809733ae04bb71c161a42e44a9fa8227abe61fa5c60ebe0",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libuuid",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libuuid-2.42.1-h5347b49_0",
        sha256 = "3f0edf1280e2f6684a986f821eaa3e123d2694a00b31b96ca0d4a4c12c129231",
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
        dist_name = "libxml2-16-2.15.3-hca6bf5a_0",
        sha256 = "3d44f737c5ae52d5af32682cc1530df433f401f8e58a7533926536244127572a",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libxml2",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libxml2-2.15.3-h49c6c72_0",
        sha256 = "3bc5551720c58591f6ea1146f7d1539c734ed1c40e7b9f5cb8cb7e900c509aba",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_libzlib",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "libzlib-1.3.2-h25fd6f3_2",
        sha256 = "55044c403570f0dc26e6364de4dc5368e5f3fc7ff103e867c487e2b5ab2bcda9",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_llvm_openmp",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "llvm-openmp-22.1.6-h4922eb0_0",
        sha256 = "e74dbafd2b420687cc913f5587050270c8f57042405b6b0d66c3a8013dc104ab",
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
        dist_name = "lz4-4.4.5-py312h3d67a73_1",
        sha256 = "e8ae9141c7afcc95555fca7ff5f91d7a84f094536715211e750569fd4bb2caa4",
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
        dist_name = "markupsafe-3.0.3-py312h8a5da7c_1",
        sha256 = "5f3aad1f3a685ed0b591faad335957dbdb1b73abfd6fc731a0d42718e0653b33",
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
        name = "conda_package_mkl_include",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "mkl-include-2023.2.0-ha770c72_50498",
        sha256 = "667e3250ec4ae0d34430533508f982cb81b9915881540cc07a182aa59e4db23c",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_mkl_service",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "mkl-service-2.4.1-py312hb582d3c_1",
        sha256 = "01ebb935c4f100e54fa82182bdc2d1554c10243ddfcd32d1e44cc8aba589aa6f",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_narwhals",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "narwhals-2.22.0-pyhcf101f3_0",
        sha256 = "ae154742b6ef2a0eecf960d5469b45fdbbf5bffa44a55bcd9cf93796690e461e",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_ncurses",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "ncurses-6.6-hdb14827_0",
        sha256 = "fc89f74bbe362fb29fa3c037697a89bec140b346a2469a90f7936d1d7ea4d8a3",
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
        name = "conda_package_networkx",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "networkx-3.6.1-pyhcf101f3_0",
        sha256 = "f6a82172afc50e54741f6f84527ef10424326611503c64e359e25a19a8e4c1c6",
        exclude = ["share/doc"],
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
        dist_name = "numpy-2.4.6-py312h33ff503_0",
        sha256 = "dfcbeadb3e7ad0da7a55a0525884ca34c19584154e13cc4159396b305d1bd445",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_onnxruntime",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "onnxruntime-1.26.0-py312hb39aaa6_0_cpu",
        sha256 = "14dd09f6435081c1318a8325cf979c97c6d78a47009bf09c851c85f46d6d2066",
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
        dist_name = "openssl-3.6.2-h35e630c_0",
        sha256 = "c0ef482280e38c71a08ad6d71448194b719630345b0c9c60744a2010e8a8e0cb",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_packaging",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "packaging-26.2-pyhc364b38_0",
        sha256 = "3906abfb6511a3bb309e39b9b1b7bc38f50a723971de2395489fd1f379255890",
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
        dist_name = "patsy-1.0.2-pyhcf101f3_0",
        sha256 = "9678f4745e6b82b36fab9657a19665081862268cb079cf9acf878ab2c4fadee9",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_pillow",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "pillow-12.2.0-py312h50c33e8_0",
        sha256 = "fa291f8915114733dc1df9f1627b8c63c517217c1eee1a6ede2ceb5e368cf27a",
        conda_repo = name,
        archive_type = "conda",
    )
    conda_package_repository(
        name = "conda_package_plotly",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "plotly-6.6.0-pyhd8ed1ab_0",
        sha256 = "c418d325359fc7a0074cea7f081ef1bce26e114d2da8a0154c5d27ecc87a08e7",
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
        name = "conda_package_polars",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "polars-1.41.2-pyh58ad624_0",
        sha256 = "7200a9b1c48fe83efa8f5a5fc35d6066a76c28cbd57cbea2f875aa6ead747ae9",
        archive_type = "conda",
        extra_deps = ["pyarrow", "polars-runtime-compat"],
        exclude_deps = ["polars-runtime-32"],
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_polars_runtime_compat",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "polars-runtime-compat-1.41.2-py310hcbd6021_0",
        sha256 = "d0deb81810480aa57500d94647927728c5612e353145fba15314e2e5edbed83d",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_protobuf",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "protobuf-5.28.3-py312h2ec8cdc_0",
        sha256 = "acb2e0ee948e3941f8ed191cb77f654e06538638aed8ccd71cbc78a15242ebbb",
        archive_type = "conda",
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
        name = "conda_package_pyfaidx",
        base_urls = ["https://conda.anaconda.org/bioconda/noarch"],
        dist_name = "pyfaidx-0.9.0.4-pyhdfd78af_0",
        sha256 = "d90a12928b785905f1004579fbc81979458adac45951df32dc1adb9cfe5bff89",
        conda_repo = name,
        archive_type = "conda",
    )
    conda_package_repository(
        name = "conda_package_pysam",
        base_urls = ["https://conda.anaconda.org/bioconda/linux-64"],
        dist_name = "pysam-0.24.0-py312hf5ad864_0",
        sha256 = "2ca1a2f1f41e6ecebdb91810d19a90062a43cc77f63624f2a3ba48c6b3fe6431",
        conda_repo = name,
        archive_type = "conda",
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
        dist_name = "pytables-3.11.1-py312h920095a_3",
        sha256 = "0e886edf85cb5119eb95a6117b47992912ad66f711de4a4dbf305f3abd7c667a",
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
        dist_name = "python-3.12.13-hd63d673_0_cpython",
        sha256 = "a44655c1c3e1d43ed8704890a91e12afd68130414ea2c0872e154e5633a13d7e",
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
        name = "conda_package_python_flatbuffers",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "python-flatbuffers-25.9.23-pyh1e1bc0e_0",
        sha256 = "1933243d9f839bb7bc6c9b75481b6445f50f8b727eb926086e8a2a347fb0467d",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_python_gil",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "python-gil-3.12.13-hd8ed1ab_0",
        sha256 = "97327b9509ae3aae28d27217a5d7bd31aff0ab61a02041e9c6f98c11d8a53b29",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_python_tzdata",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "python-tzdata-2026.2-pyhd8ed1ab_0",
        sha256 = "e943f9c15a6bdba2e1b9f423ab913b3f6b02197b0ef9f8e6b7464d78b59965b9",
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
        dist_name = "referencing-0.37.0-pyhcf101f3_0",
        sha256 = "0577eedfb347ff94d0f2fa6c052c502989b028216996b45c7f21236f25864414",
        exclude = ["site-packages/referencing/tests"],
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_requests",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "requests-2.34.2-pyhcf101f3_0",
        sha256 = "1715246b19c9f85ee022933b4845f2fc14ac9184981b7b7d9b728bec8e9588da",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_rpds_py",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "rpds-py-2026.5.1-py312h192e038_0",
        sha256 = "bc4a5045fd79e68392fb0661c698303c16e88b83d50626c2bc49c403555e900d",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_s2n",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "s2n-1.7.3-hc5a330e_0",
        sha256 = "150a0a5254e8b15ad737549721c7d13406cd96432f3f446e07073dbd98bb2491",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_scikit_image",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "scikit-image-0.26.0-np2py312h4ae17e4_0",
        sha256 = "581a2228e6963b0707562f519ff68d6c97fad44711af56d3dbeb4a7377939cce",
        patch_args = ["-p1"],
        patches = ["@cellranger//third-party/conda:scikit-image-0.26.0-nocache.patch"],
        exclude = [
            "lib/python*/site-packages/skimage/data/__pycache__/_fetchers.cpython-*.pyc",
            "lib/python*/site-packages/skimage/*/tests",
            "lib/python*/site-packages/skimage/*/*/tests",
        ],
        exclude_deps = ["_openmp_mutex"],
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_scikit_learn",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "scikit-learn-1.9.0-np2py312h3226591_0",
        sha256 = "8c0cd0326b5a17ddcf189fc4f119bf6871b7853595c088075847c484a3ed567e",
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
        dist_name = "scipy-1.17.1-py312h54fa4ab_1",
        sha256 = "d5ac05ad45c0d48731eb189c2cbb2bb99f0e3cb7e1acaad373cb2f1f2597fc75",
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
        dist_name = "setuptools-82.0.1-pyh332efcf_0",
        sha256 = "82088a6e4daa33329a30bc26dc19a98c7c1d3f05c0f73ce9845d4eab4924e9e1",
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
        dist_name = "snappy-1.2.2-h03e3b7b_1",
        sha256 = "48f3f6a76c34b2cfe80de9ce7f2283ecb55d5ed47367ba91e8bb8104e12b8f11",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_statsmodels",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "statsmodels-0.14.6-py312h4f23490_0",
        sha256 = "0c61eccf3f71b9812da8ced747b1f22bafd6f66f9a64abe06bbe147a03b7322e",
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
        dist_name = "tbb-2021.13.0-hb700be7_6",
        sha256 = "b3960446ee3d63106807c7fd7681ee9b19bc416aabd199e2a4b96fe89714c3bb",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_tbb_devel",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "tbb-devel-2021.13.0-hab88423_6",
        sha256 = "880e6b2a9a6bf4b36c58335cfd6d4d1d37409529d5101dc325283a70177db884",
        archive_type = "conda",
        exclude = ["share/doc"],
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
        name = "conda_package_tifffile",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "tifffile-2026.6.1-pyhd8ed1ab_0",
        sha256 = "7d534bbfd0b93001296076c6634f2bdb5511a0f418ecbe310c2eed71a6e66f21",
        exclude_deps = ["imagecodecs"],
        conda_repo = name,
        archive_type = "conda",
    )
    conda_package_repository(
        name = "conda_package_tk",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "tk-8.6.13-noxft_h366c992_103",
        sha256 = "cafeec44494f842ffeca27e9c8b0c27ed714f93ac77ddadc6aaf726b5554ebac",
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
        dist_name = "tzdata-2025c-hc9c84f9_1",
        sha256 = "1d30098909076af33a35017eed6f2953af1c769e273a0626a04722ac4acaba3c",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_urllib3",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "urllib3-2.7.0-pyhd8ed1ab_0",
        sha256 = "feff959a816f7988a0893201aa9727bbb7ee1e9cec2c4f0428269b489eb93fb4",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_xorg_libxau",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "xorg-libxau-1.0.12-hb03c661_1",
        sha256 = "6bc6ab7a90a5d8ac94c7e300cc10beb0500eeba4b99822768ca2f2ef356f731b",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_xorg_libxdmcp",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "xorg-libxdmcp-1.1.5-hb03c661_1",
        sha256 = "25d255fb2eef929d21ff660a0c687d38a6d2ccfbcbf0cc6aa738b12af6e9d142",
        conda_repo = name,
        archive_type = "conda",
    )
    conda_package_repository(
        name = "conda_package_xz",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "xz-5.8.3-ha02ee65_0",
        sha256 = "2553fd3ec0a1020b2ca05ca10b0036a596cb0d4bf3645922fcf69dacce0e6679",
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
        archive_type = "conda",
        exclude_deps = ["xz-gpl-tools"],
    )
    conda_package_repository(
        name = "conda_package_xz_tools",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "xz-tools-5.8.3-hb03c661_0",
        sha256 = "162ebd76803464b8c8ebc7d45df32edf0ec717b3bf369a437ae3b0254f22dc2e",
        exclude = ["share"],
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_zipp",
        base_urls = ["https://conda.anaconda.org/conda-forge/noarch"],
        dist_name = "zipp-4.1.0-pyhcf101f3_0",
        sha256 = "210bd31c22bb88f5e2a167df24c95bb5f152b2ada7502f9b8c49d1f5366db423",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_zlib",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "zlib-1.3.2-h25fd6f3_2",
        sha256 = "245c9ee8d688e23661b95e3c6dd7272ca936fabc03d423cdb3cdee1bbcf9f2f2",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_zlib_ng",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "zlib-ng-2.3.3-hceb46e0_1",
        sha256 = "ea4e50c465d70236408cb0bfe0115609fd14db1adcd8bd30d8918e0291f8a75f",
        archive_type = "conda",
        conda_repo = name,
    )
    conda_package_repository(
        name = "conda_package_zstd",
        base_urls = ["https://conda.anaconda.org/conda-forge/linux-64"],
        dist_name = "zstd-1.5.7-hb78ec9c_6",
        sha256 = "68f0206ca6e98fea941e5717cec780ed2873ffabc0e1ed34428c061e2c6268c7",
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
            "cpython",
            "cython",
            "docopt",
            "eigen",
            "h2",
            "h5py",
            "hdf5",
            "hpack",
            "htslib",
            "hyperframe",
            "icu",
            "idna",
            "imageio",
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
            "networkx",
            "numexpr",
            "numpy",
            "onnxruntime",
            "openjpeg",
            "openssl",
            "packaging",
            "pandas",
            "patsy",
            "pillow",
            "plotly",
            "polars",
            "protobuf",
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
            "s2n",
            "scipy",
            "setuptools",
            "six",
            "snappy",
            "statsmodels",
            "tbb",
            "threadpoolctl",
            "tifffile",
            "tk",
            "typing_extensions",
            "tzdata",
            "urllib3",
            "xz",
            "zipp",
            "zlib",
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
