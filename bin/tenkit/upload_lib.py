#!/usr/bin/env python3
#
# Copyright (c) 2016 10x Genomics, Inc. All rights reserved.
"""Diagnostic tar file uploader."""

from __future__ import annotations

import json
import math
import os
import ssl
import sys
import time
import xml.etree.ElementTree as ET
from http.client import HTTPSConnection
from typing import AnyStr, BinaryIO
from urllib.parse import urlencode, urlparse

import certifi

SUPPORT_URL = "support.10xgenomics.com"
GET_UPLOAD_INFO_FAIL_MSG = f"""Could not contact https://{SUPPORT_URL} (%s)
   Please make sure you have Internet connectivity and try again,
   or contact support@10xgenomics.com for help."""


def _parse_args(argv: list[str]) -> tuple[str, str]:
    """Parse command line."""
    usage_message = """Usage:
        {} upload <your_email> <file>""".format(
        os.getenv("TENX_PRODUCT", "")
    )

    if len(argv) != 3:
        print(usage_message)
        sys.exit(1)
    email, filename = argv[1:]
    return email, filename


def _verify_file(filename: str) -> tuple[str, int]:
    """Verify file input."""
    try:
        if not os.path.isfile(filename):
            print("Error: Please specify a file.")
            sys.exit(1)
        fname = os.path.basename(filename)
        size = os.path.getsize(filename)
    except Exception as e:  # pylint: disable=broad-except
        print("Error: " + str(e))
        sys.exit(1)
    return fname, size


def get_cert_paths() -> ssl.DefaultVerifyPaths:
    """Return the location of the TLS certificates.

    Search order: environment variables, host system, then bundled certificates

    Environment variables:
      SSL_CERT_FILE: path to a single file of certificates
      SSL_CERT_DIR: path to a directory of certificates
    """
    ssl_paths = ssl.get_default_verify_paths()
    cafile = next(
        (
            x
            for x in (
                ssl_paths.cafile,
                "/etc/pki/tls/cert.pem",  # Amazon Linux, CentOS, Fedora
                "/etc/ssl/ca-bundle.pem",  # OpenSUSE
                "/etc/ssl/certs/ca-certificates.crt",  # Arch Linux, Debian, Ubuntu
                certifi.where(),
            )
            if x and os.access(x, os.R_OK)
        ),
        None,
    )
    capath = next(
        (
            x
            for x in (
                ssl_paths.capath,
                "/etc/pki/tls/certs",  # Amazon Linux, CentOS, Fedora
                "/etc/ssl/certs",  # Arch Linux, Debian, OpenSUSE, Ubuntu
            )
            if x and os.access(x, os.R_OK)
        ),
        None,
    )
    return ssl.DefaultVerifyPaths(
        cafile=cafile,
        capath=capath,
        openssl_cafile_env=ssl_paths.openssl_cafile_env,
        openssl_cafile=ssl_paths.openssl_cafile,
        openssl_capath_env=ssl_paths.openssl_capath_env,
        openssl_capath=ssl_paths.openssl_capath,
    )


def _get_upload_info(email: str, fname: str, size: int, context: ssl.SSLContext) -> str:
    """Retrieve upload information."""
    try:
        print("Getting upload information...")
        sys.stdout.flush()

        conn = HTTPSConnection(SUPPORT_URL, context=context)
        conn.request(
            "GET",
            "/uploadto10x.json?" + urlencode({"email": email, "fname": fname, "size": size}),
        )
        resp = conn.getresponse()
        if resp.status != 200:
            print(GET_UPLOAD_INFO_FAIL_MSG % (str(resp.status) + " " + resp.reason))
            sys.exit(1)

        # Relay any permanent errors or deprecations to the user.
        upload = json.loads(resp.read())
        if "error" in upload:
            print("\nError: " + upload["error"])
            sys.exit(1)

        # Get upload url
        upurl = upload["url"]

        conn.close()
    except Exception as e:  # pylint: disable=broad-except
        print("\nError: " + GET_UPLOAD_INFO_FAIL_MSG % e)
        sys.exit(1)
    return upurl


def pn(n: int) -> str:
    """Pretty print number."""
    return format(n, ",d")


class FileReader(BinaryIO):
    def __init__(self, f: BinaryIO, bsize: int):
        self.f = f
        # Track progress
        self.bcount = 0
        self.bsize = bsize
        self.start = time.monotonic()

    def read(self, csize: int) -> bytes:
        """Calculate progress."""
        self.bcount = min(self.bcount + csize, self.bsize)
        elap = float(time.monotonic() - self.start)
        byps = float(self.bcount) / elap
        mbps = float(self.bcount * 8) / elap / 1000000.0
        pctg = int(float(self.bcount) / float(self.bsize) * 100.0)
        etas = math.ceil((self.bsize - self.bcount) / byps)
        etam = int(etas / 60)
        etas = int(etas % 60)

        # Update progress display
        sys.stdout.write("\r")
        sys.stdout.write(
            "{} [{}] {} {}  eta {} {}  ".format(
                ("%d%%" % pctg).rjust(4),  # percentage
                (("=" * (pctg // 3)) + ">").ljust(34),  # progress bar
                pn(self.bcount).rjust(len(pn(self.bsize))),  # bytes sent
                (f"{mbps:.2f}Mb/s").rjust(10),  # bitrate
                ("%dm" % etam),  # eta min
                ("%ds" % etas).rjust(3),  # eta sec
            )
        )
        sys.stdout.flush()

        # Return file chunk
        return self.f.read(csize)

    def __enter__(self) -> BinaryIO:
        raise NotImplementedError

    def write(self, __s: AnyStr) -> int:
        raise NotImplementedError


def _upload(filename: str, size: int, upurl: str, context: ssl.SSLContext):
    """Upload file to destination."""
    try:
        print("Uploading file...")
        print("Size: " + pn(size) + " bytes")
        sys.stdout.flush()

        url = urlparse(upurl)
        conn = HTTPSConnection(url.netloc, context=context)
        with open(filename, "rb") as fp:
            conn.request(
                "PUT",
                f"{url.path}?{url.query}",
                body=FileReader(fp, size),
                headers={"Content-Length": str(size)},
            )
            resp = conn.getresponse()
            contents = resp.read()

        # Return any errors from upload destination
        if len(contents) > 0:
            print("\nError:")
            root = ET.fromstring(contents)
            for e in root:
                if e.tag not in ("RequestId", "HostId"):
                    print(f"   {e.tag}: {e.text}")
        else:
            print("\nUpload complete!")
    except Exception as e:  # pylint: disable=broad-except
        print("\nError: " + str(e))
        sys.exit(1)


def main(argv: list[str]):
    email, filename = _parse_args(argv)
    fname, size = _verify_file(filename)

    context = ssl.create_default_context()
    cert_paths = get_cert_paths()
    context.load_verify_locations(cafile=cert_paths.cafile, capath=cert_paths.capath)

    upurl = _get_upload_info(email, fname, size, context)
    _upload(filename, size, upurl, context)


if __name__ == "__main__":
    main(sys.argv)
