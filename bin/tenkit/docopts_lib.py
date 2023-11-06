#!/usr/bin/env python3
#
# Copyright (c) 2016 10x Genomics, Inc. All rights reserved.
#

"""Shell interface for docopt, the CLI description language.

Usage:
  docopts [options] -h <msg> : [<argv>...]

Options:
  -h <msg>, --help=<msg>        The help message in docopt format.
                                If - is given, read the help message from
                                standard input.
                                If no argument is given, print docopts's own
                                help message and quit.
  -V <msg>, --version=<msg>     A version message.
                                If - is given, read the version message from
                                standard input.  If the help message is also
                                read from standard input, it is read first.
                                If no argument is given, print docopts's own
                                version message and quit.
  -O, --options-first           Disallow interspersing options and positional
                                arguments: all arguments starting from the
                                first one that does not begin with a dash will
                                be treated as positional arguments.
  -H, --no-help                 Don't handle --help and --version specially.
  -A <name>                     Export the arguments as a Bash 4.x associative
                                array called <name>.
  -s <str>, --separator=<str>   The string to use to separate the help message
                                from the version message when both are given
                                via standard input. [default: ----]
"""

from __future__ import annotations

__version__ = """docopts 0.6.1+fix
Copyright (C) 2013 Vladimir Keleshev, Lari Rasku.
License MIT <http://opensource.org/licenses/MIT>.
This is free software: you are free to change and redistribute it.
There is NO WARRANTY, to the extent permitted by law.

"""

import re
import sys
from io import StringIO

from docopt import DocoptExit, DocoptLanguageError, docopt


# helper functions
def shellquote(s):
    return "'" + s.replace("'", r"'\''") + "'"


def isbashidentifier(s):
    return re.match(r"^([A-Za-z]|[A-Za-z_][0-9A-Za-z_]+)$", s)


def to_bash(obj):
    if obj is None:
        return ""
    elif isinstance(obj, str):
        return shellquote(obj)
    elif isinstance(obj, bool):
        return "true" if obj else "false"
    elif isinstance(obj, int | float):
        return str(obj)
    elif isinstance(obj, list):
        return "(" + " ".join(shellquote(x) for x in obj) + ")"
    raise ValueError(obj)


def name_mangle(elem):
    if elem in ("-", "--"):
        return None
    elif re.match(r"^<.*>$", elem):
        var = elem[1:-1]
    elif re.match(r"^-[^-]$", elem):
        var = elem[1]
    elif re.match(r"^--.+$", elem):
        var = elem[2:]
    else:
        var = elem
    var = var.replace("-", "_")
    if not isbashidentifier(var):
        raise ValueError(elem)
    return var


def main():
    # parse docopts's own arguments
    try:
        args = docopt(__doc__, help=False, options_first=True)
    except DocoptExit as ex:
        message: str = ex.args[0]
        if message.startswith(("-h", "--help")):
            print(__doc__.strip())
            sys.exit()
        if message.startswith(("-V", "--version")):
            print(__version__.strip())
            sys.exit()
        else:
            sys.exit(message)

    argv = args["<argv>"]
    doc = args["--help"]
    version = args["--version"]
    options_first = args["--options-first"]
    help_text = not args["--no-help"]
    name = args["-A"]
    separator = args["--separator"]

    if doc == "-" and version == "-":
        doc, version = (page.strip() for page in sys.stdin.read().split(separator, 1))
    elif doc == "-":
        doc = sys.stdin.read().strip()
    elif version == "-":
        version = sys.stdin.read().strip()

    # parse options or abort if there is an error in docopt
    # temporarily redirect stdout to a StringIO so we can catch docopt()
    # output on --help and --version
    stdout = sys.stdout
    exit_message = None
    try:
        sys.stdout = StringIO()
        args = docopt(doc, argv, help_text, version, options_first)
    except DocoptLanguageError as ex:
        # invalid docstring by user
        sys.exit(f"{sys.argv[0]}: invalid doc argument: {ex}")
    except DocoptExit as ex:
        # invoked with invalid arguments
        exit_message = f"echo {shellquote(str(ex))} >&2\nexit 64"
    except SystemExit:
        # --help or --version found and --no-help was not given
        exit_message = f"echo -n {shellquote(sys.stdout.getvalue())}\nexit 0"
    finally:
        # restore stdout to normal and quit if a docopt parse error happened
        sys.stdout.close()
        sys.stdout = stdout
        if exit_message:
            print(exit_message)
            sys.exit()

    if name is not None:
        if not isbashidentifier(name):
            sys.exit(f"{sys.argv[0]}: not a valid Bash identifier: {name}")
        # fake nested Bash arrays for repeatable arguments with values
        arrays = {elem: value for elem, value in args.items() if isinstance(value, list)}
        for elem, value in arrays.items():
            del args[elem]
            args[elem + ",#"] = len(value)
            args.update(("%s,%d" % (elem, i), v) for i, v in enumerate(value))
        print(f"declare -A {name}")
        for elem, value in args.items():
            print(f"{name}[{shellquote(elem)}]={to_bash(value)}")
    else:
        try:
            variables = dict(
                zip((name_mangle(n) for n in args.keys()), (to_bash(v) for v in args.values()))
            )
        except ValueError as ex:
            sys.exit(f"{sys.argv[0]}: name could not be mangled into a valid Bash identifier: {ex}")
        else:
            variables.pop(None, None)
            args.pop("-", None)
            args.pop("--", None)
        if len(variables) < len(args):
            sys.exit(f"{sys.argv[0]}: two or more elements have identically mangled names")
        for var, value in variables.items():
            print(f"{var}={value}")


if __name__ == "__main__":
    main()
