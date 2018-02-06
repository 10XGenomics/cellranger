#
# Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
#
# Environment setup for Cell Ranger
# Source this file before running.
#

# Determine path to this script; resolve symlinks
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do
    DIR="$( cd -P "$( dirname "$SOURCE" )" > /dev/null && pwd )"
    SOURCE="$(readlink "$SOURCE")"
    [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE"
done
DIR="$( cd -P "$( dirname "$SOURCE" )" > /dev/null && pwd )"

export PATH="$DIR/bin:$DIR/lib/bin:$PATH"
export MROPATH="$DIR/mro:$MROPATH"
export PYTHONPATH="$DIR/lib/python:$PYTHONPATH"
