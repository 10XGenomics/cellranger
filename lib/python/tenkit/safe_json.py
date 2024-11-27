#!/usr/bin/env python
#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#

"""Methods for safely encoding values to json.

The json standard does not permit encoding NaN, but Python will still happily
do it, which can cause problems for other programs.  This module contains code
to fix those values up, as well as convert numpy objects to things which encode
properly in json.
"""

from __future__ import annotations

import dataclasses
import json
import math
from collections.abc import Callable, Generator, Iterable, Iterator, Mapping
from pathlib import PurePath
from typing import IO, Any, TypeVar

import numpy as np

_T0 = TypeVar("_T0", str, float, int, bool, None)  # pylint: disable=invalid-name


NAN_STRING: str = "NaN"
POS_INF_STRING: str = "inf"
NEG_INF_STRING: str = "-inf"
_NEG_INF: float = float("-Inf")
_POS_INF: float = float("+Inf")


def desanitize_value(x: _T0) -> float | _T0:
    """Converts back a special numeric value encoded by json_sanitize.

    Args:
        x: a string encoded value or float

    Returns:
        a float or the original value
    """
    if x == NAN_STRING:
        return float(NAN_STRING)
    elif x == POS_INF_STRING:
        return _POS_INF
    elif x == NEG_INF_STRING:
        return _NEG_INF
    else:
        return x


def _sanitize_key(k: Any) -> str:
    if isinstance(k, str):
        return k
    elif isinstance(k, bytes):
        return k.decode()
    elif isinstance(k, bool):
        return "true" if k else "false"
    else:
        return str(k)


def _sanitize_scalar(
    data: bytes | PurePath | float | np.integer | np.floating,
) -> None | str | int | float:
    if isinstance(data, str | int | bool):
        return data
    elif isinstance(data, bytes):
        return data.decode()
    elif isinstance(data, PurePath):
        return str(data)
    elif isinstance(data, np.integer):
        return int(data)
    elif isinstance(data, float | np.floating):
        # Handle floats specially
        if math.isnan(data):
            return NAN_STRING
        elif data == _POS_INF:
            return POS_INF_STRING
        elif data == _NEG_INF:
            return NEG_INF_STRING
        return float(data)
    return data


def _dict_items(
    data,
    # Adapted from the standard library lib/python3.8/json/encoder.py
    ## HACK: hand-optimized bytecode; turn globals into locals
    # pylint: disable=redefined-builtin
    dict: type[dict] = dict,
    fields: Callable[[Any], Iterable[dataclasses.Field[Any]]] = dataclasses.fields,
    getattr=getattr,
    hasattr=hasattr,
    is_dataclass=dataclasses.is_dataclass,
) -> Iterable[tuple[Any, Any]] | None:
    if isinstance(data, dict):
        return data.items()
    elif is_dataclass(data):
        # Don't just use dataclasses.asdict because that does a deep copy.
        return ((field.name, getattr(data, field.name)) for field in fields(data))
    elif hasattr(data, "keys"):
        if hasattr(data, "iteritems"):
            # Not a dict, but close enough.
            return data.iteritems()
        elif hasattr(data, "items"):
            # Not a dict, but close enough.
            return data.items()
        elif hasattr(data, "__getitem__"):
            # Not a dict, but has keys and key lookup
            return ((key, data[key]) for key in data.keys())
    return None


def json_sanitize(data) -> None | dict[str, Any] | int | str | float | list[Any]:
    """Convert data into a form that will correctly serialize to json.

    Yuck yuck yuck yuck yuck!
    The default JSON encoders do bad things if you try to encode
    NaN/Infinity/-Infinity as JSON.  This code takes a nested data structure
    composed of: atoms, dicts, or array-like iteratables and finds all of the
    not-really-a-number floats and converts them to an appropriate string for
    subsequent jsonification.

    Note:
        Don't call this and then pass the result to `json.dump`.  It's
        considerably more efficient to use `dump_numpy`, which avoids
        constructing the entire sanitized object in memory.
    """
    if _is_safe(data):
        # Don't construct a new object if we don't have to.
        return data
    # This really doesn't make me happy. How many cases we we have to test?
    if isinstance(data, bytes | PurePath | float | np.floating | np.integer):
        return _sanitize_scalar(data)
    elif (dict_items := _dict_items(data)) is not None:
        return {_sanitize_key(k): json_sanitize(value) for k, value in dict_items}
    elif hasattr(data, "__iter__"):
        # Anything else that looks like a list. N
        return [json_sanitize(item) for item in data]
    elif hasattr(data, "shape") and data.shape == ():
        # Numpy 0-d array
        return _sanitize_scalar(data.item())
    else:
        return data


def safe_jsonify(data: Any, pretty: bool = False, **kwargs) -> str:
    """Dump an object to a string as json, after sanitizing it.

    Note:
        Rather than writing the result of this call to a file, it's considerably
        more efficient to use `dump_numpy`, which doesn't build the entire
        string up in memory.

    Args:
        data: The object to encode.
        pretty: If true, set various defaults for pretty-printed json.
        **kwargs: Additional arguments to pass to the json encoder.
    """
    kwargs.setdefault("allow_nan", False)
    if pretty:
        _pretty_kwargs(kwargs)
    return json.dumps(data, cls=NumpyAwareJSONEncoder, **kwargs)


def _is_safe(data: Any) -> bool:
    """Returns True if we can use the default encoder.

    This allows us to short-circuit all of the custom logic to use the native
    implementation instead.

    In the case of json_sanitize, this can save quite a bit of ram copying
    objects which don't need to be copied.  In the case of NumpyAwareJSONEncoder,
    it's less important but still means using the native C implementation.

    Args:
        data: The data to check

    Returns:
        True if the object can be safely encoded with the native encoder.
    """
    if data is None or isinstance(data, str | int | bool):
        return True
    elif isinstance(data, float):
        return math.isfinite(data)
    elif isinstance(data, dict):
        return all(_is_safe(key) and _is_safe(value) for key, value in data.items())
    elif isinstance(data, list | tuple):
        return all(_is_safe(value) for value in data)

    return False


def _make_iterencode(
    enc: Callable[[Any, int], Iterable[str]],
    indent: str | None,
    sort_keys: bool,
    item_separator: str,
    key_separator: str,
    # Adapted from the standard library lib/python3.8/json/encoder.py
    ## HACK: hand-optimized bytecode; turn globals into locals
    # pylint: disable=redefined-builtin,too-many-arguments
    float: type[float] = float,
    int: type[int] = int,
    isinstance=isinstance,
    str: type[str] = str,
    bool: type[bool] = bool,
    bytes: type[bytes] = bytes,
    len=len,
    hasattr=hasattr,
    sorted=sorted,
    _is_safe: Callable[[Any], bool] = _is_safe,
    _sanitize_key: Callable[[Any], str] = _sanitize_key,
    _sanitize_scalar: Callable[[Any], None | str | int | float] = _sanitize_scalar,
    _dict_items: Callable[[Any], Iterable[tuple[Any, Any]] | None] = _dict_items,
) -> Callable[[Any, int], Generator[str, None, None]]:
    # pylint: disable=too-many-locals
    float64 = np.float64
    integer = np.integer
    generic = np.generic
    ndarray = np.ndarray

    def _iterencode_iterable(obj: Iterable[Any], depth: int) -> Generator[str, None, None]:
        if hasattr(obj, "__len__") and not len(obj):
            yield "[]"
            return
        depth += 1
        yield "["
        if indent:
            newline_indent = "\n" + indent * depth
            yield newline_indent
        else:
            newline_indent = None
        first = True
        for value in obj:
            if first:
                first = False
            else:
                yield item_separator
                if newline_indent:
                    yield newline_indent
            yield from _iterencode(value, depth)
        if newline_indent is not None:
            yield newline_indent[: -len(indent)]
        yield "]"

    def _iterencode_dict(obj: Iterable[tuple[Any, Any]], depth: int) -> Generator[str, None, None]:
        if hasattr(obj, "__len__") and not len(obj):
            yield "{}"
            return
        depth += 1
        yield "{"
        if indent:
            newline_indent = "\n" + indent * depth
            yield newline_indent
        else:
            newline_indent = None
        first = True
        if sort_keys:
            obj = sorted((_sanitize_key(k), v) for k, v in obj)
        for key, value in obj:
            if first:
                first = False
            else:
                yield item_separator
                if newline_indent:
                    yield newline_indent
            yield from enc(_sanitize_key(key), depth)
            yield key_separator
            yield from _iterencode(value, depth)
        if newline_indent is not None:
            yield newline_indent[: -len(indent)]
        yield "}"

    def _iterencode(obj: Any, depth: int) -> Generator[str, None, None]:
        if obj is None or (_is_safe(obj) if not indent else isinstance(obj, (str, int, bool))):
            yield from enc(obj, depth)
        elif isinstance(obj, bytes | PurePath | float | float64 | integer):
            yield from enc(_sanitize_scalar(obj), depth)
        elif (dict_items := _dict_items(obj)) is not None:
            yield from _iterencode_dict(dict_items, depth)
        elif isinstance(obj, ndarray) and obj.ndim == 0 or isinstance(obj, generic):
            yield from enc(_sanitize_scalar(obj.item()), depth)
        elif hasattr(obj, "__iter__"):
            # Anything else that looks like a list.
            yield from _iterencode_iterable(obj, depth)
        else:
            yield from enc(obj, depth)

    return _iterencode


class NumpyAwareJSONEncoder(json.JSONEncoder):
    """This encoder will convert 1D np.ndarrays to lists.

    For other numpy types, uses obj.item() to extract a python scalar.
    """

    def default(self, o: Any) -> Any:  # pylint: disable=method-hidden
        """Convert a 1D np.ndarray into something json-serializable.

        Returns:
            list, or a single-element array or matrix into its scalar value.
        """
        if _is_safe(o):
            return o
        elif isinstance(o, bytes | PurePath | float | np.floating | np.integer):
            return _sanitize_scalar(o)
        if isinstance(o, np.ndarray):
            if o.ndim >= 1:
                return o.tolist()
            else:
                return _sanitize_scalar(np.asscalar(o))
        elif isinstance(o, np.generic):
            return self.default(o.item())
        elif isinstance(o, Mapping):
            return {_sanitize_key(k): self.default(v) for k, v in o.items()}
        elif hasattr(o, "keys"):
            # Dictionary-like case
            return {_sanitize_key(k): self.default(o[k]) for k in o.keys()}
        elif isinstance(o, Iterable):
            # Anything else that looks like a list.
            return [self.default(item) for item in o]
        return json.JSONEncoder.default(self, o)

    def iterencode(self, o: Any, _one_shot: bool = False) -> Iterator[str]:
        """Encode the given object and yield each string as available.

        For example::

            for chunk in JSONEncoder().iterencode(bigobject):
                mysocket.write(chunk)
        """
        if _is_safe(o):
            return super().iterencode(o, _one_shot)

        return _make_iterencode(
            json.encoder.c_make_encoder(
                {} if self.check_circular else None,
                self.default,
                (
                    json.encoder.c_encode_basestring_ascii
                    if self.ensure_ascii
                    else json.encoder.c_encode_basestring
                ),
                indent := (
                    self.indent
                    if self.indent is None or isinstance(self.indent, str)
                    else " " * self.indent
                ),
                self.key_separator,
                self.item_separator,
                self.sort_keys,
                self.skipkeys,
                self.allow_nan,
            ),
            indent,
            self.sort_keys,
            self.item_separator,
            self.key_separator,
        )(o, 0)


def _pretty_kwargs(kwargs: dict[str, Any]):
    """Update kwargs with defaults for "pretty printing" json."""
    kwargs.setdefault("indent", 4)
    kwargs.setdefault("sort_keys", True)
    kwargs.setdefault("separators", (",", ":"))


def dump_numpy(
    data: Any, fp: IO[str], *, pretty: bool = False, **kwargs  # pylint: disable=invalid-name
) -> None:
    """Dump object to json, converting numpy objects to reasonable JSON.

    Also handles converting bytes to str, dict-like objects such as h5 tables
    into objects, and non-list/tuple iterables into lists.
    """
    kwargs.setdefault("allow_nan", False)
    if pretty:
        _pretty_kwargs(kwargs)
    json.dump(data, fp, cls=NumpyAwareJSONEncoder, **kwargs)
