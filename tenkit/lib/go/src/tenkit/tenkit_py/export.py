#!/usr/bin/env python

"""Generates go source code for any constants defined in the given python
modules."""

import optparse
import os.path
import subprocess
import sys

def main(argv):
    """Process the given python modules (from PYTHONPATH) into go source."""
    parser = optparse.OptionParser(
        usage='Usage: %prog <python_module> [<python_module2>...]')
    parser.add_option('-p', '--package', dest='package',
                      type='string', default=os.getenv('GOPACKAGE', 'main'),
                      help='The go package name.')
    parser.add_option('-f', '--gofmt', dest='gofmt', default=False,
                      action='store_true',
                      help='Pass the output through gofmt')
    parser.add_option('-o', '--output', dest='output', default=None,
                      type='string', action='store',
                      help='Write the output to the given file, instead of stdout.')
    options, argv = parser.parse_args(argv)
    if len(argv) < 2:
        parser.print_usage(sys.stderr)
        return 1
    if options.output:
        with open(options.output, 'w') as writer:
            return write_package(writer, options.package, argv[1:], options.gofmt)
    else:
        return write_package(sys.stdout, options.package, argv[1:], options.gofmt)


_BASE_TYPES = ['string', 'int', 'int64', 'float64', 'bool']


def format_package(writer, package, modules):
    """Pipe write_package through gofmt."""
    gofmt = subprocess.Popen(['gofmt'],
                             stdout=writer,
                             bufsize=-1,
                             stdin=subprocess.PIPE)
    gofmt_writer = gofmt.stdin
    try:
        write_package(gofmt_writer, package, modules)
        gofmt_writer.flush()
    finally:
        gofmt_writer.close()
    return_code = gofmt.wait()
    if return_code != 0:
        sys.stderr.write(
            'gofmt failed with return code %d' % return_code)
    return return_code


def write_package(writer, package, modules, gofmt=False):
    """Write the go source package to the given writer."""
    if gofmt:
        return format_package(writer, package, modules)
    write_header(writer, package, modules)
    write_modules(writer, modules)
    return 0


def write_header(writer, package, modules):
    """Write the go package header."""
    writer.write('// automatically generated from python modules ')
    writer.write(','.join(modules))
    writer.write('\n\n')
    writer.write('package %s\n\n' % package)


def write_modules(writer, modules):
    """Write the go code for each modules."""
    for mod in modules:
        output_package(writer, mod)


def get_type(val, loc):
    """Get the go type for the given value.

    If the value is not one with a good representation in Go, such as a
    list with inconsistent type, returns None.

    """
    if isinstance(val, str):
        return _string_type(val, loc)
    elif isinstance(val, int):
        return 'int'
    elif isinstance(val, long):
        return 'int64'
    elif isinstance(val, float):
        return 'float64'
    elif isinstance(val, bool):
        return 'bool'
    return _get_complex_type(val, loc)

def _string_type(val, loc):
    """Returns 'string' unless the value is a paths based on the runtime path
    of the module, since those will not be accurite."""
    if not val.startswith(loc):
        return 'string'
    return None


def _get_complex_type(val, loc):
    """Recursively gets teh go type name for the given value for sets, lists,
    or maps."""
    if isinstance(val, list):
        eltype = list(set([get_type(v, loc) for v in val]))
        if len(eltype) == 1:
            typename = eltype[0]
            if typename:
                return '[%d]%s' % (len(val), typename)
    elif isinstance(val, set):
        eltype = list(set([get_type(v, loc) for v in val]))
        if len(eltype) == 1:
            typename = eltype[0]
            if typename in _BASE_TYPES:
                return 'map[%s]struct{}' % typename
    elif isinstance(val, dict):
        keytype = list(set([get_type(v, loc) for v in val.keys()]))
        valtype = list(set([get_type(v, loc) for v in val.values()]))
        if len(keytype) == 1 and len(valtype) == 1:
            key_type = keytype[0]
            val_type = valtype[0]
            if key_type in _BASE_TYPES and val_type:
                return 'map[%s]%s' % (key_type, val_type)
    return None


def render(writer, val, tname, indent):
    """Render a variable, recursively."""
    if tname.startswith('map['):
        writer.write(tname)
        if not val:
            writer.write('{}')
            return
        writer.write('{\n')
        key_type = tname[len('map['):tname.index(']')]
        val_type = tname[tname.index(']') + 1:]
        ind = indent + '\t'
        if val_type == 'struct{}':
            _render_list(writer, val, key_type, ind, suffix=': struct{}{},\n')
        else:
            for key, value in val.items():
                writer.write(ind)
                render(writer, key, key_type, ind)
                writer.write(': ')
                render(writer, value, val_type, ind)
                writer.write(',\n')
        writer.write(indent)
        writer.write('}')
    elif tname.startswith('['):
        writer.write(tname)
        if not val:
            writer.write('{}')
            return
        writer.write('{\n')
        val_type = tname[tname.index(']') + 1:]
        _render_list(writer, val, val_type, indent + '\t')
        writer.write(indent)
        writer.write('}')
    elif tname == 'string':
        writer.write('`')
        writer.write(val)
        writer.write('`')
    else:
        writer.write(str(val))


def _render_list(writer, val, val_type, indent, suffix=',\n'):
    """Render the elements of a list."""
    for value in val:
        writer.write(indent)
        render(writer, value, val_type, indent)
        writer.write(suffix)


def output_package(writer, mod):
    """Dump all applicable constants from the given module."""
    _module = __import__(mod, level=0)
    mod_path = mod.split('.')
    for part in mod_path[1:]:
        _module = _module.__dict__[part]
    loc = os.path.dirname(_module.__file__)
    writer.write('// autogenerated from ')
    writer.write(_module.__name__)
    writer.write('\n\n')
    for name, val in sorted(_module.__dict__.items()):
        if name.startswith('_'):
            continue
        type_name = get_type(val, loc)
        if not type_name:
            continue
        if 'map[' in type_name or '[]' in type_name:
            writer.write('var ')
        else:
            writer.write('const ')
        writer.write(name)
        writer.write(' = ')
        render(writer, val, type_name, '')
        writer.write('\n')


if __name__ == '__main__':
    sys.exit(main(sys.argv))
