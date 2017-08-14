r"""Legacy support for the SDF python module.

The syntax and behaviour of the SDF module has changed substantially with
version 2.0.0 of the module and existing scripts will no longer work
without modification.

This python module provides a backwards compatibility layer so that old
scripts can be made to work without any modification. To use it, add the
following lines to the top of your script:

try:
    import sdf_legacy as sdf
except ImportError:
    import sdf
"""

import os.path
import sdf

__version__ = "2.0.0"
_module_name = "sdf_legacy"
_sdf_version = __version__


def _error_message():
    raise ImportError(r"Your sdf python module is too old for this version "
                      "of " + _module_name + ".\n"
                      "Either upgrade to sdf python " + _sdf_version +
                      " or newer, or downgrade " + _module_name)


def _check_validity():
    if not hasattr(sdf, "__version__"):
        return _error_message()
    our_version = map(int, _sdf_version.split("."))
    lib_version = map(int, sdf.__version__.split("."))
    # Just check that major version number matches
    if our_version[0] != lib_version[0]:
        return _error_message()

_check_validity()


class SDFObject:
    filename = None
    convert = 0

    def __init__(self, filename, convert):
        self.filename = filename
        self.convert = convert

    def read(self, *args, **kwargs):
        """Reads the SDF file associated with this object and returns a
        dictionary of NumPy arrays containing the file data."""
        return read(self.filename, convert=self.convert, *args, **kwargs)


def SDF(filename, convert=0):
    """Tests for the existence of a file and returns a handle object"""

    # Test for file existence and raise an exception if it isn't readable
    fd = os.open(filename, os.O_RDONLY)
    os.close(fd)
    return SDFObject(filename, convert)


def read(*args, **kwargs):
    """Reads an SDF file and returns a dictionary of NumPy arrays containing
       the file data."""

    data = sdf.read(*args, **kwargs)

    sdfdict = {}
    for key, value in data.__dict__.items():
        if hasattr(value, "name"):
            if hasattr(value, "data"):
                sdfdict[value.name] = value.data
            else:
                sdfdict[value.name] = value
        else:
            if hasattr(value, "data"):
                sdfdict[key] = value.data
            else:
                sdfdict[key] = value

        # Add separate grid variables
        t = type(value)
        if t == sdf.BlockPlainMesh or t == sdf.BlockPointMesh \
                or t == sdf.BlockLagrangianMesh:
            base = value.name
            if value.id == "grid":
                base += "_node"
            elif value.id == "grid_mid":
                base = base[:-4]
            for n in range(len(value.dims)):
                newkey = base + '/' + value.labels[n]
                sdfdict[newkey] = value.data[n]

        # Add material block
        if t == sdf.BlockStitchedMaterial:
            sdfdict["Materials"] = value.material_names

    return sdfdict
