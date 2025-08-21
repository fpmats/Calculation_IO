# fmt: off

# Copyright (C) 2010, Jesper Friis
# originally from https://gitlab.com/ase/ase/-/blob/master/ase/spacegroup/xtal.py
# modified by Jing T.
# (see accompanying license files for details).

"""
A module for ASE for simple creation of crystalline structures from
knowledge of the space group.

"""

from typing import Any, Dict

import numpy as np
from scipy import spatial
import ase
from ase.geometry import cellpar_to_cell
from ase.symbols import string2symbols
import warnings

class SpacegroupError(Exception):
    """Base exception for the spacegroup module."""

class SpacegroupValueError(SpacegroupError):
    """Raised when arguments have invalid value."""


__all__ = ['crystal']


def crystal(symbols, basis, spacegroup, 
            cell, symprec=0.001, onduplicates='warn', 
            pbc = [True,True,True], size= (1,1,1)) -> ase.Atoms:
    """Create an Atoms instance for a conventional unit cell of a
    space group.

    Parameters:

    symbols : str | sequence of str | sequence of Atom | Atoms
        Element symbols of the unique sites.  Can either be a string
        formula or a sequence of element symbols. E.g. ('Na', 'Cl')
        and 'NaCl' are equivalent.  Can also be given as a sequence of
        Atom objects or an Atoms object.
    basis : array of scaled coordinates
        Positions of the unique sites corresponding to symbols given
        either as scaled positions or through an atoms instance.  
    spacegroup : int | string | Spacegroup instance
        Space group given either as its number in International Tables
        or as its Hermann-Mauguin symbol.
    cell : 3x3 matrix
        Unit cell vectors.
    onduplicates : 'keep' | 'replace' | 'warn' | 'error'
        Action if `basis` contain symmetry-equivalent positions:
            'keep'    - ignore additional symmetry-equivalent positions
            'replace' - replace
            'warn'    - like 'keep', but issue an UserWarning
            'error'   - raises a SpacegroupValueError
    pbc : one or three bools
        Periodic boundary conditions flags.  Examples: True,
        False, 0, 1, (1, 1, 0), (True, False, False).  Default
        is True.
    size : 3 positive integers
        How many times the conventional unit cell should be repeated
        in each direction.
    symprec : float
        Minimum "distance" betweed two sites in scaled coordinates
        before they are counted as the same site.

    Keyword arguments:
    """
    
    if (
            not isinstance(symbols, str) and
            hasattr(symbols, '__getitem__') and
            len(symbols) > 0 and
            isinstance(symbols[0], ase.Atom)):
        symbols = ase.Atoms(symbols)
    if isinstance(symbols, ase.Atoms):
        basis = symbols
        symbols = basis.get_chemical_symbols()
    if isinstance(basis, ase.Atoms):
        basis_coords = basis.get_scaled_positions()
        if symbols is None:
            symbols = basis.get_chemical_symbols()
    else:
        basis_coords = np.array(basis, dtype=float, ndmin=2)

    sites, kinds = spacegroup.equivalent_sites(basis_coords,
                                       onduplicates=onduplicates,
                                       symprec=symprec)

    # this is needed to handle deuterium masses
    symbols = parse_symbols(symbols)
    symbols = [symbols[i] for i in kinds]
    
    atoms = ase.Atoms(symbols,
                      scaled_positions=sites,
                      cell=cell,
                      pbc=pbc)

    if isinstance(basis, ase.Atoms):
        for name in basis.arrays:
            if not atoms.has(name):
                array = basis.get_array(name)
                atoms.new_array(name, [array[i] for i in kinds],
                                dtype=array.dtype, shape=array.shape[1:])

    if kinds:
        atoms.new_array('spacegroup_kinds', np.asarray(kinds, dtype=int))

    if size != (1, 1, 1):
        atoms = atoms.repeat(size)
    return atoms


def parse_symbols(symbols):
    """Return `sumbols` as a sequence of element symbols."""
    if isinstance(symbols, str):
        symbols = string2symbols(symbols)
    return symbols