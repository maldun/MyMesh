# MyMesh Module - API for easier Salome smesh usage
# Basics.py: Smesh helper functions and classes for basic operations 
#
# Copyright (C) 2013  Stefan Reiterer - maldun.finsterschreck@gmail.com
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

from __future__ import print_function

import salome
import geompy
import GEOM
import smesh
import SMESH

from smesh import GetFilter, EDGE, FACE, VOLUME, FT_LinearOrQuadratic, Geom_TRIANGLE, Geom_QUADRANGLE

from numpy import array, ndarray, arange
from numpy import float64 as data_type

from Types import Element, Tria3, Quad4

def find_mesh(descriptive_string):
    """
    Help function to find mesh in a study and converts it to smesh mesh
    """
    mesh = salome.myStudy.FindObject(descriptive_string).GetObject()
    return smesh.Mesh(mesh,descriptive_string)


def apply_linear_elements(mesh,elem_ids):
    """
    Give the elements in the list the right container classes
    """
    elem_list = []
    for elem in elems:
        if Mesh.GetElementGeomType(elem) == Geom_TRIANGLE:
            elem_list += [Tria3(elem)]
        elif Mesh.GetElementGeomType(elem) == Geom_QUADRANGLE:
            elem_list += [Quad4(elem)]
    
    return elem_list
