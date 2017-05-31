#!/usr/bin/python
# -*- coding: utf-8 -*-

# MyMesh Module - API for easier Salome smesh usage
# Basics.py: Smesh helper functions and classes for basic operations 
#
# Copyright (C) 2015  Stefan Reiterer - stefan.harald.reiterer@gmail.com
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

# import salome
# import geompy
# import GEOM
# import smesh
# import SMESH

import salome
salome.salome_init()
import GEOM
from salome.geom import geomBuilder
geompy = geomBuilder.New(salome.myStudy)

import SMESH, SALOMEDS
from salome.smesh import smeshBuilder
smesh =  smeshBuilder.New(salome.myStudy)


#from smesh import GetFilter
GetFilter = smesh.GetFilter
from SMESH import Entity_Triangle, Entity_Quadrangle
from SMESH import EDGE, FACE, VOLUME, FT_LinearOrQuadratic, Geom_TRIANGLE, Geom_QUADRANGLE

from numpy import array, ndarray, arange, cross, inner, zeros
from numpy.linalg import norm
from numpy import float64 as data_type
from numpy import arccos, tan, pi

from Types import Element, Tria3, Quad4

DIMENSION = 3

def find_mesh(descriptive_string):
    """
    Help function to find mesh in a study and converts it to smesh mesh
    """
    mesh = salome.myStudy.FindObject(descriptive_string).GetObject()
    return smesh.Mesh(mesh) #,descriptive_string)

def load_mesh_from_file(filename):
    """
    Help function wich loads a mesh from a file.
    """
    return smesh.CreateMeshesFromMED(filename)[0][0]


def apply_linear_elements(mesh,elem_ids):
    """
    Give the elements in the list the right container classes
    """
    elem_list = []
    for elem in elem_ids:
        if mesh.GetElementGeomType(elem) == Entity_Triangle:
            elem_list += [Tria3(mesh,elem)]
        elif mesh.GetElementGeomType(elem) == Entity_Quadrangle:
            elem_list += [Quad4(mesh,elem)]

    return elem_list

def compute_voroni_area_of_triangle(w1,w2,l1,l2):
    """
    computes the part of the triangle for the Voroni area
    descriped in [1]

            x_i
          +
          |\
          | \
          |  \
        l1|   \l2
          |    \
          |     \
          |w1  w2\
          +-------+
        x_j  l3  x_{j+1}

    """
    # Check if triangle is obtuse in x_i
    if w1+w2 < pi/2.0:
        # Then get Area(T)/2
        return norm(cross(l1,l2))/4.0
    
    if w1 > pi/2.0 or w2 > pi/2:
        # Then get Area(T)/4
        return norm(cross(l1,l2))/8.0
            

    #Else use formula on page 9 in [1]
    return ((1/tan(w1))*inner(l2,l2) + (1/tan(w2))*inner(l1,l1))/8.0

def compute_gravity_center(mesh,group=None,element_ids = []):
    u"""
    Computes the center of gravity of a face mesh (or of a face group), by
    the discrete formula

        (Σₑ area(e)·Sₑ)/(Σₑ area(e)).
    """
    # filter linear and triangle elements
    filter_linear_tri = GetFilter(FACE, FT_LinearOrQuadratic, Geom_TRIANGLE)
    filter_linear_quad = GetFilter(FACE, FT_LinearOrQuadratic, Geom_QUADRANGLE)
    
    ids_tri = mesh.GetIdsFromFilter(filter_linear_tri)
    ids_quad = mesh.GetIdsFromFilter(filter_linear_quad)

    if group is None and element_ids == []:
        tria3 = [Tria3(mesh,id_tri) for id_tri in ids_tri]
        quad4 = [Tria3(mesh,id_tri) for id_tri in ids_tri]
    elif group is None and element_ids != []:
        group_ids = set(element_ids)
        tria3 = [Tria3(mesh,id_tri) for id_tri in ids_tri if id_tri in group_ids]
        quad4 = [Tria3(mesh,id_tri) for id_tri in ids_tri if id_tri in group_ids]
    else:
        group_ids = set(group.GetIDs())
        tria3 = [Tria3(mesh,id_tri) for id_tri in ids_tri if id_tri in group_ids]
        quad4 = [Tria3(mesh,id_tri) for id_tri in ids_tri if id_tri in group_ids]

    elems = tria3 + quad4
    S = zeros(DIMENSION)
    A = 0
    for e in elems:
        a = e.computeArea(store=False)
        S += a*e.computeGravityCenter()
        A += a

    return S/A
        
