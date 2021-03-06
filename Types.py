#!/usr/bin/python
# -*- coding: utf-8 -*-

# MyMesh Module - API for easier Salome smesh usage
# Basics.py: Smesh helper functions and classes for basic operations 
#
# Copyright (C) 2013  Stefan Reiterer - stefan.harald.reiterer@gmail.com
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

"""
REFERENCES:

-    [1] 'Discrete Differential-Geomety Operators for Triangulated 2-Manifolds'
          by Meyer, Desbrun, Schroeder and Barr. (2002) Visualization and mathematics, 3(2), 52-58.
"""

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
from SMESH import EDGE, FACE, VOLUME, FT_LinearOrQuadratic, Geom_TRIANGLE, Geom_QUADRANGLE

from numpy import array, ndarray, arange, cross, zeros, inner, append
from numpy import sum, apply_along_axis, copy
from numpy.linalg import norm
from numpy import float64 as data_type
from numpy import arccos, tan, pi
import numpy as np

# Nr. of dimension
DIMENSION=3

class Element(object):
    """
    Help class to handle elements of a mesh.
    The Intention is to use this class to add more
    OO programming into the smesh Module.
    smesh is a well designed module, much better
    than geompy. Nevertheless we will use
    This objects to enable better OO and functional
    programming capabilities. But we will try to minimize
    the use of self made tools as good as we can.

    Checks are avoided due to the use of filters etc. so be careful.
    """
    def __init__(self,mesh,id_nr):
        """
        An element is defined within it's mesh. We don't want to
        create a new Element, only a reference.
        For that reason we hold the id number and the mesh which belong to
        this element.
        """
        
        self.mesh = mesh
        self.idNr = id_nr  
        
    def getIdNr(self):
        return self.idNr

    def getMesh(self):
        return self.Mesh

class FaceElement(Element):
    """
    Template class for faces
    """
    def __init__(self,mesh,id_nr):
        super(FaceElement,self).__init__(mesh,id_nr)
        self.normals = {}
        self.voroni_area = None 

    def getNodes(self):
        return self.mesh.GetElemNodes(self.getIdNr())

    def computeNormal(self,node,store = True):
        """
        Template method for normal computation.
        If store is True, the normal will be stored
        in a dict where the node is the key.
        """
        raise NotImplementedError("Error: Not implemented!")

    def getNormals(self):
        return self.normals

    def getNormal(self,node):
        return self.normals[node]
        
    def computeArea(self,store = True):
        raise NotImplementedError("Error: Not implemented!")

    def getArea(self):
        return self.area

    def getCurvatureVector(self,node,voroni=False):
        raise NotImplementedError("Error: Not implemented!")

    def computeGravityCenter(self):
        pass


    
class Tria3(FaceElement):
    """
    Simple linear triangles with 3 corners
    """
    
    def _computeNormalOp(self):
        """
        Computes the normal of the triangle.
        The returned vector is not normalized!
        This function only executes the computation.
        The key value for node
        """

        nodes = self.getNodes()
        p0 = self.mesh.GetNodeXYZ(nodes[0])
        p1 = self.mesh.GetNodeXYZ(nodes[1])
        p2 = self.mesh.GetNodeXYZ(nodes[-1])
        
        vec1 = array(p1) - array(p0)
        vec2 = array(p2) - array(p0)
        vec_normal = cross(vec1, vec2)
        return vec_normal, norm(vec_normal)

    def computeNormal(self,node,store = True):
        """
        It stores the vector in the object per default.
        """
        if store:
            self.normals[node] = self._computeNormalOp()[0]
        else:
            return self._computeNormalOp()[0]

    def computeArea(self,store = True):
        """
        It stores the vector in the object per default.
        """
        if store:
            self.area = self._computeNormalOp()[1]/2.0
        else:
            return self._computeNormalOp()[1]/2.0

    def _computeVoroniArea(self,w1,w2,l1,l2):
        """
        computes the part of the triangle for the Voroni area
        descriped in [1]
        """
        from Tools import compute_voroni_area_of_triangle
        return compute_voroni_area_of_triangle(w1,w2,l1,l2)
        
    def computeCurvatureVector(self,node,voroni=False):
        """
        Computes the vector required for the mean curvature normal formula.
        in a triangle. There are the following formuli: (x_i is node; x_j, x_{j+1} are
        the other nodes)). We use the description in [1].

        \cos(\alpha_{ij}) = < x_i - x_{j+1}, x_j - x_{j+1}>/norms
        \cos(\beta_{ij+1}) = < x_i - x_j, x_{j+1} - x_j>/norms

        return \cot(\alpha_{ij})(x_i - x_j) + \cot(\beta_{ij+1}) (x_i - x_{j+1})

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
        nodes = self.getNodes()
        index_node = nodes.index(node)

        x_i = array(self.mesh.GetNodeXYZ(node))
        x_j = array(self.mesh.GetNodeXYZ(nodes[index_node-1]))
        x_jp = array(self.mesh.GetNodeXYZ(nodes[index_node-2]))

        l1 = x_i - x_j
        l2 = x_i - x_jp
        l3 = x_j-x_jp
        
        w1 = arccos(inner(l1,-l3)/(norm(l1)*norm(l3)))
        w2 = arccos(inner(l2,l3)/(norm(l2)*norm(l3)))
        
        if voroni:
            return (1/tan(w1))*l2 + (1/tan(w2))*l1, self._computeVoroniArea(w1,w2,l1,l2)
        
        return (1/tan(w1))*l2 + (1/tan(w2))*l1

    def computeGravityCenter(self):
        """
        Computes the center of gravity (from the zero vector) by the well known forumla

        math::
        
          (x_i + x_{i+1} + x_{xi+2})/3
          
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
        nodes = self.getNodes()
        coord_matrix = array([self.mesh.GetNodeXYZ(node) for node in nodes])
        return apply_along_axis(sum,0,coord_matrix)/3.0
    
class Quad4(FaceElement):
    """
    Simple linear quadrangles with 4 corners
    """
    
    def _computeNormalOp(self):
        """
        Computes the normal of the triangle.
        The returned vector is not normalized
        This function only executes the computation
        """

        nodes = self.mesh.GetElemNodes(self.getIdNr())
        p0 = self.mesh.GetNodeXYZ(nodes[0])
        p1 = self.mesh.GetNodeXYZ(nodes[1])
        p2 = self.mesh.GetNodeXYZ(nodes[-1])
        
        vec1 = array(p1) - array(p0)
        vec2 = array(p2) - array(p0)
        vec_normal = cross(vec1, vec2)
        return vec_normal, norm(vec_normal)

    def computeNormal(self,node,store = True):
        """
        It stores the vector in the object per default.
        """
        if store:
            self.normals[node] = self._computeNormalOp()[0]
        else:

            return self._computeNormalOp()[0]

    def computeArea(self,store = True):
        """
        It stores the vector in the object per default.
        """
        if store:
            self.area = self._computeNormalOp()[1]
        else:
            return self._computeNormalOp()[1]

    def _computeVoroniArea(self,w1,w2,w3,w4,l1,l2,l3):
        """
        Computes the the Voroni area of the quadrangle by splitting it into two
        triangles like described in the computeCurvatureVector method.
        """
        from Tools import compute_voroni_area_of_triangle
        voroni_area = compute_voroni_area_of_triangle(w1,w2,l1,l2)
        voroni_area += compute_voroni_area_of_triangle(w3,w4,l2,l3)
        return voroni_area
        
    def computeCurvatureVector(self,node,voroni=False):
        """
        Computes the vector required for the mean curvature normal formula.
        in a quadrangle. There are the following formuli: (x_i is node; x_j, x_{j+1} and x_{j+2} 
        are the other nodes). We interprete the quadrangle as a halfed triangle to derive the 
        formulas. In order to avoid obtuse triangles we distinguish between 2 cases. If 
        the angle in x_i is larger than the angles in the other 2 corners we split the quadrangle 
        this way:
        x_i    l3  x_{j+2}
           +-------+
           |\    w4|
           | \     |
           |  \l2  |
         l1|   \   |l5
           |    \w3|
           |     \ |
           |w1  w2\|
           +-------+
        x_j    l4  x_{j+1}

        Else we forget about the lower half and split this way:
        x_i    l3  x_{j+2}
           +-------+
           |    v2/:
           |     / :
           |  l6/  :
         l1|   /   :l5
           |v1/    :
           | /     :
           |/      :
           +.......+
        x_j    l4  x_{j+1}

        This handling is mathematically correct, since we could view it as a local mesh 
        refinement in order to get a better approximation in a certain area.
        """
        nodes = self.getNodes()
        index_node = nodes.index(node)

        x_i = array(self.mesh.GetNodeXYZ(node))
        x_j = array(self.mesh.GetNodeXYZ(nodes[index_node-1]))
        x_jp = array(self.mesh.GetNodeXYZ(nodes[index_node-2]))
        x_jp2 = array(self.mesh.GetNodeXYZ(nodes[index_node-3]))

        l1 = x_i - x_j
        l2 = x_i - x_jp
        l3 = x_i - x_jp2
        l4 = x_j-x_jp
        l5 = x_jp2 - x_jp 
        
        w1 = arccos(inner(l1,-l4)/(norm(l1)*norm(l4)))
        w2 = arccos(inner(l2,l4)/(norm(l2)*norm(l4)))
        w3 = arccos(inner(l2,l5)/(norm(l2)*norm(l5)))
        w4 = arccos(inner(l3,-l5)/(norm(l3)*norm(l5)))

        # first case: 
        if (w2 + w3) > pi/2.0 or (w1+w2+w3+w4) < 3.0*pi/2.0:

            if voroni:
                return ((1/tan(w1))+ (1/tan(w4)))*l2 + (1/tan(w2))*l1 + (1/tan(w3))*l3,\
                  self._computeVoroniArea(w1,w2,w3,w4,l1,l2,l3)

            return ((1/tan(w1))+ (1/tan(w4)))*l2 + (1/tan(w2))*l1 + (1/tan(w3))*l3  
        # second case
        else:
            l6 = x_jp2 - x_j
            v1 = arccos(inner(l1,l6)/(norm(l1)*norm(l6)))
            v2 = arccos(inner(l3,-l6)/(norm(l3)*norm(l6)))

            if voroni:
                from Tools import compute_voroni_area_of_triangle
                return (1/tan(v1))*l3 + (1/tan(v2))*l1,\
                  compute_voroni_area_of_triangle(v1,v2,l1,l3)
            
            return (1/tan(v1))*l3 + (1/tan(v2))*l1


    def computeGravityCenter(self):
        """
        Computes the center of gravity of a quadrangle.
        in a quadrangle. We have the quadrangle (x_j, x_{j+1} x_{j+2} and x_{j+3}) .
        We interprete the quadrangle as two halfed triangles to derive the 
        center. first we compute the center of gravity for the two triangles, and
        then we get center by taking the weighted mean.
        
        x_j l3 x_{j+3}
           +-------+
           |\    w4|
           | \     |
           |  \l2  |
         l1|   \   |l5
           |    \w3|
           |     \ |
           |w1  w2\|
           +-------+
        x_{j+1} l4 x_{j+2}

        """
        nodes = self.getNodes()

        # the strange arangement is historically
        x_j = array(self.mesh.GetNodeXYZ(nodes[0]))
        x_jp1 = array(self.mesh.GetNodeXYZ(nodes[1]))
        x_jp2 = array(self.mesh.GetNodeXYZ(nodes[2]))
        x_jp3 = array(self.mesh.GetNodeXYZ(nodes[3]))
        
        S1 = (x_j + x_jp1 + x_jp2)/3.0
        S2 = (x_j + x_jp2 + x_jp3)/3.0
        area1 = norm(cross(x_j-x_jp1,x_jp2-x_j))
        area2 = norm(cross(x_j-x_jp3,x_jp2-x_jp3))
        return (S1*area1+S2*area2)/(area1+area2)
        
class VectorField(object):
    """
    Class that represents a vector field which is applied on the nodes of 
    mesh. Thus is holds also the information about the mesh, on which it works,
    nad on which group it may be restricted
    """

    def __init__(self,mesh, scalar = 1.0, restricted_group=None):
        
        self.mesh = mesh
        self.scalar = scalar
        self.restricted_group = restricted_group

        self._applied_extrusions = 0

    def getVectorOnNode(self,node_id):
        return self.computeVectorOnNode(node_id)*self.getScalar()

    def computeVectorOnNode(self,node_id):
        raise NotImplementedError('Error: This method is a stub')

    def getScalar(self): return self.scalar

    def setScalar(self,scalar): self.scalar = scalar

    def getRestricedGroup(self): return self.restricted_group

    def setRestricedGroup(self,restricted_group): self.restricted_group = restricted_group
    
    def scalarMultiplication(self,scalar):
        """
        Implements multiplication with scalar.
        """
        mult_scalar = self.getScalar()*scalar
        self.setScalar(mult_scalar)

    def __rmul__(self,scalar): 
        
        self.scalarMultiplication(scalar)
        return self

    def updateNode(self,node_id,new_node): 
        """
        Applies necessary update operations.
        """
        pass

    def _processLookupTables(self,lookup_table):
        """
        Stub method for processing steps with
        lookup tables.
        """
    
    def applyVectorOnNode(self,node_id, mesh = None):
        """
        Apply the vector field on a given node. This method creates the
        displaced node in the given mesh, and returns the index of the new node.
        The node is created in the current mesh per default, but it is possible 
        to add it to a different mesh.

        Arguments:
        - `self`: 
        - `node_id`: integer with the Id of the node we want to apply the vector.
        """
        
        if mesh is None:
            mesh = self.mesh 

        node_vec = array(self.mesh.GetNodeXYZ(node_id))
        translated_node_vec = node_vec + self.getVectorOnNode(node_id)
        new_node = mesh.AddNode(*translated_node_vec.tolist())

        self.updateNode(node_id,new_node)
        
        return new_node

    def applyVectorFieldOnFace(self,face,mesh=None, table = {}):
        """
        Applies the Vector field on a given surface. The Face is either saved in
        the current mesh (default) or can be saved in a different mesh

        Arguments:
        - `self`:
        - `face`: The face to translate
        - `mesh`: the mesh where the face should belong to. Default is self.mesh
        - `table`: returns a dict with the corespondence table of the nodes.

        returns: The id of the new face
        """
        if mesh is None:
            mesh = self.mesh

        elem_nodes = self.mesh.GetElemNodes(face)
        new_nodes = []
        for node in elem_nodes:
            try:
                new_nodes += [table[node]]
            except KeyError:
                new_nodes += [self.applyVectorOnNode(node,mesh)]

        new_face = mesh.AddFace(new_nodes)

        return new_face, new_nodes
        

    def computeSurfaceExtrusion(self,group=None, edge_groups = [], face_groups = []):
        """
        This method applies the vector field on a surface and creates
        a translated one. This method makes the computation only.
        The method surfaceExtrusion also adds the groups to the mesh.

        Arguments:
        - `self`: 
        - `mesh`: Optional smesh.Mesh instance. Per defualt it is self.mesh
        - `group`: Optional group of elements on which we apply the vector field.
        - `table`: Variable if lookup table should be returned for further steps.
        """
        

        mesh = self.mesh

        if group is None:
            faces = self.mesh.GetElementsByType(FACE)
        else:
            faces = group.GetIDs()

        lookup_table = {}
        new_face_ids = []
        new_vol_ids = []

        for face in faces:

            face_id, new_nodes = self.applyVectorFieldOnFace(face,mesh,lookup_table)
            new_face_ids += [face_id]

            # update lookup table
            elem_nodes = self.mesh.GetElemNodes(face)
            new_lookup = [[elem_nodes[i],new_nodes[i]] for i in range(len(elem_nodes))]
            lookup_table.update(new_lookup)
            # now add volume elements
            new_vol_ids += [mesh.AddVolume(elem_nodes + new_nodes)]


        if edge_groups:
            edge_groups_faces = []
            new_edge_groups = []
            for edge_group in edge_groups:
                edge_group_faces = []
                new_edges = []
                for edge in edge_group.GetIDs():
                    edge_nodes = mesh.GetElemNodes(edge)
                    new_edge_nodes = [lookup_table[node] for node in edge_nodes]
                    new_edges += [mesh.AddEdge(new_edge_nodes)]
                    
                new_edge_groups += [new_edges]
                edge_groups_faces += [edge_group_faces]
            
        if face_groups:
            new_face_groups = []

            for face_group in face_groups:
                face_group_ids = face_group.GetIDs()
            
                new_ids = [new_face_ids[position] for position in xrange(len(faces)) if faces[position] in face_group_ids]
                new_face_groups += [new_ids]
        
        if edge_groups and not face_groups:
            return new_face_ids, new_vol_ids, new_edge_groups, lookup_table
        elif not edge_groups and face_groups:
            return new_face_ids, new_vol_ids, new_face_groups, lookup_table
        elif edge_groups and face_groups:
            return new_face_ids, new_vol_ids, new_edge_groups, new_face_groups, lookup_table
        else:
            return new_face_ids, new_vol_ids, lookup_table
    
    def extrudeSurface(self,group=None, edge_groups = [], face_groups = [], create_boundary_elements = True):
        """
        This method applies the vector field on a surface and creates
        a translated one. 

        Arguments:
        - `self`: 
        - `mesh`: Optional smesh.Mesh instance. Per defualt it is self.mesh
        - `group`: Optional group of elements on which we apply the vector field.
        - `table`: Variable if lookup table should be returned for further steps.
        """

        mesh = self.mesh

        self._applied_extrusions += 1 

        if edge_groups and not face_groups:
            new_face_ids, new_vol_ids, new_edge_groups, lookup_table = self.computeSurfaceExtrusion(group=group, edge_groups = edge_groups)
        elif not edge_groups and face_groups:
            new_face_ids, new_vol_ids, new_face_groups, lookup_table = self.computeSurfaceExtrusion(group=group, face_groups = face_groups)
        elif edge_groups and face_groups:
            new_face_ids, new_vol_ids, new_edge_groups, new_face_groups, lookup_table = self.computeSurfaceExtrusion(group=group, edge_groups = edge_groups, face_groups = face_groups)
        else:
            new_face_ids, new_vol_ids, lookup_table = self.computeSurfaceExtrusion(group=group, edge_groups = edge_groups)
        # add face and volume group    
        if group:
            face_group = mesh.MakeGroupByIds(group.GetName()+'_extruded_faces' + str(self._applied_extrusions),FACE,new_face_ids)
            vol_group = mesh.MakeGroupByIds(group.GetName()+'_extruded_volumes' + str(self._applied_extrusions),VOLUME,new_vol_ids)
        else:
            face_group = mesh.MakeGroupByIds(mesh.GetName() + '_extruded_faces' + str(self._applied_extrusions),FACE,new_face_ids)
            vol_group = mesh.MakeGroupByIds(mesh.GetName() + '_extruded_volumes' + str(self._applied_extrusions),VOLUME,new_vol_ids)
        
        if create_boundary_elements:
            bnd_faces = mesh.MakeBoundaryMesh(vol_group)
        if edge_groups:
            
            #bnd_faces = mesh.MakeGroupByIds("boundary_faces",FACE,bnd_faces)
            #bnd_nodes = bnd_faces.GetNodeIDs()

            new_salome_edge_groups = []
            new_salome_edge_face_groups = []
            for i in range(len(edge_groups)):
                new_salome_edge_groups += [mesh.MakeGroupByIds(edge_groups[i].GetName()+'_extruded' + str(self._applied_extrusions),EDGE,new_edge_groups[i])]
                
                if create_boundary_elements:
                    #find faces belonging to new created edge groups
                    new_edge_group_faces = []
                    old_edges = edge_groups[i].GetIDs()
                    new_edges = new_edge_groups[i]
                    for nr_edge in range(len(old_edges)):
                    
                        nodes_edge = mesh.GetElemNodes(old_edges[nr_edge])
                        nodes_new_edge = mesh.GetElemNodes(new_edges[nr_edge])
                        new_edge_group_faces += [mesh.FindElementByNodes(nodes_edge + nodes_new_edge)]
                    
                    new_salome_edge_face_groups += [mesh.MakeGroupByIds(edge_groups[i].GetName()+'_extruded_faces' + str(self._applied_extrusions),FACE,new_edge_group_faces)]

        if face_groups:
            new_salome_face_groups = []
            for i in range(len(face_groups)):
                new_salome_face_groups += [mesh.MakeGroupByIds(face_groups[i].GetName()+'_extruded' + str(self._applied_extrusions),FACE,new_face_groups[i])]

        if edge_groups and not face_groups:
            salome.sg.updateObjBrowser(0)
            return face_group, vol_group, new_salome_edge_groups, new_salome_edge_face_groups, [], bnd_faces, lookup_table 
        elif not edge_groups and face_groups:
            salome.sg.updateObjBrowser(0)
            return face_group, vol_group, [], [], new_salome_face_groups, bnd_faces, lookup_table
        elif edge_groups and face_groups:
            salome.sg.updateObjBrowser(0)
            return face_group, vol_group, new_salome_edge_groups, new_salome_edge_face_groups, new_salome_face_groups, bnd_faces, lookup_table
        else:
            salome.sg.updateObjBrowser(0)
            return face_group, vol_group, [], [], [], bnd_faces,lookup_table
        
    def extrudeSurfaceTimes(self,k,group=None, edge_groups = [], face_groups = []):
        """
        This method applies the vector field on a surface and creates
        a translated one over k steps.

        Arguments:
        - `self`: 
        - `k`: number of extrusions, or list of extrusion thicknesses
        - `mesh`: Optional smesh.Mesh instance. Per defualt it is self.mesh
        - `group`: Optional group of elements on which we apply the vector field.
        - `table`: Variable if lookup table should be returned for further steps.
        """
        if isinstance(k,list):
            thicknesses = k
            k = len(k)
        else:
            thicknesses = []
        
        face_groups_rest = [None]*k
        extruded_face_groups = [None]*k
        vol_groups = [None]*k
        extruded_edge_groups = [None]*k
        extruded_edge_surface_groups = [None]*k
        bnd_face_groups = [None]*k
        lookup_tables = [None]*k
        
        face_groups_rest[-1] = group
        extruded_edge_groups[-1] = edge_groups
        extruded_face_groups[-1] = face_groups

        original_rst_group = self.getRestricedGroup()

        for i in range(k):
            if thicknesses:
                self.setScalar(thicknesses[i])


            face_groups_rest[i], vol_groups[i], extruded_edge_groups[i], extruded_edge_surface_groups[i], extruded_face_groups[i], bnd_face_groups[i], lookup_tables[i] = self.extrudeSurface(group = face_groups_rest[i-1], edge_groups = extruded_edge_groups[i-1],face_groups = extruded_face_groups[i-1])
            
            self.setRestricedGroup(face_groups_rest[i])
            self._processLookupTables(lookup_tables)
                
        self.setRestricedGroup(original_rst_group)
        return face_groups_rest, vol_groups, extruded_edge_groups, extruded_edge_surface_groups, extruded_face_groups, bnd_face_groups, lookup_tables
            

    def applyVectorFieldOnSurface(self,mesh=None,group=None):
        """
        This method applies the vector field on a surface and creates
        a translated one. Optional the new surface can be stored in a new mesh.

        Arguments:
        - `self`: 
        - `mesh`: Optional smesh.Mesh instance. Per defualt it is self.mesh
        - `group`: Optional group of elements on which we apply the vector field.
        """
        
        if mesh is None:
            mesh = self.mesh

        if group is None:
            nodes = self.mesh.GetNodesId()
            faces = self.mesh.GetElementsByType(FACE)
        else:
            nodes = group.GetNodeIDs()
            faces = group.GetIDs()

        lookup_table = {}
        new_face_ids = []
        for face in faces:

            face_id, new_nodes = self.applyVectorFieldOnFace(face,mesh,lookup_table)
            new_face_ids += [face_id]

            # update lookup table
            elem_nodes = self.mesh.GetElemNodes(face)
            new_lookup = [[elem_nodes[i],new_nodes[i]] for i in range(len(elem_nodes))]
            lookup_table.update(new_lookup)

        return new_face_ids, lookup_table

    def computeNewPosition(self,node_id):
        """
        This method applies the vector field on a node and computes the position
        of a new node moved by the vector field.

        Arguments:
        - `self`: 
        - `node_id`: id of the node to be translated.
        """
        vector = self.computeVectorOnNode(node_id)
        coords = array(self.mesh.GetNodeXYZ(node_id))
        coords += vector
        return coords
        
    def moveNode(self,node_id):
        """
        This method applies the vector field on a node and moves it without creating
        a new node.

        Arguments:
        - `self`: 
        - `node_id`: id of the node to be translated.
        """
        coords = self.computeNewPosition(node_id)
        self.mesh.MoveNode(node_id,*coords)
        
    def moveSurface(self,group=None):
        """
        This method applies the vector field on a surface and creates
        a translated one. Optional the new surface can be stored in a new mesh.

        Arguments:
        - `self`: 
        - `mesh`: Optional smesh.Mesh instance. Per defualt it is self.mesh
        - `group`: Optional group of elements on which we apply the vector field.
        """
        
        if group is None:
            nodes = self.mesh.GetNodesId()
        else:
            nodes = group.GetNodeIDs()

        # To avoid wrong computations the vectors have to be computed
        # before the nodes are moved 
        new_node_positions = [[node] + self.computeNewPosition(node).tolist() for node in nodes]
        for position in new_node_positions: 
            self.mesh.MoveNode(*position)
        
class NormalVectorField(VectorField):

    """
    Class to compute the normal vector field on a mesh.
    The current version only works on surfaces, and on
    linear elements. Hopefully this will change.
    """
    def __init__(self,mesh,scalar = 1.0, restricted_group=None,fast = True):
        
        super(NormalVectorField,self).__init__(mesh,scalar,restricted_group)

        # filter linear and triangle elements
        filter_linear_tri = GetFilter(FACE, FT_LinearOrQuadratic, Geom_TRIANGLE)
        filter_linear_quad = GetFilter(FACE, FT_LinearOrQuadratic, Geom_QUADRANGLE)

        ids_tri = mesh.GetIdsFromFilter(filter_linear_tri)
        ids_quad = mesh.GetIdsFromFilter(filter_linear_quad)

        self.tria3 = [Tria3(mesh,id_tri) for id_tri in ids_tri]
        self.quad4 = [Tria3(mesh,id_tri) for id_tri in ids_tri]

        if sum((len(self.tria3),len(self.quad4))) is 0:
            raise NotImplementedError("Error: Type of element not implemented!")

        self.fast = fast
        
        if not fast:
            bound_filter = smesh.GetFilter(SMESH.EDGE, SMESH.FT_FreeBorders)
            bound_ids = mesh.GetIdsFromFilter(bound_filter)
            bound_nodes = [mesh.GetElemNodes(edge)[0] for edge in bound_ids]
            bound_nodes += [mesh.GetElemNodes(edge)[1] for edge in bound_ids]

            self.bound_nodes = set(bound_nodes)
        
        if restricted_group is not None:
            self.rst_group = set((self.getRestricedGroup()).GetIDs())
        else:
            self.rst_group = None

    def switchComputationMode(self): self.fast = not self.fast
            
    def setRestricedGroup(self,restricted_group):
        self.restricted_group = restricted_group
        if restricted_group is not None:
            self.rst_group = set((self.getRestricedGroup()).GetIDs())

    def meanNormalFormula(self,elems,node_id):
        """
        Compute the mean normal, whith the formulas

        math:: 

            m(x) = (\sum_{\Delta \in Faces(x)} ||n|| n_{\Delta}(x))/\#Faces(x)
            n(x) = m(x)/||m(x)||
        """

        result = zeros(DIMENSION)

        for elem in elems:
            result += elem.computeNormal(node_id,store=False)

        result /= len(elems)
        return result/norm(result)


    def meanCurvatureNormalFormula(self,elems,node_id,voroni=False):
        """
        Compute an averaged mean curvature normal
        due to the formula

        math::

            \frac{1}{2} \sum_{j \in N_1(i)} (\cot \alpha_{ij} + \cot \beta_{ij})(x_j - x_i)

        where x_j are the neighbouring nodes of x_i, and alpha_{ij}, beta_{ij} are the angles
        oposite to the edge (x_i,x_j) supposing we have a triangle mesh. If quad elements appear
        we cut the quadrngles in a half.
        (see [1] for more details) 
        """

        normal = zeros(DIMENSION)
        if voroni is False:
            for elem in elems:
                normal += elem.computeCurvatureVector(node_id)

        else:
            voroni_area = 0
            for elem in elems:
                vec, area = elem.computeCurvatureVector(node_id,voroni=True)
                normal += vec
                voroni_area += area 

            normal *= 1.0/(2*voroni_area)

                
        return normal
        
    
    def computeVectorOnNode(self,node_id):
        """
        We compute the normal in one node,
        with help of the mean curvature normal formula.
        If the curvature is too small we apply the 
        mean normal formula.
        """
        Mesh = self.mesh

        elems = Mesh.GetNodeInverseElements(node_id)

        rst_group = self.rst_group
        if rst_group is not None:
            elems = [elem for elem in elems if elem in rst_group]

        from Tools import apply_linear_elements
        elems = apply_linear_elements(Mesh,elems)

        if not self.fast:
            if not node_id in self.bound_nodes:
                result = self.meanCurvatureNormalFormula(elems,node_id)
                if norm(result) >= 1e-2:
                    return result/norm(result)
        
        return self.meanNormalFormula(elems,node_id)

        
        
class MeanCurvatureNormal(NormalVectorField):
    """
    Computes the mean curvature normal of a given mesh like
    proposed in [1].
    """
    def computeVectorOnNode(self,node_id):
        """
        We compute the normal in one node,
        with help of the mean curvature formula
        formula.
        """
        Mesh = self.mesh

        elems = Mesh.GetNodeInverseElements(node_id)

        rst_group = self.rst_group
        if rst_group is not None:
            elems = [elem for elem in elems if elem in rst_group]

        from Tools import apply_linear_elements
        elems = apply_linear_elements(Mesh,elems)
        
        return self.meanCurvatureNormalFormula(elems,node_id,voroni=True)

class SalomeNormalField(VectorField):
    """
    This normal vector field reverse engineers the salome extrusion
    method for better results, whithout much effort.
    The vectors are computed once and stored in a dict to have a
    fast mapping.
    """

    def __init__(self,mesh, scalar = 1.0, restricted_group=None,ByAverageNormal=False):

        super(SalomeNormalField,self).__init__(mesh,scalar,restricted_group)
        self.vec_dict = dict([])
        self.eps = 1e-6
        self.prefix = '_N_'
        self.ByAverageNormal = ByAverageNormal
        
        self.updateField(self.restricted_group)
        
    def updateField(self,face_group):
        """
        Updates the dict from a given surface group
        """
        new_mesh = self.prepareMesh(face_group)
        new_mesh, node_grps = self.createNodeGrps(new_mesh)
        new_mesh, top_node_grps = self.preExtrusion(new_mesh,self.ByAverageNormal)
        self._getNodeCoords(new_mesh,node_grps,top_node_grps)
        new_mesh.Clear()
        del new_mesh
        
    def prepareMesh(self,face_group=None):
        """
        Necessary preparation steps.
        """
        if face_group is None: 
            new_mesh = smesh.CopyMesh(self.mesh,'help_mesh',toKeepIDs = True)
        else:
            new_mesh = smesh.CopyMesh(face_group,'help_mesh',toKeepIDs = True)

        return new_mesh

    def createNodeGrps(self,new_mesh):
        """
        Generates the node groups to have a proper mapping between top and below.
        """
        nodes = new_mesh.GetNodesId()
        node_grps = [new_mesh.MakeGroupByIds(self.prefix+str(node), SMESH.NODE, [node]) for node in nodes]
        return new_mesh, node_grps

    def preExtrusion(self,new_mesh,ByAverageNormal=False):
        """
        Executes Salome Extrusion to get the top nodes.
        """
        top_node_grps = new_mesh.ExtrusionByNormal(new_mesh,1.0,1,ByAverageNormal,MakeGroups=True)
        return new_mesh, top_node_grps

    def _computeVector(self,vec1,vec2):
        """
        Computes the vector between 2 points vec1 and vec2.
        """
        vec1 = np.array(vec1,dtype=float)
        vec2 = np.array(vec2,dtype=float)

        vec = vec2-vec1
        norm = np.linalg.norm(vec)
        if norm < self.eps:
            raise ValueError("Error: Vectors are nearly the same!")
        vec /= norm
        return vec
    
    def _getNodeCoords(self,new_mesh,node_grps,top_node_grps):
        """
        Get the node coordinates of the point and his duplicate on the upper surface,
        and computes the vector to store it in the dict.
        """
        top_names = [grp.GetName() for grp in top_node_grps]
        for grp in node_grps:
            node_id = grp.GetIDs()[0]
            vec1 = new_mesh.GetNodeXYZ(node_id)
            index = top_names.index(grp.GetName()+'_top')
            top_grp = top_node_grps[index]
            node_id2 = top_grp.GetIDs()[0]
            vec2 = new_mesh.GetNodeXYZ(node_id2)

            vec = self._computeVector(vec1,vec2)
            self.vec_dict.update([[node_id,vec]])

    def extrudeSurface(self,group=None, edge_groups = [], face_groups = [], create_boundary_elements = True):
        """
        Slight modification for the surface extrusion: After a new surface is generated the field has to
        be updated.
        """
        result = super(SalomeNormalField,self).extrudeSurface(group,edge_groups,face_groups,create_boundary_elements)
        self.updateField(result[0])
        return result
        
    def computeVectorOnNode(self,node_id):
        return self.vec_dict[node_id]
            
            

    
class GroupDependentNormalVectorField(NormalVectorField):
    """
    The multiplication factor of this normal vector field depends on the group
    on which the current node lies.
    """
    def __init__(self,mesh, scalar = 1.0, special_groups = [], special_scalar_factors =[],
                 restricted_group=None):

        if len(special_groups) != len(special_scalar_factors):
            raise ValueError("Error: Number of Groups does not match number of scalar factors!")

        self.special_groups_list = []
        self.special_groups = []
        self.special_scalar_factors = []
        self.updateSpecialGroups(special_groups,special_scalar_factors)

        super(GroupDependentNormalVectorField,self).__init__(mesh,scalar=scalar,
                                                             restricted_group = restricted_group)

    def computeVectorOnNode(self,node_id):
        """
        We compute the normal in one node,
        with help of the mean normal formula.
        """
        
        result = super(GroupDependentNormalVectorField,self).computeVectorOnNode(node_id)
        for i in range(len(self.special_scalar_factors)):
            if node_id in self.special_groups[i]:
                result *= self.special_scalar_factors[i]

        return result

    def updateSpecialGroups(self,new_groups,new_factors):
        """
        As name suggests, adds new groups to the lists of special groups and adds the new factors
        """
        self.special_groups += [set(group.GetNodeIDs()) for group in new_groups]
        self.special_groups_list += new_groups

        self.special_scalar_factors += new_factors
    
    def extrudeSurface(self,group=None, edge_groups = [], face_groups = [], create_boundary_elements = True):
        """
        Slight modification for the surface extrusion: After a new surface is generated the field has to
        be updated.
        """
        result = super(NormalVectorField,self).extrudeSurface(group,edge_groups,face_groups+self.special_groups_list,create_boundary_elements)
        new_special_face_groups = result[4][len(face_groups):]
        self.updateSpecialGroups(new_special_face_groups,self.special_scalar_factors)
        return result
 

class SalomeGroupDependentNormalVectorField(SalomeNormalField):
    """
    The multiplication factor of this normal vector field depends on the group
    on which the current node lies. New Implementation with salome version
    """
    def __init__(self,mesh, scalar = 1.0, special_groups = [], special_scalar_factors =[],
                 restricted_group=None,ByAverageNormal=False):

        if len(special_groups) != len(special_scalar_factors):
            raise ValueError("Error: Number of Groups does not match number of scalar factors!")
            
        self.special_groups_list = []
        self.special_groups = []
        self.special_scalar_factors = []
        self.updateSpecialGroups(special_groups,special_scalar_factors)

        super(SalomeGroupDependentNormalVectorField,self).__init__(mesh,scalar=scalar,
                                                                   restricted_group = restricted_group,
                                                                   ByAverageNormal = ByAverageNormal)

        
    def updateSpecialGroups(self,new_groups,new_factors):
        """
        As name suggests, adds new groups to the lists of special groups and adds the new factors
        """
        self.special_groups += [set(group.GetNodeIDs()) for group in new_groups]
        self.special_groups_list += new_groups

        self.special_scalar_factors += new_factors


    def computeVectorOnNode(self,node_id):
        """
        We compute the normal in one node,
        with help of the mean normal formula.
        """
        
        result = super(SalomeGroupDependentNormalVectorField,self).computeVectorOnNode(node_id)
        for i in range(len(self.special_scalar_factors)):
            if node_id in self.special_groups[i]:
                result *= self.special_scalar_factors[i]

        return result

    def extrudeSurface(self,group=None, edge_groups = [], face_groups = [], create_boundary_elements = True):
        """
        Slight modification for the surface extrusion: After a new surface is generated the field has to
        be updated.
        """
        result = super(SalomeGroupDependentNormalVectorField,self).extrudeSurface(group,edge_groups,
                                                                                  face_groups+self.special_groups_list,
                                                                                  create_boundary_elements)
        new_special_face_groups = result[4][len(face_groups):]
        self.updateSpecialGroups(new_special_face_groups,self.special_scalar_factors)
        return result

        
############################################################################################
class MultiLayerVectorField(VectorField):
    """
    Vectorfields for which are dependent on the nr of layers, and
    the current Layer
    """
    def __init__(self,mesh, scalar = 1.0, restricted_group=None):
        """
        Init method has to set the default value for layers to one.
        """        
        self._setNrLayers(1)
        super(MultiLayerVectorField,self).__init__(mesh, scalar, restricted_group)

        
    def _setNrLayers(self,nr_layers):
        """
        Sets the parameter nr_layers.
        nr_layers is a help parameter,which is not set per default.

        Arguments:
        - `self`:
        - `nr_layers`
        """
        self.nr_layers = nr_layers

    def _setAppliedExtrusion(self,applied_extrusions):
        """
        Sets the number of the applied extrusion layers, 
        where the vector field lives currently on.

        Arguments:
        - `self`:
        - `layer`: nr of current layer
        """
        self._applied_extrusions = applied_extrusions


    def _getLayer(self,layer):
        """
        Sets the number of the current layer, where
        the vector field lives currently on.

        Arguments:
        - `self`:
        - `layer`: nr of current layer
        """
        return self.current_layer

    def _preProcess(k,group,edge_groups,face_groups):
        """
        Possible preprocessing steps for layer extrusion.

        Arguments:
        - `self`: 
        - `k`: number of extrusions, or list of extrusion thicknesses
        - `mesh`: Optional smesh.Mesh instance. Per defualt it is self.mesh
        - `groups`: Optional group of elements on which we apply the vector field.
        - `edge_groups`: Groups of edges.
        - `face_groups`: Groups of faces.
        """
        pass

    def _postProcess(k,group,edge_groups,face_groups):
        """
        Possible preprocessing steps for layer extrusion.

        Arguments:
        - `self`: 
        - `k`: number of extrusions, or list of extrusion thicknesses
        - `mesh`: Optional smesh.Mesh instance. Per defualt it is self.mesh
        - `groups`: Optional group of elements on which we apply the vector field.
        - `edge_groups`: Groups of edges.
        - `face_groups`: Groups of faces.
        """
        pass

    
    def extrudeSurfaceTimes(self,k,group=None, edge_groups = [], face_groups = []):
        """
        This method applies the vector field on a surface and creates
        a translated one over k steps.

        Arguments:
        - `self`: 
        - `k`: number of extrusions, or list of extrusion thicknesses
        - `mesh`: Optional smesh.Mesh instance. Per defualt it is self.mesh
        - `groups`: Optional group of elements on which we apply the vector field.
        - `edge_groups`: Groups of edges.
        - `face_groups`: Groups of faces.
        """
        raise NotImplementedError("Error: Multiple extrusion not tested for any subclass yet!")
        self._setNrLayers(k)
        self._setAppliedExtrusion(0)
        # preperations steps
        self._preProcess(k,group,edge_groups,face_groups)
        
        result = super(MultiLayerVectorField, self).extrudeSurfaceTimes(k,group=group, 
                                                            edge_groups = edge_groups, 
                                                            face_groups = face_groups)
        # set back to normal behaviour
        self._setNrLayers(1)
        # make postprocessing steps
        self._postProcess(k,group,edge_groups,face_groups)
        return result


class PlaneProjectionVectorField(MultiLayerVectorField):
    """
    A layer dependent VectorField which projects a
    mesh surface onto a plane which is given by a
    local coordinate system Q = array([u,v,w]) which is an
    orthonormal matrix, where u,v form the local
    xy-plane, and w is the local z-Axis, 
    and an origin O.
    Further it will be assumed, that the plane
    lies completely on one side of the surface,
    and that there is a minimal distance d between
    surface and plane. 

    To optimize computation time, vectors of given groups will be calculated beforehand
    and then looked up, instead computed everytime.
    """
    def __init__(self,mesh,O,Q,d,signum = None, scalar = 1.0, restricted_group=None):
        """
        Arguments:
        - `self`: 
        - `O`: origin given as an numpy array or an tuple/list of 3 real numbers.
        - `Q`: An orthonormal matrix provided either as numpy array or as list of vectors.
        - `d`: Real number which represents the minimal distance between the currrent surface and the plane.
        - `signum`: Sign which represents the side where the current surface should lie. If signum is None the side of the plane is not considered.
        """

        
        self.O = array(O)
        self.Q = array(Q)
        if d < 0:
            raise ValueError("Error: Distance measure is negativ!")
        self.d = d
        if not signum is None:
            if not (signum is 1) or not (signum is -1):
                raise ValueError(u"Error: Signum is not ±1!")
        self.signum = signum
        self._current_table= {}
        
        from MyMath.Types import GeometricTransformation 
        self.inv_trafo = GeometricTransformation(self.Q,self.O)
        self.trafo = self.inv_trafo.inv()

        self._vectors = None
        super(PlaneProjectionVectorField,self).__init__(mesh, scalar = scalar, 
                                                restricted_group=restricted_group)


    def _preProcess(k,group,edge_groups,face_groups):
        """
        Possible preprocessing steps for layer extrusion.
        Computes the vectors beforehand.

        Arguments:
        - `self`: 
        - `k`: number of extrusions, or list of extrusion thicknesses
        - `mesh`: Optional smesh.Mesh instance. Per defualt it is self.mesh
        - `groups`: Optional group of elements on which we apply the vector field.
        - `edge_groups`: Groups of edges.
        - `face_groups`: Groups of faces.
        """
        self._vectors = self.computeProjections(group)

    def _postProcess(k,group,edge_groups,face_groups):
        """
        Possible preprocessing steps for layer extrusion.
        Resets _vectors to None

        Arguments:
        - `self`: 
        - `k`: number of extrusions, or list of extrusion thicknesses
        - `mesh`: Optional smesh.Mesh instance. Per defualt it is self.mesh
        - `groups`: Optional group of elements on which we apply the vector field.
        - `edge_groups`: Groups of edges.
        - `face_groups`: Groups of faces.
        """
        self._vectors = None

        
    def _processLookupTables(self,lookup_tables):
        """
        Takes the latest lookup tables to find the original node.
        - `self`: 
        - `lookup_table`: dict to process. 
        """
        latest_table = lookup_tables[self._applied_extrusions]
        keys = self._current_table.keys()
        new_table = [[latest_table[key],self._current_table[key]] for key in keys]
        self._current_table = dict(new_table)

    def getNodeVectors(self,group=None):
        """
        Collects all the vectors from the nodes for transformation
        """

        if group is None:
            node_ids = self.mesh.GetNodesId()
        else:
            node_ids = group.GetNodeIDs()

        self._internal_ids = dict([[node_ids[i],i] for i in range(len(node_ids))])
        vectors = array([self.mesh.GetNodeXYZ(node_id) for node_id in node_ids])
        return vectors.transpose()

    def _makeChecks(self,trafo_vecs):
        """
        Checks if the distance criterias aren't violated.
        Takes the latest lookup tables to find the original node.
        - `self`: 
        - `trafo_vecs`: transformed vectors to check. 
        """
        # check if we have single vector
        if len(trafo_vecs.shape) is 1:
            to_check = array([trafo_vecs[-1]])
        else:
            to_check = trafo_vecs[-1,:]
            
        if not self.signum is None: 
            if self.signunm > 0:
                if not all(to_check >= self.d):
                    raise ValueError("Error: Surface does not match criterias!")

                if self.signum < 0:
                    if not all(to_check <= -self.d):
                        raise ValueError("Error: Surface does not match criterias!")
        else:
            if not all(abs(to_check) >= self.d):
                raise ValueError("Error: Surface does not match criterias!")
            
    def computeProjections(self,group=None):
        """
        Computes the vectors which have to be projeceted.
        """
        vectors = self.getNodeVectors(group)
        traf_vecs = copy(self.trafo(vectors))
        if not (self.signum is None and self.d == 0):
            self._makeChecks(traf_vecs)
        traf_vecs[-1,:] = 0.0
        proj_vecs = self.inv_trafo(traf_vecs)
        
        return  proj_vecs

    def computeSingleProjection(self,node_id):
        """
        Computes the Projection for one node

        Arguments:
        - `self`: 
        - `node_id`: id of the current node.
        """
        vector = array(self.mesh.GetNodeXYZ(node_id))
        trafo_vec = self.trafo(vector)
        if not (self.signum is None and self.d == 0):
            self._makeChecks(trafo_vec)

        trafo_vec[-1] = 0.0
        proj_vec = self.inv_trafo(trafo_vec)

        return proj_vec.reshape(DIMENSION)
        
    def distribution(self):
        """
        Distribution function for dividing edges.
        """
        return 1.0/self.nr_layers

    def moveSurface(self,group=None):
        """
        This method applies the vector field on a surface and creates
        a translated one. Optional the new surface can be stored in a new mesh.
        For The PlaneProjectionVectorField it is reasonable to preprocess 
        the computation.

        Arguments:
        - `self`: 
        - `mesh`: Optional smesh.Mesh instance. Per defualt it is self.mesh
        - `group`: Optional group of elements on which we apply the vector field.
        """
        self._vectors = self.computeProjections(group)
        super(PlaneProjectionVectorField,self).moveSurface(group)
        # reset vectors
        self._vectors = None
    
    def computeVectorOnNode(self,node_id):
        """
        The vectors are computed only at the first time,
        then lookup tables will be used to transport
        the nodes up to the plane.

        The formula is simply y = x + v/k
        where x is the current node, y is the
        new node and v is the direction.

        Every produced edge is distributed by a formula.

        Arguments:
        - `self`: 
        - `node_id`: id of the current node.
        """
        node_vec = array(self.mesh.GetNodeXYZ(node_id))
        if self._vectors is None:
            vector = (self.computeSingleProjection(node_id)*self.distribution())

        else:
            if self._current_table:
                original_id = self._current_table[node_id]
                internal_id = self._internal_ids[original_id]
            else:
                internal_id = self._internal_ids[node_id]
                
            vector = self._vectors[:,internal_id]*self.distribution()

        vector = vector.reshape(DIMENSION) #flatten vector
        return vector - node_vec
        
class FaceProjectVectorField(MultiLayerVectorField):
    """
    A Multilayer vector field which holds a face 
    provided as geometric object, and uses geompy 
    to compute the projection, similar to the
    PlaneProjectionVectorField.
    """
    def __init__(self,mesh,face,d,scalar = 1.0, restricted_group=None):
        """
        Arguments:
        - `self`: 
        - `face`: geometric object which is a face
        - `d`: Real number which represents the minimal distance between the currrent surface and the plane.
        """
        from MyGeom.Types import MyGeomObject, MyFace
        if isinstance(face,MyGeomObject) and isinstance(face,MyFace):
            self.face = face.getGeomObject()
        elif isinstance(face,GEOM._objref_GEOM_Object):
            if face.GetShapeType() == GEOM.FACE:
                self.face = face
            else:
                raise ValueError("Error: Geometric Object is not a face!")
        else:
            raise ValueError("Error: Geometric Object is not a face!")


        if d < 0:
            raise ValueError("Error: Distance measure is negativ!")
        
        self.d = d
        super(FaceProjectVectorField,self).__init__(mesh, scalar = scalar, 
                                                    restricted_group=restricted_group)

    def computeVectorOnNode(self,node_id):
        """
        The current node is projected on the face.
        The vector of the direction is returned.
        
        Arguments:
        - `self`: 
        - `node_id`: id of the current node.
        """

        coords = self.mesh.GetNodeXYZ(node_id)
        vertex = geompy.MakeVertex(*coords)
        projected_vertex = geompy.MakeProjection(vertex,self.face)
        new_coords = geompy.PointCoordinates(projected_vertex)
        
        vec =  array(new_coords) - array(coords)
        self._makeChecks(vec)
        return vec
        
    def distribution(self):
        """
        Distribution function for dividing edges.
        """
        return 1.0/self.nr_layers

    def _makeChecks(self,trafo_vec):
        """
        Checks if the distance criterias aren't violated.
        Takes the latest lookup tables to find the original node.
        - `self`: 
        - `trafo_vecs`: transformed vector to check. 
        """
        if norm(trafo_vec) < self.d:
            raise ValueError("Error: Surface does not match criterias!")
