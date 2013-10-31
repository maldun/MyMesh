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

from numpy import array, ndarray, arange, cross, zeros
from numpy.linalg import norm
from numpy import float64 as data_type

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
    def getNodes(self):
        return self.mesh.GetElemNodes(self.getIdNr())

    def computeNormal(self,store = True):
        raise NotImplementedError("Error: Not implemented!")

    def getNormal(self):
        return self.normal
        
    def computeArea(self,store = True):
        raise NotImplementedError("Error: Not implemented!")

    def getArea(self):
        return self.area


class Tria3(FaceElement):
    """
    Simple linear triangles with 3 corners
    """
    
    def computeNormalOp(self):
        """
        Computes the normal of the triangle.
        This function only executes the computation
        """

        nodes = self.getNodes()
        p0 = self.mesh.GetNodeXYZ(nodes[0])
        p1 = self.mesh.GetNodeXYZ(nodes[1])
        p2 = self.mesh.GetNodeXYZ(nodes[-1])
        
        vec1 = array(p1) - array(p0)
        vec2 = array(p2) - array(p0)
        vec_normal = cross(vec1, vec2)
        return vec_normal/norm(vec_normal), norm(vec_normal)/2.0

    def computeNormal(self,store = True):
        """
        It stores the vector in the object per default.
        """
        if store:
            self.normal = self.computeNormalOp()[0]
        else:
            return self.computeNormalOp()[0]

    def computeArea(self,store = True):
        """
        It stores the vector in the object per default.
        """
        if store:
            self.area = self.computeNormalOp()[1]
        else:
            return self.computeNormalOp()[1]

        

class Quad4(FaceElement):
    """
    Simple linear quadrangles with 4 corners
    """
    
    def computeNormalOp(self):
        """
        Computes the normal of the triangle.
        This function only executes the computation
        """

        nodes = self.mesh.GetElemNodes(self.getIdNr())
        p0 = self.mesh.GetNodeXYZ(nodes[0])
        p1 = self.mesh.GetNodeXYZ(nodes[1])
        p2 = self.mesh.GetNodeXYZ(nodes[-1])
        
        vec1 = array(p1) - array(p0)
        vec2 = array(p2) - array(p0)
        vec_normal = cross(vec1, vec2)
        return vec_normal/norm(vec_normal), norm(vec_normal)

    def computeNormal(self,store = True):
        """
        It stores the vector in the object per default.
        """
        if store:
            self.normal = self.computeNormalOp()[0]
        else:
            return self.computeNormalOp()[0]

    def computeArea(self,store = True):
        """
        It stores the vector in the object per default.
        """
        if store:
            self.area = self.computeNormalOp()[1]
        else:
            return self.computeNormalOp()[1]



class VectorField(object):
    """
    Class that represents a vector field which is applied on the nodes of 
    mesh. Thus is holds also the information about the mesh, on which it works.
    """

    def __init__(self,mesh, scalar = 1.0):
        
        self.mesh = mesh
        self.scalar = scalar

    def getVectorOnNode(self,node_id):
        return self.computeVectorOnNode(node_id)*self.getScalar()

    def computeVectorOnNode(self,node_id):
        raise NotImplementedError('Error: This method is a stub')

    def getScalar(self): return self.scalar

    def setScalar(self,scalar): self.scalar = scalar
    
    def scalarMultiplication(self,scalar):
        """
        Implements multiplication with scalar.
        """
        mult_scalar = self.getScalar()*scalar
        self.setScalar(mult_scalar)

    def __rmul__(self,scalar): 
        
        self.scalarMultiplication(scalar)
        return self

    def applyVectorOnNode(self,node_id, mesh = None):
        """
        Apply the vector field on a given node. This method creates the
        displaced node in the given mesh, and returns the index of the new node.
        The node is created in the current mesh per defualt, but it is possible 
        to add it to a different mesh.

        Arguments:
        - `self`: 
        - `node_id`: integer with the Id of the node we want to apply the vector.
        """
        
        if mesh is None:
            mesh = self.mesh 

        node_vec = self.mesh.GetNodeXYZ(node_id)
        translated_node_vec = node_vec + self.getVectorOnNode(node_id)
        new_node = mesh.AddNode(*translated_node_vec.tolist())
        
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

        # update table
        table_out = table.copy() # copy table to avoid overiding in memory
        table_out.update([[elem_nodes[i],new_nodes[i]] for i in range(len(elem_nodes))])

        return new_face, table_out
        

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

            face_id, new_lookup = self.applyVectorFieldOnFace(face,mesh,new_nodes = True)
            new_face_ids += [face_id]
            lookup_table
            

class NormalVectorField(VectorField):

    """
    Class to compute the normal vector field on a mesh.
    The current version only works on surfaces, and on
    linear elements. Hopefully this will change.
    """
    def __init__(self,mesh):
        
        super(NormalVectorField,self).__init__(mesh)

        # filter linear and triangle elements
        filter_linear_tri = GetFilter(FACE, FT_LinearOrQuadratic, Geom_TRIANGLE)
        filter_linear_quad = GetFilter(FACE, FT_LinearOrQuadratic, Geom_QUADRANGLE)

        ids_tri = mesh.GetIdsFromFilter(filter_linear_tri)
        ids_quad = mesh.GetIdsFromFilter(filter_linear_quad)

        self.tria3 = [Tria3(mesh,id_tri) for id_tri in ids_tri]
        self.quad4 = [Tria3(mesh,id_tri) for id_tri in ids_tri]

    def meanNormalFormula(self,elems):
        """
        Compute the mean normal, whith the formulas

        math:: 

            m(x) = (\sum_{\Delta \in Faces(x)} n_{\Delta}(x))/\#Faces(x)
            n(x) = m(x)/||m(x)||
        """

        result = zeros(3)

        for elem in elems:
            result += elem.computeNormal(store=False)

        result /= len(elems)
        return result/norm(result)

    def computeVectorOnNode(self,node_id):
        """
        We compute the normal in one node,
        (currently) with help of the mean value
        formula.
        """
        Mesh = self.mesh

        elems = Mesh.GetNodeInverseElements(node_id)

        from Tools import apply_linear_elements
        elems = apply_linear_elements(Mesh,elems)
        
        return self.meanNormalFormula(elems)

        
        
        
        
        
        
    
