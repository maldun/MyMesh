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
    def __init__(self,mesh,id_nr):
        super(FaceElement,self).__init__(mesh,id_nr)
        self.normals = {}

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

    def getCurvatureVector(self,node):
        raise NotImplementedError("Error: Not implemented!")


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

        node_vec = array(self.mesh.GetNodeXYZ(node_id))
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

        return new_face, new_nodes
        

    def computeSurfaceExtrusion(self,group=None, edge_groups = []):
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
            
            return new_face_ids, new_vol_ids, new_edge_groups, lookup_table
        
        #else:
        
        return new_face_ids, new_vol_ids, lookup_table
    
    def extrudeSurface(self,group=None, edge_groups = [], create_boundary_elements = True):
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

        if edge_groups:
            new_face_ids, new_vol_ids, new_edge_groups, lookup_table = self.computeSurfaceExtrusion(group=group, edge_groups = edge_groups)
        
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

            salome.sg.updateObjBrowser(0)
            return face_group, vol_group, new_salome_edge_groups, new_salome_edge_face_groups, bnd_faces, lookup_table 

        else:
            salome.sg.updateObjBrowser(0)
            return face_group, vol_group, [], [], bnd_faces,lookup_table


    def extrudeSurfaceTimes(self,k,group=None, edge_groups = []):
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
        
        face_groups = [None]*k
        vol_groups = [None]*k
        extruded_edge_groups = [None]*k
        extruded_edge_surface_groups = [None]*k
        bnd_face_groups = [None]*k
        lookup_tables = [None]*k
        
        face_groups[-1] = group
        extruded_edge_groups[-1] = edge_groups

        original_rst_group = self.getRestricedGroup()

        for i in range(k):
            if thicknesses:
                self.setScalar(thicknesses[i])

            face_groups[i], vol_groups[i], extruded_edge_groups[i], extruded_edge_surface_groups[i], bnd_face_groups[i], lookup_tables[i] = self.extrudeSurface(face_groups[i-1],extruded_edge_groups[i-1])
            self.setRestricedGroup(face_groups[i])
            
                
        self.setRestricedGroup(original_rst_group)
        return face_groups, vol_groups, extruded_edge_groups, extruded_edge_surface_groups, bnd_face_groups, lookup_tables
            

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


class NormalVectorField(VectorField):

    """
    Class to compute the normal vector field on a mesh.
    The current version only works on surfaces, and on
    linear elements. Hopefully this will change.
    """
    def __init__(self,mesh,scalar = 1.0, restricted_group=None):
        
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
            

        if restricted_group is not None:
            self.rst_group = set((self.getRestricedGroup()).GetIDs())
        else:
            self.rst_group = None

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

        result = zeros(3)

        for elem in elems:
            result += elem.computeNormal(node_id,store=False)

        result /= len(elems)
        return result/norm(result)


    def meanNormalCurvatureFormula(self,elems,node_id):
        """
        Compute an averaged mean curvature normal
        due to the formula

        math::

            \frac{1}{2} \sum_{j \in N_1(i)} (\cot \alpha_{ij} + \cot \beta_{ij})(x_j - x_i)

        where x_j are the neighbouring nodes of x_i, and alpha_{ij}, beta_{ij} are the angles
        oposite to the edge (x_i,x_j) supposing we have a triangle mesh. If quad elements appear
        we cut the quadrngles in a half.
        (see  'Discrete Differential-Geomety Operators for Triangulated 2-Manifolds'
         by Meyer, Desbrun, Schroeder and Barr for more details) 
        """

        normal = zeros(3)
        for elem in elems:
            normal += elem.getCurvatureVector(node_id)

        return normal
    
    def computeVectorOnNode(self,node_id):
        """
        We compute the normal in one node,
        (currently) with help of the mean value
        formula.
        """
        Mesh = self.mesh

        elems = Mesh.GetNodeInverseElements(node_id)

        rst_group = self.rst_group
        if rst_group is not None:
            elems = [elem for elem in elems if elem in rst_group]

        from Tools import apply_linear_elements
        elems = apply_linear_elements(Mesh,elems)
        
        return self.meanNormalFormula(elems,node_id)

        
        
class MeanCurvatureNormal(NormalVectorField):
    """
    Computes the mean curvature normal of a given mesh like
    proposed in the paper 'Discrete Differential-Geomety Operators for Triangulated 2-Manifolds'
    by Meyer, Desbrun, Schroeder and Barr. 
    """
    pass
