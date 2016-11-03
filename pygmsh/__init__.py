""" Module containing an expanded python gmsh class"""
from __future__ import print_function

import vtk
from vtk.util import numpy_support
import numpy
import struct
import copy

class LineDict(object):
    def __init__(self,*args,**kwargs):
        self._data=dict(*args,**kwargs)
    def __getitem__(self,i):
        return self._data[tuple(sorted(i))]
    def __len__(self):
        return len(self._data)
    def setdefault(self,k,v):
        if tuple(k) in self._data:
            return self._data[tuple(k)]
        ks = list(k)
        ks.reverse()
        if tuple(ks) in self._data:
            return -self._data[tuple(ks)]
        ## Otherwise add it
        ks=tuple(k)
        self._data[ks]=v
        return self._data[ks]
    def items(self):
        return self._data.items()

class Node(object):
    def __init__(self, x=None):
        """Initialise node from data if available"""
        self.vertices = numpy.empty(3)
        if x:
            self.vertices[:len(x)] = x[:]

    def __getitem__(self,i):
        return self.vertices[i]

    def __len__(self):
        return len(self.vertices)

    def __hash__(self):
        return hash(self.vertices)

    def __eq__(self, x):
        return all(self.vertices == numpy.array(x))

class Element(object):

    def __init__(self, etype, tags=(0,0), nodes=()):
        self.etype = etype
        self.tags = tags
        self.nodes = nodes

        self._data=[etype, tags, nodes]

    def __len__(self):
        return len(self.nodes)

    def __eq__(self,x):
        return self.nodes == x

    def __iter__(self):
        return iter(self._data)

    def __getitem__(self,i):
        return self._data[i]
        

class GmshMesh(object):
    """This is a class for storing nodes and elements. Based on Gmsh.py

    Members:
    nodes -- A dict of the form { nodeID: [ xcoord, ycoord, zcoord] }
    elements -- A dict of the form { elemID: (type, [tags], [nodeIDs]) }

    Methods:
    read([file]) -- Parse a Gmsh version 1.0 or 2.0 mesh file
    write([file]) -- Output a Gmsh version 2.0 mesh file
    """

    def __init__(self, filename=None):
        """Initialise Gmsh data structure"""
        self.nodes = {}
        self.elements = {}
        self.physical_points = set()
        self.physical_lines = set()
        self.physical_surfaces = set()
        self.physical_volumes = set()
        self.filename = filename
        if self.filename:
            self.read()

    def __add__(self, mesh):
        """ Get union of meshes through concatenation."""
        if type(mesh) is not GmshMesh:
            raise TypeError
        
        out = copy.deepcopy(self)
        n_nodes = len(out.nodes)
        n_ele = len(out.elements)

        for id, node in mesh.nodes.items():
            out.nodes[n_nodes+id]=Node(node)

        for id, ele in mesh.elements.items():
            out.elements[n_ele+id]=ele

        return out

    def __iadd__(self,mesh):
        """ Get in place union of meshes through concatenation."""
        if type(mesh) is not GmshMesh:
            raise TypeError

        n_nodes = len(self.nodes)
        n_ele = len(self.elements)

        for id, node in mesh.nodes.items():
            self.nodes[n_nodes+id]=Node(node)

        for id, ele in mesh.elements.items():
            self.elements[n_ele+id]=Element(*ele)

        return self

    def insert_node(self, node_id, pos=None):
        """Insert node into mesh data structure."""
        if node_id in self.nodes.keys():
            raise KeyError
        self.nodes[node_id] = Node(pos)

    def insert_element(self, element_id, etype, tags, nodes):
        """Insert element into mesh data structure."""
        if element_id in self.elements.keys():
            raise KeyError
        self.elements[element_id] = Element(etype, tags, nodes)

        self.__add_physical_id__(etype,tags[-1])

    def nodecount(self):
        """Return number of nodes in the mesh."""
        return len(self.nodes)

    def elementcount(self):
        """Return number of elements in the mesh."""
        return len(self.elements)

    def transform(self, func):
        """Apply a transformation to the mesh nodes."""

        for k, v in self.nodes.items():
            self.nodes[k] = Node(func(v))


    def reset(self):
        """Reinitialise Gmsh data structure"""
        self.nodes = {}
        self.elements = {}

    def read(self, mshfile=None):
        """Read a Gmsh .msh file.
        
        Reads Gmsh format 1.0 and 2.0 mesh files, storing the nodes and
        elements in the appropriate dicts.
        """

        if not mshfile:
            mshfile = open(self.filename,'r')

        readmode = 0
        print('Reading %s'%mshfile.name)
        line='a'
        while line:
            line=mshfile.readline()
            line = line.strip()
            if line.startswith('$'):
                if line == '$NOD' or line == '$Nodes':
                    readmode = 1
                elif line == '$ELM':
                    readmode = 2
                elif line == '$Elements':
                    readmode = 3
                elif line == '$MeshFormat':
                    readmode = 4
                else:
                    readmode = 0
            elif readmode:
                columns = line.split()
                if readmode == 4:
                    if len(columns)==3:
                        vno,ftype,dsize=(float(columns[0]),
                                         int(columns[1]),
                                         int(columns[2]))
                        print(('ASCII','Binary')[ftype]+' format')
                    else:
                        endian=struct.unpack('i',columns[0])
                if readmode == 1:
                    # Version 1.0 or 2.0 Nodes
                    try:
                        if ftype==0 and len(columns)==4:
                            self.nodes[int(columns[0])] = Node(map(float, columns[1:]))
                        elif ftype==1:
                            nnods=int(columns[0])
                            for N in range(nnods):
                                data=mshfile.read(4+3*dsize)
                                i,x,y,z=struct.unpack('=i3d',data)
                                self.nodes[i]=Node((x,y,z))
                            mshfile.read(1)
                    except ValueError:
                        print('Node format error: '+line, ERROR)
                        readmode = 0
                elif ftype==0 and  readmode > 1 and len(columns) > 5:
                    # Version 1.0 or 2.0 Elements 
                    try:
                        columns = map(int, columns)
                    except ValueError:
                        print('Element format error: '+line, ERROR)
                        readmode = 0
                    else:
                        (id, type) = columns[0:2]
                        if readmode == 2:
                            # Version 1.0 Elements
                            tags = columns[2:4]
                            nodes = columns[5:]
                        else:
                            # Version 2.0 Elements
                            ntags = columns[2]
                            tags = columns[3:3+ntags]
                            nodes = columns[3+ntags:]
                        self.elements[id] = Element(type, tags, nodes)
                        self.__add_physical_id__(type,tags[-1])
                elif readmode == 3 and ftype==1:
                    tdict={1:2,2:3,3:4,4:4,5:5,6:6,7:5,8:3,9:6,10:9,11:10}
                    try:
                        neles=int(columns[0])
                        k=0
                        while k<neles:
                            etype,ntype,ntags=struct.unpack('=3i',mshfile.read(3*4))
                            k+=ntype
                            for j in range(ntype):
                                mysize=1+ntags+tdict[etype]
                                data=struct.unpack('=%di'%mysize,
                                                   mshfile.read(4*mysize))
                                self.elements[data[0]]=Element(etype,
                                                        data[1:1+ntags],
                                                        data[1+ntags:])
                                self.__add_physical_id__(etype,data[ntags])
                    except:
                        raise
                    mshfile.read(1)
                            
        print('  %d Nodes'%len(self.nodes))
        print('  %d Elements'%len(self.elements))

        mshfile.close()

    def __add_physical_id__(self,etype,tag):
        """Insert new physical id into data structure."""

        physical_dict = {15:self.physical_points,
                         1:self.physical_lines,
                         2:self.physical_surfaces,
                         4:self.physical_volumes}

        physical_dict[etype].add(tag)
        

    def write_ascii(self, filename=None):
        """Dump the mesh out to a Gmsh 2.0 msh file."""

        if not filename:
            filename = self.filename

        mshfile = open(filename, 'w')

        print('$MeshFormat\n2.0 0 8\n$EndMeshFormat', file=mshfile)
        print('$Nodes\n%d'%len(self.nodes), file=mshfile)
        for node_id, coord in self.nodes.items():
            print(node_id,' ',' '.join([str(c) for c in  coord]), sep="", file=mshfile)
        print('$EndNodes',file=mshfile)
        print('$Elements\n%d'%len(self.elements),file=mshfile)
        for ele_id, elem in self.elements.items():
            (ele_type, tags, nodes) = tuple(elem)
            print(ele_id,' ',ele_type,' ',len(tags),' ',
                  ' '.join([str(c) for c in tags]),' ',
                  ' '.join([str(c) for c in nodes]), sep="", file=mshfile)
        print('$EndElements',file=mshfile)

    def write_binary(self, filename=None):
        """Dump the mesh out to a Gmsh 2.0 msh file."""

        if not filename:
            filename = self.filename

        mshfile = open(filename, 'wr')

        mshfile.write("$MeshFormat\n2.2 1 8\n")
        mshfile.write(struct.pack('@i',1))
        mshfile.write("\n$EndMeshFormat\n")
        mshfile.write("$Nodes\n%d\n"%(len(self.nodes)))
        for node_id, coord in self.nodes.items():
            mshfile.write(struct.pack('@i',node_id))
            mshfile.write(struct.pack('@3d',*coord))
        mshfile.write("\n$EndNodes\n")
        mshfile.write("$Elements\n%d\n"%(len(self.elements)))
        for ele_id, elem in self.elements.items():
            (ele_type, tags, nodes) = tuple(elem)
            mshfile.write(struct.pack('@i',ele_type))
            mshfile.write(struct.pack('@i',1))
            mshfile.write(struct.pack('@i',len(tags)))
            mshfile.write(struct.pack('@i',ele_id))
            for c in tags:
                mshfile.write(struct.pack('@i',c))
            for c in nodes:
                mshfile.write(struct.pack('@i',c))
        mshfile.write("\n$EndElements\n")
                      
        mshfile.close()

    def write_geo(self, filename = None, use_ids=True, use_physicals=True,
                  compound_surfaces=[]):
        """Dump the mesh out to a Gmsh .geo geometry file."""

        string=""

        if compound_surfaces:
            string += 'Mesh.RemeshAlgorithm=1;\n'

        for k,x in self.nodes.items():
            string += 'Point(%d) = {%f,%f,%f};\n'%(k,x[0],x[1],x[2])

        line_no = 1;

        lines = LineDict()
        surfaces = []
        line_ids = {}
        surface_ids = {}

        for k,l in self.elements.items():
            if l[0]==1:
                if use_physicals:
                    line_ids.setdefault(l[1][0],[]).append(line_no)
                else:
                    line_ids.setdefault(l[1][1],[]).append(line_no)
                lines.setdefault(l[2][:],len(lines)+1)
            elif l[0]==2:
                if use_physicals:
                    surface_ids.setdefault(l[1][0],[]).append(k)
                else:
                    surface_ids.setdefault(l[1][1],[]).append(k)
                e1=lines.setdefault([l[2][0],l[2][1]],len(lines)+1)
                e2=lines.setdefault([l[2][1],l[2][2]],len(lines)+1)
                e3=lines.setdefault([l[2][2],l[2][0]],len(lines)+1)
                surfaces.append([e1,e2,e3])

        for line, line_id in sorted(lines.items(),key=lambda x:x[1]):
            string += "Line(%d) = {%d,%d};\n"%(line_id,line[0],line[1])

        for surface_id, surface in enumerate(surfaces):
            string += "Line Loop(%d) = {%d,%d,%d};"%tuple([surface_id+1]+surface)
            string += "Plane Surface(%d) = %d;\n"%(surface_id+1,surface_id+1)

        if use_ids:
            for k, v in line_ids.items():
                string += "Physical Line(%d) = {%s};\n"%(k,",".join(map(str,v)))
            for k, v in surface_ids.items():
                if k in compound_surfaces:
                    string += "Compound Surface(%d) = {%s};\n"%(len(surfaces)+k+1,",".join(map(str,v)))
                    string += "Physical Surface(%d) = {%s};\n"%(k,len(surfaces)+k+1)
                else:
                    string += "Physical Surface(%d) = {%s};\n"%(k,",".join(map(str,v)))

        geofile = open(filename, 'w')
        geofile.write(string)
        geofile.close()

    def write_simple_geo(self, filename = None, use_ids=True, use_physicals=True):
        """Dump the mesh out to a ,geo file after joining coplanar surfaces."""

        
        string = ""

        lines = LineDict()
        surfaces = []
        edges = {}
        normals = []
        line_ids = {}
        vertices = set()
        surface_ids = {}
        open_surf ={}


        for k,l in self.elements.items():
            if l[0]==1:
                if use_physicals:
                    line_ids.setdefault(l[1][-1],[]).append(line_no)
                else:
                    line_ids.setdefault(l[1][1],[]).append(line_no)
                line_ids.setdefault(l[1][-1],[]).append(line_no)
                lines.setdefault(l[2][:],len(lines)+1)
            elif l[0]==2:
                if use_physicals:
                    surface_ids.setdefault(l[1][-1],[]).append(k)
                else:
                    surface_ids.setdefault(l[1][1],[]).append(k)
                e1=lines.setdefault([l[2][0],l[2][1]],len(lines)+1)
                e2=lines.setdefault([l[2][1],l[2][2]],len(lines)+1)
                e3=lines.setdefault([l[2][2],l[2][0]],len(lines)+1)
                surfaces.append([e1,e2,e3])

                sid=len(normals)
                edges.setdefault(abs(e1),[]).append(sid)
                edges.setdefault(abs(e2),[]).append(sid)
                edges.setdefault(abs(e3),[]).append(sid)

                n =numpy.cross(self.nodes[l[2][1]].vertices-self.nodes[l[2][0]].vertices,
                               self.nodes[l[2][2]].vertices-self.nodes[l[2][0]].vertices)
                n=n/numpy.sqrt(sum(n**2))

                normals.append(n)

        element_map = range(len(normals))
        rebuild = set()

        for edge, eles in sorted(edges.items()):

            ## special case of a plane touching itself
            if len(eles)<2:
                continue

            oeles=None

            
            while oeles!=eles and len(eles)>1:
                oeles=eles
                eles=[element_map[eles[0]],element_map[eles[1]]]

            eles = sorted(eles)

            if 1.0-abs(numpy.dot(normals[eles[0]],normals[eles[1]]))<1.0e-8:
                if edge in surfaces[eles[0]]:
                    I = surfaces[eles[0]].index(edge), 1
                else:
                    I = surfaces[eles[0]].index(-edge), -1
                if eles[1]!=eles[0]:
                    if edge in surfaces[eles[1]]:
                        J = surfaces[eles[1]].index(edge), 1
                    else:
                        J = surfaces[eles[1]].index(-edge), -1
                    element_map[eles[1]]=eles[0]
                    surfaces[eles[0]] = surfaces[eles[0]][I[0]+1:]+surfaces[eles[0]][0:I[0]]
                    surfaces[eles[1]] = surfaces[eles[1]][J[0]+1:]+surfaces[eles[1]][:J[0]]
                    if I[1] == J[1]:
                        surfaces[eles[1]].reverse()
                    surfaces[eles[0]].extend(surfaces[eles[1]])
                    surfaces[eles[1]] = None
                    rebuild.discard(eles[1])
                    edges.pop(edge)
                else:
                    while edge in surfaces[eles[0]]:
                            surfaces[eles[0]].remove(edge)
                    while -edge in surfaces[eles[0]]:
                            surfaces[eles[0]].remove(-edge)
                    edges.pop(edge)
                    open_surf.setdefault(eles[0],[]).append(edge)
                    rebuild.add(eles[0])

        rev_lin={}
        for line, line_id in sorted(lines.items(),key=lambda x:x[1]):
            if line_id in edges:
                vertices.add(line[0])
                vertices.add(line[1])
                rev_lin[line_id]=line[1]
                rev_lin[-line_id]=line[0]

        for k,x in self.nodes.items():
            if k in vertices:
                string += 'Point(%d) = {%f,%f,%f};\n'%(k,x[0],x[1],x[2])


        for line, line_id in sorted(lines.items(),key=lambda x:x[1]):
            if line_id in edges:
                string += "Line(%d) = {%d,%d};\n"%(line_id,line[0],line[1])

        for surface_id, surface in enumerate(surfaces):
            if element_map[surface_id]==surface_id:
                if surface_id in open_surf and False:
                    
                    cuts =[]
                    p=surface[-1]
                    for I, l in enumerate(surface):
                        if rev_lin[-l] != rev_lin[p]:
                            cuts.append(I)
                        p = l

                    print(cuts)
                
                    string += "Line Loop(%d) = {%s};"%(surface_id+1,
                                                   ",".join(map(str,surface[:I+1])))
                    if I<len(surface):
                        string += "Line Loop(%d) = {%s};"%(len(surfaces)+surface_id+1,
                                                   ",".join(map(str,surface[I:])))
                        string += "Plane Surface(%d) ={%s,%s};\n"%(surface_id+1,
                                                                  len(surfaces)+surface_id+1, surface_id+1)
                    else:
                        string += "Plane Surface(%d) = %s;\n"%(surface_id+1,surface_id+1)
                else:
                    string += "Line Loop(%d) = {%s};"%(surface_id+1,
                                                   ",".join(map(str,surface)))
                    string += "Plane Surface(%d) = %s;\n"%(surface_id+1,surface_id+1)

        if use_ids:
            for k, v in line_ids.items():
                string += "Physical Line(%d) = {%s};\n"%(k,",".join(map(str,v)))
            for k, v in surface_ids.items():
                string += "Physical Surface(%d) = {%s};\n"%(k,",".join(map(str,v)))

        journalfile = open(filename, 'w')
        journalfile.write(string)
        journalfile.close()


    def write_simple_journal(self, filename = None, use_ids=True, use_physicals=True):
        """Dump the mesh out to a ,geo file after joining coplanar surfaces."""

        string = "undo off\n"
        string += "reset\n"
        string += "set echo off\n"
        string += "set journal off\n"

        lines = LineDict()
        surfaces = []
        edges = {}
        normals = []
        line_ids = {}
        vertices = set()
        surface_ids = {}
        open_surf ={}


        for k,l in self.elements.items():
            if l[0]==1:
                if use_physicals:
                    line_ids.setdefault(l[1][-1],[]).append(line_no)
                else:
                    line_ids.setdefault(l[1][1],[]).append(line_no)
                line_ids.setdefault(l[1][-1],[]).append(line_no)
                lines.setdefault(l[2][:],len(lines)+1)
            elif l[0]==2:
                if use_physicals:
                    surface_ids.setdefault(l[1][-1],[]).append(k)
                else:
                    surface_ids.setdefault(l[1][1],[]).append(k)
                e1=lines.setdefault([l[2][0],l[2][1]],len(lines)+1)
                e2=lines.setdefault([l[2][1],l[2][2]],len(lines)+1)
                e3=lines.setdefault([l[2][2],l[2][0]],len(lines)+1)
                surfaces.append([e1,e2,e3])

                sid=len(normals)
                edges.setdefault(abs(e1),[]).append(sid)
                edges.setdefault(abs(e2),[]).append(sid)
                edges.setdefault(abs(e3),[]).append(sid)

                n =numpy.cross(self.nodes[l[2][1]].vertices-self.nodes[l[2][0]].vertices,
                               self.nodes[l[2][2]].vertices-self.nodes[l[2][0]].vertices)
                n=n/numpy.sqrt(sum(n**2))

                normals.append(n)

        element_map = range(len(normals))

        for edge, eles in sorted(edges.items()):

            print(edge,eles)

            oeles=None

            while oeles!=eles:
                oeles=eles
                eles=[element_map[eles[0]],element_map[eles[1]]]
                print(eles, oeles)

            eles = sorted(eles)

            print(eles,element_map[eles[0]],element_map[eles[1]])
            
            print(type(surfaces[eles[0]]),type(surfaces[eles[1]]))

            print(abs(numpy.dot(normals[eles[0]],normals[eles[1]])))
            print(1.0-abs(numpy.dot(normals[eles[0]],normals[eles[1]]))<1.0e-8)
            if 1.0-abs(numpy.dot(normals[eles[0]],normals[eles[1]]))<1.0e-8:
                if edge in surfaces[eles[0]]:
                    I = surfaces[eles[0]].index(edge), 1
                else:
                    I = surfaces[eles[0]].index(-edge), -1
                if eles[1]!=eles[0]:
                    if edge in surfaces[eles[1]]:
                        J = surfaces[eles[1]].index(edge), 1
                    else:
                        J = surfaces[eles[1]].index(-edge), -1
                    element_map[eles[1]]=eles[0]
                    surfaces[eles[0]] = surfaces[eles[0]][I[0]+1:]+surfaces[eles[0]][0:I[0]]
                    surfaces[eles[1]] = surfaces[eles[1]][J[0]+1:]+surfaces[eles[1]][:J[0]]
                    if I[1] == J[1]:
                        surfaces[eles[1]].reverse()
                    surfaces[eles[0]].extend(surfaces[eles[1]])
                    surfaces[eles[1]] = None
                    edges.pop(edge)
                else:
                    while edge in surfaces[eles[0]]:
                            surfaces[eles[0]].remove(edge)
                    while -edge in surfaces[eles[0]]:
                            surfaces[eles[0]].remove(-edge)
                    edges.pop(edge)
                    open_surf.setdefault(eles[0],[]).append(edge)
                    

            print(type(surfaces[eles[0]]),type(surfaces[eles[0]]))


        rev_lin={}
        for line, line_id in sorted(lines.items(),key=lambda x:x[1]):
            if line_id in edges:
                vertices.add(line[0])
                vertices.add(line[1])
                rev_lin[line_id]=line[1]
                rev_lin[-line_id]=line[0]

        vmap={}
        i=1
        for k,x in self.nodes.items():
            if k in vertices:
                string += "create vertex %f %f %f\n"%(x[0],x[1],x[2])
                vmap[k] = i
                i += 1

        lmap={}
        i=1
        for line, line_id in sorted(lines.items(),key=lambda x:x[1]):
            if line_id in edges:
                string += "create curve vertex %d %d\n"%(vmap[line[0]],
                                                         vmap[line[1]])
                lmap[line_id] = i
                i += 1


        emap = {}
        i = 1
        for surface_id, surface in enumerate(surfaces):
            if element_map[surface_id]==surface_id:
                string += "create surface curve %s\n"%" ".join(map(str,
                                                                   map(lambda x:lmap[abs(x)],
                                                                       surface)))
                emap[surface_id] = i
                i += 1
                    

        if use_ids:
            for k, v in line_ids.items():
                pass
            for k, v in surface_ids.items():
                V=[V for V in v if element_map[V-1]==V-1]
                V = map(lambda x: emap[x-1], V)
                string += "Sideset %d surface %s\n"%(k," ".join(map(str,V)))

        print(sum(numpy.array(element_map)==range(len(element_map))))

        journalfile = open(filename, 'w')
        journalfile.write(string)
        journalfile.close()

    def write_journal(self, filename = None, use_ids=True, use_physicals=True):
        """Dump the mesh out to a cubit .jou journal file."""

        
        string = "undo off\n"
        string += "reset\n"
        string += "set echo off\n"
        string += "set journal off\n"

        for k,x in self.nodes.items():
            string += "create vertex %f %f %f\n"%(x[0],x[1],x[2])

        lines = LineDict()
        surfaces = []
        edges = {}
        normals = []
        line_ids = {}
        surface_ids = {}


        for k,l in self.elements.items():
            if l[0]==1:
                if use_physicals:
                    line_ids.setdefault(l[1][-1],[]).append(line_no)
                else:
                    line_ids.setdefault(l[1][1],[]).append(line_no)
                line_ids.setdefault(l[1][-1],[]).append(line_no)
                lines.setdefault(l[2][:],len(lines)+1)
            elif l[0]==2:
                if use_physicals:
                    surface_ids.setdefault(l[1][-1],[]).append(k)
                else:
                    surface_ids.setdefault(l[1][1],[]).append(k)
                e1=abs(lines.setdefault([l[2][0],l[2][1]],len(lines)+1))
                e2=abs(lines.setdefault([l[2][1],l[2][2]],len(lines)+1))
                e3=abs(lines.setdefault([l[2][2],l[2][0]],len(lines)+1))
                surfaces.append([e1,e2,e3])

                sid=len(normals)
                edges.setdefault(e1,[]).append(sid)
                edges.setdefault(e2,[]).append(sid)
                edges.setdefault(e3,[]).append(sid)

                n =numpy.cross(self.nodes[l[2][1]].vertices-self.nodes[l[2][0]].vertices,
                               self.nodes[l[2][2]].vertices-self.nodes[l[2][0]].vertices)
                n=n/numpy.sqrt(sum(n**2))

                normals.append(n)

        element_map = range(len(normals))

        for edge, eles in sorted(edges.items()):

            print(edge,eles)

            oeles=None

            while oeles!=eles:
                oeles=eles
                eles=[element_map[eles[0]],element_map[eles[1]]]
                print(eles, oeles)

            eles = sorted(eles)

            print(eles,element_map[eles[0]],element_map[eles[1]])
            
            print(type(surfaces[eles[0]]),type(surfaces[eles[1]]))

            if abs(numpy.dot(normals[eles[0]],normals[eles[1]]))-1.0<1.0e-8:
                surfaces[eles[0]].remove(edge)
                if eles[1]!=eles[0]:
                    surfaces[eles[1]].remove(edge)
                    element_map[eles[1]]=eles[0]
                    surfaces[eles[0]].extend(surfaces[eles[1]])
                    surfaces[eles[1]]=None
                    edges.pop(edge)
                else:
                    surfaces[eles[0]].remove(edge)

            print(type(surfaces[eles[0]]),type(surfaces[eles[0]]))


        for line, line_id in sorted(lines.items(),key=lambda x:x[1]):
            if line_id in edges:
                string += "create curve vertex %d %d\n"%line
        for k,surface in enumerate(surfaces):
            if element_map[k]==k:
                string += "create surface curve %s\n"%" ".join(map(str,surface))
                
        if use_ids:
            for k, v in line_ids.items():
                pass
            for k, v in surface_ids.items():
                string += "Sideset %d surface %s\n"%(k," ".join(map(str,v)))

        journalfile = open(filename, 'w')
        journalfile.write(string)
        journalfile.close()

    def write_vtu(self, filename, **kwargs):
        """Output mesh as a .vtu file with given filename."""

        ugrid = self.as_vtk(**kwargs)

        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetFileName(filename)
        if vtk.VTK_MAJOR_VERSION<6:
            writer.SetInput(ugrid)
        else:
            writer.SetInputData(ugrid)
        writer.Write()

    def write_stl(self, filename, binary=False, **kwargs):
        """Output mesh as a .stl stereo lithography file with given filename."""

        ugrid = self.as_vtk(**kwargs)

        if binary:
            writer.SetFileTypeToBinary()
        else:
            writer.SetFileTypeToASCII()

        writer = vtk.vtkSTLWriter()
        writer.SetFileName(filename)
        if vtk.VTK_MAJOR_VERSION<6:
            writer.SetInput(ugrid)
        else:
            writer.SetInputData(ugrid)
        writer.Write()

    def as_vtk(self, elementary_index=0):
        """Convert to a VTK unstructured grid object, ugrid."""

        etype={15:vtk.VTK_PIXEL,
               1:vtk.VTK_LINE,
               2:vtk.VTK_TRIANGLE,
               4:vtk.VTK_TETRA}
        point_map = {};        
        
        ugrid = vtk.vtkUnstructuredGrid()

        pts = vtk.vtkPoints()


        for i, v in enumerate(self.nodes.items()):
            (node_id, nodes) = v
            
            pts.InsertNextPoint(nodes)
            point_map[node_id]=i


        ugrid.SetPoints(pts)
        ugrid.Allocate(len(self.elements))

        physical_ids = vtk.vtkIntArray()
        physical_ids.SetNumberOfComponents(1)
        physical_ids.SetNumberOfTuples(len(self.elements))
        physical_ids.SetName("PhysicalIds")

        elementary_entities = vtk.vtkIntArray()
        elementary_entities.SetNumberOfComponents(1)
        elementary_entities.SetNumberOfTuples(len(self.elements))
        elementary_entities.SetName("ElementaryEntities")

        for i, v in enumerate(self.elements.items()):

            k, ele = v
            ids = vtk.vtkIdList()

            for node in ele[2]:
                ids.InsertNextId(point_map[node])

            ugrid.InsertNextCell(etype[ele[0]], ids)
            elementary_entities.SetValue(i,ele[1][elementary_index])
            physical_ids.SetValue(i,ele[1][-1])

        ugrid.GetCellData().AddArray(elementary_entities)
        ugrid.GetCellData().AddArray(physical_ids)

        return  ugrid

    def from_vtk(self, ugrid):
        """Convert from a VTK unstructured grid object, ugrid."""

        etype={vtk.VTK_PIXEL:15,
               vtk.VTK_LINE:1,
               vtk.VTK_TRIANGLE:2,
               vtk.VTK_TETRA:4}

        self.reset()

        for i in range(ugrid.GetNumberOfPoints()):
            self.nodes[i+1] = Node(ugrid.GetPoint(i))
        

        for i in range(ugrid.GetNumberOfCells()):

            cell = ugrid.GetCell(i)
            tags = []
            ids = cell.GetPointIds()

            if ugrid.GetCellData().HasArray("ElementaryEntities"):
                tags.append(ugrid.GetCellData().GetArray("ElementaryEntities").GetValue(i))
            if ugrid.GetCellData().HasArray("PhysicalIds"):
                tags.append(ugrid.GetCellData().GetArray("PhysicalIds").GetValue(i))

            if not tags:
                tags.extend((i, i))

            self.elements[i+1] = (etype[cell.GetCellType()],tags,[ids.GetId(_)+1 for _ in range(ids.GetNumberOfIds())])
            self.__add_physical_id__(etype[cell.GetCellType()], tags[-1])

        print('  %d Nodes'%len(self.nodes))
        print('  %d Elements'%len(self.elements))

        return self

    def read_vtu(self, filename):
        """Convert from a VTU file."""
        
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(filename)
        reader.Update()

        self.from_vtu(reader.GetOutput)

    def read_stl(self, filename):
        """Convert from an STL file."""
        
        reader = vtk.vtkSTLReader()
        reader.SetFileName(filename)
        reader.Update()

        self.from_vtu(reader.GetOutput)

    def read_generic(self, filename, reader):
        """Convert from file using a given vtk file reader."""
        
        reader.SetFileName(filename)
        reader.Update()

        self.from_vtu(reader.GetOutput)


    def collapse_ed(self, tol):
        """Collapse and remove any edges smaller than tol"""

        ugrid = self.as_vtk()

        extract = vtk.vtkExtractEdges()

        extract.SetInputData(ugrid)

        extract.Update()

        edges = extract.GetOutput()

        elengths=[]

        loc = vtk.vtkPointLocator()
        loc.SetDataSet(ugrid)
        loc.BuildLocator()
        
        sqrlen=0

        while True:

            for _ in range(edges.GetNumberOfCells()):
                cell = edges.GetCell(_)
                elengths.append(cell.GetLength2())

    def collapse_edges(self, tol):
        """Collapse and remove any edges smaller than tol"""

        ugrid = self.as_vtk()

        extract = vtk.vtkExtractEdges()

        extract.SetInputData(ugrid)

        extract.Update()

        edges = extract.GetOutput()

        loc = vtk.vtkMergePoints()
        pts = vtk.vtkPoints()
        loc.InitPointInsertion(pts, ugrid.GetBounds())
 
        for _ in range(ugrid.GetNumberOfPoints()):
            loc.InsertNextPoint(ugrid.GetPoint(_))
        
        sqrlen=0

        while True:

            elengths=[]

            for _ in range(edges.GetNumberOfCells()):
                cell = edges.GetCell(_)

                pt0 = numpy.array(cell.GetPoints().GetPoint(0))
                pt1 = numpy.array(cell.GetPoints().GetPoint(1))
                elengths.append(sum((pt1-pt0)**2))

            elist = list(enumerate(elengths))
            elist.sort(key=lambda _:_[1])
        
            if sqrlen>tol**2:
                break

            for k,sqrlen in elist:
                if sqrlen == 0.0 :
                    continue
                if sqrlen>tol**2:
                    break
            
                cell = edges.GetCell(k)

                ptt0 = cell.GetPoints().GetPoint(0)
                ptt1 = cell.GetPoints().GetPoint(1)
                
                pid0 = loc.IsInsertedPoint(ptt0)
                pid1 = loc.IsInsertedPoint(ptt1)
                
                pt0 = ugrid.GetPoint(pid0)
                pt1 = ugrid.GetPoint(pid1)

                ptm = [(p0+p1)/2.0 for p0,p1 in zip(pt0,pt1)]
                
                ugrid.GetPoints().SetPoint(pid0,ptm)
                ugrid.GetPoints().SetPoint(pid1,ptm)
                loc.InsertPoint(pid0,ptm)
                loc.InsertPoint(pid1,ptm)
                edges.GetPoints().SetPoint(cell.GetPointIds().GetId(0),ptm)
                edges.GetPoints().SetPoint(cell.GetPointIds().GetId(1),ptm)

        self.coherence(ugrid)
            
        

    def coherence(self, ugrid=None, tol=1.0e-8):
        """Remove duplicate nodes and degenerate elements."""

        def merge_points(ugrid):

            mf = vtk.vtkMergePoints()
            pts = vtk.vtkPoints()
            mf.InitPointInsertion(pts,ugrid.GetBounds())

            point_map = []
            n_nonunique = 0

    

            for i in range(ugrid.GetNumberOfPoints()):
                newid = vtk.mutable(0)
                n_nonunique += mf.InsertUniquePoint(ugrid.GetPoint(i),newid)
                point_map.append(newid)
                
            return pts, point_map

        if ugrid is None:
            ugrid = self.as_vtk()

        pts, point_map = merge_points(ugrid)

        vgrid = vtk.vtkUnstructuredGrid()

        vgrid.SetPoints(pts)

        cell_map = []

        for i in range(ugrid.GetNumberOfCells()):
            cell = ugrid.GetCell(i)
            ids = vtk.vtkIdList()


            if cell.GetCellDimension()==1:
                if cell.GetLength2()<tol:
                    continue
            elif cell.GetCellDimension()==2:
                if cell.ComputeArea()<tol:
                    continue
            elif cell.GetCellDimension()==3:
                if cell.ComputeVolume(*[cell.GetPoints().GetPoint(_)
                                              for _ in range(4)])<tol:
                    continue

            cell_map.append(i)

            for j in range(cell.GetPointIds().GetNumberOfIds()):
                ids.InsertNextId(point_map[cell.GetPointIds().GetId(j)])

            vgrid.InsertNextCell(cell.GetCellType(),ids)

        cd = vgrid.GetCellData()
        icd = ugrid.GetCellData()

        cd.CopyStructure(icd)

        for i in range(icd.GetNumberOfArrays()):
            data = cd.GetArray(i)
            idata = icd.GetArray(i)

            data.SetNumberOfComponents(idata.GetNumberOfComponents())
            data.SetNumberOfTuples(vgrid.GetNumberOfCells())


            for i in range(vgrid.GetNumberOfCells()):

                data.SetTuple(i, cell_map[i], idata)        

        return self.from_vtk(vgrid)
            
            
            
        
