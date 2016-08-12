""" Module containing an expanded python gmsh class"""
from __future__ import print_function

import vtk
from vtk.util import numpy_support
import numpy
import struct



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
        self.filename = filename
        if self.filename:
            self.read()

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
                            self.nodes[int(columns[0])] = map(float, columns[1:])
                        elif ftype==1:
                            nnods=int(columns[0])
                            for N in range(nnods):
                                data=mshfile.read(4+3*dsize)
                                i,x,y,z=struct.unpack('=i3d',data)
                                self.nodes[i]=(x,y,z)
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
                        self.elements[id] = (type, tags, nodes)
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
                                self.elements[data[0]]=(etype,
                                                        data[1:1+ntags],
                                                        data[1+ntags:])
                    except:
                        raise
                    mshfile.read(1)
                            
        print('  %d Nodes'%len(self.nodes))
        print('  %d Elements'%len(self.elements))

        mshfile.close()

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
            (ele_type, tags, nodes) = elem
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
            (ele_type, tags, nodes) = elem
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

    def write_geo(self, filename = None, use_ids=True, use_physicals=True):
        """Dump the mesh out to a Gmsh .geo geometry file."""

        string=""

        for k,x in self.nodes.items():
            string += 'Point(%d) = {%f,%f,%f};\n'%(k,x[0],x[1],x[2])

        line_no = 1;
        line_ids = {}
        surface_ids = {}

        for k,l in self.elements.items():
            if l[0]==1:
                if use_physicals:
                    line_ids.setdefault(l[1][-1],[]).append(line_no)
                else:
                    line_ids.setdefault(l[1][1],[]).append(line_no)
                line_ids.setdefault(l[1][-1],[]).append(line_no)
                string += "Line(%d) = {%d,%d};\n"%(line_no,l[2][0],l[2][1])
                line_no+=1
            elif l[0]==2:
                if use_physicals:
                    surface_ids.setdefault(l[1][-1],[]).append(k)
                else:
                    surface_ids.setdefault(l[1][1],[]).append(k)
                string += "Line(%d) = {%d,%d};\n"%(line_no,l[2][0],l[2][1])
                string += "Line(%d) = {%d,%d};\n"%(line_no+1,l[2][1],l[2][2])
                string += "Line(%d) = {%d,%d};\n"%(line_no+2,l[2][2],l[2][0])
                string += "Line Loop(%d) = {%d,%d,%d};"%tuple([k]+range(line_no,line_no+3))
                string += "Plane Surface(%d) = %d;"%(k,k)
                line_no+=3
                
        if use_ids:
            for k, v in line_ids.items():
                string += "Physical Line(%d) = {%s};\n"%(k,",".join(map(str,v)))
            for k, v in surface_ids.items():
                string += "Physical Surface(%d) = {%s};\n"%(k,",".join(map(str,v)))

        geofile = open(filename, 'w')
        geofile.write(string)
        geofile.close()



    def as_vtk(self,elementary_index=0):
        """Convert to a VTK unstructured grid object, ugrid."""

        etype={1:vtk.VTK_LINE,
               2:vtk.VTK_TRIANGLE}
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

        etype={vtk.VTK_LINE:1,
               vtk.VTK_TRIANGLE:1}

        self.reset()

        for i in range(ugrid.GetNumberOfPoints()):
            self.nodes[i] = ugrid.GetPoint(i)
        

        for i in range(ugrid.GetNumberOfCells()):

            cell = ugrid.GetCell(i)
            tags = []
            ids = cell.GetPointIds()

            if ugrid.GetCellData().HasArray("ElementaryEntities"):
                tags.append(ugrid.GetCellData().GetArray("ElementaryEntities").GetValue(i))
            if ugrid.GetCellData().HasArray("PhysicalIds"):
                tags.append(ugrid.GetCellData().GetArray("PhysicalIds").GetValue(i))

            self.elements[i] = (etype[cell.GetCellType()],tags,[ids.GetId(_) for _ in range(ids.GetNumberOfIds())])


        print('  %d Nodes'%len(self.nodes))
        print('  %d Elements'%len(self.elements))

    def coherence(self):

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

        ugrid = self.as_vtk()

        pts, point_map = merge_points(ugrid)

        vgrid = vtk.vtkUnstructuredGrid()

        vgrid.SetPoints(pts)

        cell_map = []

        for i in range(ugrid.GetNumberOfCells()):
            cell = ugrid.GetCell(i)
            ids = vtk.vtkIdList()

            if cell.ComputeArea()==0.0:
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
            
            
            
        