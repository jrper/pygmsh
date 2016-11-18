"""Module containing readers, writers and other tools for VTK"""
import vtk

class XMLreader(object):
    """ Read a dolfin .xml geometry file into a vtk unstructured grid object."""

    def __init__(self):
        self._ugrid = None
        self._filename = None

    def SetFileName(self, filename):
        """ Set the file name for the reader."""
        self._filename = filename

    def Update(self):
        """ Update the reader output from the file."""
        import xml.etree.ElementTree as ET
        tree = ET.parse(self._filename)
        root = tree.getroot()
        points = root.find('mesh/vertices')
        cells = root.find('mesh/cells')

        self._ugrid = vtk.vtkUnstructuredGrid()
        self._ugrid.Allocate(0)
        pts = vtk.vtkPoints()

        for point in points:
            pos = [float(point.attrib[_]) for _ in ('x', 'y', 'z')]
            pts.InsertNextPoint(pos)

        self._ugrid.SetPoints(pts)

        vtk_cell = {'triangle':(vtk.VTK_TRIANGLE, 3),
                    'tetrahedron':(vtk.VTK_TETRA, 4)}

        for cell in cells:
            cell_type, npts = vtk_cell[cell.tag]
            val = [int(cell.attrib['v%d'%_]) for _ in range(npts)]
            self._ugrid.InsertNextCell(cell_type, npts, val)

    def GetOutput(self):
        """ Accessor for the output unstructured grid."""
        return self._ugrid

class XMLwriter(object):
    """Writer for dolfin .xml format"""

    def __init__(self):
        self._ugrid = None
        self._filename = None

    def SetFileName(self, filename):
        """ Set the file name for the writer."""
        self._filename = filename

    def SetInput(self, ugrid):
        """ Set the input data."""
        self._ugrid = ugrid

    def SetInputData(self, ugrid):
        """ Set the input data."""
        self._ugrid = ugrid

    def Write(self):
        """ Write the data to file."""
        from lxml import etree as ET

        vtk_cells = {vtk.VTK_TRIANGLE:('triangle', 2, vtk.VTK_TRIANGLE),
                     vtk.VTK_TETRA:('tetrahedron', 3, vtk.VTK_TETRA)}
        celltype = ('', 0, -1)
        cell_ids = []
        ncells = 0

        for _ in range(self._ugrid.GetNumberOfCells()):
            cell = self._ugrid.GetCell(_)
            if vtk_cells[cell.GetCellType()][1] > celltype[1]:
                celltype = vtk_cells[cell.GetCellType()]
                ncells = 1
                cell_ids = [[cell.GetPointIds().GetId(__)
                             for __ in range(cell.GetNumberOfPoints())]]
            elif vtk_cells[cell.GetCellType()][1] == celltype[1]:
                ncells += 1
                cell_ids.append([cell.GetPointIds().GetId(__)
                                 for __ in range(cell.GetNumberOfPoints())])

        tree = ET.ElementTree(element=ET.Element('dolfin', {},
                                                 nsmap={"dolfin":"http://www.fenicsproject.org"}))

        mesh = ET.Element('mesh', {'celltype':str(celltype[0]),
                                   'dim':str(celltype[1])})
        tree.getroot().append(mesh)
        vertices = ET.Element('vertices', {'size':str(self._ugrid.GetNumberOfPoints())})
        mesh.append(vertices)
        for _ in range(self._ugrid.GetNumberOfPoints()):
            pos = self._ugrid.GetPoint(_)
            vertices.append(ET.Element('vertex', {'index':str(_),
                                                  'x':'%1.16e'%pos[0],
                                                  'y':'%1.16e'%pos[1],
                                                  'z':'%1.16e'%pos[2]}))

        cells_element = ET.Element('cells', {'size':str(ncells)})
        mesh.append(cells_element)
        for _ in enumerate(cell_ids):
            data = {'index':str(_[0])}
            for __ in enumerate(_[1]):
                data['v%d'%__[0]] = str(__[1])
            ET.SubElement(cells_element, celltype[0], data)

        tree.write(self._filename, encoding="UTF-8",
                   xml_declaration=True,
                   pretty_print=True)
