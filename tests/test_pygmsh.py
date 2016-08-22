import pytest
import pygmsh

def test_blank_init():
    m = pygmsh.GmshMesh()
    assert(m)

def test_init_from_ascii_file():
    m = pygmsh.GmshMesh('tests/test.msh')
    assert(m)
    assert(m.nodecount()>0)
    assert(m.elementcount()>0)

def test_init_from_binary_file():
    m = pygmsh.GmshMesh('tests/test_bin.msh')
    assert(m)
    assert(m.nodecount()>0)
    assert(m.elementcount()>0)

def test_read_from_ascii_file():
    m = pygmsh.GmshMesh()
    f = open('tests/test.msh','r')
    m.read(f)
    assert(m)
    assert(m.nodecount()>0)
    assert(m.elementcount()>0)
    
def test_read_from_binary_file():
    m = pygmsh.GmshMesh()
    f = open('tests/test_bin.msh','r')
    m.read(f)
    assert(m)
    assert(m.nodecount()>0)
    assert(m.elementcount()>0)

def test_conversion_to_vtu():
    m = pygmsh.GmshMesh('tests/test.msh')
    ugrid = m.as_vtk()

    assert(ugrid.GetNumberOfPoints() == m.nodecount())
    assert(ugrid.GetNumberOfCells() == m.elementcount())


def test_conversion_from_vtu():

    reader = pygmsh.vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName('tests/test.vtu')
    reader.Update()
    ugrid = reader.GetOutput()

    m = pygmsh.GmshMesh()
    m.from_vtk(ugrid)

    assert(ugrid.GetNumberOfPoints() == m.nodecount())
    assert(ugrid.GetNumberOfCells() == m.elementcount())

def test_concatentation():
    m1 = pygmsh.GmshMesh('tests/test.msh')
    m2 = pygmsh.GmshMesh('tests/test_bin.msh')

    m3 = m1 + m2

    assert(m3.nodecount() == m1.nodecount() + m2.nodecount())
    assert(m3.elementcount() == m1.elementcount() + m2.elementcount())


def test_inplace_concatentation():
    m1 = pygmsh.GmshMesh('tests/test.msh')
    m2 = pygmsh.GmshMesh('tests/test_bin.msh')

    n_nodes = m1.nodecount()
    n_eles = m1.elementcount()

    m1 += m2

    assert(m1.nodecount() == n_nodes + m2.nodecount())
    assert(m1.elementcount() == n_eles + m2.elementcount())

def test_coherence():

    m1 = pygmsh.GmshMesh('tests/test.msh')
    m2 = pygmsh.GmshMesh('tests/test.msh')

    m1 += m2

    assert(m1.nodecount() == 2*m2.nodecount())

    m1.coherence()

    assert(m1.nodecount() == m2.nodecount())
#   need to think what coherence should mean here:
#    assert(m1.elementcount() == m2.elementcount())


def test_reset():
    m = pygmsh.GmshMesh('tests/test.msh')
    m.reset()

    assert(m.nodecount() == 0)
    assert(m.elementcount() == 0)

def test_element_insertion():

    m = pygmsh.GmshMesh()

    assert(m.nodecount()==0)

    m.insert_node(0,(0,0,0))

    assert(m.nodecount()==1)

def test_element_insertion():

    m = pygmsh.GmshMesh()

    assert(m.nodecount()==0)

    m.insert_element(0,1,(0,1),(0,1))

    assert(m.elementcount()==1)
    

def test_transformation():

    m = pygmsh.GmshMesh()

    m.insert_node(0,(0,0,0))
    m.insert_node(1,(1,0,0))
    m.insert_node(2,(0,1,0))

    def func(x):
        return (x[1],-x[0],0)

    m.transform(func)

    assert(m.nodes[0] == (0,0,0))
    assert(m.nodes[1] == (0,-1,0))
    assert(m.nodes[2] == (1,0,0))

