
import pytest
import pygmsh

def test_blank_init():
    m = pygmsh.GmshMesh()
    assert(m)

def test_init_from_ascii_file():
    m = pygmsh.GmshMesh('tests/test.msh')
    assert(m)

def test_init_from_binary_file():
    m = pygmsh.GmshMesh('tests/test_bin.msh')
    assert(m)

def test_read_from_ascii_file():
    m = pygmsh.GmshMesh()
    f = open('tests/test.msh','r')
    m.read(f)
    assert(m)
    
def test_read_from_binary_file():
    m = pygmsh.GmshMesh()
    f = open('tests/test_bin.msh','r')
    m.read(f)
    assert(m)

def test_vtu_conversion():
    m = pygmsh.GmshMesh('tests/test.msh')
    ugrid = m.as_vtk()

    assert(ugrid.GetNumberOfPoints() == m.nodecount())
    assert(ugrid.GetNumberOfCells() == m.elementcount())

    

    
