import pytest
from pygmsh import *

def test_blank_init():
    m = GmshMesh()
    assert(m)

def test_init_from_ascii_file():
    m = GmshMesh('pygmsh/tests/test.msh')
    assert(m)
    assert(m.nodecount()>0)
    assert(m.elementcount()>0)

def test_init_from_binary_file():
    m = GmshMesh('pygmsh/tests/test_bin.msh')
    assert(m)
    assert(m.nodecount()>0)
    assert(m.elementcount()>0)

def test_read_from_ascii_file():
    m = GmshMesh()
    f = open('pygmsh/tests/test.msh','r')
    m.read(f)
    assert(m)
    assert(m.nodecount()>0)
    assert(m.elementcount()>0)
    
def test_read_from_binary_file():
    m = GmshMesh()
    f = open('pygmsh/tests/test_bin.msh','r')
    m.read(f)
    assert(m)
    assert(m.nodecount()>0)
    assert(m.elementcount()>0)

def test_read_from_vtu():
    m = GmshMesh()
    m.read_vtu('pygmsh/tests/test.vtu')
    assert(m)
    assert(m.nodecount()>0)
    assert(m.elementcount()>0)

def test_read_from_stl():
    m = GmshMesh()
    m.read_stl('pygmsh/tests/test.stl')
    assert(m)
    assert(m.nodecount()>0)
    assert(m.elementcount()>0)

def test_read_from_xml():
    m = GmshMesh()
    m.read_generic('pygmsh/tests/test.xml', XMLreader())
    assert(m)
    assert(m.nodecount()>0)
    assert(m.elementcount()>0)

def test_write_ascii(tmpdir):
    m = GmshMesh('pygmsh/tests/test.msh')
    m.write_ascii(filename=tmpdir.join('test.msh').strpath)

def test_write_binary(tmpdir):
    m = GmshMesh('pygmsh/tests/test.msh')
    m.write_binary(filename=tmpdir.join('test.msh').strpath)

def test_write_geo(tmpdir):
    m = GmshMesh('pygmsh/tests/test.msh')
    m.write_geo(filename=tmpdir.join('test.geo').strpath)

def test_write_simple_geo(tmpdir):
    m = GmshMesh('pygmsh/tests/test.msh')
    m.write_simple_geo(filename=tmpdir.join('test.geo').strpath)

def test_write_journal(tmpdir):
    m = GmshMesh('pygmsh/tests/test.msh')
    m.write_journal(filename=tmpdir.join('test.e').strpath)

def test_write_simple_journal(tmpdir):
    m = GmshMesh('pygmsh/tests/test.msh')
    m.write_simple_journal(filename=tmpdir.join('test.e').strpath)

def test_write_vtu(tmpdir):
    m = GmshMesh('pygmsh/tests/test.msh')
    m.write_vtu(filename=tmpdir.join('test.vtu').strpath)

def test_write_stl(tmpdir):
    m = GmshMesh('pygmsh/tests/test.msh')
    m.write_stl(filename=tmpdir.join('test.stl').strpath)

def test_write_generic(tmpdir):
    m = GmshMesh('pygmsh/tests/test.msh')
    m.write_generic(filename=tmpdir.join('test.vtu').strpath, writer=XMLwriter())

def test_conversion_to_vtu():
    m = GmshMesh('pygmsh/tests/test.msh')
    ugrid = m.as_vtk()

    assert(ugrid.GetNumberOfPoints() == m.nodecount())
    assert(ugrid.GetNumberOfCells() == m.elementcount())

def test_conversion_from_vtu():

    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName('pygmsh/tests/test.vtu')
    reader.Update()
    ugrid = reader.GetOutput()

    m = GmshMesh()
    m.from_vtk(ugrid)

    assert(ugrid.GetNumberOfPoints() == m.nodecount())
    assert(ugrid.GetNumberOfCells() == m.elementcount())

def test_concatentation():
    m1 = GmshMesh('pygmsh/tests/test.msh')
    m2 = GmshMesh('pygmsh/tests/test_bin.msh')

    m3 = m1 + m2

    assert(m3.nodecount() == m1.nodecount() + m2.nodecount())
    assert(m3.elementcount() == m1.elementcount() + m2.elementcount())

def test_inplace_concatentation():
    m1 = GmshMesh('pygmsh/tests/test.msh')
    m2 = GmshMesh('pygmsh/tests/test_bin.msh')

    n_nodes = m1.nodecount()
    n_eles = m1.elementcount()

    m1 += m2

    assert(m1.nodecount() == n_nodes + m2.nodecount())
    assert(m1.elementcount() == n_eles + m2.elementcount())

def test_coherence():

    m1 = GmshMesh('pygmsh/tests/test.msh')
    m2 = GmshMesh('pygmsh/tests/test.msh')

    m1 += m2

    assert(m1.nodecount() == 2*m2.nodecount())

    m1.coherence()

    assert(m1.nodecount() == m2.nodecount())
#   need to think what coherence should mean here:
#    assert(m1.elementcount() == m2.elementcount())

def test_reset():
    m = GmshMesh('pygmsh/tests/test.msh')
    m.__reset__()

    assert(m.nodecount() == 0)
    assert(m.elementcount() == 0)

def test_node_insertion():

    m = GmshMesh()
    assert(m.nodecount()==0)
    m.__insert_node__(0,(0,0,0))
    assert(m.nodecount()==1)

def test_element_insertion():

    m = GmshMesh()
    assert(m.nodecount()==0)
    m.__insert_element__(0,1,(0,1),(0,1))
    assert(m.elementcount()==1)
    

def test_transformation():

    m = GmshMesh()

    m.__insert_node__(0,(0,0,0))
    m.__insert_node__(1,(1,0,0))
    m.__insert_node__(2,(0,1,0))

    def func(x):
        return (x[1],-x[0],0)

    m.transform(func)

    assert(m.nodes[0] == (0,0,0))
    assert(m.nodes[1] == (0,-1,0))
    assert(m.nodes[2] == (1,0,0))

