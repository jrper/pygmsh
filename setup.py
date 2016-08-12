import os
from setuptools import setup

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "pygmsh",
    version = "0.0.1",
    author = "James Percival",
    author_email = "j.percival@imperial.ac.uk",
    description = ("A utility for processing gmsh files using VTK"),
    license = "GPL v. 3",
    keywords = "gmsh vtk",
    url = "https://github.com/jrper/pygmsh",
    packages=['pygmsh', 'tests'],
    long_description=read('README.md'),
    classifiers=[
    ],
)
