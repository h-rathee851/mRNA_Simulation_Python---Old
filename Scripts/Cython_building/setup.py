from distutils.core import setup
from Cython.Build import cythonize
import numpy


setup(name='simulation_up_cy',
      ext_modules=cythonize("simulation_up_cy.pyx"), include_dirs=[numpy.get_include()])
