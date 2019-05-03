from distutils.core import setup
from Cython.Build import cythonize
import numpy


setup(name='helper_methods_cy',
      ext_modules=cythonize("helper_methods_cy.pyx"), include_dirs=[numpy.get_include()])
