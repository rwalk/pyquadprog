# -*- coding: utf-8 -*-
"""
Created on Tue Feb  3 21:15:19 2015
@author: rwalker (r_walker@zoho.com)
"""
from __future__ import division, absolute_import, print_function
from numpy.distutils.core import setup, Extension
from numpy.distutils.misc_util import Configuration
from numpy.distutils.system_info import get_info

# fortran/f2py source files
fsource = ['src/pyquadprog.pyf', 
    'src/achck.f', 
    'src/aind.f', 
    'src/solve.QP.compact.f', 
    'src/solve.QP.f', 
    'src/util.f']

quadprog_ext = Extension( name = 'pyquadprog', sources=fsource, libraries=["blas","lapack"])

if __name__ == "__main__":
    setup(name='quadprog',
      version='0.1',
      description='Python binding for Berwin Turlach\'s quadprog routine.',
      url='http://github.com/rwalk333/pyquadprog',
      author='rwalker',
      author_email='r_walker@zoho.com',
      py_modules = ['quadprog'],
      license='GPL V2',
      ext_modules=[quadprog_ext])
