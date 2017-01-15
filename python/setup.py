
from distutils.core import setup
from distutils.extension import Extension
import os.path
import sys
import numpy as np

if sys.platform == 'win32' or sys.platform == 'win64':
    print 'Windows is not a supported platform.'
    quit()

else:
    include_dirs = ['/usr/local/Cellar/python/2.7.9/Frameworks/Python.framework/Versions/2.7/lib/python2.7/../../include/python2.7',
                    '/Users/carguelles/Library/Python/2.7/lib/python/site-packages/numpy/core/include',
                    '/usr/local/include',np.get_include(),
                    '../Likelihood/',
                    '../LikelihoodFit/',
                    '.']
    libraries = ['python2.7','boost_python',
                 'SQuIDS','nuSQuIDS',
                 'gsl','gslcblas','m',
                 'hdf5','hdf5_hl','PhysTools']

    if sys.platform.startswith('linux'):
      libraries.append('cxxrt')

    library_dirs = ['/usr/local/Cellar/python/2.7.9/Frameworks/Python.framework/Versions/2.7/lib/python2.7',
                    '/usr/local/Cellar/python/2.7.9/Frameworks/Python.framework/Versions/2.7/lib/python2.7/../',
                    '/usr/local/lib',
                    '/usr/local/lib',
                    '/usr/local/Cellar/gsl/1.15/lib',
                    '/usr/local/opt/szip/lib',
                    '/usr/local/lib',
                    '/usr/local/lib']

files = ['SterileSearchPy.cpp']

setup(name = 'SterileSearchPy',
      ext_modules = [
          Extension('SterileSearchPy',files,
              library_dirs=library_dirs,
              libraries=libraries,
              include_dirs=include_dirs,
              extra_objects=["../LikelihoodFit/lbfgsb.o","../LikelihoodFit/linpack.o","../LikelihoodFit/SterileSearch.o","../LikelihoodFit/analysisWeighting.o","../Likelihoodfit/dataIO.o","../LikelihoodFit/oversizeWeight.o"],
              extra_compile_args=['-O3','-fPIC','-std=c++11','-Wno-unused-local-typedef'],
              depends=[]),
          ]
      )

