
from distutils.core import setup
from distutils.extension import Extension
import os,os.path
import sys
import numpy as np

if sys.platform == 'win32' or sys.platform == 'win64':
    print 'Windows is not a supported platform.'
    quit()

else:
    include_dirs = [os.environ['SNOTBUILDPATH']+"/include",
                    '/usr/local/include',
                    np.get_include(),
                    '../Likelihood/',
                    '../LikelihoodFit/',
                    '.']
    libraries = ['python2.7','boost_python','LeptonWeighter',
                 'SQuIDS','nuSQuIDS',
                 'gsl','gslcblas','m',
                 'hdf5','hdf5_hl','PhysTools']

    if sys.platform.startswith('linux'):
      libraries.append('supc++')#'cxxrt'

    library_dirs = [os.environ['SNOTBUILDPATH']+"/lib",
                    '/usr/local/lib']

files = ['SterileSearchPy.cpp']

setup(name = 'SterileSearchPy',
      ext_modules = [
          Extension('SterileSearchPy',files,
              library_dirs=library_dirs,
              libraries=libraries,
              include_dirs=include_dirs,
              extra_objects=["../LikelihoodFit/analysisWeighting.o","../LikelihoodFit/Event.o","../LikelihoodFit/linpack.o","../LikelihoodFit/SterileSearch.o",
                             "../LikelihoodFit/compactIO.o","../LikelihoodFit/lbfgsb.o","../LikelihoodFit/oversizeWeight.o"],
              extra_compile_args=['-O3','-fPIC','-std=c++11','-Wno-unused-local-typedef'],
              depends=[]),
          ]
      )

