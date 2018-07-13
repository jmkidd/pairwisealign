from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

setup(
    ext_modules = [
      Extension("pairwisealign_c1", ["pairwisealign_c1.pyx"],
      include_dirs = ['.',numpy.get_include()])
    ],
    cmdclass = {'build_ext': build_ext},
 )


setup(
    ext_modules = [
      Extension("pairwisealign_c2", ["pairwisealign_c2.pyx"],
      include_dirs = ['.',numpy.get_include()])
    ],
    cmdclass = {'build_ext': build_ext},
 )

setup(
    ext_modules = [
      Extension("pairwisealign_c3", ["pairwisealign_c3.pyx"],
      include_dirs = ['.',numpy.get_include()])
    ],
    cmdclass = {'build_ext': build_ext},
 )

setup(
    ext_modules = [
      Extension("pairwisealign_c4", ["pairwisealign_c4.pyx"],
      include_dirs = ['.',numpy.get_include()])
    ],
    cmdclass = {'build_ext': build_ext},
 )


setup(
    ext_modules = [
      Extension("pairwisealign_c5", ["pairwisealign_c5.pyx"],
      include_dirs = ['.',numpy.get_include()])
    ],
    cmdclass = {'build_ext': build_ext},
 )
