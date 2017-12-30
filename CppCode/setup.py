from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy

setup(ext_modules = cythonize(Extension("HGMP",\
           sources=["HGMP.pyx", "LogVar.cpp", \
                                "Node.cpp", \
                                "Edge.cpp", \
                                "HyperGraph.cpp", \
                                "util.cpp"], \
           include_dirs = ["/usr/local/include","."],\
           library_dirs = ["/usr/local/lib"],\
           language = "c++",\
           extra_compile_args=["-std=c++11","-stdlib=libc++"]\
      )))