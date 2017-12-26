from distutils.core import setup, Extension
from Cython.Build import cythonize

setup(ext_modules = cythonize(Extension("HGMP",\
           sources=["HGMP.pyx", "LogVar.cpp", \
                             "HypergraphBase.cpp", \
                             "HHMMMessagePassing.cpp"],\
           include_dirs = ["/usr/local/include","."],\
           library_dirs = ["/usr/local/lib"],\
           language = "c++",\
           extra_compile_args=["-std=c++11","-stdlib=libc++"]\
      )))