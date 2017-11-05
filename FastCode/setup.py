from distutils.core import setup, Extension
from Cython.Build import cythonize

setup(\
    name='fastVersions',\
    ext_modules = cythonize("/Users/Eddie/HHMM/FastCode/*.pyx")\
)

# from distutils.core import setup, Extension
# from Cython.Build import cythonize


# setup(ext_modules = cythonize(Extension("fastVersions",\
#            sources=["/Users/Eddie/HHMM/FastCode/*.pyx"],\
#            language = "c++",\
#            extra_compile_args=["-std=c++11"]\
#       )))
