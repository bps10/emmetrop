from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
from Cython.Distutils import build_ext
import sys
"""
ext_modules = [Extension("SchematicEye",
                     ["SchematicEye.pyx"],
                     SchematicEye='c++',
                     include_dirs=[r'.', r'C:/Cplusplus/gsl-1.15'],
                     library_dirs=[r'.', r'C:/Cplusplus/gsl-1.15'],
                     extra_objects=["SchematicEye.cc", "libgsl.la"],
                     )]
                     
setup(
      cmdclass = {'build_ext': build_ext},
      ext_modules = ext_modules
      )
""" 


if sys.platform == "darwin" :
    include_gsl_dir = "/usr/local/include/"
    lib_gsl_dir = "/usr/local/lib/"
    
elif sys.platform == "win32":    
    include_dir = r"c:\cygwin\usr\local\include"
    lib_dir = r"c:\cygwin\usr\local\lib" 
    
"""     
ext = Extension("SchematicEye", ["SchematicEye.pyx"],
    include_dirs=[include_gsl_dir],
    library_dirs=[lib_gsl_dir],
    libraries=["gsl", "Goptical"]
)

setup(ext_modules=[ext],
    cmdclass = {'build_ext': build_ext})
"""
  
setup(ext_modules = cythonize(
            "SchematicEye.pyx",                 # our Cython source
 
            sources=["SchematicEye.cc"],  # additional source file(s)
            language="c++"
            ))
