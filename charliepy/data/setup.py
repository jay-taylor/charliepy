from distutils.core import setup
from distutils.extension import Extension

try:
    from Cython.Build import cythonize
    USE_CYTHON = True
except ImportError:
    USE_CYTHON = False

if USE_CYTHON:
    extensions = [Extension("_chartabs", ["src/_chartabs.pyx"])]
    extensions = cythonize(extensions)
else:
    extensions = [Extension("_chartabs", ["src/_chartabs.c"])]

setup(
    ext_modules = extensions
)
