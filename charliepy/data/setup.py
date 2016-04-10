from distutils.core import setup
from distutils.extension import Extension

extensions = [Extension("_chartabs", ["src/_chartabs.c"])]

setup(
    ext_modules = extensions
)
