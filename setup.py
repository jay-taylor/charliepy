from setuptools import setup
from setuptools.extension import Extension

permutat = Extension(name = "charliepy.permutat",
                     sources = ["charliepy/src/permutat/permutat.c"],
                     include_dirs = ["charliepy/include/"])

cdata = Extension(name = "charliepy.data._cdata",
                  sources = ["charliepy/data/src/_cdata.c"],
                  include_dirs = ["charliepy/include/"])

setup(name = 'charliepy',
      version = '0.1',
      description = ('Tools for computing with characters of Lie '
                     'type objects.'),
      url = 'https://github.com/jay-taylor/charliepy',
      author = 'Jay Taylor',
      author_email = 'jaytaylor@math.arizona.edu',
      license = 'GNU GPL',
      packages = ['charliepy', 'charliepy.data'],
      ext_modules = [permutat, cdata],
      zip_safe = False)
