from distutils.core import setup
from distutils.core import Extension

permutat = Extension(name = "charliepy.permutat",
                     sources = ["src/permutat.c"],
                     depends = ["src/permutat_mul_templates.c",
                                "src/permutat_templates.c",
                                "src/temputils.h"])

setup(name = 'charliepy',
      version = '0.1',
      description = ('Tools for computing with characters of Lie '
                     'type objects.'),
      url = 'http://www.mathematik.uni-stuttgart.de/~geckmf/',
      author = 'Jay Taylor',
      author_email = 'j.taylor.maths@gmail.com',
      license = 'GNU GPL',
      packages = ['charliepy'],
      ext_modules = [permutat],
      zip_safe = False)
