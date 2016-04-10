# The import * commands bring the functions (specified with __all__) into the
# charliepy global namespace. In particular such functions (like Coxeter) are
# accesible via charliepy.Coxeter. When "from charliepy import *" is run in
# python these functions are then directly accesible in interactive mode.

# These modules have no dependencies on other modules in CharLiePy but do have
# external dependencies.
from . import permutat
from .permutat import *
from . import utils
from .utils import *
from . import data

# Call core first because chartab is dependent on core.
from . import core
from .core import *
from . import conjclass
from .conjclass import *
from . import chartab
from .chartab import *
