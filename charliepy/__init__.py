# The import * commands bring the functions (specified with __all__) into the
# charliepy global namespace. In particular such functions (like Coxeter) are
# accesible via charliepy.Coxeter. When "from charliepy import *" is run in
# python these functions are then directly accesible in interactive mode.

from . import core
from .core import *
from . import rootdata 
from .rootdata import *
from . import coxgrp
from .coxgrp import *
from . import coxcos
from .coxcos import *
from . import redgrp
from .redgrp import *
from . import data
