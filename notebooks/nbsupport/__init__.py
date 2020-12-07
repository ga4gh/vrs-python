"""jupyter notebook support for VRS

Within jupyter:

  import nbsupport
  from mynotebook import a_variable, a_function

  a_function(a_variable)

Reece Hart, 2020

"""

from .logging import set_log_level
from .nbimporter import enable_nb_import
from .sequences import translate_sequence_identifier
from .utils import ppo


set_log_level("INFO")
enable_nb_import()
    
