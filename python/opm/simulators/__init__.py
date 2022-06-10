# Instead of having the pybind11 extension module, e.g.
#   simulators.cpython-310-x86_64-linux-gnu.so  located
#   directly in the opm directory, we create a package (sub
#   directory) with the same name and place it there.
#   In this way we can do (if needed in the future)
#
#  from opm.simulators import BlackOilSimulator, FoamSimulator, PurePythonUtils, ...
#
#  where FoamSimulator and PurePythonUtils does not currently exists,
#  but could be possible future extensions..
#
from .simulators import BlackOilSimulator
