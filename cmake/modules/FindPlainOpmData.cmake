# - Try to find the an un-installed version of the opm-data module
# this rule will only do something if a properly installed version
# of opm-data has not been found before.
#
#  opm-data_FOUND     - the opm-data module is available
#  opm-data_DIR       - path to the toplevel directory of the plain opm-data tree
#  opm-data_INSTALLED - Specifies whether the opm-data module was
#                       installed on the system or is a plain tree
#
#  Usage:
#  find_package(PlainOpmData)

if (opm-data_FOUND)
  set(opm-data_INSTALLED TRUE)
else()
  set(opm-data_FOUND FALSE)
  set(opm-data_INSTALLED FALSE)
  set(opm-data_DIR NOTFOUND)
  if (opm-data_DIR AND EXISTS "${opm-data_DIR}/spe1/SPE1CASE2.DATA")
    set(opm-data_FOUND TRUE)
    set(opm-data_INSTALLED FALSE)
    set(opm-data_DIR "${opm-data_DIR}")
  elseif (OPM_DATA_DIR AND EXISTS "${OPM_DATA_DIR}/spe1/SPE1CASE2.DATA")
    set(opm-data_FOUND TRUE)
    set(opm-data_INSTALLED FALSE)
    set(opm-data_DIR "${OPM_DATA_DIR}")
  elseif (OPM_DATA_ROOT AND EXISTS "${OPM_DATA_ROOT}/spe1/SPE1CASE2.DATA")
    set(opm-data_FOUND TRUE)
    set(opm-data_INSTALLED FALSE)
    set(opm-data_DIR "${OPM_DATA_ROOT}")
  endif()
endif()

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  opm-data
  "opm-data could not be found. You can set opm-data_DIR to rectify this." opm-data_FOUND)

# stop the dune build system from complaining if a suitable opm-data
# instance was found
if (opm-data_DIR)
  set(PlainOpmData_FOUND 1)
endif()
