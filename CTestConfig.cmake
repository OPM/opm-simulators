# this is included after opm-core_NAME is set
set(CTEST_PROJECT_NAME "${${project}_NAME}")
set(CTEST_NIGHTLY_START_TIME "01:00:00 UTC")
set(CTEST_DROP_METHOD "http")
set(CTEST_DROP_SITE "opm-project.org")
set(CTEST_DROP_LOCATION "/CDash/submit.php?project=${${project}_NAME}")
set(CTEST_DROP_SITE_CDASH TRUE)
