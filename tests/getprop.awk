# Extract test property values from CTestTestfile.cmake
#
# User must initialise two awk variables, typically through the "-v" option,
# when invoking this script
#
#   search: Search pattern.  Typically a string such as
#      set_tests_properties($test_name
#
#   prop: Property name, for instance DIRNAME or SIMULATOR.
#
# Property value will be printed to the standard output stream.
#
# The script assumes that $1 on candidate lines is suitable for matching
# against 'search', and that $2 of the matching lines is the word
# 'PROPERTIES'.
#
# Example:
#   # Get value of SIMULATOR property in test named by shell variable
#   # $failed_test and assign this to shell variable 'binary'.
#
#   binary=$(awk -v search="set_tests_properties\\\($failed_test" \
#            -v prop="SIMULATOR" \
#            -f getprop.awk \
#            CTestTestfile.cmake)

$1 ~ search {
    for (i = 3; i <= NF; ++i) {
        if ($i == prop) {
            val = $(i + 1)
            gsub(/"/, "", val)
            print val
            exit
        }
    }
}
