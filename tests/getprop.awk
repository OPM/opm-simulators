# Extract test property values from CTestTestfile.cmake
#
# User must initialise two awk variables, typically through the "-v" option,
# when invoking this script
#
#   test: Test name.
#
#   prop: Property name, for instance DIRNAME or SIMULATOR.
#
# Property value will be printed to the standard output stream.
#
# Example:
#   # Get value of SIMULATOR property in test named by shell variable
#   # $failed_test and assign this to shell variable 'binary'.
#
#   binary=$(awk -v test="${failed_test}" \
#            -v prop="SIMULATOR" \
#            -f getprop.awk \
#            CTestTestfile.cmake)

BEGIN {
    # Validate that the user provided the required variables
    if (!test || !prop) {
        print "Error: Please provide both 'test' and 'prop' variables." > "/dev/stderr"
        print "Usage: awk -v test=\"NAME\" -v prop=\"PROP\" -f getprop.awk <file>" > "/dev/stderr"
        exit 1
    }

    # Construct exact anchor patterns for the test name
    bracket_target = "[=[" test "]=]"
    quoted_target  = "\"" test "\""

    # Unquoted test names must be followed by at least one blank/space
    unquoted_target = test " "
}

{
    # 1. Detect block start using strict literal string matches
    if (index($0, "set_tests_properties(") &&
        (index($0, bracket_target) ||
         index($0, quoted_target)  ||
         index($0, unquoted_target)))
    {
        in_target_block = 1
    }

    # 2. Extract target property value while inside the target block
    if (in_target_block && index($0, prop)) {
        # Match standard quoted string: PROPERTY "value" (e.g., SIMULATOR "flow")
        if (match($0, prop "[[:space:]]+\"[^\"]+\"")) {
            split(substr($0, RSTART, RLENGTH), parts, "\"")
            print parts[2]
            exit
        }
        # Match CMake bracket argument: PROPERTY [=[value]=] (e.g., SIMULATOR [=[flow]=])
        else if (match($0, prop "[[:space:]]+\\[=\\[[^\\]]+\\]=\\]")) {
            split(substr($0, RSTART, RLENGTH), parts, /\[=\[|\]=\]/)
            print parts[2]
            exit
        }
        # Match unquoted argument: PROPERTY value (e.g., SIMULATOR flow)
        else if (match($0, prop "[[:space:]]+[^[:space:]]+")) {
            # Split on FS--i.e., whitespace.
            split(substr($0, RSTART, RLENGTH), parts)
            print parts[2]
            exit
        }
    }

    # 3. Reset flag at the end of the command block
    if (in_target_block && index($0, ")")) {
        in_target_block = 0
    }
}
