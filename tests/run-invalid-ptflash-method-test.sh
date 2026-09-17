#!/bin/bash
# Copyright 2026 SINTEF Digital
#
# This file is part of the Open Porous Media project (OPM).
#
# OPM is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# OPM is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with OPM.  If not, see <http://www.gnu.org/licenses/>.

set -u

if [ "$#" -ne 4 ] \
    || [ "$1" != "--invalid-method" ] \
    || [ "$2" != "-e" ] \
    || [ "$4" != "--" ]; then
    echo "Usage: $0 --invalid-method -e <co2_ptflash_ecfv> --"
    exit 2
fi

test_binary="$3"
output="$("$test_binary" --flash-two-phase-method=invalid --print-parameters=0 2>&1)"
exit_code=$?

printf '%s\n' "$output"

if [ "$exit_code" -ne 1 ]; then
    echo "Expected exit code 1, got $exit_code"
    exit 1
fi

expected="Unknown two phase flash method invalid is specified"
case "$output" in
    *"$expected"*) ;;
    *)
        echo "Expected diagnostic not found: $expected"
        exit 1
        ;;
esac

exit 0
