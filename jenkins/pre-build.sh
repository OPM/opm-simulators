#! /bin/bash

./bin/genEvalSpecializations.py

if test -n "$(git diff)"; then
    echo "The generated source files have been manually edited or the "
    echo "code generator has been modified but not been run before "
    echo "proposing the branch for merging."
    exit 1
else
    echo "The generated source files have not been manually edited."
    exit 0
fi
