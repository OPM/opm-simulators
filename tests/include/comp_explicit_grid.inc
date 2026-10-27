-- This reservoir simulation deck is made available under the Open Database
-- License: http://opendatacommons.org/licenses/odbl/1.0/. Any rights in
-- individual contents of the database are licensed under the Database Contents
-- License: http://opendatacommons.org/licenses/dbcl/1.0/

-- Copyright (C) 2024, SINTEF Digital, TNO
-- Copyright (C) 2026, SINTEF Digital

------------------------------------------------------------------------
INIT

DXV
 9*10
/
DYV
 10
/
DZV
 1
/

BOX
1 9 1 1 1 1 /
TOPS
9*10
/
ENDBOX

PERMX
9*100
/

PERMY
9*100
/

PERMZ
9*100
/

PORO
9*0.1
/

MULTX
9*0
/

------------------------------------------------------------------------
