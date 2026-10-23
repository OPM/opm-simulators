-- This reservoir simulation deck is made available under the Open Database
-- License: http://opendatacommons.org/licenses/odbl/1.0/. Any rights in
-- individual contents of the database are licensed under the Database Contents
-- License: http://opendatacommons.org/licenses/dbcl/1.0/

-- Copyright (C) 2024, SINTEF Digital, TNO
-- Copyright (C) 2026, SINTEF Digital

------------------------------------------------------------------------

CNAMES
CO2
METHANE
DECANE
/

ROCK
68.9476
0
/

DENSITY
1* 1000.0 1* /


EOS
SRK /

BIC
0
0
0
/

ACF
0.22394
0.01142
0.4884
/
PCRIT
73.773
45.992
21.03
/
TCRIT
304.128
190.564
617.7
/
MW
44.00
16.04
142.28
/

VCRIT
0.09412
0.09863
0.60980
/

SGOF
-- Sg    Krg    Kro    Pcgo
   0.0   0.0    1.0    0.0
   0.1   0.1    0.9    0.0
   0.2   0.2    0.8    0.0
   0.3   0.3    0.7    0.0
   0.4   0.4    0.6    0.0
   0.5   0.5    0.5    0.0
   0.6   0.6    0.4    0.0
   0.7   0.7    0.3    0.0
   0.8   0.8    0.2    0.0
   0.9   0.9    0.1    0.0
   1.0   1.0    0.0    0.0
/

STCOND
15.0 1.0 /

RTEMP
150 /

------------------------------------------------------------------------
