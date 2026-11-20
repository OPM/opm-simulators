# tomlc17

TOML v1.1 parser used by OPM Flow to read configuration files ending in `.toml`
(see `Opm::PropertyTree`).

* Upstream: https://github.com/cktan/tomlc17
* Version: release tag `R261003` (commit `8d3766ddafa6515417f4885ae6878ef69fac401e`)
* License: MIT, see [LICENSE](LICENSE). The MIT license is compatible with the
  GPL-3.0-or-later license of OPM.

`tomlc17.h` and `tomlc17.c` are copied from the `src/` directory of the upstream
repository, as suggested by the upstream documentation. To update, replace both files
and `LICENSE` with the versions from a newer upstream release, update the version
information above, and re-apply the local modifications below unless upstream has
fixed the issues.

## Local modifications

* `tomlc17.c`: `scan_mark()` and `scan_restore()` also save and restore the
  `bracket_level` and `brace_level` nesting counters (fields added to
  `scanner_state_t`, marked "OPM"). Upstream only rewinds the cursor, so the
  lookahead after each table header decremented `bracket_level` twice, which let
  documents bypass the nesting limit that protects against stack overflow. Covered by
  the `Nesting_Limit_After_Table_Headers` test in `tests/test_propertytree.cpp`.
