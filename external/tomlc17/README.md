# tomlc17

TOML v1.1 parser used by OPM Flow to read configuration files ending in `.toml`
(see `Opm::PropertyTree`).

* Upstream: https://github.com/cktan/tomlc17
* Version: release tag `R261003` (commit `8d3766ddafa6515417f4885ae6878ef69fac401e`)
* License: MIT, see [LICENSE](LICENSE). The MIT license is compatible with the
  GPL-3.0-or-later license of OPM.

`tomlc17.h` and `tomlc17.c` are copied unmodified from the `src/` directory of the
upstream repository, as suggested by the upstream documentation. To update, replace
both files and `LICENSE` with the versions from a newer upstream release and update
the version information above.
