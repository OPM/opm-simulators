![OPM logo](https://github.com/OPM/.github/blob/main/profile/OPM%20logo%20small%20cropped.png)
# Open Porous Media Simulators
[![Home Page](https://img.shields.io/badge/Home_Page-opm--project.org-8888FF)](https://opm-project.org/)
[![License](https://img.shields.io/badge/License-GPL_V.3%2B-2ea44f)](#)
[![Documentation](https://img.shields.io/badge/docs-manual-999999)](https://opm-project.org/?page_id=955)
[![Main Article - 10.1016/j.camwa.2020.05.014](https://img.shields.io/badge/Main_Article-10.1016%2Fj.camwa.2020.05.014-FF9900)](https://doi.org/10.1016/j.camwa.2020.05.014)
[![Zenodo](https://img.shields.io/badge/Zenodo_Src-10.5281%2Fzenodo.12637570-ff2222)](https://doi.org/10.5281/zenodo.12637570)

> The `opm-simulators` repository contains the OPM Flow reservoir simulator which uses automatic differentiation and standard input-formats to simulate reservoirs of industrial complexity on real assets by commercial actors. It is cooperatively developed by multiple industry-partners and academic institutions under the GPL 3 license.
> 
> Flow is a fully implicit black-oil simulators which also supports CO2 storage, H2 storage, thermal simulations, solvents, as well as polymers. It works with the Eclipse formats, making it easy to integrate into existing workflows. The automatic differentiation is implemented in [opm-common](https://github.com/OPM/opm-common). The manual can be found here on the OPM [home page](https://opm-project.org/?page_id=955).

## Supported Platforms
- RHEL 8+
- Ubuntu 22.04 LTS and 24.04 LTS
- Debian
- Mac OS X
- WSL2

[Download binaries here](https://opm-project.org/?page_id=36)


## Building from Source
`opm-simulators` Depends on:
- [opm-common](https://github.com/OPM/opm-common)
- [opm-grid](https://github.com/OPM/opm-grid)
- [Dune and all other dependencies are listed and tracked on the home page](https://opm-project.org/?page_id=239)

[Follow the build instructions on the home page](http://opm-project.org/?page_id=36).


## In-Code Documentation
In addition to providing the [manual](https://opm-project.org/?page_id=955) we also document the source code with Doxygen, using the `make doc` command.


## Reporting Issues
Issues can either be reported here on [GitHub repository](https://github.com/OPM/opm-simulators/issues), or using the [OPM mailing list](https://opm-project.org/?page_id=358)

To help diagnose build errors, please provide a link to a build log together
with the issue description.

You can capture such a log from the build using the `script' utility, e.g.:

    LOGFILE=$(date +%Y%m%d-%H%M-)build.log ;
    cmake -E cmake_echo_color --cyan --bold "Log file: $LOGFILE" ;
    script -q $LOGFILE -c 'cmake ../opm-core -DCMAKE_BUILD_TYPE=Debug' &&
    script -q $LOGFILE -a -c 'ionice nice make -j 4 -l 3' ||
    cat CMakeCache.txt CMakeFiles/CMake*.log >> $LOGFILE

The resulting file can be uploaded to for instance gist.github.com.


## Citing
To cite OPM Flow we primarily use the paper `The Open Porous Media Flow reservoir simulator` published in the `Computers & Mathematics with Applications` (CAMWA) Journal. The correctly attribute credit to later contributors, we also want the Zenodo repository containing the source code to be cited. The BibTex for both can be found below.
<details>
<summary>BibTeX for CAMWA article</summary>
  
```bibtex
@article{OPMFLOW,
 title = {The {Open} {Porous} {Media} {Flow} reservoir simulator},
 journal = {Computers \& Mathematics with Applications},
 volume = {81},
 pages = {159-185},
 year = {2021},
 note = {Development and Application of Open-source Software for Problems with Numerical PDEs},
 issn = {0898-1221},
 doi = {https://doi.org/10.1016/j.camwa.2020.05.014},
 url = {https://www.sciencedirect.com/science/article/pii/S0898122120302182},
 author = {Atgeirr Flø Rasmussen and Tor Harald Sandve and Kai Bao and Andreas Lauser and Joakim Hove and Bård Skaflestad and Robert Klöfkorn and Markus Blatt and Alf Birger Rustad and Ove Sævareid and Knut-Andreas Lie and Andreas Thune}
}
```
</details>

<details>
<summary> BibTeX for source code (Zenodo) </summary>

```bibtex
@software{ahmed_2025_15573878,
  author       = {Ahmed, Elyes and
                  Alvestad, Jostein and
                  BAO, KAI and
                  Baxendale, David and
                  Berge, Runar Lie and
                  Berland, Håvard and
                  Blatt, Markus and
                  Bowden, Josh and
                  Bueno, Jose Eduardo and
                  Chang, Justin and
                  Egberts, Paul and
                  Fuchs, Franz Georg and
                  Hægland, Håkon and
                  Hove, Joakim and
                  Kippe, Vegard and
                  Klöfkorn, Robert and
                  Krogstad, Stein and
                  Kvarving, Arne Morten and
                  Landa Marban, David and
                  Logstein, Jan Inge and
                  Lye, Kjetil Olsen and
                  Machado, Cintia Goncalves and
                  Marchiori, Giacomo and
                  Meyer Andersen, Tobias and
                  Mykkeltvedt, Trine and
                  Nane, Razvan and
                  Nebel, Lisa Julia and
                  Nilsen, Halvor Møll and
                  Qiu, Tong Dong and
                  Qiu, Tuoling and
                  Rasmussen, Atgeirr Flø and
                  Ritorto, Antonella and
                  Rustad, Alf Birger and
                  Sandve, Tor Harald and
                  Sævareid, Ove and
                  Skaflestad, Bård and
                  Skille, Torbjørn and
                  Tveit, Svenn and
                  Verveer, Peter and
                  Tóth, Michal and
                  Goodfield, Matthew and
                  Sæternes, Erik Hide},
  title        = {OPM Flow 2025.04},
  month        = jun,
  year         = 2025,
  publisher    = {Zenodo},
  version      = {2025.04},
  doi          = {10.5281/zenodo.15573878},
  url          = {https://doi.org/10.5281/zenodo.15573878},
}
```
</details>

We try to keep track of all publications having used OPM Flow in their scientific papers [here](https://opm-project.org/?page_id=39). Reach out to the contact point mentioned on the page if you want your article there.

## Contributing

Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines on how to contribute to this project.
