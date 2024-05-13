# Configuration file for the Sphinx documentation builder.

project = "opm.simulators"
copyright = "2024 Equinor ASA"
author = "Håkon Hægland"

# Function to extract release version from dune.module file
def extract_opm_simulators_release():
    version_file_path = '../../../dune.module'
    with open(version_file_path, 'r') as file:
        for line in file:
            if line.startswith('Version:'):
                version_string = line.split(':')[1].strip()
                return version_string
    return "unknown"  # Fallback version

release = extract_opm_simulators_release()

# -- General configuration ---------------------------------------------------
import os
import sys

# For regular Python packages, the path to the package is usually added to sys.path
# here such that autodoc can find the modules. However, our Python module
#  opm.simulators.BlackOilSimulator is not generated yet. Since it is a pybind11
# extension module, it needs to be compiled as part of the opm-simulators build
# process (which requires building opm-common first and so on). The full compilation
# of opm-simulators requres time and storage resources, so we do not want to
# do this as part of the documentation build. Instead, we will use a sphinx extension
# (Python script) to generate documentation from a JSON input file. (This JSON file
# is also is also used to generate a C++ header file with docstrings. The header file
# is generated when opm-simulators is built and is then included
# in the pybind11 binding code for the BlackOilSimulator class. In this way, the
# opm.simulators.BlackOilSimulator Python module will also have access to the docstrings.)
sys.path.insert(0, os.path.abspath("../src"))
# Our sphinx extension that will use the docstrings.json file to generate documentation
extensions = ["opm_simulators_docs.sphinx_ext_docstrings"]
# Path to docstrings.json
opm_simulators_docstrings_path = os.path.abspath('../../docstrings_simulators.json')
opm_common_docstrings_path = os.path.abspath('../../docstrings_common.json')

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------
# See: https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output
html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_context = {
    "display_github": True,
    "github_user": "OPM",
    "github_repo": "opm-simulators",
    "github_version": "master",
    "conf_py_path": "/python/docs/",
}
