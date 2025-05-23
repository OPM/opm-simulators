import os
import sys
sys.path.insert(0, os.path.abspath('.'))

extensions = ['breathe']
breathe_projects = {
    'OPM': '../../doxygen/xml'
}
breathe_default_project = 'OPM'

html_theme = 'furo'
