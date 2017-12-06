
# This script generate the illustration pictures for the documentation.
#
# To run this script, you have to install paraview, see:
#
# http://www.paraview.org/paraview/resources/software.php
#
# Eventually, set up the paths (figure_path, tutorial_data_path) according to your own installation.
# (The default values should be ok.)
#
# Make sure that pvpython is in your path of executables.
#
# After all the tutorial programs have been executed, run the following
# command in the same directory:
#
#   pvpython generate_doc_figures.py Documentation/Figure
#

from paraview.simple import *
from os import remove, mkdir, curdir
from os.path import join, isdir
from sys import argv, exit

# we need at least the output directory
if len(argv) <= 1:
	exit('Synopsis: pvpython generate_doc_figures.py dest-dir [src-dir]')

figure_path = argv[1]

# default for the input directory is the current one
if len(argv) <= 2:
    tutorial_data_path = curdir
else:
    tutorial_data_path = argv[2]

collected_garbage_file = []

if not isdir(figure_path):
    mkdir(figure_path)
    
## [tutorial1]
# tutorial 1
data_file_name = join(tutorial_data_path, "tutorial1.vtu")
# grid = servermanager.sources.XMLUnstructuredGridReader(FileName = data_file_name)
grid = XMLUnstructuredGridReader(FileName = data_file_name)
collected_garbage_file.append(data_file_name)
grid.UpdatePipeline()
Show(grid)
dp = GetDisplayProperties(grid)
dp.Representation = 'Wireframe'
dp.LineWidth = 5
dp.AmbientColor = [1, 0, 0]
view = GetActiveView()
view.Background = [1, 1, 1]
camera = GetActiveCamera()
camera.SetPosition(4, -6, 5)
camera.SetViewUp(-0.19, 0.4, 0.9)
camera.SetViewAngle(30)
camera.SetFocalPoint(1.5, 1.5, 1)
Render()
WriteImage(join(figure_path, "tutorial1.png"))
Hide(grid)
## [tutorial1]

## [tutorial2]
# tutorial 2
data_file_name = join(tutorial_data_path, "tutorial2.vtu")
grid = XMLUnstructuredGridReader(FileName = data_file_name)
collected_garbage_file.append(data_file_name)
grid.UpdatePipeline()
Show(grid)
dp = GetDisplayProperties(grid)
dp.Representation = 'Surface'
dp.ColorArrayName = 'pressure'
pres = grid.CellData.GetArray(0)
pres_lookuptable = GetLookupTableForArray( "pressure", 1, RGBPoints=[pres.GetRange()[0], 1, 0, 0, pres.GetRange()[1], 0, 0, 1] )
dp.LookupTable = pres_lookuptable
view = GetActiveView()
view.Background = [1, 1, 1]
camera = GetActiveCamera()
camera.SetPosition(20, 20, 110)
camera.SetViewUp(0, 1, 0)
camera.SetViewAngle(30)
camera.SetFocalPoint(20, 20, 0.5)
Render()
WriteImage(join(figure_path, "tutorial2.png"))
Hide(grid)
## [tutorial2]

## [tutorial3]
# tutorial 3
for case in range(0,20):
    data_file_name = join(tutorial_data_path, "tutorial3-"+"%(case)03d"%{"case": case}+".vtu")
    collected_garbage_file.append(data_file_name)

cases = ["000", "005", "010", "015", "019"]
for case in cases:
    data_file_name = join(tutorial_data_path, "tutorial3-"+case+".vtu")
    grid = XMLUnstructuredGridReader(FileName = data_file_name)
    grid.UpdatePipeline()
    Show(grid)
    dp = GetDisplayProperties(grid)
    dp.Representation = 'Surface'
    dp.ColorArrayName = 'saturation'
    sat = grid.CellData.GetArray(1)
    sat_lookuptable = GetLookupTableForArray( "saturation", 1, RGBPoints=[0, 1, 0, 0, 1, 0, 0, 1])
    dp.LookupTable = sat_lookuptable
    view.Background = [1, 1, 1]
    camera = GetActiveCamera()
    camera.SetPosition(100, 100, 550)
    camera.SetViewUp(0, 1, 0)
    camera.SetViewAngle(30)
    camera.SetFocalPoint(100, 100, 5)
    Render()
    WriteImage(join(figure_path, "tutorial3-"+case+".png"))
Hide(grid)
## [tutorial3]

## [tutorial4]
# tutorial 4
for case in range(0,20):
    data_file_name = join(tutorial_data_path, "tutorial4-"+"%(case)03d"%{"case": case}+".vtu")
    collected_garbage_file.append(data_file_name)

cases = ["000", "005", "010", "015", "019"]
for case in cases:
    data_file_name = join(tutorial_data_path, "tutorial4-"+case+".vtu")
    grid = XMLUnstructuredGridReader(FileName = data_file_name)
    grid.UpdatePipeline()
    Show(grid)
    dp = GetDisplayProperties(grid)
    dp.Representation = 'Surface'
    dp.ColorArrayName = 'saturation'
    sat = grid.CellData.GetArray(1)
    sat_lookuptable = GetLookupTableForArray( "saturation", 1, RGBPoints=[0, 1, 0, 0, 1, 0, 0, 1])
    dp.LookupTable = sat_lookuptable
    view.Background = [1, 1, 1]
    camera = GetActiveCamera()
    camera.SetPosition(100, 100, 550)
    camera.SetViewUp(0, 1, 0)
    camera.SetViewAngle(30)
    camera.SetFocalPoint(100, 100, 5)
    Render()
    WriteImage(join(figure_path, "tutorial4-"+case+".png"))
Hide(grid)
## [tutorial4]

# remove temporary files
for f in collected_garbage_file:
    remove(f)

