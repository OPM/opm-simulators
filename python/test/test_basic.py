import os
import unittest
from pathlib import Path
from opm.simulators import BlackOilSimulator, GasWaterSimulator
from .pytest_common import pushd

class TestBasic(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        test_dir = Path(os.path.dirname(__file__))
        cls.data_dir_bo = test_dir.parent.joinpath("test_data/SPE1CASE1a")
        cls.data_dir_gw = test_dir.parent.joinpath("test_data/SPE1CASE2")

    # IMPORTANT: Since all the python unittests run in the same process we must be
    #  careful to not call MPI_Init() more than once.
    #  Tests are run alphabetically, so we need to make sure that
    #  the the first test calls MPI_Init(), therefore the name of the tests
    #  have a numeric label like "01" in test_01_blackoil to ensure that they
    #  are run in a given order.

    # IMPORTANT: This test must be run first since it calls MPI_Init()
    def test_01_blackoil(self):
        with pushd(self.data_dir_bo):
            sim = BlackOilSimulator(args=['--linear-solver=ilu0'], filename="SPE1CASE1.DATA")
            sim.setup_mpi(init=True, finalize=False)
            sim.step_init()
            sim.step()
            dt = sim.get_dt()
            # NOTE: the timestep size is reduced to 1 day to avoid regression failures
            # due to changes in time stepping.
            # NOTE: The timestep should be 1 day
            #  = 24 * 60 * 60 seconds = 86400 seconds
            self.assertAlmostEqual(dt, 86400., places=7, msg='value of timestep')
            vol = sim.get_cell_volumes()
            self.assertEqual(len(vol), 300, 'length of volume vector')
            # NOTE: The volume should be 1000 ft x 1000 ft x 20 ft * 0.3 (porosity) 
            #  = 600000 ft^3 = 566336.93 m^3
            self.assertAlmostEqual(vol[0], 566336.93, places=2, msg='value of volume')
            poro = sim.get_porosity()
            self.assertEqual(len(poro), 300, 'length of porosity vector')
            self.assertAlmostEqual(poro[0], 0.3, places=7, msg='value of porosity')
            poro = poro *.95
            sim.set_porosity(poro)
            sim.step()
            poro2 = sim.get_porosity()
            self.assertAlmostEqual(poro2[0], 0.285, places=7, msg='value of porosity 2')

    # IMPORTANT: This test must be run last since it calls MPI_Finalize()
    def test_99_gaswater(self):
        with pushd(self.data_dir_gw):
            sim = GasWaterSimulator(args=['--linear-solver=ilu0'], filename="SPE1CASE2_GASWATER.DATA")
            sim.setup_mpi(init=False, finalize=True)
            sim.step_init()
            sim.step()
            dt = sim.get_dt()  # 31 days = 31 * 24 * 60 * 60 = 2678400 seconds
            self.assertAlmostEqual(dt, 2678400., places=7, msg='value of timestep')
            vol = sim.get_cell_volumes()
            self.assertEqual(len(vol), 300, 'length of volume vector')
            self.assertAlmostEqual(vol[0], 566336.93, places=2, msg='value of volume')
            poro = sim.get_porosity()
            self.assertEqual(len(poro), 300, 'length of porosity vector')
            self.assertAlmostEqual(poro[0], 0.3, places=7, msg='value of porosity')
            poro = poro *.95
            sim.set_porosity(poro)
            sim.step()
            poro2 = sim.get_porosity()
            self.assertAlmostEqual(poro2[0], 0.285, places=7, msg='value of porosity 2')



