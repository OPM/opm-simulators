import os
import unittest
from pathlib import Path
from opm.simulators import BlackOilSimulator
from .pytest_common import pushd

class TestBasic(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        test_dir = Path(os.path.dirname(__file__))
        cls.data_dir = test_dir.parent.joinpath("test_data/SPE1CASE1a")


    def test_all(self):
        with pushd(self.data_dir):
            sim = BlackOilSimulator("SPE1CASE1.DATA")
            sim.step_init()
            sim.step()
            dt = sim.get_dt()
            # NOTE: The timestep should be 1 month = 31 days
            #  = 31 * 24 * 60 * 60 seconds = 2678400 seconds
            self.assertAlmostEqual(dt, 2678400.0, places=7, msg='value of timestep')
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


