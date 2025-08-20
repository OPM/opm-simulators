import os
import time
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
            sim1 = BlackOilSimulator("SPE1CASE1.DATA")
            sim2 = BlackOilSimulator("SPE1CASE1.DATA")
            sim1.setup_mpi(init=True, finalize=False)
            sim2.setup_mpi(init=False, finalize=False)
            sim1.step_init()
            sim2.step_init()
            sim1.step()
            sim2.step()
            poro = sim1.get_porosity()
            sim1.set_porosity(poro * 0.95)
            sim1.step()
            sim2.step()
            poro12 = sim1.get_porosity()
            self.assertAlmostEqual(poro12[0], 0.285, places=7, msg='value of porosity 2')
            poro22 = sim2.get_porosity()
            self.assertAlmostEqual(poro22[0], 0.300, places=7, msg='value of porosity 2')
            sim1.step_cleanup()
            sim2.step_cleanup()



