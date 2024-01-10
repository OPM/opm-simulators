import os
import unittest
from pathlib import Path
from opm.simulators import BlackOilSimulator
from .pytest_common import pushd

class TestBasic(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # NOTE: See comment in test_basic.py for the reason why we are
        #   only using a single test_all() function instead of splitting
        #   it up in multiple test functions
        test_dir = Path(os.path.dirname(__file__))
        cls.data_dir = test_dir.parent.joinpath("test_data/SPE1CASE1a")


    def test_all(self):
        with pushd(self.data_dir):
            sim = BlackOilSimulator("SPE1CASE1.DATA")
            sim.step_init()
            sim.step()
            #oil_pressure = sim.get_fluidstate_variable(name='po')
            #self.assertAlmostEqual(oil_pressure[0], 41729978.837, places=2, msg='value of oil pressure')
