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
            # NOTE: The following call should throw an exception since the simulation
            #   has not been initialized
            with self.assertRaises(RuntimeError):
                sim.get_dt()

