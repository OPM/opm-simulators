import os
import unittest
from pathlib import Path
from .pytest_common import pushd, create_black_oil_simulator

class TestBasic(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        test_dir = Path(os.path.dirname(__file__))
        cls.data_dir = test_dir.parent.joinpath("test_data/SPE1CASE1a")

    def test_all(self):
        with pushd(self.data_dir):
            sim = create_black_oil_simulator(filename="SPE1CASE1.DATA")
            # NOTE: The following call should throw an exception since the simulation
            #   has not been initialized
            with self.assertRaises(RuntimeError):
                sim.get_dt()
