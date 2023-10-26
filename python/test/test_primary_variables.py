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
            pressure = sim.get_primary_variable(variable='pressure')
            self.assertAlmostEqual(pressure[0], 41729978.837, places=2, msg='value of pressure')
            pressure_meaning = sim.get_primary_variable_meaning(
                variable='pressure')
            pressure_meaning_map = sim.get_primary_variable_meaning_map(
                variable='pressure')
            self.assertEqual(pressure_meaning[0], pressure_meaning_map["Po"])
            water_meaning = sim.get_primary_variable_meaning(
                variable='water')
            water_meaning_map = sim.get_primary_variable_meaning_map(
                variable='water')
            self.assertEqual(water_meaning[0], water_meaning_map["Sw"])
            gas_meaning = sim.get_primary_variable_meaning(
                variable='gas')
            gas_meaning_map = sim.get_primary_variable_meaning_map(
                variable='gas')
            self.assertEqual(gas_meaning[0], gas_meaning_map["Sg"])
            brine_meaning = sim.get_primary_variable_meaning(
                variable='brine')
            brine_meaning_map = sim.get_primary_variable_meaning_map(
                variable='brine')
            self.assertEqual(brine_meaning[0], brine_meaning_map["Disabled"])
