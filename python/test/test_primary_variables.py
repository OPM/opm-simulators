import os
import unittest
from pathlib import Path
from .pytest_common import pushd, create_black_oil_simulator, create_gas_water_simulator, create_onephase_simulator

class TestBasic(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # NOTE: See comment in test_basic.py for the reason why we are
        #   only using a single test_all() function instead of splitting
        #   it up in multiple test functions
        test_dir = Path(os.path.dirname(__file__))
        cls.data_dir_bo = test_dir.parent.joinpath("test_data/SPE1CASE1a")
        cls.data_dir_op = test_dir.parent.joinpath("test_data/SPE1CASE1")
        cls.data_dir_gw = test_dir.parent.joinpath("test_data/SPE1CASE2")

    # IMPORTANT: Since all the python unittests run in the same process we must be
    #  careful to not call MPI_Init() more than once.
    #  Tests are run alphabetically, so we need to make sure that
    #  the the first test calls MPI_Init(), therefore the name of the tests
    #  have a numeric label like "01" in test_01_bo to ensure that they
    #  are run in a given order.

    # IMPORTANT:This test must be run first since it calls MPI_Init()
    def test_01_blackoil(self):
        with pushd(self.data_dir_bo):
            sim = create_black_oil_simulator(filename="SPE1CASE1.DATA")
            sim.setup_mpi(True, False)
            sim.step_init()
            sim.step()
            pressure = sim.get_primary_variable(variable='pressure')
            self.assertAlmostEqual(pressure[0], 35795160.67, delta=1e4, msg='value of pressure')
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

    def test_02_onephase(self):
        with pushd(self.data_dir_op):
            sim = create_onephase_simulator("SPE1CASE1_WATER.DATA")
            sim.setup_mpi(False, False)
            sim.step_init()
            sim.step()
            pressure = sim.get_primary_variable(variable='pressure')
            self.assertAlmostEqual(pressure[0], 44780102.277570, delta=1e4, msg='value of pressure')
            pressure_meaning = sim.get_primary_variable_meaning(
                variable='pressure')
            pressure_meaning_map = sim.get_primary_variable_meaning_map(
                variable='pressure')
            self.assertEqual(pressure_meaning[0], pressure_meaning_map["Pw"])
            water_meaning = sim.get_primary_variable_meaning(
                variable='water')
            water_meaning_map = sim.get_primary_variable_meaning_map(
                variable='water')
            self.assertEqual(water_meaning[0], water_meaning_map["Disabled"])
            brine_meaning = sim.get_primary_variable_meaning(
                variable='brine')
            brine_meaning_map = sim.get_primary_variable_meaning_map(
                variable='brine')
            self.assertEqual(brine_meaning[0], brine_meaning_map["Disabled"])

    # IMPORTANT: This test must be run last since it calls MPI_Finalize()
    def test_99_gaswater(self):
        with pushd(self.data_dir_gw):
            sim = create_gas_water_simulator(filename="SPE1CASE2_GASWATER.DATA")
            sim.setup_mpi(False, True)
            sim.step_init()
            sim.step()
            pressure = sim.get_primary_variable(variable='pressure')
            self.assertAlmostEqual(pressure[0], 32410549.874824, delta=1e4, msg='value of pressure')
            pressure_meaning = sim.get_primary_variable_meaning(
                variable='pressure')
            pressure_meaning_map = sim.get_primary_variable_meaning_map(
                variable='pressure')
            self.assertEqual(pressure_meaning[0], pressure_meaning_map["Pg"])
            water_meaning = sim.get_primary_variable_meaning(
                variable='water')
            water_meaning_map = sim.get_primary_variable_meaning_map(
                variable='water')
            self.assertEqual(water_meaning[0], water_meaning_map["Sw"])
            gas_meaning = sim.get_primary_variable_meaning(
                variable='gas')
            gas_meaning_map = sim.get_primary_variable_meaning_map(
                variable='gas')
            self.assertEqual(gas_meaning[0], gas_meaning_map["Disabled"])
            brine_meaning = sim.get_primary_variable_meaning(
                variable='brine')
            brine_meaning_map = sim.get_primary_variable_meaning_map(
                variable='brine')
            self.assertEqual(brine_meaning[0], brine_meaning_map["Disabled"])
