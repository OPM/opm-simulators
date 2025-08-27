import os
import unittest
from pathlib import Path
from opm.simulators import BlackOilSimulator, GasWaterSimulator
from .pytest_common import pushd

class TestBasic(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # NOTE: See comment in test_basic.py for the reason why we are
        #   only using a single test_all() function instead of splitting
        #   it up in multiple test functions
        test_dir = Path(os.path.dirname(__file__))
        cls.data_dir_bo = test_dir.parent.joinpath("test_data/SPE1CASE1a")
        cls.data_dir_gw = test_dir.parent.joinpath("test_data/SPE1CASE2")


    def test_all(self):
        self._bo()
        self._gw()

    def _bo(self):
        with pushd(self.data_dir_bo):
            sim = BlackOilSimulator("SPE1CASE1.DATA")
            sim.setup_mpi(True, False)
            sim.step_init()
            sim.step()
            oil_pressure = sim.get_fluidstate_variable(name='po')
            self.assertAlmostEqual(oil_pressure[0], 35795160.67, delta=1e4, msg='value of oil pressure')
            gas_pressure = sim.get_fluidstate_variable(name='pg')
            self.assertAlmostEqual(gas_pressure[0], 35795160.67, delta=1e4, msg='value of gas pressure')
            water_pressure = sim.get_fluidstate_variable(name='pw')
            self.assertAlmostEqual(water_pressure[0], 35795160.67, delta=1e4, msg='value of water pressure')
            rho_w = sim.get_fluidstate_variable(name='rho_w')
            self.assertAlmostEqual(rho_w[0], 998.9822355, places=3, msg='value of water density')
            rho_g = sim.get_fluidstate_variable(name='rho_g')
            self.assertAlmostEqual(rho_g[0], 241.36955087, places=2, msg='value of gas density')
            rho_o = sim.get_fluidstate_variable(name='rho_o')
            self.assertAlmostEqual(rho_o[0], 631.78674, places=2, msg='value of oil density')
            Rs = sim.get_fluidstate_variable(name='Rs')
            self.assertAlmostEqual(Rs[0], 226.19666048, places=5, msg='value of solution gas-oil ratio')
            Rv = sim.get_fluidstate_variable(name='Rv')
            self.assertAlmostEqual(Rv[0], 0.0, places=5, msg='value of volatile gas-oil ratio')
            Sw = sim.get_fluidstate_variable(name='Sw')
            self.assertAlmostEqual(Sw[0], 0.11969486712, places=5, msg='value of water saturation')
            So = sim.get_fluidstate_variable(name='So')
            self.assertAlmostEqual(So[0], 0.825129, places=3, msg='value of oil saturation')
            Sg = sim.get_fluidstate_variable(name='Sg')
            self.assertAlmostEqual(Sg[0], 0.055138968544, places=3, msg='value of gas saturation')
            T = sim.get_fluidstate_variable(name='T')
            self.assertAlmostEqual(T[0], 288.705, places=3, msg='value of temperature')

    def _gw(self):
        with pushd(self.data_dir_gw):
            sim = GasWaterSimulator("SPE1CASE2_GASWATER.DATA")
            sim.setup_mpi(False, True)
            sim.step_init()
            sim.step()
            gas_pressure = sim.get_fluidstate_variable(name='pg')
            self.assertAlmostEqual(gas_pressure[0], 32410549.874824, delta=1e4, msg='value of gas pressure')
            water_pressure = sim.get_fluidstate_variable(name='pw')
            self.assertAlmostEqual(water_pressure[0], 32410549.874824, delta=1e4, msg='value of water pressure')
            rho_g = sim.get_fluidstate_variable(name='rho_g')
            self.assertAlmostEqual(rho_g[0], 219.613424, places=2, msg='value of gas density')
            rho_w = sim.get_fluidstate_variable(name='rho_w')
            self.assertAlmostEqual(rho_w[0], 1047.849235, places=3, msg='value of water density')
            Sg = sim.get_fluidstate_variable(name='Sg')
            self.assertAlmostEqual(Sg[0], 0.708648, places=3, msg='value of gas saturation')
            Sw = sim.get_fluidstate_variable(name='Sw')
            self.assertAlmostEqual(Sw[0], 0.291351, places=5, msg='value of water saturation')