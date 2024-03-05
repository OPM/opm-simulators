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
            oil_pressure = sim.get_fluidstate_variable(name='po')
            self.assertAlmostEqual(oil_pressure[0], 41729978.837, places=2, msg='value of oil pressure')
            gas_pressure = sim.get_fluidstate_variable(name='pg')
            self.assertAlmostEqual(gas_pressure[0], 41729978.837, places=2, msg='value of gas pressure')
            water_pressure = sim.get_fluidstate_variable(name='pw')
            self.assertAlmostEqual(water_pressure[0], 41729978.837, places=2, msg='value of water pressure')
            rho_w = sim.get_fluidstate_variable(name='rho_w')
            self.assertAlmostEqual(rho_w[0], 1001.7549054, places=6, msg='value of water density')
            rho_g = sim.get_fluidstate_variable(name='rho_g')
            self.assertAlmostEqual(rho_g[0], 275.72397867, places=7, msg='value of gas density')
            rho_o = sim.get_fluidstate_variable(name='rho_o')
            self.assertAlmostEqual(rho_o[0], 639.64061021, places=7, msg='value of oil density')
            Rs = sim.get_fluidstate_variable(name='Rs')
            self.assertAlmostEqual(Rs[0], 226.196660482, places=7, msg='value of solution gas-oil ratio')
            Rv = sim.get_fluidstate_variable(name='Rv')
            self.assertAlmostEqual(Rv[0], 0.0, places=7, msg='value of volatile gas-oil ratio')
            Sw = sim.get_fluidstate_variable(name='Sw')
            self.assertAlmostEqual(Sw[0], 0.11905577997, places=10, msg='value of water saturation')
            So = sim.get_fluidstate_variable(name='So')
            self.assertAlmostEqual(So[0], 0.56951652831, places=10, msg='value of oil saturation')
            Sg = sim.get_fluidstate_variable(name='Sg')
            self.assertAlmostEqual(Sg[0], 0.31142769170, places=10, msg='value of gas saturation')
            T = sim.get_fluidstate_variable(name='T')
            self.assertAlmostEqual(T[0], 288.705, places=3, msg='value of temperature')
