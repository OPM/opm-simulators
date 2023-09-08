import os
import unittest
from contextlib import contextmanager
import datetime as dt
from pathlib import Path
import re
from opm.simulators import BlackOilSimulator
from opm.io.parser import Parser
from opm.io.ecl_state import EclipseState
from opm.io.schedule import Schedule
from opm.io.summary import SummaryConfig
from .pytest_common import pushd

class TestBasic(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # NOTE: See comment in test_basic.py for the reason why we are
        #   only using a single test_all() function instead of splitting
        #   it up in multiple test functions
        test_dir = Path(os.path.dirname(__file__))
        cls.data_dir = test_dir.parent.joinpath("test_data/SPE1CASE1b")


    def test_all(self):
        with pushd(self.data_dir):
            self.deck  = Parser().parse('SPE1CASE1.DATA')
            state = EclipseState(self.deck)
            self.schedule = Schedule( self.deck, state )
            summary_config = SummaryConfig(self.deck, state, self.schedule)
            self.unit_system = self.deck.active_unit_system()

            self.assertTrue('PROD' in self.schedule)
            self.assertTrue('INJ'  in self.schedule)
            self.assertEqual(dt.datetime(2015, 1, 1),   self.schedule.start)
            self.assertEqual(dt.datetime(2016, 1, 1), self.schedule.end)
            self.sim = BlackOilSimulator(
                self.deck, state, self.schedule, summary_config  )
            tsteps = self.schedule.timesteps
            self.assertEqual(dt.datetime(2015, 1, 1), tsteps[0])
            last_step = len(tsteps) - 1
            self.assertEqual(dt.datetime(2016, 1, 1), tsteps[last_step])
            self.sim.step_init()
            report_step = 4
            self.sim.advance(report_step=report_step)
            well_name = "PROD"
            prod = self.schedule.get_well(well_name, 2)
            self.assertEqual(prod.status(), "OPEN")
            #schedule.shut_well("PROD", 3)
            #prod = schedule.get_well("PROD", 3)
            #self.assertEqual(prod.status(), "SHUT")
            self.subtest_modify_prod_weltarg_dynamically(well_name, report_step)
            self.sim.step()
            report_step = self.sim.current_step()
            well_name = "INJ"
            self.subtest_modify_inj_weltarg_dynamically(well_name, report_step)
            self.sim.advance(report_step=last_step)
            self.sim.step_cleanup()

    def subtest_modify_inj_weltarg_dynamically(self, well_name, report_step):
        prop = self.schedule.get_injection_properties(well_name, report_step)
        self.assertEqual(prop['surf_inj_rate'], 100000.0) # Mscf/day
        self.assertEqual(prop['resv_inj_rate'], 0.0) # rb/day
        self.assertEqual(prop['bhp_target'], 9014.0) # psi
        self.assertEqual(prop['thp_target'], 0.0)
        new_grat_target = prop['surf_inj_rate'] - 100  # stb/day
        self.update_inj_grat_target_wconinje(well_name, new_grat_target)
        self.sim.step()
        prop2 = self.schedule.get_injection_properties(well_name, report_step+1)
        self.assertEqual(prop2['surf_inj_rate'], new_grat_target)
        new_grat_target += 200
        self.update_inj_grat_target_weltarg(well_name, new_grat_target)
        self.sim.step()
        prop3 = self.schedule.get_injection_properties(well_name, report_step+2)
        self.assertEqual(prop3['surf_inj_rate'], new_grat_target)

    def subtest_modify_prod_weltarg_dynamically(self, well_name, report_step):
        prop = self.schedule.get_production_properties(well_name, report_step)
        self.assertEqual(prop['alq_value'], 0.0)
        self.assertEqual(prop['bhp_target'], 1000.0)
        self.assertEqual(prop['gas_rate'], 0.0)
        self.assertEqual(prop['liquid_rate'], 0.0)
        self.assertEqual(prop['oil_rate'], 20000.0)
        self.assertEqual(prop['resv_rate'], 0.0)
        self.assertEqual(prop['thp_target'], 0.0)
        self.assertEqual(prop['water_rate'], 0.0)
        new_oil_target = prop['oil_rate'] + 10000  # stb/day
        self.update_prod_orat_target_wconprod(well_name, new_oil_target)
        self.sim.step()
        prop2 = self.schedule.get_production_properties(well_name, report_step+1)
        self.assertEqual(prop2['oil_rate'], new_oil_target)
        new_oil_target += 1000
        self.update_prod_orat_target_weltarg(well_name, new_oil_target)
        self.sim.step()
        prop3 = self.schedule.get_production_properties(well_name, report_step+2)
        self.assertEqual(prop3['oil_rate'], new_oil_target)

    # This is an alternative to using WELTARG
    def update_inj_grat_target_wconinje(self, well_name, new_surf_flow_rate):
        data = self.deck["WCONINJE"]
        # assumes data looks like this:
        #  WCONINJE
        #  	'INJ'	'GAS'	'OPEN'	'RATE'	100000 1* 9014 /
        #  /
        # The initial rate can also be obtained from data[0][4].get_uda(0).get_double()
        data = re.sub(pattern='100000', repl=str(new_surf_flow_rate), string=str(data), count=1)
        report_step = self.sim.current_step()
        self.schedule.insert_keywords(
            data, step=report_step, unit_system=self.unit_system)

    # This is an alternative to using WCONINJE to modify injection properties
    def update_inj_grat_target_weltarg(self, well_name, net_surf_flow_rate):
        data = """
WELTARG
    '{}'  GRAT {} /
/
        """.format(well_name, net_surf_flow_rate)
        report_step = self.sim.current_step()
        self.schedule.insert_keywords(
            data, step=report_step, unit_system=self.unit_system)


    # This is an alternative to using WCONPROD to modify production properties
    def update_prod_orat_target_weltarg(self, well_name, oil_target):
        data = """
WELTARG
    '{}'  ORAT {} /
/
        """.format(well_name, oil_target)
        report_step = self.sim.current_step()
        self.schedule.insert_keywords(
            data, step=report_step, unit_system=self.unit_system)

    # This is an alternative to using WELTARG to modify production properties
    def update_prod_orat_target_wconprod(self, well_name, oil_target):
        well_status = "OPEN"
        control_mode = "ORAT"
        bhp_limit = 1000  # psia
        data = """
WCONPROD
	'{}' '{}' '{}' {} 4* {} /
/
        """.format(well_name, well_status, control_mode, oil_target, bhp_limit)
        report_step = self.sim.current_step()
        self.schedule.insert_keywords(
            data, step=report_step, unit_system=self.unit_system)
