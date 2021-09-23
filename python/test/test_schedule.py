import os
import unittest
from contextlib import contextmanager
import datetime as dt
from pathlib import Path
from opm2.simulators import BlackOilSimulator
from opm.io.parser import Parser
from opm.io.ecl_state import EclipseState
from opm.io.schedule import Schedule
from opm.io.summary import SummaryConfig

@contextmanager
def pushd(path):
    cwd = os.getcwd()
    if not os.path.isdir(path):
        os.makedirs(path)
    os.chdir(path)
    yield
    os.chdir(cwd)


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
            deck  = Parser().parse('SPE1CASE1.DATA')
            state = EclipseState(deck)
            schedule = Schedule( deck, state )
            summary_config = SummaryConfig(deck, state, schedule)
            self.assertTrue('PROD' in schedule)
            self.assertTrue('INJ'  in schedule)
            self.assertEqual(dt.datetime(2015, 1, 1),   schedule.start)
            self.assertEqual(dt.datetime(2016, 1, 1), schedule.end)
            sim = BlackOilSimulator( deck, state, schedule, summary_config  )
            sim.step_init()
            sim.step()
            prod = schedule.get_well("PROD", 2)
            self.assertEqual(prod.status(), "OPEN")
            #schedule.shut_well("PROD", 3)
            #prod = schedule.get_well("PROD", 3)
            #self.assertEqual(prod.status(), "SHUT")
            sim.step()
            sim.step()

