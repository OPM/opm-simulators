import os
import shutil
import unittest
from pathlib import Path
from opm.io.ecl import ESmry
from opm.io.parser import Parser
from opm.io.ecl_state import EclipseState
from opm.io.schedule import Schedule
from opm.io.summary import SummaryConfig
from .pytest_common import pushd, create_black_oil_simulator

class TestGetSchedule(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        test_dir = Path(os.path.dirname(__file__))
        cls.data_dir = test_dir.parent.joinpath("test_data/SPE1CASE1b")

    # IMPORTANT: Since all the tests in this file run in the same process we must be
    #  careful to not call MPI_Init() more than once. Tests are run alphabetically,
    #  so the numeric labels make sure that the first test to call step_init()
    #  is the one that calls MPI_Init(), and the last one calls MPI_Finalize().
    def test_01_before_step_init(self):
        with pushd(self.data_dir):
            sim = create_black_oil_simulator(filename="SPE1CASE1.DATA")
            with self.assertRaises(RuntimeError):
                sim.get_schedule()

    def test_02_filename_constructor(self):
        # A well shut through get_schedule() must be shut in the simulation itself
        output_dir = "02_get_schedule"
        with pushd(self.data_dir):
            sim = create_black_oil_simulator(
                filename="SPE1CASE1.DATA", args=[f"--output-dir={output_dir}"])
            sim.setup_mpi(init=True, finalize=False)
            sim.step_init()
            schedule = sim.get_schedule()
            self.assertEqual(schedule.get_well("PROD", 3).status(), "OPEN")
            schedule.shut_well("PROD", 3)
            self.assertEqual(schedule.get_well("PROD", 3).status(), "SHUT")
            # The change must be visible through a second call too
            self.assertEqual(sim.get_schedule().get_well("PROD", 3).status(), "SHUT")
            sim.advance(report_step=5)
            sim.step_cleanup()
            fopr = ESmry(f"{output_dir}/SPE1CASE1.SMSPEC")["FOPR", True]
            # FOPR at the end of report steps 1 to 5. PROD is the only producer.
            self.assertGreater(fopr[2], 0.0)
            self.assertEqual(fopr[3], 0.0)
            self.assertEqual(fopr[4], 0.0)

    def test_03_objects_constructor(self):
        # The simulator must return the Schedule object it was given, not a copy.
        # Run from a copy of the deck in a directory of its own: with an
        # EclipseState given, the simulator clears old result files and writes
        # the .PRT and .DBG files next to the deck the EclipseState was built
        # from, whatever --output-dir says, and test_schedule.py uses this deck.
        run_dir = self.data_dir / "03_get_schedule"
        run_dir.mkdir(exist_ok=True)
        shutil.copy(self.data_dir / "SPE1CASE1.DATA", run_dir)
        with pushd(run_dir):
            deck = Parser().parse("SPE1CASE1.DATA")
            state = EclipseState(deck)
            schedule = Schedule(deck, state)
            summary_config = SummaryConfig(deck, state, schedule)
            sim = create_black_oil_simulator(deck, state, schedule, summary_config)
            sim.setup_mpi(init=False, finalize=True)
            sim.step_init()
            self.assertIs(sim.get_schedule(), schedule)
            # NOTE: step_cleanup() needs at least one step() after step_init()
            sim.step()
            sim.step_cleanup()
