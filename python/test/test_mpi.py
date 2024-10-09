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

    # IMPORTANT: Tests are run alphabetically, so we need to make sure that
    #  the the first test calls MPI_Init(), therefore the name of the tests
    #  have a numeric label like "01" in test_01_mpi_init to ensure that they
    #  are run in a given order.
    def test_01_mpi_init(self):
        with pushd(self.data_dir):
            sim = BlackOilSimulator("SPE1CASE1.DATA")
            sim.setup_mpi(init=True, finalize=False)
            sim.step_init()  # This will create the OPM::Main() object which will call MPI_Init()
            assert True

    def test_02_mpi_no_init(self):
        with pushd(self.data_dir):
            sim = BlackOilSimulator("SPE1CASE1.DATA")
            sim.setup_mpi(init=False, finalize=True)
            sim.step_init()  # This will create the OPM::Main() object which will not call MPI_Init()
            # That this test runs shows that the simulator does not call
            # MPI_Init() a second time
            assert True
