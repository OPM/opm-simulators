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

    # IMPORTANT: Since all the python unittests run in the same process we must be
    #  careful to not call MPI_Init() more than once.
    #  Tests are run alphabetically, so we need to make sure that
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
