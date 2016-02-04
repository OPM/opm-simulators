/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.
  Copyright 2015 Statoil ASA.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <config.h>

#if HAVE_DYNAMIC_BOOST_TEST
#define BOOST_TEST_DYN_LINK
#endif

#define BOOST_TEST_MODULE RestartTests

#include <iostream>
#include <map>

#include <ert/ecl/ecl_rst_file.h>
#include <ert/ecl/ecl_kw.h>

#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_CASE(CompareRestartFileResults)
{
    const std::string& filename1 = boost::unit_test::framework::master_test_suite().argv[1];
    const std::string& filename2 = boost::unit_test::framework::master_test_suite().argv[2];
    int last_report_step = std::atoi(boost::unit_test::framework::master_test_suite().argv[3]);

    std::map<std::string, double> relative_diffs;
    relative_diffs["SWAT"]     = 0.000200;  //0.02 %
    relative_diffs["SGAS"]     = 0.000200;

    relative_diffs["RS"]       = 0.000010;  //0.001 %
    relative_diffs["RV"]       = 0.000010;
    relative_diffs["PRESSURE"] = 0.000010;

    ecl_file_type*  file1 =  ecl_file_open_rstblock_report_step( filename1.c_str() , last_report_step, 1);
    ecl_file_type*  file2 =  ecl_file_open_rstblock_report_step( filename2.c_str() , last_report_step, 1);

    for (auto key : {"PRESSURE", "SWAT", "SGAS", "RS", "RV"}) {
        ecl_kw_type * kw_1 =  ecl_file_iget_named_kw( file1 , key, 0);
        ecl_kw_type * kw_2 =  ecl_file_iget_named_kw( file2 , key, 0);

        bool numeric_equal = ecl_kw_numeric_equal(kw_1, kw_2, relative_diffs[key]);
        if (numeric_equal) {
            std::cout << " Restart results for " << key << " compared ok" << std::endl;
        } else {
            std::cout << " Restart results for " << key << " not ok, failing test " << std::endl;
        }

        BOOST_CHECK_EQUAL(numeric_equal, true);
    }

    ecl_file_close(file1);
    ecl_file_close(file2);
}


