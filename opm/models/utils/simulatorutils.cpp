// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/

#include <config.h>
#include <opm/models/utils/simulatorutils.hpp>

#include <opm/models/discretization/common/fvbaseparameters.hh>

#include <opm/models/utils/parametersystem.hpp>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

#if HAVE_QUAD
#include <opm/material/common/quad.hpp>
#endif

namespace Opm {

std::string humanReadableTime(double timeInSeconds, bool isAmendment)
{
    std::ostringstream oss;
    oss << std::setprecision(4);
    if (isAmendment) {
        oss << " (";
    }
    if (timeInSeconds >= 365.25 * 24 * 60 * 60) {
        int years = static_cast<int>(timeInSeconds / (365.25 * 24 * 60 * 60));
        int days = static_cast<int>((timeInSeconds - years*(365.25 * 24 * 60 * 60)) / (24 * 60 * 60));

        constexpr double accuracy = 1e-2;
        double hours =
            std::round(1.0 / accuracy *
                       (timeInSeconds
                        - years * (365.25 * 24 * 60 * 60)
                        - days * (24 * 60 * 60)) / (60 * 60))
            *accuracy;

        oss << years << " years, " << days << " days, "  << hours << " hours";
    }
    else if (timeInSeconds >= 24.0 * 60 * 60) {
        int days = static_cast<int>(timeInSeconds / (24 * 60 * 60));
        int hours = static_cast<int>((timeInSeconds - days * (24 * 60 * 60)) / (60 * 60));

        constexpr double accuracy = 1e-2;
        double minutes =
            std::round(1.0 / accuracy *
                       (timeInSeconds
                        - days * (24 * 60 * 60)
                        - hours * (60 * 60)) / 60)
            *accuracy;

        oss << days << " days, " << hours << " hours, " << minutes << " minutes";
    }
    else if (timeInSeconds >= 60.0 * 60) {
        int hours = static_cast<int>(timeInSeconds / (60 * 60));
        int minutes = static_cast<int>((timeInSeconds - hours * (60 * 60)) / 60);

        constexpr double accuracy = 1e-2;
        double seconds =
            std::round(1.0 / accuracy *
                       (timeInSeconds
                        - hours * (60 * 60)
                        - minutes * 60))
            * accuracy;

        oss << hours << " hours, " << minutes << " minutes, "  << seconds << " seconds";
    }
    else if (timeInSeconds >= 60.0) {
        int minutes = static_cast<int>(timeInSeconds / 60);

        constexpr double accuracy = 1e-3;
        double seconds =
            std::round(1.0 / accuracy *
                       (timeInSeconds
                        - minutes * 60))
            * accuracy;

        oss << minutes << " minutes, "  << seconds << " seconds";
    }
    else if (!isAmendment) {
        oss << timeInSeconds << " seconds";
    }
    else {
        return "";
    }
    if (isAmendment) {
        oss << ")";
    }

    return oss.str();
}

std::string simulatorOutputDir()
{
    std::string outputDir = Parameters::Get<Parameters::OutputDir>();

    if (outputDir.empty()) {
        outputDir = ".";
    }

    // TODO: replace this by std::filesystem once we require c++-2017
    struct stat st;
    if (::stat(outputDir.c_str(), &st) != 0)
        throw std::runtime_error("Could not access output directory '" + outputDir + "':" +
                                 strerror(errno));
    if (!S_ISDIR(st.st_mode)) {
        throw std::runtime_error("Path to output directory '" +outputDir +
                                 "' exists but is not a directory");
    }

    if (access(outputDir.c_str(), W_OK) != 0) {
        throw std::runtime_error("Output directory '" + outputDir +
                                 "' exists but is not writeable");
    }

    return outputDir;
}

template<class Scalar>
std::vector<Scalar> readTimeStepFile(const std::string& file)
{
    std::ifstream is(file);
    std::vector<Scalar> result;
    while (!is.eof()) {
        Scalar dt;
        is >> dt;
        result.push_back(dt);
    }

    return result;
}

template std::vector<double> readTimeStepFile(const std::string&);
template std::vector<float> readTimeStepFile(const std::string&);

#if HAVE_QUAD
template std::vector<quad> readTimeStepFile(const std::string&);
#endif


} // end namespace Opm
