/*
  Copyright 2020 Equinor ASA

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

/// This file includes the relevant OpenCL header(s)
/// All bda files using OpenCL declarations should include this header

#define CL_HPP_ENABLE_EXCEPTIONS
#define CL_TARGET_OPENCL_VERSION 120   // indicate OpenCL 1.2 is used
#define CL_HPP_TARGET_OPENCL_VERSION 120 // indicate OpenCL 1.2 is used
#define CL_HPP_MINIMUM_VERSION 120
#define CL_HPP_MINIMUM_OPENCL_VERSION 120
#if HAVE_OPENCL_HPP
#include <CL/opencl.hpp>
#else
#include <CL/cl2.hpp>                   // supports up to OpenCL 1.2
#endif

#include <string>


namespace Opm
{
namespace Accelerator
{

/// Translate OpenCL error codes to strings
/// Integer - String combinations are defined in CL/cl.h
/// \param[in] error     error code
std::string getErrorString(cl_int error);

} // namespace Accelerator
} // namespace Opm
