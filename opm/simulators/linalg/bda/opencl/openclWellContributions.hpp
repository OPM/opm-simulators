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

#ifndef WELLCONTRIBUTIONS_OPENCL_HEADER_INCLUDED
#define WELLCONTRIBUTIONS_OPENCL_HEADER_INCLUDED

#include <opm/simulators/linalg/bda/WellContributions.hpp>

#include <opm/simulators/linalg/bda/opencl/opencl.hpp>
#include <opm/simulators/linalg/bda/opencl/openclKernels.hpp>

#include <memory>
#include <vector>


namespace Opm
{

class WellContributionsOCL : public WellContributions
{
public:
    void setKernel(Opm::Accelerator::stdwell_apply_kernel_type *kernel_,
                   Opm::Accelerator::stdwell_apply_no_reorder_kernel_type *kernel_no_reorder_);
    void setOpenCLEnv(cl::Context *context_, cl::CommandQueue *queue_);

    /// Since the rows of the matrix are reordered, the columnindices of the matrixdata is incorrect
    /// Those indices need to be mapped via toOrder
    /// \param[in] toOrder    array with mappings
    /// \param[in] reorder    whether reordering is actually used or not
    void setReordering(int* toOrder, bool reorder);
    void apply_stdwells(cl::Buffer d_x, cl::Buffer d_y, cl::Buffer d_toOrder);
    void apply_mswells(cl::Buffer d_x, cl::Buffer d_y);
    void apply(cl::Buffer d_x, cl::Buffer d_y, cl::Buffer d_toOrder);

protected:
    /// Allocate memory for the StandardWells
    void APIalloc() override;

    void APIaddMatrix(MatrixType type, int *colIndices, double *values, unsigned int val_size) override;

    cl::Context* context;
    cl::CommandQueue* queue;
    Opm::Accelerator::stdwell_apply_kernel_type* kernel;
    Opm::Accelerator::stdwell_apply_no_reorder_kernel_type* kernel_no_reorder;
    std::vector<cl::Event> events;

    std::unique_ptr<cl::Buffer> d_Cnnzs_ocl, d_Dnnzs_ocl, d_Bnnzs_ocl;
    std::unique_ptr<cl::Buffer> d_Ccols_ocl, d_Bcols_ocl;
    std::unique_ptr<cl::Buffer> d_val_pointers_ocl;

    bool reorder = false;
    int *h_toOrder = nullptr;
    std::vector<double> h_x;
    std::vector<double> h_y;
};

} //namespace Opm

#endif
