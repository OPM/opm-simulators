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

#ifndef OPM_FPGASOLVER_BACKEND_HEADER_INCLUDED
#define OPM_FPGASOLVER_BACKEND_HEADER_INCLUDED

#include <opm/simulators/linalg/bda/BdaSolver.hpp>
#include <opm/simulators/linalg/bda/FPGABILU0.hpp>

#include <linearalgebra/ilu0bicgstab/xilinx/src/sda_app/bicgstab_solver_config.hpp>
#include <linearalgebra/ilu0bicgstab/xilinx/src/sda_app/common/opencl_lib.hpp>
#include <linearalgebra/ilu0bicgstab/xilinx/src/sda_app/common/fpga_functions_bicgstab.hpp>

namespace Opm
{
namespace Accelerator
{

/// This class implements an ilu0-bicgstab solver on FPGA
template <unsigned int block_size>
class FpgaSolverBackend : public BdaSolver<block_size>
{
    typedef BdaSolver<block_size> Base;
    typedef FPGABILU0<block_size> Preconditioner;

    using Base::N;
    using Base::Nb;
    using Base::nnz;
    using Base::nnzb;
    using Base::verbosity;
    using Base::maxit;
    using Base::tolerance;
    using Base::initialized;

private:
    double *rx = nullptr; // reordered x
    double *rb = nullptr; // reordered b
    int *fromOrder = nullptr, *toOrder = nullptr;
    bool analysis_done = false;
    bool level_scheduling = false;

    // LUMat will shallow copy rowPointers and colIndices of mat/rMat
    std::unique_ptr<BlockedMatrix> mat = nullptr;
    BlockedMatrix *rMat = nullptr;
    std::unique_ptr<Preconditioner> prec = nullptr;

    // vectors with data processed by the preconditioner (input to the kernel)
    void **processedPointers = nullptr;
    int *processedSizes = nullptr;

    unsigned int fpga_calls = 0;
    bool perf_call_enabled = true;

    // per call performance metrics
    typedef struct {
      double s_preconditioner_create = 0.0;
      double s_analysis = 0.0;
      double s_reorder = 0.0;
      double s_mem_setup = 0.0;
      double s_mem_h2d = 0.0;
      double s_kernel_exec = 0.0;
      unsigned int n_kernel_exec_cycles = 0;
      float n_kernel_exec_iters = 0.0;
      double s_mem_d2h = 0.0;
      double s_solve = 0.0;
      double s_postprocess = 0.0;
      bool converged = false;
      unsigned int converged_flags = 0;
    } perf_call_metrics_t;
    // cumulative performance metrics
    typedef struct {
      double s_initialization;
      double s_preconditioner_setup;
      double s_preconditioner_create;
      double s_preconditioner_create_min,s_preconditioner_create_max,s_preconditioner_create_avg;
      double s_analysis;
      double s_analysis_min,s_analysis_max,s_analysis_avg;
      double s_reorder;
      double s_reorder_min,s_reorder_max,s_reorder_avg;
      double s_mem_setup;
      double s_mem_setup_min,s_mem_setup_max,s_mem_setup_avg;
      double s_mem_h2d;
      double s_mem_h2d_min,s_mem_h2d_max,s_mem_h2d_avg;
      double s_kernel_exec;
      double s_kernel_exec_min,s_kernel_exec_max,s_kernel_exec_avg;
      unsigned long n_kernel_exec_cycles;
      unsigned long n_kernel_exec_cycles_min,n_kernel_exec_cycles_max,n_kernel_exec_cycles_avg;
      float n_kernel_exec_iters;
      float n_kernel_exec_iters_min,n_kernel_exec_iters_max,n_kernel_exec_iters_avg;
      double s_mem_d2h;
      double s_mem_d2h_min,s_mem_d2h_max,s_mem_d2h_avg;
      double s_solve;
      double s_solve_min,s_solve_max,s_solve_avg;
      double s_postprocess;
      double s_postprocess_min,s_postprocess_max,s_postprocess_avg;
      unsigned int n_converged;
    } perf_total_metrics_t;
    std::vector<perf_call_metrics_t> perf_call;
    perf_total_metrics_t perf_total;
    // fpga_config_bits: bit0=do_reset_debug: if 1, will reset debug flags at each state change, otherwise flags are sticky
    // fpga_config_bits: bit1=absolute_compare: if 1, will compare norm with provided precision value, otherwise it's incremental
    unsigned int fpga_config_bits = 0;
    bool fpga_disabled = false;
    bool platform_awsf1;
    unsigned int debugbufferSize;
    unsigned long int *debugBuffer = nullptr;
    unsigned int *databufferSize = nullptr;
    unsigned char *dataBuffer[RW_BUF] = {nullptr};
    unsigned int debug_outbuf_words;
    int resultsNum;
    int resultsBufferNum;
    unsigned int resultsBufferSize[RES_BUF_MAX];
    unsigned int result_offsets[6];
    unsigned int kernel_cycles, kernel_iter_run;
    double norms[4];
    unsigned char last_norm_idx;
    bool kernel_aborted, kernel_signature, kernel_overflow;
    bool kernel_noresults;
    bool kernel_wrafterend, kernel_dbgfifofull;
    bool use_residuals = false;
    bool use_LU_res = false;
    int sequence = 0;
    // TODO: these values may be sent via command line parameters
    unsigned int abort_cycles = 2000000000; // 2x10^9 @ 300MHz is around 6.6 s
    unsigned int debug_sample_rate = 65535; // max value allowed is 65535, 0 means disabled; reduce to get a finer debug dump
    int nnzValArrays_size = 0;
    int L_nnzValArrays_size = 0;
    int U_nnzValArrays_size = 0;
    // aliases to areas of the host data buffers
    long unsigned int *setupArray = nullptr;
    double **nnzValArrays  = nullptr;
    short unsigned int *columnIndexArray = nullptr;
    unsigned char *newRowOffsetArray = nullptr;
    unsigned int *PIndexArray = nullptr;
    unsigned int *colorSizesArray = nullptr;
    double **L_nnzValArrays = nullptr;
    short unsigned int *L_columnIndexArray = nullptr;
    unsigned char *L_newRowOffsetArray = nullptr;
    unsigned int *L_PIndexArray = nullptr;
    unsigned int *L_colorSizesArray = nullptr;
    double **U_nnzValArrays = nullptr;
    short unsigned int *U_columnIndexArray = nullptr;
    unsigned char *U_newRowOffsetArray = nullptr;
    unsigned int *U_PIndexArray = nullptr;
    unsigned int *U_colorSizesArray = nullptr;
    double *BLKDArray = nullptr;
    double *X1Array = nullptr, *X2Array = nullptr;
    double *R1Array = nullptr, *R2Array = nullptr;
    double *LresArray = nullptr, *UresArray = nullptr;
    double *resultsBuffer[RES_BUF_MAX] = {nullptr}; // alias for data output region
    // OpenCL variables
    cl_device_id device_id;
    cl_context context;
    cl_command_queue commands;
    cl_program program;
    cl_kernel kernel;
    cl_mem cldata[RW_BUF] = {nullptr};
    cl_mem cldebug = nullptr;
    // HW limits/configuration variables
    unsigned int hw_x_vector_elem;
    unsigned int hw_max_row_size;
    unsigned int hw_max_column_size;
    unsigned int hw_max_colors_size;
    unsigned short hw_max_nnzs_per_row;
    unsigned int hw_max_matrix_size;
    bool hw_use_uram;
    bool hw_write_ilu0_results;
    unsigned short hw_dma_data_width;
    unsigned char hw_x_vector_latency;
    unsigned char hw_add_latency;
    unsigned char hw_mult_latency;
    unsigned char hw_mult_num;
    unsigned char hw_num_read_ports;
    unsigned char hw_num_write_ports;
    unsigned short hw_reset_cycles;
    unsigned short hw_reset_settle;
    // debug
    bool reset_data_buffers = false;
    bool fill_results_buffers = false;
    int dump_data_buffers = 0; // 0=disabled, 1=binary format, 2=text format
    bool dump_results = false;
    char *data_dir = nullptr;
    char *basename = nullptr;
    unsigned short rst_assert_cycles = 0;
    unsigned short rst_settle_cycles = 0;

    /// Allocate host memory
    /// \param[in] N              number of nonzeroes, divide by dim*dim to get number of blocks
    /// \param[in] nnz            number of nonzeroes, divide by dim*dim to get number of blocks
    /// \param[in] dim            size of block
    /// \param[in] vals           array of nonzeroes, each block is stored row-wise and contiguous, contains nnz values
    /// \param[in] rows           array of rowPointers, contains N/dim+1 values
    /// \param[in] cols           array of columnIndices, contains nnz values
    void initialize(int N, int nnz, int dim, double *vals, int *rows, int *cols);

    /// Reorder the linear system so it corresponds with the coloring
    /// \param[in] vals           array of nonzeroes, each block is stored row-wise and contiguous, contains nnz values
    /// \param[in] b              input vector
    void update_system(double *vals, double *b);

    /// Analyse sparsity pattern to extract parallelism
    /// \return true iff analysis was successful
    bool analyse_matrix();

    /// Perform ilu0-decomposition
    /// \return true iff decomposition was successful
    bool create_preconditioner();

    /// Solve linear system
    /// \param[inout] res         summary of solver result
    void solve_system(BdaResult &res);

    /// Generate FPGA backend statistics
    void generate_statistics(void);

public:

    /// Construct an fpgaSolver
    /// \param[in] fpga_bitstream             FPGA bitstream file name
    /// \param[in] linear_solver_verbosity    verbosity of fpgaSolver
    /// \param[in] maxit                      maximum number of iterations for fpgaSolver
    /// \param[in] tolerance                  required relative tolerance for fpgaSolver
    /// \param[in] opencl_ilu_reorder         select either level_scheduling or graph_coloring, see ILUReorder.hpp for explanation
    FpgaSolverBackend(std::string fpga_bitstream, int linear_solver_verbosity, int maxit, double tolerance, ILUReorder opencl_ilu_reorder);

    /// Destroy an fpgaSolver, and free memory
    ~FpgaSolverBackend();

    /// Solve linear system, A*x = b, matrix A must be in blocked-CSR format
    /// \param[in] N              number of rows, divide by dim to get number of blockrows
    /// \param[in] nnz            number of nonzeroes, divide by dim*dim to get number of blocks
    /// \param[in] dim            size of block
    /// \param[in] vals           array of nonzeroes, each block is stored row-wise and contiguous, contains nnz values
    /// \param[in] rows           array of rowPointers, contains N/dim+1 values
    /// \param[in] cols           array of columnIndices, contains nnz values
    /// \param[in] b              input vector, contains N values
    /// \param[in] wellContribs   WellContributions, not used in FPGA solver because it requires them already added to matrix A
    /// \param[inout] res         summary of solver result
    /// \return                   status code
    SolverStatus solve_system(int N, int nnz, int dim, double *vals, int *rows, int *cols, double *b, WellContributions& wellContribs, BdaResult &res) override;

    /// Get result after linear solve, and peform postprocessing if necessary
    /// \param[inout] x           resulting x vector, caller must guarantee that x points to a valid array
    void get_result(double *x) override;

}; // end class fpgaSolverBackend

} // namespace Accelerator
} // namespace Opm

#endif

