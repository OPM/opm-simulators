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

#include <config.h>

#include <cmath>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/material/common/Unused.hpp>
#include <opm/common/ErrorMacros.hpp>

#include <opm/simulators/linalg/bda/FPGASolverBackend.hpp>
#include <opm/simulators/linalg/bda/FPGAUtils.hpp>
#include <opm/simulators/linalg/bda/Reorder.hpp>

// if defined, any FPGA kernel failure will terminate flow; otherwise, the FPGA
// kernel will be disabled and execution will continue using DUNE
#define FPGA_EXIT_WITH_HW_FAILURE
//#undef FPGA_EXIT_WITH_HW_FAILURE

// if defined, the function generate_statistics will create a CSV-formatted file
// with detailed statistics about the FPGA backend performance
//#define FPGA_STATISTICS_FILE_ENABLED
#undef FPGA_STATISTICS_FILE_ENABLED

namespace bda
{

using Opm::OpmLog;

template <unsigned int block_size>
FpgaSolverBackend<block_size>::FpgaSolverBackend(std::string fpga_bitstream, int verbosity_, int maxit_, double tolerance_, ILUReorder opencl_ilu_reorder) : BdaSolver<block_size>(fpga_bitstream, verbosity_, maxit_, tolerance_)
{
    int err;
    std::ostringstream oss;
    double start = second();

    // currently, only block size == 3 is supported by the FPGA backend
    assert(block_size == 3);

    if (verbosity < 1) {
        perf_call_enabled = false;
    }
    // setup bitstream name and other parameters
    if (fpga_bitstream.compare("") == 0) {
        OPM_THROW(std::logic_error, "fpgaSolver called but bitstream file has not been specified");
    }
    if (!fileExists(fpga_bitstream.c_str())) {
        OPM_THROW(std::logic_error, "fpgaSolver called but bitstream file specified does not exists or is not readable");
    }
    // -----------------------------
    // FPGA: setup the OpenCL platform
    // -----------------------------
    std::string main_kernel_name(KERNEL_NAME); // macro defined in bicgstab_solver_config.hpp
    // auto-select the proper FPGA device and create context and other CL objects
    err = setup_opencl(nullptr, &device_id, &context, &commands, &program, &kernel, main_kernel_name.c_str(), fpga_bitstream.c_str(), &platform_awsf1);
    if (err != 0) {
        oss << "Failed to setup the OpenCL device (" << err << ")";
        OPM_THROW(std::logic_error, oss.str());
    }
    oss << "Detected FPGA platform type is ";
    if (platform_awsf1) { oss << "AWS-F1."; } else { oss << "Xilinx Alveo."; }
    OpmLog::info(oss.str());
    oss.str("");
    oss.clear();
    // -----------------------------
    // FPGA: setup the debug buffer
    // -----------------------------
    // set kernel debug lines depending on an environment variable
    const char *xem = getenv("XCL_EMULATION_MODE");
    if ((xem != nullptr) && (strcmp(xem, "sw_emu") == 0 || strcmp(xem, "hw_emu") == 0)) {
        debug_outbuf_words = DEBUG_OUTBUF_WORDS_MAX_EMU;
        oss << "Detected co-simulation mode, debug_outbuf_words set to " << debug_outbuf_words << ".\n";
        OpmLog::info(oss.str());
        oss.str("");
        oss.clear();
    } else {
        // set to 2 to reduce overhead in reading back and interpreting the debug lines;
        // increase to get more debug info from the kernel
        // range is 2..DEBUG_OUTBUF_WORDS_MAX-1
        debug_outbuf_words = 2;
    }

    // host debug buffer setup
    err = fpga_setup_host_debugbuf(debug_outbuf_words, &debugBuffer, &debugbufferSize);
    if (err != 0) {
        oss << "Failed to call fpga_setup_host_debug_buffer (" << err << ")";
        OPM_THROW(std::logic_error, oss.str());
    }
    // device debug buffer setup
    err = fpga_setup_device_debugbuf(context, debugBuffer, &cldebug, debugbufferSize);
    if (err != 0) {
        oss << "Failed to call fpga_setup_device_debug_buffer (" << err << ").\n";
        OPM_THROW(std::logic_error, oss.str());
    }
    // copy debug buffer to device
    err = fpga_copy_to_device_debugbuf(commands, cldebug, debugBuffer, debugbufferSize, debug_outbuf_words);
    if (err != 0) {
        oss << "Failed to call fpga_copy_to_device_debugbuf (" << err << ").\n";
        OPM_THROW(std::logic_error, oss.str());
    }
    // ------------------------------------------------
    // FPGA: query the kernel for limits/configuration
    // ------------------------------------------------
    err = fpga_kernel_query(context, commands, kernel, cldebug,
                            debugBuffer, debug_outbuf_words,
                            rst_assert_cycles, rst_settle_cycles,
                            &hw_x_vector_elem, &hw_max_row_size,
                            &hw_max_column_size, &hw_max_colors_size,
                            &hw_max_nnzs_per_row, &hw_max_matrix_size,
                            &hw_use_uram, &hw_write_ilu0_results,
                            &hw_dma_data_width, &hw_mult_num,
                            &hw_x_vector_latency, &hw_add_latency, &hw_mult_latency,
                            &hw_num_read_ports, &hw_num_write_ports,
                            &hw_reset_cycles, &hw_reset_settle);
    if (err != 0) {
        oss << "Failed to call fpga_kernel_query (" << err << ")";
        OPM_THROW(std::logic_error, oss.str());
    }

    if (verbosity >= 1) {
        oss << "FPGA kernel limits/configuration:\n";
        oss << "  x_vector_elem=" << hw_max_colors_size << ", max_row_size=" << hw_max_nnzs_per_row << ", max_column_size=" << hw_max_matrix_size << "\n";
        oss << "  max_colors_size=" << hw_x_vector_elem << ", max_nnzs_per_row=" << hw_max_row_size << ", max_matrix_size=" << hw_max_column_size << "\n";
        oss << "  use_uram=" << hw_use_uram << ", write_ilu0_results=" << hw_write_ilu0_results << "\n";
        oss << "  dma_data_width=" << hw_dma_data_width << ", mult_num=" << (unsigned int)hw_mult_num << "\n";
        oss << "  x_vector_latency=" << (unsigned int)hw_x_vector_latency << "\n";
        oss << "  add_latency=" << (unsigned int)hw_add_latency << ", mult_latency=" << (unsigned int)hw_mult_latency << "\n";
        oss << "  num_read_ports=" << (unsigned int)hw_num_read_ports << ", num_write_ports=" << (unsigned int)hw_num_write_ports << "\n";
        oss << "  reset_cycles=" << hw_reset_cycles << ", reset_settle=" << hw_reset_settle;
        OpmLog::info(oss.str());
        oss.str("");
        oss.clear();
    }

    // check that LU results are generated by the kernel
    if (use_LU_res && !hw_write_ilu0_results) {
        OpmLog::warning("Kernel reports that LU results are not written to memory, but use_LU_res is set; disabling LU results usage");
        oss.str("");
        oss.clear();
        use_LU_res = false;
    }

    // setup preconditioner
    double start_prec = second();
    prec = std::make_unique<Preconditioner>(opencl_ilu_reorder, verbosity_, hw_max_row_size, hw_max_column_size, hw_max_nnzs_per_row, hw_max_colors_size);
    perf_total.s_preconditioner_setup = second() - start_prec;

    if (opencl_ilu_reorder == ILUReorder::LEVEL_SCHEDULING) {
        level_scheduling = true;
    }

    perf_total.s_initialization = second() - start;
} // end fpgaSolverBackend


template <unsigned int block_size>
FpgaSolverBackend<block_size>::~FpgaSolverBackend()
{
    if (verbosity >= 1) {
        generate_statistics();
    }
    delete[] rx;
    delete[] rb;
    if (nnzValArrays != nullptr) { free(nnzValArrays); }
    if (L_nnzValArrays != nullptr) { free(L_nnzValArrays); }
    if (U_nnzValArrays != nullptr) { free(U_nnzValArrays); }
    // FPGA: buffers
    free(debugBuffer);
    for (int b = 0; b < RW_BUF; b++) {
        free(dataBuffer[b]);
    }
    free(databufferSize);
    // FPGA: OpenCL objects
    if (cldebug != nullptr) { clReleaseMemObject(cldebug); }
    for (int b = 0; b < RW_BUF; b++) {
        if (cldata[b] != nullptr) {
            clReleaseMemObject(cldata[b]);
        }
    }
    clReleaseCommandQueue(commands);
    clReleaseContext(context);
    clReleaseKernel(kernel);
    clReleaseProgram(program);
    clReleaseDevice(device_id);
} // end ~fpgaSolverBackend()


// copy result to host memory
// caller must be sure that x is a valid array
template <unsigned int block_size>
void FpgaSolverBackend<block_size>::get_result(double *x_)
{
    double start = 0;

    if (perf_call_enabled) {
        start = second();
    }
    // apply to results the reordering (stored in toOrder)
    reorderBlockedVectorByPattern<block_size>(mat->Nb, rx, toOrder, x_);
    // TODO: check if it is more efficient to avoid copying resultsBuffer[0] to rx in solve_system (private)
    if (perf_call_enabled) {
        perf_call.back().s_postprocess = second() - start;
    }
} // end get_result()


template <unsigned int block_size>
SolverStatus FpgaSolverBackend<block_size>::solve_system(int N_, int nnz_, int dim, double *vals, int *rows, int *cols, double *b, WellContributions& wellContribs OPM_UNUSED, BdaResult &res)
{
    if (initialized == false) {
        initialize(N_, nnz_,  dim, vals, rows, cols);
        if (!analyse_matrix()) {
            return SolverStatus::BDA_SOLVER_ANALYSIS_FAILED;
        }
    }
    perf_call.emplace_back();
    update_system(vals, b);
    if (!create_preconditioner()) {
        return SolverStatus::BDA_SOLVER_CREATE_PRECONDITIONER_FAILED;
    }
    solve_system(res);

    if (verbosity >= 1) {
        std::ostringstream oss;
        oss << "fpgaSolverBackend::" << __func__ << " - converged: " << res.converged << \
            ", iterations: " << res.iterations << ", reduction: " << res.reduction << \
            ", conv_rate: " << res.conv_rate << ", elapsed: " << res.elapsed;
        OpmLog::info(oss.str());
    }
    return SolverStatus::BDA_SOLVER_SUCCESS;
}


template <unsigned int block_size>
void FpgaSolverBackend<block_size>::initialize(int N_, int nnz_, int dim, double *vals, int *rows, int *cols)
{
    double start = second();
    this->N = N_;
    this->nnz = nnz_;
    this->nnzb = nnz_ / block_size / block_size;
    Nb = (N + dim - 1) / dim;

    // allocate host memory for matrices and vectors
    // actual data for mat points to std::vector.data() in ISTLSolverEbos, so no alloc/free here
    mat.reset(new BlockedMatrix<block_size>(N_ / block_size, nnz_ / block_size / block_size, vals, cols, rows));

    std::ostringstream oss;
    oss << "Initializing FPGA data, matrix size: " << this->N << " blocks, nnz: " << this->nnzb << " blocks, " << \
        "block size: " << dim << ", total nnz: " << this->nnz << "\n";
    oss << "Maxit: " << maxit << std::scientific << ", tolerance: " << tolerance;
    OpmLog::info(oss.str());

    rx = new double[roundUpTo(N_, CACHELINE_BYTES / sizeof(double))];
    rb = new double[roundUpTo(N_, CACHELINE_BYTES / sizeof(double))];

    perf_total.s_initialization += second() - start;
    initialized = true;
} // end initialize()


template <unsigned int block_size>
bool FpgaSolverBackend<block_size>::analyse_matrix()
{
    std::ostringstream oss;
    int err;

    double start = second();
    bool success = prec->init(mat.get());

    if (!success) {
        OpmLog::warning("Preconditioner for FPGA solver failed to initialize");
        return success;
    }

    toOrder = prec->getToOrder();
    fromOrder = prec->getFromOrder();
    rMat = prec->getRMat();
    processedPointers = prec->getResultPointers();
    processedSizes = prec->getResultSizes();
    processedPointers[19] = rb;
    processedPointers[20] = rx;
    nnzValArrays_size = static_cast<int*>(processedPointers[5])[0];
    L_nnzValArrays_size = static_cast<int*>(processedPointers[11])[0];
    U_nnzValArrays_size = static_cast<int*>(processedPointers[17])[0];
    // -------------------------------------
    // FPGA: setup host/device data buffers
    // -------------------------------------
    // allocate memory and setup data layout
    err = fpga_setup_host_datamem(level_scheduling, fpga_config_bits,
                                  processedSizes,
                                  &setupArray,
                                  &nnzValArrays,   &nnzValArrays_size,   &columnIndexArray,   &newRowOffsetArray,   &PIndexArray,   &colorSizesArray,
                                  &L_nnzValArrays, &L_nnzValArrays_size, &L_columnIndexArray, &L_newRowOffsetArray, &L_PIndexArray, &L_colorSizesArray,
                                  &U_nnzValArrays, &U_nnzValArrays_size, &U_columnIndexArray, &U_newRowOffsetArray, &U_PIndexArray, &U_colorSizesArray,
                                  &BLKDArray, &X1Array, &R1Array,
                                  &X2Array, &R2Array,
                                  &LresArray, &UresArray,
                                  &databufferSize, dataBuffer,
                                  result_offsets, 1 /*num_nnz_arrays*/,
                                  true /*reset_data_buffers*/,  /* WARNING: leave reset_data_buffers always ENABLED to avoid data corruption! */
                                  debugbufferSize);
    if (err) {
        oss << "Failed to call fpga_setup_host_datamem (" << err << ")";
        OPM_THROW(std::logic_error, oss.str());
    }

    // results buffers setup
    if (use_LU_res) {
        resultsBufferNum = 4;
    } else {
        resultsBufferNum = 2;
    }
    if (resultsBufferNum > RES_BUF_MAX) {
        oss << "Number of results buffer (" << resultsBufferNum << ") is out of range (max " << RES_BUF_MAX << ")";
        OPM_THROW(std::logic_error, oss.str());
    }

    resultsNum = processedSizes[0]; // rowSize, invariant between system solves
    for (int i = 0; i < resultsBufferNum; i++) {
        resultsBufferSize[i] = roundUpTo(resultsNum, CACHELINE_BYTES / sizeof(double)) * sizeof(double);
    }

    // device data memory setup
    err = fpga_setup_device_datamem(context, databufferSize, dataBuffer, cldata);
    if (err != 0) {
        oss << "Failed to call fpga_setup_device_datamem (" << err << ")";
        OPM_THROW(std::logic_error, oss.str());
    }

    // ------------------------------------
    // FPGA: setup the kernel's parameters
    // ------------------------------------
    err = fpga_set_kernel_parameters(kernel, abort_cycles, debug_outbuf_words - 1, maxit,
                                     debug_sample_rate, tolerance, cldata, cldebug);
    if (err != 0) {
        oss << "Failed to call fpga_set_kernel_parameters (" << err << ")";
        OPM_THROW(std::logic_error, oss.str());
    }

    perf_total.s_analysis = second() - start;
    analysis_done = true;

    return success;
} // end analyse_matrix()


template <unsigned int block_size>
bool FpgaSolverBackend<block_size>::create_preconditioner()
{
    double start = 0;

    if (perf_call_enabled) {
        start = second();
    }
    memset(rx, 0, sizeof(double) * N);
    bool result = prec->create_preconditioner(mat.get());
    if (!result) {
        OpmLog::warning("fpgaSolverBackend: create_preconditioner failed");
    }

    if (perf_call_enabled) {
        perf_call.back().s_preconditioner_create = second() - start;
    }
    return result;
} // end create_preconditioner()


template <unsigned int block_size>
void FpgaSolverBackend<block_size>::solve_system(BdaResult &res)
{
    std::ostringstream oss;
    int err;
    double start = 0, start_total = 0;

    // ------------------------------------
    // FPGA: return immediately if FPGA is disabled
    // ------------------------------------
    if (fpga_disabled) {
        res.converged = false;
        OpmLog::warning("FPGA is disabled, fallback to SW execution");
        return;
    }

    fpga_calls++;

    if (perf_call_enabled) {
        start = second();
        start_total = start;
    }

    // check if any buffer is larger than the size set in preconditioner->init
    // TODO: add check for all other buffer sizes that may overflow?
    err = 0;
    if ( ((int *)processedPointers[5])[0]  > nnzValArrays_size ||
            ((int *)processedPointers[11])[0] > L_nnzValArrays_size ||
            ((int *)processedPointers[17])[0] > U_nnzValArrays_size ) {
        err = 1;
    }
    if (err != 0) {
        OPM_THROW(std::logic_error, "A buffer size is larger than the initial allocation in solve_system (check preconditioner init)");
    }

    // ------------------------------------
    // FPGA: copy input data to host data buffers
    // ------------------------------------
    if (perf_call_enabled) {
        start = second();
    }
    err = fpga_copy_host_datamem(
              processedPointers, processedSizes, setupArray,
              nnzValArrays,   &nnzValArrays_size,   columnIndexArray,   newRowOffsetArray,   PIndexArray,   colorSizesArray,
              L_nnzValArrays, &L_nnzValArrays_size, L_columnIndexArray, L_newRowOffsetArray, L_PIndexArray, L_colorSizesArray,
              U_nnzValArrays, &U_nnzValArrays_size, U_columnIndexArray, U_newRowOffsetArray, U_PIndexArray, U_colorSizesArray,
              BLKDArray, X1Array, R1Array, X2Array, R2Array,
              use_LU_res, LresArray, UresArray,
              databufferSize, dataBuffer,
              1 /* nnzValArrays_num */,
              reset_data_buffers, fill_results_buffers,
              dump_data_buffers, fpga_calls);
    if (perf_call_enabled) {
        perf_call.back().s_mem_setup = second() - start;
    }
    if (err != 0) {
        oss << "Failed to call fpga_copy_to_device_debugbuf (" << err << ")";
        OPM_THROW(std::logic_error, oss.str());
    }
    // ------------------------------------
    // FPGA: copy buffers to device
    // ------------------------------------
    // copy debug buffer to device
    if (perf_call_enabled) {
        start = second();
    }
    err = fpga_copy_to_device_debugbuf(commands,
                                       cldebug, debugBuffer, debugbufferSize,
                                       debug_outbuf_words);
    if (err != 0) {
        oss << "Failed to call fpga_copy_to_device_debugbuf (" << err << ")";
        OPM_THROW(std::logic_error, oss.str());
    }
    // copy data buffers to device
    err = fpga_copy_to_device_datamem(commands, RW_BUF, cldata);
    if (err != 0) {
        oss << "Failed to call fpga_copy_to_device_datamem (" << err << ")";
        OPM_THROW(std::logic_error, oss.str());
    }
    if (perf_call_enabled) {
        perf_call.back().s_mem_h2d = second() - start;
    }
    // ------------------------------------
    // FPGA: execute the kernel
    // ------------------------------------
    double time_elapsed_ms;
    if (perf_call_enabled) {
        start = second();
    }
    err = fpga_kernel_run(commands, kernel, &time_elapsed_ms);
    if (perf_call_enabled) {
        perf_call.back().s_kernel_exec = second() - start;
    }
    if (err != 0) {
        oss << "Failed to call fpga_kernel_run (" << err << ")";
        OPM_THROW(std::logic_error, oss.str());
    }
    // ----------------------------------------
    // FPGA: read back debug buffer from device
    // ----------------------------------------
    if (perf_call_enabled) {
        start = second();
    }
    err = fpga_copy_from_device_debugbuf((bool)(verbosity < 10),
                                         commands,
                                         debug_outbuf_words, debugbufferSize,
                                         cldebug, debugBuffer,
                                         abort_cycles,
                                         &kernel_cycles, &kernel_iter_run,
                                         norms, &last_norm_idx,
                                         &kernel_aborted, &kernel_signature, &kernel_overflow, &kernel_noresults,
                                         &kernel_wrafterend, &kernel_dbgfifofull);
    if (err != 0) {
        oss << "Failed to call fpga_copy_from_device_debugbuf (" << err << ")";
        OPM_THROW(std::logic_error, oss.str());
    }
    if (kernel_wrafterend) {
        OpmLog::warning("Detected recoverable FPGA error: kernel write after end");
    }
    if (kernel_dbgfifofull) {
        OpmLog::warning("Detected recoverable FPGA error: debug FIFO full");
    }
    if (kernel_aborted || kernel_signature || kernel_overflow) {
#if defined(FPGA_EXIT_WITH_HW_FAILURE)
        oss << "Detected unrecoverable FPGA error (ABRT=" << kernel_aborted << \
            ",SIG=" << kernel_signature << ",OVF=" << kernel_overflow << ")";
        OPM_THROW(std::logic_error, oss.str());
#else
        oss << "Detected unrecoverable FPGA error (ABRT=" << kernel_aborted << \
            ",SIG=" << kernel_signature << ",OVF=" << kernel_overflow << ")\n";
        oss << "Disabling FPGA kernel: execution will continue with SW kernel";
        OpmLog::warning(oss.str());
        oss.str("");
        oss.clear();
        fpga_disabled = true;
#endif
    }
    if (perf_call_enabled) {
        perf_call.back().n_kernel_exec_cycles = kernel_cycles;
    }
    // copy (back) results only if FPGA is not disabled
    if (!fpga_disabled) {
        if (kernel_noresults) {
            OpmLog::warning("FPGA kernel did not return results because the required precision is already reached");
            // rx still contains zeros from initial guess
        } else {
            // ------------------------------------
            // FPGA: read back results from device
            // ------------------------------------
            err = fpga_map_results(even(kernel_iter_run),
                                   use_residuals, use_LU_res, commands,
                                   resultsNum, resultsBufferNum, resultsBufferSize,
                                   debugbufferSize,
                                   cldata, resultsBuffer,
                                   result_offsets,
                                   dump_results, data_dir, basename, sequence);
            if (err != 0) {
                oss << "Failed to call fpga_map_results (" << err << ")";
                OPM_THROW(std::logic_error, oss.str());
            }
            // TODO: copy results buffers to reordering output buffers
            memcpy(rx, resultsBuffer[0], resultsNum * sizeof(double));
            err = fpga_unmap_results(even(kernel_iter_run),
                                     use_residuals, use_LU_res,
                                     commands, cldata, resultsBuffer);
            if (err != 0) {
                oss << "Failed to call fpga_unmap_results (" << err << ")";
                OPM_THROW(std::logic_error, oss.str());
            }
        }
    }
    // set results and update statistics (if enabled)
    if (perf_call_enabled) {
        perf_call.back().s_mem_d2h = second() - start;
    }
    float iter = ((float)kernel_iter_run / 2.0) + 0.5; // convert from half iteration int to actual iterationns
    res.iterations = (int)iter;
    res.reduction = norms[0] / norms[last_norm_idx]; // norms[0] is the initial norm
    res.conv_rate = pow(res.reduction, 1.0 / iter);
    res.elapsed = second() - start_total;
    if (perf_call_enabled) {
        perf_call.back().s_solve = res.elapsed;
        perf_call.back().n_kernel_exec_iters = iter;
    }
    // convergence depends on number of iterations reached and hw execution errors
    res.converged = true;
    if (fpga_disabled || kernel_aborted || kernel_signature || kernel_overflow || iter >= (float)maxit) {
        res.converged = false;
        if (verbosity >= 1) {
            oss << "FPGA kernel did not converge, reason: fpga_disabled=" << fpga_disabled << \
                ", kernel_aborted=" << kernel_aborted << ", kernel_signature=" << kernel_signature << \
                ", kernel_overflow=" << kernel_overflow << ", (iter>=" << maxit << ")=" << (iter >= (float)maxit);
            OpmLog::warning(oss.str());
            oss.str("");
            oss.clear();
        }
    }
    if (perf_call_enabled) {
        perf_call.back().converged = res.converged;
        perf_call.back().converged_flags = ((unsigned int)fpga_disabled) +
                                           ((unsigned int)kernel_aborted << 1) + ((unsigned int)kernel_signature << 2) +
                                           ((unsigned int)kernel_overflow << 3) + ((unsigned int)(iter >= (float)maxit) << 4);
    }
} // end solve_system()


template <unsigned int block_size>
void FpgaSolverBackend<block_size>::update_system(double *vals, double *b)
{
    double start = 0;

    mat->nnzValues = vals;
    // reorder inputs using previously found ordering (stored in fromOrder)
    if (perf_call_enabled) {
        start = second();
    }
    reorderBlockedVectorByPattern<block_size>(mat->Nb, b, fromOrder, rb);
    if (perf_call_enabled) {
        perf_call.back().s_reorder = second() - start;
    }
} // end update_system()


template <unsigned int block_size>
void FpgaSolverBackend<block_size>::generate_statistics()
{
    std::ostringstream oss;
    unsigned int conv_iter = 0, conv_ovf = 0;

    if (!perf_call_enabled || fpga_calls == 0) {
        OpmLog::warning("FPGA statistics were not collected");
        return;
    }
    std::printf("--- FPGA statistics ---\n");
    std::printf("total solver calls..........: %u\n", fpga_calls);
    std::printf("time initialization.........: %8.6f s\n", perf_total.s_initialization);
    std::printf("time preconditioner setup...: %8.6f s\n", perf_total.s_preconditioner_setup);

#if defined(FPGA_STATISTICS_FILE_ENABLED)
    // DEBUG: this can be enabled to gather all the statistics in a CSV-formatted file
    FILE *fout = fopen("fpga_statistics_details.csv", "w");
    if (fout != nullptr) {
        std::fprintf(fout, "call,preconditioner_create,analysis,reorder,mem_setup,mem_h2d,kernel_exec,kernel_cycles,kernel_iters,mem_d2h,solve,postprocess,converged\n");
    }
#endif
    unsigned int num_data_points = perf_call.size();
    for (unsigned int i = 0; i < num_data_points; i++) {
        perf_total.s_preconditioner_create += perf_call[i].s_preconditioner_create;
        if (perf_call[i].s_preconditioner_create > perf_total.s_preconditioner_create_max) { perf_total.s_preconditioner_create_max = perf_call[i].s_preconditioner_create; }
        if (perf_call[i].s_preconditioner_create < perf_total.s_preconditioner_create_min) { perf_total.s_preconditioner_create_min = perf_call[i].s_preconditioner_create; }
        perf_total.s_analysis += perf_call[i].s_analysis;
        if (perf_call[i].s_analysis > perf_total.s_analysis_max) { perf_total.s_analysis_max = perf_call[i].s_analysis; }
        if (perf_call[i].s_analysis < perf_total.s_analysis_min) { perf_total.s_analysis_min = perf_call[i].s_analysis; }
        perf_total.s_reorder += perf_call[i].s_reorder;
        if (perf_call[i].s_reorder > perf_total.s_reorder_max) { perf_total.s_reorder_max = perf_call[i].s_reorder; }
        if (perf_call[i].s_reorder < perf_total.s_reorder_min) { perf_total.s_reorder_min = perf_call[i].s_reorder; }
        perf_total.s_mem_setup += perf_call[i].s_mem_setup;
        if (perf_call[i].s_mem_setup > perf_total.s_mem_setup_max) { perf_total.s_mem_setup_max = perf_call[i].s_mem_setup; }
        if (perf_call[i].s_mem_setup < perf_total.s_mem_setup_min) { perf_total.s_mem_setup_min = perf_call[i].s_mem_setup; }
        perf_total.s_mem_h2d += perf_call[i].s_mem_h2d;
        if (perf_call[i].s_mem_h2d > perf_total.s_mem_h2d_max) { perf_total.s_mem_h2d_max = perf_call[i].s_mem_h2d; }
        if (perf_call[i].s_mem_h2d < perf_total.s_mem_h2d_min) { perf_total.s_mem_h2d_min = perf_call[i].s_mem_h2d; }
        perf_total.s_kernel_exec += perf_call[i].s_kernel_exec;
        if (perf_call[i].s_kernel_exec > perf_total.s_kernel_exec_max) { perf_total.s_kernel_exec_max = perf_call[i].s_kernel_exec; }
        if (perf_call[i].s_kernel_exec < perf_total.s_kernel_exec_min) { perf_total.s_kernel_exec_min = perf_call[i].s_kernel_exec; }
        perf_total.n_kernel_exec_cycles += (unsigned long)perf_call[i].n_kernel_exec_cycles;
        if (perf_call[i].n_kernel_exec_cycles > perf_total.n_kernel_exec_cycles_max) { perf_total.n_kernel_exec_cycles_max = perf_call[i].n_kernel_exec_cycles; }
        if (perf_call[i].n_kernel_exec_cycles < perf_total.n_kernel_exec_cycles_min) { perf_total.n_kernel_exec_cycles_min = perf_call[i].n_kernel_exec_cycles; }
        perf_total.n_kernel_exec_iters += perf_call[i].n_kernel_exec_iters;
        if (perf_call[i].n_kernel_exec_iters > perf_total.n_kernel_exec_iters_max) { perf_total.n_kernel_exec_iters_max = perf_call[i].n_kernel_exec_iters; }
        if (perf_call[i].n_kernel_exec_iters < perf_total.n_kernel_exec_iters_min) { perf_total.n_kernel_exec_iters_min = perf_call[i].n_kernel_exec_iters; }
        perf_total.s_mem_d2h += perf_call[i].s_mem_d2h;
        if (perf_call[i].s_mem_d2h > perf_total.s_mem_d2h_max) { perf_total.s_mem_d2h_max = perf_call[i].s_mem_d2h; }
        if (perf_call[i].s_mem_d2h < perf_total.s_mem_d2h_min) { perf_total.s_mem_d2h_min = perf_call[i].s_mem_d2h; }
        perf_total.s_solve += perf_call[i].s_solve;
        if (perf_call[i].s_solve > perf_total.s_solve_max) { perf_total.s_solve_max = perf_call[i].s_solve; }
        if (perf_call[i].s_solve < perf_total.s_solve_min) { perf_total.s_solve_min = perf_call[i].s_solve; }
        perf_total.s_postprocess += perf_call[i].s_postprocess;
        if (perf_call[i].s_postprocess > perf_total.s_postprocess_max) { perf_total.s_postprocess_max = perf_call[i].s_postprocess; }
        if (perf_call[i].s_postprocess < perf_total.s_postprocess_min) { perf_total.s_postprocess_min = perf_call[i].s_postprocess; }
        perf_total.n_converged += (unsigned int)perf_call[i].converged;
        if (perf_call[i].converged_flags & 1 << 4) { conv_iter += 1; }
        if (perf_call[i].converged_flags & 1 << 3) { conv_ovf += 1; }
#if defined(FPGA_STATISTICS_FILE_ENABLED)
        if (fout != nullptr) {
            std::fprintf(fout, "%d,%8.6f,%8.6f,%8.6f,%8.6f,%8.6f,%8.6f,%u,%.1f,%8.6f,%8.6f,%8.6f,%u\n",
                         i, perf_call[i].s_preconditioner_create, perf_call[i].s_analysis, perf_call[i].s_reorder,
                         perf_call[i].s_mem_setup, perf_call[i].s_mem_h2d, perf_call[i].s_kernel_exec, perf_call[i].n_kernel_exec_cycles,
                         perf_call[i].n_kernel_exec_iters, perf_call[i].s_mem_d2h, perf_call[i].s_solve, perf_call[i].s_postprocess,
                         (unsigned int)perf_call[i].converged);
        }
#endif
    }
#if defined(FPGA_STATISTICS_FILE_ENABLED)
    if (fout != nullptr) {
        fclose(fout);
    }
#endif
    perf_total.s_preconditioner_create_avg = perf_total.s_preconditioner_create / num_data_points;
    perf_total.s_analysis_avg = perf_total.s_analysis / num_data_points;
    perf_total.s_reorder_avg = perf_total.s_reorder / num_data_points;
    perf_total.s_mem_setup_avg = perf_total.s_mem_setup / num_data_points;
    perf_total.s_mem_h2d_avg = perf_total.s_mem_h2d / num_data_points;
    perf_total.s_kernel_exec_avg = perf_total.s_kernel_exec / num_data_points;
    perf_total.n_kernel_exec_cycles_avg = perf_total.n_kernel_exec_cycles / num_data_points;
    perf_total.n_kernel_exec_iters_avg = perf_total.n_kernel_exec_iters / num_data_points;
    perf_total.s_mem_d2h_avg = perf_total.s_mem_d2h / num_data_points;
    perf_total.s_solve_avg = perf_total.s_solve / num_data_points;
    perf_total.s_postprocess_avg = perf_total.s_postprocess / num_data_points;
    std::printf("time preconditioner creation: total %8.6f s, avg %8.6f s, min %8.6f s, max %8.6f s\n",
                perf_total.s_preconditioner_create, perf_total.s_preconditioner_create_avg, perf_total.s_preconditioner_create_min, perf_total.s_preconditioner_create_max);
    std::printf("time analysis...............: total %8.6f s, avg %8.6f s, min %8.6f s, max %8.6f s\n",
                perf_total.s_analysis, perf_total.s_analysis_avg, perf_total.s_analysis_min, perf_total.s_analysis_max);
    std::printf("time reorder................: total %8.6f s, avg %8.6f s, min %8.6f s, max %8.6f s\n",
                perf_total.s_reorder, perf_total.s_reorder_avg, perf_total.s_reorder_min, perf_total.s_reorder_max);
    std::printf("time memory setup...........: total %8.6f s, avg %8.6f s, min %8.6f s, max %8.6f s\n",
                perf_total.s_mem_setup, perf_total.s_mem_setup_avg, perf_total.s_mem_setup_min, perf_total.s_mem_setup_max);
    std::printf("time memory host2dev........: total %8.6f s, avg %8.6f s, min %8.6f s, max %8.6f s\n",
                perf_total.s_mem_h2d, perf_total.s_mem_h2d_avg, perf_total.s_mem_h2d_min, perf_total.s_mem_h2d_max);
    std::printf("time kernel execution.......: total %8.6f s, avg %8.6f s, min %8.6f s, max %8.6f s\n",
                perf_total.s_kernel_exec, perf_total.s_kernel_exec_avg, perf_total.s_kernel_exec_min, perf_total.s_kernel_exec_max);
    std::printf("cycles kernel execution.....: total %lu, avg %lu, min %lu, max %lu\n",
                perf_total.n_kernel_exec_cycles, perf_total.n_kernel_exec_cycles_avg, perf_total.n_kernel_exec_cycles_min, perf_total.n_kernel_exec_cycles_max);
    std::printf("iterations kernel execution.: total %.1f, avg %.1f, min %.1f, max %.1f\n",
                perf_total.n_kernel_exec_iters, perf_total.n_kernel_exec_iters_avg, perf_total.n_kernel_exec_iters_min, perf_total.n_kernel_exec_iters_max);
    std::printf("time memory dev2host........: total %8.6f s, avg %8.6f s, min %8.6f s, max %8.6f s\n",
                perf_total.s_mem_d2h, perf_total.s_mem_d2h_avg, perf_total.s_mem_d2h_min, perf_total.s_mem_d2h_max);
    std::printf("time solve..................: total %8.6f s, avg %8.6f s, min %8.6f s, max %8.6f s\n",
                perf_total.s_solve, perf_total.s_solve_avg, perf_total.s_solve_min, perf_total.s_solve_max);
    std::printf("time postprocess............: total %8.6f s, avg %8.6f s, min %8.6f s, max %8.6f s\n",
                perf_total.s_postprocess, perf_total.s_postprocess_avg, perf_total.s_postprocess_min, perf_total.s_postprocess_max);
    std::printf("converged...................: %u/%u, with iter>%d=%u, overflow=%u\n",
                perf_total.n_converged, num_data_points, maxit, conv_iter, conv_ovf);
    std::printf("-----------------------\n");
} //end generate_statistics()


#define INSTANTIATE_BDA_FUNCTIONS(n)                                                          \
template FpgaSolverBackend<n>::FpgaSolverBackend(std::string, int, int, double, ILUReorder);  \

INSTANTIATE_BDA_FUNCTIONS(1);
INSTANTIATE_BDA_FUNCTIONS(2);
INSTANTIATE_BDA_FUNCTIONS(3);
INSTANTIATE_BDA_FUNCTIONS(4);

#undef INSTANTIATE_BDA_FUNCTIONS

} //namespace bda
