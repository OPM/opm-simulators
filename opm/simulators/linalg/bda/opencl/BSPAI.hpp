#ifndef BSPAI_HPP
#define BSPAI_HPP

#include <set>
#include <mutex>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/istl/bccsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/io.hh>
#include <dune/istl/spqr.hh>
#include <dune/istl/solver.hh>

#include <opm/simulators/linalg/bda/opencl/opencl.hpp>
#include <opm/simulators/linalg/bda/opencl/Preconditioner.hpp>

namespace Opm
{
namespace Accelerator
{

class BlockedMatrix;

/// This class implements a Blocked version of the Sparse Approximate Inverse (SPAI) preconditioner.
template <unsigned int block_size>
class BSPAI : public Preconditioner<block_size>
{
    typedef Preconditioner<block_size> Base;

    using Base::N;
    using Base::Nb;
    using Base::nnz;
    using Base::nnzb;
    using Base::verbosity;
    using Base::context;
    using Base::queue;
    using Base::events;
    using Base::err;

    using DuneMat = Dune::BCRSMatrix<Dune::FieldMatrix<double, block_size, block_size> >;
    using DuneVec = Dune::BlockVector<Dune::FieldMatrix<double, block_size, block_size> >;

private:
    int fill_in = 2;

    std::vector<int> submatPointers;
    std::vector<int> submatIndices;
    std::vector<std::vector<int> > submatValsPositions;
    std::vector<int> eyeBlockIndices;
    std::vector<int> rowPointers, colIndices;
    std::vector<int> spaiColPointers, spaiRowIndices;
    std::vector<double> nnzValues, spaiNnzValues;

    std::set<int> iset, jset;

    DuneVec rhs, sol;
    std::vector<DuneMat> submat;

    Dune::InverseOperatorResult res;
    Dune::SPQR<DuneMat> solver;

    std::once_flag ocl_init;
    cl::Buffer d_spaiColPointers, d_spaiRowIndices, d_spaiNnzValues;

    void buildIJSets(int tcol);
    void gatherSubmatIndices();

public:
    BSPAI(ILUReorder opencl_ilu_reorder, int verbosity);

    ~BSPAI();

    // set OpenCL variables
    void setOpencl(std::shared_ptr<cl::Context>& context, std::shared_ptr<cl::CommandQueue>& queue) override;

    // analysis, find reordering if specified
    bool analyze_matrix(BlockedMatrix *mat) override;

    // ilu_decomposition
    bool create_preconditioner(BlockedMatrix *mat) override;

    // apply preconditioner, x = prec(y)
    void apply(const cl::Buffer& d_y, cl::Buffer& d_x) override;

    int* getToOrder() override {return 0;}
    int* getFromOrder() override {return 0;}
    BlockedMatrix* getRMat() override {return 0;}
};
} // namespace Accelerator
} // namespace Opm

#endif
