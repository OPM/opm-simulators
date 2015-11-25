/*
  Copyright 2015 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2015 Statoil AS

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
//#define BOOST_TEST_MODULE SendReceiveCommunicatorTest
#include <boost/test/unit_test.hpp>

#include <opm/autodiff/ParallelILU0.hpp>

#include <dune/istl/preconditioners.hh>

#include <dune/common/version.hh>

// MPI header
#if HAVE_MPI
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 3)

#include <dune/common/parallel/mpicollectivecommunication.hh>

#else

#include <dune/common/mpicollectivecommunication.hh>

#endif

#include <mpi.h>

#endif

#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 3)

#include <dune/common/parallel/collectivecommunication.hh>
#include <dune/common/parallel/mpihelper.hh>

#else

#include <dune/common/collectivecommunication.hh>
#include <dune/common/mpihelper.hh>

#endif


#include <random>
#include <algorithm>

void testIdMapping()
{
#if HAVE_MPI
    Dune::CollectiveCommunication<MPI_Comm>  cc(MPI_COMM_WORLD);
    int rank  = cc.rank();
    Opm::Detail::IdProcessMap map(cc);
    for( int i = 0 ; i < cc.size(); ++i)
    {
        BOOST_CHECK( (i < rank ) == map.compareWithOtherLabel(i, std::less<int>()));
        BOOST_CHECK( (i > rank ) == map.compareWithOtherLabel(i, std::greater<int>()));
    }
#endif
}

#if HAVE_MPI

typedef Dune::OwnerOverlapCopyAttributeSet GridAttributes;
typedef GridAttributes::AttributeSet GridFlag;
typedef Dune::ParallelLocalIndex<GridFlag> LocalIndex;
/*
  template<class M, class G, class L, int n>
  void setupPattern(int N, M& mat, Dune::ParallelIndexSet<G,L,n>& indices, int overlapStart, int overlapEnd,
  int start, int end);

  template<class M>
  void fillValues(int N, M& mat, int overlapStart, int overlapEnd, int start, int end);


  template<class M, class G, class L, class C, int n>
  M setupAnisotropic2d(int N, Dune::ParallelIndexSet<G,L,n>& indices,
  const Dune::CollectiveCommunication<C>& p, int *nout, typename M::block_type::value_type eps=1.0);
*/
// A custom mpi error handler for more easy debugging (e.g. using catch throw in gdb)
class MPIError {
public:
    /** @brief Constructor. */
    MPIError(std::string s, int e) : errorstring(s), errorcode(e){}
    /** @brief The error string. */
    std::string errorstring;
    /** @brief The mpi error code. */
    int errorcode;
};

void MPI_err_handler(MPI_Comm *, int *err_code, ...){
    char *err_string=new char[MPI_MAX_ERROR_STRING];
    int err_length;
    MPI_Error_string(*err_code, err_string, &err_length);
    std::string s(err_string, err_length);
    std::cerr << "An MPI Error ocurred:"<<std::endl<<s<<std::endl;
    delete[] err_string;
    throw MPIError(s, *err_code);
}

struct OneDRegionInformation
{
    OneDRegionInformation(int n, const Dune::CollectiveCommunication<MPI_Comm>& p)
        : N(p.size() * n)
    {
        std::size_t rank=p.rank();
        // partition start
        start = rank * n;
        end   = start+n;
        // Compute overlap region
        if(start>0)
            overlap_start = start - 1;
        else
            overlap_start = start;

        if(end<N)
            overlap_end = end + 1;
        else
            overlap_end = end;
        no_unknowns = overlap_end - overlap_start;
    }
    bool hasStartOverlap() const
    {
        return overlap_start < start;
    }
    bool hasEndOverlap() const
    {
        return end < overlap_end;
    }

    std::size_t N;
    std::size_t start, end;
    std::size_t overlap_start, overlap_end;
    std::size_t no_unknowns;
};


template<class M, class G, class L, int s>
void setupPattern(std::size_t y_size, M& mat, Dune::ParallelIndexSet<G,L,s>& indices,
                  const OneDRegionInformation& region)
{
    //int n = overlapEnd - overlapStart;

    typename M::CreateIterator iter = mat.createbegin();
    indices.beginResize();

    for(std::size_t j = 0; j < y_size; j++)
    {
        for(std::size_t i = region.overlap_start; i < region.overlap_end; i++, ++iter)
        {
            int global = j * region.N + i;
            GridFlag flag = GridAttributes::owner;
            bool isPublic = false;

            if( (i < region.start) || (i >= region.end) )
                flag=GridAttributes::copy;

            if( (i < region.start+1 && region.hasStartOverlap()) ||
                (i >= region.end-1 && region.hasEndOverlap()) ) {
                isPublic = true;
                indices.add(global, LocalIndex(iter.index(), flag, isPublic));
            }


            iter.insert(iter.index());

            // i direction
            if(i > region.overlap_start )
                // We have a left neighbour
                iter.insert(iter.index()-1);

            if(i < region.overlap_end-1)
                // We have a rigt neighbour
                iter.insert(iter.index()+1);

            // j direction
            if(j > 0)
                // lower neighbour
                iter.insert(iter.index() - region.no_unknowns);
            if(j < y_size-1)
                // upper neighbour
                iter.insert(iter.index() + region.no_unknowns);
        }
    }
    indices.endResize();
}

template<class M, class T>
void fillValues(M& mat, const OneDRegionInformation& region, T eps)
{
    typedef typename M::block_type Block;
    Block dval(0), bone(0), bmone(0), beps(0);

    for(typename Block::RowIterator b = dval.begin(); b !=  dval.end(); ++b)
        b->operator[](b.index())=2.0+2.0*eps;

    for(typename Block::RowIterator b=bone.begin(); b !=  bone.end(); ++b)
        b->operator[](b.index())=1.0;

    for(typename Block::RowIterator b=bmone.begin(); b !=  bmone.end(); ++b)
        b->operator[](b.index())=-1.0;

    for(typename Block::RowIterator b=beps.begin(); b !=  beps.end(); ++b)
        b->operator[](b.index())=-eps;


    typedef typename M::ColIterator ColIterator;
    typedef typename M::RowIterator RowIterator;

    for (RowIterator i = mat.begin(); i != mat.end(); ++i) {
        // calculate coordinate
        std::size_t y = i.index() / region.no_unknowns;
        std::size_t x = region.overlap_start + i.index() - y * region.no_unknowns;

        ColIterator endi = (*i).end();

        if(x<region.start || x >= region.end) {
            // overlap node is dirichlet border
            ColIterator j = (*i).begin();

            for(; j.index() < i.index(); ++j)
                *j=0;

            *j = bone;

            for(++j; j != endi; ++j)
                *j=0;
        }else{
            for(ColIterator j = (*i).begin(); j != endi; ++j)
                if(j.index() == i.index())
                    *j=dval;
                else if(j.index()+1==i.index() || j.index()-1==i.index())
                    *j=beps;
                else
                    *j=bmone;
        }
    }
}


template<class M, class G, class L, int s>
std::unique_ptr<M> setupAnisotropic2d(int n, Dune::ParallelIndexSet<G,L,s>& indices,
                                      const OneDRegionInformation& region, typename M::block_type::value_type eps)
{
    typedef M BCRSMat;

    // calculate size of local matrix in the distributed direction
    std::size_t no_unknowns = region.no_unknowns;

    std::unique_ptr<BCRSMat> mat;
    mat.reset(new BCRSMat(no_unknowns*n, no_unknowns*n, no_unknowns*n*5, BCRSMat::row_wise));

    setupPattern(n, *mat, indices, region);
    fillValues(*mat, region, eps);
    return mat;
}

struct LocalExtractor
{
    template<class B>
    void operator()(B& b, const B& global_b, std::size_t, std::size_t)
    {
        b = global_b;
    }
};

struct ParallelGlobalLinearSystemPermutation
    : public Opm::Detail::LinearSystemPermutation
{
    ParallelGlobalLinearSystemPermutation(std::size_t n, std::size_t procs)
        : Opm::Detail::LinearSystemPermutation(n*procs*n)
    {
        std::size_t N = n*procs;
        // interior nodes
        std::size_t new_index = 0;
        for(std::size_t p = 0; p < procs; ++p)
        {
            for(std::size_t j=0; j< n; ++j)
            {
                std::size_t end   = p * n + n - 1;
                std::size_t start = p * n + 1;
                if ( p == 0) start=0;
                if ( p == procs - 1 ) end = (p + 1) * n;
                for(std::size_t i = start; i < end; ++i, ++new_index)
                {
                    row_permutation_[ j * N + i ] = new_index;
                }
            }
        }
        if ( n > 1)
        {
            // upper interface
            for(std::size_t p = 0; p < procs - 1; ++p)
            {
                std::size_t i = p * n + n - 1;
                for(std::size_t j=0; j< n; ++j, ++new_index)
                {
                    row_permutation_[ j * N + i ] = new_index;
                }
            }
        }
        // lower interface
        for(std::size_t p = 1; p < procs; ++p)
        {
            std::size_t i = p*n;
            for(std::size_t j=0; j< n; ++j, ++new_index)
            {
                row_permutation_[ j * N + i ] = new_index;
            }
        }
    }
};



struct GlobalLocalChecker
{
    GlobalLocalChecker(std::string str)
        : str_(str)
    {
        MPI_Comm_rank(MPI_COMM_WORLD, & rank_);
    }
    template<class B>
    void operator()(B& b, const B& global_b, std::size_t local, std::size_t global)
    {
        bool is_same = (b - global_b).two_norm() < 1e-5;
        BOOST_CHECK_MESSAGE(is_same, rank_<<": global index "<<global
                            <<" "<<b<<" (local) != "
                            <<global_b<<" (global), local index: "<<local
                            <<" ("<<str_<<")!");
    }
    std::string str_;
    int rank_;
};

template<class V, class T>
T processValuesFromGlobalVector(V& v, const V& global_v,
                                const OneDRegionInformation& region,
                                std::size_t n, T processor)
{
    for ( std::size_t y=0; y < n; ++y )
    {
        for ( std::size_t x = region.overlap_start; x < region.overlap_end; ++x )
        {
            std::size_t local = y * region.no_unknowns + x - region.overlap_start;
            std::size_t global = y * region.N + x;
            processor(v[local], global_v[global], local, global);
        }
    }
    return processor;
}

template<int block_size>
void test_parallel_ilu0()
{
    typedef Dune::FieldMatrix<double, block_size, block_size> MatrixBlock;
    typedef Dune::BCRSMatrix<MatrixBlock> BCRSMat;
    typedef Dune::FieldVector<double, block_size> VectorBlock;
    typedef Dune::BlockVector<VectorBlock> Vector;
    typedef int GlobalId;
    typedef Dune::OwnerOverlapCopyCommunication<GlobalId> Communication;
    Communication comm(MPI_COMM_WORLD), self_comm(MPI_COMM_SELF);
    const int    N = 3;
    const int size = comm.communicator().size();
    const int global_unknowns = size * N * N;
    const OneDRegionInformation region(N, comm.communicator());
    const OneDRegionInformation global_region(size * N, self_comm.communicator());


    std::unique_ptr<BCRSMat> mat =
        setupAnisotropic2d<BCRSMat>(N, comm.indexSet(), region, 1);
    std::unique_ptr<BCRSMat> global_mat = setupAnisotropic2d<BCRSMat>(N, self_comm.indexSet(),
                                                                      global_region, 1);
    ParallelGlobalLinearSystemPermutation permutation(N, size);

    std::unique_ptr<BCRSMat> global_mat_permuted = permutation.createPermutedMatrix(*global_mat);

    comm.remoteIndices().template rebuild<false>();
    Opm::ParallelILU0<BCRSMat, Vector, Vector, Communication> pilu0(*mat, comm, 1.0);
    Dune::SeqILU0<BCRSMat, Vector, Vector> ilu0(*global_mat_permuted, 1.0);

    std::random_device rd;
    std::default_random_engine e1(rd());
    std::uniform_int_distribution<int> uniform_dist(1, 10);
    Vector global_x(global_unknowns), global_b(global_unknowns);
    Vector permuted_x(global_unknowns), permuted_b(global_unknowns);

    if( comm.communicator().rank() == 0 )
    {
        std::generate(global_b.begin(), global_b.end(),
                      [&] (){
                          return typename Vector::block_type
                              (static_cast<double>(uniform_dist(rd)));});
    }
    comm.communicator().broadcast(&(global_b[0]), global_b.size(), 0);
    global_x = 0.0;

    Vector x(mat->N()), b(mat->N());
    processValuesFromGlobalVector(b, global_b, region, N, LocalExtractor());
    processValuesFromGlobalVector(b, global_b, region, N, GlobalLocalChecker("extracted"));
    permutation.permutateOrder(global_b, permuted_b);

    ilu0.apply(permuted_x, permuted_b);
    pilu0.apply(x, b);
    permutation.permutateOrderBackwards(global_x, permuted_x);
    processValuesFromGlobalVector(x, global_x, region, N, GlobalLocalChecker("computed"));

}
void test_parallel_ilu0_1()
{
    test_parallel_ilu0<1>();
}

void test_parallel_ilu0_3()
{
    test_parallel_ilu0<3>();
}
#endif

bool init_function()
{
    boost::unit_test::framework::master_test_suite().
        add( BOOST_TEST_CASE(&testIdMapping) );
#if HAVE_MPI
    boost::unit_test::framework::master_test_suite().
        add( BOOST_TEST_CASE(&test_parallel_ilu0_1) );
    /*    boost::unit_test::framework::master_test_suite().
          add( BOOST_TEST_CASE(&test_parallel_ilu0_3) ); */
#endif
    return true;
}

int main(int argc, char** argv)
{
#if HAVE_MPI
    MPI_Init(&argc, &argv);
    int procs;
    MPI_Comm_size(MPI_COMM_WORLD, &procs);
    if ( procs == 1 )
    {
        std::cerr<<"This test should be run with at least two processes to make sense."<<std::endl;
    }
    MPI_Errhandler handler;
    MPI_Comm_create_errhandler(MPI_err_handler, &handler);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, handler);
#endif
    ::boost::unit_test::unit_test_main( &init_function, argc, argv );
#if HAVE_MPI
    MPI_Finalize();
#else
#warning "This file needs to compiled with MPI support to make any sense!"
#endif
}
