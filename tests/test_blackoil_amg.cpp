/*
  Copyright 2017 Dr. Blatt - HPC-Simulation-Software & Services

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
#if HAVE_MPI
#define BOOST_TEST_MODULE BlackoilAmgTest
#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>
#include <opm/simulators/linalg/PreconditionerFactory.hpp>
#include <opm/simulators/linalg/FlexibleSolver.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/linalg/getQuasiImpesWeights.hpp>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/unused.hh>
#include <dune/common/parallel/indexset.hh>
#include <dune/common/parallel/plocalindex.hh>
#include <dune/common/parallel/communication.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/schwarz.hh>

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

typedef Dune::OwnerOverlapCopyAttributeSet GridAttributes;
typedef GridAttributes::AttributeSet GridFlag;
typedef Dune::ParallelLocalIndex<GridFlag> LocalIndex;

template<class M, class G, class L, int n>
void setupPattern(int N, M& mat, Dune::ParallelIndexSet<G,L,n>& indices, int overlapStart, int overlapEnd,
                  int start, int end);

template<class M>
void fillValues(int N, M& mat, int overlapStart, int overlapEnd, int start, int end);


template<class M, class G, class L, class C, int n>
M setupAnisotropic2d(int N, Dune::ParallelIndexSet<G,L,n>& indices,
                     const C& p, int *nout, typename M::block_type::value_type eps=1.0);


template<class M, class G, class L, int s>
void setupPattern(int N, M& mat, Dune::ParallelIndexSet<G,L,s>& indices, int overlapStart, int overlapEnd,
                  int start, int end)
{
    int n = overlapEnd - overlapStart;

    typename M::CreateIterator iter = mat.createbegin();
    indices.beginResize();

    for(int j=0; j < N; j++)
        for(int i=overlapStart; i < overlapEnd; i++, ++iter) {
            int global = j*N+i;
            GridFlag flag = GridAttributes::owner;
            bool isPublic = false;

            if((i<start && i > 0) || (i>= end && i < N-1))
                flag=GridAttributes::copy;

            if(i<start+1 || i>= end-1) {
                isPublic = true;
                indices.add(global, LocalIndex(iter.index(), flag, isPublic));
            }


            iter.insert(iter.index());

            // i direction
            if(i > overlapStart )
                // We have a left neighbour
                iter.insert(iter.index()-1);

            if(i < overlapEnd-1)
                // We have a rigt neighbour
                iter.insert(iter.index()+1);

            // j direction
            // Overlap is a dirichlet border, discard neighbours
            if(flag != GridAttributes::copy) {
                if(j>0)
                    // lower neighbour
                    iter.insert(iter.index()-n);
                if(j<N-1)
                    // upper neighbour
                    iter.insert(iter.index()+n);
            }
        }
    indices.endResize();
}

template<class M, class T>
void fillValues([[maybe_unused]] int N, M& mat, int overlapStart, int overlapEnd, int start, int end, T eps)
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

    int n = overlapEnd-overlapStart;
    typedef typename M::ColIterator ColIterator;
    typedef typename M::RowIterator RowIterator;

    for (RowIterator i = mat.begin(); i != mat.end(); ++i) {
        // calculate coordinate
        int y = i.index() / n;
        int x = overlapStart + i.index() - y * n;

        ColIterator endi = (*i).end();

        if(x<start || x >= end) { // || x==0 || x==N-1 || y==0 || y==N-1){
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

template<class V, class G, class L, int s>
void setBoundary(V& lhs, V& rhs, const G& n, Dune::ParallelIndexSet<G,L,s>& indices)
{
    typedef typename Dune::ParallelIndexSet<G,L,s>::const_iterator Iter;
    for(Iter i=indices.begin(); i != indices.end(); ++i) {
        G x = i->global()/n;
        G y = i->global()%n;

        if(x==0 || y ==0 || x==n-1 || y==n-1) {
            //double h = 1.0 / ((double) (n-1));
            lhs[i->local()]=rhs[i->local()]=0; //((double)x)*((double)y)*h*h;
        }
    }
}

template<class V, class G>
void setBoundary(V& lhs, V& rhs, const G& N)
{
    typedef typename V::block_type Block;
    typedef typename Block::value_type T;
    for(int j=0; j < N; ++j)
        for(int i=0; i < N; i++)
            if(i==0 || j ==0 || i==N-1 || j==N-1) {
                T h = 1.0 / ((double) (N-1));
                T x, y;
                if(i==N-1)
                    x=1;
                else
                    x = ((T) i)*h;
                if(j==N-1)
                    y = 1;
                else
                    y = ((T) j)*h;

                lhs[j*N+i]=rhs[j*N+i]=0; //x*y;
            }
}

template<class M, class G, class L, class C, int s>
M setupAnisotropic2d(int N, Dune::ParallelIndexSet<G,L,s>& indices, const C& p, int *nout, typename M::block_type::value_type eps)
{
    int procs=p.size(), rank=p.rank();

    typedef M BCRSMat;

    // calculate size of local matrix in the distributed direction
    int start, end, overlapStart, overlapEnd;

    int n = N/procs; // number of unknowns per process
    int bigger = N%procs; // number of process with n+1 unknows

    // Compute owner region
    if(rank<bigger) {
        start = rank*(n+1);
        end   = start+(n+1);
    }else{
        start = bigger*(n+1) + (rank-bigger) * n;
        end   = start+n;
    }

    // Compute overlap region
    if(start>0)
        overlapStart = start - 1;
    else
        overlapStart = start;

    if(end<N)
        overlapEnd = end + 1;
    else
        overlapEnd = end;

    int noKnown = overlapEnd-overlapStart;

    *nout = noKnown;

    BCRSMat mat(noKnown*N, noKnown*N, noKnown*N*5, BCRSMat::row_wise);

    setupPattern(N, mat, indices, overlapStart, overlapEnd, start, end);
    fillValues(N, mat, overlapStart, overlapEnd, start, end, eps);

    //  Dune::printmatrix(std::cout,mat,"aniso 2d","row",9,1);

    return mat;
}


//BOOST_AUTO_TEST_CASE(runBlackoilAmgLaplace)
void runBlackoilAmgLaplace()
{
    constexpr int BS=2, N=100;
    typedef Dune::FieldMatrix<double,BS,BS> MatrixBlock;
    typedef Dune::BCRSMatrix<MatrixBlock> BCRSMat;
    typedef Dune::FieldVector<double,BS> VectorBlock;
    typedef Dune::BlockVector<VectorBlock> Vector;
    typedef int GlobalId;
    typedef Dune::OwnerOverlapCopyCommunication<GlobalId> Communication;
    typedef Dune::OverlappingSchwarzOperator<BCRSMat,Vector,Vector,Communication> Operator;

    const auto& ccomm = Dune::MPIHelper::getCommunication();

    Communication comm(ccomm);
    int n=0;
    BCRSMat mat = setupAnisotropic2d<BCRSMat>(N, comm.indexSet(), comm.communicator(), &n, 1);


    comm.remoteIndices().template rebuild<false>();

    Vector b(mat.N()), x(mat.M());

    b=0;
    x=100;
    setBoundary(x, b, N, comm.indexSet());

    Operator fop(mat, comm);
    Dune::OverlappingSchwarzScalarProduct<Vector,Communication> sp(comm);
    Dune::InverseOperatorResult r;

    using namespace std::string_literals;
    Opm::PropertyTree prm;
    prm.put("type", "amg"s);
    std::function<Vector()> weights = [&mat]() {
        return Opm::Amg::getQuasiImpesWeights<BCRSMat, Vector>(mat, 0, false);
    };
    auto amg = Opm::PreconditionerFactory<Operator, Communication>::create(fop, prm, weights, comm);

    Dune::CGSolver<Vector> amgCG(fop, sp, *amg, 10e-8, 300, (ccomm.rank()==0) ? 2 : 0);

    amgCG.apply(x,b,r);

}

bool init_unit_test_func()
{
    return true;
}

int main(int argc, char** argv)
{
    [[maybe_unused]] const auto& helper = Dune::MPIHelper::instance(argc, argv);
    boost::unit_test::unit_test_main(&init_unit_test_func,
                                     argc, argv);
}
#else
int main () { return 0; }
#endif // #if HAVE_MPI
