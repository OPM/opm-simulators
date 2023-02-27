#include<config.h>

#define BOOST_TEST_MODULE MILU0Test

#include<vector>
#include<memory>

#include<dune/istl/bcrsmatrix.hh>
#include<dune/istl/bvector.hh>
#include<dune/common/version.hh>
#include<dune/common/fmatrix.hh>
#include<dune/common/fvector.hh>
#include<opm/simulators/linalg/ParallelOverlappingILU0.hpp>

#include <opm/common/ErrorMacros.hpp>

#include <boost/test/unit_test.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION / 100000 == 1 && BOOST_VERSION / 100 % 1000 < 71
#include <boost/test/floating_point_comparison.hpp>
#else
#include <boost/test/tools/floating_point_comparison.hpp>
#endif

template<class M>
void test_milu0(M& A)
{
    auto ILU = A;
    
    std::shared_ptr<std::vector<typename M::block_type>> diagonal = nullptr;
    diagonal.reset(new std::vector<typename M::block_type>());

    Opm::detail::milu0_decomposition(ILU, diagonal.get());
#ifdef DEBUG
    if ( A.N() < 11)
    {
        Dune::printmatrix(std::cout, ILU, "ILU", "row");
        std::cout << "Diagonal: ";

        for (const auto& d : *diagonal)
            std::cout << d << " ";
        std::cout<<std::endl;
    }
#endif
    Dune::BlockVector<Dune::FieldVector<typename M::field_type, M::block_type::rows>> e(A.N()), x1(A.N()), x2(A.N()), t(A.N());
    e = 1;
    A.mv(e, x1);
    // Compute L*U*e;
    // t=Ue
    
    for ( auto irow = ILU.begin(), iend = ILU.end(); irow != iend; ++irow)
    {
        auto i = irow.index();
        (*diagonal)[i].mv(e[0], t[i]);
        auto col = irow->find(i);
        auto colend = irow->end();

        if ( col == colend )
        {
            OPM_THROW(std::logic_error,
                      "Missing diagonal entry for row " + std::to_string(irow.index()));
        }
        for (++col ;col != colend; ++col)
        {
            col->umv(e[0], t[i]);
        }
    }
            
    for ( auto irow = ILU.begin(), iend = ILU.end(); irow != iend; ++irow)
    {
        auto i = irow.index();
        x2[i] = 0;
        for (auto col = irow->begin(); col.index() < irow.index(); ++col)
        {
            col->umv(t[col.index()], x2[i]);
        }
        x2[i] += t[i];
    }
    auto diff = x2;
    diff-=x1;
    for ( std::size_t i = 0, end = A.N(); i < end; ++i)
    {
        auto point_difference = diff[i].two_norm();
#ifdef DEBUG
        std::cout<<"index "<<i<<" size "<<diff.size()<<" difference"<<point_difference<<std::endl;
#endif
        BOOST_CHECK(point_difference < 1e-12);
    }

    // Test that (LU)^-1Ae=e
    A.mv(e, x1);
#if DUNE_VERSION_GTE(DUNE_ISTL, 2, 8)
    Dune::ILU::blockILUBacksolve(ILU, x2, x1);
#else
    bilu_backsolve(ILU, x2, x1);
#endif
    diff = x2;
    diff -= e;

    for ( std::size_t i = 0, end = A.N(); i < end; ++i)
    {
#ifdef DEBUG
        auto point_difference = diff[i].two_norm();
        std::cout<<"index "<<i<<" size "<<diff.size()<<" point_difference "<<point_difference<<std::endl;
#endif
        BOOST_CHECK_CLOSE(x2[i].two_norm(), e[i].two_norm(), 1e-12);
    }
}

template<class B, class Alloc>
void setupSparsityPattern(Dune::BCRSMatrix<B,Alloc>& A, int N)
{
  typedef typename Dune::BCRSMatrix<B,Alloc> Matrix;
  A.setSize(N*N, N*N, N*N*5);
  A.setBuildMode(Matrix::row_wise);

  for (typename Dune::BCRSMatrix<B,Alloc>::CreateIterator i = A.createbegin(); i != A.createend(); ++i) {
    int x = i.index()%N; // x coordinate in the 2d field
    int y = i.index()/N;  // y coordinate in the 2d field

    if(y>0)
      // insert lower neighbour
      i.insert(i.index()-N);
    if(x>0)
      // insert left neighbour
      i.insert(i.index()-1);

    // insert diagonal value
    i.insert(i.index());

    if(x<N-1)
      //insert right neighbour
      i.insert(i.index()+1);
    if(y<N-1)
      // insert upper neighbour
      i.insert(i.index()+N);
  }
}


template<class B, class Alloc>
void setupLaplacian(Dune::BCRSMatrix<B,Alloc>& A, int N)
{
  typedef typename Dune::BCRSMatrix<B,Alloc>::field_type FieldType;

  setupSparsityPattern(A,N);

  B diagonal(static_cast<FieldType>(0)), bone(static_cast<FieldType>(0)),
  beps(static_cast<FieldType>(0));
  for(typename B::RowIterator b = diagonal.begin(); b !=  diagonal.end(); ++b)
    b->operator[](b.index())=4;


  for(typename B::RowIterator b=bone.begin(); b !=  bone.end(); ++b)
    b->operator[](b.index())=-1.0;


  for (typename Dune::BCRSMatrix<B,Alloc>::RowIterator i = A.begin(); i != A.end(); ++i) {
    int x = i.index()%N; // x coordinate in the 2d field
    int y = i.index()/N;  // y coordinate in the 2d field
    
    i->operator[](i.index())=diagonal;
    
    if(y>0)
        i->operator[](i.index()-N)=bone;

    if(y<N-1)
        i->operator[](i.index()+N)=bone;

    if(x>0)
        i->operator[](i.index()-1)=bone;

    if(x < N-1)
        i->operator[](i.index()+1)=bone;
  }
}

template<int bsize>
void test()
{
    std::size_t N = 32;
    Dune::BCRSMatrix<Dune::FieldMatrix<double, bsize, bsize> > A;
    setupLaplacian(A, N);
    test_milu0(A);
#ifdef DEBUG
    std::cout<< "Tested block size "<< bsize<<std::endl;
#endif
}

BOOST_AUTO_TEST_CASE(MILULaplace1)
{
    test<1>();
}

BOOST_AUTO_TEST_CASE(MILULaplace2)
{
    test<2>();
}
BOOST_AUTO_TEST_CASE(MILULaplace3)
{
    test<3>();
}
BOOST_AUTO_TEST_CASE(MILULaplace4)
{
    test<4>();
}
