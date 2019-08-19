// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef MATRIXMARKET_EWOMS_HH
#define MATRIXMARKET_EWOMS_HH
#include <dune/istl/matrixmarket.hh>
#include <ewoms/linear/matrixblock.hh>
namespace Dune
{


  namespace MatrixMarketImpl
  {


    template<typename T, typename A, int i, int j>
    struct mm_header_printer<BCRSMatrix<Ewoms::MatrixBlock<T,i,j>,A> >
    {
      static void print(std::ostream& os)
      {
        os<<"%%MatrixMarket matrix coordinate ";
        os<<mm_numeric_type<T>::str()<<" general"<<std::endl;
      }
    };


    template<typename T, int i, int j>
    struct mm_header_printer<Ewoms::MatrixBlock<T,i,j> >
    {
      static void print(std::ostream& os)
      {
        os<<"%%MatrixMarket matrix array ";
        os<<mm_numeric_type<T>::str()<<" general"<<std::endl;
      }
    };

    /**
     * @brief Metaprogram for writing the ISTL block
     * structure header.
     *
     * Member function mm_block_structure_header::print(os, mat) writes
     * the corresponding header to the specified ostream.
     * @tparam The type of the matrix to generate the header for.
     */

    template<typename T, typename A, int i, int j>
    struct mm_block_structure_header<BCRSMatrix<Ewoms::MatrixBlock<T,i,j>,A> >
    {
      typedef BCRSMatrix<Ewoms::MatrixBlock<T,i,j>,A> M;

      static void print(std::ostream& os, const M&)
      {
        os<<"% ISTL_STRUCT blocked ";
        os<<i<<" "<<j<<std::endl;
      }
    };


    template<typename T, int i, int j>
    struct mm_block_structure_header<Ewoms::MatrixBlock<T,i,j> >
    {
      typedef Ewoms::MatrixBlock<T,i,j> M;

      static void print(std::ostream& os, const M& m)
      {}
    };



    template<typename T, typename A, int brows, int bcols, typename D>
    void readSparseEntries(Dune::BCRSMatrix<Ewoms::MatrixBlock<T,brows,bcols>,A>& matrix,
                           std::istream& file, std::size_t entries,
                           const MMHeader& mmHeader, const D&)
    {
      typedef Dune::BCRSMatrix<Ewoms::MatrixBlock<T,brows,bcols>,A> Matrix;
      // First path
      // store entries together with column index in a speparate
      // data structure
      std::vector<std::set<IndexData<D> > > rows(matrix.N()*brows);

      for(; entries>0; --entries) {
        std::size_t row;
        IndexData<D> data;
        skipComments(file);
        file>>row;
        --row; // Index was 1 based.
        assert(row/bcols<matrix.N());
        file>>data;
        assert(data.index/bcols<matrix.M());
        rows[row].insert(data);
      }

      // TODO extend to capture the nongeneral cases.
      if(mmHeader.structure!= general)
        DUNE_THROW(Dune::NotImplemented, "Only general is supported right now!");

      // Setup the matrix sparsity pattern
      int nnz=0;
      for(typename Matrix::CreateIterator iter=matrix.createbegin();
          iter!= matrix.createend(); ++iter)
      {
        for(std::size_t brow=iter.index()*brows, browend=iter.index()*brows+brows;
            brow<browend; ++brow)
        {
          typedef typename std::set<IndexData<D> >::const_iterator Siter;
          for(Siter siter=rows[brow].begin(), send=rows[brow].end();
              siter != send; ++siter, ++nnz)
            iter.insert(siter->index/bcols);
        }
      }

      //Set the matrix values
      matrix=0;

      MatrixValuesSetter<D,brows,bcols> Setter;

      Setter(rows, matrix);
    }
  } // end namespace MatrixMarketImpl


  template<typename T, typename A, int brows, int bcols>
  void readMatrixMarket(Dune::BCRSMatrix<Ewoms::MatrixBlock<T,brows,bcols>,A>& matrix,
                        std::istream& istr)
  {
    using namespace MatrixMarketImpl;

    MMHeader header;
    if(!readMatrixMarketBanner(istr, header)) {
      std::cerr << "First line was not a correct Matrix Market banner. Using default:\n"
                << "%%MatrixMarket matrix coordinate real general"<<std::endl;
      // Go to the beginning of the file
      istr.clear() ;
      istr.seekg(0, std::ios::beg);
    }
    skipComments(istr);

    std::size_t rows, cols, entries;

    if(lineFeed(istr))
      throw MatrixMarketFormatError();

    istr >> rows;

    if(lineFeed(istr))
      throw MatrixMarketFormatError();
    istr >> cols;

    if(lineFeed(istr))
      throw MatrixMarketFormatError();

    istr >>entries;

    std::size_t nnz, blockrows, blockcols;

    std::tie(blockrows, blockcols, nnz) = calculateNNZ<brows, bcols>(rows, cols, entries, header);

    istr.ignore(std::numeric_limits<std::streamsize>::max(),'\n');


    matrix.setSize(blockrows, blockcols);
    matrix.setBuildMode(Dune::BCRSMatrix<Ewoms::MatrixBlock<T,brows,bcols>,A>::row_wise);

    if(header.type==array_type)
      DUNE_THROW(Dune::NotImplemented, "Array format currently not supported for matrices!");

    readSparseEntries(matrix, istr, entries, header, NumericWrapper<T>());
  }



  template<typename B, int i, int j, typename A>
  struct mm_multipliers<BCRSMatrix<Ewoms::MatrixBlock<B,i,j>,A> >
  {
    enum {
      rows = i,
      cols = j
    };
  };

  template<typename B, int i, int j>
  void mm_print_entry(const Ewoms::MatrixBlock<B,i,j>& entry,
                      typename Ewoms::MatrixBlock<B,i,j>::size_type rowidx,
                      typename Ewoms::MatrixBlock<B,i,j>::size_type colidx,
                      std::ostream& ostr)
  {
    typedef typename Ewoms::MatrixBlock<B,i,j>::const_iterator riterator;
    typedef typename Ewoms::MatrixBlock<B,i,j>::ConstColIterator citerator;
    for(riterator row=entry.begin(); row != entry.end(); ++row, ++rowidx) {
      int coli=colidx;
      for(citerator col = row->begin(); col != row->end(); ++col, ++coli)
        ostr<< rowidx<<" "<<coli<<" "<<*col<<std::endl;
    }
  }

  template<typename T, int n, int m>
  std::size_t countEntries(const Ewoms::MatrixBlock<T,n, m>& vector)
  {
    return n*m;
  }
  
  /** @} */
}

#endif // MATRIXMARKET_EWOMS_HH
