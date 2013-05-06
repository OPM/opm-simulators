/*
  Copyright 2013 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_AUTODIFFHELPERS_HEADER_INCLUDED
#define OPM_AUTODIFFHELPERS_HEADER_INCLUDED

#include "AutoDiffBlock.hpp"
#include <opm/core/grid.h>


// -------------------- class HelperOps --------------------

/// Contains vectors and sparse matrices that represent subsets or
/// operations on (AD or regular) vectors of data.
struct HelperOps
{
    typedef AutoDiff::ForwardBlock<double>::M M;
    typedef AutoDiff::ForwardBlock<double>::V V;

    /// A list of internal faces.
    typedef Eigen::Array<int, Eigen::Dynamic, 1> IFaces;
    IFaces internal_faces;

    /// Extract for each face the difference of its adjacent cells'values.
    M ngrad;
    /// Extract for each face the average of its adjacent cells' values.
    M caver;
    /// Extract for each cell the sum of its adjacent faces' (signed) values.
    M div;

    /// Constructs all helper vectors and matrices.
    HelperOps(const UnstructuredGrid& grid)
    {
        const int nc = grid.number_of_cells;
        const int nf = grid.number_of_faces;
        // Define some neighbourhood-derived helper arrays.
        typedef Eigen::Array<int, Eigen::Dynamic, 1> OneColInt;
        typedef Eigen::Array<bool, Eigen::Dynamic, 1> OneColBool;
        typedef Eigen::Array<int, Eigen::Dynamic, 2, Eigen::RowMajor> TwoColInt;
        typedef Eigen::Array<bool, Eigen::Dynamic, 2, Eigen::RowMajor> TwoColBool;
        TwoColInt nb = Eigen::Map<TwoColInt>(grid.face_cells, nf, 2);
        // std::cout << "nb = \n" << nb << std::endl;
        TwoColBool nbib = nb >= 0;
        OneColBool ifaces = nbib.rowwise().all();
        const int num_internal = ifaces.cast<int>().sum();
        // std::cout << num_internal << " internal faces." << std::endl;
        TwoColInt nbi(num_internal, 2);
        internal_faces.resize(num_internal);
        int fi = 0;
        for (int f = 0; f < nf; ++f) {
            if (ifaces[f]) {
                internal_faces[fi] = f;
                nbi.row(fi) = nb.row(f);
                ++fi;
            }
        }
        // std::cout << "nbi = \n" << nbi << std::endl;
        // Create matrices.
        ngrad.resize(num_internal, nc);
        caver.resize(num_internal, nc);
        typedef Eigen::Triplet<double> Tri;
        std::vector<Tri> ngrad_tri;
        std::vector<Tri> caver_tri;
        ngrad_tri.reserve(2*num_internal);
        caver_tri.reserve(2*num_internal);
        for (int i = 0; i < num_internal; ++i) {
            ngrad_tri.emplace_back(i, nbi(i,0), 1.0);
            ngrad_tri.emplace_back(i, nbi(i,1), -1.0);
            caver_tri.emplace_back(i, nbi(i,0), 0.5);
            caver_tri.emplace_back(i, nbi(i,1), 0.5);
        }
        ngrad.setFromTriplets(ngrad_tri.begin(), ngrad_tri.end());
        caver.setFromTriplets(caver_tri.begin(), caver_tri.end());
        div = ngrad.transpose();
    }
};



// -------------------- debugger output helpers --------------------


#if !defined(NDEBUG)
#include <cstdio>
#endif  // !defined(NDEBUG)

namespace {
#if !defined(NDEBUG)
    void
    printSparseMatrix(const Eigen::SparseMatrix<double>& A,
                      std::FILE*                         fp)
    {
        typedef Eigen::SparseMatrix<double>::Index Index;

        const Index osize = A.outerSize();

        for (Index k = 0; k < osize; ++k) {
            for (Eigen::SparseMatrix<double>::InnerIterator
                     i(A, k); i ; ++i) {
                std::fprintf(fp, "%lu %lu %26.18e\n",
                             static_cast<unsigned long>(i.row() + 1),
                             static_cast<unsigned long>(i.col() + 1),
                             i.value());
            }
        }
    }

    void
    printSparseMatrix(const Eigen::SparseMatrix<double>& A ,
                      const char* const                  fn)
    {
        std::FILE* fp;

        fp = std::fopen(fn, "w");
        if (fp != 0) {
            printSparseMatrix(A, fp);
        }

        std::fclose(fp);
    }
#endif  // !defined(NDEBUG)



// -------------------- upwinding helper class --------------------


    template <typename Scalar>
    class UpwindSelectorTotalFlux {
    public:
        typedef AutoDiff::ForwardBlock<Scalar> ADB;

        UpwindSelectorTotalFlux(const UnstructuredGrid& g,
                                const HelperOps&        h,
                                const typename ADB::V&  ifaceflux)
        {
            typedef HelperOps::IFaces::Index IFIndex;
            const IFIndex nif = h.internal_faces.size();

            // Define selector structure.
            typedef typename Eigen::Triplet<Scalar> Triplet;
            std::vector<Triplet> s;  s.reserve(nif);
            for (IFIndex iface = 0; iface < nif; ++iface) {
                const int f  = h.internal_faces[iface];
                const int c1 = g.face_cells[2*f + 0];
                const int c2 = g.face_cells[2*f + 1];

                assert ((c1 >= 0) && (c2 >= 0));

                // Select upwind cell.
                const int c = (ifaceflux[iface] >= 0) ? c1 : c2;

                s.push_back(Triplet(iface, c, Scalar(1)));
            }

            // Assemble explicit selector operator.
            select_.resize(nif, g.number_of_cells);
            select_.setFromTriplets(s.begin(), s.end());
        }

        // Upwind selection in absence of counter-current flow (i.e.,
        // without effects of gravity and/or capillary pressure).
        std::vector<ADB>
        select(const std::vector<ADB>& xc) const
        {

            // Apply selector.
            //
            // Absence of counter-current flow means that the same
            // selector applies to all quantities, 'x', defined per
            // cell.
            std::vector<ADB> xf;  xf.reserve(xc.size());
            for (typename std::vector<ADB>::const_iterator
                     b = xc.begin(), e = xc.end(); b != e; ++b)
            {
                xf.push_back(select_ * (*b));
            }

            return xf;
        }

    private:
        typename ADB::M select_;
    };
}


#endif // OPM_AUTODIFFHELPERS_HEADER_INCLUDED
