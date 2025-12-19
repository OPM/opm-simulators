/*
  Copyright 2025 Equinor ASA

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

#ifndef OPM_EXPORT_SYSTEM_HEADER_INCLUDED
#define OPM_EXPORT_SYSTEM_HEADER_INCLUDED

#include <cstdio>
#include <memory>

namespace Opm
{

    // Forward declaring as they will be called from exportSystem().
    template <class IstlMatrix>
    void exportSparsity(const IstlMatrix& A, const char *path=".");
    template <class IstlMatrix>
    void exportNonzeros(const IstlMatrix& A, const char *tag="", const char *path=".");
    template <class GlobalEqVector>
    void exportVector(const GlobalEqVector& x, const char *tag="", const char *name="export/x");

    /*!
     * \brief Export blocks-sparse linear system.
     */
    template <class IstlMatrix, class GlobalEqVector>
    void
    exportSystem(const IstlMatrix& jacobian,
                 const GlobalEqVector& residual,
                 const bool export_sparsity,
                 const char* tag,
                 const char* path = "export")
    {
        // export sparsity only if requested
        if (export_sparsity) {
            exportSparsity(jacobian,path);
        }

        // export matrix
        exportNonzeros(jacobian,tag,path);

        // export residual
        constexpr size_t bufsize = 256;
        char name[bufsize];
        std::snprintf(name,bufsize,"%s/r",path);
        exportVector(residual,tag,name);
    }

    /*!
     * \brief Export block vector.
     */
    template <class GlobalEqVector>
    void exportVector(const GlobalEqVector& x, const char *tag, const char *name)
    {
        // assume double precision and contiguous data
        const double *data = &x[0][0];

        constexpr size_t bufsize = 512;
        char filename[bufsize];
        std::snprintf(filename,bufsize,"%s%s.f64",name,tag);
        FILE *out =fopen(filename,"w");
        std::fwrite(data, sizeof(double), x.dim(),out);
        std::fclose(out);
    }

    /*!
     * \brief Export nonzero blocks of jacobian block-sparse matrix
     */
    template <class IstlMatrix>
    void exportNonzeros(const IstlMatrix& A, const char *tag, const char *path)
    {
        // assume double precision and contiguous data
        const double *data = &A[0][0][0][0];
        size_t dim = A[0][0].N()*A[0][0].M()*A.nonzeroes();

        constexpr size_t bufsize = 256;
        char filename[bufsize];
        std::snprintf(filename,bufsize,"%s/data%s.f64",path,tag);
        FILE *out =fopen(filename,"w");
        std::fwrite(data, sizeof(double), dim,out);
        std::fclose(out);
    }

    /*!
     * \brief Export sparsity pattern of jacobian block-sparse matrix
     */
    template <class IstlMatrix>
    void exportSparsity(const IstlMatrix& A, const char *path)
    {
        //assemble csr graph
        auto rows = std::make_unique<int[]>(A.N()+1);
        auto cols = std::make_unique<int[]>(A.nonzeroes());

        int irow=0;
        int icol=0;
        rows[0]=0;
        for(auto row=A.begin(); row!=A.end(); row++)
        {
            for(unsigned int i=0;i<row->getsize();i++)
            {
                cols[icol++]=row->getindexptr()[i];
            }
            rows[irow+1]= rows[irow]+row->getsize();
            irow++;
        }

        //export arrays
        FILE *out;
        constexpr size_t bufsize = 256;
        char filename[bufsize];

        std::snprintf(filename,bufsize,"%s/rows.i32",path);
        out=std::fopen(filename,"w");
        std::fwrite(rows.get(), sizeof(int), A.N()+1,out);
        std::fclose(out);

        std::snprintf(filename,bufsize,"%s/cols.i32",path);
        out=std::fopen(filename,"w");
        std::fwrite(cols.get(), sizeof(int), A.nonzeroes(),out);
        std::fclose(out);
    }

} // namespace Opm

#endif // OPM_EXPORT_SYSTEM_HEADER_INCLUDED
