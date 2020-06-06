/*
  Copyright 2019 SINTEF Digital, Mathematics and Cybernetics.

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

#ifndef OPM_PRESSURMATRIXHELPER_HEADER_INCLUDED
#define OPM_PRESSURMATRIXHELPER_HEADER_INCLUDED

namespace Opm{
    namespace PressureHelper{
        // this make pressure equations or jacobi of pressure equation
        // it is identical to the cpr PressureTransferPolicy.hpp  placed not in a linear solver setting
        template<class MatrixType,class PressureMatrixType,class PressureVectorType>
        PressureMatrixType makePressureMatrix(const MatrixType& fineLevelMatrix,
                                              int pressureVarIndex,
                                              const PressureVectorType& weights)
        {
            PressureMatrixType pressureMatrix(fineLevelMatrix.N(), fineLevelMatrix.M(), PressureMatrixType::row_wise);
            auto createIter = pressureMatrix.createbegin();
            
            for (const auto& row : fineLevelMatrix) {
                for (auto col = row.begin(), cend = row.end(); col != cend; ++col) {
                    createIter.insert(col.index());
                }
                ++createIter;
            }
            makePressureMatrixEntries(pressureMatrix,fineLevelMatrix,pressureVarIndex, weights);
            return pressureMatrix;
        }
        template<class MatrixType,class PressureMatrixType,class BlockVectorType>
        void makePressureMatrixEntries(PressureMatrixType& pmatrix,
                                       const MatrixType& fineMatrix,
                                       const int pressure_var_index,
                                       const BlockVectorType& weights){
            pmatrix = 0;
            auto rowCoarse = pmatrix.begin();
            for (auto row = fineMatrix.begin(), rowEnd = fineMatrix.end(); row != rowEnd; ++row, ++rowCoarse) {
                assert(row.index() == rowCoarse.index());
                auto entryCoarse = rowCoarse->begin();
                for (auto entry = row->begin(), entryEnd = row->end(); entry != entryEnd; ++entry, ++entryCoarse) {
                    assert(entry.index() == entryCoarse.index());
                    double matrix_el = 0;
                    const auto& bw = weights[row.index()];
                    for (size_t i = 0; i < bw.size(); ++i) {
                        matrix_el += (*entry)[i][pressure_var_index] * bw[i];
                    }
                    (*entryCoarse) = matrix_el;
                }
            }
            assert(rowCoarse == pmatrix.end());
        }
        
        template<class BlockVectorType, class PressureVectorType>
        void moveToPressureEqn(const BlockVectorType& fine,PressureVectorType& rhs,const BlockVectorType& weights) 
        {
            // Set coarse vector to zero
            rhs = 0;    
            auto end = fine.end(), begin = fine.begin();
            
            for (auto block = begin; block != end; ++block) {
                const auto& bw = weights[block.index()];
                double rhs_el = 0.0;
                for (size_t i = 0; i < block->size(); ++i) {
                    rhs_el += (*block)[i] * bw[i];
                }
                rhs[block - begin] = rhs_el;
            }
        }
        template<class BlockVectorType, class PressureVectorType>
        void movePressureToBlock(BlockVectorType& fine,const PressureVectorType& lhs,int pressure_var_index)
        {
            auto end = fine.end(), begin = fine.begin();
            for (auto block = begin; block != end; ++block) {
                (*block)[pressure_var_index] = lhs[block - begin];
            }
        }
    }
}
        
#endif
