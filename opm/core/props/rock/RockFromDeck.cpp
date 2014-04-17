/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

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


#include "config.h"
#include <opm/core/props/rock/RockFromDeck.hpp>
#include <opm/core/grid.h>
#include <opm/core/utility/ErrorMacros.hpp>

#include <opm/parser/eclipse/Deck/Deck.hpp>

#include <array>

namespace Opm
{

    // Helper functions
    namespace
    {
        enum PermeabilityKind { ScalarPerm, DiagonalPerm, TensorPerm, None, Invalid };

        void setScalarPermIfNeeded(std::array<int,9>& kmap,
                                   int i, int j, int k);
        PermeabilityKind fillTensor(Opm::DeckConstPtr newParserDeck,
                                    std::vector<const std::vector<double>*>& tensor,
                                    std::array<int,9>&                     kmap);

    } // anonymous namespace



    // ---- RockFromDeck methods ----


    /// Default constructor.
    RockFromDeck::RockFromDeck()
    {
    }

    void RockFromDeck::init(Opm::DeckConstPtr newParserDeck,
                            int number_of_cells, const int* global_cell,
                            const int* cart_dims)
    {
        assignPorosity(newParserDeck, number_of_cells, global_cell);
        permfield_valid_.assign(number_of_cells, false);
        const double perm_threshold = 0.0; // Maybe turn into parameter?
        assignPermeability(newParserDeck, number_of_cells, global_cell, cart_dims,
                           perm_threshold);
    }

    void RockFromDeck::assignPorosity(Opm::DeckConstPtr newParserDeck,
                                      int number_of_cells, const int* global_cell)
    {
        porosity_.assign(number_of_cells, 1.0);
        if (newParserDeck->hasKeyword("PORO")) {
            const std::vector<double>& poro = newParserDeck->getKeyword("PORO")->getSIDoubleData();
            for (int c = 0; c < int(porosity_.size()); ++c) {
                const int deck_pos = (global_cell == NULL) ? c : global_cell[c];
                assert(0 <= c && c < (int) porosity_.size());
                assert(0 <= deck_pos && deck_pos < (int) poro.size());
                porosity_[c] = poro[deck_pos];
            }
        }
    }

    void RockFromDeck::assignPermeability(Opm::DeckConstPtr newParserDeck,
                                          int number_of_cells,
                                          const int* global_cell,
                                          const int* cartdims,
                                          double perm_threshold)
    {
        const int dim              = 3;
        const int num_global_cells = cartdims[0]*cartdims[1]*cartdims[2];
        const int nc = number_of_cells;

        assert(num_global_cells > 0);

        permeability_.assign(dim * dim * nc, 0.0);

        std::vector<const std::vector<double>*> tensor;
        tensor.reserve(10);

        const std::vector<double> zero(num_global_cells, 0.0);
        tensor.push_back(&zero);

        std::array<int,9> kmap;
        PermeabilityKind pkind = fillTensor(newParserDeck, tensor, kmap);
        if (pkind == Invalid) {
            OPM_THROW(std::runtime_error, "Invalid permeability field.");
        }

        // Assign permeability values only if such values are
        // given in the input deck represented by 'newParserDeck'.  In
        // other words: Don't set any (arbitrary) default values.
        // It is infinitely better to experience a reproducible
        // crash than subtle errors resulting from a (poorly
        // chosen) default value...
        //
        if (tensor.size() > 1) {
            const int* gc = global_cell;
            int off = 0;

            for (int c = 0; c < nc; ++c, off += dim*dim) {
                // SharedPermTensor K(dim, dim, &permeability_[off]);
                int       kix  = 0;
                const int glob = (gc == NULL) ? c : gc[c];

                for (int i = 0; i < dim; ++i) {
                    for (int j = 0; j < dim; ++j, ++kix) {
                        // K(i,j) = (*tensor[kmap[kix]])[glob];
                        permeability_[off + kix] = (*tensor[kmap[kix]])[glob];
                    }
                    // K(i,i) = std::max(K(i,i), perm_threshold);
                    permeability_[off + 3*i + i] = std::max(permeability_[off + 3*i + i], perm_threshold);
                }

                permfield_valid_[c] = std::vector<unsigned char>::value_type(1);
            }
        }
    }

    namespace {
        /// @brief
        ///    Classify and verify a given permeability specification
        ///    from a structural point of view.  In particular, we
        ///    verify that there are no off-diagonal permeability
        ///    components such as @f$k_{xy}@f$ unless the
        ///    corresponding diagonal components are known as well.
        ///
        /// @param newParserDeck [in]
        ///    An Eclipse data parser capable of answering which
        ///    permeability components are present in a given input
        ///    deck.
        ///
        /// @return
        ///    An enum value with the following possible values:
        ///        ScalarPerm     only one component was given.
        ///        DiagonalPerm   more than one component given.
        ///        TensorPerm     at least one cross-component given.
        ///        None           no components given.
        ///        Invalid        invalid set of components given.
        PermeabilityKind classifyPermeability(Opm::DeckConstPtr newParserDeck)
        {
            const bool xx = newParserDeck->hasKeyword("PERMX" );
            const bool xy = newParserDeck->hasKeyword("PERMXY");
            const bool xz = newParserDeck->hasKeyword("PERMXZ");

            const bool yx = newParserDeck->hasKeyword("PERMYX");
            const bool yy = newParserDeck->hasKeyword("PERMY" );
            const bool yz = newParserDeck->hasKeyword("PERMYZ");

            const bool zx = newParserDeck->hasKeyword("PERMZX");
            const bool zy = newParserDeck->hasKeyword("PERMZY");
            const bool zz = newParserDeck->hasKeyword("PERMZ" );

            int num_cross_comp = xy + xz + yx + yz + zx + zy;
            int num_comp       = xx + yy + zz + num_cross_comp;
            PermeabilityKind retval = None;
            if (num_cross_comp > 0) {
                retval = TensorPerm;
            } else {
                if (num_comp == 1) {
                    retval = ScalarPerm;
                } else if (num_comp >= 2) {
                    retval = DiagonalPerm;
                }
            }

            bool ok = true;
            if (num_comp > 0) {
                // At least one tensor component specified on input.
                // Verify that any remaining components are OK from a
                // structural point of view.  In particular, there
                // must not be any cross-components (e.g., k_{xy})
                // unless the corresponding diagonal component (e.g.,
                // k_{xx}) is present as well...
                //
                ok =        xx || !(xy || xz || yx || zx) ;
                ok = ok && (yy || !(yx || yz || xy || zy));
                ok = ok && (zz || !(zx || zy || xz || yz));
            }
            if (!ok) {
                retval = Invalid;
            }

            return retval;
        }


        /// @brief
        ///    Copy isotropic (scalar) permeability to other diagonal
        ///    components if the latter have not (yet) been assigned a
        ///    separate value.  Specifically, this function assigns
        ///    copies of the @f$i@f$ permeability component (e.g.,
        ///    'PERMX') to the @f$j@f$ and @f$k@f$ permeability (e.g.,
        ///    'PERMY' and 'PERMZ') components if these have not
        ///    previously been assigned.
        ///
        /// @param kmap
        ///    Permeability indirection map.  In particular @code
        ///    kmap[i] @endcode is the index (an integral number in
        ///    the set [1..9]) into the permeability tensor
        ///    representation of function @code fillTensor @endcode
        ///    which represents permeability component @code i
        ///    @endcode.
        ///
        /// @param [in] i
        /// @param [in] j
        /// @param [in] k
        void setScalarPermIfNeeded(std::array<int,9>& kmap,
                                   int i, int j, int k)
        {
            if (kmap[j] == 0) { kmap[j] = kmap[i]; }
            if (kmap[k] == 0) { kmap[k] = kmap[i]; }
        }

        /// @brief
        ///   Extract pointers to appropriate tensor components from
        ///   input deck.  The permeability tensor is, generally,
        ///   @code
        ///        [ kxx  kxy  kxz ]
        ///    K = [ kyx  kyy  kyz ]
        ///        [ kzx  kzy  kzz ]
        ///   @endcode
        ///   We store these values in a linear array using natural
        ///   ordering with the column index cycling the most rapidly.
        ///   In particular we use the representation
        ///   @code
        ///        [  0    1    2    3    4    5    6    7    8  ]
        ///    K = [ kxx, kxy, kxz, kyx, kyy, kyz, kzx, kzy, kzz ]
        ///   @endcode
        ///   Moreover, we explicitly enforce symmetric tensors by
        ///   assigning
        ///   @code
        ///     3     1       6     2       7     5
        ///    kyx = kxy,    kzx = kxz,    kzy = kyz
        ///   @endcode
        ///   However, we make no attempt at enforcing positive
        ///   definite tensors.
        ///
        /// @param [in]  parser
        ///    An Eclipse data parser capable of answering which
        ///    permeability components are present in a given input
        ///    deck as well as retrieving the numerical value of each
        ///    permeability component in each grid cell.
        ///
        /// @param [out] tensor
        /// @param [out] kmap
        PermeabilityKind fillTensor(Opm::DeckConstPtr newParserDeck,
                                    std::vector<const std::vector<double>*>& tensor,
                                    std::array<int,9>&                     kmap)
        {
            PermeabilityKind kind = classifyPermeability(newParserDeck);
            if (kind == Invalid) {
                OPM_THROW(std::runtime_error, "Invalid set of permeability fields given.");
            }
            assert(tensor.size() == 1);
            for (int i = 0; i < 9; ++i) { kmap[i] = 0; }

            enum { xx, xy, xz,    // 0, 1, 2
                   yx, yy, yz,    // 3, 4, 5
                   zx, zy, zz };  // 6, 7, 8

            // -----------------------------------------------------------
            // 1st row: [kxx, kxy, kxz]
            if (newParserDeck->hasKeyword("PERMX" )) {
                kmap[xx] = tensor.size();
                tensor.push_back(&newParserDeck->getKeyword("PERMX")->getSIDoubleData());

                setScalarPermIfNeeded(kmap, xx, yy, zz);
            }
            if (newParserDeck->hasKeyword("PERMXY")) {
                kmap[xy] = kmap[yx] = tensor.size();  // Enforce symmetry.
                tensor.push_back(&newParserDeck->getKeyword("PERMXY")->getSIDoubleData());
            }
            if (newParserDeck->hasKeyword("PERMXZ")) {
                kmap[xz] = kmap[zx] = tensor.size();  // Enforce symmetry.
                tensor.push_back(&newParserDeck->getKeyword("PERMXZ")->getSIDoubleData());
            }

            // -----------------------------------------------------------
            // 2nd row: [kyx, kyy, kyz]
            if (newParserDeck->hasKeyword("PERMYX")) {
                kmap[yx] = kmap[xy] = tensor.size();  // Enforce symmetry.
                tensor.push_back(&newParserDeck->getKeyword("PERMYX")->getSIDoubleData());
            }
            if (newParserDeck->hasKeyword("PERMY" )) {
                kmap[yy] = tensor.size();
                tensor.push_back(&newParserDeck->getKeyword("PERMY")->getSIDoubleData());

                setScalarPermIfNeeded(kmap, yy, zz, xx);
            }
            if (newParserDeck->hasKeyword("PERMYZ")) {
                kmap[yz] = kmap[zy] = tensor.size();  // Enforce symmetry.
                tensor.push_back(&newParserDeck->getKeyword("PERMYZ")->getSIDoubleData());
            }

            // -----------------------------------------------------------
            // 3rd row: [kzx, kzy, kzz]
            if (newParserDeck->hasKeyword("PERMZX")) {
                kmap[zx] = kmap[xz] = tensor.size();  // Enforce symmetry.
                tensor.push_back(&newParserDeck->getKeyword("PERMZX")->getSIDoubleData());
            }
            if (newParserDeck->hasKeyword("PERMZY")) {
                kmap[zy] = kmap[yz] = tensor.size();  // Enforce symmetry.
                tensor.push_back(&newParserDeck->getKeyword("PERMZY")->getSIDoubleData());
            }
            if (newParserDeck->hasKeyword("PERMZ" )) {
                kmap[zz] = tensor.size();
                tensor.push_back(&newParserDeck->getKeyword("PERMZ")->getSIDoubleData());

                setScalarPermIfNeeded(kmap, zz, xx, yy);
            }
            return kind;
        }

    } // anonymous namespace

} // namespace Opm
