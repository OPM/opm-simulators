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
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>

#include <opm/core/utility/CompressedPropertyAccess.hpp>

#include <array>
#include <string>
#include <vector>

namespace Opm
{

    // Helper functions
    namespace
    {
        enum PermeabilityKind { ScalarPerm, DiagonalPerm, TensorPerm, None, Invalid };

        void setScalarPermIfNeeded(std::array<int,9>& kmap,
                                   int i, int j, int k);

        typedef GridPropertyAccess::ArrayPolicy::ExtractFromDeck<double> PermArray;

        struct PermTag {};

        typedef GridPropertyAccess::Compressed<PermArray, PermTag> PermComponent;

        PermComponent
        extractPermComponent(EclipseStateConstPtr ecl,
                             const std::string&   kw,
                             const int*           global_cell);

        PermeabilityKind
        fillTensor(EclipseStateConstPtr        eclState,
                   const int*                  global_cell,
                   std::vector<PermComponent>& tensor,
                   std::array<int,9>&          kmap);

    } // anonymous namespace



    // ---- RockFromDeck methods ----


    /// Default constructor.
    RockFromDeck::RockFromDeck()
    {
    }

    void RockFromDeck::init(Opm::EclipseStateConstPtr eclState,
                            int number_of_cells, const int* global_cell,
                            const int* cart_dims)
    {
        assignPorosity(eclState, number_of_cells, global_cell);
        permfield_valid_.assign(number_of_cells, false);
        const double perm_threshold = 0.0; // Maybe turn into parameter?
        assignPermeability(eclState, number_of_cells, global_cell, cart_dims,
                           perm_threshold);
    }

    void RockFromDeck::assignPorosity(Opm::EclipseStateConstPtr eclState,
                                      int number_of_cells, const int* global_cell)
    {
        typedef GridPropertyAccess::ArrayPolicy
            ::ExtractFromDeck<double> Array;

        Array poro_glob(eclState, "PORO", 1.0);
        GridPropertyAccess::Compressed<Array> poro(poro_glob, global_cell);

        porosity_.clear(); porosity_.reserve(number_of_cells);
        for (int c = 0; c < number_of_cells; ++c) {
            porosity_.push_back(poro[c]);
        }
    }

    void RockFromDeck::assignPermeability(Opm::EclipseStateConstPtr eclState,
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

        std::vector<PermComponent> tensor;
        tensor.reserve(6);

        std::array<int,9> kmap;
        PermeabilityKind pkind = fillTensor(eclState, global_cell,
                                            tensor, kmap);
        if (pkind == Invalid) {
            OPM_THROW(std::runtime_error, "Invalid permeability field.");
        }

        assert (! tensor.empty());
        {
            int off = 0;

            for (int c = 0; c < nc; ++c, off += dim*dim) {
                // SharedPermTensor K(dim, dim, &permeability_[off]);
                int kix = 0;

                for (int i = 0; i < dim; ++i) {
                    for (int j = 0; j < dim; ++j, ++kix) {
                        // Clients expect column-major (Fortran) order
                        // in "permeability_" so honour that
                        // requirement despite "tensor" being created
                        // row-major.  Note: The actual numerical
                        // values in the resulting array are the same
                        // in either order when viewed contiguously
                        // because fillTensor() enforces symmetry.
                        permeability_[off + (i + dim*j)] =
                            tensor[kmap[kix]][c];
                    }

                    // K(i,i) = std::max(K(i,i), perm_threshold);
                    double& kii = permeability_[off + i*(dim + 1)];
                    kii = std::max(kii, perm_threshold);
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
        /// @param eclState [in]
        ///    An internalized Eclipse deck from opm-parser which is
        ///    capable of answering which permeability components are
        ///    present in a given input deck.
        ///
        /// @return
        ///    An enum value with the following possible values:
        ///        ScalarPerm     only one component was given.
        ///        DiagonalPerm   more than one component given.
        ///        TensorPerm     at least one cross-component given.
        ///        None           no components given.
        ///        Invalid        invalid set of components given.
        PermeabilityKind classifyPermeability(Opm::EclipseStateConstPtr eclState)
        {
            const bool xx = eclState->hasDoubleGridProperty("PERMX" );
            const bool xy = eclState->hasDoubleGridProperty("PERMXY");
            const bool yx = xy;

            const bool yy = eclState->hasDoubleGridProperty("PERMY" );
            const bool yz = eclState->hasDoubleGridProperty("PERMYZ");
            const bool zy = yz;

            const bool zz = eclState->hasDoubleGridProperty("PERMZ" );
            const bool zx = eclState->hasDoubleGridProperty("PERMZX");
            const bool xz = zx;

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
            if (kmap[j] < 0) { kmap[j] = kmap[i]; }
            if (kmap[k] < 0) { kmap[k] = kmap[i]; }
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
        /// @param [in]  eclState
        ///    An internalized Eclipse deck object which capable of
        ///    answering which permeability components are present in
        ///    a given input deck as well as retrieving the numerical
        ///    value of each permeability component in each grid cell.
        ///
        /// @param [out] tensor
        /// @param [out] kmap
        PermeabilityKind
        fillTensor(EclipseStateConstPtr        eclState,
                   const int*                  global_cell,
                   std::vector<PermComponent>& tensor,
                   std::array<int,9>&          kmap)
        {
            PermeabilityKind kind = classifyPermeability(eclState);
            if (kind == Invalid) {
                OPM_THROW(std::runtime_error, "Invalid set of permeability fields given.");
            }

            assert (tensor.empty());

            for (int i = 0; i < 9; ++i) { kmap[i] = -1; }

            enum { xx, xy, xz,    // 0, 1, 2
                   yx, yy, yz,    // 3, 4, 5
                   zx, zy, zz };  // 6, 7, 8

            // -----------------------------------------------------------
            // 1st row: [ kxx, kxy ], kxz handled in kzx
            if (eclState->hasDoubleGridProperty("PERMX" )) {
                kmap[xx] = tensor.size();
                tensor.push_back(extractPermComponent(eclState, "PERMX", global_cell));

                setScalarPermIfNeeded(kmap, xx, yy, zz);
            }
            {
                kmap[xy] = kmap[yx] = tensor.size();  // Enforce symmetry.
                tensor.push_back(extractPermComponent(eclState, "PERMXY", global_cell));
            }

            // -----------------------------------------------------------
            // 2nd row: [ kyy, kyz ], kyx handled in kxy
            if (eclState->hasDoubleGridProperty("PERMY" )) {
                kmap[yy] = tensor.size();
                tensor.push_back(extractPermComponent(eclState, "PERMY", global_cell));

                setScalarPermIfNeeded(kmap, yy, zz, xx);
            }
            {
                kmap[yz] = kmap[zy] = tensor.size();  // Enforce symmetry.
                tensor.push_back(extractPermComponent(eclState, "PERMYZ", global_cell));
            }

            // -----------------------------------------------------------
            // 3rd row: [ kzx, kzz ], kzy handled in kyz
            {
                kmap[zx] = kmap[xz] = tensor.size();  // Enforce symmetry.
                tensor.push_back(extractPermComponent(eclState, "PERMZX", global_cell));
            }
            if (eclState->hasDoubleGridProperty("PERMZ" )) {
                kmap[zz] = tensor.size();
                tensor.push_back(extractPermComponent(eclState, "PERMZ", global_cell));

                setScalarPermIfNeeded(kmap, zz, xx, yy);
            }

            return kind;
        }

        PermComponent
        extractPermComponent(EclipseStateConstPtr ecl,
                             const std::string&   kw,
                             const int*           global_cell)
        {
            PermArray k(ecl, kw, 0.0); // return 0.0 if not present.

            return PermComponent(k, global_cell);
        }
    } // anonymous namespace

} // namespace Opm
