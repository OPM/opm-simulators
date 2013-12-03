/*
    
*/

#ifndef OPM_INCOMPPROPSADINTERFACE_HEADER_INCLUDED
#define OPM_INCOMPPROPSADINTERFACE_HEADER_INCLUDED
#include <opm/autodiff/AutoDiffBlock.hpp>

namespace Opm
{
    class IncompPropsAdInterface
    {
    public:
        virtual ~IncompPropsAdInterface();

        ////////////////////////////
        //      Rock interface    //
        ////////////////////////////

        /// \return   D, the number of spatial dimensions.
        virtual int numDimensions() const = 0;

        /// \return   N, the number of cells.
        virtual int numCells() const = 0;

        /// \return   Array of N porosity values.
        virtual const double* porosity() const = 0;

        /// \return   Array of ND^2 permeability values.
        ///           The D^2 permeability values for a cell are organized as a matrix,
        ///           which is symmetric (so ordering does not matter).
        virtual const double* permeability() const = 0;


        ////////////////////////////
        //      Fluid interface   //
        ////////////////////////////

        typedef AutoDiffBlock<double> ADB;
        typedef ADB::V V;
        typedef ADB::M M;
        typedef std::vector<int> Cells;

        /// \return   Number of active phases (also the number of components).
        virtual int numPhases() const = 0;


        // ------ Canonical named indices for each phase ------

        /// Canonical named indices for each phase.
        enum PhaseIndex { Water = 0, Oil = 1 };

        // ------ Density ------

        /// Densities of stock components at surface conditions.
        /// \return Array of 2 density values.
        virtual const double* surfaceDensity() const = 0;
    

        // ------ Viscosity ------

        /// Viscosity of stock components at surface conditions.
        /// \return Array of 2 viscosity values.
        virtual const double* viscosity() const = 0;



        // ------ Relative permeability ------

        /// Relative permeabilities for all phases.
        /// \param[in]  sw     Array of n water saturation values.
        /// \param[in]  so     Array of n oil saturation values.
        /// \param[in]  cells  Array of n cell indices to be associated with the saturation values.
        /// \return            An std::vector with 2 elements, each an array of n relperm values,
        ///                    containing krw, kro. Use PhaseIndex for indexing into the result.
        virtual
        std::vector<V> relperm(const V& sw,
                               const V& so,
                               const Cells& cells) const = 0;

        /// Relative permeabilities for all phases.
        /// \param[in]  sw     Array of n water saturation values.
        /// \param[in]  so     Array of n oil saturation values.
        /// \param[in]  cells  Array of n cell indices to be associated with the saturation values.
        /// \return            An std::vector with 2 elements, each an array of n relperm values,
        ///                    containing krw, kro. Use PhaseIndex for indexing into the result.
        virtual
        std::vector<ADB> relperm(const ADB& sw,
                                 const ADB& so,
                                 const Cells& cells) const = 0;

    };
} // namespace Opm

#endif// OPM_INCOMPPROPSADINTERFACE_HEADER_INCLUDED
