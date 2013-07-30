/* Copyright (c) 2013 Uni Research AS.
   This file is licensed under the GNU General Public License v3.0 or later. */
#ifndef OPM_INCOMPPROPERTIESSHADOW_HEADER_INCLUDED
#define OPM_INCOMPPROPERTIESSHADOW_HEADER_INCLUDED

#ifndef OPM_INCOMPPROPERTIESINTERFACE_HEADER_INCLUDED
#include <opm/core/props/IncompPropertiesInterface.hpp>
#endif /* OPM_INCOMPPROPERTIESINTERFACE_HEADER_INCLUDED */

namespace Opm
{
    /**
     * Override certain properties with values from elsewhere.
     *
     * This allows mixing of property objects from several sources,
     * such as rock and fluid properties from a file but unsaturated
     * properties from a function. Care must be taken to setup the
     * shadowing so no inconsistencies arise.
     *
     * @remark
     *     This object is mutable; if you change some properties
     *     it will affect all clients that have references to it.
     *     It is thus recommended to only use the mutable portion
     *     when constructing the object, before passing it to clients.
     *
     * @example
     * @code{.cpp}
     *   std::vector<double> poro;
     *   IncompPropertiesFromDeck fromDeck(deck, grid);
     *   simulate (IncompPropertiesShadow(fromDeck).usePorosity(poro));
     * @endcode
     */
    struct IncompPropertiesShadow : public IncompPropertiesInterface
    {
        /**
         * Shadow another set of properties. If no properties are
         * overridden, the values from the original will be used.
         */
        IncompPropertiesShadow (const IncompPropertiesInterface& original);

        /**
         * Implement all methods from the IncompPropertiesInterface.
         */
        virtual int numDimensions () const;
        virtual int numCells () const;
        virtual const double* porosity () const;
        virtual const double* permeability () const;
        virtual int numPhases () const;
        virtual const double* viscosity () const;
        virtual const double* density () const;
        virtual const double* surfaceDensity () const;
        virtual void relperm (const int n,
                              const double* s,
                              const int* cells,
                              double* kr,
                              double* dkrds) const;
        virtual void capPress (const int n,
                               const double* s,
                               const int* cells,
                               double* pc,
                               double* dpcds) const;
        virtual void satRange (const int n,
                               const int* cells,
                               double* smin,
                               double* smax) const;

        /**
         * Use a different set of porosities.
         *
         * @param poro
         *     Iterator containing new porosity values. It must contain
         *     numCells() values.
         * @return
         *     A reference to this object, so it can be used for chaining.
         * @remark
         *     This object does *not* assume ownership of the underlaying
         *     memory nor makes any copies of it. Hence, the calling code
         *     must manage the array so that it points to valid memory for
         *     the lifetime of this object.
         */
        IncompPropertiesShadow& usePorosity (const double* poro);
        IncompPropertiesShadow& usePorosity (const IncompPropertiesInterface& other);

        /**
         * Use a different set of permeabilities.
         *
         * @param perm
         *     Iterator containing new permeability values. It must contain
         *     numCells()*numDimensions()*numDimensions() values.
         * @return
         *     A reference to this object, so it can be used for chaining.
         * @remark
         *     This object does *not* assume ownership of the underlaying
         *     memory nor makes any copies of it. Hence, the calling code
         *     must manage the array so that it points to valid memory for
         *     the lifetime of this object.
         */
        IncompPropertiesShadow& usePermeability (const double* perm);
        IncompPropertiesShadow& usePermeability (const IncompPropertiesInterface& other);

        /**
         * Use a different set of viscosities.
         *
         * @param visc
         *     Iterator containing new viscosity values. It must contain
         *     numPhases() values.
         * @return
         *     A reference to this object, so it can be used for chaining.
         * @remark
         *     This object does *not* assume ownership of the underlaying
         *     memory nor makes any copies of it. Hence, the calling code
         *     must manage the array so that it points to valid memory for
         *     the lifetime of this object.
         */
        IncompPropertiesShadow& useViscosity (const double* visc);
        IncompPropertiesShadow& useViscosity (const IncompPropertiesInterface& other);

        /**
         * Use a different set of densities.
         *
         * @param dens
         *     Iterator containing new density values. It must contain
         *     numPhases() values.
         * @return
         *     A reference to this object, so it can be used for chaining.
         * @remark
         *     This object does *not* assume ownership of the underlaying
         *     memory nor makes any copies of it. Hence, the calling code
         *     must manage the array so that it points to valid memory for
         *     the lifetime of this object.
         */
        IncompPropertiesShadow& useDensity (const double* dens);
        IncompPropertiesShadow& useDensity (const IncompPropertiesInterface& other);

        /**
         * Use a different set of surface densities.
         *
         * @param surf
         *     Iterator containing new surface density values. It must
         *     contain numPhases() values.
         * @return
         *     A reference to this object, so it can be used for chaining.
         * @remark
         *     This object does *not* assume ownership of the underlaying
         *     memory nor makes any copies of it. Hence, the calling code
         *     must manage the array so that it points to valid memory for
         *     the lifetime of this object.
         */
        IncompPropertiesShadow& useSurfaceDensity (const double* surf);
        IncompPropertiesShadow& useSurfaceDensity (const IncompPropertiesInterface& other);

        /**
         * Convenience method to set both porosity and permeability.
         */
        IncompPropertiesShadow& useRockProps (const IncompPropertiesInterface& other);

        /**
         * Convenience method to set both viscosity and density.
         */
        IncompPropertiesShadow& useFluidProps (const IncompPropertiesInterface& other);

        /**
         * Convenience method to set both rock and fluid properties.
         */
        IncompPropertiesShadow& useRockAndFluidProps (const IncompPropertiesInterface& other);

    private:
        /**
         * If we haven't set a property explicitly, then retrieve
         * them from this. This is a kind of prototype inheritance,
         * hence the name of this field.
         */
        const IncompPropertiesInterface& prototype;

        /**
         * Bitfield which tells us which properties that has been
         * shadowed. The others are retrieved from the original
         * interface.
         */
        int shadowed;

        /**
         * Bits that indicates which fields that has been overridden.
         */
        static const int POROSITY        = 1 << 1;
        static const int PERMEABILITY    = 1 << 2;
        static const int VISCOSITY       = 1 << 3;
        static const int DENSITY         = 1 << 4;
        static const int SURFACE_DENSITY = 1 << 5;

        /**
         * Pointers to alternative values. These pointers should only
         * be assumed to be valid if the corresponding bit in the mask
         * is set. No management is done for the memory this points to!
         */
        const double* poro_;
        const double* perm_;
        const double* visc_;
        const double* dens_;
        const double* surf_;
    };
} /* namespace Opm */

// body of inline methods are defined here:
#include <opm/core/props/IncompPropertiesShadow_impl.hpp>

#endif /* OPM_INCOMPPROPERTIESSHADOW_HEADER_INCLUDED */
