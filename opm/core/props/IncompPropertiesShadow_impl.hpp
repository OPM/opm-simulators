/* Copyright (c) 2013 Uni Research AS.
   This file is licensed under the GNU General Public License v3.0 or later. */
#ifndef OPM_INCOMPPROPERTIESSHADOW_HEADER_INCLUDED
#error Do not include IncompPropertiesShadow_impl.hpp directly!
#endif /* OPM_INCOMPPROPERTIESSHADOW_HEADER_INCLUDED */

namespace Opm
{
    /**
     * Initialize so that all properties are retrieved from original.
     */
    inline IncompPropertiesShadow::IncompPropertiesShadow (const IncompPropertiesInterface& original)
        : prototype_ (original)
        , shadowed_ (0)
        , poro_ (0)
        , perm_ (0)
        , visc_ (0)
        , dens_ (0)
        , surf_ (0)
    {
    }

    /**
     * The format of the prototype and the shadow must be the same,
     * so these methods should always be forwarded directly.
     */
    inline int IncompPropertiesShadow::numDimensions () const
    {
        return prototype_.numDimensions();
    }

    inline int IncompPropertiesShadow::numCells () const
    {
        return prototype_.numCells();
    }

    inline int IncompPropertiesShadow::numPhases () const
    {
        return prototype_.numPhases();
    }

    /**
     * These methods are sufficiently advanced (the s parameter is a
     * non-integral index) for there not to be a trivial implementation,
     * so they are not overridden yet.
     */
    inline void IncompPropertiesShadow::relperm (const int n,
                                                 const double* s,
                                                 const int* cells,
                                                 double* kr,
                                                 double* dkrds) const
    {
        prototype_.relperm (n, s, cells, kr, dkrds);
    }

    inline void IncompPropertiesShadow::capPress (const int n,
                                                  const double* s,
                                                  const int* cells,
                                                  double* pc,
                                                  double* dpcds) const
    {
        prototype_.capPress (n, s, cells, pc, dpcds);
    }

    inline void IncompPropertiesShadow::satRange (const int n,
                                                  const int* cells,
                                                  double* smin,
                                                  double* smax) const
    {
        prototype_.satRange (n, cells, smin, smax);
    }

    /**
     * Return the new value if indicated in the bitfield, otherwise
     * use the original value from the other object.
     */
    inline const double* IncompPropertiesShadow::porosity () const
    {
        return (shadowed_ & POROSITY) ? poro_ : prototype_.porosity ();
    }

    inline const double* IncompPropertiesShadow::permeability () const
    {
        return (shadowed_ & PERMEABILITY) ? perm_ : prototype_.permeability ();
    }

    inline const double* IncompPropertiesShadow::viscosity () const
    {
        return (shadowed_ & VISCOSITY) ? visc_ : prototype_.viscosity ();
    }

    inline const double* IncompPropertiesShadow::density () const
    {
        return (shadowed_ & DENSITY) ? dens_ : prototype_.density ();
    }

    inline const double* IncompPropertiesShadow::surfaceDensity () const
    {
        return (shadowed_ & SURFACE_DENSITY) ? surf_ : prototype_.surfaceDensity ();
    }

    /**
     * Store the pointer and indicate that the new value should be used.
     */
    inline IncompPropertiesShadow& IncompPropertiesShadow::usePorosity (const double* poro)
    {
        this->poro_ = poro;
        shadowed_ |= POROSITY;
        return *this;
    }

    inline IncompPropertiesShadow& IncompPropertiesShadow::usePermeability (const double* perm)
    {
        this->perm_ = perm;
        shadowed_ |= PERMEABILITY;
        return *this;
    }

    inline IncompPropertiesShadow& IncompPropertiesShadow::useViscosity (const double* visc)
    {
        this->visc_ = visc;
        shadowed_ |= VISCOSITY;
        return *this;
    }

    inline IncompPropertiesShadow& IncompPropertiesShadow::useDensity (const double* dens)
    {
        this->dens_ = dens;
        shadowed_ |= DENSITY;
        return *this;
    }

    inline IncompPropertiesShadow& IncompPropertiesShadow::useSurfaceDensity (const double* surf)
    {
        this->surf_ = surf;
        shadowed_ |= SURFACE_DENSITY;
        return *this;
    }

    /**
     * Copy the pointer from another property interface, after checking
     * that they are compatible.
     */
    inline IncompPropertiesShadow& IncompPropertiesShadow::usePorosity (const IncompPropertiesInterface& other)
    {
        assert (prototype_.numCells() == other.numCells());
        return usePorosity (other.porosity());
    }

    inline IncompPropertiesShadow& IncompPropertiesShadow::usePermeability (const IncompPropertiesInterface& other)
    {
        assert (prototype_.numCells() == other.numCells());
        assert (prototype_.numDimensions() == other.numDimensions());
        return usePermeability (other.permeability());
    }

    inline IncompPropertiesShadow& IncompPropertiesShadow::useViscosity (const IncompPropertiesInterface& other)
    {
        assert (prototype_.numPhases() == other.numPhases());
        return useViscosity (other.viscosity());
    }

    inline IncompPropertiesShadow& IncompPropertiesShadow::useDensity (const IncompPropertiesInterface& other)
    {
        assert (prototype_.numPhases() == other.numPhases());
        return useDensity (other.density());
    }

    inline IncompPropertiesShadow& IncompPropertiesShadow::useSurfaceDensity (const IncompPropertiesInterface& other)
    {
        assert (prototype_.numPhases() == other.numPhases());
        return useSurfaceDensity (other.surfaceDensity());
    }

    /**
     * Convenience methods to set several set of properties at once.
     */
    inline IncompPropertiesShadow& IncompPropertiesShadow::useRockProps (const IncompPropertiesInterface& other)
    {
        return usePorosity (other).usePermeability (other);
    }

    inline IncompPropertiesShadow& IncompPropertiesShadow::useFluidProps (const IncompPropertiesInterface& other)
    {
        return useViscosity (other).useDensity (other).useSurfaceDensity (other);
    }

    inline IncompPropertiesShadow& IncompPropertiesShadow::useRockAndFluidProps (const IncompPropertiesInterface& other)
    {
        return useRockProps (other).useFluidProps (other);
    }
} /* namespace Opm */
