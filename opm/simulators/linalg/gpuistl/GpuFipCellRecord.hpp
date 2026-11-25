// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2026 Equinor ASA

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
*/
/*!
 * \file
 *
 * \brief Typed owner for the GPU property-evaluation to TPFA-assembly handoff.
 */
// Compact gas-water values needed by CPU FIP and regional-density calculations.
#ifndef OPM_GPU_FIP_CELL_RECORD_HPP
#define OPM_GPU_FIP_CELL_RECORD_HPP

namespace Opm::gpuistl {

template <class Scalar>
struct GpuFipCellRecord {
    Scalar porosityScalar{};
    Scalar referencePorosityScalar{};
    Scalar saltConcentration{};
    Scalar rsw{}, rvw{};
    Scalar gasRelativePermeability{};
    // The current bridge has water (canonical phase 0) and gas (phase 2).
    // Slot 0 is water and slot 1 is gas; the inactive oil phase is zero.
    Scalar saturation[2]{};
    Scalar invB[2]{};
    Scalar pressure[2]{};
    Scalar density[2]{};
    unsigned pvtRegionIndex{};

    static Scalar phaseValue(const Scalar (&values)[2], unsigned phase)
    {
        return phase == 1 ? Scalar{} : values[phase == 2 ? 1 : 0];
    }

    struct Value {
        Scalar scalar;
        Scalar value() const { return scalar; }
    };

    struct FluidState {
        const GpuFipCellRecord& cell;
        Scalar saturation(unsigned phase) const { return GpuFipCellRecord::phaseValue(cell.saturation, phase); }
        Scalar invB(unsigned phase) const { return GpuFipCellRecord::phaseValue(cell.invB, phase); }
        Scalar pressure(unsigned phase) const { return GpuFipCellRecord::phaseValue(cell.pressure, phase); }
        Scalar density(unsigned phase) const { return GpuFipCellRecord::phaseValue(cell.density, phase); }
        Scalar Rs() const { return Scalar{}; }
        Scalar Rv() const { return Scalar{}; }
        Scalar Rsw() const { return cell.rsw; }
        Scalar Rvw() const { return cell.rvw; }
        Value saltConcentration() const { return {cell.saltConcentration}; }
        unsigned pvtRegionIndex() const { return cell.pvtRegionIndex; }
    };

    FluidState fluidState() const { return {*this}; }
    Value porosity() const { return {porosityScalar}; }
    Scalar referencePorosity() const { return referencePorosityScalar; }
    Scalar relativePermeability(unsigned) const { return gasRelativePermeability; }
};

} // namespace Opm::gpuistl

#endif
