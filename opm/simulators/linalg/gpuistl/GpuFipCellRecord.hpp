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
