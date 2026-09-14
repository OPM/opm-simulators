#include <config.h>
#include <opm/simulators/wells/PerforationData.hpp>

namespace Opm {

template<class Scalar>
PerforationData<Scalar> PerforationData<Scalar>::serializationTestObject()
{
    PerforationData result;
    result.cell_index = 17;
    result.connection_transmissibility_factor = 1.25e-8;
    result.connection_d_factor = 3.5e-4;
    result.satnum_id = 23;
    result.ecl_index = 29;
    result.grid_id = 31;
    result.global_index = 37;
    return result;
}

template<class Scalar>
bool PerforationData<Scalar>::operator==(const PerforationData& rhs) const
{
    return this->cell_index == rhs.cell_index
        && this->connection_transmissibility_factor == rhs.connection_transmissibility_factor
        && this->connection_d_factor == rhs.connection_d_factor
        && this->satnum_id == rhs.satnum_id
        && this->ecl_index == rhs.ecl_index
        && this->grid_id == rhs.grid_id
        && this->global_index == rhs.global_index;
}

template struct PerforationData<double>;

#if FLOW_INSTANTIATE_FLOAT
template struct PerforationData<float>;
#endif

} // namespace Opm
