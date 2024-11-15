

#ifndef OPM_COMP_WELL_HPP
#define OPM_COMP_WELL_HPP

#include <string>

namespace Opm
{

template <typename Scalar>
class CompWell
{
public:
    CompWell(const std::string& name)
        : name(name)
    {
    }
private:
    std::string name;
};

} // end of namespace Opm
#endif // OPM_COMPOSITIONAL_WELL_MODEL_HPP