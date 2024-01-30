/*
  Copyright 2023 Inria, Bretagne–Atlantique Research Center

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
#include <config.h>
#include <opm/simulators/utils/DamarisVar.hpp>

#include <opm/common/ErrorMacros.hpp>

#include <Damaris.h>
#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <string_view>
#include <typeinfo>

namespace Opm::DamarisOutput {

DamarisVarXMLAttributes::DamarisVarXMLAttributes()
{
    // Additional data needed to complete an XML <variable> element
    type_ = "scalar"; // This is probably not needed as vector data is defined using
                      // the Layout paramter. Could be useful for cross checking
    visualizable_ = "false";
    time_varying_ = "true";
    centering_ = "zonal";
}

std::string DamarisVarXMLAttributes::ReturnXMLForVariable()
{
    std::string var_str;

    using Entry = std::pair<std::string_view, const std::string&>;
    auto addAttrib = [&var_str](const Entry& entry)
    {
        if (!entry.second.empty()) {
            var_str += fmt::format(" {}=\"{}\"", entry.first, entry.second);
        }
    };

    const auto entries = std::array{
        Entry{"layout", this->layout_},
        Entry{"mesh", this->mesh_},
        Entry{"type", this->type_},
        Entry{"visualizable", this->visualizable_},
        Entry{"unit", this->unit_},
        Entry{"time_varying", this->time_varying_},
        Entry{"centering", this->centering_},
        Entry{"store", this->store_},
        Entry{"select-mem", this->select_mem_},
        Entry{"select-file", this->select_file_},
        Entry{"select-subset", this->select_subset_}
    };

    std::for_each(entries.begin(), entries.end(), addAttrib);

    return var_str;
}

template<class T>
DamarisVar<T>::DamarisVar(int dims,
                          const std::vector<std::string>& param_names,
                          const std::string& variable_name, int rank)
    : param_names_(param_names)
    , variable_name_(variable_name)
    , rank_(rank)
{
    dam_err_ = DAMARIS_OK;

    assert(param_names_.size() == static_cast<std::size_t>(dims));
    assert(dims > 0);

    has_error_ = false;

    // Check that our template type T matches out Damaris XML <layout> type
    TestType(variable_name);
    if (hasError()) {
        printError(); // throws a runtime error, with error message from
                  // dam_err_sstr_
    }

    current_size_ = 0;
    num_params_ = param_names_.size();
    param_sizes_.resize(num_params_);
    positions_.resize(dims);

    data_ptr_ = nullptr;
    parameters_set_ = false;
    has_error_ = false;
}

template<class T>
DamarisVar<T>::DamarisVar(int dims,
                          const std::vector<std::string>& param_names,
                          const std::vector<int>& param_values,
                          const std::string& variable_name,
                          int rank)
    : param_names_(param_names)
    , variable_name_(variable_name)
    , rank_(rank)
{
    DamarisVar(dims, param_names, variable_name, rank);
    setDamarisParameterAndShmem(param_values); // Initialise the memory size in the constructor.
}

template<class T>
DamarisVar<T>::~DamarisVar()
{
    if (data_ptr_ != nullptr) {
        commitVariableDamarisShmem();
        clearVariableDamarisShmem();
    }
    if (this->hasError()) {
        printError(); // flush out any error messages
    }
}

template<class T>
void DamarisVar<T>::printError() const
{
    OPM_THROW(std::runtime_error, dam_err_str_);
}

template<class T>
std::string DamarisVar<T>::returnXMLForVariable()
{
    return fmt::format("<variable name=\"{}\" {} />", variable_name_,
                       xml_attributes_.ReturnXMLForVariable());
}

template<class T>
void DamarisVar<T>::setDamarisParameter(const std::vector<int>& paramSizeVal)
{
    assert(paramSizeVal.size() == static_cast<std::size_t>(num_params_));

    bool resbool = true;
    std::size_t total_size = 1;
    for (int varnum = 0; varnum < num_params_; varnum++) {
        param_sizes_[varnum] = paramSizeVal[varnum];
        total_size *= param_sizes_[varnum];

        dam_err_ = damaris_parameter_set(param_names_[varnum].c_str(), &paramSizeVal[varnum], sizeof(int));
        if (dam_err_ != DAMARIS_OK) {
            dam_err_str_ += fmt::format("  ERROR rank = {}: class DamarisVar : damaris_parameter_set(\"{}\""
                                        ", paramSizeVal, sizeof(int));  Damaris error = {}\n",
                                        rank_, param_names_[varnum], damaris_error_string(dam_err_));
            resbool = false;
            has_error_ = true;
        }
    }

    if (resbool == true) {
        parameterIsSet(); // sets parameters_set_ and gets the size of the
                          // variables block storage (as number of elemnts)
    }

    if (total_size > 0) {
        current_size_ = total_size;
    } else {
        dam_err_str_ += fmt::format("  ERROR rank = {}: class DamarisVar::getDataStoreBlockSize() "
                                    "The total size of the variable is  0 - please check "
                                    "input paramSizeVal array.\n", rank_);
        has_error_ = true;
    }

    if (hasError()) {
        printError();
    }
}

template<class T>
void DamarisVar<T>::setDamarisPosition(const std::vector<int64_t>& positionsVals)
{
    assert(positionsVals.size() == static_cast<std::size_t>(num_params_));

    for (int pos_dim = 0; pos_dim < num_params_; pos_dim++) {
        positions_[pos_dim] = positionsVals[pos_dim];
    }
    dam_err_ = damaris_set_position(variable_name_.c_str(), positionsVals.data());
    if (dam_err_ != DAMARIS_OK) {
        dam_err_str_ += fmt::format("  ERROR rank = {}: class DamarisVar : damaris_set_position(\"{}\""
                                    ", positionsVals);  Damaris error = {}\n",
                                    rank_, variable_name_, damaris_error_string(dam_err_));
        has_error_ = true;
    }

    if (hasError()) {
        printError();
    }
}

template<class T>
void DamarisVar<T>::setPointersToDamarisShmem()
{
    if (parameters_set_ == true) {
        // Allocate memory in the shared memory section...
        dam_err_ = damaris_alloc(variable_name_.c_str(), (void**)&data_ptr_);
        if (dam_err_ != DAMARIS_OK) {
            dam_err_str_ += fmt::format("  ERROR rank = {}: class DamarisVar : damaris_alloc(\"{}\""
                                        ", (void **) &ret_ptr), Damaris error = {}\n",
                                        rank_, variable_name_, damaris_error_string(dam_err_));
            has_error_ = true;
        }
    } else {
        dam_err_ = -1;
        dam_err_str_ += fmt::format("  ERROR rank = {}: class DamarisVar : "
                                    "setDamarisParameter() should be "
                                    "called first to define the size of the memory "
                                    "block required for variable: {}\n", rank_, variable_name_);
        has_error_ = true;
    }

    if (hasError()) {
        printError();
    }
}

template<class T>
void DamarisVar<T>::commitVariableDamarisShmem()
{
    // Signal to Damaris we are done writing data for this iteration
    dam_err_ = damaris_commit(variable_name_.c_str());
    if (dam_err_ != DAMARIS_OK) {
        dam_err_str_ += fmt::format("  ERROR rank = {}: class DamarisVar : damaris_commit(\"{}\")"
                                    ", Damaris error = {}\n",
                                    rank_, variable_name_, damaris_error_string(dam_err_));
        has_error_ = true;
    }
}

template<class T>
void DamarisVar<T>::clearVariableDamarisShmem()
{
    // Signal to Damaris it has complete charge of the memory area
    dam_err_ = damaris_clear(variable_name_.c_str());
    if (dam_err_ != DAMARIS_OK) {
        dam_err_str_ += fmt::format("  ERROR rank = {}: class DamarisVar : damaris_clear(\"{}\")"
                                    ", Damaris error = {}\n",
                                    rank_, variable_name_, damaris_error_string(dam_err_));
        has_error_ = true;
    }
    data_ptr_ = nullptr;
}

template<class T>
bool DamarisVar<T>::TestType(const std::string& variable_name)
{
    bool resbool = true;
    // This gets the type of the Damaris XML <variable>'s <layout>
    DAMARIS_TYPE_STR vartype;
    dam_err_ = damaris_get_type(variable_name.c_str(), &vartype);
    if (dam_err_ != DAMARIS_OK) {
        dam_err_str_ = fmt::format("  ERROR rank = {}: DamarisVar::DamarisVar () damaris_get_type(\"{}\""
                                   ", vartype);  Damaris error = {}\n",
                                   rank_, variable_name_, damaris_error_string(dam_err_));
        has_error_ = true;
        return false;
    }
    T test_id;
    const std::type_info& t1 = typeid(test_id);

    auto check = [&variable_name,this](auto td)
    {
        const std::type_info& t2 = typeid(td);
        if (t1 != t2) {
            formatTypeError(variable_name, t1.name(), t2.name());
            return false;
        }
        return true;
    };

    if (vartype == DAMARIS_TYPE_DOUBLE) {
        resbool = check(double{});
    } else if (vartype == DAMARIS_TYPE_FLOAT) {
        resbool = check(float{});
    } else if (vartype == DAMARIS_TYPE_CHAR) {
        resbool = check(char{});
    } else if (vartype == DAMARIS_TYPE_UCHAR) {
        using uchar = unsigned char;
        resbool = check(uchar{});
    } else if (vartype == DAMARIS_TYPE_SHORT) {
        resbool = check(short{});
    } else if (vartype == DAMARIS_TYPE_USHORT) {
        using ushort = unsigned short;
        resbool = check(ushort{});
    } else if (vartype == DAMARIS_TYPE_INT) {
        resbool = check(int{});
    } else if (vartype == DAMARIS_TYPE_UINT) {
        using uint = unsigned int;
        resbool = check(uint{});
    } else if (vartype == DAMARIS_TYPE_LONG) {
        resbool = check(long{});
    } else if (vartype == DAMARIS_TYPE_ULONG) {
        using ulong = unsigned long;
        resbool = check(ulong{});
    } else if (vartype == DAMARIS_TYPE_UNDEFINED) {
        dam_err_str_ += fmt::format("  ERROR rank = {}: DamarisVar::DamarisVar():: \"{}\""
                                    " has type DAMARIS_TYPE_UNDEFINED\n", rank_, variable_name);
        has_error_ = true;
        resbool = false;
    } else {
        dam_err_str_ += fmt::format("  ERROR rank = {}: DamarisVar::DamarisVar():: \"{}\""
                                    " is not of available type\n", rank_, variable_name);
        has_error_ = true;
        resbool = false;
    }

    return resbool;
}

template<class T>
void DamarisVar<T>::formatTypeError(const std::string& var_name,
                                    const std::string& type_name1,
                                    const std::string& type_name2)
{
    dam_err_str_ += fmt::format("  ERROR rank = {}: DamarisVar::DamarisVar() variable_name: \"{}\""
                                " The template type of Type of DamarisVar<T> in the code: {}"
                                " does not match type in XML: {}\n",
                                rank_, var_name, type_name1, type_name2);
    has_error_ = true;
}

template class DamarisVar<char>;
template class DamarisVar<double>;
template class DamarisVar<int>;

} // namespace Opm::DamarisOutput
