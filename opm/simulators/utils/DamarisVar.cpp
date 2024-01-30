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

#include <cassert>
#include <sstream>
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
    std::ostringstream var_sstr;

    var_sstr << " layout=\"" << this->layout_ << "\"";
    if (this->mesh_ != "")
        var_sstr << " mesh=\"" << this->mesh_ << "\"";
    if (this->type_ != "")
        var_sstr << " type=\"" << this->type_ << "\"";
    if (this->visualizable_ != "")
        var_sstr << " visualizable=\"" << this->visualizable_ << "\"";
    if (this->unit_ != "")
        var_sstr << " unit=\"" << this->unit_ << "\"";
    if (this->time_varying_ != "")
        var_sstr << " time_varying=\"" << this->time_varying_ << "\"";
    if (this->centering_ != "")
        var_sstr << " centering=\"" << this->centering_ << "\"";
    if (this->store_ != "")
        var_sstr << " store=\"" << this->store_ << "\"";
    if (this->script_ != "")
        var_sstr << " script=\"" << this->script_ << "\"";
    if (this->select_mem_ != "")
        var_sstr << " select-mem=\"" << this->select_mem_ << "\"";
    if (this->select_file_ != "")
        var_sstr << " select-file=\"" << this->select_file_ << "\"";
    if (this->select_subset_ != "")
        var_sstr << " select-subset=\"" << this->select_subset_ << "\"";

    return var_sstr.str();
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
    if (this->hasError())
        printError(); // flush out any error messages
}

template<class T>
void DamarisVar<T>::printError() const
{
    OPM_THROW(std::runtime_error, dam_err_sstr_.str());
}

template<class T>
std::string DamarisVar<T>::returnXMLForVariable()
{
    std::ostringstream var_sstr;

    var_sstr << "<variable "
             << " name=\"" << variable_name_ << "\"";
    var_sstr << xml_attributes_.ReturnXMLForVariable();
    var_sstr << " /> ";

    return var_sstr.str();
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
            dam_err_sstr_ << "  ERROR rank =" << rank_ << " : class DamarisVar : damaris_parameter_set(\""
                          << param_names_[varnum] << "\", paramSizeVal, sizeof(int));  Damaris error = "
                          << damaris_error_string(dam_err_) << std::endl;
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
        dam_err_sstr_ << "  ERROR rank =" << rank_ << " : class DamarisVar::getDataStoreBlockSize() "
                      << "The total size of the variable is  0 - please check "
                         "input paramSizeVal array."
                      << std::endl;
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
        dam_err_sstr_ << "  ERROR rank =" << rank_ << " : class DamarisVar : damaris_set_position(\""
                      << variable_name_
                      << "\", positionsVals);  Damaris error = " << damaris_error_string(dam_err_) << std::endl;
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
            dam_err_sstr_ << "  ERROR rank =" << rank_ << " : class DamarisVar : damaris_alloc(\""
                          << variable_name_ << "\", (void **) &ret_ptr)"
                          << ", Damaris error = " << damaris_error_string(dam_err_) << std::endl;
            has_error_ = true;
        }
    } else {
        dam_err_ = -1;
        dam_err_sstr_ << "  ERROR rank =" << rank_
                      << " : class DamarisVar : setDamarisParameter() should be "
                         "called first so as to define the size of the memory "
                         "block required for variable : "
                      << variable_name_ << std::endl;
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
        dam_err_sstr_ << "  ERROR rank =" << rank_ << " : class DamarisVar : damaris_commit(\""
                      << variable_name_ << "\")"
                      << ", Damaris error = " << damaris_error_string(dam_err_) << std::endl;
        has_error_ = true;
    }
}

template<class T>
void DamarisVar<T>::clearVariableDamarisShmem()
{
    // Signal to Damaris it has complete charge of the memory area
    dam_err_ = damaris_clear(variable_name_.c_str());
    if (dam_err_ != DAMARIS_OK) {
        dam_err_sstr_ << "  ERROR rank =" << rank_ << " : class DamarisVar : damaris_clear(\"" << variable_name_
                      << "\")"
                      << ", Damaris error = " << damaris_error_string(dam_err_) << std::endl;
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
        dam_err_sstr_ << "  ERROR rankDamarisVar::DamarisVar ()  damaris_get_type(\"" << variable_name_
                      << "\", vartype);  Damaris error = " << damaris_error_string(dam_err_) << std::endl;
        has_error_ = true;
        return false;
    }
    T test_id;
    const std::type_info& t1 = typeid(test_id);

    if (vartype == DAMARIS_TYPE_DOUBLE) {
        double td = 0.0;
        const std::type_info& t2 = typeid(td);
        if (t1 != t2) {
            formatTypeError(variable_name, t1.name(), t2.name());
            resbool = false;
        }
    } else if (vartype == DAMARIS_TYPE_FLOAT) {
        float td = 0.0f;
        const std::type_info& t2 = typeid(td);
        if (t1 != t2) {
            formatTypeError(variable_name, t1.name(), t2.name());
            resbool = false;
        }
    } else if (vartype == DAMARIS_TYPE_CHAR) {
        char td = 0;
        const std::type_info& t2 = typeid(td);
        if (t1 != t2) {
            formatTypeError(variable_name, t1.name(), t2.name());
            resbool = false;
        }
    } else if (vartype == DAMARIS_TYPE_UCHAR) {
        unsigned char td = 0;
        const std::type_info& t2 = typeid(td);
        if (t1 != t2) {
            formatTypeError(variable_name, t1.name(), t2.name());
            resbool = false;
        }
    } else if (vartype == DAMARIS_TYPE_SHORT) {
        short td = 0;
        const std::type_info& t2 = typeid(td);
        if (t1 != t2) {
            formatTypeError(variable_name, t1.name(), t2.name());
            resbool = false;
        }
    } else if (vartype == DAMARIS_TYPE_USHORT) {
        unsigned short td = 0;
        const std::type_info& t2 = typeid(td);
        if (t1 != t2) {
            formatTypeError(variable_name, t1.name(), t2.name());
            resbool = false;
        }
    } else if (vartype == DAMARIS_TYPE_INT) {
        int td = 0;
        const std::type_info& t2 = typeid(td);
        if (t1 != t2) {
            formatTypeError(variable_name, t1.name(), t2.name());
            resbool = false;
        }
    } else if (vartype == DAMARIS_TYPE_UINT) {
        unsigned int td = 0;
        const std::type_info& t2 = typeid(td);
        if (t1 != t2) {
            formatTypeError(variable_name, t1.name(), t2.name());
            resbool = false;
        }
    } else if (vartype == DAMARIS_TYPE_LONG) {
        long td = 0;
        const std::type_info& t2 = typeid(td);
        if (t1 != t2) {
            formatTypeError(variable_name, t1.name(), t2.name());
            resbool = false;
        }
    } else if (vartype == DAMARIS_TYPE_ULONG) {
        unsigned long td = 0;
        const std::type_info& t2 = typeid(td);
        if (t1 != t2) {
            formatTypeError(variable_name, t1.name(), t2.name());
            resbool = false;
        }
    } else if (vartype == DAMARIS_TYPE_UNDEFINED) {
        dam_err_sstr_ << "  ERROR rank =" << rank_ << " : DamarisVar::DamarisVar()::  \"" << variable_name
                      << "\" has type DAMARIS_TYPE_UNDEFINED" << std::endl;
        has_error_ = true;
        resbool = false;
    } else {
        dam_err_sstr_ << "  ERROR rank =" << rank_ << " : DamarisVar::DamarisVar():: \"" << variable_name
                      << "\" is not of available type " << std::endl;
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
    dam_err_sstr_ << "  ERROR rank =" << rank_ << " : DamarisVar::DamarisVar () variable_name_: \"" << var_name
                  << "\" The template type of Type of DamarisVar<T> in the code: " << type_name1
                  << " does not match type in XML:" << type_name2 << std::endl;
    has_error_ = true;
}

template class DamarisVar<char>;
template class DamarisVar<double>;
template class DamarisVar<int>;

} // namespace Opm::DamarisOutput
