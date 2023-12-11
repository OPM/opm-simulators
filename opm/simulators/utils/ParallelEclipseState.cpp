/*
  Copyright 2019 Equinor AS.

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
#include "opm/input/eclipse/EclipseState/Grid/FieldProps.hpp"
#include <config.h>
#include <opm/simulators/utils/ParallelEclipseState.hpp>

#include <opm/input/eclipse/EclipseState/Grid/FieldData.hpp>

#include <opm/common/ErrorMacros.hpp>

#include <cstddef>
#include <regex>
#include <string>

namespace {
    bool is_FIP(const std::string& keyword)
    {
        return std::regex_match(keyword, std::regex { "FIP[A-Z0-9]{1,5}" });
    }
}

namespace Opm {


ParallelFieldPropsManager::ParallelFieldPropsManager(FieldPropsManager& manager)
    : m_manager(manager)
    , m_comm(Parallel::Communication())
{
}

// EXPERIMENTAL FUNCTION TO ADD COMM AS INPUT
ParallelFieldPropsManager::ParallelFieldPropsManager(FieldPropsManager& manager, Parallel::Communication comm)
    : m_manager(manager)
    , m_comm(comm)
{
}


std::vector<int> ParallelFieldPropsManager::actnum() const
{
    if (m_comm.rank() == 0)
        return m_manager.actnum();

    return{};
}


void ParallelFieldPropsManager::reset_actnum(const std::vector<int>& actnum)
{
    if (m_comm.rank() != 0)
        OPM_THROW(std::runtime_error, "reset_actnum should only be called on root process.");
    m_manager.reset_actnum(actnum);
}


std::vector<double> ParallelFieldPropsManager::porv(bool global) const
{
    std::vector<double> global_porv;
    if (m_comm.rank() == 0)
        global_porv = m_manager.porv(true);

    std::size_t size = global_porv.size();
    m_comm.broadcast(&size, 1, 0);
    global_porv.resize(size);
    m_comm.broadcast(global_porv.data(), size, 0);
    if (global)
        return global_porv;

    std::vector<double> local_porv(this->m_activeSize());
    for (int i = 0; i < m_activeSize(); ++i)
    {
        local_porv[i] = global_porv[this->m_local2Global(i)];
    }
    return local_porv;
}


const std::vector<int>& ParallelFieldPropsManager::get_int(const std::string& keyword) const
{
    auto it = m_intProps.find(keyword);
    if (it == m_intProps.end())
    {
        // Some of the keywords might be defaulted.
        // We will let rank 0 create them and distribute them using get_global_int
        auto data = get_global_int(keyword);
        auto& local_data = const_cast<std::map<std::string, Fieldprops::FieldData<int>>&>(m_intProps)[keyword];
        local_data.data.resize(m_activeSize());
        local_data.value_status.resize(m_activeSize());

        for (int i = 0; i < m_activeSize(); ++i)
        {
            local_data.data[i] = data[m_local2Global(i)];
        }
        return local_data.data;
    }

    return it->second.data;
}

std::vector<int> ParallelFieldPropsManager::get_global_int(const std::string& keyword) const
{
    std::vector<int> result;
    int exceptionThrown{};

    if (m_comm.rank() == 0) {
        try {
            // Recall: FIP* keywords are special.  We care only about the
            // first three characters of the name following the initial
            // three-character "FIP" prefix, hence "substr(0, 6)".
            result = is_FIP(keyword)
                ? this->m_manager.get_global_int(keyword.substr(0, 6))
                : this->m_manager.get_global_int(keyword);
        }
        catch (std::exception& e) {
            exceptionThrown = 1;
            OpmLog::error("No integer property field: " + keyword + " ("+e.what()+")");
            m_comm.broadcast(&exceptionThrown, 1, 0);
            throw e;
        }
    }

    m_comm.broadcast(&exceptionThrown, 1, 0);

    if (exceptionThrown)
        OPM_THROW_NOLOG(std::runtime_error, "No integer property field: " + keyword);

    std::size_t size = result.size();
    m_comm.broadcast(&size, 1, 0);
    result.resize(size);
    m_comm.broadcast(result.data(), size, 0);

    return result;
}


const std::vector<double>& ParallelFieldPropsManager::get_double(const std::string& keyword) const
{
    auto it = m_doubleProps.find(keyword);
    if (it == m_doubleProps.end())
    {
        // Some of the keywords might be defaulted.
        // We will let rank 0 create them and distribute them using get_global_int
        auto data = get_global_double(keyword);
        auto& local_data = const_cast<std::map<std::string, Fieldprops::FieldData<double>>&>(m_doubleProps)[keyword];
        local_data.data.resize(m_activeSize());
        local_data.value_status.resize(m_activeSize());
        for (int i = 0; i < m_activeSize(); ++i)
        {
            local_data.data[i] = data[m_local2Global(i)];
        }
        return local_data.data;
    }

    return it->second.data;
}


std::vector<double> ParallelFieldPropsManager::get_global_double(const std::string& keyword) const
{
    std::vector<double> result;
    int exceptionThrown{};

    if (m_comm.rank() == 0)
    {
        try
        {
            result = m_manager.get_global_double(keyword);
        }catch(std::exception& e) {
            exceptionThrown = 1;
            OpmLog::error("No double property field: " + keyword + " ("+e.what()+")");
            m_comm.broadcast(&exceptionThrown, 1, 0);
            throw e;
        }
    }

    m_comm.broadcast(&exceptionThrown, 1, 0);

    if (exceptionThrown)
        OPM_THROW_NOLOG(std::runtime_error, "No double property field: " + keyword);

    std::size_t size = result.size();
    m_comm.broadcast(&size, 1, 0);
    result.resize(size);
    m_comm.broadcast(result.data(), size, 0);

    return result;
}

bool ParallelFieldPropsManager::tran_active(const std::string& keyword) const
{
    auto calculator = m_tran.find(keyword);
    return calculator != m_tran.end() && calculator->second.size();
}

void ParallelFieldPropsManager::apply_tran(const std::string& keyword,
                                           std::vector<double>& data) const
{
    Opm::apply_tran(m_tran, m_doubleProps, m_activeSize(), keyword, data);
}

bool ParallelFieldPropsManager::has_int(const std::string& keyword) const
{
    auto it = m_intProps.find(keyword);
    return it != m_intProps.end();
}

bool ParallelFieldPropsManager::has_double(const std::string& keyword) const
{
    auto it = m_doubleProps.find(keyword);
    return it != m_doubleProps.end();
}

std::vector<std::string> ParallelFieldPropsManager::fip_regions() const
{
    constexpr auto maxchar = std::string::size_type{6};

    std::vector<std::string> result;
    for (const auto& key : m_intProps) {
        if (Fieldprops::keywords::isFipxxx(key.first)) {
            result.push_back(key.first.substr(0, maxchar));
        }
    }
    return result;
}


ParallelEclipseState::ParallelEclipseState(Parallel::Communication comm)
    : m_fieldProps(field_props, comm)
    , m_comm(comm)
{
}


ParallelEclipseState::ParallelEclipseState(const Deck& deck)
    : EclipseState(deck)
    , m_fieldProps(field_props)
{
}

ParallelEclipseState::ParallelEclipseState(const Deck& deck, Parallel::Communication comm)
    : EclipseState(deck)
    , m_fieldProps(field_props, comm)
    , m_comm(comm)
{
}

const FieldPropsManager& ParallelEclipseState::fieldProps() const
{
    if (!m_parProps && m_comm.rank() != 0)
        OPM_THROW(std::runtime_error, "Attempt to access field properties on no-root process before switch to parallel properties");

    if (!m_parProps || m_comm.size() == 1)
        return this->EclipseState::fieldProps();

    return m_fieldProps;
}


const FieldPropsManager& ParallelEclipseState::globalFieldProps() const
{
    if (m_comm.rank() != 0)
        OPM_THROW(std::runtime_error, "Attempt to access global field properties on non-root process");
    return this->EclipseState::globalFieldProps();
}


const EclipseGrid& ParallelEclipseState::getInputGrid() const
{
    if (m_comm.rank() != 0)
        OPM_THROW(std::runtime_error, "Attempt to access eclipse grid on non-root process");
    return this->EclipseState::getInputGrid();
}


void ParallelEclipseState::switchToGlobalProps()
{
    m_parProps = false;
}


void ParallelEclipseState::switchToDistributedProps()
{
    if (m_comm.size() == 1) // No need for the parallel frontend
        return;
    m_parProps = true;
}

} // end namespace Opm
