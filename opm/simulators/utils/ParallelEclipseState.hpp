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
#ifndef PARALLEL_ECLIPSE_STATE_HPP
#define PARALLEL_ECLIPSE_STATE_HPP

#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <dune/common/parallel/mpihelper.hh>

namespace Opm {


class EclMpiSerializer;

/*! \brief Parallel frontend to the field properties.
 *
 * \details This is a parallel frontend to the mpi-unaware
 *          FieldPropsManager in opm-common. It contains
 *          process-local field properties on each process using
 *          compressed indexing.
*/

class ParallelFieldPropsManager : public FieldPropsManager {
public:
    friend class ParallelEclipseState; //!< Friend so props can be setup.
    //! \brief Friend to set up props
    template<class Grid>
    friend class PropsCentroidsDataHandle;

    //! \brief Constructor.
    //! \param manager The field property manager to wrap.
    ParallelFieldPropsManager(FieldPropsManager& manager);

    //! \brief Returns actnum vector.
    //! \details If called on non-root process an empty vector is returned
    std::vector<int> actnum() const override;

    //! \brief Reset the actnum vector.
    //! \details Can only be called on root process
    void reset_actnum(const std::vector<int>& actnum) override;

    //! \brief Returns the pore volume vector.
    std::vector<double> porv(bool global = false) const override;

    //! \brief Returns an int property using compressed indices.
    //! \param keyword Name of property
    const std::vector<int>& get_int(const std::string& keyword) const override;

    //! \brief Returns a double property using compressed indices.
    //! \param keyword Name of property
    const std::vector<double>& get_double(const std::string& keyword) const override;

    //! \brief Returns an int property using global cartesian indices.
    //! \param keyword Name of property
    //! \details The vector is broadcast from root process
    std::vector<int> get_global_int(const std::string& keyword) const override;

    //! \brief Returns a double property using global cartesian indices.
    //! \param keyword Name of property
    //! \details The vector is broadcast from root process
    std::vector<double> get_global_double(const std::string& keyword) const override;

    //! \brief Check if an integer property is available.
    //! \param keyword Name of property
    bool has_int(const std::string& keyword) const override;

    //! \brief Check if a double property is available.
    //! \param keyword Name of property
    bool has_double(const std::string& keyword) const override;

protected:
    std::map<std::string, std::vector<int>> m_intProps; //!< Map of integer properties in process-local compressed indices.
    std::map<std::string, std::vector<double>> m_doubleProps; //!< Map of double properties in process-local compressed indices.
    FieldPropsManager& m_manager; //!< Underlying field property manager (only used on root process).
    Dune::CollectiveCommunication<Dune::MPIHelper::MPICommunicator> m_comm; //!< Collective communication handler.
};


/*! \brief Parallel frontend to the EclipseState
 *
 * \details This is a parallel frontend to the mpi-unaware EclipseState in opm-common.
 *          It extends the eclipse state class with serialization support, and
 *          contains methods to switch between full global field properties,
 *          and distributed field properties for consumption in the simulator.
 *          Additionally, it has a few sanity checks to ensure that the data that
 *          is only available on the root process is not attempted to be accessed
 *          on non-root processes.
*/

class ParallelEclipseState : public EclipseState {
    //! \brief Friend to set up props
    template<class Grid>
    friend class PropsCentroidsDataHandle;
public:
    //! \brief Default constructor.
    ParallelEclipseState();

    //! \brief Construct from a deck instance.
    //! \param deck The deck to construct from
    //! \details Only called on root process
    ParallelEclipseState(const Deck& deck);

    //! \brief Switch to global field properties.
    //! \details Called on root process to use the global field properties
    void switchToGlobalProps();

    //! \brief Switch to distributed field properies.
    //! \details Called on root process to use the distributed field properties.
    //!          setupLocalProps must be called prior to this.
    void switchToDistributedProps();

    //! \brief Returns a const ref to current field properties.
    const FieldPropsManager& fieldProps() const override;

    //! \brief Returns a const ref to global field properties.
    //! \details Can only be called on root process.
    const FieldPropsManager& globalFieldProps() const override;

    //! \brief Returns a const ref to the eclipse grid.
    //! \details Can only be called on root process.
    const EclipseGrid& getInputGrid() const override;

private:
    bool m_parProps = false; //! True to use distributed properties on root process
    ParallelFieldPropsManager m_fieldProps; //!< The parallel field properties
};


} // end namespace Opm
#endif // PARALLEL_ECLIPSE_STATE_HPP
