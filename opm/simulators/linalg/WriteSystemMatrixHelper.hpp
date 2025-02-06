/*
  Copyright 2019 SINTEF Digital, Mathematics and Cybernetics.

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

#ifndef OPM_WRITESYSTEMMATRIXHELPER_HEADER_INCLUDED
#define OPM_WRITESYSTEMMATRIXHELPER_HEADER_INCLUDED

#include <dune/istl/matrixmarket.hh>

#include <opm/simulators/linalg/MatrixMarketSpecializations.hpp>

#include <filesystem>
#include <iomanip>
#include <sstream>
#include <string>

namespace Opm::Helper {

    namespace detail {
        /// Output a linear algebra object to file in MatrixMarket format.
        ///
        /// Common implementation function.
        ///
        /// \tparam SimulatorType Packaged simulator type.  Expected to
        /// provide, among other services, an episode index, a simulation
        /// time, and a physical model.
        ///
        /// \tparam LinalgObjectType Object type for a linear algebra
        /// object.  Expected to be a coefficient matrix type or a linear
        /// system vector type from the Dune ISTL.
        ///
        /// \tparam Communicator ISTL communication type.  Typically
        /// \code Dune::OwnerOverlapCopyCommunication<> \endcode or similar.
        ///
        /// \param[in] simulator Simulator object.
        ///
        /// \param[in] linalgObject Linear algebra object.  Expected to be a
        /// Dune ISTL matrix or vector.
        ///
        /// \param[in] objName Name of linear algebra object.  Incorporated
        /// into the name of the MatrixMarket output file.
        ///
        /// \param[in] comm Dune ISTL communication object.  Unused in
        /// builds without MPI support.  You may use a nullptr value in a
        /// sequential run of an MPI-enabled build to signify that parallel
        /// output is unneeded.
        template <class SimulatorType, class LinalgObjectType, class Communicator>
        void writeMatrixMarketObject(const SimulatorType&                 simulator,
                                     const LinalgObjectType&              linalgObject,
                                     const std::string&                   objName,
                                     [[maybe_unused]] const Communicator* comm)
        {
            namespace fs = std::filesystem;

            const auto output_dir = fs::path {simulator.problem().outputDir()} / "reports";
            if (! fs::exists(output_dir)) {
                fs::create_directory(output_dir);
            }

            // Combine and return.
            std::ostringstream oss;
            oss << "prob_"  << simulator.episodeIndex()
                << "_time_" << std::setprecision(15) << std::setw(12) << std::setfill('0') << simulator.time()
                << "_nit_"  << simulator.model().newtonMethod().numIterations()
                << '_'      << objName << "_istl";

            const auto filename = (output_dir / oss.str()).generic_string();

#if HAVE_MPI
            if (comm != nullptr) {
                Dune::storeMatrixMarket(linalgObject, filename, *comm);
            }
            else
#endif // HAVE_MPI
            {
                Dune::storeMatrixMarket(linalgObject, filename + ".mm");
            }
        }
    } // namespace detail

    /// Output a Dune ISTL vector to file in MatrixMarket format.
    ///
    /// \tparam SimulatorType Packaged simulator type.  Expected to
    /// provide, among other services, an episode index, a simulation
    /// time, and a physical model.
    ///
    /// \tparam VectorType Object type for a linear system vector object.
    /// Typically a linear system vector type from the Dune ISTL.
    ///
    /// \tparam Communicator ISTL communication type.  Typically
    /// \code Dune::OwnerOverlapCopyCommunication<> \endcode or similar.
    ///
    /// \param[in] simulator Simulator object.
    ///
    /// \param[in] rhs Linear system right-hand side vector.  Expected to be
    /// a Dune ISTL vector.
    ///
    /// \param[in] sysName Name of linear system/physical model.
    /// Incorporated into the name of the MatrixMarket output file.
    ///
    /// \param[in] comm Dune ISTL communication object.  Unused in builds
    /// without MPI support.  You may use a nullptr value in a sequential
    /// run of an MPI-enabled build to signify that parallel output is
    /// unneeded.
    template <class SimulatorType, class VectorType, class Communicator>
    void writeVector(const SimulatorType& simulator,
                     const VectorType&    rhs,
                     const std::string&   sysName,
                     const Communicator*  comm)
    {
        detail::writeMatrixMarketObject(simulator, rhs, sysName + "_vector", comm);
    }

    /// Output a Dune ISTL matrix to file in MatrixMarket format.
    ///
    /// \tparam SimulatorType Packaged simulator type.  Expected to
    /// provide, among other services, an episode index, a simulation
    /// time, and a physical model.
    ///
    /// \tparam MatrixType Object type for a linear system matrix object.
    /// Typically a linear system matrix type from the Dune ISTL.
    ///
    /// \tparam Communicator ISTL communication type.  Typically
    /// \code Dune::OwnerOverlapCopyCommunication<> \endcode or similar.
    ///
    /// \param[in] simulator Simulator object.
    ///
    /// \param[in] matrix Linear system coefficient matrix.  Expected to be
    /// a Dune ISTL (sparse) matrix.
    ///
    /// \param[in] sysName Name of linear system/physical model.
    /// Incorporated into the name of the MatrixMarket output file.
    ///
    /// \param[in] comm Dune ISTL communication object.  Unused in builds
    /// without MPI support.  You may use a nullptr value in a sequential
    /// run of an MPI-enabled build to signify that parallel output is
    /// unneeded.
    template <class SimulatorType, class MatrixType, class Communicator>
    void writeMatrix(const SimulatorType& simulator,
                     const MatrixType&    matrix,
                     const std::string&   sysName,
                     const Communicator*  comm)
    {
        detail::writeMatrixMarketObject(simulator, matrix, sysName + "_matrix", comm);
    }

    /// Output a Dune ISTL linear system to files in MatrixMarket format.
    ///
    /// This function will create one output file for the coefficient matrix
    /// and another output file for the system right-hand side vector.
    ///
    /// \tparam SimulatorType Packaged simulator type.  Expected to
    /// provide, among other services, an episode index, a simulation
    /// time, and a physical model.
    ///
    /// \tparam MatrixType Object type for a linear system matrix object.
    /// Typically a linear system matrix type from the Dune ISTL.
    ///
    /// \tparam VectorType Object type for a linear system vector object.
    /// Typically a linear system vector type from the Dune ISTL.
    ///
    /// \tparam Communicator ISTL communication type.  Typically
    /// \code Dune::OwnerOverlapCopyCommunication<> \endcode or similar.
    ///
    /// \param[in] simulator Simulator object.
    ///
    /// \param[in] matrix Linear system coefficient matrix.  Expected to be
    /// a Dune ISTL (sparse) matrix.
    ///
    /// \param[in] rhs Linear system right-hand side vector.  Expected to be
    /// a Dune ISTL vector.
    ///
    /// \param[in] sysName Name of linear system/physical model.
    /// Incorporated into the name of the MatrixMarket output file.
    ///
    /// \param[in] comm Dune ISTL communication object.  Unused in builds
    /// without MPI support.  You may use a nullptr value in a sequential
    /// run of an MPI-enabled build to signify that parallel output is
    /// unneeded.
    template <class SimulatorType, class MatrixType, class VectorType, class Communicator>
    void writeSystem(const SimulatorType& simulator,
                     const MatrixType&    matrix,
                     const VectorType&    rhs,
                     const std::string&   sysName,
                     const Communicator*  comm)
    {
        writeMatrix(simulator, matrix, sysName, comm);
        writeVector(simulator, rhs, sysName, comm);
    }

    /// Output a Dune ISTL linear system from linearising a set of flow
    /// equations to files in MatrixMarket format.
    ///
    /// This function will create one output file for the coefficient matrix
    /// and another output file for the system right-hand side vector.
    ///
    /// \tparam SimulatorType Packaged simulator type.  Expected to
    /// provide, among other services, an episode index, a simulation
    /// time, and a physical model.
    ///
    /// \tparam MatrixType Object type for a linear system matrix object.
    /// Typically a linear system matrix type from the Dune ISTL.
    ///
    /// \tparam VectorType Object type for a linear system vector object.
    /// Typically a linear system vector type from the Dune ISTL.
    ///
    /// \tparam Communicator ISTL communication type.  Typically
    /// \code Dune::OwnerOverlapCopyCommunication<> \endcode or similar.
    ///
    /// \param[in] simulator Simulator object.
    ///
    /// \param[in] matrix Linear system coefficient matrix.  Expected to be
    /// a Dune ISTL (sparse) matrix.
    ///
    /// \param[in] rhs Linear system right-hand side vector.  Expected to be
    /// a Dune ISTL vector.
    ///
    /// \param[in] comm Dune ISTL communication object.  Unused in builds
    /// without MPI support.  You may use a nullptr value in a sequential
    /// run of an MPI-enabled build to signify that parallel output is
    /// unneeded.
    template <class SimulatorType, class MatrixType, class VectorType, class Communicator>
    void writeSystem(const SimulatorType& simulator,
                     const MatrixType&    matrix,
                     const VectorType&    rhs,
                     const Communicator*  comm)
    {
        writeSystem(simulator, matrix, rhs, "flow_", comm);
    }

    /// Output a Dune ISTL linear system from linearising a set of
    /// geo-mechanical equations to files in MatrixMarket format.
    ///
    /// This function will create one output file for the coefficient matrix
    /// and another output file for the system right-hand side vector.
    ///
    /// \tparam SimulatorType Packaged simulator type.  Expected to
    /// provide, among other services, an episode index, a simulation
    /// time, and a physical model.
    ///
    /// \tparam MatrixType Object type for a linear system matrix object.
    /// Typically a linear system matrix type from the Dune ISTL.
    ///
    /// \tparam VectorType Object type for a linear system vector object.
    /// Typically a linear system vector type from the Dune ISTL.
    ///
    /// \tparam Communicator ISTL communication type.  Typically
    /// \code Dune::OwnerOverlapCopyCommunication<> \endcode or similar.
    ///
    /// \param[in] simulator Simulator object.
    ///
    /// \param[in] matrix Linear system coefficient matrix.  Expected to be
    /// a Dune ISTL (sparse) matrix.
    ///
    /// \param[in] rhs Linear system right-hand side vector.  Expected to be
    /// a Dune ISTL vector.
    ///
    /// \param[in] comm Dune ISTL communication object.  Unused in builds
    /// without MPI support.  You may use a nullptr value in a sequential
    /// run of an MPI-enabled build to signify that parallel output is
    /// unneeded.
    template <class SimulatorType, class MatrixType, class VectorType, class Communicator>
    void writeMechSystem(const SimulatorType& simulator,
                         const MatrixType&    matrix,
                         const VectorType&    rhs,
                         const Communicator*  comm)
    {
        writeSystem(simulator, matrix, rhs, "mech_", comm);
    }

} // namespace Opm::Helper

#endif // OPM_WRITESYSTEMMATRIXHELPER_HEADER_INCLUDED
