// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2022 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2023 Inria, Bretagne–Atlantique Research Center

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \copydoc Opm::DamarisWriter
 */
#ifndef OPM_DAMARIS_WRITER_HPP
#define OPM_DAMARIS_WRITER_HPP

#include <dune/grid/common/partitionset.hh>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/simulators/flow/countGlobalCells.hpp>
#include <opm/simulators/flow/DamarisParameters.hpp>
#include <opm/simulators/flow/EclGenericWriter.hpp>
#include <opm/simulators/flow/FlowBaseVanguard.hpp>
#include <opm/simulators/flow/OutputBlackoilModule.hpp>
#include <opm/simulators/utils/DamarisVar.hpp>
#include <opm/simulators/utils/DamarisKeywords.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/utils/GridDataOutput.hpp>
#include <opm/simulators/utils/ParallelSerialization.hpp>

#include <fmt/format.h>

#include <algorithm>
#include <memory>
#include <numeric>
#include <string>
#include <vector>
#include <unordered_set>

namespace Opm {

namespace DamarisOutput {

int endIteration();
int setParameter(const char* field, int value);
int setPosition(const char* field, int64_t pos);
int write(const char* field, const void* data);
int setupWritingPars(Parallel::Communication comm,
                     const int n_elements_local_grid,
                     std::vector<unsigned long long>& elements_rank_offsets);
void handleError(const int dam_err, Parallel::Communication comm, const std::string& message);
}

/*!
 * \ingroup EclBlackOilSimulator
 *
 * \brief Collects necessary output values and pass them to Damaris server processes.
 *
 * Currently only passing through PRESSURE, GLOBAL_CELL_INDEX and MPI_RANK information.
 * This class now passes through the 3D mesh information to Damaris to enable
 * in situ visualization via Paraview or Ascent. And developed so that variables specified
 * through the Eclipse input deck will be available to Damaris.
 */


template <class TypeTag>
class DamarisWriter : public EclGenericWriter<GetPropType<TypeTag, Properties::Grid>,
                                              GetPropType<TypeTag, Properties::EquilGrid>,
                                              GetPropType<TypeTag, Properties::GridView>,
                                              GetPropType<TypeTag, Properties::ElementMapper>,
                                              GetPropType<TypeTag, Properties::Scalar>>
{
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using EquilGrid = GetPropType<TypeTag, Properties::EquilGrid>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementMapper = GetPropType<TypeTag, Properties::ElementMapper>;

    using BaseType = EclGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>;
    using DamarisVarInt = DamarisOutput::DamarisVar<int>;
    using DamarisVarChar = DamarisOutput::DamarisVar<char>;
    using DamarisVarDbl = DamarisOutput::DamarisVar<double>;

public:
    static void registerParameters()
    {
        registerDamarisParameters();
    }

    // The Simulator object should preferably have been const - the
    // only reason that is not the case is due to the SummaryState
    // object owned deep down by the vanguard.
    explicit DamarisWriter(Simulator& simulator)
        : BaseType(simulator.vanguard().schedule(),
                   simulator.vanguard().eclState(),
                   simulator.vanguard().summaryConfig(),
                   simulator.vanguard().grid(),
                   ((simulator.vanguard().grid().comm().rank() == 0)
                    ? &simulator.vanguard().equilGrid()
                    : nullptr),
                   simulator.vanguard().gridView(),
                   simulator.vanguard().cartesianIndexMapper(),
                   ((simulator.vanguard().grid().comm().rank() == 0)
                    ? &simulator.vanguard().equilCartesianIndexMapper()
                    : nullptr),
                   false, false)
        , simulator_(simulator)
    {
        this->damarisUpdate_ = true ;

        this->rank_ = this->simulator_.vanguard().grid().comm().rank() ;
        this->nranks_ = this->simulator_.vanguard().grid().comm().size();

        this->elements_rank_offsets_.resize(this->nranks_);

        // Get the size of the unique vector elements (excludes the shared 'ghost' elements)
        //
        // Might possibly use
        //
        //     detail::countLocalInteriorCellsGridView(this->simulator_.gridView())
        //
        // from countGlobalCells.hpp instead of calling std::distance() directly.
        {
            const auto& gridView = this->simulator_.gridView();
            const auto& interior_elements = elements(gridView, Dune::Partitions::interior);

            this->numElements_ = std::distance(interior_elements.begin(), interior_elements.end());
        }

        if (this->nranks_ > 1) {
            auto smryCfg = (this->rank_ == 0)
                ? this->eclIO_->finalSummaryConfig()
                : SummaryConfig{};

            eclBroadcast(this->simulator_.vanguard().grid().comm(), smryCfg);

            this->damarisOutputModule_ = std::make_unique<OutputBlackOilModule<TypeTag>>
                (simulator, smryCfg, this->collectOnIORank_);
        }
        else {
            this->damarisOutputModule_ = std::make_unique<OutputBlackOilModule<TypeTag>>
                (simulator, this->eclIO_->finalSummaryConfig(), this->collectOnIORank_);
        }

        wanted_vars_set_ = DamarisOutput::getSetOfIncludedVariables();
    }

    /*!
     * \brief Writes localCellData through to Damaris servers. Sets up the unstructured mesh which is passed to Damaris.
     */
    void writeOutput(data::Solution& localCellData , bool isSubStep)
    {
        OPM_TIMEBLOCK(writeOutput);
        const int reportStepNum = simulator_.episodeIndex() + 1;
        const auto& cc = simulator_.vanguard().grid().comm();

        // added this as localCellData was not being written
        if (!isSubStep)
            this->damarisOutputModule_->invalidateLocalData() ;
        this->prepareLocalCellData(isSubStep, reportStepNum);
        this->damarisOutputModule_->outputErrorLog(cc);

        // The damarisWriter is not outputing well or aquifer data (yet)
        auto localWellData = simulator_.problem().wellModel().wellData(); // data::Well

        if (! isSubStep)
        {
            if (localCellData.size() == 0) {
                this->damarisOutputModule_->assignToSolution(localCellData);
            }

            // add cell data to perforations for Rft output
            this->damarisOutputModule_->addRftDataToWells(localWellData, reportStepNum, cc);

            // On first call and if the mesh and variable size change then set damarisUpdate_ to true
            if (damarisUpdate_ == true) {
                // Sets the damaris parameter values "n_elements_local" and "n_elements_total"
                // which define sizes of the Damaris variables, per-rank and globally (over all ranks).
                // Also sets the offsets to where a ranks array data sits within the global array.
                // This is usefull for HDF5 output and for defining distributed arrays in Dask.
                dam_err_ = DamarisOutput::setupWritingPars(cc, numElements_, elements_rank_offsets_);

                // sets positions and data for non-time-varying variables MPI_RANK and GLOBAL_CELL_INDEX
                this->setGlobalIndexForDamaris() ;

                // Set the geometry data for the mesh model.
                // this function writes the mesh data directly to Damaris shared memory using Opm::DamarisOutput::DamarisVar objects.
                this->writeDamarisGridOutput() ;

                // Currently by default we assume a static mesh grid (the geometry unchanging through the simulation)
                // Set damarisUpdate_ to true if we want to update the geometry sent to Damaris
                this->damarisUpdate_ = false;
            }


            int cell_data_written = 0 ;
            // Call damaris_write() for all available variables
            for ( const auto& damVar : localCellData )
            {
               // std::map<std::string, data::CellData>
              const std::string& name = damVar.first ;
              // If the wanted_vars_set_ set is empty then we default to passing through
              // all the variables/
              if (wanted_vars_set_.count(name) || wanted_vars_set_.empty())
              {
                  const data::CellData&  dataCol = damVar.second ;
                  OpmLog::debug(fmt::format(fmt::runtime("Name of Damaris Variable       : ( rank:{})  name: {}  "),  rank_, name));

                  // Call damaris_set_position() for all available variables
                  // There is an assumption that all variables are the same size, with the same offset.
                  // see initDamarisTemplateXmlFile.cpp for the Damaris XML descriptions.
                  dam_err_ = DamarisOutput::setPosition(name.c_str(), this->elements_rank_offsets_[rank_]);

                  // It does not seem I can test for what type of data is present (double or int)
                  // in the std::variant within the data::CellData, so I will use a try catch block.
                  try {
                    if (dataCol.data<double>().size() >= static_cast<std::vector<double>::size_type>(this->numElements_)) {
                        dam_err_ = DamarisOutput::write(name.c_str(), dataCol.data<double>().data()) ;
                    } else {
                        OpmLog::info(fmt::format(fmt::runtime("( rank:{}) The variable \"{}\" was found to be of a different size {} (not {})."),  rank_, name, dataCol.data<double>().size(), this->numElements_ ));
                    }
                  }
                  catch (std::bad_variant_access const& ex) {
                    // Not a std::vector<double>, must be a std::vector<int>
                    if (dataCol.data<int>().size() >= static_cast<std::vector<int>::size_type>(this->numElements_)) {
                        dam_err_ = DamarisOutput::write(name.c_str(), dataCol.data<int>().data()) ;
                    } else {
                        OpmLog::info(fmt::format(fmt::runtime("( rank:{}) The variable \"{}\" was found to be of a different size {} (not {})."),  rank_, name, dataCol.data<int>().size(), this->numElements_ ));
                    }
                  }
                  ++cell_data_written ;
              }
            }
            DamarisOutput::handleError(dam_err_, cc, "setPosition() and write() for available variables");

            if (!cell_data_written) {
                  OpmLog::info(fmt::format(fmt::runtime("( rank:{}) No simulation data written to the Damaris server - check --damaris-limit-variables command line option (if used) has valid variable name(s) and that the Damaris XML file contains variable names that are available in your simulation."),  rank_));
            } else {
                  OpmLog::debug(fmt::format(fmt::runtime("( rank:{}) {} Damaris Variables written to the Damaris servers"),  rank_, cell_data_written));
            }

           /*
            Code for when we want to pass to Damaris the single cell 'block' data variables
            auto mybloc = damarisOutputModule_->getBlockData() ;
            for ( auto damVar : mybloc ) {
               // std::map<std::string, data::CellData>
              const std::string name = std::get<0>(damVar.first) ;
              const int part = std::get<1>(damVar.first) ;
              double  dataCol = damVar.second ;
              std::cout << "Name of Damaris Block Varaiable : (" << rank_ << ")  "  << name  << "  part : " << part << "  Value : "  << dataCol <<  std::endl ;
            }

            dam_err_ =  DamarisOutput::endIteration();
            */
            if (!this->damarisOutputModule_->getFluidPressure().empty())
            {
                dam_err_ =  DamarisOutput::endIteration();
            }
            DamarisOutput::handleError(dam_err_, cc, "endIteration()");

         } // end of ! isSubstep
    }

private:
    int dam_err_ ;
    int rank_  ;
    int nranks_ ;
    int numElements_ ;  ///<  size of the unique vector elements
    std::unordered_set<std::string> wanted_vars_set_ ;

    Simulator& simulator_;
    std::unique_ptr<OutputBlackOilModule<TypeTag>> damarisOutputModule_;
    std::vector<unsigned long long> elements_rank_offsets_ ;
    bool damarisUpdate_ = false;  ///< Whenever this is true writeOutput() will set up Damaris mesh information and offsets of model fields

    static bool enableDamarisOutput_()
    {
        static bool enable = Parameters::Get<Parameters::EnableDamarisOutput>();
        return enable;
    }

    void setGlobalIndexForDamaris ()
    {
        const auto& cc = simulator_.vanguard().grid().comm();
        // Use damaris_set_position to set the offset in the global size of the array.
        // This is used so that output functionality (e.g. HDF5Store) knows the global offsets of
        // the data of the ranks data.
        dam_err_ = DamarisOutput::setPosition("GLOBAL_CELL_INDEX", elements_rank_offsets_[rank_]);
        DamarisOutput::handleError(dam_err_, cc, "setPosition() for GLOBAL_CELL_INDEX");

        // This is an example of writing to the Damaris shared memory directly (i.e. we allocate the
        // variable directly in the shared memory region and do not use damaris_write() to copy data there.
        // The shared memory is given back to Damaris when the DamarisVarInt goes out of scope.
        // N.B. MPI_RANK is only saved to HDF5 if --damaris-save-mesh-to-hdf=true is specified
        DamarisVarInt mpi_rank_var(1, {"n_elements_local"}, "MPI_RANK", rank_);
        mpi_rank_var.setDamarisPosition({static_cast<int64_t>(elements_rank_offsets_[rank_])});

        // GLOBAL_CELL_INDEX is used to reorder variable data when writing to disk
        // This is enabled using select-file="GLOBAL_CELL_INDEX" in the <variable> XML tag
        if (this->collectOnIORank_.isParallel()) {
            const std::vector<int>& local_to_global =
                this->collectOnIORank_.localIdxToGlobalIdxMapping();
            dam_err_ = DamarisOutput::write("GLOBAL_CELL_INDEX", local_to_global.data());
        } else {
            std::vector<int> local_to_global_filled ;
            local_to_global_filled.resize(this->numElements_) ;
            std::iota(local_to_global_filled.begin(), local_to_global_filled.end(), 0);
            dam_err_ = DamarisOutput::write("GLOBAL_CELL_INDEX", local_to_global_filled.data());
        }
        DamarisOutput::handleError(dam_err_, cc, "write() for GLOBAL_CELL_INDEX");

        mpi_rank_var.setDamarisParameterAndShmem( std::vector{this->numElements_ } ) ;
        // Fill the created memory area
        std::fill(mpi_rank_var.data(), mpi_rank_var.data() + numElements_, rank_);

        // Pass the output directory string through as a Damaris variable so that
        // Python code (as an example) can use the path as required.
        const auto& outputDir = simulator_.vanguard().eclState().cfg().io().getOutputDir();
        if (outputDir.size() > 0) {
            dam_err_ = DamarisOutput::setParameter("path_string_length", outputDir.size()) ;
            DamarisOutput::handleError(dam_err_, cc, "setParameter() for path_string_length");
            dam_err_ = DamarisOutput::write("OUTPUTDIR", outputDir.c_str());
            DamarisOutput::handleError(dam_err_, cc, "write() for OUTPUTDIR");
        }
    }

    void writeDamarisGridOutput()
    {
        const auto& gridView = simulator_.gridView();
        GridDataOutput::SimMeshDataAccessor geomData(gridView, Dune::Partitions::interior) ;

        try {
            const bool hasPolyCells = geomData.polyhedralCellPresent() ;
            if ( hasPolyCells ) {
                OpmLog::error(fmt::format(fmt::runtime("ERORR: rank {} The DUNE geometry grid has polyhedral elements - These elements are currently not supported."), rank_ ));
            }

            // This is the template XML model for x,y,z coordinates defined in initDamarisXmlFile.cpp which is used to
            // build the internally generated Damaris XML configuration file.
            // <parameter name="n_coords_local"     type="int" value="1" />
            // <parameter name="n_coords_global"    type="int" value="1" comment="only needed if we need to write to HDF5 in Collective mode"/>
            // <layout    name="n_coords_layout"    type="double" dimensions="n_coords_local"   comment="For the individual x, y and z coordinates of the mesh vertices"  />
            // <group name="coordset/coords/values">
            //     <variable name="x"    layout="n_coords_layout"  type="scalar"  visualizable="false"  unit="m"   script="PythonConduitTest" time-varying="false" />
            //     <variable name="y"    layout="n_coords_layout"  type="scalar"  visualizable="false"  unit="m"   script="PythonConduitTest" time-varying="false" />
            //     <variable name="z"    layout="n_coords_layout"  type="scalar"  visualizable="false"  unit="m"   script="PythonConduitTest" time-varying="false" />
            // </group>

            DamarisVarDbl  var_x(1, {"n_coords_local"}, "coordset/coords/values/x", rank_) ;
            // N.B. We have not set any position/offset values (using DamarisVar::SetDamarisPosition).
            // They are not needed for mesh data as each process has a local geometric model.
            // However, HDF5 collective and Dask arrays cannot be used for this data.
            var_x.setDamarisParameterAndShmem( std::vector{ geomData.getNVertices() } ) ;

            DamarisVarDbl var_y(1, {"n_coords_local"}, "coordset/coords/values/y", rank_) ;
            var_y.setDamarisParameterAndShmem( std::vector{ geomData.getNVertices() } ) ;

            DamarisVarDbl  var_z(1, {"n_coords_local"}, "coordset/coords/values/z", rank_) ;
            var_z.setDamarisParameterAndShmem( std::vector{ geomData.getNVertices() } ) ;

            // Now we can use the shared memory area that Damaris has allocated and use it to write the x,y,z coordinates
            if ( geomData.writeGridPoints(var_x, var_y, var_z) < 0)
                 DUNE_THROW(Dune::IOError, geomData.getError()  );

            //  This is the template XML model for connectivity, offsets and types, as defined in initDamarisXmlFile.cpp which is used to
            //  build the internally generated Damaris XML configuration file.
            // <parameter name="n_connectivity_ph"        type="int"  value="1" />
            // <layout    name="n_connections_layout_ph"  type="int"  dimensions="n_connectivity_ph"   comment="Layout for connectivities "  />
            // <parameter name="n_offsets_types_ph"       type="int"  value="1" />
            // <layout    name="n_offsets_layout_ph"      type="int"  dimensions="n_offsets_types_ph+1"  comment="Layout for the offsets_ph"  />
            // <layout    name="n_types_layout_ph"        type="char" dimensions="n_offsets_types_ph"  comment="Layout for the types_ph "  />
            // <group name="topologies/topo/elements">
            //     <variable name="connectivity" layout="n_connections_layout_ph"  type="scalar"  visualizable="false"  unit=""   script="PythonConduitTest" time-varying="false" />
            //     <variable name="offsets"      layout="n_offsets_layout_ph"    type="scalar"  visualizable="false"  unit=""   script="PythonConduitTest" time-varying="false" />
            //     <variable name="types"        layout="n_types_layout_ph"    type="scalar"  visualizable="false"  unit=""   script="PythonConduitTest" time-varying="false" />
            // </group>

            DamarisVarInt var_connectivity(1, {"n_connectivity_ph"},
                                           "topologies/topo/elements/connectivity", rank_) ;
            var_connectivity.setDamarisParameterAndShmem(std::vector{ geomData.getNCorners()}) ;
            DamarisVarInt  var_offsets(1, {"n_offsets_types_ph"},
                                      "topologies/topo/elements/offsets", rank_) ;
            var_offsets.setDamarisParameterAndShmem(std::vector{ geomData.getNCells()+1}) ;
            DamarisVarChar  var_types(1, {"n_offsets_types_ph"},
                                     "topologies/topo/elements/types", rank_) ;
            var_types.setDamarisParameterAndShmem(std::vector{ geomData.getNCells()}) ;

            // Copy the mesh data from the Dune grid
            long i = 0 ;
            GridDataOutput::ConnectivityVertexOrder vtkorder = GridDataOutput::VTK ;

            i = geomData.writeConnectivity(var_connectivity, vtkorder) ;
            if ( i  != geomData.getNCorners())
                 DUNE_THROW(Dune::IOError, geomData.getError());

            i = geomData.writeOffsetsCells(var_offsets);
            if ( i != geomData.getNCells()+1)
                 DUNE_THROW(Dune::IOError,geomData.getError());

            i = geomData.writeCellTypes(var_types) ;
            if ( i != geomData.getNCells())
                 DUNE_THROW(Dune::IOError,geomData.getError());
        }
        catch (std::exception& e)
        {
            OpmLog::error(e.what());
        }
    }

    void prepareLocalCellData(const bool isSubStep,
                              const int  reportStepNum)
    {
        OPM_TIMEBLOCK(prepareLocalCellData);
        if (damarisOutputModule_->localDataValid()) {
            return;
        }

        const auto& gridView = simulator_.vanguard().gridView();
        const int num_interior = detail::
            countLocalInteriorCellsGridView(gridView);
        const bool log = this->collectOnIORank_.isIORank();

        damarisOutputModule_->allocBuffers(num_interior, reportStepNum,
                                      isSubStep, log, /*isRestart*/ false);

        ElementContext elemCtx(simulator_);
        OPM_BEGIN_PARALLEL_TRY_CATCH();
        {
        OPM_TIMEBLOCK(prepareCellBasedData);
        damarisOutputModule_->setupExtractors(isSubStep, reportStepNum);
        for (const auto& elem : elements(gridView, Dune::Partitions::interior)) {
            elemCtx.updatePrimaryStencil(elem);
            elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

            damarisOutputModule_->processElement(elemCtx);
        }
        damarisOutputModule_->clearExtractors();
        }
        {
        OPM_TIMEBLOCK(prepareBlockData);
        for (const auto& elem : elements(gridView, Dune::Partitions::interior)) {
            elemCtx.updatePrimaryStencil(elem);
            elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
            damarisOutputModule_->processElementBlockData(elemCtx);
        }
        }
        {
        OPM_TIMEBLOCK(prepareFluidInPlace);
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int dofIdx=0; dofIdx < num_interior; ++dofIdx){
                const auto& intQuants = *(simulator_.model().cachedIntensiveQuantities(dofIdx, /*timeIdx=*/0));
                const auto totVolume = simulator_.model().dofTotalVolume(dofIdx);
                damarisOutputModule_->updateFluidInPlace(dofIdx, intQuants, totVolume);
        }
        }
        damarisOutputModule_->validateLocalData();
        OPM_END_PARALLEL_TRY_CATCH("DamarisWriter::prepareLocalCellData() failed: ", simulator_.vanguard().grid().comm());
    }

};

} // namespace Opm

#endif // OPM_DAMARIS_WRITER_HPP
