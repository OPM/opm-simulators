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
#ifndef EWOMS_DAMARIS_WRITER_HH
#define EWOMS_DAMARIS_WRITER_HH

#include <dune/grid/common/partitionset.hh>

#include <ebos/damaris_properties.hh>
#include <ebos/eclbasevanguard.hh>
#include <ebos/eclgenericwriter.hh>
#include <ebos/ecloutputblackoilmodule.hh>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/simulators/flow/countGlobalCells.hpp>
#include <opm/simulators/utils/DamarisVar.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/utils/GridDataOutput.hpp>

#include <fmt/format.h>

#include <algorithm>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

namespace Opm {

namespace DamarisOutput {

int endIteration(int rank);
int setParameter(const char* field, int rank, int value);
int setPosition(const char* field, int rank, int64_t pos);
int write(const char* field, int rank, const void* data);
int setupWritingPars(Parallel::Communication comm,
                     const int n_elements_local_grid,
                     std::vector<unsigned long long>& elements_rank_offsets);
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
        EWOMS_REGISTER_PARAM(TypeTag, bool, DamarisOutputHdfCollective,
                             "Write output via Damaris using parallel HDF5 to get single file and dataset per timestep instead of one per Damaris \n \
                                                   core with multiple datasets.");
        EWOMS_REGISTER_PARAM(TypeTag, bool, DamarisSaveToHdf,
                             "Set to false to prevent output to HDF5. Uses collective output by default or set --enable-damaris-collective=false to\n \
                                                   use file per core (file per Damaris server).");
        EWOMS_REGISTER_PARAM(TypeTag, bool, DamarisSaveMeshToHdf,
                             "Saves the mesh data to the HDF5 file (1st iteration only). Will set  --damaris-output-hdf-collective to false \n \
                                                   so will use file per core (file per Damaris server) output (global sizes and offset values \n \
                                                   of mesh variables are not being provided as yet).");
        EWOMS_REGISTER_PARAM(TypeTag, std::string, DamarisPythonScript,
                             "Set to the path and filename of a Python script to run on Damaris server resources with access to OPM flow data.");
        EWOMS_REGISTER_PARAM(TypeTag, std::string, DamarisPythonParaviewScript,
                             "Set to the path and filename of a Paraview Python script to run on Paraview Catalyst (1 or 2) on Damaris server \n \
                                                  resources with access to OPM flow data.");
        EWOMS_REGISTER_PARAM(TypeTag, std::string, DamarisSimName,
                             "The name of the simulation to be used by Damaris. If empty (the default) then Damaris uses \"opm-sim-<random-number>\". \n \
                                                  This name is used for the Damaris HDF5 file name prefix. Make unique if writing to the same output directory.");
        EWOMS_REGISTER_PARAM(TypeTag, std::string, DamarisLogLevel,
                             "The log level for the Damaris logging system (boost log based). \n \
                                                  Levels are: [trace, debug, info, warning, error, fatal]. Currently debug and info are useful. ");
        EWOMS_REGISTER_PARAM(TypeTag, std::string, DamarisDaskFile,
                             "The name of a Dask json configuration file (if using Dask for processing).");                                         
                                                 
        EWOMS_REGISTER_PARAM(TypeTag, int, DamarisDedicatedCores,
                             "Set the number of dedicated cores (MPI processes) that should be used for Damaris processing (per node). \n \
                                                  Must divide evenly into the number of simulation ranks (client ranks).");
        EWOMS_REGISTER_PARAM(TypeTag, int, DamarisDedicatedNodes,
                             "Set the number of dedicated nodes (full nodes) that should be used for Damaris processing (per simulation). \n \
                                                  Must divide evenly into the number of simulation nodes.");
        EWOMS_REGISTER_PARAM(TypeTag, long, DamarisSharedMemorySizeBytes,
                             "Set the size of the shared memory buffer used for IPC between the simulation and the Damaris resources. \n \
                                                  Needs to hold all the variables published, possibly over multiple simulation iterations.");
        EWOMS_REGISTER_PARAM(TypeTag, std::string, DamarisSharedMemoryName,
                             "The name of the shared memory area to be used by Damaris for the current. If empty (the default) then Damaris uses \"opm-damaris-<random-string>\". \n \
                                                  This name should be unique if multiple simulations are running on the same node/server as it is used for the Damaris shmem name and by the Python Dask \n \
                                                  library to locate sections of variables.");
    }

    // The Simulator object should preferably have been const - the
    // only reason that is not the case is due to the SummaryState
    // object owned deep down by the vanguard.
    DamarisWriter(Simulator& simulator)
        : BaseType(simulator.vanguard().schedule(),
                   simulator.vanguard().eclState(),
                   simulator.vanguard().summaryConfig(),
                   simulator.vanguard().grid(),
                   simulator.vanguard().grid().comm().rank() == 0 ? &simulator.vanguard().equilGrid() : nullptr,
                   simulator.vanguard().gridView(),
                   simulator.vanguard().cartesianIndexMapper(),
                   simulator.vanguard().grid().comm().rank() == 0 ? &simulator.vanguard().equilCartesianIndexMapper() : nullptr,
                   false, false)
        , simulator_(simulator)
    {
        this->damarisUpdate_ = true ;

        rank_ = simulator_.vanguard().grid().comm().rank() ;
        nranks_ = simulator_.vanguard().grid().comm().size();
        
        const auto& gridView = simulator_.gridView();
        const auto& interior_elements = elements(gridView, Dune::Partitions::interior);
        // Get the size of the unique vector elements (excludes the shared 'ghost' elements)
        numElements_ = std::distance(interior_elements.begin(), interior_elements.end());

        this->elements_rank_offsets_.resize(nranks_) ;
        this->damarisOutputModule_ = std::make_unique<EclOutputBlackOilModule<TypeTag>>(simulator, this->collectToIORank_);
    }

    /*!
     * \brief Writes localCellData through to Damaris servers. Sets up the unstructured mesh which is passed to Damaris.
     */
    void writeOutput(data::Solution& localCellData , bool isSubStep)
    {
        OPM_TIMEBLOCK(writeOutput);
        const int reportStepNum = simulator_.episodeIndex() + 1;

        // added this as localCellData was not being written
        if (!isSubStep)
            this->damarisOutputModule_->invalidateLocalData() ;  
        this->prepareLocalCellData(isSubStep, reportStepNum);
        this->damarisOutputModule_->outputErrorLog(simulator_.gridView().comm());

        // The damarisWriter is not outputing well or aquifer data (yet)
        auto localWellData = simulator_.problem().wellModel().wellData(); // data::Well

        if (! isSubStep) 
        {
            if (localCellData.size() == 0) {
                this->damarisOutputModule_->assignToSolution(localCellData);
            }

            // add cell data to perforations for Rft output
            this->damarisOutputModule_->addRftDataToWells(localWellData, reportStepNum);
            
            // On first call and if the mesh and variable size change then set damarisUpdate_ to true
            if (damarisUpdate_ == true) {
                // Sets the damaris parameter values "n_elements_local" and "n_elements_total" 
                // which define sizes of the Damaris variables, per-rank and globally (over all ranks).
                // Also sets the offsets to where a ranks array data sits within the global array. 
                // This is usefull for HDF5 output and for defining distributed arrays in Dask.
                dam_err_ = DamarisOutput::setupWritingPars(simulator_.vanguard().grid().comm(),
                                                           numElements_, elements_rank_offsets_);
                
                // sets data for non-time-varying variables MPI_RANK and GLOBAL_CELL_INDEX
                this->setGlobalIndexForDamaris() ; 
                
                // Set the geometry data for the mesh model.
                // this function writes the mesh data directly to Damaris shared memory using Opm::DamarisOutput::DamarisVar objects.
                this->writeDamarisGridOutput() ;
                
                // Currently by default we assume a static mesh grid (the geometry unchanging through the simulation)
                // Set damarisUpdate_ to true if we want to update the geometry sent to Damaris 
                this->damarisUpdate_ = false; 
            }

            if (this->damarisOutputModule_->getPRESSURE_ptr() != nullptr) 
            {
                dam_err_ = DamarisOutput::setPosition("PRESSURE", rank_,
                                                      this->elements_rank_offsets_[rank_]);
                dam_err_ = DamarisOutput::write("PRESSURE", rank_,
                                                this->damarisOutputModule_->getPRESSURE_ptr());

                dam_err_ =  DamarisOutput::endIteration(rank_);
            }
         } // end of ! isSubstep
    }

private:
    int dam_err_ ;
    int rank_  ;       
    int nranks_ ;
    int numElements_ ;  ///<  size of the unique vector elements
    
    Simulator& simulator_;
    std::unique_ptr<EclOutputBlackOilModule<TypeTag>> damarisOutputModule_;
    std::vector<unsigned long long> elements_rank_offsets_ ;
    bool damarisUpdate_ = false;  ///< Whenever this is true writeOutput() will set up Damaris mesh information and offsets of model fields

    static bool enableDamarisOutput_()
    { 
        return EWOMS_GET_PARAM(TypeTag, bool, EnableDamarisOutput); 
    }

    void setGlobalIndexForDamaris () 
    {
        // GLOBAL_CELL_INDEX is used to reorder variable data when writing to disk 
        // This is enabled using select-file="GLOBAL_CELL_INDEX" in the <variable> XML tag
        if (this->collectToIORank_.isParallel()) {
            const std::vector<int>& local_to_global =
                this->collectToIORank_.localIdxToGlobalIdxMapping();
            dam_err_ = DamarisOutput::write("GLOBAL_CELL_INDEX", rank_, local_to_global.data());
        } else {
            std::vector<int> local_to_global_filled ;
            local_to_global_filled.resize(this->numElements_) ;
            std::iota(local_to_global_filled.begin(), local_to_global_filled.end(), 0);
            dam_err_ = DamarisOutput::write("GLOBAL_CELL_INDEX", rank_, local_to_global_filled.data());
        }

        // This is an example of writing to the Damaris shared memory directly (i.e. not using 
        // damaris_write() to copy data there)
        // We will add the MPI rank value directly into shared memory using the DamarisVar 
        // wrapper of the C based Damaris API.
        // The shared memory is given back to Damaris when the DamarisVarInt goes out of scope.
        DamarisVarInt mpi_rank_var_test(1, {"n_elements_local"},  "MPI_RANK", rank_);
        mpi_rank_var_test.setDamarisParameterAndShmem( {this->numElements_ } ) ;
        // Fill the created memory area
        std::fill(mpi_rank_var_test.data(), mpi_rank_var_test.data() + numElements_, rank_);
    }

    void writeDamarisGridOutput()
    {
        const auto& gridView = simulator_.gridView();
        GridDataOutput::SimMeshDataAccessor geomData(gridView, Dune::Partitions::interior) ;

        try {
            const bool hasPolyCells = geomData.polyhedralCellPresent() ;
            if ( hasPolyCells ) {
                OpmLog::error(fmt::format("ERORR: rank {} The DUNE geometry grid has polyhedral elements - These elements are currently not supported.", rank_ ));
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
            var_x.setDamarisParameterAndShmem( { geomData.getNVertices() } ) ;
             
            DamarisVarDbl var_y(1, {"n_coords_local"}, "coordset/coords/values/y", rank_) ;
            var_y.setDamarisParameterAndShmem( { geomData.getNVertices() } ) ;
             
            DamarisVarDbl  var_z(1, {"n_coords_local"}, "coordset/coords/values/z", rank_) ;
            var_z.setDamarisParameterAndShmem( { geomData.getNVertices() } ) ;
            
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
            var_connectivity.setDamarisParameterAndShmem({ geomData.getNCorners()}) ;
            DamarisVarInt  var_offsets(1, {"n_offsets_types_ph"},
                                      "topologies/topo/elements/offsets", rank_) ;
            var_offsets.setDamarisParameterAndShmem({ geomData.getNCells()+1}) ;
            DamarisVarChar  var_types(1, {"n_offsets_types_ph"},
                                     "topologies/topo/elements/types", rank_) ;
            var_types.setDamarisParameterAndShmem({ geomData.getNCells()}) ;

            // Copy the mesh data from the Durne grid
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
        const bool log = this->collectToIORank_.isIORank();

        damarisOutputModule_->allocBuffers(num_interior, reportStepNum,
                                      isSubStep, log, /*isRestart*/ false);

        ElementContext elemCtx(simulator_);
        OPM_BEGIN_PARALLEL_TRY_CATCH();
        {
        OPM_TIMEBLOCK(prepareCellBasedData);
        for (const auto& elem : elements(gridView, Dune::Partitions::interior)) {
            elemCtx.updatePrimaryStencil(elem);
            elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

            damarisOutputModule_->processElement(elemCtx);
        }
        }
        if(!simulator_.model().linearizer().getFlowsInfo().empty()){
            OPM_TIMEBLOCK(prepareFlowsData);
            for (const auto& elem : elements(gridView, Dune::Partitions::interior)) {
                elemCtx.updatePrimaryStencil(elem);
                elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
                damarisOutputModule_->processElementFlows(elemCtx);
            }
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

#endif
