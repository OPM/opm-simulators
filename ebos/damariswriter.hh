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

#include <ebos/collecttoiorank.hh>
#include <ebos/eclbasevanguard.hh>
#include <ebos/eclgenericwriter.hh>
#include <ebos/ecloutputblackoilmodule.hh>

#include <opm/input/eclipse/Units/UnitSystem.hpp>

#include <opm/output/eclipse/RestartValue.hpp>

#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/utils/ParallelRestart.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <limits>
#include <stdexcept>
#include <string>

#include <omp.h>

#include <fmt/format.h>

#include <Damaris.h>
#include <opm/simulators/utils/GridDataOutput.hpp>
#include <opm/simulators/utils/DamarisVar.hpp>

namespace Opm::Properties {

template<class TypeTag, class MyTypeTag>
struct EnableDamarisOutput {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct EnableDamarisOutputCollective {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct DamarisSaveToHdf {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct DamarisPythonScript {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct DamarisPythonParaviewScript {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct DamarisSimName {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct DamarisDedicatedCores {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct DamarisDedicatedNodes {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct DamarisSharedMemorySizeBytes {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct DamarisLogLevel {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct DamarisDaskFile {
    using type = UndefinedProperty;
};
} // namespace Opm::Properties

namespace Opm {

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
    using DamarisVarInt = Opm::DamarisOutput::DamarisVar<int> ;
    using DamarisVarChar = Opm::DamarisOutput::DamarisVar<char> ;
    using DamarisVarDbl = Opm::DamarisOutput::DamarisVar<double>  ;

public:
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, bool, EnableDamarisOutputCollective,
                             "Write output via Damaris using parallel HDF5 to get single file per timestep instead of one per Damaris core.");
        EWOMS_REGISTER_PARAM(TypeTag, bool, DamarisSaveToHdf,
                             "Set to false to prevent output to HDF5. Uses collective output by default or set --enable-damaris-collective=false to\n \
                                                  use file per core (file per Damaris server).");
        EWOMS_REGISTER_PARAM(TypeTag, std::string, DamarisPythonScript,
                             "Set to the path and filename of a Python script to run on Damaris server resources with access to OPM flow data.");
        EWOMS_REGISTER_PARAM(TypeTag, std::string, DamarisPythonParaviewScript,
                             "Set to the path and filename of a Paraview Python script to run on Paraview Catalyst (1 or 2) on Damaris server \n \
                                                 resources with access to OPM flow data.");
        EWOMS_REGISTER_PARAM(TypeTag, std::string, DamarisSimName,
                             "The name of the simulation to be used by Damaris. If empty (the default) then Damaris uses \"opm-sim-<magic_number>\". \n \
                                                 This name should preferably be unique as it is used for the Damaris shmem name and by the Python Dask \n \
                                                 library to locate sections of variables.");
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

    ~DamarisWriter()
    { }

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
                this->setupDamarisWritingPars(simulator_.vanguard().grid().comm(), numElements_, elements_rank_offsets_);
                
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
                int64_t temp_int64_t[1];
                temp_int64_t[0] = static_cast<int64_t>(this->elements_rank_offsets_[rank_]);
                dam_err_ = damaris_set_position("PRESSURE", temp_int64_t);
                if (dam_err_ != DAMARIS_OK && rank_ == 0) {
                    OpmLog::error(fmt::format("damariswriter::writeOutput()       : ( rank:{})"
                                              "damaris_set_position(PRESSURE, ...), Damaris Error: {}  ",
                                              rank_, damaris_error_string(dam_err_) ));
                }

                dam_err_ = damaris_write("PRESSURE", (void*)this->damarisOutputModule_->getPRESSURE_ptr());
                if (dam_err_ != DAMARIS_OK) {
                   OpmLog::error(fmt::format("damariswriter::writeOutput()       : ( rank:{}) "
                                             "damaris_write(PRESSURE, ...), Damaris Error: {}  ",
                                             rank_, damaris_error_string(dam_err_) ));
                }

                dam_err_ =  damaris_end_iteration();
                if (dam_err_ != DAMARIS_OK) {
                    OpmLog::error(fmt::format("damariswriter::writeOutput()       : ( rank:{}) "
                                              "damaris_end_iteration(), Damaris Error: {}  ",
                                              rank_, damaris_error_string(dam_err_) ));
                }
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
        if ( this->collectToIORank_.isParallel() ){
            const std::vector<int>& local_to_global =  this->collectToIORank_.localIdxToGlobalIdxMapping(); 
            dam_err_ = damaris_write("GLOBAL_CELL_INDEX", local_to_global.data());
        } else {
            std::vector<int> local_to_global_filled ;
            local_to_global_filled.resize(this->numElements_) ;
            for (int i = 0 ; i < this->numElements_ ; i++)
            {
                local_to_global_filled[i] = i ;
            }
            dam_err_ = damaris_write("GLOBAL_CELL_INDEX", local_to_global_filled.data());
        }

        if (dam_err_ != DAMARIS_OK) {
            OpmLog::error(fmt::format("damariswriter::writeOutput()       :"
                                      "( rank:{}) damaris_write(GLOBAL_CELL_INDEX, ...), Damaris Error: {}  ",  
                                      rank_, damaris_error_string(dam_err_) ));
        }

        // This is an example of writing to the Damaris shared memory directly (i.e. not using damaris_write() to copy data there)
        // We will add the MPI rank value directly into shared memory using the DamarisVar wrapper of the C based Damaris API
        // The shared memory is given back to Damaris on object deletion - i.e. when the unique_ptr goes out of scope.
        //auto mpi_rank_var = std::make_unique<Opm::DamarisOutput::DamarisVar<int>>(
        //    1, {std::string("n_elements_local")}, std::string("MPI_RANK"), rank_)) ; 
        // std::unique_ptr<Opm::DamarisOutput::DamarisVar<int>>  
        std::unique_ptr<DamarisVarInt> mpi_rank_var( new DamarisVarInt(1, 
                                                                       {std::string("n_elements_local")}, 
                                                                       std::string("MPI_RANK"), rank_) ) ;
        // N.B. we have not set any offset values, so HDF5 collective and Dask arrays cannot be used.
        mpi_rank_var->setDamarisParameterAndShmem( {this->numElements_ } ) ;
        // Fill the created memory area
        for (int i = 0 ; i < this->numElements_; i++ )
        {
            mpi_rank_var->data()[i] = rank_ ;  // write the rank vaue to the shared memory area.
        }
    }

    void setupDamarisWritingPars(Parallel::Communication comm, const int n_elements_local_grid, std::vector<unsigned long long>& elements_rank_offsets)
    {
        // one for each rank -- to be gathered from each client rank
        std::vector<unsigned long long> elements_rank_sizes(nranks_); 
        // n_elements_local_grid should be the full model size
        const unsigned long long n_elements_local = n_elements_local_grid;

        // This gets the n_elements_local from all ranks and copies them to a std::vector of all the values on all ranks
        // (elements_rank_sizes[]).
        comm.allgather(&n_elements_local, 1, elements_rank_sizes.data());
        elements_rank_offsets[0] = 0ULL;
        // This scan makes the offsets to the start of each ranks grid section if each local grid data was concatenated (in
        // rank order)
        for (int t1 = 1; t1 < nranks_; t1++) {
            elements_rank_offsets[t1] = elements_rank_offsets[t1 - 1] + elements_rank_sizes[t1 - 1];
        }

        // find the global/total size
        unsigned long long n_elements_global_max = elements_rank_offsets[nranks_ - 1];
        n_elements_global_max += elements_rank_sizes[nranks_ - 1]; // add the last ranks size to the already accumulated offset values

        if (rank_ == 0) {
            OpmLog::debug(fmt::format("In setupDamarisWritingPars(): n_elements_global_max = {}", n_elements_global_max));
        }

        // Set the paramater so that the Damaris servers can allocate the correct amount of memory for the variabe
        // Damaris parameters only support int data types. This will limit models to be under size of 2^32-1 elements
        // ToDo: Do we need to check that local ranks are 0 based ?
        int temp_int = static_cast<int>(elements_rank_sizes[rank_]);
        dam_err_ = damaris_parameter_set("n_elements_local", &temp_int, sizeof(int));
        if (dam_err_ != DAMARIS_OK) {
            OpmLog::error("( rank:" + std::to_string(rank_)+") Damaris library produced an error result for "
                          "damaris_parameter_set(\"n_elements_local\", &temp_int, sizeof(int));");
        }
        // Damaris parameters only support int data types. This will limit models to be under size of 2^32-1 elements
        // ToDo: Do we need to check that n_elements_global_max will fit in a C int type (INT_MAX)
        if( n_elements_global_max <= std::numeric_limits<int>::max() ) {
            temp_int = static_cast<int>(n_elements_global_max);
            dam_err_ = damaris_parameter_set("n_elements_total", &temp_int, sizeof(int));
            if (dam_err_ != DAMARIS_OK) {
                OpmLog::error("( rank:" + std::to_string(rank_)+") Damaris library produced an error result for "
                              "damaris_parameter_set(\"n_elements_total\", &temp_int, sizeof(int));");
            }
        } else {
            OpmLog::error(fmt::format("( rank:{} ) The size of the global array ({}) is"
                                      "greater than what a Damaris paramater type supports ({}).  ", 
                                      rank_, n_elements_global_max, std::numeric_limits<int>::max() ));
            // assert( n_elements_global_max <= std::numeric_limits<int>::max() ) ;
            OPM_THROW(std::runtime_error, "setupDamarisWritingPars() n_elements_global_max > std::numeric_limits<int>::max() " + std::to_string(dam_err_));
        }

        // Use damaris_set_position to set the offset in the global size of the array.
        // This is used so that output functionality (e.g. HDF5Store) knows global offsets of the data of the ranks
        int64_t temp_int64_t[1];
        temp_int64_t[0] = static_cast<int64_t>(elements_rank_offsets[rank_]);
        dam_err_ = damaris_set_position("PRESSURE", temp_int64_t);
        if (dam_err_ != DAMARIS_OK) {
            OpmLog::error("( rank:" + std::to_string(rank_)+") Damaris library produced an error result for "
                          "damaris_set_position(\"PRESSURE\", temp_int64_t);");
        }
        dam_err_ = damaris_set_position("GLOBAL_CELL_INDEX", temp_int64_t);
        if (dam_err_ != DAMARIS_OK) {
            OpmLog::error("( rank:" + std::to_string(rank_)+") Damaris library produced an error result for "
                          "damaris_set_position(\"GLOBAL_CELL_INDEX\", temp_int64_t);");
        }

        //auto mpi_rank_var = std::make_unique<Opm::DamarisOutput::DamarisVar<int>>(
        //    1, {std::string("n_elements_local")}, std::string("MPI_RANK"), rank_)) ; 
            
        // std::unique_ptr<Opm::DamarisOutput::DamarisVar<int>>  
        // mpi_rank_var(new Opm::DamarisOutput::DamarisVar<int>(1, {std::string("n_elements_local")}, std::string("MPI_RANK"), rank_)) ;
        std::unique_ptr<DamarisVarInt> mpi_rank_var( new DamarisVarInt(1, 
                                                                       {std::string("n_elements_local")}, 
                                                                       std::string("MPI_RANK"), rank_) ) ;
        mpi_rank_var->setDamarisPosition({*temp_int64_t}) ;

    }
    


    void writeDamarisGridOutput( void )
    {
        const auto& gridView = simulator_.gridView();
        Opm::GridDataOutput::SimMeshDataAccessor geomData(gridView, Dune::Partitions::interior) ;

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

            std::unique_ptr<Opm::DamarisOutput::DamarisVar<double>>  var_x(new Opm::DamarisOutput::DamarisVar<double>(1, {std::string("n_coords_local")}, std::string("coordset/coords/values/x"), rank_)) ; 
            // N.B. We have not set any position/offset values (using DamarisVar::SetDamarisPosition). 
            // They are not needed for mesh data as each process has a local geometric model. 
            // However, HDF5 collective and Dask arrays cannot be used for this data.
            var_x->setDamarisParameterAndShmem( { geomData.getNVertices() } ) ;
             
            std::unique_ptr<Opm::DamarisOutput::DamarisVar<double>>  var_y(new Opm::DamarisOutput::DamarisVar<double>(1, {std::string("n_coords_local")}, std::string("coordset/coords/values/y"), rank_)) ; 
            var_y->parameterIsSet() ;
            var_y->setPointersToDamarisShmem() ;
             
            std::unique_ptr<Opm::DamarisOutput::DamarisVar<double>>  var_z(new Opm::DamarisOutput::DamarisVar<double>(1, {std::string("n_coords_local")}, std::string("coordset/coords/values/z"), rank_)) ; 
            var_z->parameterIsSet() ;
            var_z->setPointersToDamarisShmem() ;
            
            // Now we can return the memory that Damaris has allocated in shmem and use it to write the X,y,z coordinates
            double itime, ftime, exec_time;
            itime = omp_get_wtime();
            if ( geomData.writeGridPoints(*var_x,*var_y,*var_z) < 0)
                 DUNE_THROW(Dune::IOError, geomData.getError()  );
             
            //if ( geomData.writeGridPoints(var_x->data(),var_y->data(),var_z->data()) < 0)
            //     DUNE_THROW(Dune::IOError, geomData.getError()  );

            ftime = omp_get_wtime();
            exec_time = ftime - itime;
            // OpmLog::info("\n\nTime taken geomData.writeGridPoints(): is  " + std::to_string(exec_time) ) ;
            std::cout << "\n\n rank_: " << rank_ << " Time taken geomData.writeGridPoints(): is  " + std::to_string(exec_time)  << std::endl ;

            
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

            std::unique_ptr<DamarisVarInt>  var_connectivity(new DamarisVarInt(1, {std::string("n_connectivity_ph")}, std::string("topologies/topo/elements/connectivity"), rank_)) ;
            var_connectivity->setDamarisParameterAndShmem({ geomData.getNCorners()}) ;
            std::unique_ptr<DamarisVarInt>  var_offsets(new DamarisVarInt(1, {std::string("n_offsets_types_ph")}, std::string("topologies/topo/elements/offsets"), rank_)) ;
            var_offsets->setDamarisParameterAndShmem({ geomData.getNCells()}) ;
            std::unique_ptr<DamarisVarChar>  var_types(new DamarisVarChar(1, {std::string("n_offsets_types_ph")}, std::string("topologies/topo/elements/types"), rank_)) ;
            var_types->parameterIsSet() ;
            var_types->setPointersToDamarisShmem() ;

            // Copy the mesh data from the Durne grid
            long i = 0 ;
            Opm::GridDataOutput::ConnectivityVertexOrder vtkorder = Opm::GridDataOutput::VTK ;
            
            i = geomData.writeConnectivity(*var_connectivity, vtkorder) ;
            if ( i  != geomData.getNCorners())
                 DUNE_THROW(Dune::IOError, geomData.getError());

            i = geomData.writeOffsetsCells(*var_offsets);
            if ( i != geomData.getNCells()+1)
                 DUNE_THROW(Dune::IOError,geomData.getError());

            i = geomData.writeCellTypes(*var_types) ;
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
        const int numElements = gridView.size(/*codim=*/0);
        const bool log = this->collectToIORank_.isIORank();

        damarisOutputModule_->allocBuffers(numElements, reportStepNum,
                                      isSubStep, log, /*isRestart*/ false);

        ElementContext elemCtx(simulator_);
        OPM_BEGIN_PARALLEL_TRY_CATCH();
        {
        OPM_TIMEBLOCK(prepareCellBasedData);
        for (const auto& elem : elements(gridView)) {
            elemCtx.updatePrimaryStencil(elem);
            elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

            damarisOutputModule_->processElement(elemCtx);
        }
        }
        if(!simulator_.model().linearizer().getFlowsInfo().empty()){
            OPM_TIMEBLOCK(prepareFlowsData);
            for (const auto& elem : elements(gridView)) {
                elemCtx.updatePrimaryStencil(elem);
                elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
                damarisOutputModule_->processElementFlows(elemCtx);
            }
        }
        {
        OPM_TIMEBLOCK(prepareBlockData);
        for (const auto& elem : elements(gridView)) {
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
        for (int dofIdx=0; dofIdx < numElements; ++dofIdx){
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
