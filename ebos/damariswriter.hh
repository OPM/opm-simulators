// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
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

#ifdef HAVE_DAMARIS
#include <ebos/collecttoiorank.hh>
#include <ebos/damarisgenericwriter.hh>
#include <ebos/ecloutputblackoilmodule.hh>

#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/utils/ParallelRestart.hpp>

#include <opm/input/eclipse/Units/UnitSystem.hpp>

#include <opm/output/eclipse/RestartValue.hpp>

#include <dune/grid/common/partitionset.hh>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <limits>
#include <stdexcept>
#include <string>


#include <opm/simulators/utils/DamarisOutputModule.hpp>
#include <opm/simulators/utils/GridDataOutput.hpp>
#include <damaris/util/DamarisGeometryData.hpp>


namespace Opm::Properties {

// #ifdef HAVE_DAMARIS
template<class TypeTag, class MyTypeTag>
struct EnableDamarisOutput {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct EnableDamarisOutputCollective {
    using type = UndefinedProperty;
};
// #endif

} // namespace Opm::Properties

namespace Opm {

// namespace Action { class State; }
// class EclipseIO;
// class UDQState;

/*!
 * \ingroup EclBlackOilSimulator
 *
 * \brief Collects necessary output values and pass it to opm-output.
 *
 * Caveats:
 * - The only DUNE grid which is currently supported is Dune::CpGrid
 *   from the OPM module "opm-grid". Using another grid won't
 *   fail at compile time but you will provoke a fatal exception
 * - This class requires to use the black oil model with the element
 *   centered finite volume discretization.
 */
 
 
template <class TypeTag>
class DamarisWriter : public DamarisGenericWriter<GetPropType<TypeTag, Properties::Grid>,
                                          GetPropType<TypeTag, Properties::EquilGrid>,
                                          GetPropType<TypeTag, Properties::GridView>,
                                          GetPropType<TypeTag, Properties::ElementMapper>,
                                          GetPropType<TypeTag, Properties::Scalar>>
{
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Vanguard = GetPropType<TypeTag, Properties::Vanguard>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using EquilGrid = GetPropType<TypeTag, Properties::EquilGrid>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementMapper = GetPropType<TypeTag, Properties::ElementMapper>;
    using ElementIterator = typename GridView::template Codim<0>::Iterator;
    using BaseType = DamarisGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>;
    
    typedef Dune::MultipleCodimMultipleGeomTypeMapper< GridView > VertexMapper;

    enum { enableEnergy = getPropValue<TypeTag, Properties::EnableEnergy>() };
    enum { enableTemperature = getPropValue<TypeTag, Properties::EnableTemperature>() };
    enum { enableSolvent = getPropValue<TypeTag, Properties::EnableSolvent>() };

public:
    static void registerParameters()
    {
        // EclOutputBlackOilModule<TypeTag>::registerParameters();


        EWOMS_REGISTER_PARAM(TypeTag, bool, EnableDamarisOutputCollective,
                             "Write output via Damaris using parallel HDF5 to get single file per timestep instead of one per Damaris core.");


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
                   EWOMS_GET_PARAM(TypeTag, bool, EnableAsyncEclOutput), EWOMS_GET_PARAM(TypeTag, bool, EnableEsmry))
        , simulator_(simulator)
    {
        //    this->damarisUpdate_ = enableDamarisOutput_();
        this->damarisUpdate_ = true ;
        
        rank_ = simulator_.vanguard().grid().comm().rank() ;
        const auto& gridView = simulator_.gridView();
        const auto& interior_elements = elements(gridView, Dune::Partitions::interior);
        // Get the size of the unique vector elements (excludes the shared 'ghost' elements)
        numElements_ = std::distance(interior_elements.begin(), interior_elements.end());

        this->damarisOutputModule_ = std::make_unique<EclOutputBlackOilModule<TypeTag>>(simulator, this->collectToIORank_);
    }

    ~DamarisWriter()
    { }

    const EquilGrid& globalGrid() const
    {
        return simulator_.vanguard().equilGrid();
    }

    /*!
     * \brief collect and pass data and pass it to eclIO writer
     */
    void evalSummaryState(bool isSubStep)
    {
        OPM_TIMEBLOCK(evalSummaryState);
        const int reportStepNum = simulator_.episodeIndex() + 1;
        /*
          The summary data is not evaluated for timestep 0, that is
          implemented with a:

             if (time_step == 0)
                 return;

          check somewhere in the summary code. When the summary code was
          split in separate methods Summary::eval() and
          Summary::add_timestep() it was necessary to pull this test out
          here to ensure that the well and group related keywords in the
          restart file, like XWEL and XGRP were "correct" also in the
          initial report step.

          "Correct" in this context means unchanged behavior, might very
          well be more correct to actually remove this if test.
        */
        if (reportStepNum == 0)
            return;

        const Scalar curTime = simulator_.time() + simulator_.timeStepSize();
        const Scalar totalCpuTime =
            simulator_.executionTimer().realTimeElapsed() +
            simulator_.setupTimer().realTimeElapsed() +
            simulator_.vanguard().setupTime();

        const auto localWellData            = simulator_.problem().wellModel().wellData();
        const auto localGroupAndNetworkData = simulator_.problem().wellModel()
            .groupAndNetworkData(reportStepNum);

        const auto localAquiferData = simulator_.problem().aquiferModel().aquiferData();
        const auto localWellTestState = simulator_.problem().wellModel().wellTestState();
        this->prepareLocalCellData(isSubStep, reportStepNum);

        if (this->damarisOutputModule_->needInterfaceFluxes(isSubStep)) {
            this->captureLocalFluxData();
        }
         /*
            For ECL ouput this section used: 
            this->collectToIORank_.collect({},
            ... some loging of iteration numbers etc.
            this->evalSummary(...)
        */
    }

    void writeOutput(data::Solution& localCellData , bool isSubStep)
    {
        OPM_TIMEBLOCK(writeOutput);
        const int reportStepNum = simulator_.episodeIndex() + 1;
        std::cout << "INFO: damarisWriter_->writeOutput() prepareLocalCellData being called " << std::endl ;
        
        if (!isSubStep)
            this->damarisOutputModule_->invalidateLocalData() ;  // <jcb> added this as localCellData was not being written
        this->prepareLocalCellData(isSubStep, reportStepNum);
        this->damarisOutputModule_->outputErrorLog(simulator_.gridView().comm());

        // output using eclWriter if enabled
       
        //    The damarisWriter is not outputing well or aquifer data (yet)
            auto localWellData = simulator_.problem().wellModel().wellData(); // data::Well
            auto localAquiferData = simulator_.problem().aquiferModel().aquiferData();
            auto localWellTestState = simulator_.problem().wellModel().wellTestState();
            auto flowsn = this->damarisOutputModule_->getFlowsn();
            const bool isFlowsn = this->damarisOutputModule_->hasFlowsn();
            auto floresn = this->damarisOutputModule_->getFloresn();
            const bool isFloresn = this->damarisOutputModule_->hasFloresn();
        
        
        // auto localWellData = simulator_.problem().wellModel().wellData(); // data::Well
        // data::Solution localCellData = {};
        if (! isSubStep) {
            
            if (localCellData.size() == 0) {
                this->damarisOutputModule_->assignToSolution(localCellData);
            }

            // add cell data to perforations for Rft output
            this->damarisOutputModule_->addRftDataToWells(localWellData, reportStepNum);
        } // end of ! isSubstep

        this->writeDamarisGridOutput(isSubStep) ;
        // int rank =simulator_.vanguard().grid().comm().rank() ;
        
        // if (rank == 0) {
        
        for ( auto damVar : localCellData ) {
           // std::map<std::string, data::CellData>
          const std::string name = damVar.first ;
          data::CellData  dataCol = damVar.second ;
          std::cout << "Name of Damaris Varaiable       : (" << rank_ << ")  "  << name << "  Size : "  << dataCol.data.size() <<  std::endl ;  // dataCol.data().size()
          
          /*if (name == "PRESSURE") {
              if (dataCol.data.data() != nullptr) {
                  damaris_write(name.c_str(), (void*) dataCol.data.data() ) ;
              }
                    //damaris_write("PRESSURE", (void*)this->damarisOutputModule_->getPRESSURE_ptr());
          }*/
        }
        
        this->writeDamarisCellDataOutput(localCellData, isSubStep) ;
        
        auto mybloc = damarisOutputModule_->getBlockData() ;
        for ( auto damVar : mybloc ) {
           // std::map<std::string, data::CellData>
          const std::string name = std::get<0>(damVar.first) ;
          const int part = std::get<1>(damVar.first) ;
          double  dataCol = damVar.second ;
          std::cout << "Name of Damaris Block Varaiable : (" << rank_ << ")  "  << name  << "  part : " << part << "  Value : "  << dataCol <<  std::endl ;  // dataCol.data().size()
        }
        
        if (! isSubStep)
        {
            //damaris_end_iteration();
        }
        /*
            For ECL ouput this section used: 
            this->collectToIORank_.collect()
            this->doWriteOutput(...)
        */
    }
    


    void endRestart()
    {}

    const EclOutputBlackOilModule<TypeTag>& eclOutputModule() const
    { return *damarisOutputModule_; }

    EclOutputBlackOilModule<TypeTag>& mutableEclOutputModule() const
    { return *damarisOutputModule_; }

    Scalar restartTimeStepSize() const
    { return restartTimeStepSize_; }

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(*damarisOutputModule_);
    }

private:

    int damaris_err_ ;
    int rank_  ;//  = simulator_.vanguard().grid().comm().rank() ;
    int numElements_ ;
    
    static bool enableDamarisOutput_()
    { 
        return EWOMS_GET_PARAM(TypeTag, bool, EnableDamarisOutput); 
    }


     void writeDamarisCellDataOutput(data::Solution& localCellData, bool isSubStep)
    {
        static int firstcall = 0 ;
        if (!isSubStep) {
            // Output the PRESSURE field
            if (this->damarisOutputModule_->getPRESSURE_ptr() != nullptr) {
                damaris_write("PRESSURE", (void*)this->damarisOutputModule_->getPRESSURE_ptr());
               // damaris_end_iteration();
               
               // Only call this once, and only when the Damaris data is valid
               if (firstcall == 0) {
                   SetGlobalIndexForDamaris() ;
                    firstcall = 1 ;
               }
               
               damaris_err_ =  damaris_end_iteration();
               if (damaris_err_ != DAMARIS_OK) {
                    std::cerr << "ERROR rank =" << rank_ << " : damariswriter::writeDamarisCellDataOutput() : damaris_end_iteration()" 
                    << ", Damaris error = " <<  damaris_error_string(damaris_err_) << std::endl ;
               }
            }
        }
    }
    
     void SetGlobalIndexForDamaris () 
    {
        // GLOBAL_CELL_INDEX is used to reorder variable data when writing to disk 
        // This is enabled using select-file="GLOBAL_CELL_INDEX" in the <variable> XML tag
        if ( this->collectToIORank_.isParallel() ){
            const std::vector<int>& local_to_global =  this->collectToIORank_.localIdxToGlobalIdxMapping(); 
            damaris_err_ = damaris_write("GLOBAL_CELL_INDEX", local_to_global.data());
        } else {
            std::vector<int> local_to_global_filled ;
            local_to_global_filled.resize(this->numElements_) ;
            for (int i = 0 ; i < this->numElements_ ; i++)
            {
                local_to_global_filled[i] = i ;
            }
            damaris_err_ = damaris_write("GLOBAL_CELL_INDEX", local_to_global_filled.data());
        }
        
        if (damaris_err_ != DAMARIS_OK) {
            std::cerr << "ERROR rank =" << rank_ << " : eclwrite::writeOutput() : damaris_write(\"GLOBAL_CELL_INDEX\", local_to_global.data())" 
                  << ", Damaris error = " <<  damaris_error_string(damaris_err_) << std::endl ;
        }
    }
    
    
    void writeDamarisGridOutput(bool isSubStep)
    {
        // const int rank = simulator_.vanguard().grid().comm().rank() ;
        
        // N.B. damarisUpdate_ should be set to true if at any time the model geometry changes
        if (this->damarisUpdate_) {
            const auto& gridView = simulator_.gridView();
            const auto& interior_elements = elements(gridView, Dune::Partitions::interior);
            // Get the size of the unique vector elements (excludes the shared 'ghost' elements)
            const int numElements = std::distance(interior_elements.begin(), interior_elements.end());
            // Sets the damaris parameter values which then defines sizes of the arrays per-rank and globally.
            // Also sets the offsets to where a ranks array data sits within the global array. 
            // This is usefull for HDF5 output and for defining distributed arrays in Dask.
            Opm::DamarisOutput::setupDamarisWritingPars(simulator_.vanguard().grid().comm(), numElements);
            
            Opm::GridDataOutput::SimMeshDataAccessor geomData(gridView, Dune::Partitions::interior) ; // N.B. we cannot reuse the same object using a different partition.
            try {
               
                damaris::model::vertex_data_structure vertex_structure = damaris::model::VERTEX_SEPARATE_X_Y_Z ;  // define this as we know it works with Ascent
                damaris::model::DamarisGeometryData damarisMeshVars(vertex_structure, geomData.getNVertices(), 
                                                            geomData.getNCells(), geomData.getNCorners(), rank_) ;
               
                const bool hasPolyCells = geomData.polyhedralCellPresent() ;
                if ( hasPolyCells ) {
                    std::cout << "The DUNE geometry grid has polyhedral elements - currently not supported by Damaris " << std::endl ;
                } 
                damarisMeshVars.set_hasPolyhedralCells(hasPolyCells) ;
               
                /* This is our template XML model for x,y,z coordinates
                  <parameter name="n_coords_local"     type="int" value="1" />
                  <parameter name="n_coords_global"    type="int" value="1" comment="only needed if we need to write to HDF5 in Collective mode/>
                  <layout    name="n_coords_layout"    type="double" dimensions="n_coords_local"   comment="For the individual x, y and z coordinates of the mesh vertices, these values are referenced in the topologies/topo/subelements/connectivity_pg data"  />
                  <group name="coordset/coords/values"> 
                      <variable name="x"    layout="n_coords_layout"  type="scalar"  visualizable="false"  unit="m"   script="PythonConduitTest" time-varying="false" />
                      <variable name="y"    layout="n_coords_layout"  type="scalar"  visualizable="false"  unit="m"   script="PythonConduitTest" time-varying="false" />
                      <variable name="z"    layout="n_coords_layout"  type="scalar"  visualizable="false"  unit="m"   script="PythonConduitTest" time-varying="false" />
                  </group>
                */
                int xyz_coord_dims = 1 ;
                std::vector<std::string> param_names = {"n_coords_local"} ;  // a vector of strings as a variables layout may be defined by multiple parameters 
                std::string variable_x = "coordset/coords/values/x" ;  // This string must match the group/variable name of the Damaris XML file 
                std::string variable_y = "coordset/coords/values/y" ;  // This string must match the group/variable name of the Damaris XML file 
                std::string variable_z = "coordset/coords/values/z" ;  // This string must match the group/variable name of the Damaris XML file 
                
                damarisMeshVars.set_damaris_var_name_vertex_x(xyz_coord_dims,  param_names, variable_x) ;
                damarisMeshVars.set_damaris_var_name_vertex_y(xyz_coord_dims,  param_names, variable_y) ;
                damarisMeshVars.set_damaris_var_name_vertex_z(xyz_coord_dims,  param_names, variable_z) ;
                
                // Used to store Damaris parameters, which are then used to resize the 
                // shared memory region used to save the mesh data
                // there should be as many values for paramaters as names used in param_names.
                // Each name corresponds to the value at the same position in the vector.
                std::vector<int> param_vertices ;  
                std::vector<int> param_connectivity ;
                std::vector<int> param_offsets ;
                
                
                param_vertices.push_back(geomData.getNVertices() ) ; 
                // For the vertex data x, y and z arrays, SetAll_VERTEX_SEPARATE_X_Y_Z_shmem()  will set the arrays 
                // size (using the paramater value) and then allocate the shared memory region. This is where we 
                // will write the vertex data to, so that Damaris has access to it.
                damarisMeshVars.SetAll_VERTEX_SEPARATE_X_Y_Z_shmem(param_vertices) ;  
                
                // Now we can return the memory that Damaris has allocated in shmem 
                damaris::model::DamarisVar<double>* var_x =  dynamic_cast<damaris::model::DamarisVar<double>* >(damarisMeshVars.get_x()) ;
                damaris::model::DamarisVar<double>* var_y =  dynamic_cast<damaris::model::DamarisVar<double>* >(damarisMeshVars.get_y()) ;
                damaris::model::DamarisVar<double>* var_z =  dynamic_cast<damaris::model::DamarisVar<double>* >(damarisMeshVars.get_z()) ;
                
                if ( geomData.writeGridPoints(var_x->data_ptr(),var_y->data_ptr(),var_z->data_ptr()) < 0)
                     DUNE_THROW(Dune::IOError, geomData.getError()  );
                
               
                // We do not need these as the ~DamarisGeometryData destructor will call them  
                // damarisMeshVars.CommitAll_VERTEX_SEPARATE_X_Y_Z_shmem() ;
                // damarisMeshVars.ClearAll_VERTEX_SEPARATE_X_Y_Z_shmem() ;
                
                /* This is our template XML model for connectivity 
                  <parameter name="n_connectivity_ph"        type="int"  value="1" />
                  <layout    name="n_connections_layout_ph"  type="int"  dimensions="n_connectivity_ph"   comment="Layout for connectivities "  />
                  <parameter name="n_offsets_types_ph"       type="int"  value="1" />
                  <layout    name="n_offsets_layout_ph"      type="int"  dimensions="n_offsets_types_ph"  comment="Layout for the offsets_ph"  />
                  <layout    name="n_types_layout_ph"        type="char" dimensions="n_offsets_types_ph"  comment="Layout for the types_ph "  />
                  <group name="topologies/topo/elements">
                      <variable name="connectivity" layout="n_connections_layout_ph"  type="scalar"  visualizable="false"  unit=""   script="PythonConduitTest" time-varying="false" />
                      <variable name="offsets"      layout="n_offsets_layout_ph"    type="scalar"  visualizable="false"  unit=""   script="PythonConduitTest" time-varying="false" />
                      <variable name="types"        layout="n_types_layout_ph"    type="scalar"  visualizable="false"  unit=""   script="PythonConduitTest" time-varying="false" />
                  </group>
                */
               
                param_names[0] = "n_connectivity_ph" ; 
                std::string varname = std::string("topologies/topo/elements/connectivity") ;  // This string must match the group/variable name of the Damaris XML file 
                damarisMeshVars.set_damaris_var_name_connectivity(xyz_coord_dims, param_names , varname ) ;
                param_names[0] = "n_offsets_types_ph" ; 
                varname = std::string("topologies/topo/elements/offsets") ;                   // This string must match the group/variable name of the Damaris XML file 
                damarisMeshVars.set_damaris_var_name_offsets(xyz_coord_dims, param_names, varname ) ;
                param_names[0] = "n_offsets_types_ph" ; 
                varname = std::string("topologies/topo/elements/types") ;                     // This string must match the group/variable name of the Damaris XML file 
                damarisMeshVars.set_damaris_var_name_types(xyz_coord_dims, param_names, varname ) ;
                
                // Here we retrieve the DamarisVar objects. We need to *match the type* of the data as set in the Damaris XML file
                damaris::model::DamarisVar<int>* var_connectivity =  dynamic_cast<damaris::model::DamarisVar<int>* >(damarisMeshVars.get_connectivity()) ;
                damaris::model::DamarisVar<int>* var_offsets      =  dynamic_cast<damaris::model::DamarisVar<int>* >(damarisMeshVars.get_offsets()) ;
                damaris::model::DamarisVar<char>* var_types       =  dynamic_cast<damaris::model::DamarisVar<char>* >(damarisMeshVars.get_types()) ;
                
                // Set the Damaris shared memory for each variable needed to store mesh data
                param_connectivity.push_back( geomData.getNCorners() ) ;
                var_connectivity->SetDamarisParameter(param_connectivity) ;
                var_connectivity->SetPointersToDamarisShmem() ;
               
                param_offsets.push_back(geomData.getNCells() ) ;
                var_offsets->SetDamarisParameter(param_offsets) ;
                var_offsets->SetPointersToDamarisShmem() ;
                var_types->SetDamarisParameter(param_offsets) ;
                var_types->SetPointersToDamarisShmem() ;
                
                // Copy the mesh data from the Durne grid
                long i = 0 ;
                Opm::GridDataOutput::ConnectivityVertexOrder vtkorder = Opm::GridDataOutput::VTK ;
                
                i = geomData.writeConnectivity(var_connectivity->data_ptr(), vtkorder) ;
                if ( i  != geomData.getNCorners())
                     DUNE_THROW(Dune::IOError, geomData.getError() );
                
                i = geomData.writeOffsetsCells(var_offsets->data_ptr()) ;
                if ( i != geomData.getNCells()+1)
                     DUNE_THROW(Dune::IOError,geomData.getError() );
                
                i = geomData.writeCellTypes(var_types->data_ptr()) ;
                if ( i != geomData.getNCells())
                     DUNE_THROW(Dune::IOError,geomData.getError() );
                //  Commit and clear damaris functions are called when the object goes out of scope (in the destructor)
            }
            catch (std::exception& e) 
            {
                std :: cout << e.what() << std::endl;
            }


            //this-SetGlobalIndexForDamaris() ;
              
          
            // Currently by default we assume static grid (unchanging through the simulation)
            // Set damarisUpdate_ to true if we want to update the geometry to sent to Damaris 
            this->damarisUpdate_ = false; 
        }// end of if (this->damarisUpdate_)
    }

   

    void prepareLocalCellData(const bool isSubStep,
                              const int  reportStepNum)
    {
        OPM_TIMEBLOCK(prepareLocalCellData);
        if (damarisOutputModule_->localDataValid()) {
            return;
        }
        std::cout << "INFO: damarisWriter_->prepareLocalCellData()  got past localDataValid() " << std::endl ;

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


    void captureLocalFluxData()
    {
        OPM_TIMEBLOCK(captureLocalData);
        const auto& gridView = this->simulator_.vanguard().gridView();
        const auto timeIdx = 0u;

        auto elemCtx = ElementContext { this->simulator_ };

        const auto elemMapper = ElementMapper { gridView, Dune::mcmgElementLayout() };
        const auto activeIndex = [&elemMapper](const Element& e)
        {
            return elemMapper.index(e);
        };

        const auto cartesianIndex = [this](const int elemIndex)
        {
            return this->cartMapper_.cartesianIndex(elemIndex);
        };

        this->damarisOutputModule_->initializeFluxData();

        OPM_BEGIN_PARALLEL_TRY_CATCH();

        for (const auto& elem : elements(gridView, Dune::Partitions::interiorBorder)) {
            elemCtx.updateStencil(elem);
            elemCtx.updateIntensiveQuantities(timeIdx);
            elemCtx.updateExtensiveQuantities(timeIdx);

            this->damarisOutputModule_->processFluxes(elemCtx, activeIndex, cartesianIndex);
        }

        OPM_END_PARALLEL_TRY_CATCH("DamarisWriter::captureLocalFluxData() failed: ",
                                   this->simulator_.vanguard().grid().comm())

        this->damarisOutputModule_->finalizeFluxData();
    }
    

    Simulator& simulator_;
    std::unique_ptr<EclOutputBlackOilModule<TypeTag>> damarisOutputModule_;
    Scalar restartTimeStepSize_;


    bool damarisUpdate_ = false;  ///< Whenever this is true writeOutput() will set up Damaris mesh information and offsets of model fields

};
} // namespace Opm
#endif  // HAVE_DAMARIS

#endif
