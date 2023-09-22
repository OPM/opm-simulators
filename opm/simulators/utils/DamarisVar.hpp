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

#ifndef DAMARIS_MODEL_DATA_HPP
#define DAMARIS_MODEL_DATA_HPP

#include <string>
#include <iostream>
#include <sstream>
#include <cassert>
#include <typeinfo>

#include <Damaris.h>


/*
    File: DamarisVar.hpp
    Author: Joshua Bowden, Inria
    Date: 06/02/2023
    The DamarisVar class can be used to allocate memory in the Damaris shared memory region and a user can supply
    the data required for the variable to pass to Damaris.
*/


namespace Opm
{
    namespace DamarisOutput
    {
    /**
    *  This class contains the extra elements that need to be part of a Damaris <variable> type. 
    *  They are simple string values that may reference other XML elements (and could be checked for existence etc.)
    */
    class DamarisVarXMLAttributes {
            std::string layout_ ;         //!< Reference string to the XML attribute layout being used to describe the shape of the variable. This is a required attribute.
            std::string mesh_ ;           //!< Reference string to the XML attribute mesh element - the mesh is used to define the spatial layout of data and is used by visualization backends to generate 2D/3D model images
            std::string type_ ;           //!< Reference string to the XML attribute type of data - "scalar" or "vector" (others tensor maybe). TODO: check if this attribute is used by the Damaris library anywhere.
            std::string visualizable_ ;   //!< Reference string to the XML attribute property that data can be sent to vis backends - "true" | "false"
            std::string unit_ ;           //!< Reference string to the XML attribute element denoting unit of the data
            std::string time_varying_ ;   //!< Reference string to the XML attribute to indicate if data changes over iterations - "true" | "false"
            std::string centering_ ;      //!< Reference string to the XML attribute to indicate where data aligns on a mesh - "zonal" | "nodal"
            std::string store_ ;          //!< Reference string to the XML attribute to indicate if data should be passed to I/O store (e.g. to HDF5 plugin)
            std::string script_ ;         //!< Reference string to the XML attribute to indicate if data should be published as Python NumPy data
            std::string select_mem_ ;      //!< Reference string to the XML attribute select. The referenced variables data is used as indices to  select dat from memory to reorder output in the collective HDF5 data writer (Damaris version 1.8+)
            std::string select_file_ ;     //!< Reference string to the XML attribute select. The referenced variables data is used as indices to select positions in the dataset file to reorder output in the collective HDF5 data writer (Damaris version 1.8+)
            std::string select_subset_ ;   //!< Reference string to the XML attribute select. Used to specify the output dataset shape and how much data each rank contributes to it and the global offsets to the ranks data (Damaris version 1.8+)

        public:
            DamarisVarXMLAttributes(){
                // Additional data needed to complete an XML <variable> element
                layout_ = "" ;
                mesh_ = "" ;
                type_ = "scalar" ;  // This is probably not needed as vector data is defined using the Layout paramter. Could be useful for cross checking
                visualizable_ = "false";
                unit_ = "" ;
                time_varying_ = "true";
                centering_ = "zonal" ;
                store_ = "" ;
                script_ = "" ;
                select_mem_ = "" ;
            }

            /**
            * Creates the XML representation of the variable from the available strings
            */
            std::string ReturnXMLForVariable ( void )
            {
                std::ostringstream  var_sstr ;

                var_sstr << " layout=\"" << this->layout_ << "\"" ;
                if (this->mesh_ != "") var_sstr  << " mesh=\"" << this->mesh_ << "\"" ;
                if (this->type_ != "") var_sstr  << " type=\"" << this->type_ << "\"" ;
                if (this->visualizable_ != "") var_sstr  << " visualizable=\"" << this->visualizable_ << "\"" ;
                if (this->unit_ != "") var_sstr  << " unit=\"" << this->unit_ << "\"" ;
                if (this->time_varying_ != "") var_sstr  << " time_varying=\"" << this->time_varying_ << "\"" ;
                if (this->centering_ != "") var_sstr  << " centering=\"" << this->centering_ << "\"" ;
                if (this->store_ != "") var_sstr  << " store=\"" << this->store_ << "\"" ;
                if (this->script_ != "") var_sstr  << " script=\"" << this->script_  << "\"" ;
                if (this->select_mem_ != "") var_sstr  << " select-mem=\"" << this->select_mem_ << "\"" ;
                if (this->select_file_ != "") var_sstr  << " select-file=\"" << this->select_file_ << "\"" ;
                if (this->select_subset_ != "") var_sstr  << " select-subset=\"" << this->select_subset_ << "\"" ;

                return (var_sstr.str()) ;
            }
    } ;

    class DamarisVarBase { 
        public:
            // DamarisVarBase(int dims, std::vector<std::string>& param_names, std::string& variable_name, int rank=0)  = 0  ; 
            virtual ~DamarisVarBase( void )  {}  ;

            virtual void PrintError ( void ) = 0  ;
            virtual bool HasError( void ) = 0  ;
            // virtual void SetDamarisParameterAndShmem( std::vector<int>& paramSizeVal )  = 0  ;
            virtual void SetDamarisParameterAndShmem( std::vector<int> paramSizeVal )  = 0  ;
            virtual void SetDamarisParameter( std::vector<int>& paramSizeVal )  = 0  ;
            virtual void SetDamarisPosition(  std::vector<int64_t> positionsVals ) = 0  ;
            virtual void SetPointersToDamarisShmem( void ) = 0  ;
            virtual void CommitVariableDamarisShmem( void ) = 0 ;
            virtual void ClearVariableDamarisShmem( void ) = 0  ;
            // virtual void * data_ptr( void ) = 0  ;
            virtual std::string & variable_name( void ) = 0  ;



    };  // class DamarisVar

    /**
    *  class to store a Damaris variable representation for the XML file (can be used with /ref class DamarisKeywords).
    *
    *  It is thought that the details stored in the object can be used to pass into an XML generation function e.g. DamarisKeywords
    *
    *  Should be instantiated with something like the following:
    *  // The following needs to be defined in the Damaris XML file: 
    *  // <parameter name="n_elements_mpi_local"     type="int" value="1" />
    *  // <layout    name="mpi_layout"    type="int" dimensions="n_elements_mpi_local"   comment="MPI elements layout"  />
    *  // <variable name="MPI_RANK"  layout="mpi_layout"   type="scalar"  visualizable="true" mesh="unstructured_mesh" unit="rank"  centering="zonal"
    *  //                                                   store="#" time-varying="false"  script="_PYTHON_XML_NAME_" comment="The cells MPI rank"/>
    *  damaris::model::DamarisVar<int>  dam_var = new damaris::model::DamarisVar<int>(1, {std::string("n_connectivity_ph")}, std::string("topologies/topo/elements/connectivity"), rank_) ;
    *  dam_var->SetDamarisParameterAndShmem( { geomData.getNCorners() } ) ;
    *
    *  int * shmem_mpi_ptr = dam_var->data_ptr() ;
    *  // Fill the created memory area
    *  for (int i = 0 ; i < ; i++ )
    *  {
    *      shmem_mpi_ptr[i] = rank_ ;
    *  }
    *  delete dam_var ;  // this tells Damaris that the shared memory that it supplied is at its disposal. It will print error messages too.
    *
    */
    template <typename T>
    class DamarisVar :
         public DamarisVarBase {
            int   dims_ ;
            int   num_params_ ;       //!< Each paramater name string will need a value and they are set in SetDamarisParameter()
            int  * param_sizes_ ;     //!< The value for any paramaters that are being used to describe the size of the variables data array
            int64_t * positions_ ;    //!< The offsets into the array that the data in the Variable starts from for this rank.
            int   rank_ ;             //!< Rank of process - used for error reporting.
            bool  paramaters_set_ ;   //!< set to true after SetDamarisParameter() is call to ensure the variable has correct size for memory allocation in SetPointersToDamarisShmem()

            std::vector<std::string> param_names_ ;  //!< Contains one paramater name for each paramater that a variable depends on (via it's Layout)
            std::string variable_name_ ;         //!< Reference string to the XML attribute name of the variable.

            int dam_err_ ;                       //!<  Set to != DAMARIS_OK if a Daamris error was returned by a Damaris API function call
            bool has_error_ ;
            std::ostringstream  dam_err_sstr_ ;  //!< Use dam_err_sstr.str() to return an error string describing detected error

            DamarisVarXMLAttributes xml_attributes_ ; //!< The extra elements that need to be part of a Damaris <variable> type. They are simple string values that may reference other XML elements (and could be checked for existence etc.)

            T * data_ptr_ ;                      //!< This pointer will be mapped to the Damaris shared memory area for the variable in the SetPointersToDamarisShmem() method. The type T will match the Layout type

        public:

            /**
            * Constructor - sets private data values and dos not initialise the shared memory area.
            *
            * Two usages:
            *         Example XML definition:
            *          <parameter name="my_param_name1"     type="int" value="1" />
            *          <parameter name="my_param_name2"     type="int" value="1" />
            *          <layout    name="my_layout"    type="int" dimensions="my_param_name1,my_param_name2"   comment="This is a 2D variable"  />
            *          <variable name="MYVARNAME"  layout="my_layout"  visualizable="true"/>
            *
            *          1/ The variable's layout needs to be initialised via parameters :
            *          // Create the DamarisVar object:
            *          damaris::model::DamarisVar<int>  MYVARNAME_2d = new damaris::model::DamarisVar<int>(2, {std::string("my_param_name1"), std::string("my_param_name2")}, std::string("MYVARNAME"), rank_) ;
            *           // Set the paramter sizes
            *          MYVARNAME_2d->SetDamarisParameterAndShmem( {25, 100 } } ;  // sets the paramaters  (here, my_param_name1 == 25 and my_param_name2 == 100)
            *          // Get a pointer to the memeory and use it
            *          T * mymemory = MYVARNAME_2d->data_ptr() ;
            *          ... write data to mymemory ....
            *          delete MYVARNAME_2d ;
            * or,
            *  2/  The variable's layout has been initialised via parameters in another variable (i.e. "my_param_name1" and "my_param_name2" have been previously set in the code)
            *          // Create the DamarisVar object:
            *          damaris::model::DamarisVar<int>  MYVARNAME_2d = new damaris::model::DamarisVar<int>(2, {std::string("my_param_name1"), std::string("my_param_name2")}, std::string("MYVARNAME"), rank_) ;
            *          // explicitly state that the paramater values have been set somewhere else in the code previously.
            *          MYVARNAME_2d->ParameterIsSet() ;
            *          MYVARNAME_2d->SetPointersToDamarisShmem() 
            *          // Get a pointer to the memeory and use it
            *          T * mymemory = MYVARNAME_2d->data_ptr() ;
            *          ... write data to mymemory ....
            *          delete MYVARNAME_2d ;
            * 
            *  /param [IN] dims           Used to check that the inputs to SetDamarisPosition() have the same number of values - one value for each dimension
            *  /param [IN] param_names    The name the Damaris paramaters. These names (in typical use) control a Damaris variables size (names are defined in the Damaris XML file).
            *  /param [IN] variable_name  The name of the Damaris variable (defined in the Dmaaris XML file)
            *  /param [IN] rank           The rank of the process. Used for error output.
            */
            DamarisVar(int dims, std::vector<std::string> param_names, std::string variable_name, int rank) 
            : dims_( dims ),
              param_names_( param_names ),
              variable_name_( variable_name ),
              rank_(rank)
            {
                dam_err_ = DAMARIS_OK ;
                
                assert( param_names_.size() == dims ) ;
                assert( dims > 0 ) ;
                
                // Check that our template type T matches out Damaris XML <layout> type
                if ( !TestType(variable_name) ) {
                    std::exit(-1) ;
                }

                num_params_ = param_names_.size() ; 
                param_sizes_ = new int(num_params_) ;
                positions_ = new int64_t(dims) ;

                data_ptr_ = nullptr ;
                paramaters_set_ = false ;
                has_error_ = false ;
            }

            /**
            *  Constructor - Sets private data values and also initialises the Damaris shared memory area for writing (and reading) 
            *                by specifying the values for the variables parameters . i.e. makes the data_ptr() 
            * 
            *  Example use:
            *         Example XML definition:
            *          <parameter name="my_param_name1"     type="int" value="1" />
            *          <parameter name="my_param_name2"     type="int" value="1" />
            *          <layout    name="my_layout"    type="int" dimensions="my_param_name1,my_param_name2"   comment="This is a 2D variable"  />
            *          <variable name="MYVARNAME"  layout="my_layout"  visualizable="true"/>
            *  // The paramaters are intialized in the constructor code
            *  damaris::model::DamarisVar<int>  MYVARNAME_2d = new damaris::model::DamarisVar<int>(2, {std::string("my_param_name1"), std::string("my_param_name2")}, {100, 25},  std::string("MYVARNAME"), rank_) ;
            *  T * mymemory = MYVARNAME_2d->data_ptr() ;
            *  ... write data to mymemory ....
            *  delete MYVARNAME_2d ;
            * 
            *  /param [IN] dims           Used to check that the inputs to SetDamarisPosition() have the same number of values - one value for each dimension
            *  /param [IN] param_names    The name the Damaris paramaters. These names (in typical use) control a Damaris variables size (names are defined in the Damaris XML file).
            *  /param [IN] param_values   The values of the paramaters - this defines how much memory we will have access to in the shared memory area (on the current and ongoing iterations, until later modified to new values)
            *  /param [IN] variable_name  The name of the Damaris variable (defined in the Dmaaris XML file)
            *  /param [IN] rank           The rank of the process. Used for error output.
            */
            DamarisVar(int dims, std::vector<std::string> param_names, std::vector<int> param_values, std::string variable_name, int rank) 
            : dims_( dims ),
              param_names_( param_names ),
              variable_name_( variable_name ),
              rank_(rank)
            {
                dam_err_ = DAMARIS_OK ;
                
                assert( param_names_.size() == dims ) ;
                assert( dims > 0 ) ;
                
                // Check that our template type T matches out Damaris XML <layout> type
                if ( !TestType(variable_name) ) {
                    std::exit(-1) ;
                }

                num_params_ = param_names_.size() ; 
                param_sizes_ = new int(num_params_) ;
                positions_ = new int64_t(dims) ;

                data_ptr_ = nullptr ;
                paramaters_set_ = false ;
                has_error_ = false ;
                
                SetDamarisParameterAndShmem( param_values ) ;  // Initialise the memory size in the constructor.
            }

            ~DamarisVar( void ) {
                  delete [] param_sizes_ ;
                  delete [] positions_ ;
                  if (data_ptr_ != nullptr) 
                  {
                      CommitVariableDamarisShmem() ; 
                      ClearVariableDamarisShmem() ; 
                  }
                  if (this->HasError()) 
                      PrintError() ; // flush out any error messages
            }

            /**
            *  Method to check that the template paramater T is the same as the requested type for the variable in the XML file
            */
            bool TestType(std::string variable_name)
            {
                bool resbool = true ;
                // This gets the type of the Damaris XML <variable>'s <layout>
                DAMARIS_TYPE_STR vartype ;
                dam_err_ = damaris_get_type(variable_name.c_str(), &vartype) ;
                if (dam_err_ != DAMARIS_OK) {
                    dam_err_sstr_ << "  ERROR rankDamarisVar::DamarisVar ()  damaris_get_type(\"" 
                                     << variable_name_ << "\", vartype);  Damaris error = " <<  damaris_error_string(dam_err_) << std::endl ;
                    has_error_ = true ;
                    return( false ) ;
                }
                T test_id;
                const std::type_info& t1 = typeid(test_id) ;
                
                if (vartype == DAMARIS_TYPE_DOUBLE)
                {
                    double td = 0.0 ;
                    const std::type_info& t2 = typeid(td) ;
                    if (t1 != t2) {
                        OutputErrorAndAssert(variable_name, t1.name(), t2.name()) ;
                        resbool = false ;
                    }
                } 
                else if (vartype == DAMARIS_TYPE_FLOAT) 
                {
                    float td = 0.0f ;
                    const std::type_info& t2 = typeid(td) ;
                    if (t1 != t2) {
                        OutputErrorAndAssert(variable_name, t1.name(), t2.name()) ;
                        resbool = false ;
                    }
                } 
                else if (vartype == DAMARIS_TYPE_CHAR)
                {
                    char td = 0;
                    const std::type_info& t2 = typeid(td) ;
                    if (t1 != t2) {
                        OutputErrorAndAssert(variable_name, t1.name(), t2.name()) ;
                        resbool = false ;
                    }
                }
                else if (vartype == DAMARIS_TYPE_UCHAR)
                {
                    unsigned char td = 0;
                    const std::type_info& t2 = typeid(td) ;
                    if (t1 != t2) {
                        OutputErrorAndAssert(variable_name, t1.name(), t2.name()) ;
                        resbool = false ;
                    }
                }
                else if (vartype == DAMARIS_TYPE_SHORT)
                {
                    short td = 0;
                    const std::type_info& t2 = typeid(td) ;
                    if (t1 != t2) {
                        OutputErrorAndAssert(variable_name, t1.name(), t2.name()) ;
                        resbool = false ;
                    }
                }
                else if (vartype == DAMARIS_TYPE_USHORT)
                {
                    unsigned short td = 0;
                    const std::type_info& t2 = typeid(td) ;
                    if (t1 != t2) {
                        OutputErrorAndAssert(variable_name, t1.name(), t2.name()) ;
                        resbool = false ;
                    }
                }
                else if (vartype == DAMARIS_TYPE_INT)
                {
                    int td = 0;
                    const std::type_info& t2 = typeid(td) ;
                    if (t1 != t2) {
                        OutputErrorAndAssert(variable_name, t1.name(), t2.name()) ;
                        resbool = false ;
                    }
                }
                else if (vartype == DAMARIS_TYPE_UINT)
                {
                    unsigned int td = 0;
                    const std::type_info& t2 = typeid(td) ;
                    if (t1 != t2) {
                        OutputErrorAndAssert(variable_name, t1.name(), t2.name()) ;
                        resbool = false ;
                    }
                }
                else if (vartype == DAMARIS_TYPE_LONG)
                {
                    long td = 0;
                    const std::type_info& t2 = typeid(td) ;
                    if (t1 != t2) {
                        OutputErrorAndAssert(variable_name, t1.name(), t2.name()) ;
                        resbool = false ;
                    }
                }
                else if (vartype == DAMARIS_TYPE_ULONG)
                {
                    unsigned long td = 0;
                    const std::type_info& t2 = typeid(td) ;
                    if (t1 != t2) {
                        OutputErrorAndAssert(variable_name, t1.name(), t2.name()) ;
                        resbool = false ;
                    }
                }
                else if (vartype == DAMARIS_TYPE_UNDEFINED)
                {
                    std::cerr << "ERROR rank =" << rank_ << " : DamarisVar::DamarisVar()::  \"" << variable_name << "\" has type DAMARIS_TYPE_UNDEFINED"<< std::endl ;
                    resbool = false ;
                } 
                else 
                {
                    std::cerr << "ERROR rank =" << rank_ << " : DamarisVar::DamarisVar():: \"" << variable_name << "\" is not of available type "<< std::endl ;
                    resbool = false ;
                }
                
                return resbool ;
            }

            /**
            *  Allow a user to indicate that the Damaris variable has allocated a size -
             * This method is usefull as a single paramater can control one or more layouts and
             * a single layout can describe the size of multiple <variable> elements.
             * i.e. The current variable has had it's paramater(s) set through via another variable.
            */
            void ParameterIsSet() 
            {
                paramaters_set_ = true ;
            }

            void PrintError ( void )
            {
                std::cerr << dam_err_sstr_.str() ;
            }

            bool HasError( void )
            {
               return (has_error_) ;
            }

            /**
            *  Returns the data pointer to shared memory, or nullptr if it has not been allocated
            */
            T * data_ptr( void ) 
            {
               return (data_ptr_) ;
            }

            std::string & variable_name( void )
            {
               return (variable_name_) ;
            }

             /**
            * Creates the XML representation of the variable from the available strings
            */
            std::string ReturnXMLForVariable ( void )
            {
                std::ostringstream  var_sstr ;

                var_sstr << "<variable " << " name=\"" << variable_name_ << "\"" ;
                var_sstr << xml_attributes_.ReturnXMLForVariable() ;
                var_sstr << " /> " ;

                return var_sstr.str() ;
            }

            /**
            *  Method to set the Damaris paramater values and set the shmem region \ref data_ptr_
            *
            *  /param [IN] paramSizeVal : A vector of values to set the Damaris paramters to. One element per param_names_ string
            *
            *
            */
            void SetDamarisParameterAndShmem( std::vector<int> paramSizeVal ) 
            {
                this->SetDamarisParameter( paramSizeVal )  ;
                this->SetPointersToDamarisShmem() ;
            }

            /**
            *  Method to set the Damaris paramater values.
            *
            *  /param [IN] paramSizeVal : An pointer to a value or array of values to set. One element per param_names_ string
            *
            *  /implicit                : Implicitly uses the array of paramater names: \ref param_names_
            */
            void SetDamarisParameter( std::vector<int>& paramSizeVal ) 
            {
                assert(paramSizeVal.size() == num_params_) ;

                bool resbool = true ;
                for (int varnum = 0 ; varnum < num_params_ ; varnum++)
                {
                    param_sizes_[varnum] = paramSizeVal[varnum] ;

                    dam_err_ = damaris_parameter_set(param_names_[varnum].c_str(), &paramSizeVal[varnum], sizeof(int));
                    if (dam_err_ != DAMARIS_OK) {
                        dam_err_sstr_ << "  ERROR rank =" << rank_ << " : class DamarisVar : damaris_parameter_set(\"" << param_names_[varnum] 
                                          << "\", paramSizeVal, sizeof(int));  Damaris error = " <<  damaris_error_string(dam_err_) << std::endl ;
                        resbool = false ;
                        has_error_ = true ;
                    }
                }
                if (resbool == true) 
                    paramaters_set_ = true ;
            }

            /**
            *  Method to set the Damaris position values.
            *
            *  /param [IN] positionsVals : An pointer to a value or array of values to set as the offset into the array. 
            *                              One element per dimension (one value for each dim_)
            *
            *  /implicit                 : Implicitly uses the variable name: \ref variable_name_
            */
            void SetDamarisPosition( std::vector<int64_t> positionsVals ) 
            {
                assert(positionsVals.size() == dims_) ;

                for (int pos_dim = 0 ; pos_dim < dims_ ; pos_dim++)
                {
                    positions_[pos_dim] = positionsVals[pos_dim] ;
                }
                dam_err_ = damaris_set_position(variable_name_.c_str(), positionsVals.data());
                if (dam_err_ != DAMARIS_OK) {
                    dam_err_sstr_ << "  ERROR rank =" << rank_ << " : class DamarisVar : damaris_set_position(\"" 
                                     << variable_name_ << "\", positionsVals);  Damaris error = " <<  damaris_error_string(dam_err_) << std::endl ;
                    has_error_ = true ;
                }
            }

            /**
            *  Method to set the Damaris paramater values.
            *
            *  /implicit                : Implicitly uses the Damaris variable name string  \ref variable_name_
            *  /implicit                : Implicitly uses the class data element : \ref data_ptr_
            */
            void SetPointersToDamarisShmem( void )
            {
                if (paramaters_set_ == true )
                {
                    // Allocate memory in the shared memory section... 
                    dam_err_ = damaris_alloc(variable_name_.c_str(), (void **) &data_ptr_) ;
                    if (dam_err_ != DAMARIS_OK) {
                        dam_err_sstr_ << "  ERROR rank =" << rank_ << " : class DamarisVar : damaris_alloc(\"" 
                                      << variable_name_ <<"\", (void **) &ret_ptr)" << ", Damaris error = " <<  damaris_error_string(dam_err_) << std::endl ;
                        has_error_ = true ;
                    }
                } else {
                    dam_err_ = -1 ;
                    dam_err_sstr_ << "ERROR rank =" << rank_  << " : class DamarisVar : SetDamarisParameter() should be called first so as to define the size of the memory block required for variable : " << variable_name_ << std::endl ;
                    has_error_ = true ;
                }
            }

            void SetPointersToDamarisShmem( T ** ptr_in )
            {
                if (paramaters_set_ == true )
                {
                    T * temp_ptr ;
                    // Allocate memory in the shared memory section... 
                    dam_err_ = damaris_alloc(variable_name_.c_str(), (void **) &temp_ptr) ;
                    if (dam_err_ != DAMARIS_OK) {
                        dam_err_sstr_ << "  ERROR rank =" << rank_ << " : class DamarisVar : damaris_alloc(\"" 
                                          << variable_name_ <<"\", (void **) &ret_ptr)" << ", Damaris error = " <<  damaris_error_string(dam_err_) << std::endl ;
                        has_error_ = true ;
                    }

                    *ptr_in = temp_ptr ;
                    data_ptr_ = temp_ptr ;
                } else {
                    dam_err_ = -1 ;
                    dam_err_sstr_ << "  ERROR rank =" << rank_  << " : class DamarisVar : SetDamarisParameter() should be called first so as to define the size of the memory block required for variable : " << variable_name_ << std::endl ;
                    has_error_ = true ;
                }
            }

            /**
            *  Method to commit the memory of the data written to the Damaris variable - 
            *  Indicates that we will not write any more data to \ref data_ptr_
            *
            *  /implicit                : Implicitly uses the variable name string  \ref variable_name_
            */
            void CommitVariableDamarisShmem( void )
            {
                // Signal to Damaris we are done writing data for this iteration
                dam_err_ = damaris_commit (variable_name_.c_str()) ;
                if (dam_err_ != DAMARIS_OK) {
                    dam_err_sstr_ << "  ERROR rank =" << rank_ << " : class DamarisVar : damaris_commit(\"" 
                                      << variable_name_ <<"\")"  << ", Damaris error = " <<  damaris_error_string(dam_err_) << std::endl ;
                    has_error_ = true ;
                }
            }

            /**
            *  Method to release the memory of the data written to the Damaris variable - 
            *  Indicates that Damaris may take control of the shared memory area that was used for the variable \ref data_ptr_
            *
            *  /implicit                : Implicitly uses the variable name string  \ref variable_name_
            */
            void ClearVariableDamarisShmem( void )
            {
                // Signal to Damaris it has complete charge of the memory area
                dam_err_ = damaris_clear(variable_name_.c_str()) ;
                if (dam_err_ != DAMARIS_OK) {
                    dam_err_sstr_ << "  ERROR rank =" << rank_ << " : class DamarisVar : damaris_clear(\"" 
                                      << variable_name_ << "\")"  << ", Damaris error = " <<  damaris_error_string(dam_err_) << std::endl ;
                    has_error_ = true ;
                }
                data_ptr_ = nullptr ;
            }
        private:
            void OutputErrorAndAssert(std::string& var_name, std::string type_name1, std::string type_name2)
            {
                dam_err_sstr_ << "ERROR rank =" << rank_ << " : DamarisVar::DamarisVar () variable_name_: \"" << var_name 
                    << "\" The template type of Type of DamarisVar<T> in the code: " << type_name1 << " does not match type in XML (float)" << std::endl ;
                PrintError() ;
                assert( type_name1 == type_name2  ) ;
            }
    };  // class DamarisVar

    } // namespace DamarisOutput

} // namespace OPM

#endif 