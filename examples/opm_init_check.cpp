/*
  Copyright 2014 Statoil ASA

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

#include <iostream>
#include <fstream>
#include <memory>
#include <vector>
#include <string>

#include <ert/ecl/ecl_file.h>
#include <ert/ecl/ecl_kw.h>
#include <ert/ecl/ecl_grid.h>
#include <ert/ecl/ecl_nnc_export.h>

#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Grid/TransMult.hpp>

#include <opm/core/grid.h>
#include <opm/core/grid/GridManager.hpp>
#include <opm/core/wells/WellsManager.hpp>
#include <opm/core/props/BlackoilPropertiesFromDeck.hpp>


#include <opm/autodiff/GeoProps.hpp>
#include <opm/autodiff/SimulatorFullyImplicitBlackoil.hpp>
#include <opm/autodiff/BlackoilPropsAdInterface.hpp>
#include <opm/autodiff/BlackoilPropsAdFromDeck.hpp>


using namespace Opm;




class CellTrans {
public: 
    CellTrans(const std::tuple<int,int,int>& ijk) :
        m_ijk(ijk) 
    {
        
    }
    
    void update(const std::tuple<int,int,int>& ijk , double trans) {
        auto iter = m_trans.find( ijk );
        if (iter == m_trans.end())
            m_trans.insert( std::pair<std::tuple<int,int,int> , double>(ijk , trans));
        else
            iter->second *= trans;
    }
    
    int numConnections() const {
        return m_trans.size();
    }
    

    static void jsonDumpIJK(std::ostream& os , const std::tuple<int,int,int>& ijk ) {
        os << "[" << std::get<0>(ijk) << "," << std::get<1>(ijk) << "," << std::get<2>(ijk) << "]";
    }

    
    void jsonDump(std::ostream& os , bool first) {
        if (!first)
            os <<",";

        
        os << "[";
        jsonDumpIJK(os , m_ijk );
        os << ",[";
        {
            size_t count = 0;
            for (auto transIter = m_trans.begin(); transIter != m_trans.end(); ++transIter) {
                std::tuple<int,int,int> ijk = transIter->first;
                double t = transIter->second;
            
                os << "[";
                jsonDumpIJK( os ,  ijk );
                os << "," << t << "]";
                count++;
                
                if (count < m_trans.size())
                    os << ",";
                else
                    os << "]";
            }
        }
        os << "]" << std::endl;
    }


    std::tuple<int,int,int> m_ijk;
    std::map<std::tuple<int,int,int> , double> m_trans;
};
    



class TransGraph {
public:
    TransGraph(int nx , int ny , int nz) :
        m_nx(nx),
        m_ny(ny),
        m_nz(nz) 
    {
        m_transVector.resize( nx*ny*nz );
    }


    void update(int global_index1 , int global_index2 , double trans) {
        int size = m_nx * m_ny * m_nz;
        if ((global_index1 >= 0)   && 
            (global_index2 >= 0)   &&
            (global_index1 < size) && 
            (global_index2 < size)) {
            
            size_t g1 = std::min( global_index1 , global_index2 );
            size_t g2 = std::max( global_index1 , global_index2 );

            std::shared_ptr<CellTrans> cellTrans = m_transVector[g1];
            if (!cellTrans) {
                cellTrans = std::make_shared<CellTrans>( getIJK(g1) );
                m_transVector[g1] = cellTrans;
            }
            cellTrans->update( getIJK(g2) , trans );
            
        }
    }


    int activeCells() const {
        int count = 0;
        for (size_t g= 0; g < m_transVector.size(); g++) {
            std::shared_ptr<CellTrans> cellTrans = m_transVector[g];
            if (cellTrans)
                count++;
        }
        return count;
    }


    int activeConnections() const {
        int count = 0;
        for (size_t g= 0; g < m_transVector.size(); g++) {
            std::shared_ptr<CellTrans> cellTrans = m_transVector[g];
            if (cellTrans)
                count += cellTrans->numConnections();
        }
        return count;
    }

    

    std::tuple<int,int,int> getIJK(int g) const {
        int k = g / (m_nx * m_ny);
        int j = (g - k*m_nx*m_ny) / m_nx;
        int i = g - k*m_nx*m_ny - j*m_nx;

        return std::tuple<int,int,int>(i,j,k);
    }


    void jsonDump(const std::string& outputFile) {
        std::ofstream os;
        bool first = true;
        os.open(outputFile.c_str());
        
        os << "{ \"dims\" : [" << m_nx << "," << m_ny << "," << m_nz << "]" << " , \"graph\":[";
        for (size_t g= 0; g < m_transVector.size(); g++) {
            std::shared_ptr<CellTrans> cellTrans = m_transVector[g];
            if (cellTrans) {
                cellTrans->jsonDump(os , first);
                first = false;
            }
        }
        os << "]}" << std::endl;
        os.close();
    }



    int m_nx;
    int m_ny;
    int m_nz;
    std::vector<std::shared_ptr<CellTrans> > m_transVector;
};


/*****************************************************************/

void initOPMTrans(TransGraph& opmTrans , DeckConstPtr deck , std::shared_ptr<const EclipseState> eclipseState) {
    std::shared_ptr<GridManager> grid = std::make_shared<GridManager>( eclipseState->getEclipseGrid(), eclipseState->getDoubleGridProperty("PORV")->getData());
    const struct UnstructuredGrid * cGrid = grid->c_grid();
    std::shared_ptr<BlackoilPropsAdInterface> props;

    props.reset(new BlackoilPropsAdFromDeck(deck, eclipseState, *grid->c_grid()));
    DerivedGeology geology(*grid->c_grid() , *props, eclipseState);
    const double * opm_trans_data = geology.transmissibility().data();
    double SIconversion = Opm::unit::cubic(Opm::unit::meter) * Opm::unit::day * Opm::unit::barsa / (Opm::prefix::centi * Opm::unit::Poise);
    
    {
        for (int face_index = 0; face_index < cGrid->number_of_faces; face_index++ ) {
            int global_index1 = cGrid->global_cell[ cGrid->face_cells[2*face_index] ];
            int global_index2 = cGrid->global_cell[ cGrid->face_cells[2*face_index + 1] ];

            opmTrans.update( global_index1 , global_index2 , opm_trans_data[ face_index ] * SIconversion );
        }
    }
}


void initEclipseTrans(TransGraph& eclipseTrans , const ecl_grid_type * ecl_grid , const ecl_file_type * ecl_init) {
    int nx = ecl_grid_get_nx( ecl_grid );
    int ny = ecl_grid_get_ny( ecl_grid );
    int nz = ecl_grid_get_nz( ecl_grid );

    if (ecl_file_has_kw( ecl_init , "TRANX")) {
        ecl_kw_type * tranx_kw = ecl_file_iget_named_kw( ecl_init , "TRANX" , 0 );
        ecl_kw_type * trany_kw = ecl_file_iget_named_kw( ecl_init , "TRANY" , 0 );
        ecl_kw_type * tranz_kw = ecl_file_iget_named_kw( ecl_init , "TRANZ" , 0 );
        for (int k=0; k < nz; k++) {
            for (int j= 0; j < ny; j++) {
                for (int i=0; i < nx; i++) {
                    if (ecl_grid_cell_active3( ecl_grid , i , j , k )) {
                        size_t g1 = ecl_grid_get_global_index3( ecl_grid , i , j , k );
                        int a = ecl_grid_get_active_index1( ecl_grid , g1 );
                        if (a >= 0) {
                            if (i < (nx - 1) && ecl_grid_cell_active3( ecl_grid , i + 1 , j , k)) {
                                size_t g2 = ecl_grid_get_global_index3( ecl_grid , i + 1, j , k );
                                eclipseTrans.update( g1 , g2     , ecl_kw_iget_float( tranx_kw , a ));                            
                            }
                            
                            
                            if (j < (ny - 1) && ecl_grid_cell_active3( ecl_grid , i , j + 1, k)) {
                                size_t g2 = ecl_grid_get_global_index3( ecl_grid , i , j + 1, k );
                                eclipseTrans.update( g1 , g2     , ecl_kw_iget_float( trany_kw , a ));                            
                            }
                            
                            
                            if (k < (nz - 1) && ecl_grid_cell_active3( ecl_grid , i , j , k + 1)) {
                                size_t g2 = ecl_grid_get_global_index3( ecl_grid , i , j , k + 1 );
                                eclipseTrans.update( g1 , g2     , ecl_kw_iget_float( tranz_kw , a ));                            
                            }
                        }
                    }
                }
            }
        }
    } else
        std::cerr << "Init file does not have TRAN[XYZ] keywords" << std::endl;
    
    if (ecl_file_has_kw( ecl_init , "TRANX-")) {
        ecl_kw_type * tranxm_kw = ecl_file_iget_named_kw( ecl_init , "TRANX-" , 0 );
        ecl_kw_type * tranym_kw = ecl_file_iget_named_kw( ecl_init , "TRANY-" , 0 );
        ecl_kw_type * tranzm_kw = ecl_file_iget_named_kw( ecl_init , "TRANZ-" , 0 );
        for (int k=0; k < nz; k++) {
            for (int j= 0; j < ny; j++) {
                for (int i=0; i < nx; i++) {
                    if (ecl_grid_cell_active3( ecl_grid , i , j , k )) {
                        size_t g1 = ecl_grid_get_global_index3( ecl_grid , i , j , k );
                        int a = ecl_grid_get_active_index1( ecl_grid , g1 );

                        if (a >= 0) {
                            if (i > 0 && ecl_grid_cell_active3( ecl_grid , i - 1 , j , k)) {
                                size_t g2 = ecl_grid_get_global_index3( ecl_grid , i - 1, j , k );
                                eclipseTrans.update( g1 , g2     , ecl_kw_iget_float( tranxm_kw , a ));                            
                            }
                            
                            
                            if (j > 0 && ecl_grid_cell_active3( ecl_grid , i , j - 1, k)) {
                                size_t g2 = ecl_grid_get_global_index3( ecl_grid , i , j - 1, k );
                                eclipseTrans.update( g1 , g2     , ecl_kw_iget_float( tranym_kw , a ));                            
                            }
                            
                            
                            if (k > 0 && ecl_grid_cell_active3( ecl_grid , i , j , k - 1)) {
                                size_t g2 = ecl_grid_get_global_index3( ecl_grid , i , j , k - 1 );
                                eclipseTrans.update( g1 , g2     , ecl_kw_iget_float( tranzm_kw , a ));                            
                            }
                        }
                    }
                }
            }
        }
    }

    // NNC
    {
        size_t num_nnc = static_cast<size_t>( ecl_nnc_export_get_size( ecl_grid ));
        std::vector<ecl_nnc_type> nnc(num_nnc);
        
        ecl_nnc_export( ecl_grid , ecl_init , nnc.data());
        for (auto nnc_iter = nnc.begin(); nnc_iter != nnc.end(); ++nnc_iter) 
            eclipseTrans.update( nnc_iter->global_index1 , nnc_iter->global_index2 , nnc_iter->trans );
    }
}



void dump_transGraph( DeckConstPtr deck , std::shared_ptr<const EclipseState> eclipseState , const ecl_grid_type * ecl_grid , const ecl_file_type * ecl_init , size_t verbosity) {
    int nx = ecl_grid_get_nx( ecl_grid );
    int ny = ecl_grid_get_ny( ecl_grid );
    int nz = ecl_grid_get_nz( ecl_grid );
    TransGraph opmTrans(nx , ny , nz );
    TransGraph eclipseTrans( nx , ny , nz);
    
    initOPMTrans( opmTrans , deck , eclipseState );
    initEclipseTrans( eclipseTrans , ecl_grid , ecl_init );
    opmTrans.jsonDump("opm_trans.json");
    eclipseTrans.jsonDump("eclipse_trans.json");
}



int main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "The opm_init_check program needs three arguments:" << std::endl <<  std::endl;;
        std::cerr << "    ECLIPSE.DATA  ECLIPSE.INIT   ECLIPSE.EGRID" << std::endl << std::endl;
        std::cerr << "Where the ECLIPSE.INIT and ECLIPSE.EGRID are existing binary files";
        exit(1);
    }

    std::string input_file = argv[1];
    std::string init_file = argv[2];
    std::string grid_file = argv[3];
    
    ParserPtr parser(new Parser());

    std::cout << "Parsing input file ............: " << input_file << std::endl;
    DeckConstPtr deck = parser->parseFile(input_file, false);
    std::shared_ptr<EclipseState> state = std::make_shared<EclipseState>( deck );
    
    std::cout << "Loading eclipse INIT file .....: " << init_file << std::endl;
    ecl_file_type * ecl_init = ecl_file_open( init_file.c_str() , 0 );

    std::cout << "Loading eclipse EGRID file ....: " << grid_file << std::endl;
    ecl_grid_type * ecl_grid = ecl_grid_alloc( grid_file.c_str() );

    dump_transGraph( deck , state , ecl_grid , ecl_init , 3);
    
    ecl_file_close( ecl_init );
    ecl_grid_free( ecl_grid );
    return 0;
}

