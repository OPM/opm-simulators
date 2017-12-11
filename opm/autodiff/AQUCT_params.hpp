#ifndef OPM_AQUCT_HEADER_INCLUDED
#define OPM_AQUCT_HEADER_INCLUDED


    struct AQUCT_params{
        
            // Aquifer ID
            int aquiferID;
            // Table IDs
            int inftableID, pvttableID;
            // Perforation cell id
            int cell_id;
            // Variables constants
            double  phi_aq , //aquifer porosity
                    d0,   //aquifer datum depth
                    C_t , //total compressibility
                    r_o , //aquifer inner radius
                    k_a , //aquifer permeability
                    c1, // 0.008527 (METRIC, PVT-M); 0.006328 (FIELD); 3.6 (LAB)
                    h , //aquifer thickness
                    theta , //angle subtended by the aquifer boundary
                    c2 ; //6.283 (METRIC, PVT-M); 1.1191 (FIELD); 6.283 (LAB).
                    std::vector<double> td, pi;
    

    };




#endif            