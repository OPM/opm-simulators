/**/

#include <opm/autodiff/polymer/IncompPropsAdBasic.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/autodiff/AutoDiffBlock.hpp>
#include <opm/autodiff/AutoDiffHelpers.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <iostream>

namespace Opm
{
    IncompPropsAdBasic::IncompPropsAdBasic(const parameter::ParameterGroup& param,
                                           const int dim,
                                           const int num_cells)
    {
        double poro = param.getDefault("porosity", 1.0);
        using namespace Opm::unit;
        using namespace Opm::prefix;
        double perm = param.getDefault("permeability", 100) * milli * darcy;
        rock_.init(dim, num_cells, poro, perm);
        pvt_.init(param);
        satprops_.init(param);
        if (pvt_.numPhases() != satprops_.numPhases()) {
            OPM_THROW(std::runtime_error, "IncompPropsAdBasic::IncompPropsAdBasic() - Inconsistent number of phases in pvt data ("
                  << pvt_.numPhases() << ") and saturation-dependent function data (" << satprops_.numPhases() << ").");
        }
        viscosity_.resize(pvt_.numPhases());
        pvt_.mu(1, 0, 0, &viscosity_[0]);
    }
    IncompPropsAdBasic::IncompPropsAdBasic(const int num_phases,
                                                 const  SaturationPropsBasic::RelPermFunc& relpermfunc,
                                                 const std::vector<double>&  rho,
                                                 const std::vector<double>& mu,
                                                 const double por, //porosity
                                                 const double perm,
                                                 const int dim,
                                                 const int num_cells)
    {
        rock_.init(dim, num_cells, por, perm);
        pvt_.init(num_phases, rho, mu);
        satprops_.init(num_phases, relpermfunc);
        if (pvt_.numPhases() != satprops_.numPhases()) {
            OPM_THROW(std::runtime_error, "IncompPropsAdBasic::IncompPropsAdBasic() - Inconsistent number of phases in pvt data ("
                  << pvt_.numPhases() << ") and saturation-dependent function data (" << satprops_.numPhases() << ").");
        }
        viscosity_.resize(pvt_.numPhases());
        pvt_.mu(1, 0, 0, &viscosity_[0]);
    }
    IncompPropsAdBasic::~IncompPropsAdBasic()
    {
    }

    /// \return   D, the number of spatial dimensions.
    int IncompPropsAdBasic::numDimensions() const
    {
        return rock_.numDimensions();
    }

    /// \return   N, the number of cells.
    int IncompPropsAdBasic::numCells() const
    {
        return rock_.numCells();
    }

    /// \return   Array of N porosity values.
    const double* IncompPropsAdBasic::porosity() const
    {
        return rock_.porosity();
    }

    /// \return   Array of ND^2 permeability values.
    ///           The D^2 permeability values for a cell are organized as a matrix,
    ///           which is symmetric (so ordering does not matter).
    const double* IncompPropsAdBasic::permeability() const
    {
        return rock_.permeability();
    }


    // ---- Fluid interface ----

    /// \return   P, the number of phases (also the number of components).
    int IncompPropsAdBasic::numPhases() const
    {
        return pvt_.numPhases();
    }

    /// \return Array of P viscosity values.
    const double* IncompPropsAdBasic::viscosity() const
    {
        return &viscosity_[0];
    }
    /// \return Array of P density values.
    const double* IncompPropsAdBasic::density() const
    {
        return pvt_.surfaceDensities();
    }

    const double* IncompPropsAdBasic::surfaceDensity() const
    {
        return pvt_.surfaceDensities();
    }
    typedef IncompPropsAdBasic::ADB ADB;
    typedef IncompPropsAdBasic::V V;
    typedef Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Block;
    typedef std::vector<int> Cells;
    std::vector<V> IncompPropsAdBasic::relperm(const V& sw,
                                                  const V& so,
                                                  const Cells& cells) const
    {
        const int n = cells.size();
        const int np = numPhases();
        Block s_all(n, np);
        assert(sw.size() == n && so.size() == n);
        s_all.col(0) = sw;
        s_all.col(1) = so;
        Block kr(n, np);
//        satprops_.relperm(n, s_all.data(), cells.data(), kr.data(), 0);
        satprops_.relperm(n, s_all.data(), kr.data(), 0);

        std::vector<V> relperms;
        relperms.reserve(2);
        for (int phase = 0; phase < 2; ++phase) {
                relperms.emplace_back(kr.col(phase));
        }
        return relperms;
    }

    std::vector<ADB>    IncompPropsAdBasic::relperm(const ADB& sw,
                                                    const ADB& so,
                                                    const Cells& cells) const
    {
        const int n = cells.size();
        const int np = numPhases();
        Block s_all(n, np);
        assert(sw.size() == n && so.size() == n);
        s_all.col(0) = sw.value();
        s_all.col(1) = so.value();
        Block kr(n, np);
        Block dkr(n, np*np);
//        satprops_.relperm(n, s_all.data(), cells.data(), kr.data(), dkr.data());
        satprops_.relperm(n, s_all.data(), kr.data(), dkr.data());
        const int num_blocks = so.numBlocks();
        std::vector<ADB> relperms;
        relperms.reserve(2);
        typedef const ADB* ADBPtr;
        ADBPtr s[2] = { &sw, &so };
        for (int phase1 = 0; phase1 < 2; ++phase1) {
            const int phase1_pos = phase1;
            std::vector<ADB::M> jacs(num_blocks);
            for (int block = 0; block < num_blocks; ++block) {
                jacs[block] = ADB::M(n, s[phase1]->derivative()[block].cols());
            }
            for (int phase2 = 0; phase2 < 2; ++phase2) {
                const int phase2_pos = phase2;
                // Assemble dkr1/ds2.
                const int column = phase1_pos + np*phase2_pos; // Recall: Fortran ordering from props_.relperm()
                ADB::M dkr1_ds2_diag = spdiag(dkr.col(column));
                for (int block = 0; block < num_blocks; ++block) {
                    jacs[block] += dkr1_ds2_diag * s[phase2]->derivative()[block];
                }
            }
                relperms.emplace_back(ADB::function(kr.col(phase1_pos), jacs));
        }
        return relperms;
    }
}
