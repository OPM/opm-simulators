// ========================================================================
// Basically all the header files from flow.cpp
#include "config.h"

#include <flow/flow_ebos_blackoil.hpp>
#include <flow/flow_ebos_gasoil.hpp>
#include <flow/flow_ebos_oilwater.hpp>
#include <flow/flow_ebos_solvent.hpp>
#include <flow/flow_ebos_polymer.hpp>
#include <flow/flow_ebos_energy.hpp>
#include <flow/flow_ebos_oilwater_polymer.hpp>
#include <flow/flow_ebos_oilwater_polymer_injectivity.hpp>

#include <opm/simulators/flow/SimulatorFullyImplicitBlackoilEbos.hpp>
#include <opm/simulators/flow/FlowMainEbos.hpp>
#include <opm/simulators/utils/moduleVersion.hpp>
#include <ewoms/common/propertysystem.hh>
#include <ewoms/common/parametersystem.hh>
#include <opm/simulators/flow/MissingFeatures.hpp>
#include <opm/material/common/ResetLocale.hpp>

#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Parser/ParseContext.hpp>
#include <opm/parser/eclipse/Parser/ErrorGuard.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/checkDeck.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>

#include <opm/parser/eclipse/EclipseState/Schedule/ArrayDimChecker.hpp>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#else
#include <dune/common/parallel/mpihelper.hh>
#endif

// ========================================================================

#include <pybind11/pybind11.h>
#include <string>
#include <iostream>

namespace py = pybind11;

BEGIN_PROPERTIES

// this is a dummy type tag that is used to setup the parameters before the actual
// simulator.
NEW_TYPE_TAG(FlowEarlyBird, INHERITS_FROM(EclFlowProblem));

END_PROPERTIES

namespace detail
{
    boost::filesystem::path simulationCaseName( const std::string& casename ) {
        namespace fs = boost::filesystem;

        const auto exists = []( const fs::path& f ) -> bool {
            if( !fs::exists( f ) ) return false;

            if( fs::is_regular_file( f ) ) return true;

            return fs::is_symlink( f )
            && fs::is_regular_file( fs::read_symlink( f ) );
        };

        auto simcase = fs::path( casename );

        if( exists( simcase ) ) {
            return simcase;
        }

        for( const auto& ext : { std::string("data"), std::string("DATA") } ) {
            if( exists( simcase.replace_extension( ext ) ) ) {
                return simcase;
            }
        }

        throw std::invalid_argument( "Cannot find input case " + casename );
    }


    // This function is an extreme special case, if the program has been invoked
    // *exactly* as:
    //
    //    flow   --version
    //
    // the call is intercepted by this function which will print "flow $version"
    // on stdout and exit(0).
    void handleVersionCmdLine(int argc, char** argv) {
        for ( int i = 1; i < argc; ++i )
        {
            if (std::strcmp(argv[i], "--version") == 0) {
                std::cout << "flow " << Opm::moduleVersionName() << std::endl;
                std::exit(EXIT_SUCCESS);
            }
        }
    }
}

class BlackOilSimulator
{
public:
    BlackOilSimulator()
    {
        argc_ = 2;
        argv_ = new char*[2];
        argv_[0] = new char[200];
        char argv0[] = "flow";
        std::strcpy(argv_[0], argv0);
        argv_[1] = new char[200];
    }

    ~BlackOilSimulator() 
    {
        delete[] argv_[0];
        delete[] argv_[1];
        delete[] argv_;
    }

    void setDeckFilename( char* deckFilename )
    {
        std::strcpy(argv_[1], deckFilename);
    }

    void setDeck( const Opm::Deck& deck )
    {
        deck_ = std::make_shared< Opm::Deck >(deck);
        // deck_ = std::shared_ptr< Opm::Deck >( const_cast< Opm::Deck* > (&deck), [](Opm::Deck*){} );
    }

    void setEclipseState( const Opm::EclipseState& eclipseState )
    {
        eclipseState_ = std::make_shared< Opm::EclipseState >(eclipseState);
        // eclipseState_ = std::shared_ptr< Opm::EclipseState >( const_cast< Opm::EclipseState* > (&eclipseState), [](Opm::EclipseState*){} );
    }

    void setSchedule( const Opm::Schedule& schedule )
    {
        schedule_ = std::make_shared< Opm::Schedule >(schedule);
        // schedule_ = std::shared_ptr< Opm::Schedule >( const_cast< Opm::Schedule* > (&schedule), [](Opm::Schedule*){} );
    }

    void setSummaryConfig( const Opm::SummaryConfig& summaryConfig )
    {
        summaryConfig_ = std::make_shared< Opm::SummaryConfig >(summaryConfig);
        // summaryConfig_ = std::shared_ptr< Opm::SummaryConfig >( const_cast< Opm::SummaryConfig* > (&summaryConfig), [](Opm::SummaryConfig*){} );
    }

    int run()
    {
        int argc = argc_;
        char **argv = argv_;

        Dune::Timer externalSetupTimer;
        externalSetupTimer.start();

        detail::handleVersionCmdLine(argc, argv);
        // MPI setup.
#if HAVE_DUNE_FEM
        Dune::Fem::MPIManager::initialize(argc, argv);
        int mpiRank = Dune::Fem::MPIManager::rank();
#else
        // the design of the plain dune MPIHelper class is quite flawed: there is no way to
        // get the instance without having the argc and argv parameters available and it is
        // not possible to determine the MPI rank and size without an instance. (IOW: the
        // rank() and size() methods are supposed to be static.)
        const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);
        int mpiRank = mpiHelper.rank();
#endif

        // we always want to use the default locale, and thus spare us the trouble
        // with incorrect locale settings.
        Opm::resetLocale();

        // this is a work-around for a catch 22: we do not know what code path to use without
        // parsing the deck, but we don't know the deck without having access to the
        // parameters and this requires to know the type tag to be used. To solve this, we
        // use a type tag just for parsing the parameters before we instantiate the actual
        // simulator object. (Which parses the parameters again, but since this is done in an
        // identical manner it does not matter.)
        typedef TTAG(FlowEarlyBird) PreTypeTag;
        typedef GET_PROP_TYPE(PreTypeTag, Problem) PreProblem;

        PreProblem::setBriefDescription("Flow, an advanced reservoir simulator for ECL-decks provided by the Open Porous Media project.");

        int status = Opm::FlowMainEbos<PreTypeTag>::setupParameters_(argc, argv);
        if (status != 0)
            // if setupParameters_ returns a value smaller than 0, there was no error, but
            // the program should abort. This is the case e.g. for the --help and the
            // --print-properties parameters.
            return (status >= 0)?status:0;

        bool outputCout = false;
        if (mpiRank == 0)
            outputCout = EWOMS_GET_PARAM(PreTypeTag, bool, EnableTerminalOutput);

//         std::string deckFilename = EWOMS_GET_PARAM(PreTypeTag, std::string, EclDeckFileName);
//         typedef typename GET_PROP_TYPE(PreTypeTag, Vanguard) PreVanguard;
//         try {
//             deckFilename = PreVanguard::canonicalDeckPath(deckFilename).string();
//         }
//         catch (const std::exception& e) {
//             if ( mpiRank == 0 )
//                 std::cerr << "Exception received: " << e.what() << ". Try '--help' for a usage description.\n";
// #if HAVE_MPI
//             MPI_Finalize();
// #endif
//             return 1;
//         }

        if (outputCout) {
            Opm::FlowMainEbos<PreTypeTag>::printBanner();
        }

        const auto& phases = Opm::Runspec(*deck_).phases();

        // Blackoil case
        if( phases.size() == 3 ) {
            Opm::flowEbosBlackoilSetDeck(externalSetupTimer.elapsed(), *deck_, *eclipseState_, *schedule_, *summaryConfig_);
            return Opm::flowEbosBlackoilMain(argc, argv);
        }
        else
        {
            if (outputCout)
                std::cerr << "Only blackoil configuration is supported!" << std::endl;
            return EXIT_FAILURE;
        }
        
        return EXIT_SUCCESS;
    }

private:
    int argc_;
    char **argv_;

    std::shared_ptr<Opm::Deck>          deck_;
    std::shared_ptr<Opm::EclipseState>  eclipseState_;
    std::shared_ptr<Opm::Schedule>      schedule_;
    std::shared_ptr<Opm::SummaryConfig> summaryConfig_;
};

PYBIND11_MODULE(simulators, m)
{
    py::class_<BlackOilSimulator>(m, "BlackOilSimulator")
        .def(py::init<>())
        .def("run", &BlackOilSimulator::run)
        .def("setDeck", &BlackOilSimulator::setDeck)
        .def("setDeckFilename", &BlackOilSimulator::setDeckFilename)
        .def("setEclipseState", &BlackOilSimulator::setEclipseState)
        .def("setSchedule", &BlackOilSimulator::setSchedule)
        .def("setSummaryConfig", &BlackOilSimulator::setSummaryConfig);
}
