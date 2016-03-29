/*
  Copyright 2014 IRIS AS

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

#ifndef OPM_BACKUPRESTORE_HEADER_INCLUDED
#define OPM_BACKUPRESTORE_HEADER_INCLUDED

#include <iostream>
#include <opm/common/data/SimulationDataContainer.hpp>


#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/autodiff/WellStateFullyImplicitBlackoil.hpp>

namespace Opm {

    namespace {
        template <class T>
        void writeValue( std::ostream& stream, const T& value )
        {
            stream.write( (const char *) &value, sizeof( T ) );
        }

        template <class T>
        void readValue( std::istream& stream, T& value )
        {
            stream.read( (char *) &value, sizeof( T ) );
        }

        // write any container that provides a size method and a data method
        template <class Vector>
        void writeContainer( std::ostream& stream, const Vector& vector)
        {
            typedef typename Vector :: value_type T;
            unsigned int size = vector.size() * sizeof( T );
            writeValue( stream, size );
            if( size > 0 ) {
                stream.write( (const char *) vector.data(), size );
            }
        }

        // write map where the key is a std::string
        // and the value is a std::vector or std::array
        template <class Map>
        void writeMap( std::ostream& stream, const Map& m)
        {
            const unsigned int mapEntries = m.size();
            writeValue( stream, mapEntries );
            if( mapEntries > 0 )
            {
                const auto& end = m.end();
                for(auto it = m.begin(); it != end; ++it )
                {
                    // write key (assume that key is a container)
                    writeContainer( stream, (*it).first );
                    // write value (assume that value is a container)
                    writeContainer( stream, (*it).second );
                }
            }
        }

        template <class Container>
        void resizeContainer( Container& container, size_t size )
        {
            container.resize( size );
        }

        template <class T, size_t n>
        void resizeContainer( std::array<T, n>& /* a */, size_t size )
        {
            static_cast<void>(size);
            assert( int(size) == int(n) );
        }

        template <class Container>
        void readData(std::istream& stream, Container& container, size_t datasize)
        {
            stream.read( reinterpret_cast<char*>( container.data() ), datasize );
        }

        //We need to be careful with string, because string.data() returns something that
        //"A program shall not alter any of the characters in this sequence."
        template <>
        void readData<std::string>(std::istream& stream, std::string& string, size_t datasize)
        {
            std::vector<char> raw_data(datasize);
            readData(stream, raw_data, datasize);
            string = std::string(raw_data.data(), datasize);
        }

        template <class Container>
        void readContainerImpl( std::istream& stream, Container& container, const bool adjustSize )
        {
            typedef typename Container :: value_type T;
            unsigned int dataSize = 0;
            readValue( stream, dataSize );
            if( adjustSize && dataSize > 0 ) {
                resizeContainer( container, dataSize/sizeof(T) );
            }

            if( dataSize != container.size() * sizeof( T ) )
            {
                OPM_THROW(std::logic_error,
                        "Size of stored data and simulation data does not match "
                        << dataSize << " " << (container.size() * sizeof( T )) );
            }
            if( dataSize > 0 ) {
                readData(stream, container, dataSize);
            }
        }

        template <class Container>
        void readAndResizeContainer( std::istream& stream, Container& container )
        {
            readContainerImpl( stream, container, true );
        }

        template <class Container>
        void readContainer( std::istream& stream, Container& container )
        {
            readContainerImpl( stream, container, false );
        }

        template <class Map>
        void readMap( std::istream& stream, Map& m)
        {
            m.clear();
            unsigned int mapEntries = 0;
            readValue( stream, mapEntries );
            for( unsigned int entry = 0; entry<mapEntries; ++entry )
            {
                std::pair< typename Map :: key_type, typename Map :: mapped_type > mapEntry;
                // read key
                readAndResizeContainer( stream, mapEntry.first );
                // read values
                readContainer( stream, mapEntry.second );

                // insert entry into map
                m.insert( mapEntry );
            }
        }

        enum { SimulatorStateId = 0,
               WellStateId = 1,
               BlackoilStateId = 2,
               WellStateFullyImplicitBackoilId = 3 };

        inline int objectId( const SimulationDataContainer& /* state */) {
            return SimulatorStateId;
        }
        inline int objectId( const WellState& /* state */) {
            return WellStateId;
        }
        inline int objectId( const BlackoilState& /* state */) {
            return BlackoilStateId;
        }
        inline int objectId( const WellStateFullyImplicitBlackoil& /* state */) {
            return WellStateFullyImplicitBackoilId;
        }

        template <class State>
        void checkObjectId( std::istream& in, const State& state )
        {
            // read id and compare with object
            int id = -1;
            readValue( in, id );
            if( id != objectId( state ) ) {
                OPM_THROW(std::logic_error,"backup-restore object type mismatch");
            }
        }
    }

    // SimulationDataContainer
    inline
    std::ostream& operator << (std::ostream& out, const SimulationDataContainer& state )
    {
        // write id of object to stream
        writeValue( out, objectId( state ) );

        const int numPhases = state.numPhases();
        writeValue( out, numPhases );

        // write variables
        writeContainer( out, state.pressure() );
        writeContainer( out, state.temperature() );
        writeContainer( out, state.facepressure() );
        writeContainer( out, state.faceflux() );
        writeContainer( out, state.saturation() );

        return out;
    }

    inline
    std::istream& operator >> (std::istream& in, SimulationDataContainer& state )
    {
        // check id of stored object
        checkObjectId( in, state );

        int numPhases = 0;
        readValue( in, numPhases );

        if( numPhases != (int) state.numPhases() )
            OPM_THROW(std::logic_error,"num phases wrong");

        // read variables
        readContainer( in, state.pressure() );
        readContainer( in, state.temperature() );
        readContainer( in, state.facepressure() );
        readContainer( in, state.faceflux() );
        readContainer( in, state.saturation() );

        return in;
    }

    // BlackoilState
    inline
    std::ostream& operator << (std::ostream& out, const BlackoilState& state )
    {
        // write id of object to stream
        writeValue( out, objectId( state ) );

        // backup simulator state
        const SimulationDataContainer& simstate = static_cast< const SimulationDataContainer& > (state);
        out << simstate;

        // backup additional blackoil state variables
        writeContainer( out, state.surfacevol() );
        writeContainer( out, state.gasoilratio() );
        writeContainer( out, state.rv() );

        return out;
    }

    inline
    std::istream& operator >> (std::istream& in, BlackoilState& state )
    {
        // check id of stored object
        checkObjectId( in, state );

        // restore simulator state
        SimulationDataContainer& simstate = static_cast< SimulationDataContainer& > (state);
        in >> simstate;

        // restore additional blackoil state variables
        readContainer( in, state.surfacevol() );
        readContainer( in, state.gasoilratio() );
        readContainer( in, state.rv() );

        return in;
    }

    // WellState
    inline
    std::ostream& operator << (std::ostream& out, const WellState& state )
    {
        // write id of object to stream
        writeValue( out, objectId( state ) );

        // backup well state
        writeContainer( out, state.bhp() );
        writeContainer( out, state.temperature() );
        writeContainer( out, state.wellRates() );
        writeContainer( out, state.perfRates() );
        writeContainer( out, state.perfPress() );

        return out;
    }

    inline
    std::istream& operator >> (std::istream& in, WellState& state )
    {
        // check id of stored object
        checkObjectId( in, state );

        // restore well state
        readAndResizeContainer( in, state.bhp() );
        readAndResizeContainer( in, state.temperature() );
        readAndResizeContainer( in, state.wellRates() );
        readAndResizeContainer( in, state.perfRates() );
        readAndResizeContainer( in, state.perfPress() );

        return in;
    }

    // WellStateFullyImplicitBlackoil
    inline
    std::ostream& operator << (std::ostream& out, const WellStateFullyImplicitBlackoil& state )
    {
        // write id of object to stream
        writeValue( out, objectId( state ) );

        // backup well state
        const WellState& wellState = static_cast< const WellState& > (state);
        out << wellState;

        const int numWells = state.numWells();
        writeValue( out, numWells );
        if( numWells > 0 )
        {
            const int numPhases = state.numPhases();
            writeValue( out, numPhases );

            // backup additional variables
            writeContainer( out, state.perfPhaseRates() );
            writeContainer( out, state.currentControls() );
            writeMap( out, state.wellMap() );
        }

        return out;
    }

    inline
    std::istream& operator >> (std::istream& in, WellStateFullyImplicitBlackoil& state )
    {
        // check id of stored object
        checkObjectId( in, state );

        // restore well state
        WellState& wellState = static_cast< WellState& > (state);
        in >> wellState;

        int numWells = 0;
        readValue( in, numWells );
        if( numWells != state.numWells() )
            OPM_THROW(std::logic_error,"wrong numWells");

        if( numWells > 0 )
        {
            int numPhases = 0;
            readValue( in, numPhases );
            if( numPhases != state.numPhases() )
                OPM_THROW(std::logic_error,"wrong numPhases");

            // restore additional variables
            readAndResizeContainer( in, state.perfPhaseRates()  );
            readAndResizeContainer( in, state.currentControls() );
            readMap( in, state.wellMap() );
        }

        return in;
    }

} // namespace Opm

#endif // OPM_BACKUPRESTORE_HEADER_INCLUDED
