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

#include <opm/core/simulator/SimulatorState.hpp>
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

        template <class Vector>
        void writeVector( std::ostream& stream, const Vector& vector)
        {
            typedef typename Vector :: value_type T;
            unsigned int size = vector.size() * sizeof( T );
            writeValue( stream, size );
            if( size > 0 ) {
                stream.write( (const char *) vector.data(), size );
            }
        }

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
                    // write key
                    writeVector( stream, (*it).first );
                    // write value
                    writeVector( stream, (*it).second );
                }
            }
        }

        template <class Vector>
        void resizeVector( Vector& vector, size_t size )
        {
            vector.resize( size );
        }

        template <class T, unsigned long n>
        void resizeVector( std::array<T, n>& vector, size_t size )
        {
            assert( int(size) == int(n) );
        }

        template <class Vector>
        void readVectorImpl( std::istream& stream, Vector& vector, const bool adjustSize )
        {
            typedef typename Vector :: value_type T;
            unsigned int dataSize = 0;
            readValue( stream, dataSize );
            if( adjustSize && dataSize > 0 ) {
                resizeVector( vector, dataSize/sizeof(T) );
            }

            if( dataSize != vector.size() * sizeof( T ) )
            {
                OPM_THROW(std::logic_error,"Size of stored data and simulation data does not match" << dataSize << " " << (vector.size() * sizeof( T )) );
            }
            if( dataSize > 0 ) {
                stream.read( (char *) vector.data(), dataSize );
            }
        }

        template <class Vector>
        void readAndResizeVector( std::istream& stream, Vector& vector )
        {
            readVectorImpl( stream, vector, true );
        }

        template <class Vector>
        void readVector( std::istream& stream, Vector& vector )
        {
            readVectorImpl( stream, vector, false );
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
                readAndResizeVector( stream, mapEntry.first );
                // read values
                readVector( stream, mapEntry.second );

                // insert entry into map
                m.insert( mapEntry );
            }
        }

        enum { SimulatorStateId = 0,
               WellStateId = 1,
               BlackoilStateId = 2,
               WellStateFullyImplicitBackoilId = 3 };

        inline int objectId( const SimulatorState& state ) {
            return SimulatorStateId;
        }
        inline int objectId( const WellState& state ) {
            return WellStateId;
        }
        inline int objectId( const BlackoilState& state ) {
            return BlackoilStateId;
        }
        inline int objectId( const WellStateFullyImplicitBlackoil& state ) {
            return WellStateFullyImplicitBackoilId;
        }

        template <class State>
        void checkObjectId( std::istream& in, const State& state )
        {
            // read id and compare with object
            int id = -1;
            readValue( in, id );
            if( id != objectId( state ) )
                OPM_THROW(std::logic_error,"backup-restore object type missmatch");
        }
    }

    // SimulatorState
    inline
    std::ostream& operator << (std::ostream& out, const SimulatorState& state )
    {
        // write id of object to stream
        writeValue( out, objectId( state ) );

        const int numPhases = state.numPhases();
        writeValue( out, numPhases );

        // write variables
        writeVector( out, state.pressure() );
        writeVector( out, state.temperature() );
        writeVector( out, state.facepressure() );
        writeVector( out, state.faceflux() );
        writeVector( out, state.saturation() );

        return out;
    }

    inline
    std::istream& operator >> (std::istream& in, SimulatorState& state )
    {
        // check id of stored object
        checkObjectId( in, state );

        int numPhases = 0;
        readValue( in, numPhases );

        if( numPhases != state.numPhases() )
            OPM_THROW(std::logic_error,"num phases wrong");

        // read variables
        readVector( in, state.pressure() );
        readVector( in, state.temperature() );
        readVector( in, state.facepressure() );
        readVector( in, state.faceflux() );
        readVector( in, state.saturation() );

        return in;
    }

    // BlackoilState
    inline
    std::ostream& operator << (std::ostream& out, const BlackoilState& state )
    {
        // write id of object to stream
        writeValue( out, objectId( state ) );

        // backup simulator state
        const SimulatorState& simstate = static_cast< const SimulatorState& > (state);
        out << simstate;

        // backup additional blackoil state variables
        writeVector( out, state.surfacevol() );
        writeVector( out, state.gasoilratio() );
        writeVector( out, state.rv() );

        return out;
    }

    inline
    std::istream& operator >> (std::istream& in, BlackoilState& state )
    {
        // check id of stored object
        checkObjectId( in, state );

        // backup simulator state
        SimulatorState& simstate = static_cast< SimulatorState& > (state);
        in >> simstate;

        // backup additional blackoil state variables
        readVector( in, state.surfacevol() );
        readVector( in, state.gasoilratio() );
        readVector( in, state.rv() );

        return in;
    }

    // WellState
    inline
    std::ostream& operator << (std::ostream& out, const WellState& state )
    {
        // write id of object to stream
        writeValue( out, objectId( state ) );

        // backup well state
        writeVector( out, state.bhp() );
        writeVector( out, state.temperature() );
        writeVector( out, state.wellRates() );
        writeVector( out, state.perfRates() );
        writeVector( out, state.perfPress() );

        return out;
    }

    inline
    std::istream& operator >> (std::istream& in, WellState& state )
    {
        // check id of stored object
        checkObjectId( in, state );

        // restore well state
        readAndResizeVector( in, state.bhp() );
        readAndResizeVector( in, state.temperature() );
        readAndResizeVector( in, state.wellRates() );
        readAndResizeVector( in, state.perfRates() );
        readAndResizeVector( in, state.perfPress() );

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
            writeVector( out, state.perfPhaseRates() );
            writeVector( out, state.currentControls() );
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
            readAndResizeVector( in, state.perfPhaseRates()  );
            readAndResizeVector( in, state.currentControls() );
            readMap( in, state.wellMap() );
        }

        return in;
    }

} // namespace Opm

#endif // OPM_BACKUPRESTORE_HEADER_INCLUDED
