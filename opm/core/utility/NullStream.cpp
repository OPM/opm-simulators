// Copyright (C) 2013 Uni Research AS
// This file is licensed under the GNU General Public License v3.0

#include <opm/core/utility/NullStream.hpp>
#include <ostream>
#include <streambuf>

// buffer that ignores everything
// see <http://forums.codeguru.com/showthread.php?460071-ostream-bit-bucket>
struct NullBuf : public std::streambuf {};
static NullBuf null_buf_impl;

// link the stream up to the black hole buffer
struct NullStream : public std::ostream {
    NullStream () : std::ostream (&null_buf_impl) {}
};

// create a singleton and point the reference to it
static NullStream null_stream_impl;
std::ostream& Opm::null_stream = null_stream_impl;
