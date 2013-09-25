#ifndef OPM_NULLSTREAM_HEADER_INCLUDED
#define OPM_NULLSTREAM_HEADER_INCLUDED

// Copyright (C) 2013 Uni Research AS
// This file is licensed under the GNU General Public License v3.0

#include <iosfwd>

namespace Opm {

/**
 * Output stream that ignores everything written to it.
 *
 * Use this stream if you want to disable output without having a
 * lot of conditionals in your code.
 *
 * Since the null stream has no state, there is no point in
 * instantiating your own; simply use this reference instead.
 *
 * @example
 * @code{.cpp}
 *   std::ostream& outp = (quiet ? Opm::null_stream : std::cerr);
 *   outp << "Hello, World!" << std::endl;
 * @endcode
 */
extern std::ostream& null_stream;

} /* namespace Opm */

#endif /* OPM_NULLSTREAM_HEADER_INCLUDED */
