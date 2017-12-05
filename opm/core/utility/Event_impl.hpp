// Copyright (C) 2013 Uni Research AS
// This file is licensed under the GNU General Public License v3.0

#ifndef OPM_EVENT_HEADER_INCLUDED
#error Do NOT include this file directly!
#endif /* OPM_EVENT_HEADER_INCLUDED */

template <typename T, void (T::*member)()> inline Event&
Event::add (T& t) {
    // wrap the member function in a std::function and add that
    // notice the use of ref() to avoid invoking the copy constructor
    return this->add (std::bind (member, std::ref(t)));
}
