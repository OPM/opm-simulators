#ifndef OPM_EVENT_HEADER_INCLUDED
#define OPM_EVENT_HEADER_INCLUDED

// Copyright (C) 2013 Uni Research AS
// This file is licensed under the GNU General Public License v3.0

#include <functional>
#include <list>

namespace Opm {

/// Interface to register interest in receiving notifications when a
/// certain event, such as the completion of a timestep, has happened.
struct Event {
    /// Register a callback to receive notifications from this event.
    ///
    /// \param[in] handler
    /// Function object that will be invoked when the event happens.
    ///
    /// \return
    /// The event object itself, so that multiple additions can be chained.
    ///
    /// \note
    /// The event may happen several times, and the handler will receive
    /// a notification every time.
    ///
    /// \note
    /// If a handler is added more than once, it will also be called
    /// more than once.
    virtual Event& add (std::function <void ()> handler) = 0;

    /// Convenience routine to add a member function of a class as
    /// an event handler.
    ///
    /// This allows us to have all the necessary information the handler
    /// needs put into an object, and then register this with the event.
    template <typename T, void (T::*member)()> Event& add (T& t);
};

/// Generator of event notifications.
///
/// As more than one event is possible from an object, it is expected
/// that event servers implements this functionality by aggregation and
/// provide accessors to let clients reach the various events.
///
/// You should not provide the full EventSource interface to clients,
/// as this will allow them to signal the event themselves; rather, return
/// the registration-only parent interface.
///
/// \example
/// You can add an event to your code like this:
///
/// \code{.cpp}
/// struct Foo {
///   // accessor of the event that other can register at
///   Event& completed () { return completed_; }
///
///   // something that ultimately triggers the event
///   void action () { /* ... */ completed_.signal(); }
///
/// private:
///   EventSource completed_;
/// };
/// \endcode
///
/// It could then be accessed by the client like this:
///
/// \code{.cpp}
/// struct Bar {
///   void callMe() { /* ... */ }
/// };
/// \endcode
///
/// \code{.cpp}
///   Foo foo;
///   Bar bar;
///
///   // setup the connection between the two
///   foo.completed().add<Bar, &Bar::callMe>(bar);
///
///   // set events in motion
///   foo.action();
/// \endcode
class EventSource : public Event {
public:
    virtual Event& add (std::function <void ()> handler);
    virtual void signal ();
protected:
    /// List of actual handlers that will be called
    std::list <std::function <void ()> > handlers_;
};

// inline definitions
#include "Event_impl.hpp"

} /* namespace Opm */

#endif /* OPM_EVENT_HEADER_INCLUDED */
