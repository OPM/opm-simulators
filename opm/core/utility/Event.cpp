#include <opm/core/utility/Event.hpp>

using namespace std;
using namespace Opm;

Event&
EventSource::add (std::function <void ()> handler) {
    // add handler to the back of the queue
    handlers_.push_back (handler);

    // return ourselves so we can be used in a call chain
    return *this;
}

void
EventSource::signal () {
    // loop through the list of handlers, and invoke every one of them
    // (range-based for loops are not available until GCC 4.6)
    for (auto it = handlers_.begin(); it != handlers_.end(); ++it) {
        (*it) ();
    }
}
