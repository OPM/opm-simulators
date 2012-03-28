/* 
 * File:   WellsGroup.cpp
 * Author: kjetilo
 * 
 * Created on March 27, 2012, 9:27 AM
 */

#include "WellsGroup.hpp"
namespace Opm {
AbstractWellsGroup::AbstractWellsGroup() {
}

AbstractWellsGroup::~AbstractWellsGroup() {
}

const std::string& AbstractWellsGroup::get_name() {
    return name_;
}
}
}