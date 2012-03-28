/* 
 * File:   WellsGroup.cpp
 * Author: kjetilo
 * 
 * Created on March 27, 2012, 9:27 AM
 */

#include "WellsGroup.hpp"
namespace Opm {
WellsGroupInterface::AbstractWellsGroup() {
}

WellsGroupInterface::~WellsGroupInterface() {
}

const std::string& WellsGroupInterface::name() {
    return name_;
}
}
}