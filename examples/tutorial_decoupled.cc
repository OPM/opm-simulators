// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009-2010 by Bernd Flemisch                               *
 *   Copyright (C) 2010-2011 by Benjamin Faigle                              *
 *   Copyright (C) 2008-2012 by Markus Wolff                                 *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief tutorial for the sequential two-phase model
 */
#include "config.h" /*@\label{tutorial-decoupled:include-begin}@*/

#include "tutorialproblem_decoupled.hh" /*@\label{tutorial-decoupled:include-problem-header}@*/
#include <ewoms/common/start.hh> /*@\label{tutorial-decoupled:include-end}@*/

////////////////////////
// the main function
////////////////////////
int main(int argc, char** argv)
{
    typedef TTAG(TutorialProblemDecoupled) TypeTag; /*@\label{tutorial-decoupled:set-type-tag}@*/
    return Ewoms::start<TypeTag>(argc, argv); /*@\label{tutorial-decoupled:call-start}@*/
}
