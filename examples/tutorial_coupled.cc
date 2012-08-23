// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009-2012 by Andreas Lauser                               *
 *   Copyright (C) 2010-2012 by Melanie Darcis                               *
 *   Copyright (C) 2011 by Benjamin Faigle                                   *
 *   Copyright (C) 2009-2012 by Bernd Flemisch                               *
 *   Copyright (C) 2010 by Klaus Mosthaf                                     *
 *   Copyright (C) 2010 by Philipp Nuske                                     *
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
 * \brief Main file of the tutorial for a fully coupled twophase box model.
 */
#include "config.h" /*@\label{tutorial-coupled:include-begin}@*/
#include "tutorialproblem_coupled.hh"  /*@\label{tutorial-coupled:include-problem-header}@*/
#include <dumux/common/start.hh> /*@\label{tutorial-coupled:include-end}@*/

int main(int argc, char** argv)
{
    typedef TTAG(TutorialProblemCoupled) TypeTag; /*@\label{tutorial-coupled:set-type-tag}@*/
    return Dumux::start<TypeTag>(argc, argv); /*@\label{tutorial-coupled:call-start}@*/
}
