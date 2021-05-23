/*****************************************************************************
* Multiscale Universal Interface Code Coupling Library                       *
*                                                                            *
* Copyright (C) 2019 Y. H. Tang, S. Kudo, X. Bian, Z. Li, G. E. Karniadakis  *
*                                                                            *
* This software is jointly licensed under the Apache License, Version 2.0    *
* and the GNU General Public License version 3, you may use it according     *
* to either.                                                                 *
*                                                                            *
* ** Apache License, version 2.0 **                                          *
*                                                                            *
* Licensed under the Apache License, Version 2.0 (the "License");            *
* you may not use this file except in compliance with the License.           *
* You may obtain a copy of the License at                                    *
*                                                                            *
* http://www.apache.org/licenses/LICENSE-2.0                                 *
*                                                                            *
* Unless required by applicable law or agreed to in writing, software        *
* distributed under the License is distributed on an "AS IS" BASIS,          *
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
* See the License for the specific language governing permissions and        *
* limitations under the License.                                             *
*                                                                            *
* ** GNU General Public License, version 3 **                                *
*                                                                            *
* This program is free software: you can redistribute it and/or modify       *
* it under the terms of the GNU General Public License as published by       *
* the Free Software Foundation, either version 3 of the License, or          *
* (at your option) any later version.                                        *
*                                                                            *
* This program is distributed in the hope that it will be useful,            *
* but WITHOUT ANY WARRANTY; without even the implied warranty of             *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
* GNU General Public License for more details.                               *
*                                                                            *
* You should have received a copy of the GNU General Public License          *
* along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
******************************************************************************/

/**
 * @file comm_factory.h
 * @author Y. H. Tang
 * @date 14 March 2014
 * @brief Structures and methods to create a new communicator based on chosen
 * protocols.
 */

#ifndef COMM_FACTORY_H_
#define COMM_FACTORY_H_

#include "util.h"
#include "exception.h"
#include "lib_uri.h"
#include "lib_dispatcher.h"
#include "lib_singleton.h"
#include "comm.h"

namespace mui {

struct comm_factory: public singleton<dispatcher<std::string, std::function<communicator *(const char [])> > >
{
	static communicator *create_comm( const char URI[] ) {
		if ( !instance().exist(uri(URI).protocol()) ) {
			exception_segv( "MUI Error [comm_factory.h]: Unknown communicator type ", uri(URI).protocol() );
		}
		return instance()[uri(URI).protocol()]( URI );
	}
};

}

#endif /* COMM_FACTORY_H_ */
