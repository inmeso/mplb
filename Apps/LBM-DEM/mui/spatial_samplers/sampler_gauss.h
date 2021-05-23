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
*****************************************************************************/

/**
 * @file sampler_gauss.h
 * @author Y. H. Tang
 * @date 10 February 2014
 * @brief Spatial sampler that provides a value at a point using Gaussian
 * interpolation.
 */

#ifndef MUI_SAMPLER_GAUSS_H_
#define MUI_SAMPLER_GAUSS_H_

#include "../util.h"
#include "../config.h"
#include "../sampler.h"

namespace mui {

template<typename CONFIG=default_config, typename O_TP=typename CONFIG::REAL, typename I_TP=O_TP>
class sampler_gauss {
public:
	using OTYPE      = O_TP;
	using ITYPE      = I_TP;
	using REAL       = typename CONFIG::REAL;
	using INT        = typename CONFIG::INT;
	using point_type = typename CONFIG::point_type;

	sampler_gauss( REAL r_, REAL h_ ) : r(r_), h(h_), nh(std::pow(2*PI*h,-0.5*CONFIG::D)) {}

	template<template<typename,typename> class CONTAINER>
	inline OTYPE filter( point_type focus, const CONTAINER<ITYPE,CONFIG> &data_points ) const {
		REAL  wsum = 0;
		OTYPE vsum = 0;
		for( size_t i = 0 ; i < data_points.size() ; i++ ) {
			auto d = normsq( focus - data_points[i].first );
			if ( d < r*r ) {
			    REAL w = nh * std::exp( (REAL(-0.5)/h) * d );
				vsum += data_points[i].second * w;
				wsum += w;
			}
		}
		if ( wsum ) return vsum / wsum;
		else return REAL(0.);
	}

	inline geometry::any_shape<CONFIG> support( point_type focus, REAL domain_mag ) const {
		return geometry::sphere<CONFIG>( focus, r );
	}

protected:
	REAL r;
	REAL h;
	REAL nh;
};

}

#endif /* MUI_SAMPLER_GAUSS_H_ */
