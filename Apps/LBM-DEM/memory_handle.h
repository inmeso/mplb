/**
 * Copyright 2019 United Kingdom Research and Innovation
 *
 * Authors: See AUTHORS
 *
 * Contact: [jianping.meng@stfc.ac.uk and/or jpmeng@gmail.com]
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice
 *    this list of conditions and the following disclaimer in the documentation
 *    and or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * ANDANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

/*!
 * @brief   Function for dynamic allocation and re-allocation of particle related lists
 * @author  Chrysovalantis Tsigginos
 * @details Functions for dynamic allocation of particle related lists
 */

#ifndef MEMORY_HANDLE_H_
#define MEMORY_HANDLE_H_

#include <cmath>
#include <string>
#include "type.h"
#include "flowfield.h"
#include <stdlib.h>
#include <stdio.h>

template <class T>
inline void free2d(T **&a) {

	for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {
		free(a[blockIndex]);
	}

	free(a);
}

template<class T>
inline void Reallocate(T *&a, int size, int bm) {

	a= (T *) realloc(a, size * bm * sizeof(T));

}
template<class T>
inline void allocateMemory(T *&a, int size, int bm) {
	a = (T *) malloc(size * bm * sizeof(T));
}
#endif /* APPS_LBM_DEM_MEMORY_HANDLE_H_ */
