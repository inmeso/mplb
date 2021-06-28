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

/*! @brief  Class for calculating the solid fraction for spherical particles
 * @author C. Tsigginos
 */

#ifndef POROS_SPHERICAL_H_
#define POROS_SPHERICAL_H_

#include "particle_to_grid_base.h"
#include "type.h"
#include "flowfield_host_device.h"
#include <cmath>

class PorosSpherical : public ParticleToGridBase {
	public:
	PorosSpherical(int ParticleType,int spaceDim);
	virtual int  particleShape() {return 1;}
	virtual void DefineVariables(int noElem, SizeType timestep = 0);
	virtual void ParticleProjection();
	virtual void UpdateProjection();
	virtual void MappingFunction(bool flag);
	virtual void InitializeVariables();
	virtual void ReturnParticleShape() { ops_printf("Handles spherical particles only"); };
	virtual void PrintMappingVariables();

	protected:
	static void KerSolidFracSphere(ACC<int>& id, ACC<Real>& sfp, ACC<Real>& vp,
			ACC<Real>& xAvg, const ACC<Real>& xf, const Real* xPos, const Real* Radius,
			const int* idP, const Real* vPart, const Real* omPart, const Real* dx,
			const int* nelem, const int* spacedim);

	static void KerSolidVelocitySphere(ACC<Real>& vP, const ACC<int>& id,
			const ACC<Real>& xAvg, const int* idParticle,
			const Real* xPos, const Real* radPart, const Real* velP, const Real* omP,
			const Real* dx, const Real* spacedim, const Real* noelem);

	static void KerInitializePorousSpherical(ACC<Real>&sfP, ACC<Real>& vP, ACC<Real>& xAvg,
				ACC<int>& id, const int* spacedim, const int* noelem);

	static inline OPS_FUN_PREFIX Real CalculateSolidFractionSpheres(const Real* xfl,
			const Real Ravg, const Real* xPos, const Real Rp, Real* xAvg);

	static void KerPrintPorousData(const ACC<Real>&sfP,const ACC<Real>& vP,
			const ACC<Real>& xAvg, const ACC<int>& id, const ACC<Real>& xf,
			const int* spacedim, const int* noelem);

	private:
		int sizeofData = 4;
};


#endif /* APPS_LBM_DEM_POROS_SPHERICAL_H_ */
