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

/*! @brief Base class for the implementation of the PRATI scheme
 *  @author C. Tsigginos
 **/
#ifndef PRATI_H_
#define PRATI_H_

#include "fsi_base.h"
#include "field.h"
#include "block.h"
#include "type.h"
#include <string>
#include "model.h"
#include "block_particles.h"
#include "mapping_models.h"

class Prati : public FsiBase {

	protected:
		RealField Fd;
		ParticleToGridBase* poros; //TODO: Add the porosity model
		int noElem = 1;

	public:
	Prati(Component componentUser, int spacedim, Real* forceUser, bool owned = false,
			SolFracType porosModel = Mode_None, Real gammaUser = 0.0, int nelem = 2, int ParticleType = 1);
	~Prati();
	virtual void ModelCollision(); //Dpme
	virtual void CalculateDragForce();//Done
	virtual void PostVelocityCalculation(); //Done
	virtual void MappingFunction(bool flag); //Done
	virtual void DefineVariables(SizeType timestep = 0); //DONE
    virtual void InitializeVariables(); //DONE
    virtual void WriteToHdf5(const std::string& caseName, const SizeType timeStep);

	static void KerCollisionPrati3D(ACC<Real>& fcopy, ACC<Real>& Fd, const ACC<Real> f,
			const ACC<Real>& coordinates, const ACC<Real>& nodeType,
			const ACC<Real>& rho, const ACC<Real>& U, const ACC<Real>& V,
			const ACC<Real>& W,const ACC<int> &id,const ACC<Real>& sfp,
			const ACC<Real>& vp, const Real* dt, const Real* tau,
			const int* forceFlag, const Real* force, const Real* gamma,
			const int* nElem, const int* spacedim, const int* lattIdx);

	static  void KerPratiUpdateU3D(ACC<Real>&U, const ACC<int>& nodeType, const ACC<int>& idp,
			const ACC<Real>& sfm, const ACC<Real>& vp, const Real* dt, const Real* tau,
			const Real* gammaM, const int* forceflag, const Real* force, const int* noelems,
			const int* spacedim);

	static void  KerPratiUpdateV3D(ACC<Real>&V, const ACC<int>& nodeType, const ACC<int>& idp,
			const ACC<Real>& sfm, const ACC<Real>& vp, const Real* dt, const Real* tau,
			const Real* gammaM, const int* forceflag, const Real* force, const int* noelems,
			const int* spacedim);

	static void KerPratiUpdateW3D(ACC<Real>&W, const ACC<int>& nodeType, const ACC<int>& idp,
			const ACC<Real>& sfm, const ACC<Real>& vp, const Real* dt, const Real* tau,
			const Real* gammaM, const int* forceflag, const Real* force, const int* noelems,
			const int* spacedim);

	static void KerDragPRATI(const ACC<int>& id, const ACC<Real>& sfp, const ACC<Real>& xf,
			const ACC<Real>& Fd, Real* FDp, Real* TDp, const int* idp, const Real* xp,
			const Real* dt, const Real* tau, const Real* gamma,
			const int* spacedim, const int* nelem);

	static void KerInitialize(ACC<Real>& Fd, const int* size);

	//static void KerAuxiliary3D(ACC<Real>& uf);

	/*ToDo add the rest*/
};




#endif /* APPS_LBM_DEM_PRATI_H_ */
