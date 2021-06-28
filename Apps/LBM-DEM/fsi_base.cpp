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

/*! @brief Base class for fluid-particle interaction modeling
 *  @author C. Tsigginos
 *  Discription: Base class for handling the fluid-particle interactions
 **/

#include "fsi_base.h"
#include <cmath>
FsiBase::FsiBase(Component& compoUser, int spacedim, Real* forceUser, bool owned,
		SolFracType porosUser, Real gammaUser) : compo{compoUser}{

	gamma = gammaUser;
	spaceDim = spacedim;
	porosModel =  porosUser;

	force = new Real[spaceDim];
	Real sumF = 0.0;
	for (int iDim = 0; iDim < spaceDim; iDim++) {
		force[iDim] = forceUser[iDim];
		sumF += abs(force[iDim]);
	}

	forceFlag = 0;
	if (sumF != 0.0)
		forceFlag = 1;
	collisionOwned = owned;

}

FsiBase::~FsiBase() {

	delete[] force;
}

void FsiBase::ObtainID(int* idVel, int* loop, Real& tauCompo, CollisionType& collisModel,
		int& idComponent,int& rhoId,int &thId) {

	idVel[0] = compo.uId;
	idVel[1] = compo.vId;

#if OPS_3D
	idVel[2] = compo.wId;
#endif

	loop[0] = compo.index[0];
	loop[1] = compo.index[1];
	tauCompo = compo.tauRef;
	collisModel = compo.collisionType;
	idComponent = compo.id;

	rhoId = compo.macroVars.at(Variable_Rho).id;

	if (isThermalModel == 1)
		thId = compo.macroVars.at(Variable_T).id;
	else
		thId = -1;

}
