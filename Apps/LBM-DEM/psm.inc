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

/*! @brief   Implement kernel function for Partially saturated method
 * @author C. Tsigginos
*/

#ifndef PSM_INC_
#define PSM_INC_

#include "psm.h"

void Psm::KerCollisionPSM3D(ACC<Real>& fcopy, ACC<Real>& Fd, const ACC<Real>& f,
		const ACC<Real>& coordinates, const ACC<Real>& nodeType,
		const ACC<Real>& Rho, const ACC<Real>& U, const ACC<Real>& V,
		const ACC<Real>& W,const ACC<int> &id,const ACC<Real>& sfp,
		const ACC<Real>& vp, const Real* dt, const Real* tauRef,
		const int* forceFlag, const Real* force,
		const int* nElem, const int* spacedim, const int* lattIdx) {

	int NELEM = (*nElem);

	VertexType vt = (VertexType)nodeType(0, 0, 0);
	bool collision = (vt != VertexType::ImmersedSolid);
	if (collision) {
		int spaceDim = (*spacedim);
		int NELEM = (*nElem);
		Real sumx[NELEM], sumy[NELEM], sumz[NELEM];
		Real dtOverTauPlusDt, dtOverTauF;
		Real OmegaS, OmegaS1, OmegaS2;
		Real alpha, sf;
		Real poros = 1.0;
		Real eff, porosEff;

		Real u[spaceDim], rho;
		int oppos, oppos1, noTot1;
		int noTot = 0;

		rho= Rho(0, 0, 0);
		u[0] = U(0, 0, 0);
		u[1] = V(0, 0, 0);
		u[2] = W(0, 0, 0);
		const Real T{1};
		Real tau = (*tauRef);

		porosEff = 1.0;
		for (int pIndex = 0; pIndex < NELEM; pIndex++) {
			sumx[pIndex] = 0.0;
			sumy[pIndex] = 0.0;
			sumz[pIndex] = 0.0;
			noTot+= (id(pIndex, 0, 0, 0) > -1);
			poros -=sfp(pIndex, 0, 0, 0);
		}


		for (int pIndex = 0; pIndex < noTot; pIndex++)
			porosEff -= sfp(pIndex, 0, 0) * tau / (tau - (1.0 - sfp(pIndex, 0, 0, 0)) * (*dt));

		noTot1 = noTot;
		if (noTot == 0)
			noTot1 = 1;



		dtOverTauPlusDt = porosEff * (*dt) / (tau + 0.5  * (*dt));
		dtOverTauF =  porosEff *  (*dt);
		alpha = tau/(tau + 0.5 * porosEff * (*dt));
		int numXi = lattIdx[1] - lattIdx[0] + 1;
		Real up[spaceDim * noTot1];
		Real feq[numXi], feqS[numXi * noTot1], bodyForce[numXi];

		for (int pIndex  = 0; pIndex < noTot; pIndex++) {
			up[spaceDim * pIndex] =  vp(pIndex * spaceDim, 0, 0, 0);
			up[spaceDim * pIndex + 1] = vp(pIndex * spaceDim + 1, 0, 0, 0);
			up[spaceDim * pIndex + 2] = vp(pIndex * spaceDim + 2, 0, 0, 0);
		}

		//Calculate BodyForce and equilibrium functions
		for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
			feq[xiIndex] = CalcBGKFeq(xiIndex, rho, u[0], u[1], u[2], T, 2);
			if ((*forceFlag) == 1)
				bodyForce[xiIndex] = CalcBodyForce2ndOrder(xiIndex, rho, u[0], u[1], u[2], force); //TODO ADD it to the model_kernel.inc
			else
				bodyForce[xiIndex] = 0.0;

			for (int iPar = 0; iPar < noTot; iPar++) {
				if (sfp(iPar, 0, 0, 0) > 0.0)
					feqS[iPar * numXi + xiIndex] = CalcBGKFeq(xiIndex, rho, up[spaceDim * iPar],
							up[spaceDim * iPar + 1], up[spaceDim * iPar + 2], T, 2);
				else
					feqS[iPar * numXi + xiIndex] = 0.0;
			}
		}

		for (int xiIndex1 = lattIdx[0]; xiIndex1 <= lattIdx[1]; xiIndex1++) {
			OmegaS = 0.0;
			int xiIndex = xiIndex1 - lattIdx[0];
			if (noTot > 0) {
				oppos1 = OPP[xiIndex1];
				oppos =  oppos1 - lattIdx[0];
				OmegaS1 =   (f(oppos1, 0, 0, 0) - feq[oppos]);
				OmegaS2 = 0.0;
				for (int iPar = 0; iPar < noTot; iPar++) {
					eff = tau / (tau + (1.0 - sfp(iPar, 0, 0, 0)) * (*dt));
					OmegaS2 = (feqS[iPar * numXi + xiIndex] - f(xiIndex1, 0, 0, 0));
					OmegaS += sfp(iPar, 0, 0) * eff * (OmegaS1 + OmegaS2);
					sumx[iPar] -= (OmegaS1 + OmegaS2) * XI[xiIndex1 * LATTDIM] * CS;
					sumy[iPar] -= (OmegaS1 + OmegaS2) * XI[xiIndex1 * LATTDIM + 1] * CS;
					sumz[iPar] -= (OmegaS1 + OmegaS2) * XI[xiIndex1 * LATTDIM + 2] * CS;
				}
			}
			fcopy(xiIndex1, 0, 0, 0) = f(xiIndex1, 0, 0, 0) - dtOverTauPlusDt
				* (f(xiIndex1, 0, 0, 0)-feq[xiIndex]) + dtOverTauF * bodyForce[xiIndex] + OmegaS;



#ifdef CPU
			const Real res{fcopy(xiIndex, 0, 0, 0)};
			if (std::isnan(res) || res < 0 || std::isinf(res)) {
				ops_printf(
					"Error-PRATI: Distribution function %f becomes "
					"invalid for the component %i at  the lattice "
					"%i\n",
					res, 0, xiIndex);
				assert(!(std::isnan(res) || res < 0 || std::isinf(res)));
			}
#endif  // CPU


		}


		for (int iPar = 0; iPar < noTot; iPar++) {
			Fd(spaceDim * iPar, 0, 0, 0) = sumx[iPar];
			Fd(spaceDim * iPar + 1, 0, 0, 0) = sumy[iPar];
			Fd(spaceDim * iPar + 2, 0, 0, 0) = sumz[iPar];

		}
	}



}

void Psm::KerInitialize(ACC<Real>& Fd, const int* size) {

	int total = (*size);
#ifdef OPS_3D
	for (int iPart = 0; iPart < total; iPart++)
		Fd(iPart, 0, 0, 0) = 0.0;
#endif

#ifdef OPS_2D
	for (int iPart = 0; iPart < total; iPart++)
		Fd(iPart, 0, 0) = 0.0;

#endif
}

void Psm::KerDragPSM(const ACC<int>& id, const ACC<Real>& sfp, const ACC<Real>& xf,
						const ACC<Real>& Fd, Real* FDp, Real* TDp, const int* idp, const Real* xp,
						const Real* dt, const Real* tau,
						const int* spacedim, const int* nelem) {

	int SPACEDIM = (*spacedim);
	int NELEM = (*nelem);

	Real alpha, xfl[SPACEDIM];
	int loc;
	int noTot = 0;
	Real poros = 1.0;
	Real Beta, sf;
	Real porosEff, eff;

	xfl[0] = xf(0, 0, 0, 0);
	xfl[1] = xf(1, 0, 0, 0);
	xfl[2] = xf(2, 0, 0, 0);

	for (int pIndex = 0; pIndex < NELEM; pIndex++) {
		noTot += (id(pIndex, 0, 0, 0) > -1);
		poros -= sfp(pIndex, 0, 0, 0);
	}

	for (int pIndex = 0; pIndex < noTot; pIndex++) {
		if ((*idp) == id(pIndex, 0, 0, 0)) {
			eff = (*tau) / ((*tau) + (1.0 - sfp(pIndex, 0, 0, 0)) * (*dt));
			Beta = sfp(pIndex, 0, 0, 0) * eff;
			loc = pIndex * SPACEDIM;

			xfl[0] = xf(loc , 0, 0, 0);
			xfl[1] = xf(loc + 1, 0, 0, 0);
			xfl[2] = xf(loc + 2, 0, 0, 0);


			FDp[0] += Beta * Fd(loc, 0, 0, 0);
			FDp[1] += Beta * Fd(loc + 1, 0, 0, 0);
			FDp[2] += Beta * Fd(loc + 2, 0, 0, 0);
			TDp[0] += Beta * ((xfl[1] - xp[1]) * Fd(loc+2, 0, 0, 0) - (xfl[2] - xp[2]) * Fd(loc+1,0,0, 0));
			TDp[1] += Beta * ((xfl[2] - xp[2]) * Fd(loc, 0, 0, 0) - (xfl[0] - xp[0]) * Fd(loc+2, 0, 0, 0));
			TDp[2] += Beta * ((xfl[0] - xp[0]) * Fd(loc+1, 0, 0, 0) - (xfl[1] - xp[1]) * Fd(loc, 0, 0, 0));

			break;
		}
	}


}




#endif /* APPS_LBM_DEM_PSM_INC_ */
