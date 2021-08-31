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
 * @brief   Kernel functions for prati scheme
 * @author  C. Tsigginos
 * @details Kernel functions for Prati scheme
 */


#ifndef PRATI_KERNEL_INC_
#define PRATI_KERNEL_INC_

#include "ops_lib_core.h"

#ifdef OPS_MPI
#include "ops_mpi_core.h"
#endif
#include "type.h"
#include "model.h"
#include "flowfield_host_device.h"
#include "model_host_device.h"
#include "force_fsi_host_device.h"

void KerInitializePrati(ACC<Real>& Fd, const int* size) {

	int nElems = (*size);
	for (int iDir = 0; iDir < nElems; iDir++) {
#ifdef OPS_3D
		Fd(iDir, 0, 0, 0) = 0.0;
#endif

#ifdef OPS_2D
		Fd(iDir, 0, 0) =0.0;
#endif
	}

}

#ifdef OPS_3D
void KerPratiUpdateU3D(ACC<Real>&U, const ACC<int>& nodeType, const ACC<int>& idp,
		const ACC<Real>& sfm, const ACC<Real>& vp, const ACC<Real>& xf, const Real* dt, const Real* tau,
		const Real* gammaM, const int* typeForce, const Real* forceInput, const int* sizeForce,
		const int* noelems, const int* spacedim) {

	int SPACEDIM = (*spacedim);
	int NELEM = (*noelems);
	Real poros = 1.0;
	Real sfEff = 0.0;
	Real eff;
	int noTot = 0;
	ForceType forceFlag =(ForceType)static_cast<int>(*typeForce);
	int forceSize = (*sizeForce);
	Real force[SPACEDIM], xfl[SPACEDIM];
	VertexType vt = (VertexType) nodeType(0, 0, 0);
	if (vt != VertexType::ImmersedSolid) {

		for (int iPar = 0 ; iPar < NELEM; iPar++)  {
			poros -= sfm(iPar, 0, 0, 0);
			noTot+= (idp(iPar, 0, 0, 0) > -1);
		}


		if (forceFlag != noForce) {
			for (int iDir = 0; iDir < SPACEDIM; iDir++)
				xfl[iDir] = xf(iDir, 0, 0, 0);

			calculateForce(force, forceInput, xfl, forceFlag, forceSize, SPACEDIM);
			U(0, 0, 0) += 0.5 * poros * (*dt) * force[0]; //OK
		}






		for (int iPar = 0; iPar < noTot; iPar++) {
			eff = (*tau) /((*tau) + (*gammaM) * (*dt) * (1 - sfm(iPar, 0, 0, 0)));
			U(0, 0, 0) += 0.5 * sfm(iPar, 0, 0, 0) * eff * vp(SPACEDIM * iPar, 0, 0, 0);
			sfEff += eff * sfm(iPar, 0 , 0, 0);
		}

#ifdef CPU
		if (poros < 1.0) {
#endif
			U(0, 0, 0) /=(1.0 + 0.5 * sfEff);

#ifdef CPU
		}
#endif
	}


}

void KerPratiUpdateV3D(ACC<Real>&V, const ACC<int>& nodeType, const ACC<int>& idp,
		const ACC<Real>& sfm, const ACC<Real>& vp, const ACC<Real>& xf, const Real* dt, const Real* tau,
		const Real* gammaM, const int* typeForce, const Real* forceInput, const int* sizeForce,
		const int* noelems, const int* spacedim) {

	int SPACEDIM = (*spacedim);
	int NELEM = (*noelems);
	Real poros = 1.0;
	Real eff;
	Real sfEff = 0.0;
	int noTot = 0;
	ForceType forceFlag =(ForceType)static_cast<int>(*typeForce);
	int forceSize = (*sizeForce);
	Real force[SPACEDIM], xfl[SPACEDIM];

	VertexType vt = (VertexType) nodeType(0, 0, 0);
	if (vt != VertexType::ImmersedSolid) {

		for (int iPar = 0 ; iPar < NELEM; iPar++)  {
			poros -= sfm(iPar, 0, 0, 0);
			noTot+= (idp(iPar, 0, 0, 0) > -1);
		}

		if (forceFlag != noForce) {
			for (int iDir = 0; iDir < SPACEDIM; iDir++)
				xfl[iDir] = xf(iDir, 0, 0, 0);

			calculateForce(force, forceInput, xfl, forceFlag, forceSize, SPACEDIM);

			V(0, 0, 0) += 0.5 * poros * (*dt) * force[1]; //OK-We may want to add another function
		}

		for (int iPar = 0; iPar < noTot; iPar++) {
			eff = (*tau) /((*tau) + (*gammaM) * (*dt) * (1 - sfm(iPar, 0, 0, 0)));
			sfEff += eff * sfm(iPar, 0 , 0, 0);
			V(0, 0, 0) += 0.5 * sfm(iPar, 0, 0, 0) * eff * vp(SPACEDIM * iPar + 1, 0, 0, 0);
		}

#ifdef CPU
		if (poros < 1.0) {
#endif
			Real sf = 1.0 - poros;
			V(0, 0, 0) /=(1.0 + 0.5 * sfEff);

#ifdef CPU
		}
#endif
	}


}

void KerPratiUpdateW3D(ACC<Real>&W, const ACC<int>& nodeType, const ACC<int>& idp,
		const ACC<Real>& sfm, const ACC<Real>& vp, const ACC<Real>& xf, const Real* dt, const Real* tau,
		const Real* gammaM, const int* typeForce, const Real* forceInput, const int* sizeForce,
		const int* noelems, const int* spacedim) {

	int SPACEDIM = (*spacedim);
	int NELEM = (*noelems);
	Real poros = 1.0;
	Real sfEff = 0.0;
	Real eff;
	int noTot = 0;
	ForceType forceFlag =(ForceType)static_cast<int>(*typeForce);
	int forceSize = (*sizeForce);
	Real force[SPACEDIM], xfl[SPACEDIM];

	VertexType vt = (VertexType) nodeType(0, 0, 0);
	if (vt != VertexType::ImmersedSolid) {

		for (int iPar = 0 ; iPar < NELEM; iPar++)  {
			poros -= sfm(iPar, 0, 0, 0);
			noTot+= (idp(iPar, 0, 0, 0) > -1);
		}

		if (forceFlag != noForce) {
			for (int iDir = 0; iDir < SPACEDIM; iDir++)
				xfl[iDir] = xf(iDir, 0, 0, 0);

			calculateForce(force, forceInput, xfl, forceFlag, forceSize, SPACEDIM);

			W(0, 0, 0) += 0.5 * poros * (*dt) * force[2];

		}




		for (int iPar = 0; iPar < noTot; iPar++) {
			eff = (*tau) /((*tau) + (*gammaM) * (*dt) * (1 - sfm(iPar, 0, 0, 0)));
			sfEff += eff * sfm(iPar, 0 , 0, 0);
			W(0, 0, 0) += 0.5 * sfm(iPar, 0, 0, 0) * eff * vp(SPACEDIM * iPar + 2, 0, 0, 0);

		}

#ifdef CPU
		if (poros < 1.0) {
#endif
			Real sf = 1.0 - poros;
			W(0, 0, 0) /=(1.0 + 0.5 * sfEff);
#ifdef CPU
		}
#endif
	}



}

void KerCollisionPrati3D(ACC<Real>& fcopy, ACC<Real>& Fd, const ACC<Real>& f,
		const ACC<Real>& coordinates, const ACC<Real>& nodeType,
		const ACC<Real>& Rho, const ACC<Real>& U, const ACC<Real>& V,
		const ACC<Real>& W,const ACC<int> &id,const ACC<Real>& sfp,
		const ACC<Real>& vp, const Real* dt, const Real* tauRef,
		const int* typeForce, const Real* forceVariables, const int* sizeForce,
		const Real* gamma, const int* nElem, const int* spacedim, const int* lattIdx) {

	int NELEM = (*nElem);

	VertexType vt = (VertexType)nodeType(0, 0, 0);
	bool collision = (vt != VertexType::ImmersedSolid);
	if (collision) {

		int spaceDim = (*spacedim);
		int NELEM = (*nElem);
		Real sumx[NELEM], sumy[NELEM], sumz[NELEM];
		Real dtOverTauPlusDt, dtOverTauF;
		Real OmegaS, OmegaS1, OmegaS2, OmegaS3, OmegaS4, OmegaS5, OmegaS6;
		Real force[spaceDim], xfl[spaceDim];
		Real alpha, sf;
		Real poros = 1.0;
		Real eff, eff1, porosEff;
		ForceType forceFlag =(ForceType)static_cast<int>(*typeForce);
		int forceSize = (*sizeForce);

		Real u[spaceDim];
		int oppos, oppos1, noTot1;
		int noTot = 0;
		Real rho{Rho(0, 0, 0)};
		u[0] = U(0, 0, 0);
		u[1] = V(0, 0, 0);
		u[2] = W(0, 0, 0);
		const Real T{1};
		Real tau = (*tauRef);

		for (int pIndex = 0; pIndex < NELEM; pIndex++) {
			sumx[pIndex] = 0.0;
			sumy[pIndex] = 0.0;
			sumz[pIndex] = 0.0;
			noTot+= (id(pIndex, 0, 0, 0) > -1);
			poros-=sfp(pIndex, 0, 0, 0);

		}

		noTot1 = noTot;
		if (noTot == 0)
			noTot1 = 1;

		//find porosEffective
		porosEff = 1.0;
		for (int iPar = 0; iPar < noTot; iPar++)
			porosEff -= sfp(iPar, 0, 0, 0) * tau / (tau + (*gamma) * (*dt) * ( 1- sfp(iPar, 0, 0, 0)));
		sf = 1.0 - porosEff;

		dtOverTauPlusDt = porosEff * (*dt) / (tau + 0.5 * porosEff * (*dt));
		dtOverTauF = poros * (*dt) * tau / (tau+ 0.5 * porosEff * (*dt));
		alpha = tau/(tau + 0.5 * porosEff * (*dt));
		int numXi = lattIdx[1] - lattIdx[0] + 1;

		Real up[spaceDim * noTot1];
		Real feq[numXi], feqS[numXi * noTot1], bodyForce[numXi];

		for (int pIndex  = 0; pIndex < noTot; pIndex++) {
			up[spaceDim * pIndex] =  vp(pIndex * spaceDim, 0, 0, 0);
			up[spaceDim * pIndex + 1] = vp(pIndex * spaceDim + 1, 0, 0, 0);
			up[spaceDim * pIndex + 2] = vp(pIndex * spaceDim + 2, 0, 0, 0);
		}

		if (forceFlag != noForce) {
			for (int iDir = 0; iDir < spaceDim; iDir++)
				xfl[iDir] = coordinates(iDir, 0, 0, 0);

			calculateForce(force, forceVariables, xfl, forceFlag, forceSize, spaceDim);

		}


		//Calculate BodyForce and equilibrium functions
		for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
			feq[xiIndex] = CalcBGKFeq(xiIndex, rho, u[0], u[1], u[2], T, 2);
			if (forceFlag != noForce)
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

				OmegaS1 =  alpha * (f(oppos1, 0, 0, 0) - feq[oppos]);  //OK
				OmegaS2 = alpha * (feq[xiIndex]- f(xiIndex1, 0, 0, 0)); //OK
				OmegaS4 = 0.5 * alpha * poros * (*dt) * (bodyForce[oppos]-bodyForce[xiIndex]); //OK

				OmegaS5 = 0.0;
				for (int iPar = 0; iPar < noTot; iPar++) {
					eff = tau /(tau + (*gamma) * (*dt) * (1.0 - sfp(iPar, 0, 0, 0)));
					OmegaS5+= 0.5 * sfp(iPar, 0, 0, 0) * eff *
							alpha * (feqS[iPar * numXi + oppos]- feq[oppos]
							+ feqS[iPar * numXi + xiIndex] - feq[xiIndex]);
				}

				for (int iPar = 0; iPar < noTot; iPar++) {
					eff = tau /(tau + (*gamma) * (*dt) * (1.0 - sfp(iPar, 0, 0, 0)));
					OmegaS3 = feqS[iPar * numXi + xiIndex] - feq[xiIndex];
					OmegaS6 = 0.0;
					if (noTot>1) {
						for (int jPar = 0; jPar < noTot; jPar++) {
							if (iPar != jPar) {
								eff1 = tau /(tau + (*gamma) * (*dt) * (1.0 - sfp(jPar, 0, 0, 0)));
								OmegaS6 +=1.0 * sfp(jPar, 0, 0, 0)* eff1 * alpha *
										(feqS[iPar*NUMXI+xiIndex] - feqS[jPar * NUMXI + xiIndex]);
							}
						}
					}
					OmegaS += (OmegaS1 + OmegaS2 + OmegaS3 + OmegaS4 + OmegaS5 + OmegaS6) *
							sfp(iPar, 0, 0, 0) * eff * alpha/ (1 + (sf)*alpha);
					sumx[iPar] = sumx[iPar] - (OmegaS1 + OmegaS2 + OmegaS3 + OmegaS4 +
							OmegaS5 + OmegaS6) * XI[xiIndex1 * LATTDIM] * CS;
					sumy[iPar] = sumy[iPar] - (OmegaS1 + OmegaS2 + OmegaS3 + OmegaS4 +
							OmegaS5 + OmegaS6) * XI[xiIndex1 * LATTDIM + 1] * CS;
					sumz[iPar] = sumz[iPar] - (OmegaS1 + OmegaS2 + OmegaS3 + OmegaS4 +
							OmegaS5 + OmegaS6) * XI[xiIndex1 * LATTDIM + 2] * CS;
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

void KerDragPrati3D(const ACC<int>& id, const ACC<Real>& sfp, const ACC<Real>& xf,
		const ACC<Real>& Fd, Real* FDp, Real* TDp, const int* idp, const Real* xp,
		const Real* dt, const Real* tau, const Real* gamma,
		const int* spacedim, const int* nelem) {

	int SPACEDIM = (*spacedim);
	int NELEM = (*nelem);

	Real alpha, xfl[SPACEDIM];
	int loc;
	int noTot = 0;
	Real poros = 1.0;
	Real Beta, sf;
	Real porosEff, eff;



	for (int pIndex = 0; pIndex < NELEM; pIndex++) {
		noTot += (id(pIndex, 0, 0, 0) > -1);
		poros -= sfp(pIndex, 0, 0, 0);
	}

	porosEff = 1.0;
	for (int pIndex = 0; pIndex < noTot; pIndex++) {
		eff = (*tau) / ((*tau) + (*gamma) * (*dt) * (1.0 - sfp(pIndex, 0, 0, 0)));
		porosEff -= sfp(pIndex, 0, 0, 0) * eff;
	}

	sf = 1.0 - porosEff;

	alpha = (*tau)/((*tau) + 0.5 * porosEff * (*dt));


	for (int pIndex = 0; pIndex < noTot; pIndex++) {
		if ((*idp) == id(pIndex, 0, 0, 0)) {
			eff = (*tau) / ((*tau) + (*gamma) * (*dt) * (1.0 - sfp(pIndex, 0, 0, 0)));
			Beta = sfp(pIndex, 0, 0, 0) * eff /(1.0 + (sf) * alpha);
			loc = pIndex * SPACEDIM;

			xfl[0] = xf(loc, 0, 0, 0);
			xfl[1] = xf(loc + 1, 0, 0, 0);
			xfl[2] = xf(loc + 2, 0, 0, 0);
			FDp[0] += Beta * Fd(loc, 0, 0, 0);
			FDp[1] += Beta * Fd(loc + 1, 0, 0, 0);
			FDp[2] += Beta * Fd(loc + 2, 0, 0, 0);
			TDp[0] += Beta * ((xfl[1] - xp[1]) * Fd(loc+2, 0, 0, 0) - (xfl[2] - xp[2]) * Fd(loc+1,0,0, 0));
			TDp[1] += Beta * ((xfl[2] - xp[2]) * Fd(loc, 0, 0, 0) - (xfl[0] - xp[0]) * Fd(loc+2, 0, 0, 0));
			TDp[2] += Beta * ((xfl[0] - xp[0]) * Fd(loc+1, 0, 0, 0) - (xfl[1] - xp[1]) * Fd(loc, 0, 0, 0));

//	ops_printf("rank %d: FDp = [%e %e %e]\n", ops_get_proc(), FDp[0], FDp[1], FDp[2]);
			break;
		}
	}


}
#endif
//2D Functions
#ifdef OPS_2D
void KerPratiUpdateU2D(ACC<Real>&U, const ACC<int>& nodeType, const ACC<int>& idp,
		const ACC<Real>& sfm, const ACC<Real>& vp, const ACC<Real>& xf, const Real* dt, const Real* tau,
		const Real* gammaM, const int* typeForce, const Real* forceInput, const int* sizeForce,
		const int* noelems, const int* spacedim) {

	int SPACEDIM = (*spacedim);
	int NELEM = (*noelems);
	Real poros = 1.0;
	Real sfEff = 0.0;
	Real eff;
	int noTot = 0;
	ForceType forceFlag =(ForceType)static_cast<int>(*typeForce);
	int forceSize = (*sizeForce);
	Real force[SPACEDIM], xfl[SPACEDIM];
	VertexType vt = (VertexType) nodeType(0, 0);
	if (vt != VertexType::ImmersedSolid) {

		for (int iPar = 0 ; iPar < NELEM; iPar++)  {
			poros -= sfm(iPar, 0, 0);
			noTot+= (idp(iPar, 0, 0) > -1);
		}


		if (forceFlag != noForce) {
			for (int iDir = 0; iDir < SPACEDIM; iDir++)
				xfl[iDir] = xf(iDir, 0, 0);

			calculateForce(force, forceInput, xfl, forceFlag, forceSize, SPACEDIM);
			U(0, 0) += 0.5 * poros * (*dt) * force[0]; //OK
		}






		for (int iPar = 0; iPar < noTot; iPar++) {
			eff = (*tau) /((*tau) + (*gammaM) * (*dt) * (1 - sfm(iPar, 0, 0)));
			U(0, 0) += 0.5 * sfm(iPar, 0, 0) * eff * vp(SPACEDIM * iPar, 0, 0);
			sfEff += eff * sfm(iPar, 0 , 0);
		}

#ifdef CPU
		if (poros < 1.0) {
#endif
			U(0, 0) /=(1.0 + 0.5 * sfEff);

#ifdef CPU
		}
#endif
	}


}

void KerPratiUpdateV2D(ACC<Real>&V, const ACC<int>& nodeType, const ACC<int>& idp,
		const ACC<Real>& sfm, const ACC<Real>& vp, const ACC<Real>& xf, const Real* dt, const Real* tau,
		const Real* gammaM, const int* typeForce, const Real* forceInput, const int* sizeForce,
		const int* noelems, const int* spacedim) {

	int SPACEDIM = (*spacedim);
	int NELEM = (*noelems);
	Real poros = 1.0;
	Real eff;
	Real sfEff = 0.0;
	int noTot = 0;
	ForceType forceFlag =(ForceType)static_cast<int>(*typeForce);
	int forceSize = (*sizeForce);
	Real force[SPACEDIM], xfl[SPACEDIM];

	VertexType vt = (VertexType) nodeType(0, 0);
	if (vt != VertexType::ImmersedSolid) {

		for (int iPar = 0 ; iPar < NELEM; iPar++)  {
			poros -= sfm(iPar, 0, 0);
			noTot+= (idp(iPar, 0, 0) > -1);
		}

		if (forceFlag != noForce) {
			for (int iDir = 0; iDir < SPACEDIM; iDir++)
				xfl[iDir] = xf(iDir, 0, 0);

			calculateForce(force, forceInput, xfl, forceFlag, forceSize, SPACEDIM);

			V(0, 0) += 0.5 * poros * (*dt) * force[1]; //OK-We may want to add another function
		}

		for (int iPar = 0; iPar < noTot; iPar++) {
			eff = (*tau) /((*tau) + (*gammaM) * (*dt) * (1 - sfm(iPar, 0, 0)));
			sfEff += eff * sfm(iPar, 0 , 0);
			V(0, 0) += 0.5 * sfm(iPar, 0, 0) * eff * vp(SPACEDIM * iPar + 1, 0, 0);
		}

#ifdef CPU
		if (poros < 1.0) {
#endif
			Real sf = 1.0 - poros;
			V(0, 0) /=(1.0 + 0.5 * sfEff);

#ifdef CPU
		}
#endif
	}


}


void KerCollisionPrati(ACC<Real>& fcopy, ACC<Real>& Fd, const ACC<Real>& f,
		const ACC<Real>& coordinates, const ACC<Real>& nodeType,
		const ACC<Real>& Rho, const ACC<Real>& U, const ACC<Real>& V,
		const ACC<int> &id,const ACC<Real>& sfp,
		const ACC<Real>& vp, const Real* dt, const Real* tauRef,
		const int* typeForce, const Real* forceVariables, const int* sizeForce,
		const Real* gamma, const int* nElem, const int* spacedim, const int* lattIdx) {

	int NELEM = (*nElem);

	VertexType vt = (VertexType)nodeType(0, 0);
	bool collision = (vt != VertexType::ImmersedSolid);
	if (collision) {

		int spaceDim = (*spacedim);
		int NELEM = (*nElem);
		Real sumx[NELEM], sumy[NELEM];
		Real dtOverTauPlusDt, dtOverTauF;
		Real OmegaS, OmegaS1, OmegaS2, OmegaS3, OmegaS4, OmegaS5, OmegaS6;
		Real force[spaceDim], xfl[spaceDim];
		Real alpha, sf;
		Real poros = 1.0;
		Real eff, eff1, porosEff;
		ForceType forceFlag =(ForceType)static_cast<int>(*typeForce);
		int forceSize = (*sizeForce);

		Real u[spaceDim];
		int oppos, oppos1, noTot1;
		int noTot = 0;
		Real rho{Rho(0, 0)};
		u[0] = U(0, 0);
		u[1] = V(0, 0);

		const Real T{1};
		Real tau = (*tauRef);

		for (int pIndex = 0; pIndex < NELEM; pIndex++) {
			sumx[pIndex] = 0.0;
			sumy[pIndex] = 0.0;
			noTot+= (id(pIndex, 0, 0) > -1);
			poros-=sfp(pIndex, 0, 0);

		}

		noTot1 = noTot;
		if (noTot == 0)
			noTot1 = 1;

		//find porosEffective
		porosEff = 1.0;
		for (int iPar = 0; iPar < noTot; iPar++)
			porosEff -= sfp(iPar, 0, 0) * tau / (tau + (*gamma) * (*dt) * ( 1- sfp(iPar, 0, 0)));
		sf = 1.0 - porosEff;

		dtOverTauPlusDt = porosEff * (*dt) / (tau + 0.5 * porosEff * (*dt));
		dtOverTauF = poros * (*dt) * tau / (tau+ 0.5 * porosEff * (*dt));
		alpha = tau/(tau + 0.5 * porosEff * (*dt));
		int numXi = lattIdx[1] - lattIdx[0] + 1;

		Real up[spaceDim * noTot1];
		Real feq[numXi], feqS[numXi * noTot1], bodyForce[numXi];

		for (int pIndex  = 0; pIndex < noTot; pIndex++) {
			up[spaceDim * pIndex] =  vp(pIndex * spaceDim, 0, 0);
			up[spaceDim * pIndex + 1] = vp(pIndex * spaceDim + 1, 0, 0);
		}

		if (forceFlag != noForce) {
			for (int iDir = 0; iDir < spaceDim; iDir++)
				xfl[iDir] = coordinates(iDir, 0, 0);

			calculateForce(force, forceVariables, xfl, forceFlag, forceSize, spaceDim);

		}


		//Calculate BodyForce and equilibrium functions
		for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
			feq[xiIndex] = CalcBGKFeq(xiIndex, rho, u[0], u[1], T, 2);
			if (forceFlag != noForce)
				bodyForce[xiIndex] = CalcBodyForce2ndOrder(xiIndex, rho, u[0], u[1], force);
			else
				bodyForce[xiIndex] = 0.0;

			for (int iPar = 0; iPar < noTot; iPar++) {
				if (sfp(iPar, 0, 0) > 0.0)
					feqS[iPar * numXi + xiIndex] = CalcBGKFeq(xiIndex, rho, up[spaceDim * iPar],
							up[spaceDim * iPar + 1], T, 2);
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

				OmegaS1 =  alpha * (f(oppos1, 0, 0) - feq[oppos]);  //OK
				OmegaS2 = alpha * (feq[xiIndex]- f(xiIndex1, 0, 0)); //OK
				OmegaS4 = 0.5 * alpha * poros * (*dt) * (bodyForce[oppos]-bodyForce[xiIndex]); //OK

				OmegaS5 = 0.0;
				for (int iPar = 0; iPar < noTot; iPar++) {
					eff = tau /(tau + (*gamma) * (*dt) * (1.0 - sfp(iPar, 0, 0)));
					OmegaS5+= 0.5 * sfp(iPar, 0, 0) * eff *
							alpha * (feqS[iPar * numXi + oppos]- feq[oppos]
							+ feqS[iPar * numXi + xiIndex] - feq[xiIndex]);
				}

				for (int iPar = 0; iPar < noTot; iPar++) {
					eff = tau /(tau + (*gamma) * (*dt) * (1.0 - sfp(iPar, 0, 0)));
					OmegaS3 = feqS[iPar * numXi + xiIndex] - feq[xiIndex];
					OmegaS6 = 0.0;
					if (noTot>1) {
						for (int jPar = 0; jPar < noTot; jPar++) {
							if (iPar != jPar) {
								eff1 = tau /(tau + (*gamma) * (*dt) * (1.0 - sfp(jPar, 0, 0)));
								OmegaS6 +=1.0 * sfp(jPar, 0, 0)* eff1 * alpha *
										(feqS[iPar*NUMXI+xiIndex] - feqS[jPar * NUMXI + xiIndex]);
							}
						}
					}
					OmegaS += (OmegaS1 + OmegaS2 + OmegaS3 + OmegaS4 + OmegaS5 + OmegaS6) *
							sfp(iPar, 0, 0) * eff * alpha/ (1 + (sf)*alpha);
					sumx[iPar] = sumx[iPar] - (OmegaS1 + OmegaS2 + OmegaS3 + OmegaS4 +
							OmegaS5 + OmegaS6) * XI[xiIndex1 * LATTDIM] * CS;
					sumy[iPar] = sumy[iPar] - (OmegaS1 + OmegaS2 + OmegaS3 + OmegaS4 +
							OmegaS5 + OmegaS6) * XI[xiIndex1 * LATTDIM + 1] * CS;
				}
			}
			fcopy(xiIndex1, 0, 0) = f(xiIndex1, 0, 0) - dtOverTauPlusDt
					* (f(xiIndex1, 0, 0)-feq[xiIndex]) + dtOverTauF * bodyForce[xiIndex] + OmegaS;
#ifdef CPU
			const Real res{fcopy(xiIndex, 0, 0)};
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
			Fd(spaceDim * iPar, 0, 0) = sumx[iPar];
			Fd(spaceDim * iPar + 1, 0, 0) = sumy[iPar];


		}
	}

}

void KerDragPrati(const ACC<int>& id, const ACC<Real>& sfp, const ACC<Real>& xf,
		const ACC<Real>& Fd, Real* FDp, Real* TDp, const int* idp, const Real* xp,
		const Real* dt, const Real* tau, const Real* gamma,
		const int* spacedim, const int* nelem) {

	int SPACEDIM = (*spacedim);
	int NELEM = (*nelem);

	Real alpha, xfl[SPACEDIM];
	int loc;
	int noTot = 0;
	Real poros = 1.0;
	Real Beta, sf;
	Real porosEff, eff;



	for (int pIndex = 0; pIndex < NELEM; pIndex++) {
		noTot += (id(pIndex, 0, 0) > -1);
		poros -= sfp(pIndex, 0, 0);
	}

	porosEff = 1.0;
	for (int pIndex = 0; pIndex < noTot; pIndex++) {
		eff = (*tau) / ((*tau) + (*gamma) * (*dt) * (1.0 - sfp(pIndex, 0, 0)));
		porosEff -= sfp(pIndex, 0, 0) * eff;
	}

	sf = 1.0 - porosEff;

	alpha = (*tau)/((*tau) + 0.5 * porosEff * (*dt));


	for (int pIndex = 0; pIndex < noTot; pIndex++) {
		if ((*idp) == id(pIndex, 0, 0)) {
			eff = (*tau) / ((*tau) + (*gamma) * (*dt) * (1.0 - sfp(pIndex, 0, 0)));
			Beta = sfp(pIndex, 0, 0) * eff /(1.0 + (sf) * alpha);
			loc = pIndex * SPACEDIM;

			xfl[0] = xf(loc, 0, 0);
			xfl[1] = xf(loc + 1, 0, 0);

			FDp[0] += Beta * Fd(loc, 0, 0);
			FDp[1] += Beta * Fd(loc + 1, 0, 0);

			TDp[0] += Beta * ((xfl[0] - xp[0]) * Fd(loc+1, 0, 0) - (xfl[1] - xp[1]) * Fd(loc, 0, 0));

//	ops_printf("rank %d: FDp = [%e %e %e]\n", ops_get_proc(), FDp[0], FDp[1], FDp[2]);
			break;
		}
	}


}
#endif
#endif /* APPS_LBM_DEM_NOP_PRATI_KERNEL_INC_ */
