/*
 * psm_kernel.h
 *
 *  Created on: 13 Jan 2020
 *      Author: valantis
 */

#ifndef PSM_KERNEL_H_
#define PSM_KERNEL_H_

#include "psm.h"
#include  <stdio.h>
#include <math.h>
void KerSetInitialMacroVars(ACC<Real>& macros, const ACC<Real>& coords, const Real* initialMacros) {

#ifdef OPS_2D
	macros(0, 0, 0) = initialMacros[0];
	macros(1, 0, 0) = initialMacros[1];
	macros(2, 0, 0) = initialMacros[2];
#endif


#ifdef OPS_3D
	macros(0, 0, 0, 0) = initialMacros[0];
	macros(1, 0, 0, 0) = initialMacros[1];
	macros(2, 0, 0, 0) = initialMacros[2];
	macros(3, 0, 0, 0) = initialMacros[3];
#endif


}

void KerInitialize(ACC<Real>& sfp, ACC<Real>& vp, ACC<int>& id, ACC<Real>& Fd) {
	int sizet = NELEM * SPACEDIM;

#ifdef OPS_2D
	for (int elem = 0; elem < NELEM; elem++) {
		sfp(elem, 0, 0) = 0.0;
		id(elem, 0, 0) = -1;
	}

	for (int elem = 0; elem < sizet; elem++) {
		vp(elem, 0, 0) = 0.0;
		Fd(elem, 0, 0) = 0.0;
	}
#endif

#ifdef OPS_3D
	for (int elem = 0; elem < NELEM; elem++) {
			sfp(elem, 0, 0, 0) = 0.0;
			id(elem, 0, 0, 0) = -1;
	}

	for (int elem = 0; elem < sizet; elem++) {
		vp(elem, 0, 0, 0) = 0.0;
		Fd(elem, 0, 0, 0) = 0.0;
	}
#endif

}




void KerCalcMacroVarsForceFluidSolid(ACC<Real>& macroVars, const ACC<int>& nodeType, const ACC<Real>& sfp, const Real* dt, const Real* force) {
#ifdef OPS_2D
	VertexType vt = (VertexType) nodeType[OPS_ACC1(0,0)];
	if (vt !=  VertexType::ImmersedSolid) {
		Real poros = 1.0;
		for (int iPar = 0; iPar < NELEM; iPar++)
			poros -= sfp(iPar, 0, 0);


		for (int m = 0; m < SPAisnanCEDIM; m++)
			macroVars(m+1, 0, 0) += 0.5 * poros * force[m] * (*dt);

	}
#endif

#ifdef OPS_3D
	VertexType vt = (VertexType) nodeType(0, 0, 0);

	if (vt != VertexType::ImmersedSolid) {
		Real poros = 1.0;
		for (int iPar = 0; iPar < NELEM; iPar++)
			poros -= sfp(iPar, 0, 0, 0);

		for (int m = 0; m < LATTDIM; m++)
			macroVars(m+1, 0, 0, 0) += 0.5 * poros * force[m] * (*dt);

	}
#endif

}

//2D functions
#ifdef OPS_2D
void KerSolVelocityUpdateMP2D(ACC<Real>& vp, const ACC<Real>& xf, const ACC<int>& id, const int* idp, const Real* xp, const Real* up,const Real* tau) {
	  double xfl[2];
	  double vf[2];
	  int noTot = 0;


	  xfl[0]= xf(0, 0, 0);
	  xfl[1] = xf(1, 0, 0);

	  for (int pIndex = 0; pIndex < NELEM; pIndex++)
		  noTot += id(pIndex, 0, 0>-1);

	  if (noTot > 0) {
		  for (int pIndex = 0; pIndex < noTot; pIndex++) {
			  if ((*idp)==id(pIndex, 0, 0)) {
				  EstimateVelocity2D(vf, up, xp, xfl);
				  vp(SPACEDIM * pIndex, 0, 0) = vf[0];
				  vp(SPACEDIM * pIndex+1, 0, 0) = vf[1];
				  break;
			  }
		  }


	  }


}


void KerEvaluateSolidFractionMP2D(ACC<Real>& sfp, ACC<Real>& vp, ACC<int>& id, const ACC<Real>& xf, const int* idp, const Real* xp,
								  const Real* velP, const Real* dx, const int* intOrder) {

	Real xfl[SPACEDIM];
	Real xe, ye, dd;
	int nEdges, noIndex;
	int tot = 0;
	Real vf[SPACEDIM];
	int s[16] = {1,0, 0, 1, -1,0, 0, -1, 1, 1, -1, 1,-1,-1,1,-1}; //control points
	nEdges = 8;
	int no = 0;

	xfl[0] = xf(0, 0, 0);
	xfl[1] = xf(1, 0, 0);

	for (int edgeIndex = 0; edgeIndex < nEdges; edgeIndex++) {
		xe = xfl[0] + 0.5 * (*dx) * static_cast<double>(s[2*edgeIndex]);
		ye = xfl[1] + 0.5 * (*dx) * static_cast<double>(s[2*edgeIndex+1]);
		dd = (xp[0] - xe) * (xp[0] - xe) + (xp[1] - ye) * (xp[1] - ye);
		if (dd <= xp[2] * xp[2])
			tot+=1;
	}

	if (tot == nEdges) {
		no += 1;
		noIndex = no-1;
		id(noIndex,0,0) = *idp;
		sfp(noIndex, 0, 0)= 1.0;
		vp(SPACEDIM * noIndex, 0, 0) = velP[0] - (xfl[1] - xp[1]) * velP[2];
		vp(SPACEDIM * noIndex + 1, 0, 0) = velP[1] + (xfl[0] - xp[0]) * velP[2];
	}
	else if (tot > 0) { //Partially saturated node
		no += 1;
		noIndex = no-1;
		id(noIndex,0,0) = *idp;

		sfp(noIndex, 0, 0) =  ComputeCellSF(xfl, xp, (*dx), *intOrder);

		EstimateVelocity2D(vf, velP, xp, xfl);
		vp(2 * noIndex, 0, 0) = vf[0];
		vp(2 * noIndex + 1, 0, 0) = vf[1];
	}

}



void KerPRE2D(ACC<Real>& fcopy, ACC<Real>& Fd,const ACC<Real>& f, const ACC<Real>& macros,const ACC<int>& nodeType,
					 const ACC<int>& id, const ACC<Real>& sfp, const Real* tau, const Real* dt, const int* forceFlag, const Real* force) {

	VertexType vt = (VertexType) nodeType(0, 0);
	bool collisionRequired = (vt != VertexType::ImmersedSolid);

	if (collisionRequired) {
		Real dtOverTauPlusDt, dtOverTauF;
		Real sumx[NELEM], sumy[NELEM];
		Real OmegaS, OmegaS1, OmegaS2, OmegaS3;
		Real  sf;
		Real poros = 1.0;
		Real u[SPACEDIM], rho;
		int oppos;
		int noTot = 0;

		for (int pIndex = 0; pIndex < NELEM; pIndex++) {
			sumx[pIndex] = 0.0;
			sumy[pIndex] = 0.0;
			noTot+= (id(pIndex, 0, 0) > -1);
			poros-=sf(pIndex, 0, 0);
		}

		sf = 1.0 - poros;


		int noTot1;
		if (noTot > 0)
			noTot1 = noTot;
		else
			noTot1 = 1;

		dtOverTauPlusDt = poros * (*dt) / ((*tau) + 0.5 * poros * (*dt));
		dtOverTauF = dtOverTauPlusDt * (*tau);
		invDtOverTauF = (*tau ) / ((*tau) + 0.5 * poros * (*dt))


		Real up[SPACEDIM * noTot1];
		Real feq[NUMXI], feqS[NUMXI * noTot1], bodyForce[NUMXI];

		for (int pIndex = 0; pIndex < noTot; pIndex++) {
			up[SPACEDIM * pIndex] = vp(pIndex * SPACEDIM, 0, 0);
			up[SPACEDIM * pIndex + 1] = vp(pIndex * SPACEDIM +1, 0, 0);

		}

		rho = macros(0, 0, 0);
		u[0] = macros(1, 0, 0);
		u[1] = macros(2, 0, 0);


		//Calculate BodyForces and equilibrium functions
		for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
			feq[xiIndex] = CalcBGKFeq(xiIndex, rho, u[0], u[1], 1, 2);
			//printf("xiIndex = %d: feq = %12.9e,u = %12.9e v= %12.9e\n",xiIndex, feq[xiIndex], u[0], u[1]);
			if ((*forceFlag) == 1)
				bodyForce[xiIndex] = CalcBodyForce2ndOrder(xiIndex, rho, u, force);
			else
				bodyForce[xiIndex] = 0.0;

			for (int iPar = 0; iPar < noTot; iPar++)
				feqS[iPar * NUMXI + xiIndex] = CalcBGKFeq(xiIndex, rho, up[SPACEDIM * iPar], up[SPACEDIM * iPar + 1], 1, 2);

		}
		
		for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
			OmegaS = 0.0;
			if (noTot > 0) {
				oppos= OPP[xiIndex];
				OmegaS1 =  (f(oppos, 0, 0) - feq[oppos]);


				OmegaS3 = 0.5 * porosEff * (*dt) * (bodyForce[oppos]-bodyForce[xiIndex]);


				for (int iPar = 0; iPar < noTot; iPar++) {


					OmegaS2 = feqS[iPar * NUXMI xiIndex] - f(xiIndex, 0, 0);


					OmegaS += (OmegaS1 + OmegaS2 + OmegaS3) * sfp(iPar, 0, 0);

					sumx[iPar] = sumx[iPar] - (OmegaS1 + OmegaS2  + OmegaS3) * XI[xiIndex * LATTDIM] * CS;
					sumy[iPar] = sumy[iPar] - (OmegaS1 + OmegaS2  + OmegaS3) * XI[xiIndex * LATTDIM + 1] * CS;
				}

			}
			//OmegaS = 0.0;
			fcopy(xiIndex, 0, 0) = f(xiIndex, 0, 0) - dtOverTauPlusDt * (f(xiIndex, 0, 0)-feq[xiIndex])
								   + dtOverTauF * bodyForce[xiIndex] + OmegaS * invDtOverTauF;

		}

		for (int iPar = 0; iPar < noTot; iPar++) {
			Fd(2*iPar, 0, 0) = sumx[iPar];
			Fd(2*iPar+1, 0, 0) = sumy[iPar];
		}
	}
}

void KerDragForcePRE2D(const ACC<int>&id, const ACC<Real>& sfp, const ACC<Real>&xf, const ACC<Real>&Fd, const int* idp, const Real* xp,
					 Real* FDp, const Real* dt, const Real* tau) {

	Real alpha;
	int loc;
	int noTot = 0;
	Real poros = 1.0;
	Real Beta, sf;




	for (int pIndex = 0; pIndex < NELEM; pIndex++)
		noTot += (id(pIndex, 0, 0) > -1);

	for (int pIndex = 0; pIndex < noTot; pIndex++)
		poros -= sfp(pIndex, 0, 0);



	sf = 1.0 - poros;

	for (int pIndex = 0; pIndex < noTot; pIndex++) {
		if ((*idp) == id(pIndex, 0, 0)) {
			Beta = sfp(pIndex, 0, 0) * (*tau) / ((*tau) + 0.5 * poros * (*dt));
			loc = pIndex * SPACEDIM;
			FDp[0] += Beta * Fd(loc, 0, 0);
			FDp[1] += Beta * Fd(loc + 1, 0, 0);
			FDp[2] +=Beta * ((xf(0,0,0)-xp[0]) * Fd(loc+1,0,0) -
							 (xf(1,0,0) - xp[1]) * Fd(loc,0,0));

			break;
		}
	}


}
#endif


//3D kernels


void KerEvaluateSolidFractionMP3D(ACC<int>& id, ACC<Real>& sfp, ACC<Real>& vp, const ACC<Real>& xf, const int* idp, const Real* xp,
								  const Real* velP, const Real* dx) {

	Real vf[SPACEDIM] , xfl[SPACEDIM];
	Real xe, ye, ze;
	Real dd;
	Real dx1, dy1, dz1;
	int ip, jp, kp;
	int status, nEdges, noIndex;
	int itot = 0;
	int sx[26] = { 1,  1,  1,  1, -1, -1, -1, -1,  1, -1,  0,  0,  0,  0,  1, -1, 1, -1,  0,  0,  0,  0,  1,  1, -1, -1};
	int sy[26] = { 1,  1, -1, -1,  1,  1, -1, -1,  0,  0,  1, -1,  0,  0,  1,  1,-1, -1,  1, -1,  1, -1,  0,  0,  0,  0};
	int sz[26] = { 1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0,  1, -1,  0,  0, 0,  0,  1,  1, -1, -1,  1, -1,  1, -1};


	int Npoints = 9;
	nEdges = 26;
	dx1 = (*dx);
	int no = 0;

	xfl[0] = xf(0, 0, 0, 0);
	xfl[1] = xf(1, 0, 0, 0);
	xfl[2] = xf(2, 0, 0, 0);

	for (int edgeIndex = 0; edgeIndex < nEdges; edgeIndex++) {
		xe = xfl[0] + 0.5 * (*dx) * static_cast<Real>(sx[edgeIndex]);
		ye = xfl[1] + 0.5 * (*dx) * static_cast<Real>(sy[edgeIndex]);
		ze = xfl[2] + 0.5 * (*dx) * static_cast<Real>(sz[edgeIndex]);

		dd = (xe - xp[0]) * (xe - xp[0]) + (ye - xp[1]) * (ye - xp[1]) + (ze - xp[2]) * (ze - xp[2]);

		if (dd <= xp[3] * xp[3])
			itot += 1;
	}

	if (itot == nEdges) {
		no += 1;
		noIndex = no - 1;
		id(noIndex, 0, 0, 0) = (*idp);
		sfp(noIndex, 0, 0, 0) = 1.0;
		dx1 = xfl[0] - xp[0];
		dy1 = xfl[1] - xp[1];
		dz1 = xfl[2] - xp[2];
		vp(SPACEDIM * noIndex, 0, 0, 0) = velP[0] + velP[4] * dz1 - velP[5] * dy1;
		vp(SPACEDIM * noIndex +1, 0, 0, 0) = velP[1] + velP[5] * dx1 - velP[3] * dz1;
		vp(SPACEDIM * noIndex +2, 0, 0, 0) = velP[2] + velP[3] * dy1 - velP[4] * dx1;
	}
	else if (itot > 0) {
		no += 1;
		noIndex = no - 1;
		id(noIndex, 0, 0, 0) = (*idp);
		sfp(noIndex, 0, 0, 0) = ComputeSolidFraction3D(xfl,xp, dx1, Npoints);
		EstimateVelocity3D(vf, velP, xp, xfl); //To be replaced by an average concept
		vp(SPACEDIM * noIndex, 0, 0, 0) = vf[0];
		vp(SPACEDIM * noIndex +1, 0, 0, 0) = vf[1];
		vp(SPACEDIM * noIndex +2, 0, 0, 0) = vf[2];
	}


#ifdef CPU
        const Real res{sfp(0, 0, 0, 0)};
        if (std::isnan(res) || res < 0 || res > 1.0) {
        ops_printf(
                "Error! Solid fraction %f not in [0, 1]\n");
        assert(!(std::isnan(res) || res < 0 || res > 1.0));
            }
#endif  // CPU

}


void KerSolVelocityUpdateMP3D(ACC<Real>& vp, const ACC<int>& id, const ACC<Real>& xf, const int* idp, const Real* xp, const Real* up) {

	Real xfl[3];
	Real vf[3];
	int noTot = 0;

	for (int iPar = 0; iPar < NELEM; iPar++) {
		noTot += (id(iPar, 0, 0, 0)>-1);
	}

	xfl[0] = xf(0, 0, 0, 0);
	xfl[1] = xf(1, 0, 0, 0);
	xfl[2] = xf(2, 0, 0, 0);

	if (noTot > 0) {
		for (int pIndex = 0; pIndex < noTot; pIndex++) {
			if ((*idp)==id(pIndex, 0, 0, 0)) {
				EstimateVelocity3D(vf, up, xp, xfl);
				vp(SPACEDIM * pIndex, 0, 0, 0) = vf[0];
				vp(SPACEDIM * pIndex + 1, 0, 0, 0) = vf[1];
				vp(SPACEDIM * pIndex + 2, 0, 0, 0) = vf[2];
				break;
			}
		}
	}

}



void  KerPRE3D(ACC<Real>& fcopy, ACC<Real>& Fd,const ACC<Real>& f, const ACC<Real>& macros,const ACC<int>& nodeType,
			   const ACC<int>& id, const ACC<Real>& sfp, const ACC<Real>& vp, const Real* tau, const Real* dt,
			   const int* forceFlag, const Real* force) {

	VertexType vt = (VertexType) nodeType(0, 0, 0);
	bool collisionRequired = (vt != VertexType::ImmersedSolid);

	if (collisionRequired) {
		Real dtOverTauPlusDt, dtOverTauF, invDtOverTau;
		Real sumx[NELEM], sumy[NELEM], sumz[NELEM];
		Real OmegaS, OmegaS1, OmegaS2, OmegaS3;
		Real alpha, sf;
		Real poros = 1.0;
		Real eff, porosEff;


		Real u[SPACEDIM], rho;
		int oppos;
		int noTot = 0;



		for (int pIndex = 0; pIndex < NELEM; pIndex++) {
			sumx[pIndex] = 0.0;
			sumy[pIndex] = 0.0;
			sumz[pIndex] = 0.0;
			noTot+= (id(pIndex, 0, 0, 0) > -1);
			poros-=sfp(pIndex, 0, 0, 0);

		}

		int noTot1;
		if (noTot > 0)
			noTot1 = noTot;
		else
			noTot1 = 1;

		sf = 1.0 - poros;



		dtOverTauPlusDt = poros * (*dt) / ((*tau) + 0.5 * poros * (*dt));
		dtOverTauF = poros * (*dt) * (*tau) / ((*tau) + 0.5 * poros * (*dt));
		invDtOverTau = (*tau) / ((*tau) + 0.5 * poros * (*dt));


		//if (poros < 1)
		//	printf("Rank %d: sf = %f alpha = %f eff = %f\n", ops_get_proc(), 1-poros, alpha, eff);

		Real up[SPACEDIM * noTot1];
		Real feq[NUMXI], feqS[NUMXI * noTot1], bodyForce[NUMXI];

		for (int pIndex = 0; pIndex < noTot; pIndex++) {
			up[SPACEDIM * pIndex] = vp(pIndex * SPACEDIM, 0, 0, 0);
			up[SPACEDIM * pIndex + 1] = vp(pIndex * SPACEDIM + 1, 0, 0, 0);
			up[SPACEDIM * pIndex + 2] = vp(pIndex * SPACEDIM + 2, 0, 0, 0);
		}

		rho = macros(0, 0, 0, 0);
		u[0] = macros(1, 0, 0, 0);
		u[1] = macros(2, 0, 0, 0);
		u[2] = macros(3, 0, 0, 0);


		//Calculate BodyForces and equilibrium functions
		for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
			feq[xiIndex] = CalcBGKFeq(xiIndex, rho, u[0], u[1], u[2], 1, 2);
			//printf("xiIndex = %d: feq = %12.9e,u = %12.9e v= %12.9e\n",xiIndex, feq[xiIndex], u[0], u[1]);
			if ((*forceFlag) == 1)
				bodyForce[xiIndex] = CalcBodyForce2ndOrder(xiIndex, rho, u, force);
			else
				bodyForce[xiIndex] = 0.0;

			for (int iPar = 0; iPar < noTot; iPar++)
				if (sfp(iPar, 0, 0, 0) > 0)
					feqS[iPar * NUMXI + xiIndex] = CalcBGKFeq(xiIndex, rho, up[SPACEDIM * iPar], up[SPACEDIM * iPar + 1],up[SPACEDIM * iPar + 2], 1, 2);
				else
					feqS[iPar * NUMXI + xiIndex] = 0.0;
		}

		for (int xiIndex = 0; xiIndex < NUMXI; xiIndex++) {
			OmegaS = 0.0;
			if (noTot > 0) {
				oppos= OPP[xiIndex];
				OmegaS1 = (f(oppos, 0, 0, 0) - feq[oppos]);
				OmegaS2 = 0.5 * poros * (*dt) * (bodyForce[oppos]-bodyForce[xiIndex]);


				for (int iPar = 0; iPar < noTot; iPar++) {

					OmegaS2 = feqS(iPar * NUMXI + xiIndex, 0, 0, 0) - f(xiIndex, 0, 0, 0);

					OmegaS += (OmegaS1 + OmegaS2 + OmegaS3) * sfp(iPar, 0, 0);

#if DebugLevel >= 3
					printf("Rank %d: xiIndex %d OmegaS = %f\n", ops_get_proc(), xiIndex, OmegaS);
#endif
					sumx[iPar] = sumx[iPar] - (OmegaS1 + OmegaS2 + OmegaS3) * XI[xiIndex * LATTDIM] * CS;
					sumy[iPar] = sumy[iPar] - (OmegaS1 + OmegaS2 + OmegaS3) * XI[xiIndex * LATTDIM + 1] * CS;
					sumz[iPar] = sumz[iPar] - (OmegaS1 + OmegaS2 + OmegaS3) * XI[xiIndex * LATTDIM + 2] * CS;
				}

			}
			//OmegaS = 0.0;
			fcopy(xiIndex, 0, 0, 0) = f(xiIndex, 0, 0, 0) - dtOverTauPlusDt * (f(xiIndex, 0, 0, 0)-feq[xiIndex])
								   + dtOverTauF * bodyForce[xiIndex] + invDtOverTau * OmegaS;


#ifdef CPU
            const Real res{fcopy(xiIndex, 0, 0, 0)};
            if (std::isnan(res) || res < 0 || std::isinf(res)) {
                ops_printf(
                    "Error-PSM: Distribution function %f becomes "
                    "invalid for the component %i at  the lattice "
                    "%i\n",
                    res, 0, xiIndex);
                assert(!(std::isnan(res) || res < 0 || std::isinf(res)));
            }
#endif  // CPU

		}

		for (int iPar = 0; iPar < noTot; iPar++) {
			Fd(SPACEDIM * iPar, 0, 0, 0) = sumx[iPar];
			Fd(SPACEDIM * iPar + 1, 0, 0, 0) = sumy[iPar];
			Fd(SPACEDIM * iPar + 2, 0, 0, 0) = sumz[iPar];
		}
	}
}

void KerDragPRE3D(const ACC<int>& id, const ACC<Real>& sfp, const ACC<Real>& xf, const ACC<Real>& Fd,
					   Real* FDp, const int* idp, const Real* xp, const Real* dt, const Real* tau, const Real* gamma) {

	Real alpha, xfl[SPACEDIM];
	int loc;
	int noTot = 0;
	Real poros = 1.0;
	Real Beta, sf;
	Real porosEff, eff;
	xfl[0] = xf(0, 0, 0, 0);
	xfl[1] = xf(1, 0, 0, 0);
	xfl[2] = xf(2, 0, 0, 0);

	for (int pIndex = 0; pIndex < NELEM; pIndex++)
		noTot += (id(pIndex, 0, 0, 0) > -1);


	for (int pIndex = 0; pIndex < noTot; pIndex++)
		poros -= sfp(pIndex, 0, 0, 0);
	sf = 1.0 - poros;



	for (int pIndex = 0; pIndex < noTot; pIndex++) {
		if ((*idp) == id(pIndex, 0, 0, 0)) {
			Beta = sfp(pIndex, 0, 0, 0) * (*tau) / ((*tau) + 0.5 * poros * (*dt));
			loc = pIndex * SPACEDIM;

			FDp[0] += Beta * Fd(loc, 0, 0, 0);
			FDp[1] += Beta * Fd(loc + 1, 0, 0, 0);
			FDp[2] += Beta * Fd(loc + 2, 0, 0, 0);
			FDp[3] += Beta * ((xfl[1] - xp[1]) * Fd(loc+2, 0, 0, 0) - (xfl[2] - xp[2]) * Fd(loc+1,0,0, 0));
			FDp[4] += Beta * ((xfl[2] - xp[2]) * Fd(loc, 0, 0, 0) - (xfl[0] - xp[0]) * Fd(loc+2, 0, 0, 0));
			FDp[5] += Beta * ((xfl[0] - xp[0]) * Fd(loc+1, 0, 0, 0) - (xfl[1] - xp[1]) * Fd(loc, 0, 0, 0));

		//	ops_printf("rank %d: FDp = [%e %e %e]\n", ops_get_proc(), FDp[0], FDp[1], FDp[2]);
			break;
		}
	}

}
void KerStreamPeriodic3D(ACC<Real> & f, const ACC<Real>& fStage,
                 const ACC<int>& nodeType, const ACC<int>& geometry) {
#ifdef OPS_3D
    VertexGeometryType vg = (VertexGeometryType)geometry(0, 0, 0);
    for (int compoIndex = 0; compoIndex < NUMCOMPONENTS; compoIndex++) {
        VertexType vt = (VertexType)nodeType(compoIndex, 0, 0, 0);
        for (int xiIndex = COMPOINDEX[2 * compoIndex];
             xiIndex <= COMPOINDEX[2 * compoIndex + 1]; xiIndex++) {
            int cx = (int)XI[xiIndex * LATTDIM];
            int cy = (int)XI[xiIndex * LATTDIM + 1];
            int cz = (int)XI[xiIndex * LATTDIM + 2];

            if (vt == VertexType::Fluid) {
                f(xiIndex, 0, 0, 0) = fStage(xiIndex, -cx, -cy, -cz);
            }

            if (vt == VertexType::Periodic)
            	f(xiIndex, 0, 0, 0) = fStage(xiIndex, -cx, -cy, -cz);

            if (vt != VertexType::ImmersedSolid && vt != VertexType::Fluid && vt!= VertexType::Periodic) {
                // TODO te be determined if necessary
                bool streamRequired{true};
                if (streamRequired) {
                    if ((cx == 0) && (cy == 0) && (cz == 0)) {
                        f(xiIndex, 0, 0, 0) = fStage(xiIndex, 0, 0, 0);
                        continue;
                    }
                }

                switch (vg) {
                        // faces six types
                    case VG_IP:
                        // (cx=0 means stream is implemented at i=0,so here we
                        //  disable the step at boundary)
                        if (streamRequired) {
                            if (cx <= 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        } else {
                            if (cx < 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        }
                        break;
                    case VG_IM:
                        if (streamRequired) {
                            if (cx >= 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        } else {
                            if (cx > 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        }
                        break;
                    case VG_JP:
                        if (streamRequired) {
                            if (cy <= 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        } else {
                            if (cy < 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        }
                        break;
                    case VG_JM:
                        if (streamRequired) {
                            if (cy >= 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        } else {
                            if (cy > 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        }
                        break;
                    case VG_KP:
                        // (cx=0 means stream is implemented at i=0,so here we
                        //  disable the step at boundary)
                        if (streamRequired) {
                            if (cz <= 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        } else {
                            if (cz < 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        }
                        break;
                    case VG_KM:
                        if (streamRequired) {
                            if (cz >= 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }

                        } else {
                            if (cz > 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        }
                        break;
                    // faces six types end
                    // 12 edges
                    case VG_IPJP_I:
                        if (streamRequired) {
                            if (cy <= 0 && cx <= 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        } else {
                            if (cy < 0 && cx < 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        }
                        break;
                    case VG_IPJM_I:
                        if (streamRequired) {
                            if (cy >= 0 && cx <= 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        } else {
                            if (cy > 0 && cx < 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        }
                        break;
                    case VG_IMJP_I:
                        if (streamRequired) {
                            if (cy <= 0 && cx >= 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        } else {
                            if (cy < 0 && cx > 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        }
                        break;
                    case VG_IMJM_I:
                        if (streamRequired) {
                            if (cy >= 0 && cx >= 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        } else {
                            if (cy > 0 && cx > 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        }
                        break;
                    // k
                    case VG_IPKP_I:
                        if (streamRequired) {
                            if (cz <= 0 && cx <= 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        } else {
                            if (cz < 0 && cx < 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        }
                        break;
                    case VG_IPKM_I:
                        if (streamRequired) {
                            if (cz >= 0 && cx <= 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        } else {
                            if (cz > 0 && cx < 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        }
                        break;
                    case VG_IMKP_I:
                        if (streamRequired) {
                            if (cz <= 0 && cx >= 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        } else {
                            if (cz < 0 && cx > 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        }
                        break;
                    case VG_IMKM_I:
                        if (streamRequired) {
                            if (cz >= 0 && cx >= 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        } else {
                            if (cz > 0 && cx > 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        }
                        break;
                    case VG_JPKP_I:
                        if (streamRequired) {
                            if (cz <= 0 && cy <= 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        } else {
                            if (cz < 0 && cy < 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        }
                        break;
                    case VG_JPKM_I:
                        if (streamRequired) {
                            if (cz >= 0 && cy <= 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        } else {
                            if (cz > 0 && cy < 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        }
                        break;
                    case VG_JMKP_I:
                        if (streamRequired) {
                            if (cz <= 0 && cy >= 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        } else {
                            if (cz < 0 && cy > 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        }
                        break;
                    case VG_JMKM_I:
                        if (streamRequired) {
                            if (cz >= 0 && cy >= 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        } else {
                            if (cz > 0 && cy > 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        }
                        break;
                    // K
                    // k_out
                    case VG_IPJP_O:
                        if (streamRequired) {
                            if (cy <= 0 || cx <= 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        } else {
                            if (cy < 0 || cx < 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        }
                        break;
                    case VG_IPJM_O:
                        if (streamRequired) {
                            if (cy >= 0 || cx <= 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        } else {
                            if (cy > 0 || cx < 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        }
                        break;
                    case VG_IMJP_O:
                        if (streamRequired) {
                            if (cy <= 0 || cx >= 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        } else {
                            if (cy < 0 || cx > 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        }
                        break;
                    case VG_IMJM_O:
                        if (streamRequired) {
                            if (cy >= 0 || cx >= 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        } else {
                            if (cy > 0 || cx > 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        }
                        break;
                        // IK
                    case VG_IPKP_O:
                        if (streamRequired) {
                            if (cz <= 0 || cx <= 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        } else {
                            if (cz < 0 || cx < 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        }
                        break;
                    case VG_IPKM_O:
                        if (streamRequired) {
                            if (cz >= 0 || cx <= 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        } else {
                            if (cz > 0 || cx < 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        }
                        break;
                    case VG_IMKP_O:
                        if (streamRequired) {
                            if (cz <= 0 || cx >= 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        } else {
                            if (cz < 0 || cx > 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        }
                        break;
                    case VG_IMKM_O:
                        if (streamRequired) {
                            if (cz >= 0 || cx >= 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        } else {
                            if (cz > 0 || cx > 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        }
                        break;
                        // JK
                    case VG_JPKP_O:
                        if (streamRequired) {
                            if (cz <= 0 || cy <= 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        } else {
                            if (cz < 0 || cy < 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        }
                        break;
                    case VG_JPKM_O:
                        if (streamRequired) {
                            if (cz >= 0 || cy <= 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        } else {
                            if (cz > 0 || cy < 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        }
                        break;
                    case VG_JMKP_O:
                        if (streamRequired) {
                            if (cz <= 0 || cy >= 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        } else {
                            if (cz < 0 || cy > 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        }
                        break;
                    case VG_JMKM_O:
                        if (streamRequired) {
                            if (cz >= 0 || cy >= 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        } else {
                            if (cz > 0 || cy > 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        }
                        break;
                    // k_out end
                    // 12 edges end
                    // 8 corners
                    // inner corners
                    case VG_IPJPKP_I:
                        if (streamRequired) {
                            if (cx <= 0 && cy <= 0 && cz <= 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        } else {
                            if (cx < 0 && cy < 0 && cz < 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        }
                        break;
                    case VG_IPJPKM_I:
                        if (streamRequired) {
                            if (cx <= 0 && cy <= 0 && cz >= 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        } else {
                            if (cx < 0 && cy < 0 && cz > 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        }
                        break;
                    case VG_IPJMKP_I:
                        if (streamRequired) {
                            if (cx <= 0 && cy >= 0 && cz <= 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        } else {
                            if (cx < 0 && cy > 0 && cz < 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        }
                        break;
                    case VG_IPJMKM_I:
                        if (streamRequired) {
                            if (cx <= 0 && cy >= 0 && cz >= 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        } else {
                            if (cx < 0 && cy > 0 && cz > 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        }
                        break;
                    case VG_IMJPKP_I:
                        if (streamRequired) {
                            if (cx >= 0 && cy <= 0 && cz <= 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        } else {
                            if (cx > 0 && cy < 0 && cz < 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        }
                        break;
                    case VG_IMJPKM_I:
                        if (streamRequired) {
                            if (cx >= 0 && cy <= 0 && cz >= 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        } else {
                            if (cx > 0 && cy < 0 && cz > 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        }
                        break;
                    case VG_IMJMKP_I:
                        if (streamRequired) {
                            if (cx >= 0 && cy >= 0 && cz <= 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        } else {
                            if (cx > 0 && cy > 0 && cz < 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        }
                        break;
                    case VG_IMJMKM_I:
                        if (streamRequired) {
                            if (cx >= 0 && cy >= 0 && cz >= 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        } else {
                            if (cx > 0 && cy > 0 && cz > 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        }
                        break;
                    // out corner
                    case VG_IPJPKP_O:
                        if (streamRequired) {
                            if (cx <= 0 || cy <= 0 || cz <= 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        } else {
                            if (cx < 0 || cy < 0 || cz < 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        }
                        break;
                    case VG_IPJPKM_O:
                        if (streamRequired) {
                            if (cx <= 0 || cy <= 0 || cz >= 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        } else {
                            if (cx < 0 || cy < 0 || cz > 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        }
                        break;
                    case VG_IPJMKP_O:
                        if (streamRequired) {
                            if (cx <= 0 || cy >= 0 || cz <= 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        } else {
                            if (cx < 0 || cy > 0 || cz < 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        }
                        break;
                    case VG_IPJMKM_O:
                        if (streamRequired) {
                            if (cx <= 0 || cy >= 0 || cz >= 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        } else {
                            if (cx < 0 || cy > 0 || cz > 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        }
                        break;
                    case VG_IMJPKP_O:
                        if (streamRequired) {
                            if (cx >= 0 || cy <= 0 || cz <= 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        } else {
                            if (cx > 0 || cy < 0 || cz < 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        }
                        break;
                    case VG_IMJPKM_O:
                        if (streamRequired) {
                            if (cx >= 0 || cy <= 0 || cz >= 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        } else {
                            if (cx > 0 || cy < 0 || cz > 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        }
                        break;
                    case VG_IMJMKP_O:
                        if (streamRequired) {
                            if (cx >= 0 || cy >= 0 || cz <= 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        } else {
                            if (cx > 0 || cy > 0 || cz < 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        }
                        break;
                    case VG_IMJMKM_O:
                        if (streamRequired) {
                            if (cx >= 0 || cy >= 0 || cz >= 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        } else {
                            if (cx > 0 || cy > 0 || cz > 0) {
                                f(xiIndex, 0, 0, 0) =
                                    fStage(xiIndex, -cx, -cy, -cz);
                            }
                        }
                        break;
                    default:
                        break;
                }
            }
        }
    }
#endif  // OPS_3D
}
#endif /* APPS_LBM_DEM_PSM_KERNEL_H_ */
