/*
 * poros_grid.h
 *
 *  Created on: May 29, 2021
 *      Author: jpd38567
 */

#ifndef POROS_GRID_H_
#define POROS_GRID_H_

#include "particle_to_grid_base.h"


class PorosGrid : public ParticleToGridBase {
	public:
	PorosGrid(int ParticleType,int spaceDim);
	virtual int  particleShape() {return 1;}
	virtual void DefineVariables(int noElem, SizeType timestep = 0);
	virtual void ParticleProjection();
	virtual void UpdateProjection();
	virtual void MappingFunction(bool flag);
	virtual void InitializeVariables();
	virtual void ReturnParticleShape() { ops_printf("Class handles abritrary shape particles\n"); }

	protected:
	static void KerSolidFracGridSphere(ACC<int>& idp, ACC<Real>& sfp, ACC<Real>& uPar,
				ACC<Real>& xAvg,const ACC<Real>& xF, const Real* xPart, const Real* Radius,
				const int* idParticle, const Real* vParticle, const Real* omParticle,
				const Real* dx, const int* Ngrid, const int* noelem, const int* spacedim);

	static void KerSolidVelocityUpdate(ACC<Real>& vP, const ACC<int>& id,
			const ACC<Real>& xAvg, const int* idParticle,
			const Real* xPos, const Real* radPart, const Real* velP, const Real* omP,
			const Real* dx, const Real* spacedim, const Real* noelem);


	static void KerInitializePorousGrid(ACC<Real>&sfP, ACC<Real>& vP, ACC<Real>& xAvg,
			ACC<int>& id, const int* spacedim, const int* noelem);

	static inline OPS_FUN_PREFIX Real ComputeCircularSolidFractionGrid(const Real* xfl,const Real* xPart,
			const Real Rp, Real* xavg, const Real dx, const int noGrid);



	private:
		int sizeofData = 4;
		int Npoints = 9;
};





#endif /* APPS_LBM_DEM_POROS_GRID_H_ */
