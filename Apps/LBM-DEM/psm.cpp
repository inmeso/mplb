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
 * @brief   Wrap functions for the implementation of the PSM scheme
 * @author  Chrysovalantis Tsigginos
 * @details Define wrap functions for implmenting the PSM scheme
 */

#include "psm.h"

#define PI 3.1415926535897932384626433832795028841971693993751058209

Field<Real> g_Fd{"Fd"};
Field<Real> g_sfm{"sfp"};
Field<int> g_id{"id"};
Field<Real> g_vp{"vp"};

Real force[3];
int forceFlag;
int NELEM{2}; //Max. number of particles intersecting a computational cell
Real gammaX;

void SetInitialMacroVars(Real* initialMacros) {
	int size = SPACEDIM + 1;
	for (int blockIdx = 0; blockIdx < BlockNum(); blockIdx++) {
		int *iterRng = BlockIterRng(blockIdx, IterRngWhole());
		ops_par_loop(KerSetInitialMacroVars, "KerSetInitialMacroVars", g_Block[blockIdx], SPACEDIM, iterRng,
					 ops_arg_dat(g_MacroVars[blockIdx], NUMMACROVAR, LOCALSTENCIL, "double", OPS_RW),
					 ops_arg_dat(g_CoordinateXYZ[blockIdx], SPACEDIM, LOCALSTENCIL, "double", OPS_READ),
					 ops_arg_gbl(initialMacros, size, "double", OPS_READ));

	}
}

void DefinePSMVariables(SizeType timeStep, Real gamma) {

	//Define PSM model variables
	NELEM = 2;
	gammaX = gamma;

	int sizeVp = NELEM * SPACEDIM;
	int sizeFd = NELEM * SPACEDIM;

	ops_printf("MP-LBM: PSM gamma paramater is %f\n",gammaX);

	ops_decl_const("NELEM", 1, "int", &NELEM);

	//Setting up PSM variables

	g_sfm.SetDataDim(NELEM);
	if (timeStep == 0)
		g_sfm.CreateFieldFromScratch();
	else
		g_sfm.CreateFieldFromCheckPoint(timeStep);


	g_id.SetDataDim(NELEM);
	if (timeStep == 0)
		g_id.CreateFieldFromScratch();
	else
		g_id.CreateFieldFromCheckPoint(timeStep);

	g_vp.SetDataDim(sizeVp);
	g_vp.CreateFieldFromScratch();

	g_Fd.SetDataDim(sizeFd);
	if (timeStep == 0)
		g_Fd.CreateFieldFromScratch();
	else
		g_Fd.CreateFieldFromCheckPoint(timeStep);


}

void SetupBodyForces() {

	force[0] = 0.0;
	force[1] = 0.0;
#ifdef OPS_3D
	force[2] = 0.0;
#endif

	forceFlag = 0;
	for (int index = 0; index < SPACEDIM; index++)
		if (abs(force[index])> 0.000000001)
			forceFlag += 1;

	if (forceFlag > 1)
		forceFlag = 1;
}

void InitializePSMLists() {

	int size = NELEM * SPACEDIM;


	for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {
#ifdef OPS_MPI
		sub_block_list sb = OPS_sub_block_list[blockIndex];
		if (!sb->owned) continue;
#endif
		int* iterRng = BlockIterRng(blockIndex, IterRngWhole());
		ops_par_loop(KerInitialize, "KerInitialize", g_Block[blockIndex], SPACEDIM, iterRng,
					 ops_arg_dat(g_sfm[blockIndex], NELEM, LOCALSTENCIL, "double", OPS_WRITE),
					 ops_arg_dat(g_vp[blockIndex], size, LOCALSTENCIL, "double", OPS_WRITE),
					 ops_arg_dat(g_id[blockIndex], NELEM, LOCALSTENCIL, "int", OPS_WRITE),
					 ops_arg_dat(g_Fd[blockIndex], size , LOCALSTENCIL, "double", OPS_WRITE));
	}

}

void UpdateMacroVarsForceFluidSolid() {

	for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {
		int *iterRng = BlockIterRng(blockIndex, IterRngWhole());
		ops_par_loop(KerCalcMacroVarsForceFluidSolid,"KerCalcMacroVarsForceFluidSolid", g_Block[blockIndex],
					 SPACEDIM, iterRng,
					 ops_arg_dat(g_MacroVars[blockIndex], NUMMACROVAR, LOCALSTENCIL, "Real", OPS_RW),
					 ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int", OPS_READ),
					 ops_arg_dat(g_sfm[blockIndex], NELEM, LOCALSTENCIL, "double", OPS_READ),
					 ops_arg_gbl(pTimeStep(), 1, "double", OPS_READ),
					 ops_arg_gbl(force, LATTDIM, "double", OPS_READ));

	}
}


//2D functions
#ifdef OPS_2D
void EstimateSolidFractionMP2D() {
	int iminl, imaxl, jminl, jmaxl, kminl, kmaxl;
	double xPar[3];
	double velP[3];
	int intOrder  = 4;
	int idp;
	int size = NELEM * SPACEDIM;

	for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {
#ifdef OPS_MPI
		sub_block_list sb = OPS_sub_block_list[blockIndex];
		if (!sb->owned) continue;
#endif
		for (int particleIndex = 0; particleIndex < Nparticles[blockIndex]; particleIndex++) {
			 xPar[0] = xp[blockIndex][particleIndex];
			 xPar[1] = yp[blockIndex][particleIndex];
			 xPar[2] = Radius[blockIndex][particleIndex];

			 velP[0] = up[blockIndex][particleIndex];
			 velP[1] = vp[blockIndex][particleIndex];
			 velP[2] = omegaZ[blockIndex][particleIndex];

			 StenList[0] = ParticleList[blockIndex][2 * SPACEDIM * particleIndex];
			 StenList[1] = ParticleList[blockIndex][2 * SPACEDIM * particleIndex + 1];
			 StenList[2] = ParticleList[blockIndex][2 * SPACEDIM * particleIndex + 2];
			 StenList[3] = ParticleList[blockIndex][2 * SPACEDIM * particleIndex + 3];

			 idp = particleIndex;

			 ops_par_loop(KerEvaluateSolidFractionMP2D,"KerEvaluateSolidFractionMP2D", g_Block[blockIndex], SPACEDIM, StenList,
					 	  ops_arg_gbl(g_sfm[blockIndex], NELEM, LOCALSTENCIL, "int", OPS_RW),
						  ops_arg_dat(g_vp[blockIndex], size, LOCALSTENCIL, "double", OPS_RW),
						  ops_arg_dat(g_id[blockIndex], NELEM, LOCALSTENCIL, "int",OPS_RW),
						  ops_arg_dat(g_CoordinateXYZ[blockIndex], SPACEDIM, LOCALSTENCIL, "int", OPS_READ),
						  ops_arg_gbl(&idp, 1, "int", OPS_READ),
						  UpdateMacroVars3D				  ops_arg_gbl(xPar, 3, "double", OPS_READ),
						  ops_arg_gbl(velP, 3, "double", OPS_READ),
						  ops_arg_gbl(pdxLBM(), 1, "double", OPS_READ),
						  ops_arg_gbl(&intOrder,1,"int", OPS_READ));





		}
	}
}


void UpdateSolidVelocity2D() {

	int stenIndex = 2 * SPACEDIM;
	double vel_p[3], x_par[3];
	int size_v = NELEM*SPACEDIM;
	int idp;


	for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {
#ifdef OPS_MPI
		sub_block_list sb = OPS_sub_block_list[blockIndex];
		if (!sb->owned) continue;
#endif
		for (int particleIndex = 0; particleIndex < Nparticles[blockIndex]; particleIndex++) {
			vel_p[0] = up[blockIndex][particleIndex];
			vel_p[1] = vp[blockIndex][particleIndex];
			vel_p[2] = omegaZ[blockIndex][particleIndex];
				//printf("Particle velocity %12.9f,  %12.9f, %12.9f\n",vel_p[0],vel_p[1],vel_p[2]);
			x_par[0] = xp[blockIndex][particleIndex];
			x_par[1] = yp[blockIndex][particleIndex];
			x_par[2] = Radius[blockIndex][particleIndex];
			idp = particleIndex;

			StenList[0] = ParticleList[blockIndex][stenIndex*particleIndex];
			StenList[1] = ParticleList[blockIndex][stenIndex*particleIndex+1];
			StenList[2] = ParticleList[blockIndex][stenIndex*particleIndex+2];
			StenList[3] = ParticleList[blockIndex][stenIndex*particleIndex+3];

			ops_par_loop(KerSolVelocityUpdateMP2D, "KerSolVelocityUpdateMP2D",g_Block[blockIndex], SPACEDIM, StenList,
						 ops_arg_dat(g_vp[blockIndex], size_v, LOCALSTENCIL, "double", OPS_RW)),
						 ops_arg_dat(g_CoordinateXYZ[blockIndex], SPACEDIM, LOCALSTENCIL, "double", OPS_READ),
						 ops_arg_dat(g_id[blockIndex], NELEM, LOCALSTENCIL,"int", OPS_READ),
						 ops_arg_gbl(&idp, 1, "int", OPS_READ),
						 ops_arg_gbl(x_par, 3, "double", OPS_READ),
						 UpdateMacroVars3D				 ops_arg_gbl(vel_p, 3, "double", OPS_READ));





			}

		}


}


void SolidFluidCollision2D() {

	int size = NELEM * SPACEDIM;

	for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {
		int* iterRng = BlockIterRng(blockIndex, IterRngWhole());
		const Real tau{TauRef()[0]};
		ops_par_loop(KerPRE2D, "KerPRE2D",  g_Block[blockIndex], SPACEDIM, iterRng,
					ops_arg_dat(g_fStage[blockIndex], NUMXI, LOCALSTENCIL, "double", OPS_WRITE),
					ops_arg_dat(g_Fd[blockIndex], size, LOCALSTENCIL, "double", OPS_WRITE),
					ops_arg_dat(g_f[blockIndex], NUMXI, LOCALSTENCIL, "double", OPS_READ),
					ops_arg_dat(g_MacroVars[blockIndex], NUMMACROVAR, LOCALSTENCIL, "double", OPS_READ),
					ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int", OPS_READ),
					ops_arg_dat(g_id[blockIndex], NELEM, LOCALSTENCIL, "int", OPS_READ),
					ops_arg_dat(g_sfm[blockIndex], NELEM, LOCALSTENCIL, "double", OPS_READ),
					ops_arg_dat(g_vp[blockIndex], size, LOCALSTENCIL, "double", OPS_READ),
					ops_arg_gbl(&tau, 1, "double", OPS_READ),
					ops_arg_gbl(pTimeStep(), 1, "double", OPS_READ),
					ops_arg_gbl(&forceFlag, 1, "int", OPS_READ),
					ops_arg_gbl(force, SPACEDIM, "double", OPS_READ));

	}
}

void EstimateDragForce2ndOrder2D() {

	double xpos[SPACEDIM];
	double Fdp[3];
	int size = NELEM * SPACEDIM;
	int idp, nlocal;
	Real dx = dxLBM();
	const Real tau{TauRef()[0]};

	for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {
#ifdef OPS_MPI
		sub_block_list sb = OPS_sub_block_list[blockIndex];
		if (!sb->owned) continue;
#endif

		nlocal = Nparticles[blockIndex];
		for (int particleIndex = 0; particleIndex < nlocal; particleIndex++) {


			xpos[0] = xp[blockIndex][particleIndex];
			xpos[1] = yp[blockIndex][particleIndex];
			Fdp[0] = 0.0;
			Fdp[1] = 0.0;
			Fdp[2] = 0.0;
			idp = particleIndex;
			StenList[0] = ParticleList[blockIndex][4*particleIndex];
			StenList[1] = ParticleList[blockIndex][4*particleIndex+1];
			StenList[2] = ParticleList[blockIndex][4*particleIndex+2];
			StenList[3] = ParticleList[blockIndex][4*particleIndex+3];

			ops_par_loop(KerDragForcePRE2D, "KerDragForcePRE2D", g_Block[blockIndex], SPACEDIM, StenList,
					     ops_arg_dat(g_id[blockIndex], NELEM, LOCALSTENCIL, "int", OPS_READ),
						 ops_arg_dat(g_sfm[blockIndex],NELEM, LOCALSTENCIL, "double", OPS_READ),
						 ops_arg_dat(g_CoordinateXYZ[blockIndex], LATTDIM, LOCALSTENCIL, "double", OPS_READ),
						 ops_arg_dat(g_Fd[blockIndex],size,LOCALSTENCIL,"double", OPS_READ),
					     ops_arg_gbl(&idp, 1, "int", OPS_READ),
						 ops_arg_gbl(xpos, SPACEDIM, "double", OPS_READ),
						 ops_arg_gbl(Fdp, 3, "double", OPS_READ),
						 ops_arg_gbl(pTimeStep(),1,"double", OPS_READ),
						 ops_arg_gbl(&tau, 1, "double", OPS_READ));

			FDrag[blockIndex][3*particleIndex] +=  Fdp[0] * dx * dx;
			FDrag[blockIndex][3*particleIndex+1] += Fdp[1] * dx * dx;
			FDrag[blockIndex][3*particleIndex+2] += Fdp[2] * dx * dx;

		}

	}
}

//Auxiliary sf calculation functions

void GauUpdateMacroVars3DssIntegrationTable(double *t, double* W, const int n_int) {
	if (n_int==0) {

		t[0] = -1.;
		t[1] =  1.;
		W[1] = 1./2.;
		W[0] = 1./2.;
	}
	else {

		if (n_int==1) {
			t[0] = 0.0;
			W[0] = 2.0;
		}
		else if (n_int == 2) {
			t[0] = -0.577350269189626;
			t[1] = -0.577350269189626;

			W[1] = 1.0;
			W[0] = 1.0;
		}
		else if (n_int == 3) {
			t[0] = -0.774596669241483;
			t[1] = 0.0;
			t[2] =  0.774596669241483;

			W[0] = 0.55555555555556;
			W[1] = 0.88888888888889;
			W[2] = 0.55555555555556;
		}
		else if (n_int == 4) {
			t[0] = -0.861136311594053;
			t[1] = -0.339981043584856;
			t[2] =  0.339981043584856;
			t[3] =  0.861136311594053;

			W[0] =  0.347854845137454;
			W[1] =  0.652145154862546;
			W[2] =  0.652145154862546;
			W[3] =  0.347854845137454;
 		}
		else if (n_int == 5) {
			t[0] = -0.906179845938664;
			t[1] = -0.538469310105683;
			t[2] =  0.0;
			t[3] =  0.538469310105683;
			t[4] =  0.906179845938664;

			UpdateMacroVars3D	W[0] =  0.236926885056189;
			W[1] =  0.478628670499366;
			W[2] =  0.568888888888889;
			W[3] =  0.478628670499366;
			W[4] =  0.236926885056189;
		}
		else if (n_int == 6) { //Copy paste from table
			t[0] = -0.932469514203152;
			t[1] = -0.661209386466265;
			t[2] = -0.238619186083197;
			t[3] =  0.238619186083197;
			t[4] =  0.661209386466265;
			t[5] =  0.932469514203152;

			W[0] =  0.171324492379170;
			W[1] =  0.360761573048139;
			W[2] =  0.467913934572691;
			W[3] =  0.467913934572691;
			W[4] =  0.360761573048139;
			W[5] =  0.171324492379170;
		}
		else if (n_int == 7) { //copy paste from table
			t[0] = -0.949107912342759;
			t[1] = -0.741531185599394;
			t[2] = -0.405845151377397;
			t[3] =  0.0;
			t[4] =  0.405845151377397;
			t[5] =  0.741531185599394;
			t[6] =  0.949107912342759;

			W[0] = 0.129484966168870;
			W[1] = 0.279705391489277;
			W[2] = 0.381830050505119;
			W[3] = 0.417959183673469;
			W[4] = 0.381830050505119;
			W[5] = 0.279705391489277;
			W[6] = 0.129484966168870;
		}
		else if (n_int == 8) { //copy-paste from table
			t[0] = -0.960289856497536;
			t[1] = -0.796666477413627;
			t[2] = -0.525532409916329;
			t[3] = -0.183434642495650;
			t[4] =  0.183434642495650;
			t[5] =  0.525532409916329;
			t[6] =  0.796666477413627;
			t[7] =  0.960289856497536;

			W[0] = 0.101228536290376;
			W[1] = 0.222381034453374;
			W[2] = 0.313706645877887;
			W[3] = 0.362683783378362;
			W[4] = 0.362683783378362;
			W[5] = 0.313706645877887;
			W[6] = 0.222381034453374;
			W[7] = 0.101228536290376;
		}
	}

}

double sign(double x) {
	if (x==0)
		return 1;

	return x/fabs(x);
}

double ReturnRoot(double xmin,double xmax, int target, double prod, double Rp, double ym) {

	double xm;
	double am;
	//printf("I entered in return root with prod %f\n",prod);
	//printf("xmin = %f, xmax = %f, target = %d, prod = %f, Rp = %f, ym = %f\n", xmin, xmax, target, prod, Rp, ym);
	//Two cases for fine discretization
	if (prod > 0.0)  {// xmin * xmax > 0 : 1 root {
		am = Rp * Rp - ym * ym;

		//printf("am = %12.10e\n",am);
		if (am < 0.0)
			am = -am;
		return xm = sign(xmin) * (sqrt(am));
	}
	else if (prod == 0.0) {
		//printf("I have a prod equal to zero\n");
		am = Rp * Rp - ym * ym;
		//printf("am = %12.9e\n",am);
		if (am < 0.0)
			am = fabs(am);
		return xm = sign(xmin) * sqrt(am);
	}

	else { //Different sign
		int a = 0;
		int b = 0;

		xm = sqrt(Rp * Rp - ym * ym);
		//printf("Rp = %f, ym =%f, xm=%f\n",Rp,ym,xm);
		if ((xm >= xmin) && (xm <= xmax))
			a = 1;

		if ((-xm >= xmin) && (-xm <= xmax))
			b = 1;



		if ((a+b)==1) {
			if (a==1)
				return xm;
			if (b==1)
				return -xm;
		}
		//else if (a+b==2) {

		//}

	}
	return 0;
}

Real EvaluateIntegral(Real x, Real ym, Real Rp, Real dx,int n_int) {

    Real *t;
    Real *W;
    if ((n_int < 0) && (n_int > 8)) {
    	printf("Integration order not consistent!!\n Use value between 0 and 8\nProgram terminating\n");
        exit (EXIT_FAILURE);
    }
    if (n_int==0) {
    	t = new Real[2];
    	W = new Real[2];
    }
    else {
    	t = new Real[n_int];
    	W = new Real[n_int];
    }
    GaussIntegrationTable(t,W,n_int);
    if (n_int==0)
    	n_int=2;

    Real integral = 0.0;
    Real ax, bx, ft, xt;


    Real xmin, xmax;
    Real ymin, ymax;
    Real x1 = 0.0;
    Real x2 = 0.0;
    Real prod;
    //find corners f(x)
    xmin = x - 0.5 * dx;
    xmax = x + 0.5 * dx;

    //NOTE: The if structure is used to define the variables in cases of non-fuinction definition
    if (Rp*Rp - xmin * xmin > 0)
       	ymin = sqrt(Rp * Rp - xmin * xmin) - ym;
       else
       	ymin = -100.;
       if (Rp*Rp - xmax * xmax > 0 )
       	ymax = sqrt(Rp * Rp - xmax * xmax) - ym;
       else
       	ymax = -100.;
       prod = xmin * xmax;

    //Assume fine discretization

    //Update xmin and xmax
       if (ymin <= 0) //These solutions are for circle
          	xmin = ReturnRoot(xmin,xmax, -1, prod, Rp, ym);
          if (ymax <= 0)
              xmax = ReturnRoot(xmin,xmax, 1, prod,Rp,ym);

              //xmin = -sqrt(Rp * Rp - ym * ym); //Correct it or not

          if (ymin >= dx) { //modify xmin
          	x2 = -xmin;
          	xmin = ReturnRoot(xmin,xmax, -1, prod, Rp, ym+dx);
          	x2 = x2 + xmin;
          	x2 = dx * x2;
          }



          	//xmax = sqrt(Rp * Rp - ym * ym);

          if (ymax >= dx) {//modify xmax
          	x1 = xmax;
          	xmax = ReturnRoot(xmin, xmax, 1, prod, Rp, ym+dx);
          	x1 = x1 - xmax;
          	x1 = dx * x1;
          }

    //Transform the integral into \int_{-1}^{1}
    ax = (xmax - xmin) * 0.5;
    bx = (xmax + xmin) * 0.5;
    for (int iIndex = 0; iIndex < n_int; iIndex++) {
        xt = bx + ax * t[iIndex];
        ft = sqrt(Rp * Rp - xt * xt) - ym;
        integral += ft * W[iIndex];
    }
   // printf("I evaluate the integral\n");
    integral=integral*(xmax-xmin)*0.5;
    delete[] t;
    delete[] W;
    return integral + x1 + x2;
}



Real  ComputeCellSF(Real* xf,Real* xp,const Real dx,const int intOrder) {

	Real vol = dx * dx;
	Real xl, yl, ryx;
	xl = xf[0] - xp[0];
	yl = xf[1] - xp[1];
	Real Rp = xp[2];

	if (xl==0)
		ryx = 1000000;
	else
		ryx = fabs(yl/xl);

	if ((ryx <= 1.0) && (xl < 0))
		return EvaluateIntegral(yl ,xl - 0.5 * dx, Rp, dx, intOrder) / vol;
	else if ((ryx <= 1.0 ) && (xl < 0))
		return EvaluateIntegral(yl , -(xl + 0.5 * dx), Rp, dx, intOrdUpdateMacroVars3Der) / vol;
	else if ((ryx > 1.0) && (yl >=0)) { //y>0 y/x>1
		return EvaluateIntegral(xl ,(yl - 0.5 *dx), Rp, dx, intOrder) / vol;
	}
	else if ((ryx > 1.0) && (yl < 0)) {
		return EvaluateIntegral(xl , -(yl + 0.5 *dx), Rp, dx, intOrder) / vol;
	}

	return 0;

}

void EstimateVelocity2D(Real* vf,const Real* vel_p,const Real* xp,Real *xfl) {

	double dx = xfl[0] - xp[0];
	double dy = xfl[1] - xp[1];
	double nx, ny, dsqrt;

	double dd = dx * dx + dy * dy;

	if (dd > xp[2] * xp[2]) { //Projection in particle surface
		dsqrt = sqrt(dd);

		nx = dx / dsqrt;
		ny = dy / dsqrt;
		dx = xp[2] * nx;
		dy = xp[2] * ny;
	}

	vf[0] = vel_p[0] - vel_p[2] * dy;
	vf[1] = vel_p[1] + vel_p[2] * dx
}#include "evolution3d.h"




void StreamCollisionPSM2D() {

	int checkDistance;

	//Mapping particles into fluid domain
	checkDistance = CheckDistance();
	if (checkDistance == 1)  {
		UpdateOldParticleLocation();
		ParticleMapping();

		InitializePSMLists();
		EstimateSolidFractionMP2D();
	}
	else
		UpdateSolidVelocity2D();


	//Update fluid velocities
	UpdateMacroVars();

	if (forceFlag == 1)
		UpdateMacroVarsForceFluidSolid();

	CorrectFluidVelocities();

	CopyBlockEnvelopDistribution2D(g_fStage, g_f);

	SolidFluidCollision2D();

	Stream();

	//TransferHalos(HaloGroups()); We do not identified a function

	ImplementBoundary();

	EstimateDragForce2ndOrder();



}
#endif

void EstimateSolidFractionMP3D() {

	int iminl, imaxl, jminl, jmaxl, kminl, kmaxl;
	double xPar[4];
	double velP[6];
	int idp, nlocal;
	int size = NELEM * SPACEDIM;
	int sPar = SPACEDIM + 1;
	for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {
#ifdef OPS_MPI
		sub_block_list sb = OPS_sub_block_list[blockIndex];
		if (!sb->owned) continue;
#endif
		nlocal = Nparticles[blockIndex] + Nperiodic[blockIndex];
		for (int particleIndex = 0; particleIndex < nlocal; particleIndex++) {
			 xPar[0] = xp[blockIndex][particleIndex];
			 xPar[1] = yp[blockIndex][particleIndex];
			 xPar[2] = zp[blockIndex][particleIndex];
			 xPar[3] = Radius[blockIndex][particleIndex];

			 velP[0] = up[blockIndex][particleIndex];
			 velP[1] = vp[blockIndex][particleIndex];
			 velP[2] = wp[blockIndex][particleIndex];
			 velP[3] = omegaX[blockIndex][particleIndex];
			 velP[4] = omegaY[blockIndex][particleIndex];
			 velP[5] = omegaZ[blockIndex][particleIndex];

			 StenList[0] = ParticleList[blockIndex][2 * SPACEDIM * particleIndex];
			 StenList[1] = ParticleList[blockIndex][2 * SPACEDIM * particleIndex + 1];
			 StenList[2] = ParticleList[blockIndex][2 * SPACEDIM * particleIndex + 2];
			 StenList[3] = ParticleList[blockIndex][2 * SPACEDIM * particleIndex + 3];
			 StenList[4] = ParticleList[blockIndex][2 * SPACEDIM * particleIndex + 4];
			 StenList[5] = ParticleList[blockIndex][2 * SPACEDIM * particleIndex + 5];

			 idp = particleIndex;
			// ops_printf("Particle %d, Stencil: [%d %d %d %d %d %d]\n",idp, StenList[0], StenList[1], StenList[2],
			//		 	 StenList[3], StenList[4], StenList[5]);
			 ops_par_loop(KerEvaluateSolidFractionMP3D,"KerEvaluateSolidFractionMP3D", g_Block[blockIndex], SPACEDIM, StenList,
					 	  ops_arg_dat(g_id[blockIndex], NELEM, LOCALSTENCIL, "int",OPS_RW),
					 	  ops_arg_dat(g_sfm[blockIndex], NELEM, LOCALSTENCIL, "double", OPS_RW),
					 	  ops_arg_dat(g_vp[blockIndex], SPACEDIM * NELEM, LOCALSTENCIL, "double", OPS_RW),
						  ops_arg_dat(g_CoordinateXYZ[blockIndex], SPACEDIM, LOCALSTENCIL, "int", OPS_READ),
						  ops_arg_gbl(&idp, 1, "int", OPS_READ),
						  ops_arg_gbl(xPar, sPar, "double", OPS_READ),
						  ops_arg_gbl(velP, 2 * SPACEDIM, "double", OPS_READ),
						  ops_arg_gbl(pdxLBM(), 1, "double", OPS_READ));
		}
	}
}

void UpdateSolidVelocity3D() {
	int ix, iy;
	double vel_p[6], x_par[4];
	int size = NELEM * SPACEDIM;
	int idp;
	int stenIndex = 2 * SPACEDIM;
	int nlocal;
	for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {
#ifdef OPS_MPI
		sub_block_list sb = OPS_sub_block_list[blockIndex];
		if (!sb->owned) continue;
#endif
			nlocal = Nparticles[blockIndex] + Nperiodic[blockIndex];
		for (int particleIndex = 0; particleIndex < nlocal; particleIndex++) {
			vel_p[0] = up[blockIndex][particleIndex];
			vel_p[1] = vp[blockIndex][particleIndex];
			vel_p[2] = wp[blockIndex][particleIndex];
			vel_p[3] = omegaX[blockIndex][particleIndex];
			vel_p[4] = omegaY[blockIndex][particleIndex];
			vel_p[5] = omegaZ[blockIndex][particleIndex];

			x_par[0] = xp[blockIndex][particleIndex];
			x_par[1] = yp[blockIndex][particleIndex];
			x_par[2] = zp[blockIndex][particleIndex];
			x_par[3] = Radius[blockIndex][particleIndex];
			idp = particleIndex;

			StenList[0] = ParticleList[blockIndex][stenIndex * particleIndex];
			StenList[1] = ParticleList[blockIndex][stenIndex * particleIndex+1];
			StenList[2] = ParticleList[blockIndex][stenIndex * particleIndex+2];
			StenList[3] = ParticleList[blockIndex][stenIndex * particleIndex+3];
			StenList[4] = ParticleList[blockIndex][stenIndex * particleIndex+4];
			StenList[5] = ParticleList[blockIndex][stenIndex * particleIndex+5];

			ops_par_loop(KerSolVelocityUpdateMP3D,"KerSolVelocityUpdateMP3D",g_Block[blockIndex], SPACEDIM, StenList,
					     ops_arg_dat(g_vp[blockIndex], size, LOCALSTENCIL, "double", OPS_RW),
						 ops_arg_dat(g_id[blockIndex], NELEM, LOCALSTENCIL,"int", OPS_READ),
						 ops_arg_dat(g_CoordinateXYZ[blockIndex], SPACEDIM, LOCALSTENCIL, "double", OPS_READ),
						 ops_arg_gbl(&idp, 1, "int", OPS_READ),
						 ops_arg_gbl(x_par, 4, "double", OPS_READ),
						 ops_arg_gbl(vel_p, 6, "double", OPS_READ));
		}
	}

}

void SolidFluidCollision3D() {

	int size = NELEM * SPACEDIM;
	const Real tau{TauRef()[0]};
	for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {
		int* iterRng = BlockIterRng(blockIndex, IterRngWhole());

		ops_par_loop(KerPRE3D, "KerPRE3D",  g_Block[blockIndex], SPACEDIM, iterRng,
					ops_arg_dat(g_fStage[blockIndex], NUMXI, LOCALSTENCIL, "double", OPS_WRITE),
					ops_arg_dat(g_Fd[blockIndex], size, LOCALSTENCIL, "double", OPS_WRITE),
					ops_arg_dat(g_f[blockIndex], NUMXI, LOCALSTENCIL, "double", OPS_READ),
					ops_arg_dat(g_MacroVars[blockIndex], NUMMACROVAR, LOCALSTENCIL, "double", OPS_READ),
					ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int", OPS_READ),
					ops_arg_dat(g_id[blockIndex], NELEM, LOCALSTENCIL, "int", OPS_READ),
					ops_arg_dat(g_sfm[blockIndex],NELEM, LOCALSTENCIL, "double", OPS_READ),
					ops_arg_dat(g_vp[blockIndex], size, LOCALSTENCIL, "double", OPS_READ),
					ops_arg_gbl(&tau, 1, "double", OPS_READ),
					ops_arg_gbl(pTimeStep(), 1, "double", OPS_READ),
					ops_arg_gbl(&forceFlag, 1, "int", OPS_READ),
					ops_arg_gbl(force, SPACEDIM, "double", OPS_READ));

	}


}


void Stream3DPeriodic() {
    for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {
        int* iterRng = BlockIterRng(blockIndex, IterRngWhole());
        ops_par_loop(
        		KerStreamPeriodic3D, "KerStreamPeriodic3D", g_Block[blockIndex], SPACEDIM, iterRng,
            ops_arg_dat(g_f[blockIndex], NUMXI, LOCALSTENCIL, "double", OPS_RW),
            ops_arg_dat(g_fStage[blockIndex], NUMXI, ONEPTLATTICESTENCIL,
                        "double", OPS_READ),
            ops_arg_dat(g_NodeType[blockIndex], NUMCOMPONENTS, LOCALSTENCIL,
                        "int", OPS_READ),
            ops_arg_dat(g_GeometryProperty[blockIndex], 1, LOCALSTENCIL, "int",
                        OPS_READ));
    }
}

void EstimateDragForce2ndOrder3D() {

	double xpos[3];
	int sizeFdp = 2 * SPACEDIM;
	double Fdp[sizeFdp];
	const Real tau{TauRef()[0]};

	int sizeFd = SPACEDIM * NELEM; //Different for 2D and 3D
	int idp;
	double dx =dxLBM();

	for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {
#ifdef OPS_MPI
		sub_block_list sb = OPS_sub_block_list[blockIndex];
		if (!sb->owned) continue;
#endif
		int nlocal = Nparticles[blockIndex] + Nperiodic[blockIndex];
		for (int particleIndex = 0; particleIndex < nlocal; particleIndex++) {
			xpos[0] = xp[blockIndex][particleIndex];
			xpos[1] = yp[blockIndex][particleIndex];
			xpos[2] = zp[blockIndex][particleIndex];
			for (int idir = 0; idir < sizeFdp; idir++)
				Fdp[idir] = 0.0;

			idp = particleIndex;
			StenList[0] = ParticleList[blockIndex][2 * SPACEDIM * particleIndex];
			StenList[1] = ParticleList[blockIndex][2 * SPACEDIM * particleIndex + 1];
			StenList[2] = ParticleList[blockIndex][2 * SPACEDIM * particleIndex + 2];
			StenList[3] = ParticleList[blockIndex][2 * SPACEDIM * particleIndex + 3];
			StenList[4] = ParticleList[blockIndex][2 * SPACEDIM * particleIndex + 4];
			StenList[5] = ParticleList[blockIndex][2 * SPACEDIM * particleIndex + 5];

			ops_par_loop(KerDragPRE3D, "KerDragPRE3D", g_Block[blockIndex], SPACEDIM, StenList,
					 	 ops_arg_dat(g_id[blockIndex], NELEM, LOCALSTENCIL, "int", OPS_READ),
						 ops_arg_dat(g_sfm[blockIndex], NELEM, LOCALSTENCIL, "double", OPS_READ),
						 ops_arg_dat(g_CoordinateXYZ[blockIndex], SPACEDIM, LOCALSTENCIL, "double", OPS_READ),
						 ops_arg_dat(g_Fd[blockIndex], sizeFd, LOCALSTENCIL, "double", OPS_READ),
						 ops_arg_gbl(Fdp, sizeFdp, "double", OPS_READ),
						 ops_arg_gbl(&idp, 1, "int", OPS_READ),
						 ops_arg_gbl(xpos, SPACEDIM, "double", OPS_READ),
						 ops_arg_gbl(pTimeStep(), 1, "double", OPS_READ),
						 ops_arg_gbl(&tau, 1, "double", OPS_READ));

			for (int idir = 0; idir < sizeFdp; idir++)
				FDrag[blockIndex][sizeFd * particleIndex + idir] += Fdp[idir] * dx * dx * dx;
		}
	}


}

//Auxiliary 3D functions

Real CalculateSolidFractionSpheres(Real* xf,Real Rc,Real* xp,Real Rp) {

	Real d, x, sf, h1, h2, V1, V2, Vol;

	d = (xp[0] - xf[0]) * (xp[0] - xf[0]) + (xp[1] - xf[1]) * (xp[1] - xf[1]) + (xp[2] - xf[2]) * (xp[2] - xf[2]);

	d = sqrt(d);

	if (d < Rc + Rp) {

		if (d <=fabs(Rp-Rc))
			return 1.0;

		x = (d * d - Rp * Rp + Rc * Rc)/(2.0 * d);

		h1 = Rc - x;
		h2 = Rp - d + x;
		V1 = (1.0/3.0) * PI * h1 * h1 * (3 * Rc - h1);
		V2 = (1.0/3.0) * PI * h2 * h2 * (3 * Rp - h2);
		Vol = (4.0/3.0) * PI * Rc * Rc * Rc;
		return (V1 + V2) / Vol;

	}

	return 0.0;
}

Real  ComputeSolidFraction3D(const Real* xfl,const Real *xp, Real dx,int nmc) {
	int Nx = nmc + 2;
	int Ny = nmc + 2;
	int Nz = nmc + 2;
	int inter = 0;
	Real dxp, dyp, dzp, dd;
	Real xmin, ymin, zmin;
	Real xtrial, ytrial, ztrial;
	Real dxc = dx/static_cast<Real>(Nx-1);

	xmin = xfl[0] - 0.5 * dx;
	ymin = xfl[1] - 0.5 * dx;
	zmin = xfl[2] - 0.5 * dx;

	for (int ip = 0; ip < Nx; ip ++) {
		xtrial = xmin + static_cast<Real>(ip) * dxc;
		dxp = (xtrial - xp[0]) * (xtrial - xp[0]);
		for (int jp = 0; jp < Ny; jp++) {
			ytrial = ymin + static_cast<Real>(jp) * dxc;
			dyp = (ytrial - xp[1]) * (ytrial - xp[1]);

			for (int kp = 0; kp < Nz; kp++) {
				ztrial = zmin + static_cast<Real>(kp) * dxc;
				dzp = (ztrial - xp[2]) * (ztrial - xp[2]);
				dd = dxp + dyp + dzp;
				if (dd <= xp[3] * xp[3])
					inter += 1;
			}
		}
	}

	return  static_cast<Real>(inter)/static_cast<Real>(Nx * Ny * Nz);
}

void EstimateVelocity3D(Real* vf,const Real* vel_p,const Real* xp, Real *xfl) {

	double dx = xfl[0] - xp[0];
	double dy = xfl[1] - xp[1];
	double dz = xfl[2] - xp[2];
	double nx, ny, nz, dsqrt;
	double dd = dx * dx + dy * dy + dz * dz;
	if (dd > xp[3] * xp[3]) {//Projection in the surface
		dsqrt = sqrt(dd);
		nx = dx / dsqrt;
		ny = dy / dsqrt;
		nz = dz / dsqrt;
		dx = xp[3] * nx;
		dy = xp[3] * ny;
		dz = xp[3] * nz;
	}

	//finding velocities
	vf[0] = vel_p[0] + vel_p[4] * dz - vel_p[5] * dy;
	vf[1] = vel_p[1] + vel_p[5] * dx - vel_p[3] * dz;
	vf[2] = vel_p[2] + vel_p[3] * dy - vel_p[4] * dx;
}

Real CalcBodyForce2ndOrder(const int xiIndex,const Real rho,const Real* u, const Real*  accel) {

	Real sum, sum1;

	sum = 0.0;

	for (int iIndex = 0; iIndex < SPACEDIM; iIndex++) {
		sum1 = XI[xiIndex * SPACEDIM + iIndex] * CS;
		for (int jIndex = 0; jIndex < SPACEDIM; jIndex++) {
			sum1 += XI[SPACEDIM * xiIndex + iIndex] * CS * XI[SPACEDIM * xiIndex + jIndex] * CS * u[jIndex];
			if (iIndex == jIndex)
				sum1 += -1.0 * u[jIndex];

		}
		sum += sum1 * accel[iIndex];
	}

	return sum * WEIGHTS[xiIndex] * rho;
}


void StreamCollisionPSM3D(int flag) {

	int checkDistance;

	checkDistance = CheckDistance();
	if (flag == 1)
		checkDistance = 1;

	if (checkDistance == 1)  {

		if (periodicFlag == 1) {

			PeriodicExchange();

			InitializePeriodicDragForce();


		}

		UpdateOldParticleLocation();

		ParticleMapping();
		InitializePSMLists();

		EstimateSolidFractionMP3D();

	}
	else {
		if (periodicFlag == 1)
			ForwardComm();


		UpdateSolidVelocity3D();
	}
	//ops_printf("Pass solid mapping-Ready to update Fluid Vels\n");


	double time = 0.0;

	UpdateMacroVars3D();

	//ops_printf("Pass macro Calc\n");

	if (forceFlag == 1)
		UpdateMacroVarsForceFluidSolid();


	//ops_printf("Pass Fluid velocities\n");


	CopyBlockEnvelopDistribution3D(g_fStage, g_f);

	//ops_printf("Pass Copy Distribution\n");


	SolidFluidCollision3D();

//ops_printf("Pass Stream Collision 3D\n");

	if (periodicFlag == 1)
		PeriodicHaloTransfer();


	if (periodicFlag == 1)
		Stream3DPeriodic();
	else
		Stream3D();



	//ops_printf("Pass streaming\n");
	if (periodicFlag == 0)
		TransferHalos(HaloGroups());


	//ops_printf("Pass Halo Transfer\n");
	if (periodicFlag == 1)
		ImplementBoundary3DPeriodic();
	else
		ImplementBoundary3D();



  //  ops_printf("Pass boundary condition\n");



    EstimateDragForce2ndOrder3D();
    if (periodicFlag == 1)
    	ReverseComm();

   // ops_printf("Pass Drag Force\n");


}

void WritePSMVariablesToHdf5(const long TimeStep) {

	std::string caseName =  CaseName();

	for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {
		 std::string blockName("Block_");
		 std::string label(std::to_string(blockIndex));
		 std::string time(std::to_string(TimeStep));

		 blockName += (label + "_" + time);
		 std::string fileName = caseName + "_" + blockName + ".h5";

		 ops_fetch_block_hdf5_file(g_Block[blockIndex], fileName.c_str());

		 ops_fetch_dat_hdf5_file(g_Fd[blockIndex], fileName.c_str());
		 ops_fetch_dat_hdf5_file(g_vp[blockIndex], fileName.c_str());
		 ops_fetch_dat_hdf5_file(g_sfm[blockIndex], fileName.c_str());
		 ops_fetch_dat_hdf5_file(g_id[blockIndex], fileName.c_str());
	}



}

#include "psm_kernel.h"



