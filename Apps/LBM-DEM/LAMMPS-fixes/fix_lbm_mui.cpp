/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   C. Tsigginos STFC
   Copyright 2019- STFC Daresbury Laboratory, Daresbury, UK

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_lbm_mui.h"
#include "../mui/mui.h" //to be verified
#include "atom.h"
#include "memory.h"
#include "force.h"
#include "update.h"
#include "domain.h"
#include "domain.h"
#include "error.h"
#include <string>



using namespace LAMMPS_NS;
using namespace FixConst;

FixLbmMui::FixLbmMui(LAMMPS *lmp, int narg, char **arg) :
   Fix(lmp, narg, arg) {

	int periodicFlag = 0;
	double dx;
	//create interface

	//In the future modified accordingly. Based on keys, But now
	//want to debug.

	if (domain->dimension == 2) {
		interface2d = new mui::uniface2d("mpi://DEM/ifs");

	}
	else if (domain->dimension == 3) {
		interface3d = new mui::uniface3d("mpi://DEM/ifs");
		fprintf(screen, "3D interface is created\n");
	}

	//Pass the only argument
	if (narg < 5) error->all(FLERR,"Illegal fix lbm_mui command");

	int iarg = 3;

	 while (iarg < narg) {
	    if (strcmp(arg[iarg],"dtLBM") == 0) {
	    	dtLBM = atof(arg[iarg+1]);
	    	iarg += 2;
	    } else error->all(FLERR,"Illegal fix nve/sphere command");
	  }

	 if (dtLBM <= 0.0)
	 		error->all(FLERR,"Fix lbm/mui/iter negative timestep");

	 timeInit = update->ntimestep;

	 alpha1  =  update ->dt / dtLBM;

	 //setting up nevery
	 if (alpha1 < 0.0000001)
		 alpha1 = 1;


	 if (alpha1 <= 1) {
		 nevery = static_cast<int>(1.0/alpha1);
		 update->dt = dtLBM/ static_cast<double>(nevery);
		 update->niter = 1;
	 }
	 else if (alpha1 > 1.0) {
		 nevery = 1;
		 update->niter = static_cast<int>(alpha1); //DEM subcycling scheme-OPS will handle the
	 //transformation to integer
		 update->dt = dtLBM * static_cast<double>(update->niter);
	 }
	 alpha1  =  update ->dt / dtLBM;
	 flagSetup = false;

	 Fd = NULL; //To be verified
	 Md = NULL;
	 //setting up calling back
	 grow_arrays(atom->nmax);
	 atom->add_callback(0);
	 atom->add_callback(1);

	// size_peratom_cols = 6;
	// array_flag = 1;
	// peratom_flag  = 1;
	// peratom_freq = 1;

	 //modify the DEM domain
	 if (domain->dimension==3)
		 if ((domain->periodicity[0]==1) || (domain->periodicity[1]==1)|| (domain->periodicity[2]==1)) {
			 periodicFlag = 1;
		 }
	else if (domain->dimension==2)
		if ((domain->periodicity[0]==1) || (domain->periodicity[1]==1))
			periodicFlag = 1;

	 if (periodicFlag == 1) {
		 dx = dtLBM * sqrt(3);
		 if (domain->periodicity[0]==1) {
			 domain->boxlo[0] -=0.5 * dx;
			 domain->boxhi[0] +=0.5 * dx;
		 }

		 if (domain->periodicity[1]==1) {
			 domain->boxlo[1] -=0.5 * dx;
			 domain->boxhi[1] +=0.5 * dx;
		 }

		 if (domain->dimension == 3) {
			 if (domain->periodicity[2] == 1) {
				 domain->boxlo[2] -= 0.5 * dx;
				 domain->boxhi[2] += 0.5 * dx;
			 }
		 }

		 domain->set_initial_box();
		 domain->set_global_box();
		 domain->set_local_box();
		 domain->print_box("  ");
	 }
}

FixLbmMui::~FixLbmMui() {

	atom->delete_callback(id,0);
	atom->delete_callback(id,1);

	memory->destroy(Fd);
	memory->destroy(Md);

	if (domain->dimension==2)
		delete interface2d;
	else if (domain->dimension==3)
		delete interface3d;
	//Check if I need to destroy the MUI interface
}

int FixLbmMui::setmask() {

	 int mask = 0;
	 mask |= PRE_FORCE;
	 mask |= POST_FORCE;
	 mask |= END_OF_STEP;

	 return mask;
}
void FixLbmMui::init() {

	 if (!atom->sphere_flag)
	    error->all(FLERR,"Fix lbm/mui requires atom style sphere");
}
void FixLbmMui::setup(int vflag) {



     flagSetup = false;
     if (domain->dimension==2)
    	 send_MUI_data();
     else if (domain->dimension ==3)
    	 send_MUI_data3D();



     if (domain->dimension==2)
    	 recv_MUI_data();
     else if (domain->dimension == 3) {
    	 recv_MUI_data3D();
     }

     flagSetup = true;

	 if (strstr(update->integrate_style,"verlet"))
	    post_force(vflag);


}
void FixLbmMui::pre_force(int vflag) {

	//printf("Timestep %d\n",update->ntimestep);
	if (update->ntimestep % nevery == 0) {
		if (domain->dimension == 2)
			send_MUI_data();
		else if (domain->dimension == 3)
			send_MUI_data3D();
	}


}
void FixLbmMui::post_force(int vflag) {


	//fprintf(screen,"Timestep %d\n",update->ntimestep);
	double **f = atom->f;
	double **torque = atom-> torque;
	int nlocal = atom->nlocal;
	int *mask = atom->mask;





	for (int i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit) {
			f[i][0] += Fd[i][0];
			f[i][1] += Fd[i][1];
			f[i][2] += Fd[i][2];

			torque[i][0] += Md[i][0];
			torque[i][1] += Md[i][1];
			torque[i][2] += Md[i][2];
		}

	}
}

void FixLbmMui::send_MUI_data3D() {

	double sublox, subhix, subloy, subhiy, subloz, subhiz;
	int timeCur, timeNext, data1;
	int nlocal = atom->nlocal;
	double **x = atom->x;
	double **v = atom->v;
	double **omega = atom->omega;
	double *radius = atom->radius;
	int *mask = atom->mask;

	sublox = domain->sublo[0];
	subhix = domain->subhi[0];

	subloy = domain->sublo[1];
	subhiy = domain->subhi[1];


	subloz = domain->sublo[2];
	subhiz = domain->subhi[2];

	mui::geometry::box3d send_recv_region({sublox, subloy, subloz}, {subhix, subhiy, subhiz});

	if (flagSetup == false) {
		fprintf(screen, "FlagSetup is false\n");
		if (domain->box_change == 1) {
			timeCur = 0;
			timeNext = update->nsteps + 1;
		}
		else {
			timeCur = 0;
			timeNext = update->nsteps + 10000;
		}

		interface3d->push("alpha",alpha1);
		interface3d->push("Nsteps",update->nsteps);

		interface3d->push("xPerFlag", domain->periodicity[0]);
		interface3d->push("yPerFlag", domain->periodicity[1]);
		interface3d->push("zPerFlag", domain->periodicity[2]);

		if ((domain->periodicity[0] == 1) || (domain->periodicity[1]==1) || (domain->periodicity[2]==1)) {
					interface3d->push("xper", domain->prd[0]);
					interface3d->push("yper", domain->prd[1]);
					interface3d->push("zper", domain->prd[2]);

					double cutoffx = 0.5 * comm->cutghost[0] + 0.5 * neighbor->skin;
					double cutoffy = 0.5 * comm->cutghost[1] + 0.5 * neighbor->skin;
					double cutoffz = 0.5 * comm->cutghost[2] + 0.5 * neighbor->skin;
					interface3d->push("xCutOff", cutoffx);
					interface3d->push("yCutOff", cutoffy);
					interface3d->push("zCutOff", cutoffz);

		}

		interface3d->announce_send_span(timeCur, timeNext, send_recv_region);
		interface3d->announce_recv_span(timeCur, timeNext, send_recv_region);


	}
	else {
		timeCur = update->ntimestep - timeInit;


		if (domain->box_change==1) {
			timeNext = (update-> ntimestep + 1) - timeInit;
			interface3d->announce_send_span(timeCur, timeNext, send_recv_region);
			interface3d->announce_recv_span(timeCur, timeNext, send_recv_region);
		}

	}

	//sending data to LBM solver
	//fprintf(screen, "LIGGGHTS: Number of particles of %d\n",nlocal);
	for (int i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit) {

				interface3d->push("radius",{x[i][0],x[i][1],x[i][2]},radius[i]);
				interface3d->push("u",{x[i][0],x[i][1],x[i][2]},v[i][0]);
				interface3d->push("v",{x[i][0],x[i][1],x[i][2]},v[i][1]);
				interface3d->push("w",{x[i][0],x[i][1],x[i][2]},v[i][2]);
				interface3d->push("ox",{x[i][0],x[i][1],x[i][2]},omega[i][0]);
				interface3d->push("oy",{x[i][0],x[i][1],x[i][2]},omega[i][1]);
				interface3d->push("oz",{x[i][0],x[i][1],x[i][2]},omega[i][2]);



		}
	}
	//fprintf(screen,"LIGGGHTS: Data are send at time %d\n",timeCur);
    data1 = interface3d->commit(timeCur);
   // fprintf(screen,"Data that I send: %d\n",data1);
}

void FixLbmMui::send_MUI_data() {

	double sublox, subhix, subloy, subhiy, subloz, subhiz;
	int timeCur, timeNext, data1;
	int nlocal = atom->nlocal;
	double **x = atom->x;
	double **v = atom->v;
	double **omega = atom->omega;
	double *radius = atom->radius;
	int *mask = atom->mask;

	//setting_up communication domain


	mui::geometry::box2d send_recv_region({sublox, subloy},{subhix,subhiy});

	if (flagSetup == false) {
		if (domain->box_change==1) {
			timeCur = 0;
			timeNext  = update->nsteps + 1;
		}
		else {
			timeCur = 0;
			timeNext = update->nsteps + 10000;
		}
		interface2d->push("alpha",alpha1);
		interface2d->push("Nsteps",update->nsteps); //Flag Setup-Not general

		interface2d->push("xPerFlag", domain->periodicity[0]);
		interface2d->push("yPerFlag", domain->periodicity[1]);

		if ((domain->periodicity[0] == 1) || (domain->periodicity[1]==1)) {
			interface2d->push("xper", domain->prd[0]);
			interface2d->push("yper", domain->prd[1]);

			double cutoffx = 0.5 * comm->cutghost[0] + 0.5 * neighbor->skin;
			double cutoffy = 0.5 * comm->cutghost[1] + 0.5 * neighbor->skin;
			interface2d->push("xCutOff", cutoffx);
			interface2d->push("yCutOff", cutoffy);

		}

		interface2d->announce_send_span(timeCur, timeNext, send_recv_region);
		interface2d->announce_recv_span(timeCur, timeNext, send_recv_region);
		fprintf(screen, "LIGGGHTS: I send alpha and Nsteps\n");
	}
	else {
		timeCur = update->ntimestep - timeInit;


		if (domain->box_change==1) {
			timeNext = (update-> ntimestep + 1) - timeInit;
			interface2d->announce_send_span(timeCur, timeNext, send_recv_region);
			interface2d->announce_recv_span(timeCur, timeNext, send_recv_region);
		}

	}



	//Sending data to LBM solver
	for (int i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit) {
			if (domain->dimension==2) {
				interface2d->push("radius",{x[i][0],x[i][1]},radius[i]);
				interface2d->push("u",{x[i][0],x[i][1]},v[i][0]);
				interface2d->push("v",{x[i][0],x[i][1]},v[i][1]);
				interface2d->push("oz",{x[i][0],x[i][1]},omega[i][2]);
			}
		}
	}

	data1 = interface2d->commit(timeCur); //replace with integer timestep for safety


}
void FixLbmMui::recv_MUI_data3D() {
	int timeCur;
	int nlocal = atom->nlocal;
	double **x = atom->x;
	int *mask = atom->mask;

	double dt = maxdt(update->dt, dtLBM);
	if (flagSetup == false)
		timeCur = 0;
	else
		timeCur = update->ntimestep - timeInit;

	mui::chrono_sampler_exact3d timeS;
	mui::sampler_dem3d<double> spatial;

	for (int i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit) {
			Fd[i][0] = 0.0;
			Fd[i][1] = 0.0;
			Fd[i][2] = 0.0;
			Md[i][0] = 0.0;
			Md[i][1] = 0.0;
			Md[i][2] = 0.0;

			auto FxTemp = interface3d->fetch("Fdx", {x[i][0],x[i][1],x[i][2]}, timeCur, spatial, timeS);
			auto FyTemp = interface3d->fetch("Fdy", {x[i][0],x[i][1],x[i][2]}, timeCur, spatial, timeS);
			auto FzTemp = interface3d->fetch("Fdz", {x[i][0],x[i][1],x[i][2]}, timeCur, spatial, timeS);
			auto MxTemp = interface3d->fetch("Mdx", {x[i][0],x[i][1],x[i][2]}, timeCur, spatial, timeS);
			auto MyTemp = interface3d->fetch("Mdy", {x[i][0],x[i][1],x[i][2]}, timeCur, spatial, timeS);
			auto MzTemp = interface3d->fetch("Mdz", {x[i][0],x[i][1],x[i][2]}, timeCur, spatial, timeS);

			Fd[i][0] += FxTemp;
			Fd[i][1] += FyTemp;
			Fd[i][2] += FzTemp;
			Md[i][0] += MxTemp;
			Md[i][1] += MyTemp;
			Md[i][2] += MzTemp;

		}
	}
	//printf("LIGGGHTS: I received particle data at timestep %d",timeCur);
	interface3d->forget(timeCur);
}
void FixLbmMui::recv_MUI_data() {

	int timeCur;
	//setting up the names of the receiving lists
	std::string Fx1 = "Fdx";
	std::string Fy1 = "Fdy";
	std::string Mz1 = "Mdz";

	int nlocal = atom->nlocal;
	double **x = atom->x;
	int *mask = atom->mask;

	double dt = maxdt(update->dt, dtLBM);
	//setting up time span for data Receiving
	if (flagSetup == false) {
		timeCur = 0;
	}
	else {
		timeCur = update->ntimestep - timeInit;
	}
	//Receiving timestep
	mui::chrono_sampler_exact2d   timeS;
	mui::sampler_dem2d<double> spatial;


	for (int i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit) {
			Fd[i][0] = 0.0;
			Fd[i][1] = 0.0;
			Fd[i][2] = 0.0;
			Md[i][0] = 0.0;
			Md[i][1] = 0.0;
			Md[i][2] = 0.0;


			auto FxTemp = interface2d->fetch("Fdx",{x[i][0],x[i][1]}, timeCur, spatial, timeS);
			auto FyTemp = interface2d->fetch("Fdy",{x[i][0],x[i][1]}, timeCur, spatial, timeS);
			auto MzTemp = interface2d->fetch("Mdz",{x[i][0],x[i][1]}, timeCur, spatial, timeS);

			Fd[i][0] += FxTemp;
			Fd[i][1] += FyTemp;
			Md[i][2] += MzTemp;


		}

	}

	interface2d->forget(timeCur);
	printf("LIGGGHTS: I received my drag forces for t=%d\n",timeCur);


}

void  FixLbmMui::end_of_step() {

	//send_MUI_data();
	if (domain->dimension == 2)
		recv_MUI_data();
	else if(domain->dimension == 3)
		recv_MUI_data3D();
}


double FixLbmMui::memory_usage() {

	double bytes;
	bytes = 6.0 * atom->nmax * sizeof(double);

	return bytes;

}

void FixLbmMui::grow_arrays(int nmax) {

	memory->grow(Fd, nmax, 3, "lbm_mui:Fd");
	memory->grow(Md, nmax, 3, "lbm_mui:Md");
}

int FixLbmMui::pack_exchange(int i, double *buf) {

	int m = 0;

	buf[m++] = Fd[i][0];
	buf[m++] = Fd[i][1];
	buf[m++] = Fd[i][2];
	buf[m++] = Md[i][0];
	buf[m++] = Md[i][1];
	buf[m++] = Md[i][2];

	return m;
}

int FixLbmMui::unpack_exchange(int nlocal, double *buf) {

	int m = 0;
	Fd[nlocal][0] = buf[m++];
	Fd[nlocal][1] = buf[m++];
	Fd[nlocal][2] = buf[m++];

	Md[nlocal][0] = buf[m++];
	Md[nlocal][1] = buf[m++];
	Md[nlocal][2] = buf[m++];

	return m;
}

void FixLbmMui::set_arrays(int i) {

	Fd[i][0] = 0.0;
	Fd[i][1] = 0.0;
	Fd[i][2] = 0.0;

	Md[i][0] = 0.0;
	Md[i][1] = 0.0;
	Md[i][2] = 0.0;
}

void FixLbmMui::copy_arrays(int i, int j) {

	Fd[j][0] = Fd[i][0];
	Fd[j][1] = Fd[i][1];
	Fd[j][2] = Fd[i][2];

	Md[j][0] = Md[i][0];
	Md[j][1] = Md[i][1];
	Md[j][2] = Md[i][2];

}

double FixLbmMui::maxdt(double a, double b) {

	if (a>=b)
		return a;

	return b;
}
