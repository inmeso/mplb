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
#ifdef FIX_CLASS

FixStyle(lbm3d/mui,FixLbm3DMui);

#else

#ifndef FIX_LBM3D_MUI_H_
#define FIX_LBM3D_MUI_H_

#include "fix.h"
#include "../mui/mui.h" //Define location

namespace LAMMPS_NS {

class FixLbm3DMui : public Fix {
public:

	FixLbm3DMui(class LAMMPS*, int, char **);
	~FixLbm3DMui();

	virtual int setmask();
	virtual double memory_usage();
	virtual void grow_arrays(int );
	virtual void copy_arrays(int, int);
	virtual void set_arrays(int);
	virtual int pack_exchange(int, double *);
	virtual int unpack_exchange(int, double *);

	virtual void end_of_step();
	virtual void post_force(int);
	virtual void pre_force(int);

	virtual void init();
	virtual void setup(int);

private:
	void send_MUI_data();
	void recv_MUI_data();




	double maxdt(double, double);
	double dtLBM; //Timestep of LBM, actually it is fixed.
	double **Fd; //Drag force
	double **Md; //Drag moment

	mui::uniface3d *interface;
	int timeInit; //Initial time step that the fix is constructed.
	bool flagSetup;
	double alpha1; //Timestep ratio.
	int boxExpand;
};
}


#endif
#endif /* FIX_LBM3D_MUI_H_ */
