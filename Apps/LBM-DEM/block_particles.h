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
 * @brief   Class for handling particles at various blocks
 * @author  C. Tsigginos
 */
#ifndef BLOCK_PARTICLES_H_
#define BLOCK_PARTICLES_H_

#include <vector>
#include "particle.h"
#include "type.h"
#include "block.h"
#include "scheme.h"
#include "flowfield.h"
#include "ops_lib_core.h"
#ifdef OPS_MPI
#include "ops_mpi_core.h"
#endif
//#include "block_particles_wrapper.h"

class BlockParticles {

	public:
		int NParticles; //Local particle owned by the block
		int NPeriodic; //Periodic particle inserted in the box
		int Nmax; //Maximum allowed particle in the block
		bool owned; //Flag for defining if the given process owns the block
		std::vector<Particle> particleList; //PArticle data
		BlockParticles(int spacedim, Real dx1, Real cutoff, Block myBlock);
		~BlockParticles();
		int checkDistanceBlock();
		void initializeDragForce();
		void evaluateDragForce(Real dt);
		void updateParticle(Real* xpos, Real radius, std::vector<Real> shape);
		void insertParticle(Real* xpos, Real radius, std::vector<Real> shape,
							Real* uPart, Real* omegaT);
		void updateParticleVelocities(Real* uPart,Real* omegaT);
		void findStencil();
		void FindBoxLocalBound();
		static void  KerCarBound(const ACC<Real>& xf, Real* xb,const int* spacedim);
		void getLocalBound(Real* xBound);
		void setGlobalBound(Real* xBound);
		void getOwnership(bool flag);
	private:
		int spaceDim;
		int cutOff;
		Real* xBoundGlobal;
		Real* xBoundLocal;
		Real dx;
		int *Nf;
		Block* currentBlock;
		bool SphereParallepipedIntersection(Real* xpos,Real radius);
};



#endif
