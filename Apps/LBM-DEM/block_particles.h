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
#include "dem_data.h"
#include <string>
#include <new>
#ifdef OPS_MPI
#include "ops_mpi_core.h"
#endif
//#include "block_particles_wrapper.h"

class BlockParticles {

	public:
		int NParticles; 		//number of particles owned by the block
		int NPeriodic; 			//Periodic particle inserted in the block
		int Nmax;			   //Maximum particles inserted in the given block
		bool owned; 		   //Block is owned or not
		std::vector<Particle> particleList;  	//Particle data
		BlockParticles(int spacedim, Real dx1, Real cutoff,const Block& myBlock,
				ParticleShapeDiscriptor particleType,
				int nInputExtraVariables = 0, int nOutputExtraVariables = 0);
		~BlockParticles();
		int CheckDistanceBlock();
		void InitializeDragForce();
		void EvaluateDragForce(Real dt);
		void ExtractDragForce(Real* Fd, Real* Td,int iParticle);
		void ExtractPositions(Real* xPos,int iParticle);
		int UpdateParticle(Real* xpos, Real radius, std::vector<Real> shape);
		int InsertParticle(Real* xpos, Real radius, std::vector<Real> shape,
							Real* uPart, Real* omegaT, std::vector<Real>& inputData);
		void UpdateParticleVelocities(Real* uPart,Real* omegaT);
		void FindStencil();
		static void  KerCarBound(const ACC<Real>& xf, Real* xb,const int* spacedim);
		void GetLocalBound(Real* xBound);
		void SetGlobalBound(Real* xBound);
		void GetOwnership();
		ParticleShapeDiscriptor GetParticleShape() { return particleShape;};
		void UpdateOldParticlePosition();
		void ExtractBound(Real* xMin, Real* xMax);
		void GetAdditionalOutputVariables(std::vector<Real>& output,int iParticle);
		void GetAdditionalInputVariables(std::vector<Real>& inputVariables,int idParticle);
		bool OwnedStatus() { return owned;};
		bool HasExtraInput() { return hasExtraInputVariables;};
		bool HasExtraOutput() {return hasExtraOutputVariables;};
		void SetNfLocal(int* Ndata);
		void SetLocalBound(Real* xBound);


		inline void ClearParticles() {NParticles = 0;
					NPeriodic = 0;};
		const Block& GetBlock() { return currentBlock; };
	private:
		int spaceDim;			//Spatial size
		Real cutOff;			//Distance for not updating the particle mapping
		Real xBoundGlobal[6];   //Region of block
		Real xBoundLocal[6];    //Region  of block owned by this rank
		Real dx;                //Grid size
		int Nf[6];              //Part of the block in grid nodes
		int NReserve; //Maximum size that resrved for particleList
		const Block& currentBlock;   //The block to which the particle list is mapped
		ParticleShapeDiscriptor particleShape; //Type of particles
		bool SphereParallepipedIntersection(Real* xpos,Real radius);
		bool hasExtraInputVariables;		//Additional input parameters (i.e. T)
		bool hasExtraOutputVariables;       //Output parameters other than (Fd, Td)
		int extraInputSize;
		int extraOutputSize;
};


using ListBlockParticles = std::map<int, BlockParticles>;
extern ListBlockParticles BlockParticleList;
#endif
