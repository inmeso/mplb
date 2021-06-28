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

/*! @brief  Classes for handling and storing different particle shapes.
 * @author C. Tsigginos
 */

#ifndef PARTICLE_SHAPE_H_
#define PARTICLE_SHAPE_H_

#include "type.h"
#include <string>
#include <vector>
#include "dem_data.h"
class ParticleShape{

	public:
		Real Rparticle;
		ParticleShapeDiscriptor typeParticle;
		ParticleShape(Real Rp, ParticleShapeDiscriptor particleShape);
		virtual ~ParticleShape() { };
		virtual void Rotate() { };
		virtual void UpdateShape(Real radius, std::vector<Real> shape) { Rparticle = radius; };
		virtual Real GetArea() { return 4.0/3.0 * PI * Rparticle * Rparticle * Rparticle;};
		virtual Real GetSurface() {return 4.0 * PI * Rparticle * Rparticle;}
		virtual Real GetEquivalentRadius() {return Rparticle; };


};

class ParticleShapeQuadratic : public ParticleShape {

	public:
		std::vector<Real> particleParameters;
		ParticleShapeQuadratic(Real Rp, ParticleShapeDiscriptor particleShape,
				std::vector<Real> shape);
		virtual void Rotate() { };
		virtual void UpdateShape(Real radius, std::vector<Real> shape) { };
		virtual Real GetArea() {return 0;};
		virtual Real GetSurface() {return 0;};
		virtual Real GetEquivalentRadius() {return 0;};

};

class ParticleShapeMesh : public ParticleShape {

	public:
		std::vector<Real> gridPoints;
		ParticleShapeMesh(Real Rp, ParticleShapeDiscriptor particleShape, std::vector<Real> shape);
		virtual void UpdateShape(Real radius, std::vector<Real> shape) { };
		virtual void Rotate() { };
		virtual Real GetArea() {return 0;}
		virtual Real GetSurface() {return 0;};
		virtual Real GetEquivalentRadius() {return 0;};


};




#endif /* APPS_LBM_DEM_PARTICLE_SHAPE_H_ */
