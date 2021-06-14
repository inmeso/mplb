/*
 * particle_shape.cpp
 *
 *  Created on: May 6, 2021
 *      Author: jpd38567
 */



#include "particle_shape.h"


ParticleShape::ParticleShape(Real Rp, ParticleShapeDiscriptor particleShape) {

	Rparticle = Rp;
	typeParticle = particleShape;

}

ParticleShapeQuadratic::ParticleShapeQuadratic(Real Rp, ParticleShapeDiscriptor particleShape,
		std::vector<Real> shape) :	ParticleShape(Rp, particleShape) {

	particleParameters = shape;

}
