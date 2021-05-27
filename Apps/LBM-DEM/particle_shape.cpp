/*
 * particle_shape.cpp
 *
 *  Created on: May 6, 2021
 *      Author: jpd38567
 */



#include "particle_shape.h"


ParticleShape::ParticleShape(Real Rp, std::string name) {

	Rparticle = Rp;
	typeParticle = name;

}

void ParticleShape::updateShape(Real radius, std::vector<Real> shape) {

	Rparticle = radius;
}

Real ParticleShape::getArea() {


	return 4.0/3.0 * PI * Rparticle * Rparticle * Rparticle;

}

Real ParticleShape::getSurface() {

	return 4.0 * PI * Rparticle * Rparticle;
}

Real ParticleShape::getEquivalentRadius() {

	return Rparticle;
}

ParticleShapeQuadratic::ParticleShapeQuadratic(Real Rp, std::string name,
		std::vector<Real> shape) :	ParticleShape(Rp, name) {

	particleParameters = shape;

}
