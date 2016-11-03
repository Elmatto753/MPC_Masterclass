#ifndef __SPH_PARTICLE_H__
#define __SPH_PARTICLE_H__

#include "sph_type.h"

class Particle
{
public:
	uint id;
	vec3 pos;
	vec3 vel;

	vec3 acc;
	vec3 ev;

	float dens;
	float pres;
	float mass;

	float surf_norm;

	Particle *next;
};


#endif
