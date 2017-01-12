#ifndef __SPH_PARTICLE_H__
#define __SPH_PARTICLE_H__

#include <vector>
#include "sph_type.h"

class Particle
{
public:
	uint id;
	vec3 pos;
	vec3 vel = vec3(0.0f, 0.0f, 0.0f);

	vec3 acc;
	vec3 ev;
	vec3 colour;

	float restDensity;
	float interp_density;
	float pressure;
	float restMass;
	float viscosity;
	uint phase;
	vec3 conv_mom_change;
	vec3 div_viscosity_tensor;
	std::vector<vec3> pressureGradientk;
	vec3 pressureGradientm;
	std::vector<vec3> diffusionGradient;

	float sum_mass_density = 0.0f;
	vec3 sum_mass_pressure = vec3(0.0f, 0.0f, 0.0f);
	vec3 sum_mass_diffusion = vec3(0.0f, 0.0f, 0.0f);

	float surf_norm;

	std::vector<vec3> driftVelocity;
	std::vector<float> volumeFraction;
	std::vector<float> massFraction;

	Particle *next;
};


#endif
