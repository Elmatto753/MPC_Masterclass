/** File:		sph_system.h
 ** Author:		Dongli Zhang
 ** Contact:	dongli.zhang0129@gmail.com
 **
 ** Copyright (C) Dongli Zhang 2013
 **
 ** This program is free software;  you can redistribute it and/or modify
 ** it under the terms of the GNU General Public License as published by
 ** the Free Software Foundation; either version 2 of the License, or
 ** (at your option) any later version.
 **
 ** This program is distributed in the hope that it will be useful,
 ** but WITHOUT ANY WARRANTY;  without even the implied warranty of
 ** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See
 ** the GNU General Public License for more details.
 **
 ** You should have received a copy of the GNU General Public License
 ** along with this program;  if not, write to the Free Software 
 ** Foundation, 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

#ifndef __SPHSYSTEM_H__
#define __SPHSYSTEM_H__

#include "sph_type.h"
#include "sph_phase.h"
#include "sph_particle.h"
#include <math.h>



class SPHSystem
{
public:
	uint max_particle;
	uint num_particle;

	float kernel;
	float mass;

	vec3 world_size;
	float cell_size;
	uint3 grid_size;
	uint tot_cell;

	vec3 gravity;
	float wall_damping;
	float rest_density;

	float gas_constant;
	float viscosity;
	float time_step;
	float surf_norm;
	float surf_coe;

	float poly6_value;
	float spiky_value;
	float visco_value;

	float grad_poly6;
	float lplc_poly6;

	float kernel_2;
	float self_dens;
	float self_lplc_color;

	vec3 prevVel;
	vec3 drift_velocity;
	vec3 acceleration;


	Phase *phase1;
	Phase *phase2;

	Particle *mem;
	Particle **cell;

	std::vector<Phase *> phaseList;

	uint num_phases;

	uint sys_running;

public:
	SPHSystem();
	~SPHSystem();
	void animation();
	void init_system();
	void add_particle(Phase *phase, vec3 pos, vec3 vel);

private:
	void build_table();
	void comp_dens_pres();
	void comp_force_adv();
	void advection();
	void adv_vol_frac();
	void calc_acceleration();

	float interp_mix_dens();
	float poly6(Particle *p, Particle *np);
	vec3 spiky(Particle *p, Particle *np);
	float visco(Particle *p, Particle *np);
	float length(vec3 vector);
	vec3 normalize(vec3 vector);
	void driftVelocity();

private:
	int3 calc_cell_pos(vec3 p);
	uint calc_cell_hash(int3 cell_pos);
};

#endif
