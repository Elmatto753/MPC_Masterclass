/** File:		sph_system.cpp
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

#include "sph_system.h"
#include "sph_header.h"
#include <iostream>

SPHSystem::SPHSystem()
{
	max_particle=300000;
	num_particle=0;

	kernel=0.04f;
	mass=0.02f;

	world_size.m_x=0.2f;
	world_size.m_y=0.2f;
	world_size.m_z=0.2f;
	cell_size=0.005;
	grid_size.m_x=(uint)ceil(world_size.m_x/cell_size);
	grid_size.m_y=(uint)ceil(world_size.m_y/cell_size);
	grid_size.m_z=(uint)ceil(world_size.m_z/cell_size);
	tot_cell=grid_size.m_x*grid_size.m_y*grid_size.m_z;

	gravity.m_x=0.0f;
	gravity.m_y=-6.8f;
	gravity.m_z=0.0f;
	wall_damping=-0.5f;
	rest_density=1000.0f;
	gas_constant=1.0f;
	viscosity=6.5f;
	time_step=0.003f;
	surf_norm=6.0f;
	surf_coe=0.1f;

	poly6_value=315.0f/(64.0f * PI * pow(kernel, 9));;
	spiky_value=-45.0f/(PI * pow(kernel, 6));
	visco_value=45.0f/(PI * pow(kernel, 6));

	grad_poly6=-945/(32 * PI * pow(kernel, 9));
	lplc_poly6=-945/(8 * PI * pow(kernel, 9));

	kernel_2=kernel*kernel;
	self_dens=mass*poly6_value*pow(kernel, 6);
	self_lplc_color=lplc_poly6*mass*kernel_2*(0-3/4*kernel_2);

	mem=(Particle *)malloc(sizeof(Particle)*max_particle);

	cell=(Particle **)malloc(sizeof(Particle *)*tot_cell);

	phase1 = new Phase;
	phaseList.push_back(phase1);
	num_phases++;
	phase2 = new Phase;
	phaseList.push_back(phase2);
	num_phases++;
	phase1->phaseColour = vec3(1.0f, 0.0f, 0.0f);
	phase2->phaseColour = vec3(0.0f, 0.0f, 1.0f);
	phase1->restDensity = 1000.0f;
	phase2->restDensity = 500.0f;
	phase1->restMass = 0.01f;
	phase2->restMass = 0.05f;

	sys_running=0;

	printf("Initialize SPH:\n");
	printf("World Width : %f\n", world_size.m_x);
	printf("World Height: %f\n", world_size.m_y);
	printf("World Length: %f\n", world_size.m_z);
	printf("Cell Size  : %f\n", cell_size);
	printf("Grid Width : %u\n", grid_size.m_x);
	printf("Grid Height: %u\n", grid_size.m_y);
	printf("Grid Length: %u\n", grid_size.m_z);
	printf("Total Cell : %u\n", tot_cell);
	printf("Poly6 Kernel: %f\n", poly6_value);
	printf("Spiky Kernel: %f\n", spiky_value);
	printf("Visco Kernel: %f\n", visco_value);
	printf("Self Density: %f\n", self_dens);
}

SPHSystem::~SPHSystem()
{
	free(mem);
	free(cell);
}

void SPHSystem::animation()
{
	if(sys_running == 0)
	{
		return;
	}

	build_table();

	driftVelocity();

	adv_vol_frac();

	comp_dens_pres();

	comp_force_adv();

	advection();

}

void SPHSystem::init_system()
{
	vec3 pos;
	vec3 vel;

	vel.m_x=0.0f;
	vel.m_y=0.0f;
	vel.m_z=0.0f;

	for(pos.m_x=world_size.m_x*0.0f; pos.m_x<world_size.m_x*0.2f; pos.m_x+=(kernel*0.1f))
	{
		for(pos.m_y=world_size.m_y*0.0f; pos.m_y<world_size.m_y*0.2f; pos.m_y+=(kernel*0.1f))
		{
			for(pos.m_z=world_size.m_z*0.0f; pos.m_z<world_size.m_z*0.2f; pos.m_z+=(kernel*0.1f))
			{
				add_particle(phase1, pos, vel);
			}
		}
	}

	for(pos.m_x=world_size.m_x*0.0f; pos.m_x<world_size.m_x*0.2f; pos.m_x+=(kernel*0.1f))
	{
		for(pos.m_y=world_size.m_y*0.0f; pos.m_y<world_size.m_y*0.2f; pos.m_y+=(kernel*0.1f))
		{
			for(pos.m_z=world_size.m_z*0.0f; pos.m_z<world_size.m_z*0.2f; pos.m_z+=(kernel*0.1f))
			{

				add_particle(phase2, pos+vec3(0.1f, 0.0f, 0.0f), vel);
			}
		}
	}

	std::cout<<phase1->particleList.size() << "\n";
	std::cout<< phaseList.size() << "\n";
	std::cout<< phaseList[0]->particleList.size() << "\n";;
	std::cout<<"Init Particle: "<<num_particle<<"\n";
}

void SPHSystem::add_particle(Phase *phase, vec3 pos, vec3 vel)
{
//  if(phase->numParticle != 0)
//  {
	Particle *p=&(mem[num_particle]);

	//std::cout<<num_particle<<"\n";

	p->id=phase->numParticle;

	p->pos=pos;
	p->vel=vel;

	p->acc.m_x=0.0f;
	p->acc.m_y=0.0f;
	p->acc.m_z=0.0f;
	p->ev.m_x=0.0f;
	p->ev.m_y=0.0f;
	p->ev.m_z=0.0f;

	p->restDensity=phase->restDensity;
	p->pressure=0.0f;
	p->restMass=phase->restMass;
	p->colour = phase->phaseColour;
	p->interp_density = 0.01f;

	for(uint i = 0; i<phaseList.size(); i++)
	{
	  p->driftVelocity.push_back(vec3(0.0f, 0.0f, 0.0f));
	}

	p->next=NULL;

	for(uint i = 0; i<num_phases; i++)
	{
//	  std::cout<<"i = "<<num_phases<<"\n";
	  if(phaseList[i] != phase)
	  {
	    p->volumeFraction.push_back(0.01f);
	  }

	  else
	  {
	    p->volumeFraction.push_back(0.99f);
	  }


	}

	phase->particleList.push_back(p);
//	Particle test = (phase->particleList[phase->numParticle]);
//	std::cout<<"stuff: "<<test.pos.m_x<<"\n";
//  }
//  else
//  {
//    std::cerr<<"hi";
//    Particle *o = new Particle;
//    phase->particleList.push_back(o);
//    o->id=phase->numParticle;

//    o->pos=pos;
//    o->vel=vel;

//    o->acc.m_x=0.0f;
//    o->acc.m_y=0.0f;
//    o->acc.m_z=0.0f;
//    o->ev.m_x=0.0f;
//    o->ev.m_y=0.0f;
//    o->ev.m_z=0.0f;

//    o->restDensity=rest_density;
//    o->pressure=0.0f;
//    o->restMass=0.01f;
//    o->colour = phase->phaseColour;

//    o->next=NULL;
//    phase->particleList.push_back(o);


//  }


	phase->numParticle++;
	num_particle++;
	//std::cout<<"num: "<<num_particle<<"\n";
}

void SPHSystem::build_table()
{
	Particle *p;
	uint hash;

	for(uint i=0; i<tot_cell; i++)
	{
		cell[i]=NULL;
	}

	for(uint i=0; i<phaseList.size(); i++)
	{
	  for(uint j=0; j<phaseList[i]->numParticle; j++)
	  {
	    //std::cerr<<i<<"\n";
		  p=phaseList[i]->particleList[j];
		  //std::cerr<<p->pos.m_x<<"\n";
		  hash=calc_cell_hash(calc_cell_pos(p->pos));

                  if(cell[hash] == NULL)
                  {
                          p->next=NULL;
                          cell[hash]=p;
                  }
                  else
                  {
                          p->next=cell[hash];
                          cell[hash]=p;
                  }
          }
        }

}

void SPHSystem::comp_dens_pres()
{
	Particle *p;
	Particle *np;

	int3 cell_pos;
	int3 near_pos;
	uint hash;

	vec3 rel_pos;
	float r2;

        for(uint i = 0; i<phaseList.size(); i++)
        {
          for(uint j = 0; j<phaseList[i]->particleList.size(); j++)
          {
                p=phaseList[i]->particleList[j];
		cell_pos=calc_cell_pos(p->pos);

		//p->restDensity=0.0f;
//		p->pressure=0.0f;

		for(int x=-1; x<=1; x++)
		{
			for(int y=-1; y<=1; y++)
			{
				for(int z=-1; z<=1; z++)
				{
					near_pos.m_x=cell_pos.m_x+x;
					near_pos.m_y=cell_pos.m_y+y;
					near_pos.m_z=cell_pos.m_z+z;
					hash=calc_cell_hash(near_pos);

					if(hash == 0xffffffff)
					{
						continue;
					}

					np=cell[hash];
					int neighbornum = 0;
					while(np != NULL)
					{
						rel_pos=p->pos-np->pos;
						r2=rel_pos.m_x*rel_pos.m_x+rel_pos.m_y*rel_pos.m_y+rel_pos.m_z*rel_pos.m_z;

						if(r2<INF || r2>=kernel_2)
						{
							np=np->next;
							continue;
						}

						p->restDensity=p->restDensity + p->restMass * poly6(p, np);

						p->restDensity=p->restDensity+self_dens;
						p->interp_density += p->restMass * poly6(p, np);
						p->pressure=(pow(p->interp_density / p->restDensity, 7) - 1) * ((gas_constant * p->restDensity) / 7);
						//std::cout<<p->pressure<<"\n";

						np=np->next;

					}
				}
			}
		}



	  }

	}
}

void SPHSystem::comp_force_adv()
{
	Particle *p;
	Particle *np;

	int3 cell_pos;
	int3 near_pos;
	uint hash;

	vec3 rel_pos;
	vec3 rel_vel;

	float r2;
	float r;
	float kernel_r;
	float V;

	float pres_kernel;
	float visc_kernel;
	float temp_force;

	vec3 grad_color;
	float lplc_color;

	for(uint i=0; i<phaseList.size(); i++)
	{
	  for(uint j=0; j<phaseList[i]->numParticle; j++)
	  {
		  p=phaseList[i]->particleList[j];
		  cell_pos=calc_cell_pos(p->pos);

                  p->acc.m_x=0.0f;
                  p->acc.m_y=0.0f;
                  p->acc.m_z=0.0f;

                  grad_color.m_x=0.0f;
                  grad_color.m_y=0.0f;
                  grad_color.m_z=0.0f;
                  lplc_color=0.0f;

                  for(int x=-1; x<=1; x++)
                  {
                          for(int y=-1; y<=1; y++)
                          {
                                  for(int z=-1; z<=1; z++)
                                  {
                                          near_pos.m_x=cell_pos.m_x+x;
                                          near_pos.m_y=cell_pos.m_y+y;
                                          near_pos.m_z=cell_pos.m_z+z;
                                          hash=calc_cell_hash(near_pos);

                                          if(hash == 0xffffffff)
                                          {
                                                  continue;
                                          }

                                          np=cell[hash];
                                          while(np != NULL)
                                          {
                                                  rel_pos= p->pos-np->pos;
                                                  //rel_pos.m_y=p->pos.m_y-np->pos.m_y;
                                                  //rel_pos.m_z=p->pos.m_z-np->pos.m_z;
                                                  r2=rel_pos.m_x*rel_pos.m_x+rel_pos.m_y*rel_pos.m_y+rel_pos.m_z*rel_pos.m_z;

                                                  if(r2 < kernel_2 && r2 > INF)
                                                  {
                                                          r=sqrt(r2);
                                                          V=mass/np->restDensity/2;
                                                          kernel_r=kernel-r;

                                                          pres_kernel=spiky_value * kernel_r * kernel_r;

                                                          temp_force=V * (p->pressure+np->pressure) * pres_kernel;
                                                          p->acc=p->acc-rel_pos*temp_force/r;



                                                          rel_vel=np->ev-p->ev;

                                                          visc_kernel=visco_value*(kernel-r);
                                                          temp_force=V * viscosity * visc_kernel;
                                                          p->acc=p->acc + rel_vel*temp_force;


                                                          float temp=(-1) * grad_poly6 * V * pow(kernel_2-r2, 2);
                                                          grad_color += rel_pos * temp;
                                                          lplc_color += lplc_poly6 * V * (kernel_2-r2) * (r2-3/4*(kernel_2-r2));
                                                  }

                                                  np=np->next;
                                          }
                                  }
                          }
                  }

                  lplc_color+=self_lplc_color/p->restDensity;
                  p->surf_norm=sqrt(grad_color.m_x*grad_color.m_x+grad_color.m_y*grad_color.m_y+grad_color.m_z*grad_color.m_z);

                  if(p->surf_norm > surf_norm)
                  {
                          p->acc.m_x+=surf_coe * lplc_color * grad_color.m_x / p->surf_norm;
                          p->acc.m_y+=surf_coe * lplc_color * grad_color.m_y / p->surf_norm;
                          p->acc.m_z+=surf_coe * lplc_color * grad_color.m_z / p->surf_norm;


                  }

          }
        }

}

void SPHSystem::advection()
{
	Particle *p;
	for(uint i=0; i<phaseList.size(); i++)
	{
	  for(uint j=0; j<phaseList[i]->numParticle; j++)
	  {
		  p=phaseList[i]->particleList[j];


                  p->vel=p->vel+p->acc*time_step/-(p->interp_density/100)+gravity*time_step;
  //		p->vel.m_y=p->vel.m_y+p->acc.m_y*time_step/p->restDensity+gravity.m_y*time_step;
  //		p->vel.m_z=p->vel.m_z+p->acc.m_z*time_step/p->restDensity+gravity.m_z*time_step;

                  p->pos=p->pos+p->vel*time_step;
  //		p->pos.m_y=p->pos.m_y+p->vel.m_y*time_step;
  //		p->pos.m_z=p->pos.m_z+p->vel.m_z*time_step;

                  if(p->pos.m_x >= world_size.m_x-BOUNDARY)
                  {
                          p->vel.m_x=p->vel.m_x*wall_damping;
                          p->pos.m_x=world_size.m_x-BOUNDARY;
                  }

                  if(p->pos.m_x < 0.0f)
                  {
                          p->vel.m_x=p->vel.m_x*wall_damping;
                          p->pos.m_x=0.0f;
                  }

                  if(p->pos.m_y >= world_size.m_y-BOUNDARY)
                  {
                          p->vel.m_y=p->vel.m_y*wall_damping;
                          p->pos.m_y=world_size.m_y-BOUNDARY;
                  }

                  if(p->pos.m_y < 0.0f)
                  {
                          p->vel.m_y=p->vel.m_y*wall_damping;
                          p->pos.m_y=0.0f;
                  }

                  if(p->pos.m_z >= world_size.m_z-BOUNDARY)
                  {
                          p->vel.m_z=p->vel.m_z*wall_damping;
                          p->pos.m_z=world_size.m_z-BOUNDARY;
                  }

                  if(p->pos.m_z < 0.0f)
                  {
                          p->vel.m_z=p->vel.m_z*wall_damping;
                          p->pos.m_z=0.0f;
                  }

                  p->ev=(p->ev+p->vel)/2;
  //		p->ev.m_y=(p->ev.m_y+p->vel.m_y)/2;
  //		p->ev.m_z=(p->ev.m_z+p->vel.m_z)/2;
          }
        }
}

void SPHSystem::adv_vol_frac()
{
  Particle *p;
  for(uint i =0; i<phaseList.size(); i++)
  {
    for(uint j = 0; j<phaseList[i]->particleList.size(); j++)
    {
      p = phaseList[i]->particleList[j];
      p->colour = vec3(0.0f, 0.0f, 0.0f);

      p->volumeFraction[i] += p->vel.dot(p->vel) * -p->volumeFraction[i] - (sqrt(p->driftVelocity[i].dot(p->driftVelocity[i])) * p->volumeFraction[i]) - p->vel.dot(p->vel) * p->volumeFraction[i];
    }
  }

  for(uint i =0; i<phaseList.size(); i++)
  {
    for(uint j = 0; j<phaseList[i]->particleList.size(); j++)
    {
      p = phaseList[i]->particleList[j];

      float sum_vol_frac = 0.0f;
      for(uint k = 0; k<phaseList.size(); k++)
      {
        sum_vol_frac += p->volumeFraction[k];
      }
      if(sum_vol_frac != 1.0f)
      {
        float divisor = 1.0f/sum_vol_frac;
        for(uint l = 0; l<phaseList.size(); l++)
        {
          p->volumeFraction[l] *= divisor;
//          std::cout<<p->volumeFraction[l]<<"\n";

        }
      }
      for(uint m = 0; m<phaseList.size(); m++)
      {
        p->colour += phaseList[m]->phaseColour * p->volumeFraction[m];
      }
    }
  }

}

int3 SPHSystem::calc_cell_pos(vec3 p)
{
	int3 cell_pos;
	cell_pos.m_x = int(floor((p.m_x) / cell_size));
	cell_pos.m_y = int(floor((p.m_y) / cell_size));
	cell_pos.m_z = int(floor((p.m_z) / cell_size));

    return cell_pos;
}

uint SPHSystem::calc_cell_hash(int3 cell_pos)
{
	if(cell_pos.m_x<0 || cell_pos.m_x>=(int)grid_size.m_x || cell_pos.m_y<0 || cell_pos.m_y>=(int)grid_size.m_y || cell_pos.m_z<0 || cell_pos.m_z>=(int)grid_size.m_z)
	{
		return (uint)0xffffffff;
	}

        cell_pos.m_x = cell_pos.m_x & (grid_size.m_x-1);
        cell_pos.m_y = cell_pos.m_y & (grid_size.m_y-1);
        cell_pos.m_z = cell_pos.m_z & (grid_size.m_z-1);

        return ((uint)(cell_pos.m_z))*grid_size.m_y*grid_size.m_x + ((uint)(cell_pos.m_y))*grid_size.m_x + (uint)(cell_pos.m_x);
}

float SPHSystem::poly6(Particle *p, Particle *np)
{
  vec3 r = p->pos - np->pos;
  float r2 = r.dot(r);
  if(sqrt(r2)>kernel)
  {
    return 0;
  }
  else
  {
    return poly6_value * pow(kernel_2 - r2, 3);
  }
}

vec3 SPHSystem::spiky(Particle *p, Particle *np)
{
  vec3 r = p->pos - np->pos;
  float rlen= length(r);

  if(rlen>kernel)
  {
    return vec3(0.0f, 0.0f, 0.0f);
  }
  else
  {
    vec3 normal = r/rlen;
    return normal * spiky_value * pow(kernel - rlen, 2) ;
  }

}

float SPHSystem::visco(Particle *p, Particle *np)
{
  vec3 r = p->pos - np->pos;
  float rlen = length(r);

  if(rlen>kernel)
  {
    return 0.0f;
  }
  else
  {
    return visco_value * (kernel - rlen);
  }
}

float SPHSystem::length(vec3 vector)
{
  float X = vector.m_x;
  float Y = vector.m_y;
  float Z = vector.m_z;
  return sqrt( X * X + Y * Y + Z * Z );
}

vec3 SPHSystem::normalize(vec3 vector)
{
  vec3 normal;
  float vecLength = length(vector);

  if(vecLength != 0)
  {
    normal.m_x = vector.m_x / vecLength;
    normal.m_y = vector.m_y / vecLength;
    normal.m_z = vector.m_z / vecLength;
  }

  return normal;

}

void SPHSystem::driftVelocity()
{
  Particle *p;
  Particle *np;
//  vec3 pressureGradient;
  int3 cell_pos;
  int3 near_pos;
  uint hash;

  for(uint i=0; i<phaseList.size(); i++)
  {
    for(uint j=0; j<phaseList[i]->numParticle; j++)
    {
      p=phaseList[i]->particleList[j];
      p->viscosity = 0.0f;
      float sum_mass_dens = 0;
      vec3 sum_drift_vel = vec3(0.0f, 0.0f, 0.0f);

      for(int x=-1; x<=1; x++)
      {
        for(int y=-1; y<=1; y++)
        {
          for(int z=-1; z<=1; z++)
          {
            near_pos.m_x=cell_pos.m_x+x;
            near_pos.m_y=cell_pos.m_y+y;
            near_pos.m_z=cell_pos.m_z+z;
            hash=calc_cell_hash(near_pos);

            if(hash == 0xffffffff)
            {
              continue;
            }

            np=cell[hash];

            while(np != NULL)
            {

              p->pressureGradientk += spiky(p, np) * p->restMass/p->restDensity * (np->pressure - p->pressure);
              p->diffusionGradient += spiky(p, np) * p->restMass/p->restDensity * (np->volumeFraction[i] - p->volumeFraction[i]);

              p->pressureGradientm += spiky(p, np) * np->restMass * (np->pressure + p->pressure)/(2*np->interp_density);

              sum_mass_dens += np->restMass/np->interp_density;

              for(uint k = 0; k<phaseList.size(); k++)
              {
                sum_drift_vel += (np->driftVelocity[k] * p->volumeFraction[k] * (np->driftVelocity[k].dot(spiky(p, np))) + p->driftVelocity[k] * p->volumeFraction[k] * (p->driftVelocity[k].dot(spiky(p, np))) ) * p->restDensity;
                np->viscosity += np->volumeFraction[k] * phaseList[k]->phaseViscosity;
              }

              p->div_viscosity_tensor += (np->vel - p->vel) * (np->restMass/p->interp_density) * (p->viscosity + np->viscosity) * ((vec3(np->pos - p->pos).dot(spiky(p, np)))/pow(length(vec3(np->pos - p->pos)), 2));
            }
          }
        }
      }
      p->conv_mom_change = sum_drift_vel * - sum_mass_dens;
    }
  }

  for(uint i=0; i<phaseList.size(); i++)
  {
    for(uint j=0; j<phaseList[i]->numParticle; j++)
    {
      p=phaseList[i]->particleList[j];

      for(uint k = 0; k<phaseList.size(); k++)
      {
        p->massFraction.push_back((p->volumeFraction[k] * p->restDensity)/p->interp_density);
        p->sum_mass_density += p->massFraction[k] * p->restDensity;

        p->acc = (p->pressureGradientm/p->restDensity) - (p->conv_mom_change/p->restDensity) - (p->div_viscosity_tensor/p->restDensity);
      }

    }

  }


  for(uint i=0; i<phaseList.size(); i++)
  {
    for(uint j=0; j<phaseList[i]->numParticle; j++)
    {
      p=phaseList[i]->particleList[j];

      p->driftVelocity.push_back(p->acc * 1 * (p->restDensity - p->sum_mass_density));

    //pressureGradient += p->mass * ()
    }

  }
 // drift_velocity = 1 * (rest_density - mass * rest_density) * acceleration - 1 * (normalize(p->pres) );
  prevVel = drift_velocity;
}


