/** File:		sph_system.cpp
 ** Author:		Dongli Zhang
 ** Contact:	dongli.m_zhang0129@gmail.com
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
#include <cassert>

SPHSystem::SPHSystem()
{
	max_particle=300000;
	num_particle=0;

	kernel=0.04f;
//	mass=0.02f;

	world_size.m_x=0.64f;
	world_size.m_y=0.64f;
	world_size.m_z=0.64f;
	cell_size=kernel;
	grid_size.m_x=(uint)ceil(world_size.m_x/cell_size);
	grid_size.m_y=(uint)ceil(world_size.m_y/cell_size);
	grid_size.m_z=(uint)ceil(world_size.m_z/cell_size);
	tot_cell=grid_size.m_x*grid_size.m_y*grid_size.m_z;

	gravity.m_x=0.0f;
	gravity.m_y=-9.8f;
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
	phase2->restDensity = 1000.0f;
	phase1->restMass = 0.01f;
	phase2->restMass = 0.01f;
	phase1->phaseViscosity = 2.5f;
	phase2->phaseViscosity = 2.5f;

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

	// Paper step 1
	driftVelocity();

	// Paper step 2 & 3
//	adv_vol_frac();

	// Paper step 4
	calc_acceleration();

//	comp_dens_pres();

//	comp_force_adv();

	// Paper step 5
	advection();

}

void SPHSystem::init_system()
{
	vec3 pos;
	vec3 vel;

	vel.m_x=0.0f;
	vel.m_y=0.0f;
	vel.m_z=0.0f;

//	for(pos.m_x=0.0f; pos.m_x<0.64f*0.1f; pos.m_x+=(0.05f*0.2f))
//	{
//		for(pos.m_y=0.0f; pos.m_y<0.64f*0.1f; pos.m_y+=(0.05f*0.2f))
//		{
//			for(pos.m_z=0.0f; pos.m_z<0.64f*0.1f; pos.m_z+=(0.05f*0.2f))
//			{
//				add_particle(phase1, pos+vec3(0.0f, 0.01f, 0.0f), vel);
//			}
//		}
//	}

//	for(pos.m_x=0.0f; pos.m_x<0.64f*0.1f; pos.m_x+=(0.05f*0.2f))
//	{
//		for(pos.m_y=0.0f; pos.m_y<0.64f*0.1f; pos.m_y+=(0.05f*0.2f))
//		{
//			for(pos.m_z=0.0f; pos.m_z<0.64f*0.1f; pos.m_z+=(0.05f*0.2f))
//			{

//				add_particle(phase2, pos+vec3(0.08f, 0.01f, 0.0f), vel);
//			}
//		}
//	}

	for(pos.m_x=world_size.m_x*0.0f; pos.m_x<world_size.m_x*0.3f; pos.m_x+=(0.04*world_size.m_x))
	{
		for(pos.m_y=world_size.m_y*0.0f; pos.m_y<world_size.m_y*0.3f; pos.m_y+=(0.04*world_size.m_y))
		{
			for(pos.m_z=world_size.m_z*0.0f; pos.m_z<world_size.m_z*0.3f; pos.m_z+=(0.04*world_size.m_z))
			{
				add_particle(phase1, pos+vec3(0.0f, 0.01f, 0.0f), vel);
			}
		}
	}

	for(pos.m_x=world_size.m_x*0.0f; pos.m_x<world_size.m_x*0.3f; pos.m_x+=(0.04*world_size.m_x))
	{
		for(pos.m_y=world_size.m_y*0.0f; pos.m_y<world_size.m_y*0.3f; pos.m_y+=(0.04*world_size.m_y))
		{
			for(pos.m_z=world_size.m_z*0.0f; pos.m_z<world_size.m_z*0.3f; pos.m_z+=(0.04*world_size.m_z))
			{
				add_particle(phase2, pos+vec3(world_size.m_x/2.f, 0.01f, 0.0f), vel);
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

	p->acc = gravity;
	p->ev.m_x=0.0f;
	p->ev.m_y=0.0f;
	p->ev.m_z=0.0f;

	p->restDensity=phase->restDensity;
	p->pressure=0.0f;
	p->restMass=phase->restMass;
	p->colour = phase->phaseColour;
	p->interp_density = 0.01f;
	if(phase == phase1)
	{
	  p->phase = 0;
	}

	if(phase == phase2)
	{
	  p->phase = 1;
	}

	for(uint i = 0; i<phaseList.size(); i++)
	{
	  p->driftVelocity.push_back(vec3(0.0f, 0.0f, 0.0f));
	  p->massFraction.push_back(0.0f);
	  p->diffusionGradient.push_back(vec3(0.0f, 0.0f, 0.0f));
	  p->pressureGradientk.push_back(vec3(0.0f, 0.0f, 0.0f));
	}

	p->next=NULL;

	for(uint i = 0; i<num_phases; i++)
	{
//	  std::cout<<"i = "<<num_phases<<"\n";
	  if(phaseList[i] != phase)
	  {
	    p->volumeFraction.push_back(0.f);
	  }

	  else
	  {
	    p->volumeFraction.push_back(1.f);
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


                  if(hash == 4294967295)
                  {
                    hash = 0;
                  }


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

					np=p->next;
					int neighbornum = 0;
					while(np != NULL)
					{
						rel_pos=np->pos-p->pos;
						r2=rel_pos.m_x*rel_pos.m_x+rel_pos.m_y*rel_pos.m_y+rel_pos.m_z*rel_pos.m_z;

						if(r2<INF || r2>=kernel_2)
						{
							np=np->next;
							continue;
						}

						//p->restDensity=p->restDensity + p->restMass * poly6(p, np);

						//p->interp_density += p->restMass * poly6(p, np);
						//std::cout<<p->pressure<<"\n";

						np=np->next;

					}
				}
			}
		}

		for(uint k=0; k<phaseList.size(); k++)
		{
		  p->restDensity += p->volumeFraction[k] * phaseList[k]->restDensity;

		}
		//p->pressure=(pow(p->interp_density / p->restDensity, 7) - 1) * ((gas_constant * p->restDensity) / 7);
		//p->pressure = 1 * (p->interp_density - p->restDensity);


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

                                          np=p->next;
                                          while(np != NULL)
                                          {
                                                  rel_pos= p->pos-np->pos;
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
                                                          temp_force=V * p->viscosity * visc_kernel;
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


//                  p->vel=p->vel+p->acc*time_step/(p->interp_density)+gravity*time_step;
//                  p->vel=p->vel+p->acc*time_step;

                  vec3 oldpos = p->pos;
                  p->pos=p->pos+p->vel*time_step + p->acc*time_step*time_step*0.5;
                  p->vel= (p->pos - oldpos) / time_step;

                  if(p->pos.m_x >= world_size.m_x-BOUNDARY)
                  {
                          p->vel.m_x=p->vel.m_x*wall_damping;
                          p->pos.m_x=world_size.m_x-BOUNDARY - 0.00f;
                  }

                  if(p->pos.m_x < 0.0f)
                  {
                          p->vel.m_x=p->vel.m_x*wall_damping;
                          p->pos.m_x=0.00f;
                  }

                  if(p->pos.m_y >= world_size.m_y-BOUNDARY)
                  {
                          p->vel.m_y=p->vel.m_y*wall_damping;
                          p->pos.m_y=world_size.m_y-BOUNDARY - 0.00f;
                  }

                  if(p->pos.m_y < 0.0f)
                  {
                          p->vel.m_y=p->vel.m_y*wall_damping;
                          p->pos.m_y=0.00f;
                  }

                  if(p->pos.m_z >= world_size.m_z-BOUNDARY)
                  {
                          p->vel.m_z=p->vel.m_z*wall_damping;
                          p->pos.m_z=world_size.m_z-BOUNDARY - 0.00f;
                  }

                  if(p->pos.m_z < 0.0f)
                  {
                          p->vel.m_z=p->vel.m_z*wall_damping;
                          p->pos.m_z=0.00f;
                  }

                  p->ev=(p->ev+p->vel)/2;
          }
        }
}

void SPHSystem::adv_vol_frac()
{
  Particle *p;
  Particle *np;
  int3 cell_pos;
  int3 near_pos;
  uint hash;
  float pressureChange = 0.0f;
  for(uint i =0; i<phaseList.size(); i++)
  {
    for(uint j = 0; j<phaseList[i]->particleList.size(); j++)
    {
      p = phaseList[i]->particleList[j];
      cell_pos=calc_cell_pos(p->pos);
//      p->colour = vec3(0.0f, 0.0f, 0.0f);

      for(uint k = 0; k < phaseList.size(); k++)
      {
        //p->volumeFraction[k] += (length(p->vel) * -p->volumeFraction[k] - (length((p->driftVelocity[k]) * p->volumeFraction[k])) - length(p->vel) * p->volumeFraction[k]) / time_step;
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

              np = cell[hash];
              while(np != NULL)
              {
                if(np == p)
                {
                  np = np->next;
                  continue;
                }
                if(np->interp_density == 0.0f)
                {
                  np = np->next;
                  continue;
                }
                p->volumeFraction[k] += ((np->restMass)/(np->interp_density)) * ((np->volumeFraction[k] + p->volumeFraction[k])/2) * (np->vel - p->vel).dot(spiky(p, np));
                np = np->next;
              }
            }
          }
        }
      }
    }
  }

  for(uint i =0; i<phaseList.size(); i++)
  {
    for(uint j = 0; j<phaseList[i]->particleList.size(); j++)
    {
      p = phaseList[i]->particleList[j];
      p->colour = vec3(0.0f, 0.0f, 0.0f);
      float sum_vol_frac = 0.0f;
      std::vector<float> volFracChange;
      for(uint x = 0; x<phaseList.size(); x++)
      {
        volFracChange.push_back(0.0f);
      }
      pressureChange = 0.0f;
      for(uint k = 0; k<phaseList.size(); k++)
      {
        sum_vol_frac += p->volumeFraction[k];
      }
      //std::cout<<sum_vol_frac<<"\n";
      if(sum_vol_frac != 1.0f && sum_vol_frac != 0.0f)
      {
        float divisor = 1.0f/sum_vol_frac;
        for(uint l = 0; l<phaseList.size(); l++)
        {
          float oldVolFrac = p->volumeFraction[l];
          p->volumeFraction[l] *= divisor;
          //p->pressure *= divisor;
          volFracChange[l] = p->volumeFraction[l] - oldVolFrac;

//          std::cout<<p->volumeFraction[l]<<"\n";

        }
      }
      sum_vol_frac = 0.0f;
      for(uint k = 0; k<phaseList.size(); k++)
      {
        sum_vol_frac += p->volumeFraction[k];
      }
      if(sum_vol_frac <= 0.000001f || isnan(sum_vol_frac) == true)
      {
        p->volumeFraction[i] = 1.0f;
        for(uint l = 0; l< phaseList.size(); l++)
        {
          if( l != i)
          {
            p->volumeFraction[l] = 0.0f;
          }
        }
      }
      for(uint m = 0; m<phaseList.size(); m++)
      {
        p->pressure += (-gas_constant) * phaseList[m]->restDensity * volFracChange[m];
        //std::cout<<"phase "<<m<<" = "<<p->volumeFraction[m]<<"\n";
        p->colour.m_x += phaseList[m]->phaseColour.m_x * (p->volumeFraction[m]);
        p->colour.m_y += phaseList[m]->phaseColour.m_y * (p->volumeFraction[m]);
        p->colour.m_z += phaseList[m]->phaseColour.m_z * (p->volumeFraction[m]);
        if(p->colour.m_x > 1.0f)
        {
          p->colour.m_x = 1.0f;
        }
        if(p->colour.m_z > 1.0f)
        {
          p->colour.m_z = 1.0f;
        }
        p->colour.m_y = 0.0f;
      }
//      p->pressure = p->pressure + pressureChange;
    }
  }

}

int3 SPHSystem::calc_cell_pos(vec3 p)
{
	int3 cell_pos;
	cell_pos.m_x = int(floor((p.m_x) / cell_size));
	cell_pos.m_y = int(floor((p.m_y) / cell_size));
	cell_pos.m_z = int(floor((p.m_z) / cell_size));

	//std::cout<<cell_pos.m_x<<" "<< cell_pos.m_y<<" "<<cell_pos.m_z<<"\n";

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
    vec3 normal = normalize(r) * spiky_value * pow(kernel - rlen, 2);
//    std::cout<<normal.m_x<<" "<<normal.m_y<<" "<<normal.m_z<<"\n";
    return normal;
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
  float len = sqrt( X * X + Y * Y + Z * Z );
  if(isnan(len) == true)
  {
    return 0.0f;
  }
  else
  {
    return len;
  }
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

  else
  {
    normal = vec3(0.0f, 0.0f, 0.0f);
  }

  return normal;

}

void SPHSystem::driftVelocity()
{
  Particle *p;
  Particle *np;
  int3 cell_pos;
  int3 near_pos;
  uint hash;

  for(uint i=0; i<phaseList.size(); i++)
  {
    for(uint j=0; j<phaseList[i]->numParticle; j++)
    {
      p = phaseList[i]->particleList[j];

      p->sum_mass_density = 0.0f;
      p->sum_mass_diffusion = vec3(0.0f, 0.0f, 0.0f);
      p->sum_mass_pressure = vec3(0.0f, 0.0f, 0.0f);
      p->interp_density = 0.0f;
      cell_pos=calc_cell_pos(p->pos);
      float massFraction[]  {0.0f};

      for(uint k=0; k<phaseList.size(); k++)
      {
        p->pressureGradientk[k] = vec3(0.0f, 0.0f, 0.0f);
        p->diffusionGradient[k] = vec3(0.0f, 0.0f, 0.0f);
        massFraction[k] = (p->volumeFraction[k] * phaseList[k]->restDensity)/p->restDensity;
      }

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
              p->interp_density += p->restMass * poly6(p, np);

              np = np->next;
            }

            np =cell[hash];

            while(np != NULL)
            {
              vec3 spike = spiky(p, np);
              for(uint a=0; a<phaseList.size(); a++)
              {
                p->pressureGradientk[a] += spike * np->restMass/np->interp_density * (np->pressure - p->pressure);
                p->diffusionGradient[a] += spike * np->restMass/np->interp_density * (np->volumeFraction[i] - p->volumeFraction[i]);
              }

              np = np->next;
            }
          }
        }
      }
      p->pressure = 1.0f * (p->interp_density - p->restDensity);

      for(uint l=0; l<phaseList.size(); l++)
      {
        p->sum_mass_density += massFraction[l] * phaseList[l]->restDensity;
        p->sum_mass_pressure += p->pressureGradientk[l] * massFraction[l];
        if(p->volumeFraction[l] > 0.0f)
        {
          p->sum_mass_diffusion += (p->diffusionGradient[l] / p->volumeFraction[l]) * massFraction[l];
        }
      }

      float inertiaConst = 0.00000001;
      float diffusionConst = 0.0001;
      for(uint m=0; m<phaseList.size(); m++)
      {

        p->driftVelocity[m] = ((gravity - p->acc) * inertiaConst * (phaseList[m]->restDensity -  p->sum_mass_density) -
                              (p->pressureGradientk[m] - p->sum_mass_pressure) * inertiaConst);
        if(p->volumeFraction[m] > 0.0f)
        {
          p->driftVelocity[m] -= ((p->diffusionGradient[m] / p->volumeFraction[m]) - p->sum_mass_diffusion) * diffusionConst;
        }
      }

    }
  }

  for(uint i=0; i<phaseList.size(); i++)
  {
    for(uint j=0; j<phaseList[i]->numParticle; j++)
    {
      p = phaseList[i]->particleList[j];
      vec3 newVel(0.0f, 0.0f, 0.0f);

      for(uint k=0; k<phaseList.size(); k++)
      {
        newVel += (p->driftVelocity[k] + p->vel) * p->volumeFraction[k] * phaseList[k]->restDensity;
      }
      p->vel = newVel * (1.0f / p->restDensity);
    }
  }




//  vec3 pressureGradient;
//  int3 cell_pos;
//  int3 near_pos;
//  uint hash;

//  for(uint i=0; i<phaseList.size(); i++)
//  {
//    for(uint j=0; j<phaseList[i]->numParticle; j++)
//    {
//      p=phaseList[i]->particleList[j];
//      p->acc = vec3(0.0f, 0.0f, 0.0f);
//      p->pressureGradientk = vec3(0.0f, 0.0f, 0.0f);
//      p->diffusionGradient = vec3(0.0f, 0.0f, 0.0f);
//      p->pressureGradientm = vec3(0.0f, 0.0f, 0.0f);
//      cell_pos=calc_cell_pos(p->pos);
//      p->viscosity = 0.0f;
//      float sum_mass_dens = 0;
//      vec3 sum_drift_vel = vec3(0.0f, 0.0f, 0.0f);
//      for(uint l = 0; l < phaseList.size(); l++)
//      {
//        p->viscosity += p->volumeFraction[l] * phaseList[l]->phaseViscosity;
//      }

//      for(int x=-1; x<=1; x++)
//      {
//        for(int y=-1; y<=1; y++)
//        {
//          for(int z=-1; z<=1; z++)
//          {
//            near_pos.m_x=cell_pos.m_x+x;
//            near_pos.m_y=cell_pos.m_y+y;
//            near_pos.m_z=cell_pos.m_z+z;
//            hash=calc_cell_hash(near_pos);

//            if(hash == 0xffffffff)
//            {
//              continue;
//            }

//            np=p->next;

//            while(np != NULL)
//            {
//              vec3 spike = spiky(p, np);
//              p->pressureGradientk += spike * np->restMass/np->interp_density * (np->pressure - p->pressure);
//              p->diffusionGradient += spike * np->restMass/np->interp_density * (np->volumeFraction[i] - p->volumeFraction[i]);

//              p->pressureGradientm += spike * np->restMass * (np->pressure + p->pressure)/(2*np->interp_density);

//              sum_mass_dens += np->restMass/np->interp_density;

//              for(uint k = 0; k<phaseList.size(); k++)
//              {
//                sum_drift_vel += (np->driftVelocity[k] * p->volumeFraction[k] * (np->driftVelocity[k].dot(spike)) + p->driftVelocity[k] * p->volumeFraction[k] * (p->driftVelocity[k].dot(spike)) ) * p->restDensity;
//              }

//              p->div_viscosity_tensor += (np->vel - p->vel) * (np->restMass/p->interp_density) * (p->viscosity + np->viscosity) * ((vec3(np->pos - p->pos).dot(spike))/pow(length(vec3(np->pos - p->pos)), 2));

//              np=np->next;
//            }
//          }
//        }
//      }
//      p->conv_mom_change = sum_drift_vel * - sum_mass_dens;
//    }
//  }

//  for(uint i=0; i<phaseList.size(); i++)
//  {
//    for(uint j=0; j<phaseList[i]->numParticle; j++)
//    {
//      p=phaseList[i]->particleList[j];

//      for(uint k = 0; k<phaseList.size(); k++)
//      {
//        p->massFraction[k] = ((p->volumeFraction[k] * p->restDensity)/p->interp_density);
//        p->sum_mass_density += p->massFraction[k] * p->restDensity;
//        p->sum_mass_pressure += p->pressureGradientk * p->massFraction[k] ;

//        p->acc = (p->pressureGradientm/p->restDensity) - (p->conv_mom_change/p->restDensity) - (p->div_viscosity_tensor/p->restDensity);
//      }

//    }

//  }


//  for(uint i=0; i<phaseList.size(); i++)
//  {
//    for(uint j=0; j<phaseList[i]->numParticle; j++)
//    {
//      p=phaseList[i]->particleList[j];

//      for(uint k = 0; k<phaseList.size(); k++)
//      {

//        p->driftVelocity[k] = (p->acc * 1 * (phaseList[k]->restDensity - p->sum_mass_density) - (p->pressureGradientk - p->sum_mass_pressure) * 1);
//        std::cout<<p->driftVelocity[k].m_x<<" "<<p->driftVelocity[k].m_y<<" "<< p->driftVelocity[k].m_z<<"\n";

//      }

//    //pressureGradient += p->mass * ()
//    }

//  }
// // drift_velocity = 1 * (rest_density - mass * rest_density) * acceleration - 1 * (normalize(p->pres) );
}

void SPHSystem::calc_acceleration()
{
  Particle *p;
  Particle *np;
  int3 cell_pos;
  int3 near_pos;
  uint hash;

  for(uint i=0; i<phaseList.size(); i++)
  {
    for(uint j=0; j<phaseList[i]->numParticle; j++)
    {
      p = phaseList[i]->particleList[j];
      p->pressureGradientm = vec3(0.0f, 0.0f, 0.0f);
      p->div_viscosity_tensor = vec3(0.0f, 0.0f, 0.0f);
      p->conv_mom_change = vec3(0.0f, 0.0f, 0.0f);
      p->viscosity = 0.0f;
      cell_pos=calc_cell_pos(p->pos);

      vec3 sum_drift_vel = vec3(0.0f, 0.0f, 0.0f);

      for(uint k=0; k<phaseList.size(); k++)
      {
        p->viscosity += p->volumeFraction[k] * phaseList[k]->phaseViscosity;
      }

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

            np = cell[hash];
            while(np != NULL)
            {
              if(np == p)
              {
                np = np->next;
                continue;
              }
              vec3 dist = np->pos - p->pos;
              if(length(dist) > 0.00001f)
              {
                sum_drift_vel = vec3(0.0f, 0.0f, 0.0f);
                p->pressureGradientm += spiky(p, np) * np->restMass * ((p->pressure + np->pressure)/(2 * np->interp_density));
                p->div_viscosity_tensor += (np->vel - p->vel) * (np->restMass/np->interp_density) * (p->viscosity + np->viscosity) * ((np->pos - p->pos).dot(spiky(p, np))/((np->pos - p->pos).dot(np->pos - p->pos)));

                for(uint k = 0; k<phaseList.size(); k++)
                {
                  sum_drift_vel += (np->driftVelocity[k] * np->volumeFraction[k] * (np->driftVelocity[k].dot(spiky(p, np))) + p->driftVelocity[k] * p->volumeFraction[k] * (p->driftVelocity[k].dot(spiky(p, np))) ) * phaseList[k]->restDensity;
                }
                p->conv_mom_change -=  sum_drift_vel * (np->restMass/np->interp_density);
              }

              np = np->next;
            }
          }
        }
      }
      p->acc = (p->pressureGradientm*-1)/p->restDensity + gravity + p->div_viscosity_tensor/p->restDensity + (p->conv_mom_change*-1)/p->restDensity;
//      std::cout<< p->acc.m_x<<" "<< p->acc.m_y<<" "<< p->acc.m_z<<"\n";
    }
  }
}

