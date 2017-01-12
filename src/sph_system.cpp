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

  world_size.m_x=0.32f;
  world_size.m_y=0.32f;
  world_size.m_z=0.32f;
  cell_size=0.02f;
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
  time_step=0.003f;
  surf_norm=6.0f;
  surf_coe=0.1f;

  //Precomputed kernel values
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

  //Creating and assigning values to phases
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
  phase1->restMass = 0.02f;
  phase2->restMass = 0.01f;
  phase1->phaseViscosity = 15.0f;
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

  //Seperate out into function steps

  //Build up hash table of particle positions
  build_table();

  // Paper step 1
  driftVelocity();

  // Paper step 2 & 3 (enabling currently breaks the simulation, but shows colour mixing)
//  adv_vol_frac();

  // Paper step 4
  calc_acceleration();

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

  //Adds particles in a grid in the world
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

  std::cout<<"Number of phases:"<<phaseList.size() << "\n";
  std::cout<<"Init Particle: "<<num_particle<<"\n";
  }

void SPHSystem::add_particle(Phase *phase, vec3 pos, vec3 vel)
{

  Particle *p=&(mem[num_particle]);

  //Set up starting values for particle
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

  //Ensure correct amount of values in vector arrays
  for(uint i = 0; i<phaseList.size(); i++)
  {
    p->driftVelocity.push_back(vec3(0.0f, 0.0f, 0.0f));
    p->massFraction.push_back(0.0f);
    p->diffusionGradient.push_back(vec3(0.0f, 0.0f, 0.0f));
    p->pressureGradientk.push_back(vec3(0.0f, 0.0f, 0.0f));
  }

  p->next=NULL;

  //Initial volume fractions
  for(uint i = 0; i<num_phases; i++)
  {
    //The phases the particle is not part of
    if(phaseList[i] != phase)
    {
      p->volumeFraction.push_back(0.0f);
    }

    //The particle's phase
    else
    {
      p->volumeFraction.push_back(1.0f);
    }
  }
  //Push the particle into this phase
  phase->particleList.push_back(p);

  phase->numParticle++;
  num_particle++;
}

void SPHSystem::build_table()
{
  Particle *p;
  uint hash;

  for(uint i=0; i<tot_cell; i++)
  {
    cell[i]=NULL;
  }

//Loop through the particles in each phase
  for(uint i=0; i<phaseList.size(); i++)
  {
    for(uint j=0; j<phaseList[i]->numParticle; j++)
    {
      p=phaseList[i]->particleList[j];

      hash=calc_cell_hash(calc_cell_pos(p->pos));

      //In some cases the hash can go below zero, hence this number
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

void SPHSystem::advection()
{
  Particle *p;
  for(uint i=0; i<phaseList.size(); i++)
  {
    for(uint j=0; j<phaseList[i]->numParticle; j++)
    {
      p=phaseList[i]->particleList[j];

      //Incrementing position based on the particle's acceleration, velocity and the time step
      vec3 oldpos = p->pos;
      p->pos=p->pos+p->vel*time_step + p->acc*time_step*time_step*0.5;
      p->vel= (p->pos - oldpos) / time_step;

      //Ensuring the particle does not go out of bounds
      if(p->pos.m_x >= world_size.m_x-BOUNDARY)
      {
              p->vel.m_x=p->vel.m_x*wall_damping;
              p->pos.m_x=world_size.m_x-BOUNDARY - 0.001f;
      }

      if(p->pos.m_x < 0.0f)
      {
              p->vel.m_x=p->vel.m_x*wall_damping;
              p->pos.m_x=0.001f;
      }

      if(p->pos.m_y >= world_size.m_y-BOUNDARY)
      {
              p->vel.m_y=p->vel.m_y*wall_damping;
              p->pos.m_y=world_size.m_y-BOUNDARY - 0.001f;
      }

      if(p->pos.m_y < 0.0f)
      {
              p->vel.m_y=p->vel.m_y*wall_damping;
              p->pos.m_y=0.001f;
      }

      if(p->pos.m_z >= world_size.m_z-BOUNDARY)
      {
              p->vel.m_z=p->vel.m_z*wall_damping;
              p->pos.m_z=world_size.m_z-BOUNDARY - 0.001f;
      }

      if(p->pos.m_z < 0.0f)
      {
              p->vel.m_z=p->vel.m_z*wall_damping;
              p->pos.m_z=0.001f;
      }

      p->ev=(p->ev+p->vel)/2;
    }
  }
}

//Volume fraction advection
void SPHSystem::adv_vol_frac()
{
  Particle *p;
  Particle *np;
  int3 cell_pos;
  int3 near_pos;
  uint hash;

  for(uint i =0; i<phaseList.size(); i++)
  {
    for(uint j = 0; j<phaseList[i]->particleList.size(); j++)
    {
      p = phaseList[i]->particleList[j];
      cell_pos=calc_cell_pos(p->pos);

      //Loop through each phase
      for(uint k = 0; k < phaseList.size(); k++)
      {
        //Loop through all neighbours
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
                //Prevent dividing by zero
                if(np->interp_density == 0.0f)
                {
                  np = np->next;
                  continue;
                }
                // Particle's volume fraction for this phase
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
      //Reset colour to prevent using old values
      p->colour = vec3(0.0f, 0.0f, 0.0f);
      //Used to check if a particle's volume fractions add up to 1
      float sum_vol_frac = 0.0f;
      //Change in each volume fraction
      std::vector<float> volFracChange;
      for(uint x = 0; x<phaseList.size(); x++)
      {
        volFracChange.push_back(0.0f);
      }

      for(uint k = 0; k<phaseList.size(); k++)
      {
        sum_vol_frac += p->volumeFraction[k];
      }
      //std::cout<<sum_vol_frac<<"\n";
      if(sum_vol_frac != 1.0f && sum_vol_frac != 0.0f)
      {
        //Readjust so volume fractions do add to 1
        float divisor = 1.0f/sum_vol_frac;
        for(uint l = 0; l<phaseList.size(); l++)
        {
          float oldVolFrac = p->volumeFraction[l];
          p->volumeFraction[l] *= divisor;
          volFracChange[l] = p->volumeFraction[l] - oldVolFrac;
        }
      }
      //Recalculate volume fraction total
      sum_vol_frac = 0.0f;
      for(uint k = 0; k<phaseList.size(); k++)
      {
        sum_vol_frac += p->volumeFraction[k];
      }
      //Check if the sum is 0 (accounting for float precision) or if it is not-a-number
      if(sum_vol_frac <= 0.000001f || isnan(sum_vol_frac) == true)
      {
        //Reset volume fractions to original values
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
        //Adjust pressure as stated in the paper (equation 30)
        p->pressure += (((-gas_constant)*phaseList[m]->restDensity)/7) * ((7-1)*pow((p->interp_density/p->restDensity), 7) + 1) * volFracChange[m];
        //Colour particles based on their volume fractions
        //This creates the "mixing" visual effect
        p->colour.m_x += phaseList[m]->phaseColour.m_x * (p->volumeFraction[m]);
        p->colour.m_y += phaseList[m]->phaseColour.m_y * (p->volumeFraction[m]);
        p->colour.m_z += phaseList[m]->phaseColour.m_z * (p->volumeFraction[m]);
        //Ensure colour values do not exceed 1
        if(p->colour.m_x > 1.0f)
        {
          p->colour.m_x = 1.0f;
        }
        if(p->colour.m_y > 1.0f)
        {
          p->colour.m_y = 1.0f;
        }
        if(p->colour.m_z > 1.0f)
        {
          p->colour.m_z = 1.0f;
        }
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
  //Poly6 kernel calculation (for some Wij instances in the paper)
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
  //Spiky kernel calculation
  vec3 r = p->pos - np->pos;
  float rlen= length(r);

  if(rlen>kernel)
  {
    return vec3(0.0f, 0.0f, 0.0f);
  }
  else
  {
    vec3 normal = normalize(r) * spiky_value * pow(kernel - rlen, 2);
    return normal;
  }

}

float SPHSystem::visco(Particle *p, Particle *np)
{
  //Unused viscosity kernel calculation, may be useful in the future
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
  //Function to calculate a vector's length
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
  //Function to normalise a vector
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

      //Resetting variables to prevent use of old values
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
        //In some instances parts of equations need to be pre-calculated. I have given names to
        //parts which the paper does not clearly name. The massFraction is labelled Ck in the paper
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
              //Pre-calculate interpolated density
              p->interp_density += p->restMass * poly6(p, np);

              np = np->next;
            }

            np =cell[hash];

            while(np != NULL)
            {
              vec3 spike = spiky(p, np);
              for(uint a=0; a<phaseList.size(); a++)
              {
                //More pre-calculations
                if(np->interp_density != 0.0f)
                {
                  p->pressureGradientk[a] += spike * np->restMass/np->interp_density * (np->pressure - p->pressure);
                  p->diffusionGradient[a] += spike * np->restMass/np->interp_density * (np->volumeFraction[i] - p->volumeFraction[i]);
                }
              }

              np = np->next;
            }
          }
        }
      }
      p->pressure = gas_constant * (p->interp_density - p->restDensity);

      for(uint l=0; l<phaseList.size(); l++)
      {
        //The summations from the three "sections" of the drift velocity equation
        p->sum_mass_density += massFraction[l] * phaseList[l]->restDensity;
        p->sum_mass_pressure += p->pressureGradientk[l] * massFraction[l];
        if(p->volumeFraction[l] > 0.0f)
        {
          p->sum_mass_diffusion += (p->diffusionGradient[l] / p->volumeFraction[l]) * massFraction[l];
        }
      }

      //User-defined constants
      float inertiaConst = 0.00000001;
      float diffusionConst = 0.0001;
      for(uint m=0; m<phaseList.size(); m++)
      {
        //Putting together the pieces
        p->driftVelocity[m] = ((gravity - p->acc) * inertiaConst * (phaseList[m]->restDensity -  p->sum_mass_density) -
                              (p->pressureGradientk[m] - p->sum_mass_pressure) * inertiaConst);
        //Prevent division by zero
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
      //Calculating the particle's new velocity
      vec3 newVel(0.0f, 0.0f, 0.0f);

      for(uint k=0; k<phaseList.size(); k++)
      {
        newVel += (p->driftVelocity[k] + p->vel) * p->volumeFraction[k] * phaseList[k]->restDensity;
      }
      p->vel = newVel * (1.0f / p->restDensity);
    }
  }
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
      //Variable resets
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
              //Ensure the particles are not inside each other, preventing division by zero
              vec3 dist = np->pos - p->pos;
              if(length(dist) > 0.00001f)
              {
                sum_drift_vel = vec3(0.0f, 0.0f, 0.0f);
                //Pressure gradient (equation 20)
                p->pressureGradientm += spiky(p, np) * np->restMass * ((p->pressure + np->pressure)/(2 * np->interp_density));
                //Divergence of viscosity tensor
                p->div_viscosity_tensor += (np->vel - p->vel) * (np->restMass/np->interp_density) * (p->viscosity + np->viscosity) * ((dist).dot(spiky(p, np))/((dist).dot(dist)));

                for(uint k = 0; k<phaseList.size(); k++)
                {
                  //The summation utilising drift velocities in euqation 19, pre-calculated
                  sum_drift_vel += (np->driftVelocity[k] * np->volumeFraction[k] * (np->driftVelocity[k].dot(spiky(p, np))) + p->driftVelocity[k] * p->volumeFraction[k] * (p->driftVelocity[k].dot(spiky(p, np))) ) * phaseList[k]->restDensity;
                }
                //Convective momentum change (equation 19)
                p->conv_mom_change -=  sum_drift_vel * (np->restMass/np->interp_density);
              }

              np = np->next;
            }
          }
        }
      }
      //Acceleration calculation (equation 8)
      p->acc = (p->pressureGradientm*-1)/p->restDensity + gravity + p->div_viscosity_tensor/p->restDensity + (p->conv_mom_change*-1)/p->restDensity;
    }
  }
}

