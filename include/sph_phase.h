#ifndef __SPH_PHASE_H__
#define __SPH_PHASE_H__

#include <vector>
#include <iostream>
#include "sph_system.h"
#include "sph_particle.h"


class Phase
{
public:

  Phase(){}

  ~Phase(){}

  uint max_particle;
  float restDensity;
  float restMass;
  vec3 phaseColour;
  float phaseViscosity;
  float numParticle = 0;


  std::vector<Particle *> particleList;

  void addParticle(Particle *p) { particleList.push_back(p); }

};


#endif
