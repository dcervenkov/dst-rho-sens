#ifndef PRINTTREE_H_INCLUDED
#define PRINTTREE_H_INCLUDED

#include <string>
#include <vector>

#include "Particle.h"

namespace printTree {

std::string getRealName(int idhep);
void printDecayTree(std::vector<Particle> list);
void printParticle(const Particle&, std::string leader);
bool isParticleLastDaughter(const Particle& particle);

}  // namespace printTree

#endif