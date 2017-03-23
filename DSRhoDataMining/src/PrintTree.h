#ifndef PRINTTREE_H_INCLUDED
#define PRINTTREE_H_INCLUDED

#include <vector>

#include "Particle.h"

namespace printTree {

std::string getRealName(int idhep);
void printEvent(std::vector<Particle> list);
void printParticle(const Particle&, int level);

}  // namespace printTree

#endif