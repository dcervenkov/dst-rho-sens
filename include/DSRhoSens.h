#include "Particle.h"

static const int maxParticles = 300;
static std::vector< std::vector<Particle> >* events = new std::vector< std::vector<Particle> >;

int ReadEvents(char fileName[]);
bool GetRelevantParticles(int eventNo, Particle* DSD0, Particle* DSPi, Particle* RhoPi0, Particle* RhoPi);
bool GetDSRhoFromB0(Particle* B0, Particle* DS, Particle* Rho);
bool GetD0PiFromDS(Particle* DS, Particle* D0, Particle* Pi);
bool GetPi0PiFromRho(Particle* Rho, Particle* Pi0, Particle* Pi);