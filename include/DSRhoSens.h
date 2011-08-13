#include "Particle.h"

static const int maxParticles = 300;
static std::vector< std::vector<Particle> >* events = new std::vector< std::vector<Particle> >;

int ReadEvents(char fileName[]);
void PrintEvent(int evtNo);
bool GetRelevantParticles(int eventNo, Particle** DS, Particle** DSD0, Particle** DSPi, \
                          Particle** Rho, Particle** RhoPi0, Particle** RhoPi);
bool GetDSRhoFromB0(const Particle* const B0, Particle** DS, Particle** Rho);
bool GetD0PiFromDS(const Particle* const DS, Particle** D0, Particle** Pi);
bool GetPi0PiFromRho(const Particle* const Rho, Particle** Pi0, Particle** Pi);
TRotation GetRotationToZ(const Particle* const DS);