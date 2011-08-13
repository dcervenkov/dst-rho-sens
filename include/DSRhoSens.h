#include "Particle.h"

static const int maxParticles = 300;
static std::vector< std::vector<Particle> >* events = new std::vector< std::vector<Particle> >;

int ReadEvents(char fileName[]);
void PrintEvent(int evtNo);
void PrintRelevantParticles(const Particle* DS,const Particle* DSD0,const Particle* DSPi, \
                          const Particle* Rho,const Particle* RhoPi0,const Particle* RhoPi);
bool GetRelevantParticles(int eventNo, Particle** B, Particle** DS, Particle** DSD0, Particle** DSPi, \
                          Particle** Rho, Particle** RhoPi0, Particle** RhoPi);
bool GetDSRhoFromB0(const Particle* const B0, Particle** DS, Particle** Rho);
bool GetD0PiFromDS(const Particle* const DS, Particle** D0, Particle** Pi);
bool GetPi0PiFromRho(const Particle* const Rho, Particle** Pi0, Particle** Pi);
TRotation GetRotationToZ(const Particle* const DS);
void TransformHel (Particle* B0, Particle* DS, Particle* DSD0, Particle* DSPi, \
                           Particle* Rho, Particle* RhoPi0, Particle* RhoPi);
void GetAngles(Particle* a, Particle* b, double& chi, double& theta_a, double& theta_b);