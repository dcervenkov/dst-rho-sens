#include "Particle.h"

static const int maxParticles = 300;
static std::vector< std::vector<Particle> >* events = new std::vector< std::vector<Particle> >;

int ReadEvents(char fileName[]);
void PrintEvent(int evtNo);
void PrintEventOrig(int evtNo,int numParticles,int* id,int* idhep,int* mother,int* da1,int* da2,double(*p)[5],double(*v)[4]);
void PrintRelevantParticles(const Particle* DS,const Particle* DSD0,const Particle* DSPi, \
                            const Particle* Rho,const Particle* RhoPi0,const Particle* RhoPi);
bool GetRelevantParticles(int eventNo, Particle** B, Particle** DS, Particle** DSD0, Particle** DSPi, \
                          Particle** Rho, Particle** RhoPi0, Particle** RhoPi);
bool GetDSRhoFromB0(const Particle* const B0, Particle** DS, Particle** Rho);
bool GetD0PiFromDS(const Particle* const DS, Particle** D0, Particle** Pi);
bool GetPi0PiFromRho(const Particle* const Rho, Particle** Pi0, Particle** Pi);
TRotation GetRotationToZ(const Particle* const DS);
TRotation GetZRotationToX(const Particle* const RhoPi);
void TransformHel(Particle* B0, Particle* DS, Particle* DSD0, Particle* DSPi, \
                  Particle* Rho, Particle* RhoPi0, Particle* RhoPi);
void TransformTrans(Particle* B0, Particle* DS, Particle* DSD0, Particle* DSPi, \
                    Particle* Rho, Particle* RhoPi0, Particle* RhoPi);
void GetAnglesHel(Particle* a, Particle* b, double& chi, double& theta_a, double& theta_b);
void GetAnglesTrans(Particle* a, Particle* b, double& theta_t, double& phi_t, double& theta_b);
void Analyze();
