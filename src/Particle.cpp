#include "Particle.h"
#include "ASSERT.h"

Particle::Particle()
{
    m_idhep = 0;
    for(int i = 0; i < 4; i++)
    {
        m_p[i] = 0;
        m_v[i] = 0;
    }
    m_m = -1;
    m_mother = 0;
    for(int i = 0; i < maxDaughters; i++)
        m_daughter[i] = 0;
    m_numDaughters = 0;
}

Particle::~Particle()
{
    //dtor
}

int Particle::GetIdhep()
{
    return m_idhep;
}

void Particle::SetIdhep(int val)
{
    m_idhep = val;
}

double Particle::GetP(int i)
{
    ASSERT(i >= 0 && i < 4);
    return m_p[i];
}

void Particle::SetP(int i, double val)
{
    ASSERT(i >= 0 && i < 4);
    m_p[i] = val;
}

double Particle::GetV(int i)
{
    ASSERT(i >= 0 && i < 4);
    return m_v[i];
}

void Particle::SetV(int i, double val)
{
    ASSERT(i >= 0 && i < 4);
    m_v[i] = val;
}

double Particle::GetM()
{
    return m_m;
}

void Particle::SetM(double val)
{
    m_m = val;
}

Particle* Particle::GetMother()
{
    return m_mother;
}

void Particle::SetMother(Particle& val)
{
    m_mother = &val;
}

Particle* Particle::GetDaughter(int i)
{
    ASSERT(i < maxDaughters);
    return m_daughter[i];
}

void Particle::SetDaughter(int i, Particle& val)
{
    ASSERT(i < maxDaughters);
    m_daughter[i] = &val;
    m_numDaughters++;
}

int Particle::GetNumDaughters()
{
    return m_numDaughters;
}