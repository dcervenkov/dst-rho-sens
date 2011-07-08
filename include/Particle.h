#ifndef PARTICLE_H_INCLUDED
#define PARTICLE_H_INCLUDED

static const int maxDaughters = 10;

class Particle
{
    public:
        Particle();
        ~Particle();
        int GetIdhep();
        void SetIdhep(int val);
        double GetP(int i);
        void SetP(int i, double val);
        double GetV(int i);
        void SetV(int i, double val);
        double GetM();
        void SetM(double val);
        Particle* GetMother();
        void SetMother(Particle& val);
        Particle* GetDaughter(int i);
        void SetDaughter(int i, Particle& val);
        int GetNumDaughters();

    private:
        int m_idhep;
        double m_p[4];
        double m_v[4];
        double m_m;
        Particle* m_mother;
        Particle* m_daughter[maxDaughters];
        int m_numDaughters;
};

#endif