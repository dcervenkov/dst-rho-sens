#ifndef PARTICLE_H_INCLUDED
#define PARTICLE_H_INCLUDED

static const int maxDaughters = 10;

class Particle
{
    public:
        Particle();
        ~Particle();
        int GetIdhep() const;
        void SetIdhep(int val);
        double GetP(int i) const;
        void SetP(int i, double val);
        double GetV(int i) const;
        void SetV(int i, double val);
        double GetM() const;
        void SetM(double val);
        Particle* GetMother() const;
        void SetMother(Particle& val);
        Particle* GetDaughter(int i) const;
        void SetDaughter(int i, Particle& val);
        int GetNumDaughters() const;

    private:
        int m_idhep;
        /// p = (p_x,p_y,p_z,E)
        double m_p[4];
        /// v = (x,y,z,t)
        double m_v[4];
        double m_m;
        Particle* m_mother;
        Particle* m_daughter[maxDaughters];
        int m_numDaughters;
};

#endif