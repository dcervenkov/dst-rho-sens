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
        Particle* GetMother() const;
        void SetMother(Particle& val);
        Particle* GetDaughter(int i) const;
        void SetDaughter(int i, Particle& val);
        int GetNumDaughters() const;
        void RotateP(TRotation rot);
        void BoostP(TVector3 beta);
        TVector3 GetBoost();
        double GetPhi();
        double GetTheta();

    private:
        int m_idhep;
        /// p = (p_x,p_y,p_z,E)
        TLorentzVector m_p;
        /// v = (x,y,z,t)
        TLorentzVector m_v;
        Particle* m_mother;
        Particle* m_daughter[maxDaughters];
        int m_numDaughters;
};

#endif