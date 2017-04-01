#ifndef CONSTANTS_H_INCLUDED
#define CONSTANTS_H_INCLUDED

static const double PI = 3.141592;

static const int B0_IDHEP = 511;
static const int DS_IDHEP = 413;
static const int RHO_IDHEP = 213;
static const int D0_IDHEP = 421;
static const int PI0_IDHEP = 111;
static const int PI_IDHEP = 211;
static const int PHOTON_IDHEP = 22;

static const Double_t Btau = 1.53439;  // [ps]
static const Double_t Bdm = 0.507;     // [h-bar (ps)^-1]

static const Double_t c = 0.29979;  // [mm/ps]

static const Double_t Bgamma = 1 / Btau / c;  // 2.17242; // [h-bar (mm/c)^-1]
static const Double_t Bfreq = Bdm / c;  // 1.69117; //[h-bar/c^2 (mm/c)^-1] number taken directly from evtgen (in debug mode)
#endif
