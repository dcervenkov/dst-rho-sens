#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom1.h"
#include "TH1F.h"

#include <stdio.h>

static const int maxParticles = 300;
static const char fileName[] = "data/3events.root";