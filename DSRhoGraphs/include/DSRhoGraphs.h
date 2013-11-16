#include "TH2D.h"

static const double PI = 3.141592;

std::vector<TH2D*> CreateMaps(const RooDataSet* const dataset);
void CreateGeneralPlots(RooDataSet* const dataset, const ObservablesCollection c);
