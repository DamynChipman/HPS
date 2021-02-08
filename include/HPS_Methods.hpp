#ifndef HPS_METHODS_HPP_
#define HPS_METHODS_HPP_

#include "HPS_Base.hpp"
#include "HPS_Patch.hpp"

namespace hps {

// Helper functions
int computeNumberOfPatches(int N_levels);
int computeGlobalID(int level, int level_ID);
int computeNextLevelID(int level_ID);
int computePrevLevelID(int level_ID);

// HPS Functions
Vector<Patch> setupHPS(CellGrid<double, 2> parent_grid, int N_levels);
void buildStage(Vector<Patch> patches, int N_levels, double lambda);
void solveStage();



}

#endif // HPS_METHODS_HPP_
