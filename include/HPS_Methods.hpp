#ifndef HPS_METHODS_HPP_
#define HPS_METHODS_HPP_

#include "HPS_Base.hpp"
#include "HPS_Patch.hpp"
#include "HPS_PatchSolver.hpp"
#include "HPS_Merge.hpp"

namespace hps {

// Globals
static int PARENT_ID = 0;

// Helper functions
int computeNumberOfPatches(int N_levels);
int computeGlobalID(int level, int level_ID);
int computeLevelID(int global_ID, int level);
int computeNextLevelID(int level_ID);
int computePrevLevelID(int level_ID);

// HPS Functions
Vector<Patch> setupHPS(CellGrid<double, 2> parent_grid, int N_levels);
void buildStage(Vector<Patch>& patches, int N_levels, PoissonProblem poisson, double lambda);
void solveStage(Vector<Patch>& patches, int N_levels, PoissonProblem poisson, Vector<double> g_parent, double lambda);



}

#endif // HPS_METHODS_HPP_
