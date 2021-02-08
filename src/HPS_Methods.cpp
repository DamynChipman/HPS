#include "HPS_Methods.hpp"

namespace hps {

/**
 * Computes the number of patches for a given level
 * @param  levels Number of levels
 * @return        Number of patches
 */
int computeNumberOfPatches(int levels) {
    int N_patches = 0;
    for (int n = 0; n < levels; n++) {
        N_patches += pow(2, n);
    }
    return N_patches;
}

/**
 * Computes the global ID for a patch in the patch tree, given the `level` and the `level_ID`.
 * @param  level    Level in patch tree
 * @param  level_ID ID within the patch tree level
 * @return          Global ID in patch tree
 */
int computeGlobalID(int level, int level_ID) {
    return (pow(2, level) - 1) + level_ID;
}

/**
 * Computes the next level's `level_ID` from the current level
 * @param  level_ID ID within the patch tree level
 * @return          Next level's `level_ID`
 */
int computeNextLevelID(int level_ID) {
    return 2*level_ID;
}

/**
 * Computes the previous level's `level_ID` from the current level
 * @param level_ID ID within the patch tree level
 * @return         Previous level's `level_ID`
 */
int computePrevLevelID(int level_ID) {
    return level_ID / 2;
}

Vector<Patch> setupHPS(CellGrid<double, 2> parent_grid, int N_levels) {

    int N_patches = computeNumberOfPatches(N_levels);
    Vector<Patch> patches(N_patches);
    Patch parent_patch(parent_grid, 0, 0, true);
    patches[0] = parent_patch;

    Patch* patch_pointer;

    for (int level = 0; level < N_levels-1; level++) {
        for (int level_ID = 0; level_ID < pow(2, level); level_ID++) {
            int ID = computeGlobalID(level, level_ID);
            int first_ID = computeGlobalID(level+1, computeNextLevelID(level_ID));
            int second_ID = computeGlobalID(level+1, computeNextLevelID(level_ID)) + 1;

            patch_pointer = &patches[ID];
            std::pair<Patch, Patch> patch_pair = patch_pointer->split(first_ID, second_ID);
            patches[first_ID] = patch_pair.first;
            patches[second_ID] = patch_pair.second;
        }
    }

    return patches;

}

void buildStage(Vector<Patch> patches, int N_levels, double lambda) {

    int N_patches = computeNumberOfPatches(N_levels);
    for (int tau = N_patches-1; tau >= 0; tau--) {

        if (patches[tau].is_leaf) {

            // Build local DtN operator
            patches[tau].T = patches[tau].buildDirichletToNeumannMatrix(lambda);

        }
        // else {
        //
        // }

    }

}


// void solveStage();

} // NAMESPACE hps
