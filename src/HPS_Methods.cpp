#include "HPS_Methods.hpp"
using namespace std;
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
 * Computes the level ID for a patch from it's global ID and level in tree.
 * @param  global_ID Global ID in patch tree
 * @param  level     Level in patch tree
 * @return           ID within the patch tree level
 */
int computeLevelID(int global_ID, int level) {
    return global_ID - (pow(2, level) - 1);
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

std::pair<int, int> computeChildren(int ID, int level) {
    int level_ID = computeLevelID(ID, level);
    int first_child_level_ID = level_ID * 2;
    int second_child_level_ID = first_child_level_ID + 1;
    int first_child_ID = computeGlobalID(level+1, first_child_level_ID);
    int second_child_ID = computeGlobalID(level+1, second_child_level_ID);
    return std::pair<int, int>(first_child_ID, second_child_ID);
}

int computeParent(int ID, int level) {
    int level_ID = computeLevelID(ID, level);
    int parent_level_ID = level_ID / 2;
    return computeGlobalID(level-1, parent_level_ID);
}

Vector<Patch> setupHPS(CellGrid<double, 2> parent_grid, int N_levels) {

    int N_patches = computeNumberOfPatches(N_levels);
    Vector<Patch> patches(N_patches);
    Patch parent_patch(parent_grid);
    patches[0] = parent_patch;

    // Patch* patch_pointer;

    for (int level = 0; level < N_levels-1; level++) {
        for (int level_ID = 0; level_ID < pow(2, level); level_ID++) {
            int ID = computeGlobalID(level, level_ID);
            int first_ID = computeGlobalID(level+1, computeNextLevelID(level_ID));
            int second_ID = computeGlobalID(level+1, computeNextLevelID(level_ID)) + 1;

            // patch_pointer = &patches[ID];
            std::pair<Patch, Patch> patch_pair = patches[ID].split(first_ID, second_ID);
            patches[first_ID] = patch_pair.first;
            patches[second_ID] = patch_pair.second;
        }
    }

    return patches;

}

void buildStage(Vector<Patch>& patches, int N_levels, PoissonProblem poisson, double lambda) {

    int N_patches = computeNumberOfPatches(N_levels);
    for (int tau = N_patches-1; tau >= 0; tau--) {

        if (patches[tau].is_leaf) {

            // Build local DtN operator
            cout << "building local T for tau = " << tau << endl;
            patches[tau].buildT(lambda);

            // Compute local particular solution
            int N_leaf = patches[tau].grid.N_pts[X];
            patches[tau].u = Matrix<double>(N_leaf, N_leaf);
            patches[tau].f = Matrix<double>(N_leaf, N_leaf);
            for (int i = 0; i < N_leaf; i++) {
                for (int j = 0; j < N_leaf; j++) {
                    patches[tau].f.at(i,j) = poisson.f(patches[tau].grid.point(X, i), patches[tau].grid.point(Y, j));
                }
            }
            Vector<double> g_zero_leaf(4*patches[tau].grid.N_pts[X], 0.0);
            patches[tau].u = mapSolution(patches[tau].grid, g_zero_leaf, patches[tau].f, lambda);
            patches[tau].fhat = mapDirichletToNeumann(patches[tau].grid, g_zero_leaf, patches[tau].f, lambda);

        }
        else {

            // Get IDs for children of tau
            std::pair<int, int> children = computeChildren(tau, patches[tau].level);
            int alpha = children.first;
            int beta = children.second;
            cout << "merging alpha = " << alpha << " and beta = " << beta << " to form tau = " << tau << endl;

            // Check if vertical or horizontal merge and perform merge
            if ((patches[tau].level+1) % 2 == 0) {
                // Horizontal merge
                merge_values merged = mergeHorizontal2(
                    patches[alpha].T,
                    patches[beta].T,
                    patches[alpha].fhat,
                    patches[beta].fhat,
                    N_levels,
                    patches[tau].level+1
                );

                // Store data in patch
                patches[tau].T = merged.T;
                patches[tau].S = merged.S;
                patches[tau].fhat = merged.fhat;
            }
            else {
                // Vertical merge
                merge_values merged = mergeVertical2(
                    patches[alpha].T,
                    patches[beta].T,
                    patches[alpha].fhat,
                    patches[beta].fhat,
                    N_levels,
                    patches[tau].level+1
                );

                // Store data in patch
                patches[tau].T = merged.T;
                patches[tau].S = merged.S;
                patches[tau].fhat = merged.fhat;
            }
        }
    }
}


void solveStage(Vector<Patch>& patches, int N_levels, PoissonProblem poisson, Vector<double> g_parent, double lambda) {

    // Set sizes for later use
    int N_patches = computeNumberOfPatches(N_levels);
    int N_leaf = patches[PARENT_ID].grid.N_pts[X];
    int N_boundary = g_parent.size() / 4;

    // Downwards pass: Apply solution operator on patch boundaries to get interior solution
    patches[PARENT_ID].g = g_parent;
    for (int tau = 0; tau < N_patches; tau++) {
        if (patches[tau].is_leaf) {
            cout << "applying S to leaf patch tau = " << tau << endl;

            // Map Dirichlet data to solution on leaf
            Matrix<double> u_homogeneous = mapSolution(patches[tau].grid, patches[tau].g, patches[tau].f, lambda);

            // Add particular solution from build stage to homogeneous solution
            patches[tau].u = patches[tau].u + u_homogeneous;
        }
        else {
            cout << "applying S to parent patch tau = " << tau << endl;

            // Apply solution operator to exterior Dirichlet data to get interior Dirchlet data
            Vector<double> g_interior = patches[tau].S * patches[tau].g; // TODO: Need to add particular solution to g_interior

            // Create children Dirichlet data vectors
            //    Get children IDs
            std::pair<int, int> children = computeChildren(tau, patches[tau].level);
            int alpha = children.first;
            int beta = children.second;

            //    Horizontal vs. Vertical split
            if (patches[tau].level % 2 == 0) {
                // Square to two horizontal rectangles
                //    Get sizes for alpha and beta
                int N_half_boundary = N_boundary / 2;
                int N_alpha = 2*N_boundary + 2*N_half_boundary;
                int N_beta = N_alpha;

                //    Extract parent boundary data
                Vector<double> g0 = patches[tau].g.extract(0*N_boundary, N_boundary);
                Vector<double> g1 = patches[tau].g.extract(1*N_boundary, N_boundary);
                Vector<double> g2 = patches[tau].g.extract(2*N_boundary, N_boundary);
                Vector<double> g3 = patches[tau].g.extract(3*N_boundary, N_boundary);

                //     Create boundary data for alpha
                Vector<double> g_alpha(N_alpha);
                Vector<double> g_alpha0 = g0.extract(0, N_half_boundary);
                Vector<double> g_alpha1 = g1.extract(0, N_half_boundary);
                Vector<double> g_alpha2 = g2;
                Vector<double> g_alpha3 = g_interior;
                g_alpha.intract(0*N_half_boundary, g_alpha0);
                g_alpha.intract(1*N_half_boundary, g_alpha1);
                g_alpha.intract(2*N_half_boundary, g_alpha2);
                g_alpha.intract(4*N_half_boundary, g_alpha3);
                patches[alpha].g = g_alpha;

                //    Create boundary data for beta
                Vector<double> g_beta(N_beta);
                Vector<double> g_beta0 = g0.extract(N_half_boundary, N_half_boundary);
                Vector<double> g_beta1 = g1.extract(N_half_boundary, N_half_boundary);
                Vector<double> g_beta2 = g_interior;
                Vector<double> g_beta3 = g3;
                g_beta.intract(0*N_half_boundary, g_beta0);
                g_beta.intract(1*N_half_boundary, g_beta1);
                g_beta.intract(2*N_half_boundary, g_beta2);
                g_beta.intract(4*N_half_boundary, g_beta3);
                patches[beta].g = g_beta;
            }
            else {
                // Rectangle to two squares
                //    Get sizes for alpha and beta
                int N_half_boundary = N_boundary / 2;
                int N_alpha = 4*N_half_boundary;
                int N_beta = N_alpha;

                //    Extract parent boundary data
                Vector<double> g0 = patches[tau].g.extract(0*N_half_boundary, N_half_boundary);
                Vector<double> g1 = patches[tau].g.extract(1*N_half_boundary, N_half_boundary);
                Vector<double> g2 = patches[tau].g.extract(2*N_half_boundary, N_boundary);
                Vector<double> g3 = patches[tau].g.extract(4*N_half_boundary, N_boundary);

                //    Create boundary data for alpha
                Vector<double> g_alpha(N_alpha);
                Vector<double> g_alpha0 = g0;
                Vector<double> g_alpha1 = g_interior;
                Vector<double> g_alpha2 = g2.extract(0, N_half_boundary);
                Vector<double> g_alpha3 = g3.extract(0, N_half_boundary);
                g_alpha.intract(0*N_half_boundary, g_alpha0);
                g_alpha.intract(1*N_half_boundary, g_alpha1);
                g_alpha.intract(2*N_half_boundary, g_alpha2);
                g_alpha.intract(3*N_half_boundary, g_alpha3);
                patches[alpha].g = g_alpha;

                //    Create boundary data for beta
                Vector<double> g_beta(N_beta);
                Vector<double> g_beta0 = g_interior;
                Vector<double> g_beta1 = g1;
                Vector<double> g_beta2 = g2.extract(N_half_boundary, N_half_boundary);
                Vector<double> g_beta3 = g3.extract(N_half_boundary, N_half_boundary);
                g_beta.intract(0*N_half_boundary, g_beta0);
                g_beta.intract(1*N_half_boundary, g_beta1);
                g_beta.intract(2*N_half_boundary, g_beta2);
                g_beta.intract(3*N_half_boundary, g_beta3);
                patches[beta].g = g_beta;

                // Update boundary size
                if (patches[tau].level % 2 == 0) {
                    N_boundary = N_boundary / 2;
                } // TODO: Fix for greater than 4-to-1 patch!
            }

        }


    }

}

} // NAMESPACE hps
