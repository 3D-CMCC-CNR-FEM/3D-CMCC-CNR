/* biomass.h */
#ifndef BIOMASS_H_
#define BIOMASS_H_

#include "matrix.h"

void live_total_wood_age   ( const age_t *const a, species_t *const s );

void average_tree_pools    ( cell_t *const c );

void tree_branch_and_bark  ( cell_t *const c, const int height, const int dbh, const int age, const int species );

void abg_bgb_biomass       ( cell_t *const c, const int height, const int dbh, const int age, const int species );

void tree_biomass_remove   ( cell_t *const c, const int height, const int dbh, const int age, const int species, const int tree_remove, const int nat_man );

#endif /* BIOMASS_H */
