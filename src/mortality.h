/* mortality.h */
#ifndef MORTALITY_H_
#define MORTALITY_H_

#include "matrix.h"

int growth_efficiency_mortality        ( cell_t *const c, const int height, const int dbh, const int age, const int species );

int annual_growth_efficiency_mortality ( cell_t *const c, const int height, const int dbh, const int age, const int species );

void self_pruning                      ( cell_t *const c, const int height, const int dbh, const int age, const int species, const double old, const double current );

void self_thinning_mortality           ( cell_t *const c, const int layer, const int year );

void age_mortality                     ( cell_t *const c, const int height, const int dbh, const int age, const int species );

void stochastic_mortality              ( cell_t *const c, const int height, const int dbh, const int age, const int species );

void mortality                         (cell_t *const c, const int height, const int dbh, const int age, const int species, const int tree_remove);

void pruning_daily					   (matrix_t* m, int cell_index, pruning_t* p);

void saplings_mortality                (cell_t *const c, const int height, const int dbh, const int age, const int species);

#endif /* MORTALITY_H_ */
