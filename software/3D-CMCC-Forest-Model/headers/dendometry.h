/* dendometry.h */
#ifndef DENDOMETRY_H_
#define DENDOMETRY_H_

#include "matrix.h"

void dendrometry_old(cell_t *const c, const int layer, const int height, const int dbh, const int age, const int species, const int year );

void dendrometry(cell_t *const c, const int layer, const int height, const int dbh, const int age, const int species, const meteo_daily_t *const meteo_daily, const int month, const int year );

void annual_tree_increment ( cell_t *const c, const int height, const int dbh, const int age, const int species, const int year );

void annual_gross_tree_increment ( cell_t *const c, const int height, const int dbh, const int age, const int species, const int year );


void annual_minimum_reserve (species_t *const s);

#endif /* DENDOMETRY_H_ */
