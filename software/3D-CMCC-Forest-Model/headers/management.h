/*
 * management.h
 *
 *  Created on: 18/ago/2016
 *      Author: alessio-cmcc
 */

#ifndef MANAGEMENT_H_
#define MANAGEMENT_H_

#include "matrix.h"

typedef struct {
	int *thinning_years;
	int thinning_years_count;
	int *harvesting_years;
	int harvesting_years_count;
        int *thinning_intensity;     //added by ddalmo may 2021
        int thinning_intensity_count; 
        int *regeneration_years;     //added by ddalmo august 2021
        int regeneration_years_count; 
} management_t;

int forest_management (cell_t *const c, const int day, const int month, const int year);

void choose_management(species_t *const s, const int years);

void clearcut_timber_upon_request(cell_t *const c, const int layer, const int height, const int dbh, const int age, const int species);

void thinning (cell_t *const c, const int height, const int dbh, const int age, const int species, const int year, const int THINNING_INDEX, const int counter_thinning);

void prescribed_thinning (cell_t *const c, const int height, const int dbh, const int age, const int species, const int year);

//void harvesting (cell_t *const c, const int height, const int dbh, const int age, const int species);

void clearcut_coppice(cell_t *const c, const int layer, const int height, const int dbh, const int age, const int species);

//management_t* management_load(const char* const filename);
management_t* management_load_presc(const char* const filename);  //5p6

void management_free(management_t* p);

#endif /* MANAGEMENT_H_ */
