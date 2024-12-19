/*
 * C-fruit-partitioning-allocatio.h
 *
 *  Created on: 24 giu 2016
 *      Author: alessio
 */

#ifndef SRC_C_FRUIT_PARTITIONING_ALLOCATION_H_
#define SRC_C_FRUIT_PARTITIONING_ALLOCATION_H_

#include "matrix.h"

void Daily_C_Evergreen_Partitioning_Allocation (species_t *const s, cell_t *const c, const meteo_t *const met, const int day, const int month,
		const int years, const int DaysInMonth, const int height, const int age, const int species);

#endif /* SRC_C_FRUIT_PARTITIONING_ALLOCATION_H_ */

/*
 * C-fruit-partitioning.h
 *
 * Created on: 21 sep 2022
 * Author: Vincenzo Saponaro
 */

#ifndef C_FRUIT_PARTITIONING_H_
#define C_FRUIT_PARTITIONING_H_

#include "matrix.h"

//int Fruit_to_seeds_function (cell_t *const c, const int layer, const age_t *const a, species_t *const s, const int year);
void Fruit_to_seeds_function_npp (cell_t *const c, const age_t *const a, species_t *const s, const int day, const int month, const int year);
//void Fruit_to_seeds_function_structure (cell_t *const c, const int height, const int dbh, const int age, const int species, const int day, const int month, const int year);
//void Fruit_to_seeds_function_age(cell_t *const c, const int height, const int dbh, const int age, const int species, const int day, const int month, const int year);

#endif /* C_FRUIT_PARTITIONING_H_ */
