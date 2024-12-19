/*
 * recruitment.h
 *
 *  Created on: 17 april 2023
 *      Author: Vincenzo Saponaro
 */

#ifndef RECRUITMENT_H_
#define RECRUITMENT_H_

#include "matrix.h"

//static int fill_cell_for_regeneration (cell_t *const c);

int recruitment (cell_t *const c, species_t * s, const int day, const int month, const int year);

#endif /* RECRUITMENT_H_ */
