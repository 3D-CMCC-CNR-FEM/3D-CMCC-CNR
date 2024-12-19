/*
 * regeneration.h
 *
 *  Created on: 15 set 2016
 *      Author: alessio
 */

#ifndef HEADERS_REGENERATION_H_
#define HEADERS_REGENERATION_H_

#include "matrix.h"

int germination (cell_t *const c, const meteo_daily_t *const meteo_daily, species_t *const s, const int day, const int month, const int year);

int establishment (cell_t *const c, const meteo_daily_t *const meteo_daily, species_t *const s, const int day, const int month, const int year);

#endif /* HEADERS_REGENERATION_H_ */
