/*
 * canopy_radiation_sw_band.h
 *
 *  Created on: 05 lug 2016
 *      Author: alessio
 */

#ifndef HEADERS_CANOPY_RADIATION_SW_BAND_H_
#define HEADERS_CANOPY_RADIATION_SW_BAND_H_

#include "matrix.h"

void canopy_sw_band_abs_trans_refl_radiation(cell_t *const , const int , const int , const int , const int , const int , const meteo_daily_t *const,
		double, double, double, double, double, double, double, double, double);

void canopy_radiation_sw_band(cell_t *const c, const int layer, const int height, const int dbh, const int age, const int species,  meteo_daily_t *const meteo_daily);

#endif /* HEADERS_CANOPY_RADIATION_SW_BAND_H_ */
