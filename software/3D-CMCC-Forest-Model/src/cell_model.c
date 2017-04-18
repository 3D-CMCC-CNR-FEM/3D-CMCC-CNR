/*
 * cell_model.c
 *
 *  Created on: 31 ago 2016
 *      Author: alessio
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "matrix.h"
#include "constants.h"
#include "logger.h"
#include "canopy_radiation_lw_band.h"
#include "canopy_radiation_sw_band.h"
#include "evapotranspiration.h"
#include "heat_fluxes.h"
#include "fluxes.h"
#include "check_balance.h"

extern logger_t* g_debug_log;

int Cell_model_daily (matrix_t *const m, const int cell, const int day, const int month, const int year)
{
	cell_t *c;
	meteo_daily_t *meteo_daily;
	yos_t *meteo_yearly;

	/* assign shortcuts */
	c = &m->cells[cell];
	meteo_daily  = &m->cells[cell].years[year].m[month].d[day];
	meteo_yearly = &m->cells[cell].years[year];

	/* check parameters */
	assert(m);

	logger (g_debug_log, "\n*********CELL_MODEL_DAILY*********\n");

	/* OVERALL CELL */

	/* compute cell level carbon fluxes */
	//carbon_fluxes        ( c );

	/* compute cell level evapotranspiration */
	evapotranspiration   ( c );

	/* compute cell level latent heat fluxes */
	latent_heat_fluxes   ( c, meteo_daily );

	/* compute cell level sensible heat fluxes */
	sensible_heat_fluxes ( c, meteo_daily );

	/* compute cell level water fluxes */
	water_fluxes         ( c, meteo_daily );

	/*******************************************************************************************************/

	/* CHECK FOR BALANCE CLOSURE */

	/* CHECK FOR RADIATIVE FLUX BALANCE CLOSURE */
	/* 1 */ if ( ! check_cell_radiation_flux_balance ( c, meteo_daily ) ) return 0;

	/* CHECK FOR CARBON FLUX BALANCE CLOSURE */
	/* 2 */ if ( ! check_cell_carbon_flux_balance    ( c ) ) return 0;

	/* CHECK FOR NITROGEN FLUX BALANCE CLOSURE */
	/* 3 */ //fixme if ( ! check_cell_nitrogen_flux_balance    ( c ) ) return 0;

	/* CHECK FOR CARBON MASS BALANCE CLOSURE */
	/* 4 */ if ( ! check_cell_carbon_mass_balance    ( c ) ) return 0;

	/* CHECK FOR WATER FLUX BALANCE CLOSURE */
	/* 5 */ if ( ! check_cell_water_flux_balance     ( c, meteo_daily ) ) return 0;

	/* CHECK FOR NITROGEN MASS BALANCE CLOSURE */
	/* 6 */ //fixme if ( ! check_cell_nitrogen_mass_balance    ( c, meteo_yearly ) ) return 0;

	/*******************************************************************************************************/

	/* ok */
	return 1;

}
