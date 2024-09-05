/*
 * litter_model.c
 *
 *  Created on: 22 mar 2017
 *      Author: alessio
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "matrix.h"
#include "common.h"
#include "constants.h"
#include "logger.h"
#include "settings.h"
#include "het_respiration.h"
#include "litter_carbon_balance.h"
#include "litter_nitrogen_balance.h"
#include "decomposition.h"
#include "check_balance.h"
#include "littering.h"

extern settings_t* g_settings;

int Litter_model_daily (matrix_t *const m, const int cell, const int day, const int month, const int year)
{
	/* shortcuts */
	cell_t *c;
	meteo_daily_t *meteo_daily;

	/* assign shortcuts */
	c = &m->cells[cell];
	meteo_daily = &m->cells[cell].years[year].m[month].d[day];

	/* check parameters */
	assert(m);

	/* if spinup on than does littering */
	if ( g_settings->spinup )
	{
		spinup_littering ( c );
	}

	/* compute litter decomposition */
	litter_decomposition             ( c, meteo_daily );

	/* compute heterotrophic respiration from litter pool */
	litter_heterotrophic_respiration ( c, meteo_daily );

	/* compute soil carbon balance */
	litter_carbon_balance            ( c, year );

	/* compute soil nitrogen balance */
	// litter_nitrogen_balance          ( c, year );


	/*******************************************************************************************************/

	/* CHECK FOR BALANCE CLOSURE */

	/* CHECK FOR CARBON FLUX BALANCE CLOSURE */
	/* 1 */ if ( ! check_litter_carbon_flux_balance    ( c ) ) return 0;

	/* CHECK FOR CARBON MASS BALANCE CLOSURE */
	/* 2 */ if ( ! check_litter_carbon_mass_balance    ( c ) ) return 0;

	/* CHECK FOR NITROGEN FLUX BALANCE CLOSURE */
	// /* 1 */ if ( ! check_litter_nitrogen_flux_balance  ( c ) ) return 0;

	/* CHECK FOR NITROGEN MASS BALANCE CLOSURE */
	// /* 2 */ if ( ! check_litter_nitrogen_mass_balance  ( c ) ) return 0;

	/*******************************************************************************************************/

	/* ok */
	return 1;
}
