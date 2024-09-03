/*
 * CN-balance.c
 *
 *  Created on: 13 apr 2017
 *      Author: alessio
 */

// NOTE : CURRENTLY NOT USED

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "common.h"
#include "constants.h"
#include "settings.h"
#include "biomass.h"

extern settings_t* g_settings;
//extern logger_t* g_debug_log;

void carbon_balance (cell_t *const c, const int height, const int dbh, const int age, const int species)
{

	// note: currently not used.

	species_t *s;
	s = &c->heights[height].dbhs[dbh].ages[age].species[species];

 
	/*** update class level carbon mass pools ***/
	s->value[LEAF_C]     -= (s->value[C_LEAF_TO_LITR]  + s->value[C_LEAF_TO_RESERVE]);
	s->value[FROOT_C]    -= (s->value[C_FROOT_TO_LITR] + s->value[C_FROOT_TO_RESERVE]);
	s->value[STEM_C]     -= (s->value[C_STEM_TO_CWD]);
	s->value[CROOT_C]    -= (s->value[C_CROOT_TO_CWD]);
	s->value[BRANCH_C]   -= (s->value[C_BRANCH_TO_CWD]);
	s->value[RESERVE_C]  -= (s->value[C_RESERVE_TO_LITR]);
	s->value[FRUIT_C]    -= (s->value[C_FRUIT_TO_CWD]);
	s->value[LITR_C]     += (s->value[C_LEAF_TO_LITR]     +
			s->value[C_FROOT_TO_LITR]);
	s->value[CWD_C]      += (s->value[C_STEM_TO_CWD] +
			s->value[C_CROOT_TO_CWD]                 +
			s->value[C_BRANCH_TO_CWD]                +
					s->value[C_RESERVE_TO_LITR]               +
			s->value[C_FRUIT_TO_CWD]);

	/* check */
	CHECK_CONDITION ( s->value[LEAF_C],     < , ZERO );
	CHECK_CONDITION ( s->value[FROOT_C],    < , ZERO );
	CHECK_CONDITION ( s->value[STEM_C],     < , ZERO );
	CHECK_CONDITION ( s->value[BRANCH_C],   < , ZERO );
	CHECK_CONDITION ( s->value[CROOT_C],    < , ZERO );
	CHECK_CONDITION ( s->value[FRUIT_C],    < , ZERO );
	CHECK_CONDITION ( s->value[LITR_C],     < , ZERO );
	CHECK_CONDITION ( s->value[CWD_C],      < , ZERO );

	/*** update cell level carbon pools ***/

	/* tree */
	c->leaf_carbon              -= (s->value[C_LEAF_TO_LITR]  * 1e6 / g_settings->sizeCell);
	c->froot_carbon             -= (s->value[C_FROOT_TO_LITR] * 1e6 / g_settings->sizeCell);
	c->stem_carbon              -= (s->value[C_STEM_TO_CWD]   * 1e6 / g_settings->sizeCell);
	c->branch_carbon            -= (s->value[C_CROOT_TO_CWD]  * 1e6 / g_settings->sizeCell);
	c->croot_carbon             -= (s->value[C_BRANCH_TO_CWD] * 1e6 / g_settings->sizeCell);
	c->reserve_carbon           -= (s->value[C_BRANCH_TO_CWD] * 1e6 / g_settings->sizeCell);
	c->fruit_carbon             -= (s->value[C_FRUIT_TO_CWD]  * 1e6 / g_settings->sizeCell);

	/* check */
	CHECK_CONDITION ( c->leaf_carbon,    < , ZERO );
	CHECK_CONDITION ( c->froot_carbon,   < , ZERO );
	CHECK_CONDITION ( c->stem_carbon,    < , ZERO );
	CHECK_CONDITION ( c->branch_carbon,  < , ZERO );
	CHECK_CONDITION ( c->croot_carbon,   < , ZERO );
	CHECK_CONDITION ( c->fruit_carbon,   < , ZERO );
	CHECK_CONDITION ( c->litrC,          < , ZERO );
	CHECK_CONDITION ( c->cwd_C,          < , ZERO );
	/***************************************************************************************/
}

