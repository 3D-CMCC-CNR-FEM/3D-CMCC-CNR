/*
 * littering.c
 *
 *  Created on: 12 apr 2017
 *      Author: alessio
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "common.h"
#include "leaf_fall.h"
#include "constants.h"
#include "settings.h"
#include "logger.h"

extern settings_t* g_settings;
extern logger_t* g_debug_log;

void littering ( cell_t *const c, species_t *const s )
{
	logger(g_debug_log, "\nLITTERING\n");

	/*********************************************************************************************************************************/
	/*** deadwood carbon ***/

	/*** compute class-level deadwood carbon fluxes (tC/sizecell/day) ****/
	s->value[CWD_TO_LITRC]         = s->value[C_TO_CWD];
	s->value[CWD_TO_LITR2C]        = s->value[C_TO_CWD]       * s->value[DEADWOOD_USCEL_FRAC];
	s->value[CWD_TO_LITR3C]        = s->value[C_TO_CWD]       * s->value[DEADWOOD_SCEL_FRAC];
	s->value[CWD_TO_LITR4C]        = s->value[C_TO_CWD]       * s->value[DEADWOOD_LIGN_FRAC];
	/* check */
	CHECK_CONDITION ( s->value[CWD_TO_LITR2C] + s->value[CWD_TO_LITR3C] + s->value[CWD_TO_LITR4C] , == , s->value[CWD_TO_LITRC] + eps );
    //printf("c->value[CWD_TO_LITRC]           = %g tC/cell/day\n", (s->value[CWD_TO_LITRC] * 1e6 / g_settings->sizeCell));

	/*** compute cell-level cwd carbon pools (gC/m2) ***/
	c->cwd_C                      += s->value[CWD_TO_LITRC]   * 1e6 / g_settings->sizeCell;
	c->cwd_2C                     += s->value[CWD_TO_LITR2C]  * 1e6 / g_settings->sizeCell;
	c->cwd_3C                     += s->value[CWD_TO_LITR3C]  * 1e6 / g_settings->sizeCell;
	c->cwd_4C                     += s->value[CWD_TO_LITR4C]  * 1e6 / g_settings->sizeCell;
	/* check */
	CHECK_CONDITION ( c->cwd_2C + c->cwd_3C + c->cwd_4C , == , c->cwd_C + eps );

	/*********************************************************************************************************************************/
	/*** leaf litter carbon ***/

	/*** compute leaf litter carbon fluxes (tC/sizecell/day) ****/
	s->value[LEAF_TO_LITRC]             = s->value[C_LEAF_TO_LITR];
	s->value[LEAF_TO_LITR1C]            = s->value[C_LEAF_TO_LITR]      * s->value[LEAF_LITR_LAB_FRAC] ;
	s->value[LEAF_TO_LITR2C]            = s->value[C_LEAF_TO_LITR]      * s->value[LEAF_LITR_USCEL_FRAC];
	s->value[LEAF_TO_LITR3C]            = s->value[C_LEAF_TO_LITR]      * s->value[LEAF_LITR_SCEL_FRAC];
	s->value[LEAF_TO_LITR4C]            = s->value[C_LEAF_TO_LITR]      * s->value[LEAF_LITR_LIGN_FRAC];
	/* check */
	CHECK_CONDITION ( s->value[LEAF_TO_LITR1C] + s->value[LEAF_TO_LITR2C] + s->value[LEAF_TO_LITR3C] + s->value[LEAF_TO_LITR4C] , == , s->value[LEAF_TO_LITRC] + eps );
  //  printf("c->value[LEAF_TO_LITRC]           = %g tC/cell/day\n", (s->value[LEAF_TO_LITRC]* 1e6 / g_settings->sizeCell));

	/*** update cell-level leaf litter carbon fluxes (gC/m2/day) ****/
	c->daily_leaf_to_litrC             += s->value[LEAF_TO_LITRC]       * 1e6 / g_settings->sizeCell;
	c->daily_leaf_to_litr1C            += s->value[LEAF_TO_LITR1C]      * 1e6 / g_settings->sizeCell;
	c->daily_leaf_to_litr2C            += s->value[LEAF_TO_LITR2C]      * 1e6 / g_settings->sizeCell;
	c->daily_leaf_to_litr3C            += s->value[LEAF_TO_LITR3C]      * 1e6 / g_settings->sizeCell;
	c->daily_leaf_to_litr4C            += s->value[LEAF_TO_LITR4C]      * 1e6 / g_settings->sizeCell;
	/* check */
	CHECK_CONDITION ( c->daily_leaf_to_litr1C + c->daily_leaf_to_litr2C + c->daily_leaf_to_litr3C + c->daily_leaf_to_litr4C , == , c->daily_leaf_to_litrC + eps );
   // printf("c->daily_leaf_to_litrC           = %g tC/cell/day\n",c->daily_leaf_to_litrC );

	/*** compute cell-level leaf litter carbon pools (gC/m2) ***/
	c->leaf_litrC                      += s->value[LEAF_TO_LITRC]       * 1e6 / g_settings->sizeCell;
	c->leaf_litr1C                     += s->value[LEAF_TO_LITR1C]      * 1e6 / g_settings->sizeCell;
	c->leaf_litr2C                     += s->value[LEAF_TO_LITR2C]      * 1e6 / g_settings->sizeCell;
	c->leaf_litr3C                     += s->value[LEAF_TO_LITR3C]      * 1e6 / g_settings->sizeCell;
	c->leaf_litr4C                     += s->value[LEAF_TO_LITR4C]      * 1e6 / g_settings->sizeCell;
	/* check */
	CHECK_CONDITION ( c->leaf_litr1C + c->leaf_litr2C + c->leaf_litr3C + c->leaf_litr4C , == , c->leaf_litrC + eps );

	/*********************************************************************************************************************************/
	/*** fine root litter carbon ***/

	/*** compute fine root litter carbon fluxes (tC/sizecell/day) ****/
	s->value[FROOT_TO_LITRC]            = s->value[C_FROOT_TO_LITR];
	s->value[FROOT_TO_LITR1C]           = s->value[C_FROOT_TO_LITR]     * s->value[FROOT_LITR_LAB_FRAC] ;
	s->value[FROOT_TO_LITR2C]           = s->value[C_FROOT_TO_LITR]     * s->value[FROOT_LITR_USCEL_FRAC];
	s->value[FROOT_TO_LITR3C]           = s->value[C_FROOT_TO_LITR]     * s->value[FROOT_LITR_SCEL_FRAC];
	s->value[FROOT_TO_LITR4C]           = s->value[C_FROOT_TO_LITR]     * s->value[FROOT_LITR_LIGN_FRAC];
	/* check */
	CHECK_CONDITION ( s->value[LEAF_TO_LITR1C] + s->value[LEAF_TO_LITR2C] + s->value[LEAF_TO_LITR3C] + s->value[LEAF_TO_LITR4C], == , s->value[LEAF_TO_LITRC] + eps );

   // printf("c->value[FROOT_TO_LITRC]           = %g tC/cell/day\n", (s->value[FROOT_TO_LITRC]* 1e6 / g_settings->sizeCell));


	/*** update cell-level fine root litter carbon fluxes (gC/m2/day) ****/
	c->daily_froot_to_litrC            += s->value[FROOT_TO_LITRC]      * 1e6 / g_settings->sizeCell;
	c->daily_froot_to_litr1C           += s->value[FROOT_TO_LITR1C]     * 1e6 / g_settings->sizeCell;
	c->daily_froot_to_litr2C           += s->value[FROOT_TO_LITR2C]     * 1e6 / g_settings->sizeCell;
	c->daily_froot_to_litr3C           += s->value[FROOT_TO_LITR3C]     * 1e6 / g_settings->sizeCell;
	c->daily_froot_to_litr4C           += s->value[FROOT_TO_LITR4C]     * 1e6 / g_settings->sizeCell;
	/* check */
	CHECK_CONDITION ( c->daily_froot_to_litr1C + c->daily_froot_to_litr2C + c->daily_froot_to_litr3C + c->daily_froot_to_litr4C , == , c->daily_froot_to_litrC + eps );
  // printf("c->daily_froot_to_litrC           = %g tC/cell/day\n",c->daily_froot_to_litrC );


	/*** compute cell-level fine root litter carbon pools (gC/m2) ***/
	c->froot_litrC                     += s->value[FROOT_TO_LITRC]      * 1e6 / g_settings->sizeCell;
	c->froot_litr1C                    += s->value[FROOT_TO_LITR1C]     * 1e6 / g_settings->sizeCell;
	c->froot_litr2C                    += s->value[FROOT_TO_LITR2C]     * 1e6 / g_settings->sizeCell;
	c->froot_litr3C                    += s->value[FROOT_TO_LITR3C]     * 1e6 / g_settings->sizeCell;
	c->froot_litr4C                    += s->value[FROOT_TO_LITR4C]     * 1e6 / g_settings->sizeCell;
	/* check */
	CHECK_CONDITION ( c->froot_litr1C + c->froot_litr2C + c->froot_litr3C + c->froot_litr4C , == , c->froot_litrC + eps );

    // reserve 
/*** compute researve litter carbon fluxes (tC/sizecell/day) ****/
	s->value[RESERVE_TO_LITRC]             = s->value[C_RESERVE_TO_LITR];
	s->value[RESERVE_TO_LITR1C]            = s->value[C_RESERVE_TO_LITR] ;

//	printf("c->value[RESERVE_TO_LITRC]            = %g tC/cell/day\n", (s->value[RESERVE_TO_LITRC]* 1e6 / g_settings->sizeCell));

	/*** update cell-level leaf litter carbon fluxes (gC/m2/day) ****/
	c->daily_reserve_to_litrC             += s->value[RESERVE_TO_LITRC]       * 1e6 / g_settings->sizeCell;
	c->daily_reserve_to_litr1C            += s->value[RESERVE_TO_LITR1C]      * 1e6 / g_settings->sizeCell;

	/* check */
	//CHECK_CONDITION ( c->daily_leaf_to_litr1C + c->daily_leaf_to_litr2C + c->daily_leaf_to_litr3C + c->daily_leaf_to_litr4C , == , c->daily_leaf_to_litrC + eps );

	/*** compute cell-level leaf litter carbon pools (gC/m2) ***/
	c->reserve_litrC                      += s->value[RESERVE_TO_LITRC]       * 1e6 / g_settings->sizeCell;
	c->reserve_litr1C                     += s->value[RESERVE_TO_LITR1C]      * 1e6 / g_settings->sizeCell;

	/* check */
	//CHECK_CONDITION ( c->leaf_litr1C + c->leaf_litr2C + c->leaf_litr3C + c->leaf_litr4C , == , c->leaf_litrC + eps );


	/*********************************************************************************************************************************/
	/*** deadwood nitrogen ***/

	/*** compute class-level deadwood nitrogen fluxes (tN/sizecell/day) ****/
	s->value[CWD_TO_LITRN]         = s->value[N_TO_CWD];
	s->value[CWD_TO_LITR2N]        = s->value[N_TO_CWD]       * s->value[DEADWOOD_USCEL_FRAC];
	s->value[CWD_TO_LITR3N]        = s->value[N_TO_CWD]       * s->value[DEADWOOD_SCEL_FRAC];
	s->value[CWD_TO_LITR4N]        = s->value[N_TO_CWD]       * s->value[DEADWOOD_LIGN_FRAC];
	/* check */
	CHECK_CONDITION ( s->value[CWD_TO_LITR2N] + s->value[CWD_TO_LITR3N] + s->value[CWD_TO_LITR4N], == , s->value[CWD_TO_LITRN] + eps );

	/*** update cell-level cwd nitrogen fluxes (gN/m2/day) ****/
	c->daily_cwd_to_litrN         += s->value[CWD_TO_LITRN]   * 1e6 / g_settings->sizeCell;
	c->daily_cwd_to_litr2N        += s->value[CWD_TO_LITR2N]  * 1e6 / g_settings->sizeCell;
	c->daily_cwd_to_litr3N        += s->value[CWD_TO_LITR3N]  * 1e6 / g_settings->sizeCell;
	c->daily_cwd_to_litr4N        += s->value[CWD_TO_LITR4N]  * 1e6 / g_settings->sizeCell;
	/* check */
	CHECK_CONDITION ( c->daily_cwd_to_litr2N + c->daily_cwd_to_litr3N + c->daily_cwd_to_litr4N , == , c->daily_cwd_to_litrN + eps );

	/*** compute cell-level cwd nitrogen pools (gN/m2) ***/
	c->cwd_N                      += s->value[CWD_TO_LITRN]   * 1e6 / g_settings->sizeCell;
	c->cwd_2N                     += s->value[CWD_TO_LITR2N]  * 1e6 / g_settings->sizeCell;
	c->cwd_3N                     += s->value[CWD_TO_LITR3N]  * 1e6 / g_settings->sizeCell;
	c->cwd_4N                     += s->value[CWD_TO_LITR4N]  * 1e6 / g_settings->sizeCell;

	/* check */
	CHECK_CONDITION ( c->cwd_2N + c->cwd_3N + c->cwd_4N , == , c->cwd_N + eps );

	/*********************************************************************************************************************************/
	/*** leaf litter nitrogen ***/

	/*** compute leaf litter nitrogen fluxes (tN/sizecell/day) ****/
	s->value[LEAF_TO_LITRN]             = s->value[N_LEAF_TO_LITR];
	s->value[LEAF_TO_LITR1N]            = s->value[N_LEAF_TO_LITR]      * s->value[LEAF_LITR_LAB_FRAC] ;
	s->value[LEAF_TO_LITR2N]            = s->value[N_LEAF_TO_LITR]      * s->value[LEAF_LITR_USCEL_FRAC];
	s->value[LEAF_TO_LITR3N]            = s->value[N_LEAF_TO_LITR]      * s->value[LEAF_LITR_SCEL_FRAC];
	s->value[LEAF_TO_LITR4N]            = s->value[N_LEAF_TO_LITR]      * s->value[LEAF_LITR_LIGN_FRAC];
	/* check */
	CHECK_CONDITION ( s->value[LEAF_TO_LITR1N] + s->value[LEAF_TO_LITR2N] + s->value[LEAF_TO_LITR3N] + s->value[LEAF_TO_LITR4N], == , s->value[LEAF_TO_LITRN] + eps );

	/*** update cell-level leaf litter nitrogen fluxes (gN/m2/day) ****/
	c->daily_leaf_to_litrN             += s->value[LEAF_TO_LITRN]       * 1e6 / g_settings->sizeCell;
	c->daily_leaf_to_litr1N            += s->value[LEAF_TO_LITR1N]      * 1e6 / g_settings->sizeCell;
	c->daily_leaf_to_litr2N            += s->value[LEAF_TO_LITR2N]      * 1e6 / g_settings->sizeCell;
	c->daily_leaf_to_litr3N            += s->value[LEAF_TO_LITR3N]      * 1e6 / g_settings->sizeCell;
	c->daily_leaf_to_litr4N            += s->value[LEAF_TO_LITR4N]      * 1e6 / g_settings->sizeCell;
	/* check */
	CHECK_CONDITION ( c->daily_leaf_to_litr1N + c->daily_leaf_to_litr2N + c->daily_leaf_to_litr3N + c->daily_leaf_to_litr4N , == , c->daily_leaf_to_litrN + eps );

	/*** compute cell-level leaf litter nitrogen pools (gN/m2) ***/
	c->leaf_litrN                      += s->value[LEAF_TO_LITRN]       * 1e6 / g_settings->sizeCell;
	c->leaf_litr1N                     += s->value[LEAF_TO_LITR1N]      * 1e6 / g_settings->sizeCell;
	c->leaf_litr2N                     += s->value[LEAF_TO_LITR2N]      * 1e6 / g_settings->sizeCell;
	c->leaf_litr3N                     += s->value[LEAF_TO_LITR3N]      * 1e6 / g_settings->sizeCell;
	c->leaf_litr4N                     += s->value[LEAF_TO_LITR4N]      * 1e6 / g_settings->sizeCell;
	/* check */
	CHECK_CONDITION ( c->leaf_litr1N + c->leaf_litr2N + c->leaf_litr3N + c->leaf_litr4N , == , c->leaf_litrN + eps );

	/*********************************************************************************************************************************/
	/*** fine root litter nitrogen ***/

	/*** compute fine root litter carbon fluxes (tN/sizecell/day) ****/
	s->value[FROOT_TO_LITRN]            = s->value[N_FROOT_TO_LITR];
	s->value[FROOT_TO_LITR1N]           = s->value[N_FROOT_TO_LITR]     * s->value[FROOT_LITR_LAB_FRAC] ;
	s->value[FROOT_TO_LITR2N]           = s->value[N_FROOT_TO_LITR]     * s->value[FROOT_LITR_USCEL_FRAC];
	s->value[FROOT_TO_LITR3N]           = s->value[N_FROOT_TO_LITR]     * s->value[FROOT_LITR_SCEL_FRAC];
	s->value[FROOT_TO_LITR4N]           = s->value[N_FROOT_TO_LITR]     * s->value[FROOT_LITR_LIGN_FRAC];
	/* check */
	CHECK_CONDITION ( s->value[FROOT_TO_LITR1N] + s->value[FROOT_TO_LITR2N] + s->value[FROOT_TO_LITR3N] + s->value[FROOT_TO_LITR4N], == , s->value[FROOT_TO_LITRN] + eps );

	/*** update cell-level fine root litter nitrogen fluxes (gN/m2/day) ****/
	c->daily_froot_to_litrN            += s->value[FROOT_TO_LITRN]      * 1e6 / g_settings->sizeCell;
	c->daily_froot_to_litr1N           += s->value[FROOT_TO_LITR1N]     * 1e6 / g_settings->sizeCell;
	c->daily_froot_to_litr2N           += s->value[FROOT_TO_LITR2N]     * 1e6 / g_settings->sizeCell;
	c->daily_froot_to_litr3N           += s->value[FROOT_TO_LITR3N]     * 1e6 / g_settings->sizeCell;
	c->daily_froot_to_litr4N           += s->value[FROOT_TO_LITR4N]     * 1e6 / g_settings->sizeCell;
	/* check */
	CHECK_CONDITION ( c->daily_froot_to_litr1N + c->daily_froot_to_litr2N + c->daily_froot_to_litr3N + c->daily_froot_to_litr4N , == , c->daily_froot_to_litrN + eps );

	/*** compute cell-level fine root litter nitrogen pools (gN/m2) ***/
	c->froot_litrN                     += s->value[FROOT_TO_LITRN]      * 1e6 / g_settings->sizeCell;
	c->froot_litr1N                    += s->value[FROOT_TO_LITR1N]     * 1e6 / g_settings->sizeCell;
	c->froot_litr2N                    += s->value[FROOT_TO_LITR2N]     * 1e6 / g_settings->sizeCell;
	c->froot_litr3N                    += s->value[FROOT_TO_LITR3N]     * 1e6 / g_settings->sizeCell;
	c->froot_litr4N                    += s->value[FROOT_TO_LITR4N]     * 1e6 / g_settings->sizeCell;
	/* check */
	CHECK_CONDITION ( c->froot_litr1N + c->froot_litr2N + c->froot_litr3N + c->froot_litr4N , == , c->froot_litrN + eps );

}

void spinup_littering ( cell_t *const c )
{
	double r1;
	double temp_var;
	double deadwood_cell_frac;
	double deadwood_uscel_frac;
	double deadwood_scel_frac;
	double deadwood_lign_frac;
	double leaf_litr_lab_frac;
	double leaf_ltr_cell_frac;
	double leaf_litr_uscel_frac;
	double leaf_litr_scel_frac;
	double leaf_litr_lign_frac;
	double froot_litr_lab_frac;
	double froot_ltr_cell_frac;
	double froot_litr_uscel_frac;
	double froot_litr_scel_frac;
	double froot_litr_lign_frac;
	double annual_C_to_cwd , annual_N_to_cwd;
	double annual_C_to_leaf_litr , annual_N_to_leaf_litr;
	double annual_C_to_froot_litr, annual_N_to_froot_litr;
	double C_to_cwd, N_to_cwd;
	double C_to_leaf_litr  , N_to_leaf_litr;
	double C_to_froot_litr , N_to_froot_litr;
	double cwd_to_litrC, cwd_to_litr2C, cwd_to_litr3C, cwd_to_litr4C;
	double leaf_to_litrC, leaf_to_litr1C, leaf_to_litr2C, leaf_to_litr3C, leaf_to_litr4C;
	double froot_to_litrC, froot_to_litr1C, froot_to_litr2C, froot_to_litr3C, froot_to_litr4C;
	double cwd_to_litrN, cwd_to_litr2N, cwd_to_litr3N, cwd_to_litr4N;
	double leaf_to_litrN, leaf_to_litr1N, leaf_to_litr2N, leaf_to_litr3N, leaf_to_litr4N;
	double froot_to_litrN, froot_to_litr1N, froot_to_litr2N, froot_to_litr3N, froot_to_litr4N;

	double cn_deadwood, cn_leaf_litr, cn_froot_litr;

	/* note assume general partitioning ratios values */

	deadwood_cell_frac   = 0.75;
	deadwood_lign_frac   = 0.25;
	/* check */
	CHECK_CONDITION ( fabs ( deadwood_cell_frac + deadwood_lign_frac ) - 1. , > ,  eps );

	leaf_litr_lab_frac   = 0.2;
	leaf_ltr_cell_frac   = 0.5;
	leaf_litr_lign_frac  = 0.3;
	/* check */
	CHECK_CONDITION ( fabs ( leaf_litr_lab_frac + leaf_ltr_cell_frac + leaf_litr_lign_frac ) -1 , > , eps );

	froot_litr_lab_frac  = 0.3;
	froot_ltr_cell_frac  = 0.4;
	froot_litr_lign_frac = 0.3;
	/* check */
	CHECK_CONDITION ( fabs (froot_litr_lab_frac + froot_ltr_cell_frac + froot_litr_lign_frac) -1.0 , > , eps );

	/******************************************************************************************************************************************************************/

	/* note assume general incoming fluxes of carbon values (tC/sizeCell/year) */
	/* annual daily values for littering */
	annual_C_to_cwd          = 1.;
	annual_C_to_leaf_litr    = 1.5;
	annual_C_to_froot_litr   = 1.5;

	/* convert to daily values */
	C_to_cwd                 = annual_C_to_cwd        / 365.;
	C_to_leaf_litr           = annual_C_to_leaf_litr  / 365.;
	C_to_froot_litr          = annual_C_to_froot_litr / 365.;

	/******************************************************************************************************************************************************************/

	/* note: assume general CN ratios */
	cn_deadwood              = 600;
	cn_leaf_litr             = 70;
	cn_froot_litr            = 60;

	/* calculate Nitrogen pools */
	annual_N_to_cwd          = annual_C_to_cwd        /  cn_deadwood;
	annual_N_to_leaf_litr    = annual_C_to_leaf_litr  /  cn_leaf_litr;
	annual_N_to_froot_litr   = annual_C_to_froot_litr /  cn_froot_litr;

	/* convert to daily values */
	N_to_cwd                 = annual_N_to_cwd        / 365.;
	N_to_leaf_litr           = annual_N_to_leaf_litr  / 365.;
	N_to_froot_litr          = annual_N_to_froot_litr / 365.;

	/******************************************************************************************************************************************************************/
	/*** compute coarse woody debris litter fractions ***/

	/* partitioning coarse woody debris litter carbon cellulose into shielded and unshielded pools */
	r1 = deadwood_lign_frac / deadwood_cell_frac;

	if ( r1 <= 0.45 )
	{
		deadwood_scel_frac  = 0.;
		deadwood_uscel_frac = deadwood_cell_frac;
	}
	else if ( r1 > 0.45 && r1 < 0.7 )
	{
		temp_var = ( r1 - 0.45 ) * 3.2;
		deadwood_scel_frac  = temp_var        * deadwood_cell_frac;
		deadwood_uscel_frac = (1. - temp_var) * deadwood_cell_frac;
	}
	else
	{
		deadwood_scel_frac  = 0.8 * deadwood_cell_frac;
		deadwood_uscel_frac = 0.2 * deadwood_cell_frac;
	}
	/* check */
	CHECK_CONDITION ( fabs ( deadwood_scel_frac + deadwood_uscel_frac + deadwood_lign_frac ), >, 1 + eps);

	/******************************************************************************************************************************************************************/
	/** calculate shielded and unshielded cellulose fraction for leaf pool **/

	/*** compute leaf litter fractions ***/

	/* partitioning leaf litter carbon cellulose into shielded and unshielded pools */
	r1 = leaf_litr_lign_frac / leaf_litr_lab_frac;

	if ( r1 <= 0.45 )
	{
		leaf_litr_scel_frac   = 0.;
		leaf_litr_uscel_frac  = leaf_ltr_cell_frac;
	}
	else if ( r1 > 0.45 && r1 < 0.7 )
	{
		temp_var = ( r1 - 0.45 ) * 3.2;
		leaf_litr_scel_frac   = temp_var          * leaf_ltr_cell_frac;
		leaf_litr_uscel_frac  = ( 1. - temp_var ) * leaf_ltr_cell_frac;
	}
	else
	{
		leaf_litr_scel_frac   = 0.8 * leaf_ltr_cell_frac;
		leaf_litr_uscel_frac  = 0.2 * leaf_ltr_cell_frac;
	}
	/* check */
	CHECK_CONDITION ( fabs ( leaf_litr_lab_frac + leaf_litr_scel_frac + leaf_litr_uscel_frac + leaf_litr_lign_frac ), >, 1 + eps);

	/******************************************************************************************************************************************************************/
	/** calculate shielded and unshielded cellulose fraction for froot pool **/

	/*** compute fine root litter fractions ***/

	/* partitioning fine root litter carbon cellulose into shielded and unshielded pools */
	r1 = froot_litr_lign_frac / froot_litr_lab_frac;

	if ( r1 <= 0.45 )
	{
		froot_litr_scel_frac   = 0.;
		froot_litr_uscel_frac  = froot_ltr_cell_frac;
	}
	else if ( r1 > 0.45 && r1 < 0.7 )
	{
		temp_var = ( r1 - 0.45 ) * 3.2;
		froot_litr_scel_frac   = temp_var          * froot_ltr_cell_frac;
		froot_litr_uscel_frac  = ( 1. - temp_var ) * froot_ltr_cell_frac;
	}
	else
	{
		froot_litr_scel_frac   = 0.8 * froot_ltr_cell_frac;
		froot_litr_uscel_frac  = 0.2 * froot_ltr_cell_frac;
	}
	/* check */
	CHECK_CONDITION ( fabs ( froot_litr_lab_frac + froot_litr_scel_frac + froot_litr_uscel_frac + froot_litr_lign_frac ), >, 1 + eps);

	/******************************************************************************************************************************************************************/

	/*********************************************************************************************************************************/
	/*** deadwood carbon ***/

	/*** compute class-level deadwood carbon fluxes (tC/sizecell/day) ****/
	cwd_to_litrC                  = C_to_cwd;
	cwd_to_litr2C                 = C_to_cwd       * deadwood_uscel_frac;
	cwd_to_litr3C                 = C_to_cwd       * deadwood_scel_frac;
	cwd_to_litr4C                 = C_to_cwd       * deadwood_lign_frac;
	/* check */
	CHECK_CONDITION ( cwd_to_litr2C + cwd_to_litr3C + cwd_to_litr4C , == , cwd_to_litrC + eps );

	/*** update cell-level cwd carbon fluxes (gC/m2/day) ****/
	c->daily_cwd_to_litrC         = cwd_to_litrC    * 1e6 / g_settings->sizeCell;
	c->daily_cwd_to_litr2C        = cwd_to_litr2C   * 1e6 / g_settings->sizeCell;
	c->daily_cwd_to_litr3C        = cwd_to_litr3C   * 1e6 / g_settings->sizeCell;
	c->daily_cwd_to_litr4C        = cwd_to_litr4C   * 1e6 / g_settings->sizeCell;
	/* check */
	CHECK_CONDITION ( c->daily_cwd_to_litr2C + c->daily_cwd_to_litr3C + c->daily_cwd_to_litr4C , == , c->daily_cwd_to_litrC + eps );

	/*** compute cell-level cwd carbon pools (gC/m2) ***/
	c->cwd_C                      += c->daily_cwd_to_litrC;
	c->cwd_2C                     += c->daily_cwd_to_litr2C;
	c->cwd_3C                     += c->daily_cwd_to_litr3C;
	c->cwd_4C                     += c->daily_cwd_to_litr4C;
	/* check */
	CHECK_CONDITION ( c->cwd_2C + c->cwd_3C + c->cwd_4C , == , c->cwd_C + eps );

	/*** deadwood nitrogen ***/

	/*** compute class-level deadwood nitrogen fluxes (tN/sizecell/day) ****/
	cwd_to_litrN                  = N_to_cwd;
	cwd_to_litr2N                 = N_to_cwd       * deadwood_uscel_frac;
	cwd_to_litr3N                 = N_to_cwd       * deadwood_scel_frac;
	cwd_to_litr4N                 = N_to_cwd       * deadwood_lign_frac;
	/* check */
	CHECK_CONDITION ( cwd_to_litr2N + cwd_to_litr3N + cwd_to_litr4N , == , cwd_to_litrN + eps );

	/*** update cell-level cwd nitrogen fluxes (gN/m2/day) ****/
	c->daily_cwd_to_litrN         = cwd_to_litrN    * 1e6 / g_settings->sizeCell;
	c->daily_cwd_to_litr2N        = cwd_to_litr2N   * 1e6 / g_settings->sizeCell;
	c->daily_cwd_to_litr3N        = cwd_to_litr3N   * 1e6 / g_settings->sizeCell;
	c->daily_cwd_to_litr4N        = cwd_to_litr4N   * 1e6 / g_settings->sizeCell;
	/* check */
	CHECK_CONDITION ( c->daily_cwd_to_litr2N + c->daily_cwd_to_litr3N + c->daily_cwd_to_litr4N , == , c->daily_cwd_to_litrN + eps );

	/*** compute cell-level cwd nitrogen pools (gN/m2) ***/
	c->cwd_N                      += c->daily_cwd_to_litrN;
	c->cwd_2N                     += c->daily_cwd_to_litr2N;
	c->cwd_3N                     += c->daily_cwd_to_litr3N;
	c->cwd_4N                     += c->daily_cwd_to_litr4N;
	/* check */
	CHECK_CONDITION ( c->cwd_2N + c->cwd_3N + c->cwd_4N , == , c->cwd_N + eps );

	/*********************************************************************************************************************************/

	/*** leaf litter carbon ***/

	/*** compute leaf litter carbon fluxes (tC/sizecell/day) ****/
	leaf_to_litrC                       = C_to_leaf_litr;
	leaf_to_litr1C                      = C_to_leaf_litr      * leaf_litr_lab_frac ;
	leaf_to_litr2C                      = C_to_leaf_litr      * leaf_litr_uscel_frac;
	leaf_to_litr3C                      = C_to_leaf_litr      * leaf_litr_scel_frac;
	leaf_to_litr4C                      = C_to_leaf_litr      * leaf_litr_lign_frac;
	/* check */
	CHECK_CONDITION ( leaf_to_litr1C + leaf_to_litr2C + leaf_to_litr3C + leaf_to_litr4C , == , leaf_to_litrC + eps );

	/*** update cell-level leaf litter carbon fluxes (gC/m2/day) ****/
	c->daily_leaf_to_litrC              = leaf_to_litrC       * 1e6 / g_settings->sizeCell;
	c->daily_leaf_to_litr1C             = leaf_to_litr1C      * 1e6 / g_settings->sizeCell;
	c->daily_leaf_to_litr2C             = leaf_to_litr2C      * 1e6 / g_settings->sizeCell;
	c->daily_leaf_to_litr3C             = leaf_to_litr3C      * 1e6 / g_settings->sizeCell;
	c->daily_leaf_to_litr4C             = leaf_to_litr4C      * 1e6 / g_settings->sizeCell;
	/* check */
	CHECK_CONDITION ( c->daily_leaf_to_litr1C + c->daily_leaf_to_litr2C + c->daily_leaf_to_litr3C + c->daily_leaf_to_litr4C , == , c->daily_leaf_to_litrC + eps );

	/*** compute cell-level leaf litter carbon pools (gC/m2) ***/
	c->leaf_litrC                      += c->daily_leaf_to_litrC;
	c->leaf_litr1C                     += c->daily_leaf_to_litr1C;
	c->leaf_litr2C                     += c->daily_leaf_to_litr2C;
	c->leaf_litr3C                     += c->daily_leaf_to_litr3C;
	c->leaf_litr4C                     += c->daily_leaf_to_litr4C;
	/* check */
	CHECK_CONDITION ( c->leaf_litr1C + c->leaf_litr2C + c->leaf_litr3C + c->leaf_litr4C , == , c->leaf_litrC + eps );

	/*** leaf litter nitrogen ***/

	/*** compute leaf litter nitrogen fluxes (tN/sizecell/day) ****/
	leaf_to_litrN                       = N_to_leaf_litr;
	leaf_to_litr1N                      = N_to_leaf_litr      * leaf_litr_lab_frac ;
	leaf_to_litr2N                      = N_to_leaf_litr      * leaf_litr_uscel_frac;
	leaf_to_litr3N                      = N_to_leaf_litr      * leaf_litr_scel_frac;
	leaf_to_litr4N                      = N_to_leaf_litr      * leaf_litr_lign_frac;
	/* check */
	CHECK_CONDITION ( leaf_to_litr1N + leaf_to_litr2N + leaf_to_litr3N + leaf_to_litr4N , == , leaf_to_litrN + eps );

	/*** update cell-level leaf litter nitrogen fluxes (gN/m2/day) ****/
	c->daily_leaf_to_litrN              = leaf_to_litrN       * 1e6 / g_settings->sizeCell;
	c->daily_leaf_to_litr1N             = leaf_to_litr1N      * 1e6 / g_settings->sizeCell;
	c->daily_leaf_to_litr2N             = leaf_to_litr2N      * 1e6 / g_settings->sizeCell;
	c->daily_leaf_to_litr3N             = leaf_to_litr3N      * 1e6 / g_settings->sizeCell;
	c->daily_leaf_to_litr4N             = leaf_to_litr4N      * 1e6 / g_settings->sizeCell;
	/* check */
	CHECK_CONDITION ( c->daily_leaf_to_litr1N + c->daily_leaf_to_litr2N + c->daily_leaf_to_litr3N + c->daily_leaf_to_litr4N , == , c->daily_leaf_to_litrN + eps );

	/*** compute cell-level leaf litter nitrogen pools (gN/m2) ***/
	c->leaf_litrN                      += c->daily_leaf_to_litrN;
	c->leaf_litr1N                     += c->daily_leaf_to_litr1N;
	c->leaf_litr2N                     += c->daily_leaf_to_litr2N;
	c->leaf_litr3N                     += c->daily_leaf_to_litr3N;
	c->leaf_litr4N                     += c->daily_leaf_to_litr4N;
	/* check */
	CHECK_CONDITION ( c->leaf_litr1N + c->leaf_litr2N + c->leaf_litr3N + c->leaf_litr4N , == , c->leaf_litrN + eps );

	/*********************************************************************************************************************************/

	/*** froot litter carbon ***/

	/*** compute froot litter carbon fluxes (tC/sizecell/day) ****/
	froot_to_litrC                       = C_to_froot_litr;
	froot_to_litr1C                      = C_to_froot_litr      * froot_litr_lab_frac ;
	froot_to_litr2C                      = C_to_froot_litr      * froot_litr_uscel_frac;
	froot_to_litr3C                      = C_to_froot_litr      * froot_litr_scel_frac;
	froot_to_litr4C                      = C_to_froot_litr      * froot_litr_lign_frac;
	/* check */
	CHECK_CONDITION ( froot_to_litr1C + froot_to_litr2C + froot_to_litr3C + froot_to_litr4C , == , froot_to_litrC + eps );

	/*** update cell-level froot litter carbon fluxes (gC/m2/day) ****/
	c->daily_froot_to_litrC              = froot_to_litrC       * 1e6 / g_settings->sizeCell;
	c->daily_froot_to_litr1C             = froot_to_litr1C      * 1e6 / g_settings->sizeCell;
	c->daily_froot_to_litr2C             = froot_to_litr2C      * 1e6 / g_settings->sizeCell;
	c->daily_froot_to_litr3C             = froot_to_litr3C      * 1e6 / g_settings->sizeCell;
	c->daily_froot_to_litr4C             = froot_to_litr4C      * 1e6 / g_settings->sizeCell;
	/* check */
	CHECK_CONDITION ( c->daily_froot_to_litr1C + c->daily_froot_to_litr2C + c->daily_froot_to_litr3C + c->daily_froot_to_litr4C , == , c->daily_froot_to_litrC + eps );

	/*** compute cell-level froot litter carbon pools (gC/m2) ***/
	c->froot_litrC                      += c->daily_froot_to_litrC;
	c->froot_litr1C                     += c->daily_froot_to_litr1C;
	c->froot_litr2C                     += c->daily_froot_to_litr2C;
	c->froot_litr3C                     += c->daily_froot_to_litr3C;
	c->froot_litr4C                     += c->daily_froot_to_litr4C;
	/* check */
	CHECK_CONDITION ( c->froot_litr1C + c->froot_litr2C + c->froot_litr3C + c->froot_litr4C , == , c->froot_litrC + eps );

	/*** froot litter nitrogen ***/

	/*** compute froot litter nitrogen fluxes (tN/sizeCell/day) ****/
	froot_to_litrN                       = N_to_froot_litr;
	froot_to_litr1N                      = N_to_froot_litr      * froot_litr_lab_frac ;
	froot_to_litr2N                      = N_to_froot_litr      * froot_litr_uscel_frac;
	froot_to_litr3N                      = N_to_froot_litr      * froot_litr_scel_frac;
	froot_to_litr4N                      = N_to_froot_litr      * froot_litr_lign_frac;
	/* check */
	CHECK_CONDITION ( froot_to_litr1N + froot_to_litr2N + froot_to_litr3N + froot_to_litr4N , == , froot_to_litrN + eps );

	/*** update cell-level froot litter nitrogen fluxes (gN/m2/day) ****/
	c->daily_froot_to_litrN              = froot_to_litrN       * 1e6 / g_settings->sizeCell;
	c->daily_froot_to_litr1N             = froot_to_litr1N      * 1e6 / g_settings->sizeCell;
	c->daily_froot_to_litr2N             = froot_to_litr2N      * 1e6 / g_settings->sizeCell;
	c->daily_froot_to_litr3N             = froot_to_litr3N      * 1e6 / g_settings->sizeCell;
	c->daily_froot_to_litr4N             = froot_to_litr4N      * 1e6 / g_settings->sizeCell;
	/* check */
	CHECK_CONDITION ( c->daily_froot_to_litr1N + c->daily_froot_to_litr2N + c->daily_froot_to_litr3N + c->daily_froot_to_litr4N , == , c->daily_froot_to_litrN + eps );

	/*** compute cell-level froot litter nitrogen pools (gN/m2) ***/
	c->froot_litrN                      += c->daily_froot_to_litrN;
	c->froot_litr1N                     += c->daily_froot_to_litr1N;
	c->froot_litr2N                     += c->daily_froot_to_litr2N;
	c->froot_litr3N                     += c->daily_froot_to_litr3N;
	c->froot_litr4N                     += c->daily_froot_to_litr4N;
	/* check */
	CHECK_CONDITION ( c->froot_litr1N + c->froot_litr2N + c->froot_litr3N + c->froot_litr4N , == , c->froot_litrN + eps );

	/*********************************************************************************************************************************/



}
