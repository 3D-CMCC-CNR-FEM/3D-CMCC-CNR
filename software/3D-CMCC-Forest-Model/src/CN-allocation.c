/*
 * CN-allocation.c
 *
 * Created on: 01 lug 2016
 *   Author: alessio
 */
/* includes */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "common.h"
#include "constants.h"
#include "settings.h"
#include "logger.h"
#include "lai.h"
#include "turnover.h"
#include "dendometry.h"
#include "biomass.h"
#include "check_balance.h"
#include "mortality.h"
#include "remove_tree_class.h"

extern settings_t* g_settings;
extern logger_t* g_debug_log;


void carbon_allocation ( cell_t *const c, age_t *const a, species_t *const s, const int day, const int month, const int year )
{

	/* it allocates daily assimilated carbon for both deciduous and evergreen daily
	 * and removes the respired and dead parts */

	s->value[OLD_RESERVE_C] = s->value[RESERVE_C];

	logger(g_debug_log, "\n**CARBON ALLOCATION**\n");

	/*** removing growth respiration from carbon flux pools ***/
	s->value[C_TO_LEAF]    -= (s->value[LEAF_GROWTH_RESP]   / 1e6 * g_settings->sizeCell);
	s->value[C_TO_FROOT]   -= (s->value[FROOT_GROWTH_RESP]  / 1e6 * g_settings->sizeCell);
	s->value[C_TO_STEM]    -= (s->value[STEM_GROWTH_RESP]   / 1e6 * g_settings->sizeCell);
	s->value[C_TO_CROOT]   -= (s->value[CROOT_GROWTH_RESP]  / 1e6 * g_settings->sizeCell);
	s->value[C_TO_BRANCH]  -= (s->value[BRANCH_GROWTH_RESP] / 1e6 * g_settings->sizeCell);

	/***************************************************************************************/

	/* stem */
	s->value[STEM_SAPWOOD_C]          += s->value[C_TO_STEM] - s->value[C_STEM_SAPWOOD_TO_CWD];


	/* live stem */
	s->value[STEM_LIVEWOOD_C]          = s->value[STEM_SAPWOOD_C] * s->value[EFF_LIVE_TOTAL_WOOD_FRAC];

	/* to avoid that self-thinning mortality happens to remove too much biomass */
	if ( day && month ) s->value[YEARLY_C_TO_STEM] += s->value[C_TO_STEM];

	/***************************************************************************************/

	/* coarse root */
	s->value[CROOT_SAPWOOD_C]         += s->value[C_TO_CROOT] - s->value[C_CROOT_SAPWOOD_TO_CWD];

	/* live coarse root */
	s->value[CROOT_LIVEWOOD_C]         = s->value[CROOT_SAPWOOD_C] * s->value[EFF_LIVE_TOTAL_WOOD_FRAC];

	/* to avoid that self-thinning mortality happens to remove too much biomass */
	if ( day && month ) s->value[YEARLY_C_TO_CROOT] += s->value[C_TO_CROOT];

	/***************************************************************************************/
	/* branch */
	s->value[BRANCH_SAPWOOD_C]        += s->value[C_TO_BRANCH] - s->value[C_BRANCH_SAPWOOD_TO_CWD];

	/* live branch */
	s->value[BRANCH_LIVEWOOD_C]        = s->value[BRANCH_SAPWOOD_C] * s->value[EFF_LIVE_TOTAL_WOOD_FRAC];

	/* to avoid that self-thinning mortality happens to remove too much biomass */
	if ( day && month ) s->value[YEARLY_C_TO_BRANCH]      += s->value[C_TO_BRANCH];

	/***************************************************************************************/

	/*** removing dead pools from carbon flux pools ***/
	s->value[C_TO_LEAF]    -= s->value[C_LEAF_TO_LITR]  + s->value[C_LEAF_TO_RESERVE];
	s->value[C_TO_FROOT]   -= s->value[C_FROOT_TO_LITR] + s->value[C_FROOT_TO_RESERVE];
	s->value[C_TO_STEM]    -= s->value[C_STEM_TO_CWD];
	s->value[C_TO_CROOT]   -= s->value[C_CROOT_TO_CWD];
	s->value[C_TO_BRANCH]  -= s->value[C_BRANCH_TO_CWD];
	s->value[C_TO_FRUIT]   -= s->value[C_FRUIT_TO_CWD];
	s->value[C_TO_RESERVE] -= s->value[C_RESERVE_TO_CWD];

	/***************************************************************************************/

	/* for stocking wood */
	/* note: if here and NOT above model doesn't take into account dead parts */
	if ( s->value[C_TO_STEM]   > 0. ) s->value[YEARLY_C_TO_WOOD]     += s->value[C_TO_STEM];
	if ( s->value[C_TO_CROOT]  > 0. ) s->value[YEARLY_C_TO_WOOD]     += s->value[C_TO_BRANCH];
	if ( s->value[C_TO_BRANCH] > 0. ) s->value[YEARLY_C_TO_WOOD]     += s->value[C_TO_CROOT];
	if ( s->value[C_TO_STEM]   > 0. ) s->value[CUM_YEARLY_C_TO_WOOD] += s->value[C_TO_STEM];
	if ( s->value[C_TO_CROOT]  > 0. ) s->value[CUM_YEARLY_C_TO_WOOD] += s->value[C_TO_BRANCH];
	if ( s->value[C_TO_BRANCH] > 0. ) s->value[CUM_YEARLY_C_TO_WOOD] += s->value[C_TO_CROOT];

	/***************************************************************************************/
//
//	/* update Leaf Area Index */
//	daily_lai             ( c, a, s );
//
//	/***************************************************************************************/
//
//	/* computing leaf sun and shaded Carbon pools (tC/ha) */
//	if ( ! s->value[LEAF_C] )
//	{
//		s->value[LEAF_SUN_C]   = 0.;
//		s->value[LEAF_SHADE_C] = 0.;
//	}
//	else
//	{
//		/* compute based on SLA Leaf carbon for sun and shaded leaves */
//		s->value[LEAF_SUN_C]   = ( ( s->value[LAI_SUN_PROJ]   * ( s->value[CANOPY_COVER_EXP] * g_settings->sizeCell ) ) / s->value[SLA_SUN_PROJ] )   / 1e3;
//		s->value[LEAF_SHADE_C] = ( ( s->value[LAI_SHADE_PROJ] * ( s->value[CANOPY_COVER_EXP] * g_settings->sizeCell ) ) / s->value[SLA_SHADE_PROJ] ) / 1e3;
//	}

	/***************************************************************************************/

	/*** update class level carbon mass pools ***/
	s->value[LEAF_C]      += s->value[C_TO_LEAF];
	s->value[FROOT_C]     += s->value[C_TO_FROOT];
	s->value[STEM_C]      += s->value[C_TO_STEM];
	s->value[CROOT_C]     += s->value[C_TO_CROOT];
	s->value[BRANCH_C]    += s->value[C_TO_BRANCH];
	s->value[RESERVE_C]   += s->value[C_TO_RESERVE];
	s->value[FRUIT_C]     += s->value[C_TO_FRUIT];

	/*** update live carbon mass pools **/
	s->value[TOT_LIVEWOOD_C]           = ( s->value[STEM_LIVEWOOD_C] + s->value[CROOT_LIVEWOOD_C] + s->value[BRANCH_LIVEWOOD_C] );
	s->value[TOT_SAPWOOD_C]            = ( s->value[STEM_SAPWOOD_C]  + s->value[CROOT_SAPWOOD_C]  + s->value[BRANCH_SAPWOOD_C] );


	/*** update heartwood and dead carbon mass pools **/
	s->value[STEM_HEARTWOOD_C]         = s->value[STEM_C]   - s->value[STEM_SAPWOOD_C];
	s->value[STEM_DEADWOOD_C]          = s->value[STEM_C]   - s->value[STEM_LIVEWOOD_C];
	s->value[CROOT_HEARTWOOD_C]        = s->value[CROOT_C]  - s->value[CROOT_SAPWOOD_C];
	s->value[CROOT_DEADWOOD_C]         = s->value[CROOT_C]  - s->value[CROOT_LIVEWOOD_C];
	s->value[BRANCH_HEARTWOOD_C]       = s->value[BRANCH_C] - s->value[BRANCH_SAPWOOD_C];
	s->value[BRANCH_DEADWOOD_C]        = s->value[BRANCH_C] - s->value[BRANCH_LIVEWOOD_C];
	s->value[TOT_DEADWOOD_C]           = (s->value[STEM_DEADWOOD_C] + s->value[CROOT_DEADWOOD_C] + s->value[BRANCH_DEADWOOD_C] );

	/***************************************************************************************/

	/* check for closure */
	CHECK_CONDITION(fabs((s->value[STEM_SAPWOOD_C]   + s->value[STEM_HEARTWOOD_C])  - s->value[STEM_C]),  >,eps);
	CHECK_CONDITION(fabs((s->value[CROOT_SAPWOOD_C]  + s->value[CROOT_HEARTWOOD_C]) - s->value[CROOT_C]), >,eps);
	CHECK_CONDITION(fabs((s->value[BRANCH_SAPWOOD_C] + s->value[BRANCH_HEARTWOOD_C])- s->value[BRANCH_C]),>,eps);

	/* check */
	CHECK_CONDITION ( s->value[LEAF_C],     < , ZERO );
	CHECK_CONDITION ( s->value[FROOT_C],    < , ZERO );
	CHECK_CONDITION ( s->value[STEM_C],     < , ZERO );
	CHECK_CONDITION ( s->value[BRANCH_C],   < , ZERO );
	CHECK_CONDITION ( s->value[CROOT_C],    < , ZERO );
	CHECK_CONDITION ( s->value[FRUIT_C],    < , ZERO );
	/* sub pools */
	CHECK_CONDITION ( s->value[STEM_SAPWOOD_C],     < , ZERO );
	CHECK_CONDITION ( s->value[STEM_HEARTWOOD_C],   < , ZERO );
	CHECK_CONDITION ( s->value[CROOT_SAPWOOD_C],    < , ZERO );
	CHECK_CONDITION ( s->value[CROOT_HEARTWOOD_C],  < , ZERO );
	CHECK_CONDITION ( s->value[BRANCH_SAPWOOD_C],   < , ZERO );
	CHECK_CONDITION ( s->value[BRANCH_HEARTWOOD_C], < , ZERO );

	s->value[TOTAL_C] = s->value[LEAF_C] +
			s->value[FROOT_C]            +
			s->value[STEM_C]             +
			s->value[BRANCH_C]           +
			s->value[CROOT_C]            +
			s->value[FRUIT_C]            +
			s->value[RESERVE_C]          ;

	/* compute maximum annual amount to fruit */
	if ( s->value[C_TO_FRUIT] > 0. )
	{
		/* special case for fruit */
		s->value[MAX_FRUIT_C] += s->value[C_TO_FRUIT];
	}

	/* compute maximum and minimum annual NSC concentration */
	if ( s->value[RESERVE_C] > s->value[OLD_RESERVE_C] )
	{
		s->value[MAX_RESERVE_C]      = s->value[RESERVE_C];

		s->value[MAX_RESERVE_C_CONC] = ( s->value[MAX_RESERVE_C] / (s->value[TOT_SAPWOOD_C] * GC_GDM ) ) * 100.;
	}
	else
	{
		s->value[MIN_RESERVE_C_CONC] = ( s->value[RESERVE_C]     / (s->value[TOT_SAPWOOD_C] * GC_GDM ) ) * 100.;
	}

	/*** update cell level carbon fluxes (gC/m2/day)***/
	c->daily_leaf_carbon        += (s->value[C_TO_LEAF]    * 1e6 / g_settings->sizeCell);
	c->daily_froot_carbon       += (s->value[C_TO_FROOT]   * 1e6 / g_settings->sizeCell);
	c->daily_stem_carbon        += (s->value[C_TO_STEM]    * 1e6 / g_settings->sizeCell);
	c->daily_croot_carbon       += (s->value[C_TO_CROOT]   * 1e6 / g_settings->sizeCell);
	c->daily_branch_carbon      += (s->value[C_TO_BRANCH]  * 1e6 / g_settings->sizeCell);
	c->daily_reserve_carbon     += (s->value[C_TO_RESERVE] * 1e6 / g_settings->sizeCell);
	c->daily_fruit_carbon       += (s->value[C_TO_FRUIT]   * 1e6 / g_settings->sizeCell);

	/*** update cell level carbon pools (gC/m2)***/
	c->leaf_carbon              += (s->value[C_TO_LEAF]    * 1e6 / g_settings->sizeCell);
	c->froot_carbon             += (s->value[C_TO_FROOT]   * 1e6 / g_settings->sizeCell);
	c->stem_carbon              += (s->value[C_TO_STEM]    * 1e6 / g_settings->sizeCell);
	c->branch_carbon            += (s->value[C_TO_BRANCH]  * 1e6 / g_settings->sizeCell);
	c->croot_carbon             += (s->value[C_TO_CROOT]   * 1e6 / g_settings->sizeCell);
	c->reserve_carbon           += (s->value[C_TO_RESERVE] * 1e6 / g_settings->sizeCell);
	c->fruit_carbon             += (s->value[C_TO_FRUIT]   * 1e6 / g_settings->sizeCell);


	/* check */
	CHECK_CONDITION ( c->leaf_carbon,    < , ZERO );
	CHECK_CONDITION ( c->froot_carbon,   < , ZERO );
	CHECK_CONDITION ( c->stem_carbon,    < , ZERO );
	CHECK_CONDITION ( c->branch_carbon,  < , ZERO );
	CHECK_CONDITION ( c->croot_carbon,   < , ZERO );
	CHECK_CONDITION ( c->fruit_carbon,   < , ZERO );

}

/********************************************************************************************************************************************/

void nitrogen_allocation ( cell_t *const c, species_t *const s, const int day, const int month, const int year )
{
	double n_to_leaf = 0.;
	double n_to_froot = 0.;
	double n_to_stem = 0.;
	double n_to_croot = 0.;
	double n_to_branch = 0.;
	double n_to_fruit = 0.;

	//TODO IT SHOULD ALLOCATE NITROGEN ONCE EFFECTIVE NITROGEN AVAILABILITY IS COMPUTED

	/* it allocates Daily assimilated Nitrogen for both deciduous and evergreen and compute Nitrogen demand */

	logger(g_debug_log, "\n**NITROGEN ALLOCATION**\n");

	/*** update class level carbon Nitrogen pools ***/

	/*************************************************************************/

	/* leaf */
	s->value[N_TO_LEAF]   = s->value[C_TO_LEAF]   / s->value[CN_LEAVES];

	if ( s->value[N_TO_LEAF]   > 0. ) n_to_leaf   = s->value[N_TO_LEAF];
	else                              n_to_leaf   = 0.;

	/*************************************************************************/

	/* fine root */
	s->value[N_TO_FROOT]  = s->value[C_TO_FROOT]  / s->value[CN_FINE_ROOTS];

	if ( s->value[N_TO_FROOT]  > 0. ) n_to_froot  = s->value[N_TO_FROOT];
	else                              n_to_froot  = 0.;

	/*************************************************************************/

	/* fruit */
	s->value[N_TO_FRUIT]   = s->value[C_TO_FRUIT] / s->value[CN_LEAVES];//FIXME IT USES CN_LEAVES INSTEAD A CN_FRUITS

	if ( s->value[N_TO_FRUIT]  > 0. ) n_to_fruit  = s->value[N_TO_FRUIT];
	else                              n_to_fruit  = 0.;

	/*************************************************************************/

	/* reserve */
	/* note: special case reserve */
	s->value[N_TO_RESERVE] = s->value[N_LEAF_TO_RESERVE] +
			s->value[N_FROOT_TO_RESERVE] +
			s->value[N_CROOT_TO_RESERVE] +
			s->value[N_BRANCH_TO_RESERVE];

	/*************************************************************************/

	//note: if carbon transfer fluxes are positive than carbon and nitrogen to move are considered as "live tissues"
	//note: otherwise e.g. they need to be considered in their live and dead wood parts

	/* stem */
	if ( s->value[C_TO_STEM] > 0. )
	{
		s->value[N_TO_STEM] = s->value[C_TO_STEM] / s->value[CN_LIVEWOOD];

		n_to_stem  = s->value[N_TO_STEM];

	}
	else
	{
		s->value[N_TO_STEM] = (s->value[C_TO_STEM] * s->value[EFF_LIVE_TOTAL_WOOD_FRAC] / s->value[CN_LIVEWOOD]) +
				(s->value[C_TO_STEM] * ( 1. - s->value[EFF_LIVE_TOTAL_WOOD_FRAC]) / s->value[CN_DEADWOOD]);

		//		n_to_stem   = 0.;
	}

	/*************************************************************************/

	/* coarse root */
	if ( s->value[C_TO_CROOT] > 0. )
	{
		s->value[N_TO_CROOT] = s->value[C_TO_CROOT] / s->value[CN_LIVEWOOD];

		n_to_croot   = s->value[N_TO_CROOT];

	}
	else
	{
		s->value[N_TO_CROOT] = (s->value[C_TO_CROOT] * s->value[EFF_LIVE_TOTAL_WOOD_FRAC] / s->value[CN_LIVEWOOD]) +
				(s->value[C_TO_CROOT] * ( 1. - s->value[EFF_LIVE_TOTAL_WOOD_FRAC]) / s->value[CN_DEADWOOD]);

		//		n_to_croot   = 0.;
	}

	/*************************************************************************/

	/* branch */
	if ( s->value[C_TO_BRANCH] > 0.)
	{
		s->value[N_TO_BRANCH] = s->value[C_TO_BRANCH] / s->value[CN_LIVEWOOD];

		n_to_branch   = s->value[N_TO_BRANCH];

	}
	else
	{
		s->value[N_TO_BRANCH] = (s->value[C_TO_BRANCH] * s->value[EFF_LIVE_TOTAL_WOOD_FRAC] / s->value[CN_LIVEWOOD]) +
				(s->value[C_TO_BRANCH] * ( 1. - s->value[EFF_LIVE_TOTAL_WOOD_FRAC]) / s->value[CN_DEADWOOD]);

		n_to_branch   = 0.;
	}


	//TODO REMOVE ONCE MODIFIED
	//TODO 2 UPDATE NITROGEN POOLS

	/*** compute daily nitrogen demand ***/

	/* daily nitrogen demand */
	s->value[NPP_tN_DEMAND] = n_to_leaf + n_to_froot + n_to_stem + n_to_croot + n_to_branch + n_to_fruit /*+ n_to_reserve*/;

	/* tN/Cell/day -> gC/m2/day */
	s->value[NPP_gN_DEMAND] = s->value[NPP_tN_DEMAND] * 1e6 / g_settings->sizeCell;

	/* daily Nitrogen demand */
	s->value[TREE_N_DEMAND] = s->value[NPP_gN_DEMAND];

	//fixme
	if (s->value[NPP_gN_DEMAND] > c->soilN)
	{
		//todo back to partitioning-allocation routine and recompute both NPP in gC and NPP in gN based on the available soil nitrogen content
	}

	/* computing leaf sun and shaded Nitrogen pools (tN/ha) */
	if ( ! s->value[LEAF_C] )
	{
		s->value[LEAF_N]            = 0.;
	}
	else
	{
		s->value[LEAF_N]            = s->value[LEAF_C] / s->value[CN_LEAVES];
	}

	if ( ! s->value[FROOT_C] )
	{
		s->value[FROOT_N]           = 0.;
	}
	else
	{
		s->value[FROOT_N]           = s->value[FROOT_C]           / s->value[CN_FINE_ROOTS];
	}
	//s->value[LEAF_FALLING_N]  = s->value[LEAF_FALLING_C]        / s->value[CN_FALLING_LEAVES];
	s->value[STEM_LIVEWOOD_N]   = s->value[STEM_LIVEWOOD_C]       / s->value[CN_LIVEWOOD];
	s->value[STEM_DEADWOOD_N]   = s->value[STEM_DEADWOOD_C]       / s->value[CN_DEADWOOD];
	s->value[STEM_N]            = s->value[STEM_LIVEWOOD_N]       + s->value[STEM_DEADWOOD_N];
	s->value[CROOT_LIVEWOOD_N]  = s->value[CROOT_LIVEWOOD_C]      / s->value[CN_LIVEWOOD];
	s->value[CROOT_DEADWOOD_N]  = s->value[CROOT_DEADWOOD_C]      / s->value[CN_DEADWOOD];
	s->value[CROOT_N]           = s->value[CROOT_LIVEWOOD_N]      + s->value[CROOT_DEADWOOD_N];
	s->value[BRANCH_LIVEWOOD_N] = s->value[BRANCH_LIVEWOOD_C]     / s->value[CN_LIVEWOOD];
	s->value[BRANCH_DEADWOOD_N] = s->value[BRANCH_DEADWOOD_C]     / s->value[CN_DEADWOOD];
	s->value[BRANCH_N]          = s->value[BRANCH_LIVEWOOD_N]     + s->value[BRANCH_DEADWOOD_N];

	/* note: special case reserve */
	s->value[RESERVE_N]        += s->value[N_TO_RESERVE];

	/* check */
	CHECK_CONDITION ( s->value[LEAF_N],     < , ZERO );
	CHECK_CONDITION ( s->value[FROOT_N],    < , ZERO );
	CHECK_CONDITION ( s->value[STEM_N],     < , ZERO );
	CHECK_CONDITION ( s->value[BRANCH_N],   < , ZERO );
	CHECK_CONDITION ( s->value[CROOT_N],    < , ZERO );
	CHECK_CONDITION ( s->value[FRUIT_N],    < , ZERO );


}
