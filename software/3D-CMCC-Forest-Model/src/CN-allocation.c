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
	//s->value[C_TO_RESERVE] -= s->value[C_RESERVE_TO_CWD];
	s->value[C_TO_RESERVE] -= s->value[C_RESERVE_TO_LITR];

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

	/*** update dead carbon mass pools **/
	s->value[TOT_HEARTWOOD_C]          = (s->value[STEM_HEARTWOOD_C] + s->value[CROOT_HEARTWOOD_C] + s->value[BRANCH_HEARTWOOD_C] );
	s->value[TOT_DEADWOOD_C]           = (s->value[STEM_DEADWOOD_C]  + s->value[CROOT_DEADWOOD_C]  + s->value[BRANCH_DEADWOOD_C] );

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
		s->value[MAX_ANN_RESERVE_C]      = s->value[RESERVE_C];

		s->value[MAX_RESERVE_C_CONC] = ( s->value[MAX_ANN_RESERVE_C] / (s->value[TOT_SAPWOOD_C] * GC_GDM ) ) * 100.;
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

void carbon_allocation_new ( cell_t *const c, age_t *const a, species_t *const s, const int day, const int month, const int year )
{

	// 2024 new: 
	// it allocates daily assimilated carbon for both deciduous and evergreen daily
	// and removes the respired and dead parts of leaf,branch,fruits removed because of 
	// pruning and senescence.
	/* note: this mostly differentiates the 5.606 to 5.5 and 5p6 versions */
	// where the deadpool after mortality were removed from allocatable fluxes. 
	// however this created some numerical issues under some conditions.

    double RESERVE_C_CONC ;

	s->value[OLD_RESERVE_C] = s->value[RESERVE_C];

	logger(g_debug_log, "\n**CARBON ALLOCATION**\n");

	/*** removing growth respiration from carbon flux pools ***/
	s->value[C_TO_LEAF]    -= (s->value[LEAF_GROWTH_RESP]   / 1e6 * g_settings->sizeCell);
	s->value[C_TO_FROOT]   -= (s->value[FROOT_GROWTH_RESP]  / 1e6 * g_settings->sizeCell);
	s->value[C_TO_STEM]    -= (s->value[STEM_GROWTH_RESP]   / 1e6 * g_settings->sizeCell);
	s->value[C_TO_CROOT]   -= (s->value[CROOT_GROWTH_RESP]  / 1e6 * g_settings->sizeCell);
	s->value[C_TO_BRANCH]  -= (s->value[BRANCH_GROWTH_RESP] / 1e6 * g_settings->sizeCell);

	/***************************************************************************************/
	/***************************************************************************************/
	/* stem sapwood */
	s->value[STEM_SAPWOOD_C]          += ( s->value[C_TO_STEM]);
	/* live stem */
	s->value[STEM_LIVEWOOD_C]         += ( ( s->value[C_TO_STEM] * s->value[EFF_LIVE_TOTAL_WOOD_FRAC] ) ) ;
	/* coarse root sapwood */
	s->value[CROOT_SAPWOOD_C]         += ( s->value[C_TO_CROOT]- s->value[CROOT_C_TO_REMOVE]);
	/* live coarse root */
	s->value[CROOT_LIVEWOOD_C]        +=  (s->value[C_TO_CROOT]- s->value[CROOT_C_TO_REMOVE]) * s->value[EFF_LIVE_TOTAL_WOOD_FRAC]; 
	/* branch sapwood */
	s->value[BRANCH_SAPWOOD_C]        += ( s->value[C_TO_BRANCH] - s->value[BRANCH_C_TO_REMOVE] );
	/* live branch */
	s->value[BRANCH_LIVEWOOD_C]       += ( s->value[C_TO_BRANCH] - s->value[BRANCH_C_TO_REMOVE] ) * s->value[EFF_LIVE_TOTAL_WOOD_FRAC];

	/***************************************************************************************/
	/***************************************************************************************/

//   printf(" CN_allocation  species %s!!!\n", s->name);
   //printf(" s->value[BRANCH_C_TO_REMOVE]  %f!!!\n", s->value[BRANCH_C_TO_REMOVE]);
   //printf(" s->value[C_TO_BRANCH]  %f!!!\n",s->value[C_TO_BRANCH]);
   //printf(" s->value[BRANCH_C]  %f!!!\n",s->value[BRANCH_C]);
   //printf(" s->value[FROOT_C_TO_REMOVE]]  %f!!!\n",s->value[FROOT_C_TO_REMOVE]);
   //printf(" s->value[LEAF_C_TO_REMOVE]  %f!!!\n",s->value[LEAF_C_TO_REMOVE]);
   //printf(" s->value[FRUIT_C_TO_REMOVE]  %f!!!\n",s->value[FRUIT_C_TO_REMOVE]);

	/*** removing dead pools from carbon flux pools ***/
    // from pruning, senescence 

    s->value[C_TO_LEAF]    -= s->value[LEAF_C_TO_REMOVE]  ;
	s->value[C_TO_FROOT]   -= s->value[FROOT_C_TO_REMOVE] ;
	s->value[C_TO_CROOT]   -= s->value[CROOT_C_TO_REMOVE];
	s->value[C_TO_BRANCH]  -= s->value[BRANCH_C_TO_REMOVE];
	s->value[C_TO_FRUIT]   -= s->value[FRUIT_C_TO_REMOVE];

	/***************************************************************************************/
	/***************************************************************************************/

	/* for stocking wood */
	/* note: if here and NOT above model doesn't take into account dead parts */
	// net flux in to woody pool
	
	if ( s->value[C_TO_STEM]   > 0. ) s->value[YEARLY_C_TO_WOOD]     += s->value[C_TO_STEM];
	if ( s->value[C_TO_CROOT]  > 0. ) s->value[YEARLY_C_TO_WOOD]     += s->value[C_TO_BRANCH];
	if ( s->value[C_TO_BRANCH] > 0. ) s->value[YEARLY_C_TO_WOOD]     += s->value[C_TO_CROOT];
	if ( s->value[C_TO_STEM]   > 0. ) s->value[CUM_YEARLY_C_TO_WOOD] += s->value[C_TO_STEM];
	if ( s->value[C_TO_CROOT]  > 0. ) s->value[CUM_YEARLY_C_TO_WOOD] += s->value[C_TO_BRANCH];
	if ( s->value[C_TO_BRANCH] > 0. ) s->value[CUM_YEARLY_C_TO_WOOD] += s->value[C_TO_CROOT];

	/***************************************************************************************/
   #if 1
   //  printf(" ddalmo species %s!!!\n", s->name);
	// printf("in allocation  s->value[LEAF_C]              = %g tC/cell/day\n", (s->value[LEAF_C]* 1e6 / g_settings->sizeCell));
  //   printf("FROOT_C              = %g tC/cell/day\n", (s->value[FROOT_C]* 1e6 / g_settings->sizeCell));
	// printf("CROOT_C              = %g tC/cell/day\n", (s->value[CROOT_C]* 1e6 / g_settings->sizeCell));
	// printf("STEM_C               = %g tC/cell/day\n", (s->value[STEM_C]* 1e6 / g_settings->sizeCell));
	// printf("BRANCH_C             = %g tC/cell/day\n", (s->value[BRANCH_C]* 1e6 / g_settings->sizeCell));
	// printf("RESERVE_C            = %g tC/cell/day\n", (s->value[RESERVE_C]* 1e6 / g_settings->sizeCell));
  //  printf("FRUIT_C              = %g tC/cell/day\n", (s->value[FRUIT_C]* 1e6 / g_settings->sizeCell));

 
	// printf(" LEAF_C TO REMOVE                = %g tC/cell/day\n", (s->value[LEAF_C_TO_REMOVE]* 1e6 / g_settings->sizeCell));
    // printf("FROOT_C   TO REMOVE             = %g tC/cell/day\n", (s->value[FROOT_C_TO_REMOVE]* 1e6 / g_settings->sizeCell));
	// printf("BRANCH_C  TO REMOVE             = %g tC/cell/day\n", (s->value[BRANCH_C_TO_REMOVE]* 1e6 / g_settings->sizeCell));
	// printf("CROOT_C   TO REMOVE             = %g tC/cell/day\n", (s->value[CROOT_C_TO_REMOVE]* 1e6 / g_settings->sizeCell));
	// printf("FRUIT_C   TO REMOVE           = %g tC/cell/day\n", (s->value[FRUIT_C_TO_REMOVE]* 1e6 / g_settings->sizeCell));
	// printf("C_TO_RESERVE            = %g tC/cell/day\n",(s->value[C_TO_RESERVE]* 1e6 / g_settings->sizeCell));
   //printf("FRUIT_C              = %g tC/cell/day\n", s->value[FRUIT_C]);
#endif 

// * 1e6 / g_settings->sizeCell)

// // printf("in allocation s->value[C_TO_LEAF]             = %g tC/cell/day\n", (s->value[C_TO_LEAF]* 1e6 / g_settings->sizeCell));
// printf("in allocation s->value[C_TO_FROOT]             = %g tC/cell/day\n", (s->value[C_TO_FROOT]* 1e6 / g_settings->sizeCell));
// printf("in allocation s->value[C_TO_CROOT]             = %g tC/cell/day\n", (s->value[C_TO_CROOT]* 1e6 / g_settings->sizeCell));
// printf("in allocation s->value[C_TO_STEM]             = %g tC/cell/day\n", (s->value[C_TO_STEM]* 1e6 / g_settings->sizeCell));
// printf("in allocation s->value[C_TO_BRANCH]             = %g tC/cell/day\n", (s->value[C_TO_BRANCH]* 1e6 / g_settings->sizeCell));
// printf("in allocation s->value[C_TO_FRUIT]             = %g tC/cell/day\n", (s->value[C_TO_FRUIT]* 1e6 / g_settings->sizeCell));
// printf("in allocation s->value[C_TO_RESERVE]             = %g tC/cell/day\n", (s->value[C_TO_RESERVE]* 1e6 / g_settings->sizeCell));

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

	/*** update dead carbon mass pools **/
	s->value[TOT_DEADWOOD_C]           = (s->value[STEM_DEADWOOD_C]  + s->value[CROOT_DEADWOOD_C]  + s->value[BRANCH_DEADWOOD_C] );
	s->value[TOT_HEARTWOOD_C]          = (s->value[STEM_HEARTWOOD_C] + s->value[CROOT_HEARTWOOD_C] + s->value[BRANCH_HEARTWOOD_C] );

	/***************************************************************************************/

	/* check for closure */
	CHECK_CONDITION(fabs((s->value[STEM_SAPWOOD_C]   + s->value[STEM_HEARTWOOD_C])  - s->value[STEM_C]),  >,eps);
	CHECK_CONDITION(fabs((s->value[CROOT_SAPWOOD_C]  + s->value[CROOT_HEARTWOOD_C]) - s->value[CROOT_C]), >,eps);
	CHECK_CONDITION(fabs((s->value[BRANCH_SAPWOOD_C] + s->value[BRANCH_HEARTWOOD_C])- s->value[BRANCH_C]),>,eps);

	/* check */
	CHECK_CONDITION ( s->value[LEAF_C],     < , ZERO_C );
	CHECK_CONDITION ( s->value[FROOT_C],    < , ZERO_C );
	CHECK_CONDITION ( s->value[STEM_C],     < , ZERO_C );
	CHECK_CONDITION ( s->value[BRANCH_C],   < , ZERO_C );
	CHECK_CONDITION ( s->value[CROOT_C],    < , ZERO_C );
	CHECK_CONDITION ( s->value[FRUIT_C],    < , ZERO_C );
	
	// *********************************************************************
	
	// Note: this should not happen from the the 5.6
	// however if this happen, better to sort out why
#if 0
	if ( s->value[STEM_SAPWOOD_C] <  ZERO )
	{
		s->value[STEM_SAPWOOD_C]   = 0.;
		s->value[STEM_HEARTWOOD_C] = s->value[STEM_C];
		//printf("Warning: s->value[STEM_SAPWOOD_C] < ZERO \n");
	}
	if ( s->value[STEM_HEARTWOOD_C] < ZERO )
	{
		s->value[STEM_HEARTWOOD_C]   = 0.;
		s->value[STEM_SAPWOOD_C] = s->value[STEM_C];
		//printf("Warning: s->value[STEM_HEARTWOOD_C] < ZERO\n");
	}

	if ( s->value[CROOT_SAPWOOD_C] <  ZERO )
	{
		s->value[CROOT_SAPWOOD_C]   = 0.;
		s->value[CROOT_HEARTWOOD_C] = s->value[CROOT_C];
		//printf("Warning: s->value[CROOT_SAPWOOD_C] < ZERO\n");
	}
	if ( s->value[CROOT_HEARTWOOD_C] < ZERO )
	{	
		s->value[CROOT_HEARTWOOD_C]   = 0.;
		s->value[CROOT_SAPWOOD_C] = s->value[CROOT_C];
		//printf("Warning: s->value[CROOT_HEARTWOOD_C] < ZERO\n");

	}

	if ( s->value[BRANCH_SAPWOOD_C] <  ZERO )
	{
		s->value[BRANCH_SAPWOOD_C]   = 0.;
		s->value[BRANCH_HEARTWOOD_C] = s->value[BRANCH_C];
		//printf("Warning: s->value[BRANCH_SAPWOOD_C] < ZERO\n");
	}
	if ( s->value[BRANCH_HEARTWOOD_C] < ZERO )
	{
		s->value[BRANCH_HEARTWOOD_C]   = 0.;
		s->value[BRANCH_SAPWOOD_C] = s->value[BRANCH_C];
		//printf("Warning: s->value[BRANCH_HEARTWOOD_C] < ZERO\n" );
	}
#endif
	/*********************************************************************/
	// if the sapwood goes to 0, basically there is no Leaf, no GPP etc.
	// we put a more restrictive condition here
	
	CHECK_CONDITION ( s->value[STEM_SAPWOOD_C],     < , -ZERO );
	//CHECK_CONDITION ( s->value[STEM_HEARTWOOD_C],   < , -ZERO );
	CHECK_CONDITION ( s->value[CROOT_SAPWOOD_C],    < , -ZERO );
	//CHECK_CONDITION ( s->value[CROOT_HEARTWOOD_C],  < , -ZERO );
	CHECK_CONDITION ( s->value[BRANCH_SAPWOOD_C],   < , -ZERO );
	//CHECK_CONDITION ( s->value[BRANCH_HEARTWOOD_C], < , -ZERO );

	/* sub pools */
	CHECK_CONDITION ( s->value[STEM_SAPWOOD_C],     < , ZERO );
	CHECK_CONDITION ( s->value[STEM_HEARTWOOD_C],   < , ZERO );
	CHECK_CONDITION ( s->value[CROOT_SAPWOOD_C],    < , ZERO );
	CHECK_CONDITION ( s->value[CROOT_HEARTWOOD_C],  < , ZERO );
	CHECK_CONDITION ( s->value[BRANCH_SAPWOOD_C],   < , ZERO );
	CHECK_CONDITION ( s->value[BRANCH_HEARTWOOD_C], < , ZERO );

	/* sub pools */
	CHECK_CONDITION ( s->value[STEM_LIVEWOOD_C],    < , ZERO );
	CHECK_CONDITION ( s->value[STEM_DEADWOOD_C],    < , ZERO );
	CHECK_CONDITION ( s->value[CROOT_LIVEWOOD_C],   < , ZERO );
	CHECK_CONDITION ( s->value[CROOT_DEADWOOD_C],   < , ZERO );
	CHECK_CONDITION ( s->value[BRANCH_LIVEWOOD_C],  < , ZERO );
	CHECK_CONDITION ( s->value[BRANCH_DEADWOOD_C],  < , ZERO );

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
    // Compute Reserve-pool variables

    if ( s->value[RESERVE_C] > s->value[MAX_ANN_RESERVE_C])
	{
		s->value[MAX_ANN_RESERVE_C]               = s->value[RESERVE_C];
        s->value[TREE_MAX_ANN_RESERVE_C]          = s->value[RESERVE_C] / (double)s->counter[N_TREE];

	}

	if ( s->value[RESERVE_C] <= s->value[MIN_ANN_RESERVE_C])
	{
        s->value[MIN_ANN_RESERVE_C]               = s->value[RESERVE_C];
        s->value[TREE_MIN_ANN_RESERVE_C]          = s->value[RESERVE_C] / (double)s->counter[N_TREE]; 
	}

    // compute annual maximum and minimum of NSC concentration 
	// in relation to the actual ratio between reserve and sapwood available
       
    RESERVE_C_CONC = ( s->value[RESERVE_C] / (s->value[TOT_SAPWOOD_C] * GC_GDM ) ) * 100.;
        
    if ( RESERVE_C_CONC > s->value[MAX_RESERVE_C_CONC])
	{
        s->value[MAX_RESERVE_C_CONC] =  RESERVE_C_CONC ;
	}

	if ( RESERVE_C_CONC <= s->value[MIN_RESERVE_C_CONC])
	{
		s->value[MIN_RESERVE_C_CONC] = RESERVE_C_CONC ;
	}
		// 	printf(" CN ALLO  1 c->leaf_carbon   %f!!!\n",c->leaf_carbon  );
		// 	printf(" CN ALLO  1  c->croot_carbon   %f!!!\n",c->croot_carbon   );
		// 	printf(" CN ALLO  1 c->froot_carbon   %f!!!\n",c->froot_carbon   );
		// 	printf(" CN ALLO  1  c->stem_carbon   %f!!!\n",c->stem_carbon  );
		// 	printf(" CN ALLO  1  c->branch_carbon   %f!!!\n",c->branch_carbon  );
		// 	printf(" CN ALLO  1 c->reserve_carbon   %f!!!\n",c->reserve_carbon  );
		// 	printf(" CN ALLO  1  c->fruit_carbon   %f!!!\n",c->fruit_carbon  );

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

		// 		printf(" CN ALLO   c->leaf_carbon   %f!!!\n",c->leaf_carbon  );
		// 	printf(" CN ALLO    c->croot_carbon   %f!!!\n",c->croot_carbon   );
		// 	printf(" CN ALLO   c->froot_carbon   %f!!!\n",c->froot_carbon   );
		// 	printf(" CN ALLO    c->stem_carbon   %f!!!\n",c->stem_carbon  );
		// 	printf(" CN ALLO    c->branch_carbon   %f!!!\n",c->branch_carbon  );
		// 	printf(" CN ALLO   c->reserve_carbon   %f!!!\n",c->reserve_carbon  );
		// 	printf(" CN ALLO    c->fruit_carbon   %f!!!\n",c->fruit_carbon  );


        // Precision control for LEAF_C 
        // i.e. when LEAF_C is very low, it forced to 0
        // otherwise, if the LEAF_C value/ leaf_carbon is not exactly == 0, this leads the model to compute LAI_PROJ, and in cascade
        // all the leaf-related processes (e.g.radiation, water interception) are carried, even if the growing season is over. */

        if ( (s->value[LEAF_C]) < CRIT_PREC & (s->value[LEAF_C]) > 0. )  
         { 

			//printf(" SONO QUA PRECISION ISSUE!!!! \n");
            
            s->value[C_LEAF_TO_LITR]     += s->value[LEAF_C];
            s->value[TOTAL_C]             = s->value[TOTAL_C] - s->value[LEAF_C];
            
            // update pool
           //  c-> leaf_carbon               = c-> leaf_carbon - (s->value[LEAF_C] * 1e6 / g_settings->sizeCell);
            c-> leaf_carbon               = c-> leaf_carbon - (s->value[LEAF_C] * 1e6 / g_settings->sizeCell);
            s->value[LEAF_C]              = 0. ;
           
         } 

        if ( (s->value[FROOT_C]) < CRIT_PREC & (s->value[FROOT_C]) > 0. )  
         { 

            // TODO check if there are cases where this might lead to 
            // inconsistency when defoliation occurs 
            // i.e. if we need to force FROOT_C to 0  
            
               s->value[C_FROOT_TO_LITR]  += s->value[FROOT_C];
               s->value[TOTAL_C]           = s->value[TOTAL_C] - s->value[FROOT_C];
               c-> froot_carbon            = c-> froot_carbon - (s->value[FROOT_C]* 1e6 / g_settings->sizeCell);
               s->value[FROOT_C]           = 0. ;    
			   
         } 

        if ( (s->value[FRUIT_C]) < CRIT_PREC & (s->value[FRUIT_C]) > 0.  )  
         { 
			   s->value[C_FRUIT_TO_CWD]  += s->value[FRUIT_C];
			   s->value[C_TO_CWD]        += s->value[FRUIT_C];
               s->value[TOTAL_C]           -=  s->value[FRUIT_C];
               c-> fruit_carbon            -= (s->value[FRUIT_C]* 1e6 / g_settings->sizeCell);
               s->value[FRUIT_C]           = 0. ;     
         } 

	/* check */
	CHECK_CONDITION ( c->leaf_carbon,    < , ZERO_C );
	CHECK_CONDITION ( c->froot_carbon,   < , ZERO_C );
	CHECK_CONDITION ( c->stem_carbon,    < , ZERO_C );
	CHECK_CONDITION ( c->branch_carbon,  < , ZERO_C );
	CHECK_CONDITION ( c->croot_carbon,   < , ZERO_C );
	CHECK_CONDITION ( c->fruit_carbon,   < , ZERO_C );
	
	// 5p606
	// update tree level carbon pools
	
	
					/* compute tree average C pools */
					s->value[TREE_LEAF_C]                = (s->value[LEAF_C]             / (double)s->counter[N_TREE]);
					s->value[TREE_STEM_C]                = (s->value[STEM_C]             / (double)s->counter[N_TREE]);
					s->value[TREE_FROOT_C]               = (s->value[FROOT_C]            / (double)s->counter[N_TREE]);
					s->value[TREE_CROOT_C]               = (s->value[CROOT_C]            / (double)s->counter[N_TREE]);
					s->value[TREE_RESERVE_C]             = (s->value[RESERVE_C]          / (double)s->counter[N_TREE]);
					s->value[TREE_BRANCH_C]              = (s->value[BRANCH_C]           / (double)s->counter[N_TREE]);
					s->value[TREE_FRUIT_C]               = (s->value[FRUIT_C]            / (double)s->counter[N_TREE]);

					s->value[TREE_STEM_SAPWOOD_C]        = (s->value[STEM_SAPWOOD_C]     / (double)s->counter[N_TREE]);
					s->value[TREE_STEM_HEARTWOOD_C]      = (s->value[STEM_HEARTWOOD_C]   / (double)s->counter[N_TREE]);
					s->value[TREE_CROOT_SAPWOOD_C]       = (s->value[CROOT_SAPWOOD_C]    / (double)s->counter[N_TREE]);
					s->value[TREE_CROOT_HEARTWOOD_C]     = (s->value[CROOT_HEARTWOOD_C]  / (double)s->counter[N_TREE]);
					s->value[TREE_BRANCH_SAPWOOD_C]      = (s->value[BRANCH_SAPWOOD_C]   / (double)s->counter[N_TREE]);
					s->value[TREE_BRANCH_HEARTWOOD_C]    = (s->value[BRANCH_HEARTWOOD_C] / (double)s->counter[N_TREE]);
					s->value[TREE_SAPWOOD_C]             = (s->value[TOT_SAPWOOD_C]      / (double)s->counter[N_TREE]);
					s->value[TREE_HEARTWOOD_C]           = (s->value[TOT_HEARTWOOD_C]    / (double)s->counter[N_TREE]);

					s->value[TREE_STEM_LIVEWOOD_C]       = (s->value[STEM_LIVEWOOD_C]    / (double)s->counter[N_TREE]);
					s->value[TREE_STEM_DEADWOOD_C]       = (s->value[STEM_DEADWOOD_C]    / (double)s->counter[N_TREE]);
					s->value[TREE_CROOT_LIVEWOOD_C]      = (s->value[CROOT_LIVEWOOD_C]   / (double)s->counter[N_TREE]);
					s->value[TREE_CROOT_DEADWOOD_C]      = (s->value[CROOT_DEADWOOD_C]   / (double)s->counter[N_TREE]);
					s->value[TREE_BRANCH_LIVEWOOD_C]     = (s->value[BRANCH_LIVEWOOD_C]  / (double)s->counter[N_TREE]);
					s->value[TREE_BRANCH_DEADWOOD_C]     = (s->value[BRANCH_DEADWOOD_C]  / (double)s->counter[N_TREE]);
					s->value[TREE_TOT_LIVEWOOD_C]        = (s->value[TOT_LIVEWOOD_C]     / (double)s->counter[N_TREE]);
					s->value[TREE_TOT_DEADWOOD_C]        = (s->value[TOT_DEADWOOD_C]     / (double)s->counter[N_TREE]);

  // set to 0 
  s->value[LEAF_C_TO_REMOVE] = 0. ;
	 s->value[FROOT_C_TO_REMOVE]  = 0.;
	s->value[CROOT_C_TO_REMOVE] = 0.;
	 s->value[BRANCH_C_TO_REMOVE] = 0.;
	 s->value[FRUIT_C_TO_REMOVE] = 0.;

}
#if 0
// version for 5p5 and 5p6
void carbon_allocation_new_old ( cell_t *const c, age_t *const a, species_t *const s, const int day, const int month, const int year )
{

	/* it allocates daily assimilated carbon for both deciduous and evergreen daily
	 * and removes the respired and dead parts */
	/* note: this mostly differentiates the 5.4. to 5.5 (and subsequests) versions */

        double RESERVE_C_CONC ;

	s->value[OLD_RESERVE_C] = s->value[RESERVE_C];

	logger(g_debug_log, "\n**CARBON ALLOCATION**\n");

	/*** removing growth respiration from carbon flux pools ***/
	s->value[C_TO_LEAF]    -= (s->value[LEAF_GROWTH_RESP]   / 1e6 * g_settings->sizeCell);
	s->value[C_TO_FROOT]   -= (s->value[FROOT_GROWTH_RESP]  / 1e6 * g_settings->sizeCell);
	s->value[C_TO_STEM]    -= (s->value[STEM_GROWTH_RESP]   / 1e6 * g_settings->sizeCell);
	s->value[C_TO_CROOT]   -= (s->value[CROOT_GROWTH_RESP]  / 1e6 * g_settings->sizeCell);
	s->value[C_TO_BRANCH]  -= (s->value[BRANCH_GROWTH_RESP] / 1e6 * g_settings->sizeCell);

	/***************************************************************************************/
	/***************************************************************************************/


	/* stem sapwood */
	s->value[STEM_SAPWOOD_C]          += ( s->value[C_TO_STEM] - s->value[C_STEM_SAPWOOD_TO_CWD] );

	/* to avoid that self-thinning mortality happens to remove too much biomass */
	//if ( day && month )     // if (c->doy > 1 )  
        if (c->doy > 1 )
        {
         s->value[YEARLY_C_TO_STEM_SAPWOOD] += (s->value[C_TO_STEM]) - s->value[C_STEM_SAPWOOD_TO_CWD];
         } 
      
	/***************************************************************************************/

	/* live stem */
	s->value[STEM_LIVEWOOD_C]         += ( ( s->value[C_TO_STEM] * s->value[EFF_LIVE_TOTAL_WOOD_FRAC] ) - s->value[C_STEM_LIVEWOOD_TO_CWD] ) ;

	/* to avoid that self-thinning mortality happens to remove too much biomass */
 //        if ( day && month ) 
        if (c->doy > 1 )
        {
          s->value[YEARLY_C_TO_STEM_LIVEWOOD] += ( ( s->value[C_TO_STEM] * s->value[EFF_LIVE_TOTAL_WOOD_FRAC] ) - s->value[C_STEM_LIVEWOOD_TO_CWD]);
     
         }

	/***************************************************************************************/
	/***************************************************************************************/

	/* coarse root sapwood */
	s->value[CROOT_SAPWOOD_C]         += ( s->value[C_TO_CROOT] - s->value[C_CROOT_SAPWOOD_TO_CWD] );

	/* to avoid that self-thinning mortality happens to remove too much biomass */ 
        //if ( day && month ) 
        if (c->doy > 1 )
        {
         s->value[YEARLY_C_TO_CROOT_SAPWOOD] += (s->value[C_TO_CROOT] - s->value[C_CROOT_SAPWOOD_TO_CWD]);
        }


	/***************************************************************************************/

	/* live coarse root */
	s->value[CROOT_LIVEWOOD_C]        += ( ( s->value[C_TO_CROOT] * s->value[EFF_LIVE_TOTAL_WOOD_FRAC] ) - s->value[C_CROOT_LIVEWOOD_TO_CWD] );

	/* to avoid that self-thinning mortality happens to remove too much biomass */
        // if ( day && month ) 
        if (c->doy > 1 )
        {
         s->value[YEARLY_C_TO_CROOT_LIVEWOOD] += ( ( s->value[C_TO_CROOT] * s->value[EFF_LIVE_TOTAL_WOOD_FRAC] ) - s->value[C_CROOT_LIVEWOOD_TO_CWD] );
        }


	/***************************************************************************************/
	/***************************************************************************************/

	/* branch sapwood */
	s->value[BRANCH_SAPWOOD_C]        += ( s->value[C_TO_BRANCH] - s->value[C_BRANCH_SAPWOOD_TO_CWD] );

	/* to avoid that self-thinning mortality happens to remove too much biomass */
        //	if ( day && month )   
        if (c->doy > 1 )
        {
         s->value[YEARLY_C_TO_BRANCH_SAPWOOD]      += (s->value[C_TO_BRANCH] - s->value[C_BRANCH_SAPWOOD_TO_CWD]);     
        }
	/***************************************************************************************/

	/* live branch */
	s->value[BRANCH_LIVEWOOD_C]       += ( ( s->value[C_TO_BRANCH] * s->value[EFF_LIVE_TOTAL_WOOD_FRAC] ) - s->value[C_BRANCH_LIVEWOOD_TO_CWD] );

	/* to avoid that self-thinning mortality happens to remove too much biomass */
        //if ( day && month ) 
        if (c->doy > 1 )
        {
         s->value[YEARLY_C_TO_BRANCH_LIVEWOOD] += ( ( s->value[C_TO_BRANCH] * s->value[EFF_LIVE_TOTAL_WOOD_FRAC] ) - s->value[C_BRANCH_LIVEWOOD_TO_CWD] ); 
        }
	/***************************************************************************************/
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
	/***************************************************************************************/

	/* for stocking wood */
	/* note: if here and NOT above model doesn't take into account dead parts */
	// net flux in to woody pool
	
	if ( s->value[C_TO_STEM]   > 0. ) s->value[YEARLY_C_TO_WOOD]     += s->value[C_TO_STEM];
	if ( s->value[C_TO_CROOT]  > 0. ) s->value[YEARLY_C_TO_WOOD]     += s->value[C_TO_BRANCH];
	if ( s->value[C_TO_BRANCH] > 0. ) s->value[YEARLY_C_TO_WOOD]     += s->value[C_TO_CROOT];
	if ( s->value[C_TO_STEM]   > 0. ) s->value[CUM_YEARLY_C_TO_WOOD] += s->value[C_TO_STEM];
	if ( s->value[C_TO_CROOT]  > 0. ) s->value[CUM_YEARLY_C_TO_WOOD] += s->value[C_TO_BRANCH];
	if ( s->value[C_TO_BRANCH] > 0. ) s->value[CUM_YEARLY_C_TO_WOOD] += s->value[C_TO_CROOT];

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

	/*** update dead carbon mass pools **/
	s->value[TOT_DEADWOOD_C]           = (s->value[STEM_DEADWOOD_C]  + s->value[CROOT_DEADWOOD_C]  + s->value[BRANCH_DEADWOOD_C] );
	s->value[TOT_HEARTWOOD_C]          = (s->value[STEM_HEARTWOOD_C] + s->value[CROOT_HEARTWOOD_C] + s->value[BRANCH_HEARTWOOD_C] );

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
	
	// *********************************************************************
	// TODO check here on LEAF C
	
	// Precision control for LEAF_C 
        // i.e. when LEAF_C is very low, it forced to 0
        // otherwise, if the LEAF_C value/ leaf_carbon is not exactly == 0, this leads the model to compute LAI_PROJ, and in cascade
        // all the leaf-related processes (e.g.radiation, water interception) are carried, even if the growing season is over. */

         // if ( (s->value[LEAF_C]) < CRIT_PREC & (s->value[LEAF_C]) > 0.)    // in modo che saltiamo quando arriva in inverno che e' 0.
         //  { 
        
         ////   s->value[C_LEAF_TO_LITR]     += s->value[LEAF_C];
            //s->value[TOTAL_C]             = s->value[TOTAL_C] - s->value[LEAF_C];
            
            // update pool
            
           // c-> leaf_carbon               = c-> leaf_carbon - (s->value[LEAF_C] * 1e6 / g_settings->sizeCell);
          ////  s->value[LEAF_C]              = 0. ;
            
            // ned to fix C_TO_LEAF
           
            // TODO check if there are cases where this might lead to 
            // inconsistency when defoliation occurs 
            // i.e. if we need to force FROOT_C to 0  
            
            //   s->value[C_FROOT_TO_LITR]  += s->value[FROOT_C];
            //   s->value[TOTAL_C]           = s->value[TOTAL_C] - s->value[FROOT_C];
            //   c-> froot_carbon            = c-> froot_carbon - (s->value[FROOT_C]* 1e6 / g_settings->sizeCell);
            //   s->value[FROOT_C]           = 0. ;     
          // } 
	

	/*********************************************************************/

	//Note: this should not happen from the the 5.6
	// however if this happen, better to sort out why
	  
	if ( s->value[STEM_SAPWOOD_C] <  ZERO )
	{
		s->value[STEM_SAPWOOD_C]   = 0.;
		s->value[STEM_HEARTWOOD_C] = s->value[STEM_C];
		printf("Warning: s->value[STEM_SAPWOOD_C] < ZERO\n");
	}
	if ( s->value[STEM_HEARTWOOD_C] < ZERO )
	{
		s->value[STEM_HEARTWOOD_C]   = 0.;
		s->value[STEM_SAPWOOD_C] = s->value[STEM_C];
		printf("Warning: s->value[STEM_HEARTWOOD_C] < ZERO\n");
	}

	if ( s->value[CROOT_SAPWOOD_C] <  ZERO )
	{
		s->value[CROOT_SAPWOOD_C]   = 0.;
		s->value[CROOT_HEARTWOOD_C] = s->value[CROOT_C];
		printf("Warning: s->value[CROOT_SAPWOOD_C] < ZERO\n");
	}
	if ( s->value[CROOT_HEARTWOOD_C] < ZERO )
	{	
		s->value[CROOT_HEARTWOOD_C]   = 0.;
		s->value[CROOT_SAPWOOD_C] = s->value[CROOT_C];
		printf("Warning: s->value[CROOT_HEARTWOOD_C] < ZERO\n");

	}

	if ( s->value[BRANCH_SAPWOOD_C] <  ZERO )
	{
		s->value[BRANCH_SAPWOOD_C]   = 0.;
		s->value[BRANCH_HEARTWOOD_C] = s->value[BRANCH_C];
		printf("Warning: s->value[BRANCH_SAPWOOD_C] < ZERO\n");
	}
	if ( s->value[BRANCH_HEARTWOOD_C] < ZERO )
	{
		s->value[BRANCH_HEARTWOOD_C]   = 0.;
		s->value[BRANCH_SAPWOOD_C] = s->value[BRANCH_C];
		printf("Warning: s->value[BRANCH_HEARTWOOD_C] < ZERO\n");
	}

	/*********************************************************************/

	/* sub pools */
	CHECK_CONDITION ( s->value[STEM_SAPWOOD_C],     < , ZERO );
	CHECK_CONDITION ( s->value[STEM_HEARTWOOD_C],   < , ZERO );
	CHECK_CONDITION ( s->value[CROOT_SAPWOOD_C],    < , ZERO );
	CHECK_CONDITION ( s->value[CROOT_HEARTWOOD_C],  < , ZERO );
	CHECK_CONDITION ( s->value[BRANCH_SAPWOOD_C],   < , ZERO );
	CHECK_CONDITION ( s->value[BRANCH_HEARTWOOD_C], < , ZERO );

	/* sub pools */
	CHECK_CONDITION ( s->value[STEM_LIVEWOOD_C],    < , ZERO );
	CHECK_CONDITION ( s->value[STEM_DEADWOOD_C],    < , ZERO );
	CHECK_CONDITION ( s->value[CROOT_LIVEWOOD_C],   < , ZERO );
	CHECK_CONDITION ( s->value[CROOT_DEADWOOD_C],   < , ZERO );
	CHECK_CONDITION ( s->value[BRANCH_LIVEWOOD_C],  < , ZERO );
	CHECK_CONDITION ( s->value[BRANCH_DEADWOOD_C],  < , ZERO );

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

	

        // Compute Reserve-pool variables

        if ( s->value[RESERVE_C] > s->value[MAX_ANN_RESERVE_C])
	{
		s->value[MAX_ANN_RESERVE_C]               = s->value[RESERVE_C];
                s->value[TREE_MAX_ANN_RESERVE_C]          = s->value[RESERVE_C] / (double)s->counter[N_TREE];

	}

	if ( s->value[RESERVE_C] <= s->value[MIN_ANN_RESERVE_C])
	{
                s->value[MIN_ANN_RESERVE_C]               = s->value[RESERVE_C];
                s->value[TREE_MIN_ANN_RESERVE_C]          = s->value[RESERVE_C] / (double)s->counter[N_TREE]; 
	}

        // compute annual maximum and minimum of NSC concentration in relation to the actual ratio between reserve and sapwood available
       
        RESERVE_C_CONC = ( s->value[RESERVE_C] / (s->value[TOT_SAPWOOD_C] * GC_GDM ) ) * 100.;
        
        if ( RESERVE_C_CONC > s->value[MAX_RESERVE_C_CONC])
	{
        	s->value[MAX_RESERVE_C_CONC] =  RESERVE_C_CONC ;
	}

	if ( RESERVE_C_CONC <= s->value[MIN_RESERVE_C_CONC])
	{
		s->value[MIN_RESERVE_C_CONC] = RESERVE_C_CONC ;
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

        // Precision control for LEAF_C 
        // i.e. when LEAF_C is very low, it forced to 0
        // otherwise, if the LEAF_C value/ leaf_carbon is not exactly == 0, this leads the model to compute LAI_PROJ, and in cascade
        // all the leaf-related processes (e.g.radiation, water interception) are carried, even if the growing season is over. */

        if ( (s->value[LEAF_C]) < CRIT_PREC)  
         { 
            
            s->value[C_LEAF_TO_LITR]     += s->value[LEAF_C];
            s->value[TOTAL_C]             = s->value[TOTAL_C] - s->value[LEAF_C];
            
            // update pool
           //  c-> leaf_carbon               = c-> leaf_carbon - (s->value[LEAF_C] * 1e6 / g_settings->sizeCell);
            c-> leaf_carbon               = c-> leaf_carbon - (s->value[LEAF_C] * 1e6 / g_settings->sizeCell);
            s->value[LEAF_C]              = 0. ;
           
            // TODO check if there are cases where this might lead to 
            // inconsistency when defoliation occurs 
            // i.e. if we need to force FROOT_C to 0  
            
            //   s->value[C_FROOT_TO_LITR]  += s->value[FROOT_C];
            //   s->value[TOTAL_C]           = s->value[TOTAL_C] - s->value[FROOT_C];
            //   c-> froot_carbon            = c-> froot_carbon - (s->value[FROOT_C]* 1e6 / g_settings->sizeCell);
            //   s->value[FROOT_C]           = 0. ;     
         } 

	/* check */
	CHECK_CONDITION ( c->leaf_carbon,    < , ZERO );
	CHECK_CONDITION ( c->froot_carbon,   < , ZERO );
	CHECK_CONDITION ( c->stem_carbon,    < , ZERO );
	CHECK_CONDITION ( c->branch_carbon,  < , ZERO );
	CHECK_CONDITION ( c->croot_carbon,   < , ZERO );
	CHECK_CONDITION ( c->fruit_carbon,   < , ZERO );
	
	
	#if 1
	
	// 5p606
	// update tree level carbon pools
	
	
					/* compute tree average C pools */
					s->value[TREE_LEAF_C]                = (s->value[LEAF_C]             / (double)s->counter[N_TREE]);
					s->value[TREE_STEM_C]                = (s->value[STEM_C]             / (double)s->counter[N_TREE]);
					s->value[TREE_FROOT_C]               = (s->value[FROOT_C]            / (double)s->counter[N_TREE]);
					s->value[TREE_CROOT_C]               = (s->value[CROOT_C]            / (double)s->counter[N_TREE]);
					s->value[TREE_RESERVE_C]             = (s->value[RESERVE_C]          / (double)s->counter[N_TREE]);
					s->value[TREE_BRANCH_C]              = (s->value[BRANCH_C]           / (double)s->counter[N_TREE]);
					s->value[TREE_FRUIT_C]               = (s->value[FRUIT_C]            / (double)s->counter[N_TREE]);

					s->value[TREE_STEM_SAPWOOD_C]        = (s->value[STEM_SAPWOOD_C]     / (double)s->counter[N_TREE]);
					s->value[TREE_STEM_HEARTWOOD_C]      = (s->value[STEM_HEARTWOOD_C]   / (double)s->counter[N_TREE]);
					s->value[TREE_CROOT_SAPWOOD_C]       = (s->value[CROOT_SAPWOOD_C]    / (double)s->counter[N_TREE]);
					s->value[TREE_CROOT_HEARTWOOD_C]     = (s->value[CROOT_HEARTWOOD_C]  / (double)s->counter[N_TREE]);
					s->value[TREE_BRANCH_SAPWOOD_C]      = (s->value[BRANCH_SAPWOOD_C]   / (double)s->counter[N_TREE]);
					s->value[TREE_BRANCH_HEARTWOOD_C]    = (s->value[BRANCH_HEARTWOOD_C] / (double)s->counter[N_TREE]);
					s->value[TREE_SAPWOOD_C]             = (s->value[TOT_SAPWOOD_C]      / (double)s->counter[N_TREE]);
					s->value[TREE_HEARTWOOD_C]           = (s->value[TOT_HEARTWOOD_C]    / (double)s->counter[N_TREE]);

					s->value[TREE_STEM_LIVEWOOD_C]       = (s->value[STEM_LIVEWOOD_C]    / (double)s->counter[N_TREE]);
					s->value[TREE_STEM_DEADWOOD_C]       = (s->value[STEM_DEADWOOD_C]    / (double)s->counter[N_TREE]);
					s->value[TREE_CROOT_LIVEWOOD_C]      = (s->value[CROOT_LIVEWOOD_C]   / (double)s->counter[N_TREE]);
					s->value[TREE_CROOT_DEADWOOD_C]      = (s->value[CROOT_DEADWOOD_C]   / (double)s->counter[N_TREE]);
					s->value[TREE_BRANCH_LIVEWOOD_C]     = (s->value[BRANCH_LIVEWOOD_C]  / (double)s->counter[N_TREE]);
					s->value[TREE_BRANCH_DEADWOOD_C]     = (s->value[BRANCH_DEADWOOD_C]  / (double)s->counter[N_TREE]);
					s->value[TREE_TOT_LIVEWOOD_C]        = (s->value[TOT_LIVEWOOD_C]     / (double)s->counter[N_TREE]);
					s->value[TREE_TOT_DEADWOOD_C]        = (s->value[TOT_DEADWOOD_C]     / (double)s->counter[N_TREE]);

 
    #endif

}
#endif

/******************************************************************************************************************************/

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
	//additional condition, to reconcilie when LEAF_C is forced to be 0
	
	if ( ! s->value[LEAF_C] )
	{  
                s->value[N_LEAF_TO_LITR]     += s->value[LEAF_N];
   
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
