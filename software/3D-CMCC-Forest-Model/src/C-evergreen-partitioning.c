/*
 * C-evergreen-partitioning-allocation.c
 *
 *  Created on: 20/nov/2013
 *      Author: alessio
 */
/* includes */
#include "CN-allocation.h" 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "C-evergreen-partitioning.h"
#include "common.h"
#include "constants.h"
#include "settings.h"
#include "logger.h"
#include "lai.h"
#include "leaf_fall.h"
#include "turnover.h"
#include "dendometry.h"
#include "biomass.h"
#include "check_balance.h"
#include "mortality.h"
#include "aut_respiration.h"

extern logger_t* g_debug_log;
extern settings_t* g_settings;

#define TEST 0

/* Evergreen carbon allocation routine */
void daily_C_evergreen_partitioning (cell_t *const c, const int layer, const int height, const int dbh, const int age, const int species,
		const meteo_daily_t *const meteo_daily, const int day, const int month, const int year)
{
	double s0;
	double r0;
	double omega;
	double pS;
	double pR;
	double pL;
	double Light_trasm;
	double npp_to_alloc;

	double delta_leaf    = 0.;
	double delta_froot   = 0.;
	double delta_croot   = 0.;
	double delta_branch  = 0.;
	double delta_stem    = 0.;
	double delta_fruit   = 0.;
	double delta_reserve_alloc  = 0.;
	double delta_reserve_deple  = 0.;

	height_t *h;
	dbh_t *d;
	age_t *a;
	species_t *s;

	h = &c->heights[height];
	d = &h->dbhs[dbh];
	a = &d->ages[age];
	s = &a->species[species];

	s0    = s->value[S0CTEM];        /* parameter controlling allocation to stem (minimum ratio to stem pool) */
	r0    = s->value[R0CTEM];        /* parameter controlling allocation to roots (minimum ratio to root pools) */
	omega = s->value[OMEGA_CTEM];    /* controls the sensitivity of allocation to changes in water and light availability */

	//fixme it should takes into account above layers
	Light_trasm = exp( - s->value[K] * s->value[LAI_PROJ] );

	/* Marconi: here the allocation of biomass reserve is divided in fine root and leaves following the
	 * allocation ratio parameter between them. That because
	 * in evergreen we don't have bud burst phenology phase, and indeed there are two phenology phases;
	 * the former in which carbon is allocated in fine root and foliage, the latter in
	 * every pool except foliage and fine root */

	logger(g_debug_log, "\n**C-PARTITIONING**\n");
	logger(g_debug_log, "Carbon partitioning for evergreen\n");

	/* partitioning block using approach of:
	 *  Potter et al., 1993, Friendlingstein et al., 1999;
	 *  Schwalm & Ek, 2004; Arora and Boer 2005 */
	logger(g_debug_log, "*Partitioning ratios*\n");

	/* roots */
	pR = ( r0 + ( omega * ( 1. - s->value[F_SW] ) ) ) / ( 1. + ( omega * ( 2. - Light_trasm - s->value[F_SW] ) ) );
	//logger(g_debug_log, "Roots CTEM ratio = %g %%\n", pR * 100.);

	/* stem */
	pS = ( s0 + ( omega * ( 1. - Light_trasm ) ) ) / ( 1. + ( omega * ( 2. - Light_trasm - s->value[F_SW] ) ) );
	//logger(g_debug_log, "Stem CTEM ratio = %g %%\n", pS * 100.);

	/* reserve and leaves */
	pL = ( 1. - pS - pR );
	logger(g_debug_log, "Reserve CTEM ratio = %g %%\n", pL * 100.);
	CHECK_CONDITION( fabs ( pR + pS + pL ), >, 1 + eps );


	logger(g_debug_log, "PHENOLOGICAL PHASE = %d\n", s->phenology_phase);
	logger(g_debug_log, "LAI_PROJ = %g\n", s->value[LAI_PROJ]);
	logger(g_debug_log, "PEAK_LAI_PROJ = %g \n", s->value[PEAK_LAI_PROJ]);

	if ( g_settings->Prog_Aut_Resp )
	{
		/* assign NPP to local variables and remove the maintenance respiration */
		npp_to_alloc = s->value[GPP_tC] - s->value[TOTAL_MAINT_RESP_tC] ;
	}
	else
	{
		/* assign NPP to local variables and remove the prognostic respiration */
		npp_to_alloc = s->value[GPP_tC] * ( 1. - g_settings->Fixed_Aut_Resp_rate );
	}
	logger(g_debug_log, "s->value[GPP_tC] = %g\n", s->value[GPP_tC]);
	logger(g_debug_log, "s->value[TOTAL_MAINT_RESP_tC]  = %g\n", s->value[TOTAL_MAINT_RESP_tC]);
	logger(g_debug_log, "npp_to_alloc = %g tC/sizecell/day\n", npp_to_alloc);

       

	/* note: none carbon pool is refilled if reserve is lower than minimum */
	/* reserves have priority before all other pools!!! */
	/* "...species-specific minimum pool (of NSC) sizes of stored C may be necessary to avoid catastrophic xylem failure" 
	   Sala et al., 2012, Tree Physiology */ 


#if 0
	switch ( s->phenology_phase )
	{
	/************************************************************/
	case 1:
		/*
		 * this phenological phase happens when:
		 * - month < 6 for north hemisphere
		 * - LAI < PEAK_LAI
		 * */

		/*just a fraction of biomass reserve is used for foliage the other part is allocated to the stem (Magnani pers comm),
		 * and Barbaroux et al., 2002, the ratio is driven by the BIOME_BGC newStem:newLeaf ratio
		 * the fraction of reserve to allocate for foliage is re-computed for each of the BUD_BURST days
		 * sharing the daily remaining amount (taking into account respiration costs)of NSC */

		/*note: every day of the year leaf and fine root turnover happens, that's why every day
		 * an amount of leaf and fine root C is re-translocated to the reserve pool*/

		/* partitioning */
		if (npp_to_alloc > 0.0)
		{
			/* check if minimum reserve pool needs to be refilled */
			/* it doesn't need */
			if(s->value[RESERVE_C] >= s->value[MIN_RESERVE_C])
			{
				double max_leafC        = 0.;
				double exceeding_leafC  = 0.;
				double max_frootC       = 0.;
				double exceeding_frootC = 0.;
				double exceeding_C      = 0.;
				double pR1              = 0.;
				double pS1              = 0.;


				logger(g_debug_log, "Allocating only into leaf and fine root pools (positive NPP)\n");

				/* update carbon flux */
				s->value[C_TO_LEAF]    = npp_to_alloc * (1. - s->value[FROOT_LEAF_FRAC]);
				s->value[C_TO_FROOT]   = npp_to_alloc * s->value[FROOT_LEAF_FRAC];

				/********************************************************************************************************************************************/
				/* check for leaf C > max leaf C */
				if ( ( ( s->value[C_TO_LEAF] * ( 1. - s->value[EFF_GRPERC] ) ) + s->value[LEAF_C]) > s->value[MAX_LEAF_C])
				{
					logger(g_debug_log, "..exceeding leaf carbon to other pools\n");

					max_leafC       = s->value[MAX_LEAF_C] - s->value[LEAF_C];
					exceeding_leafC = ( ( npp_to_alloc * ( 1. - s->value[FROOT_LEAF_FRAC] )) - ( max_leafC + ( max_leafC * s->value[EFF_GRPERC] ) ) );

					s->value[C_TO_LEAF] = max_leafC + ( max_leafC * s->value[EFF_GRPERC] );

				}

				/* check for fine root C > max fine root C */
				if ( ( (s->value[C_TO_FROOT] * ( 1. - s->value[EFF_GRPERC] ) ) + s->value[FROOT_C] ) > s->value[MAX_FROOT_C])
				{
					logger(g_debug_log, "..exceeding fine root carbon to other pools\n");

					max_frootC       = s->value[MAX_FROOT_C] - s->value[FROOT_C];
					exceeding_frootC = ( ( npp_to_alloc * s->value[FROOT_LEAF_FRAC] ) - ( max_frootC + ( max_frootC * s->value[EFF_GRPERC] ) ) );

					s->value[C_TO_FROOT] = max_frootC + ( max_frootC * s->value[EFF_GRPERC] );
				}

				/* exceeding carbon */
				exceeding_C = exceeding_leafC + exceeding_frootC;

				/* compute ratios */
				pR1 = pR / ( pR + pS );
				pS1 = 1. - pR1;
				CHECK_CONDITION( fabs ( pR1 + pS1 ), >, 1 + eps );

				/* to other carbon pools */
				s->value[C_TO_CROOT]     = (exceeding_C * pR1);
				s->value[C_TO_STEM]      = (exceeding_C * pS1) * (1. - s->value[FRACBB]);
				s->value[C_TO_BRANCH]    = (exceeding_C * pS1) * s->value[FRACBB];

				/********************************************************************************************************************************************/
			}
			/* it needs */
			else
			{
				logger(g_debug_log, "Allocating only into reserve pool (low reserves, positive NPP)\n");
				s->value[C_TO_RESERVE] = npp_to_alloc;
			}
		}
		else
		{
			logger(g_debug_log, "Consuming reserve pool (negative NPP)\n");
			s->value[C_TO_RESERVE]     = npp_to_alloc;
		}
		/**********************************************************************/
		break;
	case 2:
		/*
		 * this phenological phase happens when:
		 * - month > 6 for north hemisphere
		 * - LAI > PEAK_LAI
		 */

		/* partitioning */
		if (npp_to_alloc > 0.0)
		{

			/* check if minimum reserve pool needs to be refilled */
			/* it doesn't need */
			if(s->value[RESERVE_C] >= s->value[MIN_RESERVE_C])
			{
				logger(g_debug_log, "Allocating only into Coarse root, Reserve, Stem and Branch pools (positive NPP)\n");

				/* reproduction */
				if ( g_settings->regeneration && ( a->value > s->value[SEXAGE] ) )
				{
					logger(g_debug_log, "allocating into fruit pool\n");

					s->value[C_TO_FRUIT] = npp_to_alloc * s->value[FRUIT_PERC];
					npp_to_alloc        -= s->value[C_TO_FRUIT];
				}

				s->value[C_TO_CROOT]     = (npp_to_alloc * pR);
				s->value[C_TO_RESERVE]   = (npp_to_alloc * pL);
				s->value[C_TO_STEM]      = (npp_to_alloc * pS) * ( 1. - s->value[FRACBB] );
				s->value[C_TO_BRANCH]    = (npp_to_alloc * pS) * s->value[FRACBB];


			}
			/* it needs */
			else
			{
				logger(g_debug_log, "Allocating only into reserve pool (low reserves, positive NPP)\n");
				s->value[C_TO_RESERVE]   = npp_to_alloc;
			}
		}
		else
		{
			logger(g_debug_log, "Consuming reserve pool (negative NPP)\n");
			s->value[C_TO_RESERVE]       = npp_to_alloc;
		}
		/**********************************************************************/
		break;
	}
#else

	//NEW function TEST 31 MARCH 2017
	//note: this is the simplest function for evergreen species
	//todo in case it will be definitely used then remove evergreen phenology from phenology function

	/*just a fraction of biomass reserve is used for foliage the other part is allocated to the stem (Magnani pers comm),
	 * and Barbaroux et al., 2002, the ratio is driven by the BIOME_BGC newStem:newLeaf ratio
	 * the fraction of reserve to allocate for foliage is re-computed for each of the BUD_BURST days
	 * sharing the daily remaining amount (taking into account respiration costs) of NSC */

	/*note: every day of the year leaf and fine root turnover happen, that's why every day
	 * an amount of leaf and fine root C is re-translocated to the reserve pool*/
	/*note: for patterns and NSC usage behavior in conifer see also Woodruff & Meinzer, 2011, Plant Cell and Environment */

	/* partitioning */
	if ( npp_to_alloc > 0. )
	{
		/* check if minimum reserve pool needs to be refilled */
		/* it doesn't need */
		if( s->value[RESERVE_C] >= s->value[MIN_RESERVE_C] )
		{
			double max_leafC        = 0.;
			double exceeding_leafC  = 0.;
			double max_frootC       = 0.;
			double exceeding_frootC = 0.;
			double exceeding_C      = 0.;
			double pR1              = 0.;
			double pS1              = 0.;
                     
			if ( month < JULY )
			// test jan 24, for boreal sites OPJ on the verision 5p6 clicfor
			//if ( month < SEPTEMBER ) 
			{

				/* if management doens't happen (this to avoid problems in carbon balance) */
				if ( ! s->counter[THINNING_HAPPENS] )
				{
                                    
					/* update carbon flux */
					s->value[C_TO_LEAF]      = npp_to_alloc * ( 1. - s->value[FROOT_LEAF_FRAC] );
					s->value[C_TO_FROOT]     = npp_to_alloc * s->value[FROOT_LEAF_FRAC];
                                                 
					/********************************************************************************************************************************************/
					/* check for leaf C > max leaf C */
					if ( ( ( s->value[C_TO_LEAF] * ( 1. - s->value[EFF_GRPERC] ) ) + s->value[LEAF_C] ) > s->value[MAX_LEAF_C] )
					{
						logger(g_debug_log, "..exceeding leaf carbon to other pools\n");

						max_leafC           = s->value[MAX_LEAF_C] - s->value[LEAF_C];

						exceeding_leafC     = ( ( npp_to_alloc * ( 1. - s->value[FROOT_LEAF_FRAC] ) ) - ( max_leafC + ( max_leafC * s->value[EFF_GRPERC] ) ) );

						s->value[C_TO_LEAF] = max_leafC + ( max_leafC * s->value[EFF_GRPERC] );
					}

					/* check for fine root C > max fine root C */
					if ( ( ( s->value[C_TO_FROOT] * ( 1. - s->value[EFF_GRPERC] ) ) + s->value[FROOT_C] ) > s->value[MAX_FROOT_C] )
					{
						logger(g_debug_log, "..exceeding fine root carbon to other pools\n");

						max_frootC           = s->value[MAX_FROOT_C] - s->value[FROOT_C];

						exceeding_frootC     = ( ( npp_to_alloc * s->value[FROOT_LEAF_FRAC] ) - ( max_frootC + ( max_frootC * s->value[EFF_GRPERC] ) ) );

						s->value[C_TO_FROOT] = max_frootC + ( max_frootC * s->value[EFF_GRPERC] );
                                         
					}

					/* exceeding carbon */
					exceeding_C = exceeding_leafC + exceeding_frootC;

					/* compute ratios */
					pR1 = pR / ( pR + pS );
					pS1 = 1. - pR1;
					CHECK_CONDITION( fabs ( pR1 + pS1 ), >, 1 + eps );

					/* to other carbon pools */
					s->value[C_TO_CROOT]     = (exceeding_C * pR1);
					s->value[C_TO_STEM]      = (exceeding_C * pS1) * (1. - s->value[FRACBB]);
					s->value[C_TO_BRANCH]    = (exceeding_C * pS1) * s->value[FRACBB];
                                      

				}
				else
				{
					s->value[C_TO_RESERVE] = npp_to_alloc;

				}
			}
			else
			{
				/* as in Bossel (1996), second priority is assigned to fruit */
				/* reproduction */
				if ( a->value > s->value[SEXAGE] )
				{
					logger(g_debug_log, "allocating into fruit pool\n");

					s->value[C_TO_FRUIT] = npp_to_alloc * s->value[FRUIT_PERC];
					npp_to_alloc        -= s->value[C_TO_FRUIT];
				}

				s->value[C_TO_CROOT]     = (npp_to_alloc * pR);
				s->value[C_TO_LEAF]      = (npp_to_alloc * pL) * ( 1. - s->value[FROOT_LEAF_FRAC] );
				s->value[C_TO_FROOT]     = (npp_to_alloc * pL) * s->value[FROOT_LEAF_FRAC];
				s->value[C_TO_STEM]      = (npp_to_alloc * pS) * ( 1. - s->value[FRACBB] );
				s->value[C_TO_BRANCH]    = (npp_to_alloc * pS) * s->value[FRACBB];

				/********************************************************************************************************************************************/
				/* check for leaf C > max leaf C */
				if ( ( ( s->value[C_TO_LEAF] * ( 1. - s->value[EFF_GRPERC] ) ) + s->value[LEAF_C] ) > s->value[MAX_LEAF_C] )
				{
					logger(g_debug_log, "..exceeding leaf carbon to other pools\n");

					max_leafC       = s->value[MAX_LEAF_C] - s->value[LEAF_C];
					exceeding_leafC = ( ( ( npp_to_alloc * pL ) * ( 1. - s->value[FROOT_LEAF_FRAC] ) ) - ( max_leafC + ( max_leafC * s->value[EFF_GRPERC] ) ) );

					s->value[C_TO_LEAF] = max_leafC + ( max_leafC * s->value[EFF_GRPERC] );
				}

				/* check for fine root C > max fine root C */
				if ( ( (s->value[C_TO_FROOT] * ( 1. - s->value[EFF_GRPERC] ) ) + s->value[FROOT_C] ) > s->value[MAX_FROOT_C])
				{
					logger(g_debug_log, "..exceeding fine root carbon to other pools\n");

					max_frootC       = s->value[MAX_FROOT_C] - s->value[FROOT_C];
					exceeding_frootC = ( ( ( npp_to_alloc * pL ) * s->value[FROOT_LEAF_FRAC] )  - ( max_frootC + ( max_frootC * s->value[EFF_GRPERC] ) ) );

					s->value[C_TO_FROOT] = max_frootC + ( max_frootC * s->value[EFF_GRPERC] );
				}

				/* exceeding carbon */
				exceeding_C = exceeding_leafC + exceeding_frootC;

				/* compute ratios */
				pR1 = pR / ( pR + pS );
				pS1 = 1. - pR1;
				CHECK_CONDITION( fabs ( pR1 + pS1 ), >, 1 + eps );

				/* to other carbon pools */
				s->value[C_TO_CROOT]    += (exceeding_C * pR1);
				s->value[C_TO_STEM]     += (exceeding_C * pS1) * ( 1. - s->value[FRACBB] );
				s->value[C_TO_BRANCH]   += (exceeding_C * pS1) * s->value[FRACBB];

			}

			/********************************************************************************************************************************************/
		}
		/* it needs */
		else
		{
			logger(g_debug_log, "Allocating only into reserve pool (low reserves, positive NPP)\n");
			s->value[C_TO_RESERVE] = npp_to_alloc;
		}
	}
	else
	{

		logger(g_debug_log, "Consuming reserve pool (negative NPP)\n");

		/* consuming reserve carbon pools */

		s->value[C_TO_RESERVE]     = npp_to_alloc;


#if 0 //TODO to be tested
		if ( s->value[RESERVE_C] >= s->value[MIN_RESERVE_C] )
		{
			logger(g_debug_log, "Consuming reserve pool (negative NPP)\n");

			/* consuming reserve carbon pools */

			s->value[C_TO_RESERVE]     = npp_to_alloc;
		}
		else
		{
			double leaf_froot_ratio;
			double leaf_reserve_to_remove;
			double leaf_litter_to_remove;
			double leaf_to_remove;
			double froot_reserve_to_remove;
			double froot_litter_to_remove;
			double froot_to_remove;

			printf("Defoliation (negative NPP)\n");

			logger(g_debug_log, "Defoliation (negative NPP)\n");
			/* note: see Jaquet et al., 2014 Tree Phys. */

			/* consuming reserve carbon pools */
			s->value[C_TO_RESERVE]     = npp_to_alloc;

			/* TEST special case when DEFOLIATION happens */

			leaf_froot_ratio = s->value[LEAF_C] / ( s->value[LEAF_C] + s->value[FROOT_C] );

			/* leaf and fine rooot carbon to litter and to reserve for respiration demand */
			/* to reserve pool */
			leaf_reserve_to_remove  = fabs(npp_to_alloc * leaf_froot_ratio);
			froot_reserve_to_remove = fabs(npp_to_alloc * (1. - leaf_froot_ratio));

			/* to litterfall */
			leaf_litter_to_remove   = fabs(leaf_reserve_to_remove  * (1. - C_FRAC_TO_RETRANSL));
			froot_litter_to_remove  = fabs(froot_reserve_to_remove * (1. - C_FRAC_TO_RETRANSL));

			/* overall */
			leaf_to_remove          = leaf_reserve_to_remove  + leaf_litter_to_remove;
			froot_to_remove         = froot_reserve_to_remove + froot_litter_to_remove;

			if ( ( s->value[LEAF_C] > leaf_to_remove ) && ( s->value[FROOT_C] > froot_to_remove ) )
			{
				/* leaf carbon to remove and to retranslocate for respiration demand */
				s->value[C_LEAF_TO_RESERVE]  += leaf_reserve_to_remove;
				s->value[C_LEAF_TO_LITR]     += leaf_litter_to_remove;

				/* leaf carbon to remove and to retranslocate for respiration demand */
				s->value[C_FROOT_TO_RESERVE] += froot_reserve_to_remove;
				s->value[C_FROOT_TO_LITR]    += froot_litter_to_remove;

				/* refill reserve pool with retranslocated reserve from leaf and fine root */
				s->value[C_TO_RESERVE]       += s->value[C_LEAF_TO_RESERVE] + s->value[C_FROOT_TO_RESERVE];
			}
		}
#endif
	}

	/* check for Maximum DBH */
	/* if dbh exceeds its maximum all C goes to branches */
	if ( d->value > DBH_MAX )
	{
		s->value[C_TO_BRANCH]    += s->value[C_TO_STEM];
		s->value[C_TO_STEM]       = 0.;
	}

#endif


	logger(g_debug_log, "C_TO_LEAF    = %f tC/cell\n", s->value[C_TO_LEAF]);
	logger(g_debug_log, "C_TO_FROOT   = %f tC/cell\n", s->value[C_TO_FROOT]);
	logger(g_debug_log, "C_TO_STEM    = %f tC/cell\n", s->value[C_TO_STEM]);
	logger(g_debug_log, "C_TO_CROOT   = %f tC/cell\n", s->value[C_TO_CROOT]);
	logger(g_debug_log, "C_TO_BRANCH  = %f tC/cell\n", s->value[C_TO_BRANCH]);
	logger(g_debug_log, "C_TO_RESERVE = %f tC/cell\n", s->value[C_TO_RESERVE]);

	if ( s->value[C_TO_LEAF]     > 0. ) delta_leaf     = s->value[C_TO_LEAF];
	else delta_leaf     = 0.;
	if ( s->value[C_TO_FROOT]    > 0. ) delta_froot    = s->value[C_TO_FROOT];
	else delta_froot    = 0.;
	if ( s->value[C_TO_STEM]     > 0. ) delta_stem     = s->value[C_TO_STEM];
	else delta_stem     = 0.;
	if ( s->value[C_TO_CROOT]    > 0. ) delta_croot    = s->value[C_TO_CROOT];
	else delta_croot    = 0.;
	if ( s->value[C_TO_BRANCH]   > 0. ) delta_branch   = s->value[C_TO_BRANCH];
	else delta_branch   = 0.;
	if ( s->value[C_TO_FRUIT]    > 0. ) delta_fruit    = s->value[C_TO_FRUIT];
	else delta_fruit    = 0.;
	if ( s->value[C_TO_RESERVE]  > 0. ) delta_reserve_alloc  = s->value[C_TO_RESERVE];
	else delta_reserve_deple  = fabs(s->value[C_TO_RESERVE]);


	s->value[YEARLY_RESERVE_ALLOC] += ( delta_reserve_alloc * 1e6 / g_settings->sizeCell );
	s->value[YEARLY_RESERVE_DEPLE] += ( delta_reserve_deple * 1e6 / g_settings->sizeCell );

	s->value[YEARLY_RESERVE_USAGE] += ( s->value[C_TO_RESERVE] * 1e6 / g_settings->sizeCell );

	/* biomass production */
	s->value[BP] += ( ( delta_leaf + delta_froot + delta_stem + delta_branch + delta_croot + delta_fruit ) * 1e6 / g_settings->sizeCell );

	s->value[BP] *= ( 1. - s->value[EFF_GRPERC] );

	/* leaf and fine root fall */
	leaffall_evergreen ( c, height, dbh, age, species, year );

	/* allocate daily carbon */
	live_total_wood_age ( a, s );
}

/**/

