/* includes */
#include "CN-allocation.h"  //ddalmo 2024 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "C-deciduous-partitioning.h"
#include "common.h"
#include "constants.h"
#include "settings.h"
#include "logger.h"
#include "leaf_fall.h"
#include "biomass.h"
#include "dendometry.h"
#include "turnover.h"
#include "lai.h"
#include "mortality.h"
#include "aut_respiration.h"

extern logger_t* g_debug_log;
extern settings_t* g_settings;

/* Deciduous carbon allocation routine */
void daily_C_deciduous_partitioning (cell_t *const c, const int layer, const int height, const int dbh, const int age, const int species,
		const meteo_daily_t *const meteo_daily, const int day, const int month, const int year)
{
	double s0;
	double r0;
	double omega;
	double pS;
	double pR;
	double pL;
	double Light_trasm;
	int i;

	double delta_leaf     = 0.;
	double delta_froot    = 0.;
	double delta_croot    = 0.;
	double delta_branch   = 0.;
	double delta_stem     = 0.;
	double delta_fruit    = 0.;
	double delta_reserve_alloc  = 0.;
	double delta_reserve_deple  = 0.;

	/* for check */
	double npp_to_alloc;

	height_t *h;
	dbh_t *d;
	age_t *a;
	species_t *s;

	h = &c->heights[height];
	d = &h->dbhs[dbh];
	a = &d->ages[age];
	s = &a->species[species];

	s0    = s->value[S0CTEM];        /* parameter controlling allocation to stem (minimum ratio to stem pool */
	r0    = s->value[R0CTEM];        /* parameter controlling allocation to roots (minimum ratio to root pools */
	omega = s->value[OMEGA_CTEM];    /* controls the sensitivity of allocation to changes in water and light availability */

	//FIXME it should takes into account above layers
	Light_trasm = exp ( - s->value[K] * s->value[LAI_PROJ] );

	/* note: in Biome a constant proportion (50%) (Growth:storage parameter) of NPP that goes to the c-pools is allocated
	 * to each storage_pool, i.e. each carbon pools receive just a part of NPP (50%) the remaining remain as storage
	 * and used to maintain trees when NPP is < 0 */

	logger(g_debug_log, "\n**C-PARTITIONING**\n");
	logger(g_debug_log, "Carbon partitioning for deciduous\n");

	/* partitioning block using approach of:
	 *  Potter et al., 1993, Friendlingstein et al., 1999;
	 *  Schwalm & Ek, 2004; Arora and Boer 2005 */

	/* roots */
	pR = ( r0 + ( omega * ( 1. - s->value[F_SW] ) ) ) / ( 1. + (omega * ( 2. - Light_trasm - s->value[F_SW] ) ) );
	//logger(g_debug_log, "Roots CTEM ratio = %g %%\n", pR * 100);

	/* stem */
	pS = ( s0 + ( omega * ( 1. - Light_trasm ) ) ) / ( 1. + ( omega * ( 2. - Light_trasm - s->value[F_SW] ) ) );
	//logger(g_debug_log, "Stem CTEM ratio = %g %%\n", pS * 100. );

	/* reserve and leaves */
	pL = ( 1. - pS - pR );
	//logger(g_debug_log, "Reserve CTEM ratio = %g %%\n", pL * 100. );
	CHECK_CONDITION( fabs ( pR + pS + pL ), > , 1 + eps );

	if (s->counter[VEG_DAYS] == 1)
	{
		s->counter[BUD_BURST_COUNTER] = (int)s->value[BUD_BURST];
	}
	if (s->counter[VEG_DAYS] > 1 && s->counter[VEG_DAYS] <= s->value[BUD_BURST])
	{
		s->counter[BUD_BURST_COUNTER] --;
	}
	if (s->counter[VEG_DAYS] > s->value[BUD_BURST])
	{
		s->counter[BUD_BURST_COUNTER] = 0;
	}

	/* check for respiration */
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

	/* note: none carbon pool is refilled if reserve is lower than minimum */
	/* reserves have priority before all other pools */
	/* "...species-specific minimum pool (of NSC) sizes of stored C may be necessary to avoid catastrophic xylem failure" 
	   Sala et al., 2012, Tree Physiology */ 

	switch ( s->phenology_phase )
	{
	/************************************************************/
	case 1:

		/* if management doens't happen (this to avoid problems in carbon balance) */
		if ( ! s->counter[THINNING_HAPPENS] )
		{
			/* days for budburst */
			s->counter[DAYS_FOR_BUDBURST] = s->value[BUD_BURST] - 1.;

#if 1
		
			//FIXME these should be all class dependent variables
			//somehow there's a time lag by which model reduces allocation before to finish to allocate
			if ( ! s->counter[BUD_BURST_DAY_COUNTER] )
			{
				s->counter[BUDBURST_A]       = 0;
				s->counter[BUDBURST_B]       = 0;
				s->counter[BUD_BURST_WEIGHT] = 0;

				for ( i = 0; i < s->counter[DAYS_FOR_BUDBURST]; ++i )
				{
					s->counter[BUDBURST_A]       += 1;
					s->counter[BUDBURST_B]       += s->counter[BUDBURST_A];
					s->counter[BUD_BURST_WEIGHT] += s->counter[BUDBURST_B];
				}
				s->counter[BUDBURST_A] = 0;
				s->counter[BUDBURST_B] = 0;

				s->value[TOT_C_RESERVE_TO_LEAF_BUDBURST]  = ( s->value[MAX_LEAF_C]   + ( s->value[MAX_LEAF_C]  * s->value[EFF_GRPERC] ) );
				s->value[TOT_C_RESERVE_TO_FROOT_BUDBURST] = ( s->value[MAX_FROOT_C]  + ( s->value[MAX_FROOT_C] * s->value[EFF_GRPERC] ) );
			}

			s->counter[BUDBURST_A] +=1;
			s->counter[BUDBURST_B] += s->counter[BUDBURST_A];

			s->value[C_RESERVE_TO_LEAF_BUDBURST]  = s->value[TOT_C_RESERVE_TO_LEAF_BUDBURST]  * ( s->counter[BUDBURST_B] / (double)s->counter[BUD_BURST_WEIGHT] );
			s->value[C_RESERVE_TO_FROOT_BUDBURST] = s->value[TOT_C_RESERVE_TO_FROOT_BUDBURST] * ( s->counter[BUDBURST_B] / (double)s->counter[BUD_BURST_WEIGHT] );

			s->counter[BUD_BURST_DAY_COUNTER] ++;

#else
			/* "This has recently been confirmed by Dyckmans et al. (2000)
			who showed that only 44% of carbon in leaves came from
			carbon reserves in beech trees" (Barbaroux et al., 2003) */

			/* following Friedlingstein et al.,1998 and Krinner et al.,2005 during budburst model allocates to leaf and fine root */
			/* following Campioli et al., 2008, Maillard et al., 1994, Barbaroux et al., 2003 */

			/* compute amount of leaf carbon and relative growth respiration amount */
			daily_reserve_to_leaf_budburst   = ( s->value[MAX_LEAF_C]   + ( s->value[MAX_LEAF_C]  * s->value[EFF_GRPERC] ) )  / days_for_budburst;

			/* compute amount of fine root carbon and relative growth respiration amount */
			daily_reserve_to_froot_budburst  = ( s->value[MAX_FROOT_C]  + ( s->value[MAX_FROOT_C] * s->value[EFF_GRPERC] ) )  / days_for_budburst;

#endif
			/* compute reserve needed for budburst */
			s->value[C_RESERVE_TO_BUDBURST]  = s->value[C_RESERVE_TO_LEAF_BUDBURST] + s->value[C_RESERVE_TO_FROOT_BUDBURST];

			/* update carbon flux */
			s->value[C_TO_LEAF]        = s->value[C_RESERVE_TO_LEAF_BUDBURST];
			s->value[C_TO_FROOT]       = s->value[C_RESERVE_TO_FROOT_BUDBURST];
			s->value[C_TO_RESERVE]     = npp_to_alloc - s->value[C_RESERVE_TO_BUDBURST];

			/**********************************************************************/
			/* check for leaf C > max leaf C */
			if ( ( ( s->value[C_TO_LEAF] * ( 1. - s->value[EFF_GRPERC] ) ) + s->value[LEAF_C]) > s->value[MAX_LEAF_C] )
			{
				double max_leafC;

				max_leafC = s->value[MAX_LEAF_C] - s->value[LEAF_C];

				s->value[C_TO_LEAF]     = max_leafC + ( max_leafC * s->value[EFF_GRPERC] );
				s->value[C_TO_RESERVE] += ( s->value[C_RESERVE_TO_LEAF_BUDBURST] - ( max_leafC + ( max_leafC * s->value[EFF_GRPERC] ) ) );
			}

			/* check for fine root C > max fine root C */
			if ( ( (s->value[C_TO_FROOT] * ( 1. - s->value[EFF_GRPERC] ) ) + s->value[FROOT_C] ) > s->value[MAX_FROOT_C] )
			{
				double max_frootC;

				max_frootC = s->value[MAX_FROOT_C] - s->value[FROOT_C];

				s->value[C_TO_FROOT]    = max_frootC + ( max_frootC * s->value[EFF_GRPERC] );
				s->value[C_TO_RESERVE] += ( s->value[C_RESERVE_TO_FROOT_BUDBURST]  - ( max_frootC + ( max_frootC * s->value[EFF_GRPERC] ) ) );
			}
		}
		else
		{
			logger(g_debug_log, "Allocating only into reserve pool (special case thinning happens on)\n");
			s->value[C_TO_RESERVE] = npp_to_alloc;
		}
		/**********************************************************************/

		break;

		/********************************************************************************************************************************************/
	case 2:
		logger(g_debug_log, "\n*NORMAL GROWTH*\n");
		logger(g_debug_log, "(LAI == PEAK LAI)\n");

		/* reset budburst day counter */
		s->counter[BUD_BURST_DAY_COUNTER] = 0.;

		/* see Barbaroux et al., 2002, Scartazza et al., 2013 */

		if (npp_to_alloc > 0.0)
		{
			/* check if minimum reserve pool needs to be refilled */
			/* see Litton et al., GCB, 2007 */

			/* it doesn't need */
			if( s->value[RESERVE_C] >= s->value[MIN_RESERVE_C] )
			{
				/* as in Bossel (1996), second priority is assigned to fruit */
				/* reproduction */
				if ( a->value > s->value[SEXAGE] )
				{
					logger(g_debug_log, "allocating into fruit pool\n");

					s->value[C_TO_FRUIT] = npp_to_alloc * s->value[FRUIT_PERC];
					npp_to_alloc        -= s->value[C_TO_FRUIT];
				}
				else
				{
					s->value[C_TO_FRUIT] = 0.;
				}

				logger(g_debug_log, "allocating into the three pools stem, branch, coarse root and reserve\n");

				/* allocating into c pools */
				/* note: for references see: Potter et al., 1993; Schwalm and Ek 2004 */
				s->value[C_TO_RESERVE]   = (npp_to_alloc * pL);
				s->value[C_TO_CROOT]     = (npp_to_alloc * pR);
				s->value[C_TO_STEM]      = (npp_to_alloc * pS) * ( 1. - s->value[FRACBB] );
				s->value[C_TO_BRANCH]    = (npp_to_alloc * pS) * s->value[FRACBB];

			}
			/* it needs */
			else
			{
				logger(g_debug_log, "allocating into reserve\n");

				/* allocating into c pools */
				s->value[C_TO_RESERVE]     = npp_to_alloc;
			}
		}
		else  // npp < 0
		{
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

				logger(g_debug_log, "Defoliation (negative NPP)\n");
				/* note: see Jaquet et al., 2014 Tree Phys. */

				/* consuming reserve carbon pools */
				s->value[C_TO_RESERVE]     = npp_to_alloc;

				/* special case when DEFOLIATION happens */

				leaf_froot_ratio = s->value[LEAF_C] / ( s->value[LEAF_C] + s->value[FROOT_C] );

				/* leaf and fine rooot carbon to litter and to reserve for respiration demand */
				/* to reserve pool */
				leaf_reserve_to_remove  = fabs(npp_to_alloc * leaf_froot_ratio);
				froot_reserve_to_remove = fabs(npp_to_alloc * ( 1. - leaf_froot_ratio ) );

				/* to litterfall */
				leaf_litter_to_remove   = fabs(leaf_reserve_to_remove  * ( 1. - C_FRAC_TO_RETRANSL ) );
				froot_litter_to_remove  = fabs(froot_reserve_to_remove * ( 1. - C_FRAC_TO_RETRANSL ) );

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
		}

		break;
		/**********************************************************************/
	case 3:
		logger(g_debug_log, "\n*LEAF FALL*\n");
		logger(g_debug_log, "(DayLength < MINDAYLENGTH)\n");
		logger(g_debug_log, "allocating into W reserve pool\n");

		/* allocating into c pools */
		/* including retranslocated C */

		s->value[C_TO_RESERVE] = npp_to_alloc ;

		/* leaf fall */
		leaffall_deciduous( c, height, dbh, age, species );

		break;
		/**********************************************************************/
	case 0:
		logger(g_debug_log, "\n*DORMANT SEASON*\n");
		logger(g_debug_log, "consuming reserve pool\n");

		/* allocating into c pools */
		s->value[C_TO_RESERVE] = npp_to_alloc;

		break;
		/**********************************************************************/
	}

	/* check for Maximum DBH */
	/* if dbh exceeds its maximum all C goes to branches */
	
	if ( d->value > DBH_MAX )
	{
		s->value[C_TO_BRANCH]    += s->value[C_TO_STEM];
		s->value[C_TO_STEM]       = 0.;
	}

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

	/* update live_total wood fraction based on age */
	live_total_wood_age ( a, s );
}




