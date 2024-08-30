/* biomass.c */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "biomass.h"
#include "common.h"
#include "constants.h"
#include "settings.h"
#include "logger.h"

extern logger_t* g_debug_log;
extern settings_t* g_settings;

void live_total_wood_age(const age_t *const a, species_t *const s)
{
#if 0
	/* this function update based on current tree age the amount of live:total wood ratio
	 * based on the assumption that the live wood decrease linearly increasing age */
	/* e.g. for Fagus sylvatica  base on simulation for Hesse site (age 30) and Collelongo site (age 160)*/
	// fixme values should be included in species.txt
	double max_live_total_ratio;
	double min_live_total_ratio;
	double t1;
	double t2;

	int min_age;
	int max_age;

	logger(g_debug_log, "\n*live:total wood ratio based on tree age*\n");

	max_age = (int)s->value[MAXAGE];

	//fixme values should be included in species.txt
	/* age at which live:total wood ratio is maximum */
	min_age = 25;

	max_live_total_ratio = s->value[LIVE_TOTAL_WOOD];

	//fixme values should be included in species.txt
	min_live_total_ratio = 0.03;

	t1 = max_live_total_ratio - min_live_total_ratio;
	t2 = max_age - min_age;

	s->value[EFF_LIVE_TOTAL_WOOD_FRAC] = (t1/t2)*(max_age - a->value) + min_live_total_ratio;
#else
	s->value[EFF_LIVE_TOTAL_WOOD_FRAC] = s->value[LIVE_TOTAL_WOOD];
#endif
	logger(g_debug_log, "Effective live:total wood fraction based on stand age = %g\n", s->value[EFF_LIVE_TOTAL_WOOD_FRAC]);
}


void abg_bgb_biomass(cell_t *const c, const int height, const int dbh, const int age, const int species)
{
	species_t *s;

	//ALESSIOR
	/*
	if ( ! c->heights ) return;
	if ( ! c->heights[height].dbhs ) return;
	if ( ! c->heights[height].dbhs[dbh].ages ) return;
	if ( ! c->heights[height].dbhs[dbh].ages[age].species ) return;
	 */
	s = &c->heights[height].dbhs[dbh].ages[age].species[species];

	logger(g_debug_log, "**AGB & BGB**\n");

	s->value[OLD_AGB]      = s->value[AGB];
	s->value[OLD_BGB]      = s->value[BGB];
	s->value[TREE_OLD_AGB] = s->value[TREE_AGB];
	s->value[TREE_OLD_BGB] = s->value[TREE_BGB];

	s->value[AGB]            = s->value[LEAF_C]  + s->value[STEM_C] + s->value[BRANCH_C] +s->value[FRUIT_C];
	logger(g_debug_log, "AGB              = %f tC/cell\n", s->value[AGB]);

	s->value[BGB]            = s->value[FROOT_C] + s->value[CROOT_C];
	logger(g_debug_log, "BGB              = %f tC/cell\n", s->value[BGB]);

	s->value[STANDING_WOOD]  = s->value[STEM_C] + s->value[BRANCH_C] + s->value[CROOT_C];
	logger(g_debug_log, "STANDING_WOOD    = %f tC/cell\n", s->value[STANDING_WOOD]);

	s->value[DELTA_AGB]      = s->value[AGB] - s->value[OLD_AGB];
	logger(g_debug_log, "DELTA_AGB        = %f tC/cell/year\n", s->value[DELTA_AGB]);

	s->value[DELTA_BGB]      = s->value[BGB] - s->value[OLD_BGB];
	logger(g_debug_log, "DELTA_BGB        = %f tC/cell/year\n", s->value[DELTA_BGB]);

	s->value[TREE_AGB]       = s->value[AGB] / (double)s->counter[N_TREE];
	logger(g_debug_log, "Yearly Class AGB = %f tC/tree\n", s->value[TREE_AGB]);

	s->value[TREE_BGB]       = s->value[BGB] / (double)s->counter[N_TREE];
	logger(g_debug_log, "Yearly Class BGB = %f tC/tree\n", s->value[TREE_BGB]);

	s->value[DELTA_TREE_AGB] = s->value[TREE_AGB] - s->value[TREE_OLD_AGB];
	logger(g_debug_log, "DELTA_TREE_AGB   = %f tC/tree/year\n", s->value[DELTA_TREE_AGB]);

	s->value[DELTA_TREE_BGB] = s->value[TREE_BGB] - s->value[TREE_OLD_BGB];
	logger(g_debug_log, "DELTA_TREE_BGB   = %f tC/tree/year\n", s->value[DELTA_TREE_BGB]);
}

void average_tree_pools(cell_t *const c)
{
	int height;
	int dbh;
	int age;
	int species;

	species_t *s;

	logger(g_debug_log, "\n*AVERAGE TREE POOLS*\n");
	for ( height = 0; height < c->heights_count ; ++height )
	{

		for ( dbh = 0; dbh < c->heights[height].dbhs_count; ++dbh )
		{
			for ( age = 0; age < c->heights[height].dbhs[dbh].ages_count ; ++age )
			{
				for ( species = 0; species < c->heights[height].dbhs[dbh].ages[age].species_count; ++species )
				{
					s = &c->heights[height].dbhs[dbh].ages[age].species[species];

					logger(g_debug_log,  "N_TREE = %d\n", s->counter[N_TREE]);

					/* note: it takes into account previous biomass to remove if function is called twice in a day */


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

					/* compute tree average N pools */
					s->value[TREE_LEAF_N]                = (s->value[LEAF_N]             / (double)s->counter[N_TREE]);
					s->value[TREE_STEM_N]                = (s->value[STEM_N]             / (double)s->counter[N_TREE]);
					s->value[TREE_FROOT_N]               = (s->value[FROOT_N]            / (double)s->counter[N_TREE]);
					s->value[TREE_CROOT_N]               = (s->value[CROOT_N]            / (double)s->counter[N_TREE]);
					s->value[TREE_BRANCH_N]              = (s->value[BRANCH_N]           / (double)s->counter[N_TREE]);
					s->value[TREE_RESERVE_N]             = (s->value[RESERVE_N]          / (double)s->counter[N_TREE]);
					s->value[TREE_FRUIT_N]               = (s->value[FRUIT_N]            / (double)s->counter[N_TREE]);

				}
			}
		}
	}
}


void tree_branch_and_bark (cell_t *const c, const int height, const int dbh, const int age, const int species)
{
	age_t *a;
	species_t *s;

	a = &c->heights[height].dbhs[dbh].ages[age];
	// ALESSIOR
	//if ( ! a ) return;
	s = &c->heights[height].dbhs[dbh].ages[age].species[species];

	/* compute branch and bark fractions */

	if ( !s->value[FRACBB0] ) s->value[FRACBB0] = s->value[FRACBB1];
	else s->value[FRACBB] = s->value[FRACBB1] + ( s->value[FRACBB0] - s->value[FRACBB1] )* exp( -LN2 * ( (double)a->value / s->value[TBB] ) );
}

void tree_biomass_remove (cell_t *const c, const int height, const int dbh, const int age, const int species, const int tree_remove, const int nat_man)
{
	species_t *s;
	s = &c->heights[height].dbhs[dbh].ages[age].species[species];

	logger(g_debug_log, "\n* TREE BIOMASS REMOVE *\n");

	/******************************************************************************************/

	/* carbon to litter fluxes */
	s->value[C_LEAF_TO_LITR]     += (s->value[TREE_LEAF_C]    * tree_remove);

	s->value[C_FROOT_TO_LITR]    += (s->value[TREE_FROOT_C]   * tree_remove);

	/* overall litter */
	s->value[C_TO_LITR]          += s->value[C_LEAF_TO_LITR] + s->value[C_FROOT_TO_LITR];

	/* carbon to cwd fluxes */
	s->value[C_STEM_TO_CWD]      += (s->value[TREE_STEM_C]    * tree_remove);

	s->value[C_CROOT_TO_CWD]     += (s->value[TREE_CROOT_C]   * tree_remove);

	s->value[C_BRANCH_TO_CWD]    += (s->value[TREE_BRANCH_C]  * tree_remove);

	s->value[C_RESERVE_TO_CWD]   += (s->value[TREE_RESERVE_C] * tree_remove);

	s->value[C_FRUIT_TO_CWD]     += (s->value[TREE_FRUIT_C]   * tree_remove);

	/******************************************************************************************/

	/* sapwood and heartwood */
	/* note: this mostly differentiates the 5.4. to 5.5 (and subsequests) versions */

	s->value[C_STEM_SAPWOOD_TO_CWD]     += (s->value[TREE_STEM_SAPWOOD_C]     * tree_remove);

	s->value[C_CROOT_SAPWOOD_TO_CWD]    += (s->value[TREE_CROOT_SAPWOOD_C]    * tree_remove);

	s->value[C_BRANCH_SAPWOOD_TO_CWD]   += (s->value[TREE_BRANCH_SAPWOOD_C]   * tree_remove);

	s->value[C_STEM_HEARTWOOD_TO_CWD]   += (s->value[TREE_STEM_HEARTWOOD_C]   * tree_remove);

	s->value[C_CROOT_HEARTWOOD_TO_CWD]  += (s->value[TREE_CROOT_HEARTWOOD_C]  * tree_remove);

	s->value[C_BRANCH_HEARTWOOD_TO_CWD] += (s->value[TREE_BRANCH_HEARTWOOD_C] * tree_remove);

	/******************************************************************************************/

	/* livewood and deadwood */

	s->value[C_STEM_LIVEWOOD_TO_CWD]    += (s->value[TREE_STEM_LIVEWOOD_C]     * tree_remove);

	s->value[C_CROOT_LIVEWOOD_TO_CWD]   += (s->value[TREE_CROOT_LIVEWOOD_C]    * tree_remove);

	s->value[C_BRANCH_LIVEWOOD_TO_CWD]  += (s->value[TREE_BRANCH_LIVEWOOD_C]   * tree_remove);

	s->value[C_STEM_DEADWOOD_TO_CWD]    += (s->value[TREE_STEM_DEADWOOD_C]     * tree_remove);

	s->value[C_CROOT_DEADWOOD_TO_CWD]   += (s->value[TREE_CROOT_DEADWOOD_C]    * tree_remove);

	s->value[C_BRANCH_DEADWOOD_TO_CWD]  += (s->value[TREE_BRANCH_DEADWOOD_C]   * tree_remove);

	/******************************************************************************************/
    /* update class carbon pools */

	s->value[LEAF_C]    -= (s->value[TREE_LEAF_C]    * tree_remove);
	s->value[FROOT_C]   -= (s->value[TREE_FROOT_C]   * tree_remove);

	s->value[STEM_C]      -= (s->value[TREE_STEM_C]    * tree_remove);
	s->value[CROOT_C]     -= (s->value[TREE_CROOT_C]   * tree_remove);
	s->value[BRANCH_C]    -= (s->value[TREE_BRANCH_C]  * tree_remove);

	s->value[RESERVE_C]   -= (s->value[TREE_RESERVE_C] * tree_remove);
	s->value[FRUIT_C]     -= (s->value[TREE_FRUIT_C]   * tree_remove);

	/* sapwood and heartwood */

	s->value[STEM_SAPWOOD_C]     -= (s->value[TREE_STEM_SAPWOOD_C]     * tree_remove);
	s->value[CROOT_SAPWOOD_C]    -= (s->value[TREE_CROOT_SAPWOOD_C]    * tree_remove);
	s->value[BRANCH_SAPWOOD_C]   -= (s->value[TREE_BRANCH_SAPWOOD_C]   * tree_remove);
	s->value[STEM_HEARTWOOD_C]   -= (s->value[TREE_STEM_HEARTWOOD_C]   * tree_remove);
	s->value[CROOT_HEARTWOOD_C]  -= (s->value[TREE_CROOT_HEARTWOOD_C]  * tree_remove);
	s->value[BRANCH_HEARTWOOD_C] -= (s->value[TREE_BRANCH_HEARTWOOD_C] * tree_remove);

	/* livewood and deadwood */

	s->value[STEM_LIVEWOOD_C]    -= (s->value[TREE_STEM_LIVEWOOD_C]     * tree_remove);
	s->value[CROOT_LIVEWOOD_C]   -= (s->value[TREE_CROOT_LIVEWOOD_C]    * tree_remove);
	s->value[BRANCH_LIVEWOOD_C]  -= (s->value[TREE_BRANCH_LIVEWOOD_C]   * tree_remove);
	s->value[STEM_DEADWOOD_C]    -= (s->value[TREE_STEM_DEADWOOD_C]     * tree_remove);
	s->value[CROOT_DEADWOOD_C]   -= (s->value[TREE_CROOT_DEADWOOD_C]    * tree_remove);
	s->value[BRANCH_DEADWOOD_C]  -= (s->value[TREE_BRANCH_DEADWOOD_C]   * tree_remove);


	// Aggregate the fluxes toward CWD pool 

	if ( ! nat_man )
	{
		/* natural mortality */
		/* overall cwd */

		if ( s->value[PHENOLOGY] == 0.1 || s->value[PHENOLOGY] == 0.2 )
		{
			/* deciduous */
			s->value[C_TO_CWD]                += s->value[C_STEM_TO_CWD] +
					s->value[C_CROOT_TO_CWD]       +
					s->value[C_BRANCH_TO_CWD]      +
					s->value[C_RESERVE_TO_CWD]     +
					s->value[C_FRUIT_TO_CWD]       ;
		}
		else
		{
			/* evergreen */
			s->value[C_TO_CWD]                += s->value[C_STEM_TO_CWD] +
					s->value[C_CROOT_TO_CWD]       +
					s->value[C_BRANCH_TO_CWD]      +
					s->value[C_RESERVE_TO_CWD]     +
					s->value[C_FRUIT_TO_CWD]       +
					s->value[C_LEAF_TO_LITR]       +
					s->value[C_FROOT_TO_LITR]      ;
		}
	}
	else
	{
		/* managed mortality */
		/* note: we consider that just tree stems are removed from stand during management */  //ddalmo note: NOT TRUE:see below!
		/* overall cwd */
		if ( s->value[PHENOLOGY] == 0.1 || s->value[PHENOLOGY] == 0.2 )
		{
			/* deciduous */
			s->value[C_TO_CWD]                += s->value[C_CROOT_TO_CWD]  +
					s->value[C_BRANCH_TO_CWD] +
					s->value[C_RESERVE_TO_CWD]+
					s->value[C_FRUIT_TO_CWD] ;
		}
		else
		{
			/* evergreen */
			s->value[C_TO_CWD]                += s->value[C_CROOT_TO_CWD]  +
					s->value[C_BRANCH_TO_CWD] +
					s->value[C_RESERVE_TO_CWD]+
					s->value[C_FRUIT_TO_CWD]  +
					s->value[C_LEAF_TO_LITR]  +
					s->value[C_FROOT_TO_LITR] ;
		}
	}

	if ( s->counter[HARVESTING_HAPPENS] == 1 ) 
	{      // TODO check this! as this should also be the case for thinning, if only stem is removed..and mortality in general. for N cycle too.
	
		/*** compute class-level Coarse Woody Debris carbon fluxes (tC/sizecell/day) ****/
		s->value[CWD_TO_LITRC]           = s->value[C_TO_CWD];
		s->value[CWD_TO_LITR2C]          = s->value[C_TO_CWD] * s->value[DEADWOOD_USCEL_FRAC];
		s->value[CWD_TO_LITR3C]          = s->value[C_TO_CWD] * s->value[DEADWOOD_SCEL_FRAC];
		s->value[CWD_TO_LITR4C]          = s->value[C_TO_CWD] * s->value[DEADWOOD_LIGN_FRAC];
		/* check */
		CHECK_CONDITION ( s->value[CWD_TO_LITR2C] + s->value[CWD_TO_LITR3C] + s->value[CWD_TO_LITR4C] , == , s->value[CWD_TO_LITRC] + eps );

		/*** compute class-level Coarse Woody Debris carbon pools (tC/sizecell) ****/
		s->value[CWD_LITRC]             += s->value[C_TO_CWD];
		s->value[CWD_LITR2C]            += s->value[CWD_TO_LITR2C];
		s->value[CWD_LITR3C]            += s->value[CWD_TO_LITR3C];
		s->value[CWD_LITR4C]            += s->value[CWD_TO_LITR4C];
		/* check */
		CHECK_CONDITION ( s->value[CWD_LITR2C] + s->value[CWD_LITR3C] + s->value[CWD_LITR4C] , == , s->value[CWD_LITRC] + eps );

		/*** update cell-level Coarse Woody Debris carbon fluxes (gC/m2/day) ****/
		c->daily_cwd_to_litrC           += s->value[CWD_LITRC]  * 1e6 / g_settings->sizeCell;
		c->daily_cwd_to_litr2C          += s->value[CWD_LITR2C] * 1e6 / g_settings->sizeCell;
		c->daily_cwd_to_litr3C          += s->value[CWD_LITR3C] * 1e6 / g_settings->sizeCell;
		c->daily_cwd_to_litr4C          += s->value[CWD_LITR4C] * 1e6 / g_settings->sizeCell;
		/* check */
		CHECK_CONDITION ( c->daily_cwd_to_litr2C + c->daily_cwd_to_litr3C + c->daily_cwd_to_litr4C , == , c->daily_cwd_to_litrC + eps );

		/*** compute cell-level Coarse Woody Debris carbon pools (gC/m2) ***/
		c->cwd_C                        += s->value[CWD_LITRC]  * 1e6 / g_settings->sizeCell;
		c->cwd_2C                       += s->value[CWD_LITR2C] * 1e6 / g_settings->sizeCell;
		c->cwd_3C                       += s->value[CWD_LITR3C] * 1e6 / g_settings->sizeCell;
		c->cwd_4C                       += s->value[CWD_LITR4C] * 1e6 / g_settings->sizeCell;
		/* check */
		CHECK_CONDITION ( c->cwd_2C + c->cwd_3C + c->cwd_4C , == , c->cwd_C + eps );
	}

	/******************************************************************************************/

	/* accounting for harvested/thinned wood products (HWP) */
	if ( s->counter[THINNING_HAPPENS] == 1 || s->counter[HARVESTING_HAPPENS] == 1 )   // TODO include the condition nat_man=1  
	{
		/* compute woody biomass removed (tC/ha/yr) */
		s->value[C_HWP]          += s->value[C_STEM_TO_CWD] + s->value[C_CROOT_TO_CWD] + s->value[C_BRANCH_TO_CWD] ;
		s->value[CUM_C_HWP]      += s->value[C_HWP];
		/* compute stem volume removed (m3/ha/yr) */
		s->value[VOLUME_HWP]     += s->value[TREE_VOLUME] * tree_remove;
		s->value[CUM_VOLUME_HWP] += s->value[TREE_VOLUME] * tree_remove;
	}

	/******************************************************************************************/

	/* note: special case for turnover when mortality and thinning management happen */
	/* removing biomass to NOT consider in turnover of the subsequent year */
//	s->value[YEARLY_C_TO_STEM]            -= ((s->value[YEARLY_C_TO_STEM]            / s->counter[N_TREE]) * tree_remove);
//	s->value[YEARLY_C_TO_CROOT]           -= ((s->value[YEARLY_C_TO_CROOT]           / s->counter[N_TREE]) * tree_remove);
//	s->value[YEARLY_C_TO_BRANCH]          -= ((s->value[YEARLY_C_TO_BRANCH]          / s->counter[N_TREE]) * tree_remove);
//
//	s->value[YEARLY_C_TO_STEM_SAPWOOD]    -= ((s->value[YEARLY_C_TO_STEM_SAPWOOD]    / s->counter[N_TREE]) * tree_remove);
//	s->value[YEARLY_C_TO_CROOT_SAPWOOD]   -= ((s->value[YEARLY_C_TO_CROOT_SAPWOOD]   / s->counter[N_TREE]) * tree_remove);
//	s->value[YEARLY_C_TO_BRANCH_SAPWOOD]  -= ((s->value[YEARLY_C_TO_BRANCH_SAPWOOD]  / s->counter[N_TREE]) * tree_remove);
//
//	s->value[YEARLY_C_TO_STEM_LIVEWOOD]   -= ((s->value[YEARLY_C_TO_STEM_LIVEWOOD]   / s->counter[N_TREE]) * tree_remove);
//	s->value[YEARLY_C_TO_CROOT_LIVEWOOD]  -= ((s->value[YEARLY_C_TO_CROOT_LIVEWOOD]  / s->counter[N_TREE]) * tree_remove);
//	s->value[YEARLY_C_TO_BRANCH_LIVEWOOD] -= ((s->value[YEARLY_C_TO_BRANCH_LIVEWOOD] / s->counter[N_TREE]) * tree_remove);

	/******************************************************************************************/

	/* nitrogen to litter pool */
	s->value[N_LEAF_TO_LITR]          += (s->value[TREE_LEAF_N]       * tree_remove);

	s->value[N_FROOT_TO_LITR]         += (s->value[TREE_FROOT_N]      * tree_remove);

	/* overall litter */
	s->value[N_TO_LITR]               += s->value[N_LEAF_TO_LITR] +	s->value[N_FROOT_TO_LITR];

	s->value[LITR_N]                  += s->value[N_TO_LITR];

	/* nitrogen to cwd fluxes */
	s->value[N_STEM_TO_CWD]           += (s->value[TREE_STEM_N]       * tree_remove);

	s->value[N_CROOT_TO_CWD]          += (s->value[TREE_CROOT_N]      * tree_remove);

	s->value[N_BRANCH_TO_CWD]         += (s->value[TREE_BRANCH_N]     * tree_remove);

	s->value[N_BRANCH_TO_CWD]         += (s->value[TREE_RESERVE_N]    * tree_remove);

	s->value[N_FRUIT_TO_CWD]          += (s->value[TREE_FRUIT_N]      * tree_remove);

	s->value[N_RESERVE_TO_CWD]        += (s->value[TREE_RESERVE_N]    * tree_remove);

	if ( ! nat_man )
	{
		/* natural mortality */
		/* overall cwd */

		if ( s->value[PHENOLOGY] == 0.1 || s->value[PHENOLOGY] == 0.2 )
		{
			/* deciduous */
			s->value[N_TO_CWD]                += s->value[N_STEM_TO_CWD] +
					s->value[N_CROOT_TO_CWD]  +
					s->value[N_BRANCH_TO_CWD] +
					s->value[N_RESERVE_TO_CWD]+
					s->value[N_FRUIT_TO_CWD]  ;
		}
		else
		{
			/* evergreen */
			s->value[N_TO_CWD]                += s->value[N_STEM_TO_CWD] +
					s->value[N_CROOT_TO_CWD]  +
					s->value[N_BRANCH_TO_CWD] +
					s->value[N_RESERVE_TO_CWD]+
					s->value[N_FRUIT_TO_CWD]  +
					s->value[N_LEAF_TO_LITR]       +
					s->value[N_FROOT_TO_LITR]      ;
		}
	}
	else
	{
		/* managed mortality */
		/* note: we consider that just tree stems are removed from stand during management */
		/* overall cwd */
		if ( s->value[PHENOLOGY] == 0.1 || s->value[PHENOLOGY] == 0.2 )
		{
			/* deciduous */
			s->value[N_TO_CWD]                += s->value[N_CROOT_TO_CWD]  +
					s->value[N_BRANCH_TO_CWD] +
					s->value[N_RESERVE_TO_CWD]+
					s->value[N_FRUIT_TO_CWD] ;
		}
		else
		{
			/* evergreen */
			s->value[N_TO_CWD]                += s->value[N_CROOT_TO_CWD]  +
					s->value[N_BRANCH_TO_CWD] +
					s->value[N_RESERVE_TO_CWD]+
					s->value[N_FRUIT_TO_CWD]  +
					s->value[N_LEAF_TO_LITR]       +
					s->value[N_FROOT_TO_LITR]      ;
		}
	}

	if ( s->counter[HARVESTING_HAPPENS] == 1 )
	{
		/*** compute class-level deadwood nitrogen fluxes (tN/sizecell/day) ****/
		s->value[CWD_TO_LITRN]           = s->value[N_TO_CWD];
		s->value[CWD_TO_LITR2N]          = s->value[N_TO_CWD] * s->value[DEADWOOD_USCEL_FRAC];
		s->value[CWD_TO_LITR3N]          = s->value[N_TO_CWD] * s->value[DEADWOOD_SCEL_FRAC];
		s->value[CWD_TO_LITR4N]          = s->value[N_TO_CWD] * s->value[DEADWOOD_LIGN_FRAC];
		/* check */
		CHECK_CONDITION ( s->value[CWD_TO_LITR2N] + s->value[CWD_TO_LITR3N] + s->value[CWD_TO_LITR4N] , == , s->value[CWD_TO_LITRN] + eps );

		/*** compute class-level deadwood nitrogen pools (tN/sizecell) ****/
		s->value[CWD_LITRN]             += s->value[N_TO_CWD];
		s->value[CWD_LITR2N]            += s->value[CWD_TO_LITR2N];
		s->value[CWD_LITR3N]            += s->value[CWD_TO_LITR3N];
		s->value[CWD_LITR4N]            += s->value[CWD_TO_LITR4N];
		/* check */
		CHECK_CONDITION ( s->value[CWD_LITR2N] + s->value[CWD_LITR3N] + s->value[CWD_LITR4N] , == , s->value[CWD_LITRN] + eps );

		/*** update cell-level cwd nitrogen fluxes (gN/m2/day) ****/
		c->daily_cwd_to_litrN           += s->value[CWD_LITRN]  * 1e6 / g_settings->sizeCell;
		c->daily_cwd_to_litr2N          += s->value[CWD_LITR2N] * 1e6 / g_settings->sizeCell;
		c->daily_cwd_to_litr3N          += s->value[CWD_LITR3N] * 1e6 / g_settings->sizeCell;
		c->daily_cwd_to_litr4N          += s->value[CWD_LITR4N] * 1e6 / g_settings->sizeCell;
		/* check */
		CHECK_CONDITION ( c->daily_cwd_to_litr2N + c->daily_cwd_to_litr3N + c->daily_cwd_to_litr4N , == , c->daily_cwd_to_litrN + eps );

		/*** compute cell-level cwd nitrogen pools (gN/m2) ***/
		c->cwd_N                        += s->value[CWD_LITRN]  * 1e6 / g_settings->sizeCell;
		c->cwd_2N                       += s->value[CWD_LITR2N] * 1e6 / g_settings->sizeCell;
		c->cwd_3N                       += s->value[CWD_LITR3N] * 1e6 / g_settings->sizeCell;
		c->cwd_4N                       += s->value[CWD_LITR4N] * 1e6 / g_settings->sizeCell;
		/* check */
		CHECK_CONDITION ( c->cwd_2N + c->cwd_3N + c->cwd_4N , == , c->cwd_N + eps );
	}


}
