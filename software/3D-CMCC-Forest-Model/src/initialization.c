/* initialization.c */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "constants.h"
#include "soil_settings.h"
#include "settings.h"
#include "logger.h"
#include "common.h"
#include "structure.h"
#include "biomass.h"
#include "lai.h"
#include "aut_respiration.h"
#include "test.h"

extern settings_t *g_settings;
extern logger_t* g_debug_log;
//extern soil_settings_t *g_soil_settings;

extern const char sz_err_out_of_memory[];


void initialization_forest_structure(cell_t *const c, const int day, const int month, const int year)
{
	/* note: this function from here are called just for initialization */

	logger(g_debug_log,"\n*******INITIALIZE FOREST STRUCTURE*******\n");

	/* annual forest structure */
	if ( ! annual_forest_structure ( c, year ) )
	{
		puts(sz_err_out_of_memory);
		exit(1);
	}
}

void initialization_forest_class_C (cell_t *const c, const int height, const int dbh, const int age, const int species)
{
	double r1;
	double temp_var;                                /* temporary variable */

	height_t *h;
	dbh_t *d;
	age_t *a;
	species_t *s;

	h = &c->heights[height];
	d = &c->heights[height].dbhs[dbh];
	a = &c->heights[height].dbhs[dbh].ages[age];
	s = &c->heights[height].dbhs[dbh].ages[age].species[species];

	/* check that all mandatory variables are initialized */
	CHECK_CONDITION(h->value, <=, 0);
	CHECK_CONDITION(d->value, <=, 0);
	CHECK_CONDITION(a->value, <=, 0);

	/* check parameters values */
	CHECK_CONDITION (s->value[SLA_AVG0],               <, s->value[SLA_AVG1]);
	CHECK_CONDITION (s->value[FRACBB0],                <, s->value[FRACBB1]);
	CHECK_CONDITION (s->value[GROWTHTMIN],             >, s->value[GROWTHTOPT]);
	CHECK_CONDITION (s->value[GROWTHTMIN],             >, s->value[GROWTHTMAX]);
	CHECK_CONDITION (s->value[GROWTHTOPT],             >, s->value[GROWTHTMAX]);
	CHECK_CONDITION (s->value[SWPOPEN],                <, s->value[SWPCLOSE]);
	CHECK_CONDITION (s->value[FRUIT_PERC],             >, 1.);
	CHECK_CONDITION (s->value[LEAF_FALL_FRAC_GROWING], >, 1.);
	CHECK_CONDITION (s->value[LIVEWOOD_TURNOVER],      >, 1.);
	CHECK_CONDITION ( fabs (s->value[S0CTEM] + s->value[R0CTEM] + s->value[F0CTEM]), >, 1 + eps)

	if (s->value[PHENOLOGY] == 0.1 || s->value[PHENOLOGY] == 0.2)
	{

		CHECK_CONDITION (s->value[LEAF_FROOT_TURNOVER],    !=, 1.);
	}
	else
	{
             
		CHECK_CONDITION (s->value[LEAF_FROOT_TURNOVER],     >, 1.);
	}
   

	/* compute growth respiration fraction */
	growth_respiration_frac ( a, s );

	/* compute age-related specific leaf area */
	specific_leaf_area ( a, s);

	/**************************************************************************************************/

	/* note: special case, model recomputes tree height based on its parameterizations to avoid differences between
	 * modelled and initialized tree height data */

	/* Chapman-Richards function */
	/* for references see also: R. Pilli et al. Forest Ecology and Management 237 (2006) 583–593 */
	/* note: this shouldn't be applied to saplings that are lower than 1.3 meter */
	
	h->value = DBH_ref + s->value[CRA] * pow (1. - exp ( - s->value[CRB] * d->value) , s->value[CRC]);

	if ( h->value > s->value[CRA] )
	{
		h->value = s->value[CRA];
	}
	logger(g_debug_log, "-Tree Height using Chapman-Richard function = %f m\n", h->value);

	/* check */
	CHECK_CONDITION (h->value, >, s->value[CRA] + eps )
	logger(g_debug_log, "-New Tree Height = %f m\n", h->value);

	/**************************************************************************************************/


	logger(g_debug_log,"\n*******INITIALIZE FOREST CLASS CARBON POOLS (%s)*******\n", s->name);
	logger(g_debug_log, "...checking initial biomass data for class: height %f, dbh %f, age %d, species %s...\n", h->value, d->value, a->value, s->name);

	/* compute mass density */
	s->value[MASS_DENSITY] = s->value[RHO1] + (s->value[RHO0] - s->value[RHO1]) * exp(-LN2 * ((double)a->value / s->value[TRHO]));
	logger(g_debug_log, "-Mass Density = %f tDM/m3\n", s->value[MASS_DENSITY]);

	/* compute volume */
	s->value[TREE_VOLUME] = (Pi * s->value[FORM_FACTOR] * pow((d->value / 100.) , 2.) * h->value) / 4.;
	s->value[VOLUME]      = s->value[TREE_VOLUME] * s->counter[N_TREE];

	logger(g_debug_log, "-Single tree volume = %f m3/tree\n", s->value[TREE_VOLUME]);
	logger(g_debug_log, "-Class volume       = %f m3/cell\n", s->value[VOLUME]);

	/*check for initial biomass*/
	if (s->value[STEM_DM] == 0.0 || s->value[STEM_DM] == NO_DATA)
	{
		logger(g_debug_log, "\nNo Stem Biomass Data are available for model initialization\n");
		logger(g_debug_log, "...Generating input Stem Biomass biomass data from DBH = %f cm\n", d->value);

		/* compute stem biomass from DBH */
		if (s->value[STEMCONST_P] == NO_DATA && s->value[STEMPOWER_P] == NO_DATA)
		{
			/* use generic stemconst stempower values */
			logger(g_debug_log, "..computing stem biomass from generic stempower and stemconst DBH = %f cm\n", d->value);
			if (d->value < 9)
			{
				s->value[TREE_STEM_C] = (s->value[STEMCONST] * (pow (d->value, STEMPOWER_A))) / GC_GDM / 1e3;
			}
			else if (d->value > 9 && d->value < 15)
			{
				s->value[TREE_STEM_C] = (s->value[STEMCONST] * (pow (d->value, STEMPOWER_B))) / GC_GDM / 1e3;
			}
			else
			{
				s->value[TREE_STEM_C] = (s->value[STEMCONST] * (pow (d->value, STEMPOWER_C))) / GC_GDM / 1e3;
			}
		}

		else
		{
			/* use site specific stemconst stempower values */
			logger(g_debug_log, "..computing stem biomass from site specific stempower and stemconst DBH = %f cm\n", d->value);
			s->value[TREE_STEM_C] = (s->value[STEMCONST_P] * pow (d->value, s->value[STEMPOWER_P])) / GC_GDM / 1e3;

		}
		logger(g_debug_log, "Single tree stem mass = %f tC/tree\n", s->value[TREE_STEM_C]);

		/* tree tC -> class tDM */
		s->value[STEM_DM] = s->value[TREE_STEM_C] * s->counter[N_TREE] * GC_GDM;

		/* tree tC -> class tC */
		s->value[STEM_C]  = s->value[TREE_STEM_C] * s->counter[N_TREE];
	}
	else
	{
		logger(g_debug_log, "Ok stem biomass..\n");

		s->value[STEM_C]  = s->value[STEM_DM] / GC_GDM;
	}
	logger(g_debug_log, "---Stem Biomass = %f tC/tree\n",  s->value[TREE_STEM_C]);
	logger(g_debug_log, "---Stem Biomass = %f tDM/cell\n", s->value[STEM_DM]);
	logger(g_debug_log, "---Stem Biomass = %f tC/cell\n",  s->value[STEM_C]);

	if (s->value[BRANCH_DM] == 0.0 || s->value[BRANCH_DM] == NO_DATA)
	{
		logger(g_debug_log, "\nNo Branch and Bark Data are available from model initialization\n"
				"Is the Stem biomass initial value with Branch and Bark?\n");
		if (s->value[FRACBB0] == 0)
		{
			logger(g_debug_log, "I don't have FRACBB0 = FRACBB1 \n");

			s->value[FRACBB0] = s->value[FRACBB1];
			logger(g_debug_log, "FRACBB0 = %f\n", s->value[FRACBB0]);
		}
		else
		{
			s->value[FRACBB]    = s->value[FRACBB1] + (s->value[FRACBB0] - s->value[FRACBB1]) * exp(-LN2 * (h->value / s->value[TBB]));
			s->value[BRANCH_DM] = s->value[STEM_DM] * s->value[FRACBB];
			s->value[BRANCH_C]  = s->value[STEM_C] * s->value[FRACBB];
		}
	}
	else
	{
		logger(g_debug_log, "Ok stem branch biomass..\n");
		s->value[BRANCH_C] = s->value[BRANCH_DM] / GC_GDM;
	}

	s->value[TREE_BRANCH_C]      = s->value[BRANCH_C] / (double)s->counter[N_TREE];

	logger(g_debug_log, "---Branch Biomass = %f tC/tree\n",  s->value[TREE_BRANCH_C]);
	logger(g_debug_log, "---Branch Biomass = %f tDM/cell\n", s->value[BRANCH_DM]);
	logger(g_debug_log, "---Branch Biomass = %f tC/cell\n",  s->value[BRANCH_C]);


	if( s->value[CROOT_DM] == 0 || s->value[CROOT_DM] == NO_DATA)
	{
		logger(g_debug_log, "\nNo Coarse root Biomass Data are available for model initialization \n");
		logger(g_debug_log, "...Generating input Coarse root Biomass biomass data from DBH data...\n");

		/* compute coarse root biomass using root to shoot ratio */
		s->value[CROOT_DM] = s->value[STEM_DM] * s->value[COARSE_ROOT_STEM];
		s->value[CROOT_C]  = s->value[STEM_C] * s->value[COARSE_ROOT_STEM];
	}
	else
	{
		logger(g_debug_log, "Ok coarse root biomass..\n");
		s->value[CROOT_C] = s->value[CROOT_DM] / GC_GDM;
	}

	s->value[TREE_CROOT_C]    = s->value[CROOT_C] / (double)s->counter[N_TREE];
	logger(g_debug_log, "-Coarse root biomass = %f tC/tree\n",   s->value[TREE_CROOT_C]);
	logger(g_debug_log, "-Coarse root biomass = %f tDM/cell \n", s->value[CROOT_DM]);
	logger(g_debug_log, "-Coarse root biomass = %f tC/cell \n",  s->value[CROOT_C]);

	/* total wood */
	s->value[TOT_WOOD_C]      = s->value[STEM_C] + s->value[CROOT_C] + s->value[BRANCH_C];
	s->value[TREE_TOT_WOOD_C] = s->value[TREE_STEM_C] + s->value[TREE_CROOT_C] + s->value[TREE_BRANCH_C];

	logger(g_debug_log, "-Total wood = %f tC/cell\n",s->value[TOT_WOOD_C]);
	logger(g_debug_log, "-Total wood = %f tC/tree\n",s->value[TREE_TOT_WOOD_C]);

	/* basal area */
	s->value[BASAL_AREA]          = ( ( pow ( ( d->value / 2. ), 2. ) ) * Pi );
	s->value[BASAL_AREA_m2]       = s->value[BASAL_AREA] * 0.0001;
	s->value[STAND_BASAL_AREA]    = s->value[BASAL_AREA] * s->counter[N_TREE];
	s->value[STAND_BASAL_AREA_m2] = s->value[BASAL_AREA_m2] * s->counter[N_TREE];

	logger(g_debug_log, "-Basal area                       = %f cm2/tree\n",  s->value[BASAL_AREA]);
	logger(g_debug_log, "-Basal basal_area_m2              = %f m2/tree\n",   s->value[BASAL_AREA_m2]);
	logger(g_debug_log, "-Stand level class basal cell     = %f cm2/class\n", s->value[STAND_BASAL_AREA]);
	logger(g_debug_log, "-Stand level class basal cell     = %f cm2/class\n", s->value[STAND_BASAL_AREA]);

	/* sapwood and heartwood area are based on  Vertessy  et  al.  (1995)  and  Meinzer  et  al.(2001, 2005)*/
	s->value[SAPWOOD_AREA]           = s->value[SAP_A] * pow ( d->value, s->value[SAP_B] );

	//fixme special case
	if( s->value[SAPWOOD_AREA] > s->value[BASAL_AREA] ) s->value[SAPWOOD_AREA] = s->value[BASAL_AREA];

	s->value[HEARTWOOD_AREA]         = s->value[BASAL_AREA] -  s->value[SAPWOOD_AREA];
	s->value[SAPWOOD_PERC]           = s->value[SAPWOOD_AREA] / s->value[BASAL_AREA];

	logger(g_debug_log, "-Sapwood_area                     = %f cm^2\n",       s->value[SAPWOOD_AREA]);
	logger(g_debug_log, "-Heartwood_area                   = %f cm^2\n",       s->value[HEARTWOOD_AREA]);
	logger(g_debug_log, "-Sapwood perc                     = %f %%\n",         s->value[SAPWOOD_PERC] * 100. );

	/* sapwood and heartwood biomass */
	/* stem */
	s->value[STEM_SAPWOOD_DM]         = s->value[STEM_DM] * s->value[SAPWOOD_PERC];    // FIXME sapwood_perc is simply the ration of two areas,
	                                                                                   // not really correct to apply then the relationship to mass
	                                                                                   //  (and hence volume)
	s->value[STEM_SAPWOOD_C]          = s->value[STEM_C]  * s->value[SAPWOOD_PERC];
	s->value[STEM_HEARTWOOD_DM]       = s->value[STEM_DM] - s->value[STEM_SAPWOOD_DM];
	s->value[STEM_HEARTWOOD_C]        = s->value[STEM_C]  - s->value[STEM_SAPWOOD_C];
	s->value[TREE_STEM_SAPWOOD_C]     = s->value[STEM_SAPWOOD_C] / (double)s->counter[N_TREE];
        s->value[TREE_STEM_HEARTWOOD_C]   = (s->value[STEM_C] - s->value[STEM_SAPWOOD_C]) / (double)s->counter[N_TREE];     

	/* coarse root */
	s->value[CROOT_SAPWOOD_DM]        = s->value[CROOT_DM] * s->value[SAPWOOD_PERC];
	s->value[CROOT_SAPWOOD_C]         = s->value[CROOT_C]  * s->value[SAPWOOD_PERC];
	s->value[CROOT_HEARTWOOD_DM]      = s->value[CROOT_DM] - s->value[CROOT_SAPWOOD_DM];
	s->value[CROOT_HEARTWOOD_C]       = s->value[CROOT_C]  - s->value[CROOT_SAPWOOD_C];
	s->value[TREE_CROOT_SAPWOOD_C]    = s->value[CROOT_SAPWOOD_C] / (double)s->counter[N_TREE];
        s->value[TREE_CROOT_HEARTWOOD_C]  = (s->value[CROOT_C] - s->value[CROOT_SAPWOOD_C] ) / (double)s->counter[N_TREE];

	/* branch */
	s->value[BRANCH_SAPWOOD_DM]       = s->value[BRANCH_DM] * s->value[SAPWOOD_PERC];
	s->value[BRANCH_SAPWOOD_C]        = s->value[BRANCH_C]  * s->value[SAPWOOD_PERC];
	s->value[BRANCH_HEARTWOOD_DM]     = s->value[BRANCH_DM] - s->value[BRANCH_SAPWOOD_DM];
	s->value[BRANCH_HEARTWOOD_C]      = s->value[BRANCH_C]  - s->value[BRANCH_SAPWOOD_C];
	s->value[TREE_BRANCH_SAPWOOD_C]   = s->value[BRANCH_SAPWOOD_C] / (double)s->counter[N_TREE];
        s->value[TREE_BRANCH_HEARTWOOD_C] = (s->value[BRANCH_C] - s->value[BRANCH_SAPWOOD_C]) / (double)s->counter[N_TREE];
 
	/* overall */
	s->value[TOT_SAPWOOD_DM]          = s->value[STEM_SAPWOOD_DM]   + s->value[CROOT_SAPWOOD_DM]   + s->value[BRANCH_SAPWOOD_DM];
	s->value[TOT_SAPWOOD_C]           = s->value[STEM_SAPWOOD_C]    + s->value[CROOT_SAPWOOD_C]    + s->value[BRANCH_SAPWOOD_C];
	s->value[TOT_HEARTWOOD_DM]        = s->value[STEM_HEARTWOOD_DM] + s->value[CROOT_HEARTWOOD_DM] + s->value[BRANCH_HEARTWOOD_DM];
	s->value[TOT_HEARTWOOD_C]         = s->value[STEM_HEARTWOOD_C]  + s->value[CROOT_HEARTWOOD_C]  + s->value[BRANCH_HEARTWOOD_C];
        s->value[TREE_SAPWOOD_C]          = s->value[TOT_SAPWOOD_C] / (double)s->counter[N_TREE];
        s->value[TREE_HEARTWOOD_C]        = (s->value[TOT_WOOD_C] - s->value[TOT_SAPWOOD_C]) / (double)s->counter[N_TREE];

	logger(g_debug_log, "-Sapwood stem biomass             = %f tC/cell\n",  s->value[STEM_SAPWOOD_C]);
	logger(g_debug_log, "-Heartwood stem biomass           = %f tC/cell\n",  s->value[STEM_HEARTWOOD_C]);
	logger(g_debug_log, "-Sapwood coarse root biomass      = %f tC/cell\n",  s->value[CROOT_SAPWOOD_C]);
	logger(g_debug_log, "-Heartwood coarse root biomass    = %f tC/cell\n",  s->value[CROOT_HEARTWOOD_C]);
	logger(g_debug_log, "-Sapwood branch and bark biomass  = %f tC/cell\n",  s->value[BRANCH_SAPWOOD_C]);
	logger(g_debug_log, "-Heartwood branch biomass         = %f tC/cell\n",  s->value[BRANCH_HEARTWOOD_C]);
	logger(g_debug_log, "-Total Sapwood biomass            = %f tC/cell\n",  s->value[TOT_SAPWOOD_C]);
	logger(g_debug_log, "-Total Sapwood biomass per tree   = %f tC/tree\n",  s->value[TOT_SAPWOOD_C]    / (double)s->counter[N_TREE]);
	logger(g_debug_log, "-Total Sapwood biomass per tree   = %f KgC/tree\n", (s->value[TOT_SAPWOOD_C]   / (double)s->counter[N_TREE]) * 1e3);
	logger(g_debug_log, "-Total Heartwood biomass          = %f tC/cell\n",  s->value[TOT_HEARTWOOD_C]);
	logger(g_debug_log, "-Total Heartwood biomass per tree = %f tC/tree\n",  s->value[TOT_HEARTWOOD_C]  / (double)s->counter[N_TREE]);
	logger(g_debug_log, "-Total Heartwood biomass per tree = %f KgC/tree\n", (s->value[TOT_HEARTWOOD_C] / (double)s->counter[N_TREE]) * 1e3);

	/* reserve */
	if ( s->value[RESERVE_DM] == 0 || s->value[RESERVE_DM] == NO_DATA )
	{
		logger(g_debug_log, "\nNo Reserve Biomass Data are available for model initialization \n");
		logger(g_debug_log, "...Generating input Reserve Biomass biomass data\n");

		/* note: these values are taken from: Schwalm and Ek, 2004 Ecological Modelling */
		// see if change with the ratio reported from Barbaroux et al., 2002 (using DryMatter)
		/* NSC are initialized at the minimum annual value at the beginning of simulation (1. january hence winter)
		   see Woodruff and Meinzer 2011, Plant Cell and Environment */

		/* IMPORTANT! reserve computation if not in initialized is computed from DryMatter */
		s->value[RESERVE_DM] = s->value[TOT_SAPWOOD_DM] * s->value[SAP_WRES];

		s->value[RESERVE_C]  = s->value[TOT_SAPWOOD_DM] * s->value[SAP_WRES];
	}

	s->value[TREE_RESERVE_C]    = s->value[RESERVE_C] / (double)s->counter[N_TREE];

	/* set minimum reserve pool */
	s->value[MIN_RESERVE_C]      = s->value[RESERVE_C];
	s->value[TREE_MIN_RESERVE_C] = s->value[MIN_RESERVE_C] / (double)s->counter[N_TREE];

	/*************************************************************************************************************************************************/

	/* LEAF */
	if ( ! s->value[LEAF_DM] || s->value[LEAF_DM] == NO_DATA )
	{
		/* deciduous */
		if ( s->value[PHENOLOGY] == 0.1 || s->value[PHENOLOGY] == 0.2 )
		{
			/* assuming no leaf at 1st of January */
			s->value[LEAF_DM]        = 0.;
			s->value[LEAF_C]         = 0.;
			s->value[LEAF_SUN_C]     = 0.;
			s->value[LEAF_SHADE_C]   = 0.;
			s->value[LAI_PROJ]       = 0.;
			s->value[LAI_SUN_PROJ]   = 0.;
			s->value[LAI_SHADE_PROJ] = 0.;
		}
		/* evergreen */
		else
		{
			/* compute maximum LAI at peak value */
			s->value[PEAK_LAI_PROJ]  = ( ( s->value[SAPWOOD_AREA] / 1e4 ) * s->value[SAP_LEAF]) / s->value[CROWN_AREA_PROJ];

			/* compute leaf carbon to PEAK LAI down-scaled to canopy cover */
			s->value[LEAF_C]         = ( ( s->value[PEAK_LAI_PROJ] / s->value[SLA_AVG] ) / 1e3 ) * ( s->value[CANOPY_COVER_PROJ] * g_settings->sizeCell );

			/* Calculate projected LAI */
			/* note: assume that at the beginning of simulations LAI = PEAK_LAI */
			s->value[LAI_PROJ]       = s->value[PEAK_LAI_PROJ];

			/* Calculate projected LAI for sun and shaded leaves */
			/* sun */
			s->value[LAI_SUN_PROJ]   = 1. - exp( - s->value[LAI_PROJ] );

			/* shade */
			s->value[LAI_SHADE_PROJ] = s->value[LAI_PROJ] - s->value[LAI_SUN_PROJ];
		}
	}
	else
	{
		/* if no values for LAI are available */
		if ( ! s->value[LAI_PROJ] && ( s->value[PHENOLOGY] == 1.1 || s->value[PHENOLOGY] == 1.2 ) )
		{
			/* Calculate projected LAI */
			s->value[LAI_PROJ]       = ( ( s->value[LEAF_C]  * s->value[SLA_AVG] ) * 1e3 ) / ( s->value[CANOPY_COVER_PROJ] * g_settings->sizeCell );

			/* Calculate projected LAI for sun and shaded leaves */
			/* sun */
			s->value[LAI_SUN_PROJ]   = 1. - exp( - s->value[LAI_PROJ] );

			/* shade */
			s->value[LAI_SHADE_PROJ] = s->value[LAI_PROJ] - s->value[LAI_SUN_PROJ];
		}
		else
		{
			s->value[LAI_PROJ]       = 0.;
			s->value[LAI_SUN_PROJ]   = 0.;
			s->value[LAI_SHADE_PROJ] = 0.;
		}
	}

	/*************************************************************************************/

	/* compute specific leaf area projected for sun and shaded (m2/kgC) */
	if ( s->value[LEAF_C] > 0. )
	{
		/* for total leaves */
		s->value[SLA_PROJ]           =  ( s->value[LAI_PROJ] * ( s->value[CANOPY_COVER_PROJ] * g_settings->sizeCell ) ) / ( s->value[LEAF_C] * 1e3 ) ;

		/* sun leaves */
		s->value[SLA_SUN_PROJ]       = ( ( s->value[LAI_SUN_PROJ] + ( s->value[LAI_SHADE_PROJ] / s->value[SLA_RATIO] ) ) * ( s->value[CANOPY_COVER_PROJ] * g_settings->sizeCell ) ) / ( s->value[LEAF_C] * 1e3 );

		/* shaded leaves */
		s->value[SLA_SHADE_PROJ]     = s->value[SLA_SUN_PROJ] * s->value[SLA_RATIO];

		/* compute based on SLA Leaf carbon for sun and shaded leaves */
		s->value[LEAF_SUN_C]   = ( ( s->value[LAI_SUN_PROJ]   * ( s->value[CANOPY_COVER_PROJ] * g_settings->sizeCell ) ) / s->value[SLA_SUN_PROJ] )   / 1e3;
		s->value[LEAF_SHADE_C] = s->value[LEAF_C] - s->value[LEAF_SUN_C];

		/* check */
		CHECK_CONDITION  (fabs( ( s->value[LEAF_SUN_C] + s->value[LEAF_SHADE_C] ) - s->value[LEAF_C]) , > , eps);
	}

	/*************************************************************************************/

	/* convert tC/cell to tCDM/cell */
	s->value[LEAF_DM]        = s->value[LEAF_DM] * GC_GDM;

	/* compute single tree leaf carbon amount */
	s->value[TREE_LEAF_C] = s->value[LEAF_C] / s->counter[N_TREE];

	/* compute all-sided Leaf Area */
	s->value[ALL_LAI_PROJ] = s->value[LAI_PROJ] * s->value[CANOPY_COVER_PROJ];

	/*************************************************************************************************************************************************/

	/* note: model assumes that if no fine-root biomass are available the same ratio foliage-fine roots is used */
	if (( s->value[FROOT_DM] == 0.0 || s->value[FROOT_DM] == NO_DATA)
			&& (s->value[PHENOLOGY] == 1.1 || s->value[PHENOLOGY] == 1.2))
	{
		s->value[FROOT_DM] = s->value[LEAF_DM] * s->value[FINE_ROOT_LEAF];
		s->value[FROOT_C]  = s->value[LEAF_C] * s->value[FINE_ROOT_LEAF];
	}
	else
	{
		/* assuming no fine root at 1st of January */
		s->value[FROOT_DM]     = 0.;
		s->value[FROOT_C]      = 0.;
	}

	s->value[TREE_FROOT_C]    = s->value[FROOT_C] / (double)s->counter[N_TREE];
 
	/***** INITIALIZE LITTER POOL *****/


	if ( ! c->init_litter_C || c->init_litter_C == NO_DATA )
	{
		/* compute leaf litter (assuming that at year zero litter is composed by the amount of peak lai of the previous year */
		s->value[LEAF_LITRC]  = ( ( s->value[PEAK_LAI_PROJ] / s->value[SLA_AVG] ) / 1e3 * ( s->value[CANOPY_COVER_PROJ] * g_settings->sizeCell ) ) / GC_GDM;
		s->value[FROOT_LITRC] = s->value[LEAF_LITRC] * ( 1. - s->value[FROOT_LEAF_FRAC] );

		if ( s->value[PHENOLOGY] == 1.1 || s->value[PHENOLOGY] == 1.2 )
		{
			/* note: if evergreen we assume that only a fraction of previous year leaf carbon goes to litter */
			/*(based on species leaf fall turnover rate) (not simulated) has gone to the litter carbon pool */
			s->value[LEAF_LITRC]  *= s->value[LEAF_FROOT_TURNOVER];
			s->value[FROOT_LITRC] *= s->value[LEAF_FROOT_TURNOVER];
		}
	}
	else
	{
		/* get data from soil setting file */
		s->value[LEAF_LITRC]  = c->init_litter_C * s->value[FROOT_LEAF_FRAC];
		s->value[FROOT_LITRC] = c->init_litter_C * ( 1. - s->value[FINE_ROOT_LEAF]);
	}

	/***** COMPUTE LIVE DEAD BIOMASS *****/

	live_total_wood_age( a, s );

	/* stem */
#if TEST_RESP
	//old
	s->value[STEM_LIVEWOOD_DM]     = s->value[STEM_DM] * s->value[EFF_LIVE_TOTAL_WOOD_FRAC];
	s->value[STEM_LIVEWOOD_C]      = s->value[STEM_C]  * s->value[EFF_LIVE_TOTAL_WOOD_FRAC];
	s->value[STEM_DEADWOOD_DM]     = s->value[STEM_DM] - s->value[STEM_LIVEWOOD_DM];
	s->value[STEM_DEADWOOD_C]      = s->value[STEM_C]  - s->value[STEM_LIVEWOOD_C];

	s->value[TREE_STEM_LIVE_KgDM]  = s->value[STEM_LIVEWOOD_DM] * 1e3 / (double)s->counter[N_TREE];
	s->value[TREE_STEM_LIVE_C]     = s->value[STEM_LIVEWOOD_C] / (double)s->counter[N_TREE];
	s->value[TREE_STEM_DEAD_KgDM]  = s->value[STEM_DEADWOOD_DM] * 1e3 / s->counter[N_TREE];
	s->value[TREE_STEM_DEAD_C]     = s->value[STEM_DEADWOOD_C] / (double)s->counter[N_TREE];

	/* check */
	CHECK_CONDITION(fabs((s->value[STEM_C])-(s->value[STEM_LIVEWOOD_C] + s->value[STEM_DEADWOOD_C])), >,eps);

#else
	//new 13 May 2017
	/* class */
	s->value[STEM_LIVEWOOD_DM]     = s->value[STEM_SAPWOOD_DM] * s->value[EFF_LIVE_TOTAL_WOOD_FRAC];
	s->value[STEM_LIVEWOOD_C]      = s->value[STEM_SAPWOOD_C]  * s->value[EFF_LIVE_TOTAL_WOOD_FRAC];
	s->value[STEM_DEADWOOD_DM]     = s->value[STEM_DM]         - s->value[STEM_LIVEWOOD_DM];
	s->value[STEM_DEADWOOD_C]      = s->value[STEM_C]          - s->value[STEM_LIVEWOOD_C];

	/* tree */
	s->value[TREE_STEM_LIVEWOOD_C] = s->value[STEM_LIVEWOOD_C] / (double)s->counter[N_TREE];
	s->value[TREE_STEM_DEADWOOD_C] = s->value[STEM_DEADWOOD_C] / (double)s->counter[N_TREE];

	/* check */
	CHECK_CONDITION(fabs((s->value[STEM_DM])-(s->value[STEM_LIVEWOOD_DM] + s->value[STEM_DEADWOOD_DM])), >,eps);
	CHECK_CONDITION(fabs((s->value[STEM_C]) -(s->value[STEM_LIVEWOOD_C]  + s->value[STEM_DEADWOOD_C])),  >,eps);

#endif


#if TEST_RESP
	//old
	/* coarse root */
	s->value[CROOT_LIVEWOOD_DM]     = s->value[CROOT_DM] * (s->value[EFF_LIVE_TOTAL_WOOD_FRAC]);
	s->value[CROOT_LIVEWOOD_C]      = s->value[CROOT_C] * (s->value[EFF_LIVE_TOTAL_WOOD_FRAC]);
	s->value[CROOT_DEADWOOD_DM]     = s->value[CROOT_DM] -s->value[CROOT_LIVEWOOD_DM];
	s->value[CROOT_DEADWOOD_C]      = s->value[CROOT_C] -s->value[CROOT_LIVEWOOD_C];

	s->value[TREE_CROOT_LIVE_KgDM]  = s->value[CROOT_LIVEWOOD_DM] * 1e3 / (double)s->counter[N_TREE];
	s->value[TREE_CROOT_LIVE_C]     = s->value[CROOT_LIVEWOOD_C] / (double)s->counter[N_TREE];
	s->value[TREE_CROOT_DEAD_KgDM]  = s->value[CROOT_DEADWOOD_DM] * 1e3 / (double)s->counter[N_TREE];
	s->value[TREE_CROOT_DEAD_C]     = s->value[CROOT_DEADWOOD_C] / (double)s->counter[N_TREE];

	/* check */
	CHECK_CONDITION(fabs((s->value[CROOT_C])-(s->value[CROOT_LIVEWOOD_C] + s->value[CROOT_DEADWOOD_C])), >,eps);

#else
	//new 13 May 2017
	/* class */
	s->value[CROOT_LIVEWOOD_DM]     = s->value[CROOT_SAPWOOD_DM] * s->value[EFF_LIVE_TOTAL_WOOD_FRAC];
	s->value[CROOT_LIVEWOOD_C]      = s->value[CROOT_SAPWOOD_C]  * s->value[EFF_LIVE_TOTAL_WOOD_FRAC];
	s->value[CROOT_DEADWOOD_DM]     = s->value[CROOT_DM]         - s->value[CROOT_LIVEWOOD_DM];
	s->value[CROOT_DEADWOOD_C]      = s->value[CROOT_C]          - s->value[CROOT_LIVEWOOD_C];

	/* tree */
	s->value[TREE_CROOT_LIVEWOOD_C] = s->value[CROOT_LIVEWOOD_C] / (double)s->counter[N_TREE];
	s->value[TREE_CROOT_DEADWOOD_C] = s->value[CROOT_DEADWOOD_C] / (double)s->counter[N_TREE];

	/* check */
	CHECK_CONDITION(fabs((s->value[CROOT_DM])-(s->value[CROOT_LIVEWOOD_DM] + s->value[CROOT_DEADWOOD_DM])), >,eps);
	CHECK_CONDITION(fabs((s->value[CROOT_C]) -(s->value[CROOT_LIVEWOOD_C]  + s->value[CROOT_DEADWOOD_C])),  >,eps);

#endif

	logger(g_debug_log, "-Coarse root Biomass      = %f tC/cell\n", s->value[CROOT_C]);
	logger(g_debug_log, "-Live Coarse Root Biomass = %f tC/cell\n", s->value[CROOT_LIVEWOOD_C]);
	logger(g_debug_log, "-Dead Coarse Root Biomass = %f tC/cell\n", s->value[CROOT_DEADWOOD_C]);
	logger(g_debug_log, "-Live Coarse Root Biomass = %f tC/tree\n", s->value[TREE_CROOT_LIVEWOOD_C]);
	logger(g_debug_log, "-Dead Coarse Root Biomass = %f tC/tree\n", s->value[TREE_CROOT_DEADWOOD_C]);

#if TEST_RESP
	//old
	/* branch */
	s->value[BRANCH_LIVEWOOD_DM]      = s->value[BRANCH_DM] * (s->value[EFF_LIVE_TOTAL_WOOD_FRAC]);
	s->value[BRANCH_LIVEWOOD_C]       = s->value[BRANCH_C]  * (s->value[EFF_LIVE_TOTAL_WOOD_FRAC]);
	s->value[BRANCH_DEADWOOD_DM]      = s->value[BRANCH_DM] - s->value[BRANCH_LIVEWOOD_DM];
	s->value[BRANCH_DEADWOOD_C]       = s->value[BRANCH_C]  - s->value[BRANCH_LIVEWOOD_C];

	s->value[TREE_BRANCH_LIVE_KgDM]   = s->value[BRANCH_LIVEWOOD_DM] * 1e3 / (double)s->counter[N_TREE];
	s->value[TREE_BRANCH_LIVE_C]      = s->value[BRANCH_LIVEWOOD_C] / s->counter[N_TREE];
	s->value[TREE_BRANCH_DEAD_KgDM]   = s->value[BRANCH_DEADWOOD_DM] * 1e3 / (double)s->counter[N_TREE];
	s->value[TREE_BRANCH_DEAD_C]      = s->value[BRANCH_DEADWOOD_C] / (double)s->counter[N_TREE];

	/* check */
	CHECK_CONDITION(fabs((s->value[BRANCH_C])-(s->value[BRANCH_LIVEWOOD_C] + s->value[BRANCH_DEADWOOD_C])), >,eps);

#else
	//new 13 May 2017
	/* class */
	s->value[BRANCH_LIVEWOOD_DM]      = s->value[BRANCH_SAPWOOD_DM] * s->value[EFF_LIVE_TOTAL_WOOD_FRAC];
	s->value[BRANCH_LIVEWOOD_C]       = s->value[BRANCH_SAPWOOD_C]  * s->value[EFF_LIVE_TOTAL_WOOD_FRAC];
	s->value[BRANCH_DEADWOOD_DM]      = s->value[BRANCH_DM]         - s->value[BRANCH_LIVEWOOD_DM];
	s->value[BRANCH_DEADWOOD_C]       = s->value[BRANCH_C]          - s->value[BRANCH_LIVEWOOD_C];

	/* tree */
	s->value[TREE_BRANCH_LIVEWOOD_C]  = s->value[BRANCH_LIVEWOOD_C] / (double)s->counter[N_TREE];
	s->value[TREE_BRANCH_DEADWOOD_C]  = s->value[BRANCH_DEADWOOD_C] / (double)s->counter[N_TREE];

	/* check */
	CHECK_CONDITION(fabs((s->value[BRANCH_DM])-(s->value[BRANCH_LIVEWOOD_DM] + s->value[BRANCH_DEADWOOD_DM])), >,eps);
	CHECK_CONDITION(fabs((s->value[BRANCH_C]) -(s->value[BRANCH_LIVEWOOD_C]  + s->value[BRANCH_DEADWOOD_C])),  >,eps);

#endif

	logger(g_debug_log, "-Branch Biuomass          = %f tC/cell\n", s->value[BRANCH_C]);
	logger(g_debug_log, "-Live Stem Branch Biomass = %f tC/cell\n", s->value[BRANCH_LIVEWOOD_C]);
	logger(g_debug_log, "-Dead Stem Branch Biomass = %f tC/cell\n", s->value[BRANCH_DEADWOOD_C]);
	logger(g_debug_log, "-Live Stem Branch Biomass = %f tC/tree\n", s->value[TREE_BRANCH_LIVEWOOD_C]);
	logger(g_debug_log, "-Dead Stem Branch Biomass = %f tC/tree\n", s->value[TREE_BRANCH_DEADWOOD_C]);

	/** live biomass **/
	/* class */
	s->value[TOT_LIVEWOOD_DM]       = s->value[STEM_LIVEWOOD_DM]  + s->value[CROOT_LIVEWOOD_DM] + s->value[BRANCH_LIVEWOOD_DM];
	s->value[TOT_LIVEWOOD_C]        = s->value[STEM_LIVEWOOD_C]   + s->value[CROOT_LIVEWOOD_C]  + s->value[BRANCH_LIVEWOOD_C];

	/* tree */
	s->value[TREE_TOT_LIVEWOOD_C]   = s->value[TOT_LIVEWOOD_C] / (double)s->counter[N_TREE];

	logger(g_debug_log, "-Live biomass = %f tC/cell\n", s->value[TOT_LIVEWOOD_C]);
	logger(g_debug_log, "-Live biomass = %f tC/tree\n", s->value[TREE_TOT_LIVEWOOD_C]);

	/** dead biomass **/
	/* class */
	s->value[TOT_DEADWOOD_DM]        = s->value[STEM_DEADWOOD_DM] + s->value[CROOT_DEADWOOD_DM] + s->value[BRANCH_DEADWOOD_DM];
	s->value[TOT_DEADWOOD_C]         = s->value[STEM_DEADWOOD_C]  + s->value[CROOT_DEADWOOD_C]  + s->value[BRANCH_DEADWOOD_C];

	/* tree */
	s->value[TREE_TOT_DEADWOOD_C]    = s->value[TOT_DEADWOOD_C] / (double)s->counter[N_TREE];

	logger(g_debug_log, "-Dead biomass = %f tC/cell\n", s->value[TOT_DEADWOOD_C]);
	logger(g_debug_log, "-Dead biomass = %f tC/tree\n", s->value[TREE_TOT_DEADWOOD_C]);

	logger(g_debug_log, "-Live wood vs total biomass = %f %%\n", (s->value[TOT_LIVEWOOD_C] / s->value[TOT_WOOD_C]) * 100.);
	logger(g_debug_log, "-Dead wood vs total biomass = %f %%\n", (s->value[TOT_DEADWOOD_C] / s->value[TOT_WOOD_C]) * 100.);
	logger(g_debug_log, "-Reserve                    = %f tC/cell\n", s->value[TOT_LIVEWOOD_C] * s->value[SAP_WRES] );

	/* fruit */
	s->value[FRUIT_C]      = 0.;
	s->value[TREE_FRUIT_C] = 0.;

	/* compute AGB and BGB */
	s->value[AGB]      = s->value[STEM_C] + s->value[BRANCH_C] + s->value[LEAF_C] + s->value[FRUIT_C];
	s->value[BGB]      = s->value[FROOT_C] + s->value[CROOT_C];
	s->value[TREE_AGB] = s->value[AGB] / (double)s->counter[N_TREE];
	s->value[TREE_BGB] = s->value[BGB] / (double)s->counter[N_TREE];

	logger(g_debug_log, "-Yearly Class AGB = %f tC/cell\n", s->value[AGB]);
	logger(g_debug_log, "-Yearly Class BGB = %f tC/cell\n", s->value[BGB]);
	logger(g_debug_log, "-Yearly Class AGB = %f tC/tree\n", s->value[TREE_AGB]);
	logger(g_debug_log, "-Yearly Class BGB = %f tC/tree\n", s->value[TREE_BGB]);

	/* check */
	CHECK_CONDITION ( fabs (s->value[LEAF_LITR_LAB_FRAC]  + s->value[LEAF_LITR_CEL_FRAC]  + s->value[LEAF_LITR_LIGN_FRAC]),  >, 1 + eps);
	CHECK_CONDITION ( fabs (s->value[FROOT_LITR_LAB_FRAC] + s->value[FROOT_LITR_CEL_FRAC] + s->value[FROOT_LITR_LIGN_FRAC]), >, 1 + eps);
	CHECK_CONDITION ( fabs (s->value[DEADWOOD_CEL_FRAC]   + s->value[DEADWOOD_LIGN_FRAC]), >, 1 + eps);

	/******************************************************************************************************************************************************************/
	/*** compute coarse woody debris litter fractions ***/

	/* partitioning coarse woody debris litter carbon cellulose into shielded and unshielded pools */
	r1 = s->value[DEADWOOD_LIGN_FRAC] / s->value[DEADWOOD_CEL_FRAC];

	if ( r1 <= 0.45 )
	{
		s->value[DEADWOOD_SCEL_FRAC]  = 0.;
		s->value[DEADWOOD_USCEL_FRAC] = s->value[DEADWOOD_CEL_FRAC];
	}
	else if ( r1 > 0.45 && r1 < 0.7 )
	{
		temp_var = ( r1 - 0.45 ) * 3.2;
		s->value[DEADWOOD_SCEL_FRAC]  = temp_var * s->value[DEADWOOD_CEL_FRAC];
		s->value[DEADWOOD_USCEL_FRAC] = (1. - temp_var) * s->value[DEADWOOD_CEL_FRAC];
	}
	else
	{
		s->value[DEADWOOD_SCEL_FRAC]  = 0.8 * s->value[DEADWOOD_CEL_FRAC];
		s->value[DEADWOOD_USCEL_FRAC] = 0.2 * s->value[DEADWOOD_CEL_FRAC];
	}

	/******************************************************************************************************************************************************************/
	/** calculate shielded and unshielded cellulose fraction for leaf pool **/

	/*** compute leaf litter fractions ***/

	/* partitioning leaf litter carbon cellulose into shielded and unshielded pools */
	r1 = s->value[LEAF_LITR_LIGN_FRAC] / s->value[LEAF_LITR_LAB_FRAC];

	if ( r1 <= 0.45 )
	{
		s->value[LEAF_LITR_SCEL_FRAC]   = 0.;
		s->value[LEAF_LITR_USCEL_FRAC]  = s->value[LEAF_LITR_CEL_FRAC];
	}
	else if ( r1 > 0.45 && r1 < 0.7 )
	{
		temp_var = ( r1 - 0.45 ) * 3.2;
		s->value[LEAF_LITR_SCEL_FRAC]   = temp_var * s->value[LEAF_LITR_CEL_FRAC];
		s->value[LEAF_LITR_USCEL_FRAC]  = ( 1. - temp_var ) * s->value[LEAF_LITR_CEL_FRAC];
	}
	else
	{
		s->value[LEAF_LITR_SCEL_FRAC]   = 0.8 * s->value[LEAF_LITR_CEL_FRAC];
		s->value[LEAF_LITR_USCEL_FRAC]  = 0.2 * s->value[LEAF_LITR_CEL_FRAC];
	}
	/* check */
	CHECK_CONDITION ( fabs ( s->value[LEAF_LITR_LAB_FRAC] + s->value[LEAF_LITR_SCEL_FRAC] + s->value[LEAF_LITR_USCEL_FRAC] + s->value[LEAF_LITR_LIGN_FRAC] ), >, 1 + eps);

	/******************************************************************************************************************************************************************/
	/** calculate shielded and unshielded cellulose fraction for froot pool **/

	/*** compute fine root litter fractions ***/

	/* partitioning fine root litter carbon cellulose into shielded and unshielded pools */
	r1 = s->value[FROOT_LITR_LIGN_FRAC] / s->value[FROOT_LITR_LAB_FRAC];

	if ( r1 <= 0.45 )
	{
		s->value[FROOT_LITR_SCEL_FRAC]   = 0.;
		s->value[FROOT_LITR_USCEL_FRAC]  = s->value[FROOT_LITR_CEL_FRAC];
	}
	else if ( r1 > 0.45 && r1 < 0.7 )
	{
		temp_var = ( r1 - 0.45 ) * 3.2;
		s->value[FROOT_LITR_SCEL_FRAC]   = temp_var * s->value[FROOT_LITR_CEL_FRAC];
		s->value[FROOT_LITR_USCEL_FRAC]  = ( 1. - temp_var ) * s->value[FROOT_LITR_CEL_FRAC];
	}
	else
	{
		s->value[FROOT_LITR_SCEL_FRAC]   = 0.8 * s->value[FROOT_LITR_CEL_FRAC];
		s->value[FROOT_LITR_USCEL_FRAC]  = 0.2 * s->value[FROOT_LITR_CEL_FRAC];
	}
	/* check */
	CHECK_CONDITION ( fabs ( s->value[FROOT_LITR_LAB_FRAC] + s->value[FROOT_LITR_SCEL_FRAC] + s->value[FROOT_LITR_USCEL_FRAC] + s->value[FROOT_LITR_LIGN_FRAC] ), >, 1 + eps);

	/******************************************************************************************************************************************************************/

	/* ONLY for evergreen */
	if ( s->value[PHENOLOGY] == 1.1 || s->value[PHENOLOGY] == 1.2 )
	{
		CHECK_CONDITION(s->value[LEAF_C],              <=, ZERO);
		CHECK_CONDITION(s->value[FROOT_C],             <=, ZERO);
		CHECK_CONDITION(s->value[LAI_PROJ],            <=, ZERO);
		CHECK_CONDITION(s->value[LAI_SUN_PROJ],        <=, ZERO);
		CHECK_CONDITION(s->value[LAI_SHADE_PROJ],      <=, ZERO);
	}
	CHECK_CONDITION(s->value[STEM_C],                  <=, ZERO);
	CHECK_CONDITION(s->value[BRANCH_C],                <=, ZERO);
	CHECK_CONDITION(s->value[CROOT_C],                 <=, ZERO);
	CHECK_CONDITION(s->value[TOT_WOOD_C],              <=, ZERO);
	CHECK_CONDITION(s->value[RESERVE_C],               <=, ZERO);
	CHECK_CONDITION(s->value[MIN_RESERVE_C],           <=, ZERO);
	CHECK_CONDITION(s->value[STEM_SAPWOOD_C],          <=, ZERO);
	CHECK_CONDITION(s->value[CROOT_SAPWOOD_C],         <=, ZERO);
	CHECK_CONDITION(s->value[BRANCH_SAPWOOD_C],        <=, ZERO);
	CHECK_CONDITION(s->value[STEM_HEARTWOOD_C],        <=, ZERO);
	CHECK_CONDITION(s->value[CROOT_HEARTWOOD_C],       <=, ZERO);
	CHECK_CONDITION(s->value[BRANCH_HEARTWOOD_C],      <=, ZERO);
	CHECK_CONDITION(s->value[TOT_SAPWOOD_C],           <=, ZERO);
	CHECK_CONDITION(s->value[STEM_LIVEWOOD_C],         <=, ZERO);
	CHECK_CONDITION(s->value[STEM_DEADWOOD_C],         <=, ZERO);
	CHECK_CONDITION(s->value[CROOT_LIVEWOOD_C],        <=, ZERO);
	CHECK_CONDITION(s->value[CROOT_DEADWOOD_C],        <=, ZERO);
	CHECK_CONDITION(s->value[BRANCH_LIVEWOOD_C],       <=, ZERO);
	CHECK_CONDITION(s->value[BRANCH_DEADWOOD_C],       <=, ZERO);
	CHECK_CONDITION(s->value[TOT_LIVEWOOD_C],          <=, ZERO);
	CHECK_CONDITION(s->value[TOT_DEADWOOD_C],          <=, ZERO);
	CHECK_CONDITION(s->value[BASAL_AREA],              <=, ZERO);
	CHECK_CONDITION(s->value[AGB],                     <=, ZERO);
	CHECK_CONDITION(s->value[BGB],                     <=, ZERO);
	CHECK_CONDITION(s->value[VOLUME],                  <=, ZERO);
	CHECK_CONDITION(s->value[TREE_VOLUME],             <=, ZERO);  
}

void initialization_forest_cell_C (cell_t *const c, const int height, const int dbh, const int age, const int species)
{
	species_t *s;
	s = &c->heights[height].dbhs[dbh].ages[age].species[species];

	/*** update at cell level (gC/m2)  ***/
	c->leaf_carbon              += (s->value[LEAF_C]             * 1e6 / g_settings->sizeCell);
	c->froot_carbon             += (s->value[FROOT_C]            * 1e6 / g_settings->sizeCell);
	c->stem_carbon              += (s->value[STEM_C]             * 1e6 / g_settings->sizeCell);
	c->stem_live_wood_carbon    += (s->value[STEM_LIVEWOOD_C]    * 1e6 / g_settings->sizeCell);
	c->stem_dead_wood_carbon    += (s->value[STEM_DEADWOOD_C]    * 1e6 / g_settings->sizeCell);
	c->croot_carbon             += (s->value[CROOT_C]            * 1e6 / g_settings->sizeCell);
	c->croot_live_wood_carbon   += (s->value[CROOT_LIVEWOOD_C]   * 1e6 / g_settings->sizeCell);
	c->croot_dead_wood_carbon   += (s->value[CROOT_DEADWOOD_C]   * 1e6 / g_settings->sizeCell);
	c->branch_carbon            += (s->value[BRANCH_C]           * 1e6 / g_settings->sizeCell);
	c->branch_live_wood_carbon  += (s->value[BRANCH_LIVEWOOD_C]  * 1e6 / g_settings->sizeCell);
	c->branch_dead_wood_carbon  += (s->value[BRANCH_DEADWOOD_C]  * 1e6 / g_settings->sizeCell);
	c->reserve_carbon           += (s->value[RESERVE_C]          * 1e6 / g_settings->sizeCell);
	c->fruit_carbon             += (s->value[FRUIT_C]            * 1e6 / g_settings->sizeCell);
 
	c->agb                      += ((s->value[LEAF_C] + s->value[STEM_C] + s->value[BRANCH_C] + s->value[FRUIT_C]) * 1e6 / g_settings->sizeCell);
	c->bgb                      += ((s->value[FROOT_C] + s->value[CROOT_C]) * 1e6 / g_settings->sizeCell);

	/* check */
	CHECK_CONDITION(c->leaf_carbon,                    <=, ZERO);
	CHECK_CONDITION(c->froot_carbon,                   <=, ZERO);
	CHECK_CONDITION(c->stem_carbon,                    <=, ZERO);
	CHECK_CONDITION(c->stem_live_wood_carbon,          <=, ZERO);
	CHECK_CONDITION(c->stem_dead_wood_carbon,          <=, ZERO);
	CHECK_CONDITION(c->croot_carbon,                   <=, ZERO);
	CHECK_CONDITION(c->croot_live_wood_carbon,         <=, ZERO);
	CHECK_CONDITION(c->croot_live_wood_carbon,         <=, ZERO);
	CHECK_CONDITION(c->branch_carbon,                  <=, ZERO);
	CHECK_CONDITION(c->branch_live_wood_carbon,        <=, ZERO);
	CHECK_CONDITION(c->branch_dead_wood_carbon,        <=, ZERO);
	CHECK_CONDITION(c->reserve_carbon,                 <=, ZERO);
	CHECK_CONDITION(c->fruit_carbon,                   <=, ZERO);
	CHECK_CONDITION(c->agb,                            <=, ZERO);
	CHECK_CONDITION(c->bgb,                            <=, ZERO);

}

void initialization_forest_class_N (cell_t *const c, const int height, const int dbh, const int age, const int species)
{
	species_t *s;
	s = &c->heights[height].dbhs[dbh].ages[age].species[species];

	logger(g_debug_log,"\n*******INITIALIZE FOREST CLASS NITROGEN POOLS (%s)*******\n", s->name);

	/* leaf */
	if ( ! s->value[LEAF_C] )
	{
		s->value[LEAF_N]       = 0.;
		s->value[LEAF_SUN_N]   = 0.;
		s->value[LEAF_SHADE_N] = 0.;
	}
	else
	{
		s->value[LEAF_N]       = s->value[LEAF_C]       / s->value[CN_LEAVES];
		s->value[LEAF_SUN_N]   = s->value[LEAF_SUN_C]   / s->value[CN_LEAVES];
		s->value[LEAF_SHADE_N] = s->value[LEAF_SHADE_C] / s->value[CN_LEAVES];
	}

	/* fine root */
	if ( ! s->value[FROOT_C] )
	{
		s->value[FROOT_N] = 0.;
	}
	else
	{
		s->value[FROOT_N]       = s->value[FROOT_C] / s->value[CN_FINE_ROOTS];
	}

	s->value[DEADWOODN]         = c->init_dead_C / s->value[CN_DEADWOOD];

	if ( ! c->init_litter_N || c->init_litter_N == NO_DATA )
	{
		s->value[LEAF_LITRN]    = (s->value[LEAF_LITRC]  / s->value[CN_LEAVES])    * (1. - N_FRAC_TO_RETRANSL);
		s->value[FROOT_LITRN]   = s->value[FROOT_LITRC]  / s->value[CN_FINE_ROOTS] * (1. - N_FRAC_TO_RETRANSL);
	}
	else
	{
		/* get data from soil setting file */
		s->value[LEAF_LITRN]    = c->init_litter_N * s->value[FROOT_LEAF_FRAC];
		s->value[FROOT_LITRN]   = c->init_litter_N * ( 1. - s->value[FINE_ROOT_LEAF]);
	}

	s->value[RESERVE_N]         = 0.;

	/* stem */
	s->value[STEM_LIVEWOOD_N]   = s->value[STEM_LIVEWOOD_C] / s->value[CN_LIVEWOOD];
	s->value[STEM_DEADWOOD_N]   = s->value[STEM_DEADWOOD_C] / s->value[CN_DEADWOOD];
	s->value[STEM_N]            = s->value[STEM_LIVEWOOD_N] + s->value[STEM_DEADWOOD_N];

	/* coarse root */
	s->value[CROOT_LIVEWOOD_N]  = s->value[CROOT_LIVEWOOD_C] / s->value[CN_LIVEWOOD];
	s->value[CROOT_DEADWOOD_N]  = s->value[CROOT_DEADWOOD_C] / s->value[CN_DEADWOOD];
	s->value[CROOT_N]           = s->value[CROOT_LIVEWOOD_N] + s->value[CROOT_DEADWOOD_N];

	/* branch */
	s->value[BRANCH_LIVEWOOD_N] = s->value[BRANCH_LIVEWOOD_C] / s->value[CN_LIVEWOOD];
	s->value[BRANCH_DEADWOOD_N] = s->value[BRANCH_DEADWOOD_C] / s->value[CN_DEADWOOD];
	s->value[BRANCH_N]          = s->value[BRANCH_LIVEWOOD_N] + s->value[BRANCH_DEADWOOD_N];

	/* check that all mandatory variables are initialized */
	CHECK_CONDITION(s->value[STEM_N],          <=, ZERO);
	CHECK_CONDITION(s->value[CROOT_N],         <=, ZERO);
	CHECK_CONDITION(s->value[BRANCH_N],        <=, ZERO);
	CHECK_CONDITION(fabs( ( s->value[LEAF_SUN_N] + s->value[LEAF_SHADE_N] ) - s->value[LEAF_N]) , > , eps);

	/* just for evergreen */
	if ( s->value[PHENOLOGY] == 1.1 || s->value[PHENOLOGY] == 1.2 )
	{
		CHECK_CONDITION(s->value[LEAF_N],      <=, ZERO);
		CHECK_CONDITION(s->value[FROOT_N],     <=, ZERO);
	}
}

void initialization_forest_cell_N (cell_t *const c, const int height, const int dbh, const int age, const int species)
{
	species_t *s;
	s = &c->heights[height].dbhs[dbh].ages[age].species[species];

	/*** update at cell level (gN/m2)  ***/
	c->leaf_nitrogen              += (s->value[LEAF_N]            * 1e6 / g_settings->sizeCell);
	c->froot_nitrogen             += (s->value[FROOT_N]           * 1e6 / g_settings->sizeCell);
	c->stem_nitrogen              += (s->value[STEM_N]            * 1e6 / g_settings->sizeCell);
	c->stem_live_wood_nitrogen    += (s->value[STEM_LIVEWOOD_N]   * 1e6 / g_settings->sizeCell);
	c->stem_dead_wood_nitrogen    += (s->value[STEM_DEADWOOD_N]   * 1e6 / g_settings->sizeCell);
	c->croot_nitrogen             += (s->value[CROOT_N]           * 1e6 / g_settings->sizeCell);
	c->croot_live_wood_nitrogen   += (s->value[CROOT_LIVEWOOD_N]  * 1e6 / g_settings->sizeCell);
	c->croot_dead_wood_nitrogen   += (s->value[CROOT_DEADWOOD_N]  * 1e6 / g_settings->sizeCell);
	c->branch_nitrogen            += (s->value[BRANCH_N]          * 1e6 / g_settings->sizeCell);
	c->branch_live_wood_nitrogen  += (s->value[BRANCH_LIVEWOOD_N] * 1e6 / g_settings->sizeCell);
	c->branch_dead_wood_nitrogen  += (s->value[BRANCH_DEADWOOD_N] * 1e6 / g_settings->sizeCell);
	c->reserve_nitrogen           += (s->value[RESERVE_N]         * 1e6 / g_settings->sizeCell);
	c->fruit_nitrogen             += (s->value[FRUIT_N]           * 1e6 / g_settings->sizeCell);

	/* check */
	CHECK_CONDITION(c->leaf_nitrogen,             <=, ZERO);
	CHECK_CONDITION(c->froot_nitrogen,            <=, ZERO);
	CHECK_CONDITION(c->stem_nitrogen,             <=, ZERO);
	CHECK_CONDITION(c->stem_live_wood_nitrogen,   <=, ZERO);
	CHECK_CONDITION(c->stem_dead_wood_nitrogen,   <=, ZERO);
	CHECK_CONDITION(c->croot_nitrogen,            <=, ZERO);
	CHECK_CONDITION(c->croot_live_wood_nitrogen,  <=, ZERO);
	CHECK_CONDITION(c->croot_live_wood_nitrogen,  <=, ZERO);
	CHECK_CONDITION(c->branch_nitrogen,           <=, ZERO);
	CHECK_CONDITION(c->branch_live_wood_nitrogen, <=, ZERO);
	CHECK_CONDITION(c->branch_dead_wood_nitrogen, <=, ZERO);
	CHECK_CONDITION(c->reserve_nitrogen,          <=, ZERO);
	CHECK_CONDITION(c->fruit_nitrogen,            <=, ZERO);


}

void initialization_forest_class_litter (cell_t *const c, const int height, const int dbh, const int age, const int species)
{
	double cwd_litrC;
	double cwd_litr2C;
	double cwd_litr3C;
	double cwd_litr4C;
	double leaf_litrC;
	double leaf_litr1C;
	double leaf_litr2C;
	double leaf_litr3C;
	double leaf_litr4C;
	double froot_litrC;
	double froot_litr1C;
	double froot_litr2C;
	double froot_litr3C;
	double froot_litr4C;
	double cwd_litrN;
	double cwd_litr2N;
	double cwd_litr3N;
	double cwd_litr4N;
	double leaf_litrN;
	double leaf_litr1N;
	double leaf_litr2N;
	double leaf_litr3N;
	double leaf_litr4N;
	double froot_litrN;
	double froot_litr1N;
	double froot_litr2N;
	double froot_litr3N;
	double froot_litr4N;

	species_t *s;
	s = &c->heights[height].dbhs[dbh].ages[age].species[species];

	/* if spinup is off */
	if ( ! g_settings->spinup )
	{
		/*** compute coarse woody debris carbon pools ****/
		if ( ! c->init_dead_C || c->init_dead_C == NO_DATA )
		{
			//note:this must be initialized although to a minimum value to avoid model crashes
			s->value[CWD_LITRC]  = 0.001;
			s->value[CWD_LITRC]  = 10.;  // ddalmo rough estimates for the forest Mg C /ha
			s->value[CWD_LITRC]  = 3.;  // we use this value to initialize the litter pool Mg C /ha
		}
		else
		{
			s->value[CWD_LITRC]  =   c->init_dead_C;
		}
		s->value[CWD_LITR2C]     = s->value[CWD_LITRC]    * s->value[DEADWOOD_USCEL_FRAC];
		s->value[CWD_LITR3C]     = s->value[CWD_LITRC]    * s->value[DEADWOOD_SCEL_FRAC] ;
		s->value[CWD_LITR4C]     = s->value[CWD_LITRC]    * s->value[DEADWOOD_LIGN_FRAC] ;
		/* check */
		CHECK_CONDITION ( fabs ( s->value[CWD_LITR2C] + s->value[CWD_LITR3C] + s->value[CWD_LITR4C] ), > , s->value[CWD_LITRC] + eps);

		cwd_litrC                = s->value[CWD_LITRC];
		cwd_litr2C               = s->value[CWD_LITR2C];
		cwd_litr3C               = s->value[CWD_LITR3C];
		cwd_litr4C               = s->value[CWD_LITR4C];

		/*** compute leaf litter carbon pools ****/
		s->value[LEAF_LITR1C]         = s->value[LEAF_LITRC]        * s->value[LEAF_LITR_LAB_FRAC]  ;
		s->value[LEAF_LITR2C]         = s->value[LEAF_LITRC]        * s->value[LEAF_LITR_USCEL_FRAC];
		s->value[LEAF_LITR3C]         = s->value[LEAF_LITRC]        * s->value[LEAF_LITR_SCEL_FRAC] ;
		s->value[LEAF_LITR4C]         = s->value[LEAF_LITRC]        * s->value[LEAF_LITR_LIGN_FRAC] ;
		/* check */
		CHECK_CONDITION ( fabs ( s->value[LEAF_LITR1C] + s->value[LEAF_LITR2C] + s->value[LEAF_LITR3C] + s->value[LEAF_LITR4C] ), > , s->value[LEAF_LITRC] + eps);

		leaf_litrC                    = s->value[LEAF_LITRC];
		leaf_litr1C                   = s->value[LEAF_LITR1C];
		leaf_litr2C                   = s->value[LEAF_LITR2C];
		leaf_litr3C                   = s->value[LEAF_LITR3C];
		leaf_litr4C                   = s->value[LEAF_LITR4C];

		/*** compute froot litter carbon pools ****/
		s->value[FROOT_LITR1C]       = s->value[FROOT_LITRC]       * s->value[FROOT_LITR_LAB_FRAC]  ;
		s->value[FROOT_LITR2C]       = s->value[FROOT_LITRC]       * s->value[FROOT_LITR_USCEL_FRAC];
		s->value[FROOT_LITR3C]       = s->value[FROOT_LITRC]       * s->value[FROOT_LITR_SCEL_FRAC] ;
		s->value[FROOT_LITR4C]       = s->value[FROOT_LITRC]       * s->value[FROOT_LITR_LIGN_FRAC] ;
		/* check */
		CHECK_CONDITION ( fabs ( s->value[FROOT_LITR1C] + s->value[FROOT_LITR2C] + s->value[FROOT_LITR3C] + s->value[FROOT_LITR4C]), > , s->value[FROOT_LITRC] + eps);

		froot_litrC                   = s->value[FROOT_LITRC];
		froot_litr1C                  = s->value[FROOT_LITR1C];
		froot_litr2C                  = s->value[FROOT_LITR2C];
		froot_litr3C                  = s->value[FROOT_LITR3C];
		froot_litr4C                  = s->value[FROOT_LITR4C];

	}
	else
	{
		/* spinup initialization pools */
		cwd_litrC                     = 0.001;
		cwd_litr2C                    = 0.;
		cwd_litr3C                    = 0.;
		cwd_litr4C                    = 0.;

		leaf_litrC                    = 0.001;
		leaf_litr1C                   = 0.;
		leaf_litr2C                   = 0.;
		leaf_litr3C                   = 0.;
		leaf_litr4C                   = 0.;

		froot_litrC                   = 0.001;
		froot_litr1C                  = 0.;
		froot_litr2C                  = 0.;
		froot_litr3C                  = 0.;
		froot_litr4C                  = 0.;

	}

	c->daily_cwd_to_litrC  += cwd_litrC    * 1e6 / g_settings->sizeCell;
	c->daily_cwd_to_litr2C += cwd_litr2C   * 1e6 / g_settings->sizeCell;
	c->daily_cwd_to_litr3C += cwd_litr3C   * 1e6 / g_settings->sizeCell;
	c->daily_cwd_to_litr4C += cwd_litr4C   * 1e6 / g_settings->sizeCell;
	/* check */
	CHECK_CONDITION ( fabs ( c->daily_cwd_to_litr2C + c->daily_cwd_to_litr3C + c->daily_cwd_to_litr4C ), >, c->daily_cwd_to_litrC+ eps);

	c->daily_leaf_to_litrC      += leaf_litrC        * 1e6 / g_settings->sizeCell;
	c->daily_leaf_to_litr1C     += leaf_litr1C       * 1e6 / g_settings->sizeCell;
	c->daily_leaf_to_litr2C     += leaf_litr2C       * 1e6 / g_settings->sizeCell;
	c->daily_leaf_to_litr3C     += leaf_litr3C       * 1e6 / g_settings->sizeCell;
	c->daily_leaf_to_litr4C     += leaf_litr4C       * 1e6 / g_settings->sizeCell;
	/* check */
	CHECK_CONDITION ( fabs ( c->daily_leaf_to_litr1C + c->daily_leaf_to_litr2C + c->daily_leaf_to_litr3C + c->daily_leaf_to_litr4C ), >, c->daily_leaf_to_litrC+ eps);

	c->daily_froot_to_litrC     += froot_litrC       * 1e6 / g_settings->sizeCell;
	c->daily_froot_to_litr1C    += froot_litr1C      * 1e6 / g_settings->sizeCell;
	c->daily_froot_to_litr2C    += froot_litr2C      * 1e6 / g_settings->sizeCell;
	c->daily_froot_to_litr3C    += froot_litr3C      * 1e6 / g_settings->sizeCell;
	c->daily_froot_to_litr4C    += froot_litr4C      * 1e6 / g_settings->sizeCell;
	/* check */
	CHECK_CONDITION ( fabs ( c->daily_froot_to_litr1C + c->daily_froot_to_litr2C + c->daily_froot_to_litr3C + c->daily_froot_to_litr4C ), >, c->daily_froot_to_litrC+ eps);

	/******************************************************************************************************************************************************************/

	/* if spinup is off */
	if ( ! g_settings->spinup )
	{
		/*** compute coarse woody debris nitrogen pools ****/

		s->value[CWD_LITRN]     = s->value[CWD_LITRC]    / s->value[CN_DEADWOOD];
		s->value[CWD_LITR2N]    = s->value[CWD_LITRN]    * s->value[DEADWOOD_USCEL_FRAC];
		s->value[CWD_LITR3N]    = s->value[CWD_LITRN]    * s->value[DEADWOOD_SCEL_FRAC] ;
		s->value[CWD_LITR4N]    = s->value[CWD_LITRN]    * s->value[DEADWOOD_LIGN_FRAC] ;
		/* check */
		CHECK_CONDITION ( fabs ( s->value[CWD_LITR2N] + s->value[CWD_LITR3N] + s->value[CWD_LITR4N]), >, s->value[CWD_LITRN] + eps);

		cwd_litrN                = s->value[CWD_LITRN];
		cwd_litr2N               = s->value[CWD_LITR2N];
		cwd_litr3N               = s->value[CWD_LITR3N];
		cwd_litr4N               = s->value[CWD_LITR4N];

		/*** compute leaf litter nitrogen pools ****/
		s->value[LEAF_LITR1N]        = s->value[LEAF_LITRN]        * s->value[LEAF_LITR_LAB_FRAC]  ;
		s->value[LEAF_LITR2N]        = s->value[LEAF_LITRN]        * s->value[LEAF_LITR_USCEL_FRAC];
		s->value[LEAF_LITR3N]        = s->value[LEAF_LITRN]        * s->value[LEAF_LITR_SCEL_FRAC] ;
		s->value[LEAF_LITR4N]        = s->value[LEAF_LITRN]        * s->value[LEAF_LITR_LIGN_FRAC] ;
		/* check */
		CHECK_CONDITION ( fabs ( s->value[LEAF_LITR1N] + s->value[LEAF_LITR2N] + s->value[LEAF_LITR3N] + s->value[LEAF_LITR4N] ), >, s->value[LEAF_LITRN] + eps);

		leaf_litrN                    = s->value[LEAF_LITRN];
		leaf_litr1N                   = s->value[LEAF_LITR1N];
		leaf_litr2N                   = s->value[LEAF_LITR2N];
		leaf_litr3N                   = s->value[LEAF_LITR3N];
		leaf_litr4N                   = s->value[LEAF_LITR4N];

		/*** compute froot litter nitrogen pools ****/
		s->value[FROOT_LITR1N]       = s->value[FROOT_LITRN]       * s->value[FROOT_LITR_LAB_FRAC]  ;
		s->value[FROOT_LITR2N]       = s->value[FROOT_LITRN]       * s->value[FROOT_LITR_USCEL_FRAC];
		s->value[FROOT_LITR3N]       = s->value[FROOT_LITRN]       * s->value[FROOT_LITR_SCEL_FRAC] ;
		s->value[FROOT_LITR4N]       = s->value[FROOT_LITRN]       * s->value[FROOT_LITR_LIGN_FRAC] ;
		/* check */
		CHECK_CONDITION ( fabs ( s->value[FROOT_LITR1N] + s->value[FROOT_LITR2N] + s->value[FROOT_LITR3N] + s->value[FROOT_LITR4N]), >,  s->value[FROOT_LITRN] + eps);

		froot_litrN                   = s->value[FROOT_LITRN];
		froot_litr1N                  = s->value[FROOT_LITR1N];
		froot_litr2N                  = s->value[FROOT_LITR2N];
		froot_litr3N                  = s->value[FROOT_LITR3N];
		froot_litr4N                  = s->value[FROOT_LITR4N];


	}
	else
	{
		/* spinup initialization pools */
		cwd_litrN                = 0.001;
		cwd_litr2N               = 0.;
		cwd_litr3N               = 0.;
		cwd_litr4N               = 0.;

		leaf_litrN                    = 0.001;
		leaf_litr1N                   = 0.;
		leaf_litr2N                   = 0.;
		leaf_litr3N                   = 0.;
		leaf_litr4N                   = 0.;

		froot_litrN                   = 0.001;
		froot_litr1N                  = 0.;
		froot_litr2N                  = 0.;
		froot_litr3N                  = 0.;
		froot_litr4N                  = 0.;
	}

	c->daily_cwd_to_litrN       += cwd_litrN    * 1e6 / g_settings->sizeCell;
	c->daily_cwd_to_litr2N      += cwd_litr2N   * 1e6 / g_settings->sizeCell;
	c->daily_cwd_to_litr3N      += cwd_litr3N   * 1e6 / g_settings->sizeCell;
	c->daily_cwd_to_litr4N      += cwd_litr4N   * 1e6 / g_settings->sizeCell;
	/* check */
	CHECK_CONDITION ( fabs ( c->daily_cwd_to_litr2C + c->daily_cwd_to_litr3C + c->daily_cwd_to_litr4C ), >, c->daily_cwd_to_litrC+ eps);

	c->daily_leaf_to_litrN      += leaf_litrN        * 1e6 / g_settings->sizeCell;
	c->daily_leaf_to_litr1N     += leaf_litr1N       * 1e6 / g_settings->sizeCell;
	c->daily_leaf_to_litr2N     += leaf_litr2N       * 1e6 / g_settings->sizeCell;
	c->daily_leaf_to_litr3N     += leaf_litr3N       * 1e6 / g_settings->sizeCell;
	c->daily_leaf_to_litr4N     += leaf_litr4N       * 1e6 / g_settings->sizeCell;
	/* check */
	CHECK_CONDITION ( fabs ( c->daily_leaf_to_litr1C + c->daily_leaf_to_litr2C + c->daily_leaf_to_litr3C + c->daily_leaf_to_litr4C ), >, c->daily_leaf_to_litrC+ eps);

	c->daily_froot_to_litrN     += froot_litrN      * 1e6 / g_settings->sizeCell;
	c->daily_froot_to_litr1N    += froot_litr1N      * 1e6 / g_settings->sizeCell;
	c->daily_froot_to_litr2N    += froot_litr2N      * 1e6 / g_settings->sizeCell;
	c->daily_froot_to_litr3N    += froot_litr3N      * 1e6 / g_settings->sizeCell;
	c->daily_froot_to_litr4N    += froot_litr4N     * 1e6 / g_settings->sizeCell;
	/* check */
	CHECK_CONDITION ( fabs ( c->daily_froot_to_litr1N + c->daily_froot_to_litr2N + c->daily_froot_to_litr3N + c->daily_froot_to_litr4N ), >, c->daily_froot_to_litrN+ eps);

	/******************************************************************************************************************************************************************/

}
void initialization_cell_litter_biochem ( cell_t *const c )
{
	/******************************************************************************************************************************************************************/

	/*** create cwd carbon pools fractions (gC/m2) ***/
	//fixme it should sum throughout all classes...
	c->cwd_C     = c->daily_cwd_to_litrC;
	c->cwd_2C    = c->daily_cwd_to_litr2C;
	c->cwd_3C    = c->daily_cwd_to_litr3C;
	c->cwd_4C    = c->daily_cwd_to_litr4C;
	/* check */
	CHECK_CONDITION ( fabs ( c->cwd_2C + c->cwd_3C + c->cwd_4C ) , > , c->cwd_C + eps );

	/*** create leaf litter carbon pools fractions (gC/m2) ****/
	c->leaf_litrC     = c->daily_leaf_to_litrC;
	c->leaf_litr1C    = c->daily_leaf_to_litr1C;
	c->leaf_litr2C    = c->daily_leaf_to_litr2C;
	c->leaf_litr3C    = c->daily_leaf_to_litr3C;
	c->leaf_litr4C    = c->daily_leaf_to_litr4C;
	/* check */
	CHECK_CONDITION ( fabs ( c->leaf_litr1C + c->leaf_litr2C + c->leaf_litr3C + c->leaf_litr4C ) , > , c->leaf_litrC + eps );

	/*** create froot litter carbon pools fractions (gN/m2) ****/
	c->froot_litrC    = c->daily_froot_to_litrC;
	c->froot_litr1C   = c->daily_froot_to_litr1C;
	c->froot_litr2C   = c->daily_froot_to_litr2C;
	c->froot_litr3C   = c->daily_froot_to_litr3C;
	c->froot_litr4C   = c->daily_froot_to_litr4C;
	/* check */
	CHECK_CONDITION ( fabs ( c->froot_litr1C + c->froot_litr2C + c->froot_litr3C + c->froot_litr4C ) , > , c->froot_litrC + eps );

	/* cumulate overall */
	c->litrC          = c->leaf_litrC + c->froot_litrC + c->cwd_C;

	/* cumulate carbon pools */
	c->litr1C         = c->leaf_litr1C + c->froot_litr1C;
	c->litr2C         = c->leaf_litr2C + c->froot_litr2C + c->cwd_2C;
	c->litr3C         = c->leaf_litr3C + c->froot_litr3C + c->cwd_3C;
	c->litr4C         = c->leaf_litr4C + c->froot_litr4C + c->cwd_4C;

	/******************************************************************************************************************************************************************/

	/*** create cwd nitrogen pools fractions (gN/m2) ***/
	//fixme it should sum throghout all classes...
	c->cwd_N     = c->daily_cwd_to_litrN;
	c->cwd_2N    = c->daily_cwd_to_litr2N;;
	c->cwd_3N    = c->daily_cwd_to_litr3N;;
	c->cwd_4N    = c->daily_cwd_to_litr4N;;
	/* check */
	CHECK_CONDITION ( fabs ( c->cwd_2N + c->cwd_3N + c->cwd_4N ) , > , c->cwd_N + eps );

	/*** create leaf litter nitrogen pools fractions (gN/m2) ****/
	c->leaf_litrN     = c->daily_leaf_to_litrN;
	c->leaf_litr1N    = c->daily_leaf_to_litr1N;
	c->leaf_litr2N    = c->daily_leaf_to_litr2N;
	c->leaf_litr3N    = c->daily_leaf_to_litr3N;
	c->leaf_litr4N    = c->daily_leaf_to_litr4N;
	/* check */
	CHECK_CONDITION ( fabs ( c->leaf_litr1N + c->leaf_litr2N + c->leaf_litr3N + c->leaf_litr4N ) , > , c->leaf_litrN + eps );

	/*** create froot litter nitrogen pools fractions (gN/m2) ****/
	c->froot_litrN    = c->daily_froot_to_litrN;
	c->froot_litr1N   = c->daily_froot_to_litr1N;
	c->froot_litr2N   = c->daily_froot_to_litr2N;
	c->froot_litr3N   = c->daily_froot_to_litr3N;
	c->froot_litr4N   = c->daily_froot_to_litr4N;
	/* check */
	CHECK_CONDITION ( fabs ( c->froot_litr1N + c->froot_litr2N + c->froot_litr3N + c->froot_litr4N ) , > , c->froot_litrN + eps );

	/* cumulate overall */
	c->litrN          = c->leaf_litrN + c->froot_litrN + c->cwd_N;

	/* cumulate nitrogen pools */
	c->litr1N         = c->leaf_litr1N + c->froot_litr1N;
	c->litr2N         = c->leaf_litr2N + c->froot_litr2N + c->cwd_2N;
	c->litr3N         = c->leaf_litr3N + c->froot_litr3N + c->cwd_3N;
	c->litr4N         = c->leaf_litr4N + c->froot_litr4N + c->cwd_4N;

}
void initialization_cell_soil_biochem (cell_t *const c)
{
	/* initialize soil carbon */
	if ( ! c->init_soil_C || c->init_soil_C == NO_DATA )
	{
		//c->soilC      = 0.001;
		c->soilC      = 0.;
	}
	else
	{
		c->soilC      = c->init_soil_C;
	}

    c->soil1C         = 0.;
	c->soil2C         = 0.;
	c->soil3C         = 0.;
	c->soil4C         = 0.;

	/* initialize soil nitrogen */
	if ( ! c->init_soil_N || c->init_soil_N == NO_DATA )
	{
		c->soilN      = 0.001;
	}
	else
	{
		c->soilN      = c->init_soil_N;
	}
	c->soil1N         = 0.001;
	c->soil2N         = 0.001;
	c->soil3N         = 0.001;
	c->soil4N         = 0.001;
}

void initialization_cell_soil_physic(cell_t *const c)   
{
	float acoeff;
	float bcoeff;
	float sat;
	float volumetric_wilting_point;
	float volumetric_field_capacity;
	float volumetric_saturated_hydraulic_conductivity;

	logger(g_debug_log,"\n******* INITIALIZE SOIL *******\n");

	/*soil matric potential*/
	CHECK_CONDITION(fabs((c->sand_perc + c->clay_perc + c->silt_perc) -100.0 ), >, eps);

	/* BIOME-BGC METHOD */
	/*
	For further discussion see:
	Cosby, B.J., G.M. Hornberger, R.B. Clapp, and T.R. Ginn, 1984.  A
	   statistical exploration of the relationships of soil moisture
	   characteristics to the physical properties of soils.  Water Res.
	   Res. 20:682-690.

	Saxton, K.E., W.J. Rawls, J.S. Romberger, and R.I. Papendick, 1986.
		Estimating generalized soil-water characteristics from texture.
		Soil Sci. Soc. Am. J. 50:1031-1036.
	 */

	logger(g_debug_log, "BIOME soil characteristics\n");


	/* (DIM) Clapp-Hornberger "b" parameter */
	c->soil_b = -( 3.10 + 0.157 * c->clay_perc - 0.003 * c->sand_perc); /* ok for schwalm*/
	logger(g_debug_log, "soil_b = %f (DIM)\n", c->soil_b);

	//test
	/* following Rawls et al., 1992 and Schwalm et al., 2004 */
	/*texture-dependent empirical coefficinet */
	// c->soil_b = 11.43 - (0.1034*g_soil_settings->values[SOIL_SAND_PERC) - (0.0687*0.157*g_soil_settings->values[SOIL_silt_perc);

	/* calculate the soil pressure-volume coefficients from texture data */
	/* Uses the multivariate regressions from Cosby et al., 1984 */
	/* (DIM) Soil volumetric water content at saturation */
	c->vwc_sat = ( 50.5 - 0.142 * c->sand_perc - 0.037 * c->clay_perc ) / 100.; /* ok for schwalm*/
	logger(g_debug_log, "volumetric water content at saturation (BIOME) = %f %%(vol)\n", c->vwc_sat);
	/* (MPa) soil matric potential at saturation */
	c->psi_sat = -( exp ( ( 1.54 - 0.0095 * c->sand_perc + 0.0063 * c->silt_perc ) * log ( 10.) ) * 9.8e-5 ); /* ok for schwalm*/
	logger(g_debug_log, "psi_sat = %f MPa \n", c->psi_sat);
	/* Clapp-Hornenberger function 1978 (DIM) Soil Field Capacity Volumetric Water Content at field capacity ( = -0.015 MPa) */
	c->vwc_fc = c->vwc_sat * pow ( ( -0.015 / c->psi_sat ) , ( 1. / c->soil_b ) );
	logger(g_debug_log, "volumetric water content at field capacity (BIOME) = %f %%(vol) \n", c->vwc_fc);

	/* define maximum soilwater content, for outflow calculation */
	/* converts volumetric water content (m3/m3) --> (kg/m2) */

	/* (kgH2O/m2) soil water at field capacity */
	c->soilw_fc = ( c->soil_depth / 100. ) * c->vwc_fc * 1e3;
	logger(g_debug_log, "soilw_fc BIOME (MAXASW FC BIOME)= %f (kgH2O/m2)\n", c->soilw_fc);
	//equal to MAXASW

	/* (kgH2O/m2) soil water at saturation */
	c->soilw_sat = ( c->soil_depth / 100. ) * c->vwc_sat * 1e3;
	logger(g_debug_log, "soilw_sat BIOME (MAXASW SAT BIOME)= %f (kgH2O/m2)\n", c->soilw_sat);

	/* (kgH2O/m2) maximum soil water at field capacity */
	c->max_asw_fc = c->soilw_fc;

	/* compute initialization soil water content */
	c->asw = (c->soilw_sat * g_settings->init_frac_maxasw);
	logger(g_debug_log, "Initialization ASW = %f (mm-kgH2O/m2)\n\n\n", c->asw);


	/** CENTURY METHOD **/

	/* soil data from https://www.nrel.colostate.edu/projects/century/soilCalculatorHelp.htm */
	/* following Saxton et al 1986, 2006, 2008 */
	logger(g_debug_log, "CENTURY soil characteristics\n");
	acoeff = (float)exp(-4.396 - 0.0715 * c->clay_perc - 4.88e-4 * pow(c->sand_perc, 2) - 4.285e-5 * pow( c->sand_perc, 2 ) * c->clay_perc);
	bcoeff = (float)(-3.14 - 0.00222 * pow(c->clay_perc, 2) - 3.484e-5 * pow(c->sand_perc, 2) * c->clay_perc);
	sat    = (float)(0.332 - 7.251e-4 * c->sand_perc + 0.1276 * log10(c->clay_perc));

	/* volumetric percentage for wilting point */
	volumetric_wilting_point = (float)pow((15.0/acoeff), (1.0/bcoeff));
	/* volumetric percentage for field capacity */
	volumetric_field_capacity = (float)pow((0.333/acoeff),(1.0/bcoeff));
	/* volumetric percentage for saturated hydraulic conductivity */
	volumetric_saturated_hydraulic_conductivity = (float)exp((12.012 - 0.0755 * c->sand_perc) + (-3.895 + 0.03671 * c->sand_perc - 0.1103 * c->clay_perc + 8.7546e-4 * pow(c->clay_perc, 2))/sat);
	/* bulk density g/cm3 */
	c->bulk_density = (1 - sat) * 2.65;

	/* corrections from Steve Del Grosso */
	/* volumetric percentage wilting point */
	volumetric_wilting_point += (float)(-0.15 * volumetric_wilting_point);
	logger(g_debug_log, "volumetric water content at wilting point (CENTURY) = %f %%(vol)\n", volumetric_wilting_point);
	/* (kgH2O/m2) soil water at wilting point */
	c->wilting_point = ( c->soil_depth / 100. ) * volumetric_wilting_point * 1e3;
	logger(g_debug_log, "Wilting point (CENTURY) = %f mm/m2\n", c->wilting_point);

	/* volumetric percentage field capacity */
	volumetric_field_capacity += (float)(0.07 * volumetric_field_capacity);
	logger(g_debug_log, "volumetric water content at field capacity (CENTURY) = %f %%(vol)\n", volumetric_field_capacity);
	/* (kgH2O/m2) soil water at field capacity */
	c->field_capacity = ( c->soil_depth / 100. ) * volumetric_field_capacity * 1e3;
	logger(g_debug_log, "Field capacity (CENTURY) = %f mm/m2\n", c->field_capacity);

	/* volumetric percentage saturated hydraulic conductivity */
	volumetric_saturated_hydraulic_conductivity /= 1500.;
	logger(g_debug_log, "volumetric water content at saturated hydraulic conductance (CENTURY) = %f %%(vol)\n", volumetric_saturated_hydraulic_conductivity);

	/* (kgH2O/m2) soil water at saturated hydraulic conductivity */
	c->sat_hydr_conduct = (c->soil_depth / 100.) * volumetric_saturated_hydraulic_conductivity * 1e3;
	logger(g_debug_log, "Saturated hydraulic conductivity (CENTURY) = %f mm/m2\n", c->sat_hydr_conduct);

	/* bulk density g/cm3 */
	c->bulk_density += (-0.08 * c->bulk_density);
	logger(g_debug_log, "Bulk density = %f g/cm^3\n", c->bulk_density);
	logger(g_debug_log, "***************************************************\n\n");

}

