/*
 * tree_model.c
 *
 *  Created on: 08/nov/2016
 *      Author: alessio-cmcc
 */
/* tree_model.c */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <new_forest_tree_class.h>
#include <N-assimilation.h>
#include "constants.h"
#include "common.h"
#include "print_output.h"
#include "settings.h"
#include "logger.h"
#include "matrix.h"
#include "g-function.h"
#include "structure.h"
#include "cue.h"
#include "wue.h"
#include "check_balance.h"
#include "print.h"
#include "utility.h"
#include "dendometry.h"
#include "snow.h"
#include "phenology.h"
#include "peak_lai.h"
#include "canopy_radiation_sw_band.h"
#include "canopy_radiation_lw_band.h"
#include "canopy_net_radiation.h"
#include "canopy_temperature.h"
#include "modifiers.h"
#include "N-assimilation.h"
#include "canopy_evapotranspiration.h"
#include "photosynthesis.h"
#include "aut_respiration.h"
#include "fluxes.h"
#include "C-assimilation.h"
#include "C-deciduous-partitioning.h"
#include "C-evergreen-partitioning.h"
#include "CN-allocation.h"
#include "lai.h"
#include "turnover.h"
#include "mortality.h"
#include "biomass.h"
#include "management.h"
#include "soil_evaporation.h"
#include "soil_respiration.h"
#include "evapotranspiration.h"
#include "heat_fluxes.h"
#include "soil_water_balance.h"
#include "leaf_fall.h"
#include "settings.h"
#include "tree_model.h"
#include "remove_tree_class.h"
#include "regeneration.h"
#include "decomposition.h"
#include "littering.h"
#include "CN-balance.h"

extern logger_t* g_debug_log;
extern soil_settings_t* g_soil_settings;
extern settings_t* g_settings;

//extern const char sz_err_out_of_memory[];

//extern const char *szMonth[MONTHS_COUNT];

/* Days in Months */
extern int DaysInMonth[];

/* Last cumulative days in months */
//extern int MonthLength [];
//extern int MonthLength_Leap [];

/*****************************************************************************************************************/
int Tree_model(matrix_t *const m, const int cell, const int halfhour, const int hour, const int day, const int month, const int year)
{
	static int layer;
	static int height;
	int dbh;
	static int age;
	static int species;
	//int management = 0;

	/* shortcuts */
	cell_t *c;
	tree_layer_t *l;
	height_t *h;
	dbh_t *d;
	age_t *a;
	species_t *s;
	meteo_daily_t *meteo_daily;
	meteo_annual_t *meteo_annual;
	
	meteo_daily = NULL;
	
	assert(m);

	/* assign shortcuts */
	c = &m->cells[cell];

	if ( DAILY == g_settings->time )
	{
		meteo_daily = &METEO_DAILY(m->cells[cell].years[year].m)[month].d[day];
		//meteo_daily  = &m->cells[cell].years[year].m[month].d[day];
	}
	else if ( HOURLY == g_settings->time )
	{
		// TODO
	}
	else if ( HALFHOURLY == g_settings->time )
	{
		// TODO
	}
	
	assert(meteo_daily);
	
	meteo_annual = &m->cells[cell].years[year];

	/* check parameters */
	assert( m );

	logger (g_debug_log, "\n********* TREE_MODEL_DAILY *********\n");

	/* compute single tree average biomass */
	average_tree_pools ( c );

	/****************************************************************************/

	if ( ! day && ! month )
	{
		/* management blocks */
		forest_management ( c, day, month, year );
	}

	/****************************************************************************/

	/* annual forest structure (only the year after the first) */
	if ( ( ! day && ! month && year ) && ( ! c->harvesting ) )
	{
		annual_forest_structure ( c, year );
	}

	/* daily forest structure*/
	daily_forest_structure ( c,  meteo_daily);

	/* print forest cell data */
	print_daily_forest_data ( c );

	/* prephenology */
	prephenology ( c, meteo_daily, day, month, year );

	/***********************************************************************************************/

	/* sort heights in descending order */
#ifndef USE_NEW_OUTPUT
	qsort ( c->heights, c->heights_count, sizeof (height_t), sort_by_heights_desc );
#endif

	/* loop on each cell layers starting from highest to lower */
	for ( layer = c->tree_layers_count -1 ; layer >= 0; --layer )
	{
		/* assign shortcut */
		l = &m->cells[cell].tree_layers[layer];

		logger(g_debug_log,"*****************************************************************************\n"
				"                                layer = %d                                 \n"
				"*****************************************************************************\n",layer);

		l->layer_z = layer;

		/* loop on each heights starting from highest to lower */
		for ( height = 0; height <  c->heights_count; ++height )
		{
			/* assign shortcut */
			h = &m->cells[cell].heights[height];

			/* sort dbhs in descending order */
		#ifndef USE_NEW_OUTPUT
			qsort ( h->dbhs, h->dbhs_count, sizeof (dbh_t), sort_by_dbhs_desc );
		#endif
			//ALESSIOC FIXME
			c->cell_heights_count ++;

			/* check if tree height class matches with corresponding cell layer */
			if( h->height_z == l->layer_z )
			{
				logger(g_debug_log,"*****************************************************************************\n"
						"                              height = %f                              \n"
						"*****************************************************************************\n", h->value);

				/* loop on each dbh starting from highest to lower */
				for ( dbh = 0; dbh < h->dbhs_count; ++dbh )
				{
					/* assign shortcut */
					d = &h->dbhs[dbh];

					/* sort ages in descending order */
				#ifndef USE_NEW_OUTPUT
					qsort ( d->ages, d->ages_count, sizeof (age_t), sort_by_ages_desc );
				#endif

					logger(g_debug_log,"*****************************************************************************\n"
							"                              dbh = %f                              \n"
							"*****************************************************************************\n", d->value);

					/* loop on each age class */
					for ( age = 0; age < d->ages_count; ++age )
					{
						/* assign shortcut */
						a = &m->cells[cell].heights[height].dbhs[dbh].ages[age];

						logger(g_debug_log,"*****************************************************************************\n"
								"                                  age = %d                                 \n"
								"*****************************************************************************\n", a->value);

						/* increment age after first year */
						if( !day && !month && year ) a->value += 1;

						/* loop on each species class */
						for ( species = 0; species < a->species_count; ++species )
						{
							/* assign shortcut */
							s = &m->cells[cell].heights[height].dbhs[dbh].ages[age].species[species];

							logger(g_debug_log,"*****************************************************************************\n"
									"*                              species = %s                         *\n"
									"*****************************************************************************\n", s->name);

							/**********************************/

							/* counter for class days, months and years of simulation */
							++s->counter[DOS];
							if ( ! day && ! month ) ++s->counter[MOS];
							if ( ! day && ! month ) ++s->counter[YOS];

							/**********************************/

							/* beginning of simulation (every year included the first one) */
							if ( ( g_soil_settings->values[SOIL_LAT] > 0 && c->doy == 1)  ||
									( g_soil_settings->values[SOIL_LAT] < 0 && c->doy == 180) )
							{
								/* compute annual minimum reserve for incoming year */
								annual_minimum_reserve( s );

								/* compute age-related sla */
								specific_leaf_area ( a, s);

								/* compute annual potential Maximum LAI */
								peak_lai( s, day, month, year );

								/* compute growth respiration fraction */
								growth_respiration_frac ( a, s );
							}

							/***************************************************************/

							/* print at the beginning of simulation forest class data */
							print_daily_forest_class_data ( c, layer, height, dbh, age, species );

							/* compute species-specific phenological phase */
							phenology ( c, layer, height, dbh, age, species, meteo_daily, day, month, year );

							logger(g_debug_log, "--PHYSIOLOGICAL PROCESSES LAYER %d --\n", l->layer_z);

							/*********************************************************************/

							/* short wave band */
							canopy_radiation_sw_band ( c, layer, height, dbh, age, species, meteo_daily );

							/* long wave band */
							canopy_radiation_lw_band ( c, layer, height, dbh, age, species, meteo_daily );

							/* net radiation */
							canopy_net_radiation ( c, layer, height, dbh, age, species );

							/**********************************************************************/

							/* canopy temperature */
							canopy_temperature ( c, layer, height, dbh, age, species, meteo_daily );

							/* daily modifier */
							modifiers ( c, layer, height, dbh, age, species, meteo_daily, meteo_annual );

							/* canopy water fluxes */
							canopy_evapotranspiration ( c, layer, height, dbh, age, species, meteo_daily );

							/* note: following Piao et al., 2010 */
							/* (Maint Resp)->(Growth Resp = (GPP - Maint Resp) * eff_grperc)->(NPP) */

							/* canopy carbon assimilation */
							photosynthesis ( c, layer, height, dbh, age, species, DaysInMonth[month], meteo_annual );

							/* canopy carbon assimilation (BIOME-BGC) */
							total_photosynthesis_biome ( c, height, dbh, age, species, meteo_daily, meteo_annual );

							/* maintenance respiration */
							maintenance_respiration ( c, layer, height, dbh, age, species, meteo_daily );

							/* C-N-partitioning */
							if ( s->value[PHENOLOGY] == 0.1 || s->value[PHENOLOGY] == 0.2 )
							{
								/* deciduous */
								daily_C_deciduous_partitioning ( c, layer, height, dbh, age, species, meteo_daily, day, month, year );
							}
							else
							{
								/* evergreen */
								daily_C_evergreen_partitioning ( c, layer, height, dbh, age, species, meteo_daily, day, month, year );
							}

							/* growth respiration */
							growth_respiration      ( c, layer, height, dbh, age, species );

							/* autotrophic respiration */
							autotrophic_respiration ( c, layer, height, dbh, age, species, meteo_daily );

							/* carbon fluxes */
							carbon_fluxes           ( c, height, dbh, age, species );

							/* C assimilation */
							carbon_assimilation     ( c, height, dbh, age, species );

							if ( c->doy == ( IS_LEAP_YEAR ( c->years[year].year ) ? 366 : 365) )
							{
								/* ISIMIP: exclude age mortality function when management is "var" and year < year start management */
								if ( c->years[year].year > g_settings->year_start_management && g_settings->management != MANAGEMENT_VAR)
								{
									/* Mortality based on tree Age (LPJ) */
									age_mortality ( c, height, dbh, age, species);
								}
							}

							/* allocate daily carbon */
							carbon_allocation       ( c, height, dbh, age, species, day, month, year );

							/* allocate daily nitrogen */
							nitrogen_allocation     ( c, s );

							/* MORTALITY */
							/* Mortality based on growth efficiency */
							/* note: when it happens the overall class is removed */
							if ( ! growth_efficiency_mortality ( c, height, dbh, age, species ) )
							{

								/* turnover */
								turnover ( c, a, s, day, month, year );

								/* carbon use efficiency */
								carbon_use_efficiency ( c, height, dbh, age, species, day, month, year );

								/* water use efficiency */
								water_use_efficiency  ( c, height, dbh, age, species, day, month, year );

								/* update Leaf Area Index */
								daily_lai             ( s );

								/* N assimilation */
								nitrogen_assimilation ( s );

								/* litter fluxes and pools */
								littering             ( c, s );

								/* tree level dendrometry */
								dendrometry_old       ( c, layer, height, dbh, age, species, year );

								/** END OF YEAR **/
								/* last day of the year */
								if ( c->doy == ( IS_LEAP_YEAR ( c->years[year].year ) ? 366 : 365) )
								{

									/* above ground-below ground stocks */
									abg_bgb_biomass ( c, height, dbh, age, species );

									/* annual branch and bark fraction */
									tree_branch_and_bark ( c, height, dbh, age, species );

									/* annual volume, MAI and CAI */
									annual_tree_increment ( c, height, dbh, age, species, year );
								}

								/** check for fluxes and mass balance closure at the tree class level **/

								logger(g_debug_log, "\n**TREE CLASS LEVEL BALANCE**\n");

								/* check for radiative flux balance closure */
								/* 1 */ if ( ! check_tree_class_radiation_flux_balance ( c, layer, height, dbh, age, species ) ) return 0;

								/* check for carbon flux balance closure */
								/* 2 */ if ( ! check_tree_class_carbon_flux_balance    ( c, layer, height, dbh, age, species ) ) return 0;

								/* check for nitrogen flux balance closure */
								/* 3 */ //fixme if ( ! check_tree_class_nitrogen_flux_balance  ( c, layer, height, dbh, age, species ) ) return 0;

								/* check for water flux balance closure */
								/* 4 */ if ( ! check_tree_class_water_flux_balance     ( c, layer, height, dbh, age, species ) ) return 0;

								/* check for carbon mass balance closure */
								/* 5 */ if ( ! check_tree_class_carbon_mass_balance    ( c, layer, height, dbh, age, species ) ) return 0;

								/* check for nitrogen mass balance closure */
								/* 6 */ //fixme if ( ! check_tree_class_nitrogen_mass_balance  ( c, layer, height, dbh, age, species ) ) return 0;
							}
							else
							{
								//fixme here model should remove class just after have checked that the balances are closed
								//so model has to include c fluxes that go out to litter and cwd
								if ( height >= m->cells[cell].heights_count ) goto height_end;
								if ( dbh >= m->cells[cell].heights[height].dbhs_count ) goto dbh_end;
								if ( age >= m->cells[cell].heights[height].dbhs[dbh].ages_count ) goto age_end;
								//if ( species >= m->cells[cell].heights[height].dbhs[dbh].ages[age].species_count ) goto species_end;
							}
							/****************************************************************************************************************************************/
							/****************************************************************************************************************************************/
						}
						logger(g_debug_log, "****************END OF SPECIES CLASS***************\n");
					}
				age_end:
					logger(g_debug_log, "****************END OF AGES CLASS***************\n");
				}
			dbh_end:
				logger(g_debug_log, "****************END OF DBH CLASS***************\n");
			}
		}
	height_end:
		logger(g_debug_log, "****************END OF HEIGHT CLASS***************\n");
	}
	logger(g_debug_log, "****************END OF LAYER CLASS***************\n");

	/* ok */
	return 1;
}


