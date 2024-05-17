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
#include "new_forest_tree_class.h"
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
#include "decomposition.h"
#include "littering.h"
#include "CN-balance.h"
#include "C-fruit-partitioning.h"
#include "regeneration.h"
#include "recruitment.h"
#include "canopy_cover.h"
#include "met_data.h"

extern logger_t* g_debug_log;
//extern soil_settings_t* g_soil_settings;
extern settings_t* g_settings;
extern dataset_t* g_dataset;


//extern const char sz_err_out_of_memory[];

//extern const char *szMonth[MONTHS_COUNT];

/* Days in Months */
//extern int DaysInMonth[];

/* Last cumulative days in months */
//extern int MonthLength [];
//extern int MonthLength_Leap [];

/*****************************************************************************************************************/
/*****************************************************************************************************************/
/*****************************************************************************************************************/
int Tree_model_daily (matrix_t *const m, const int cell, const int day, const int month, const int year)
{
	int layer;
	int height;
	int dbh;
	int age;
	int species;
        int row = 0;
        int year_dens_fin = 0;  // only used if MANAGEMENT == VAR or VAR1


	/* shortcuts */
	cell_t *c;
	tree_layer_t *l;
	height_t *h;
	dbh_t *d;
	age_t *a;
	species_t *s;
	meteo_daily_t *meteo_daily;
	meteo_annual_t *meteo_annual;

	/* assign shortcuts */
	c = &m->cells[cell];
	meteo_daily  = &m->cells[cell].years[year].m[month].d[day];
	meteo_annual = &m->cells[cell].years[year];

	/* check parameters */
	assert( m );

        assert(g_dataset);

	logger (g_debug_log, "\n********* TREE_MODEL_DAILY *********\n");

	/* compute single tree average biomass */
	average_tree_pools ( c );


	/****************************************************************************/

	if ( ! day && ! month )
	{
                printf("\n anno, year= %d\n",c->years[year].year);

               // ddalmo at some point: let's perform the if statement here
               //if (MANAGEMENT_VAR == g_settings->management || MANAGEMENT_ON == g_settings->management)
               // {
		    /* management blocks */
	     	forest_management ( c, day, month, year );

             // }

	}

	/****************************************************************************/

	/* annual forest structure (only the year after the first) */
	if ( ( ! day && ! month && year ) && ( ! c->harvesting ) )
	{

		annual_forest_structure ( c, year );
		//printf("\n after annual_forest_structure\n");

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
							if ( ( c->lat > 0. && c->doy == 1 ) || ( c->lat < 0. && c->doy == 180 ) )
							{

                                                                /* compute annual minimum reserve for incoming year */
								annual_minimum_reserve( s );

								/* compute annual potential Maximum LAI */
								peak_lai                ( a, s, day, month, year );

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
							canopy_net_radiation     ( c, layer, height, dbh, age, species );

							/**********************************************************************/

							/* canopy temperature */
							canopy_temperature ( c, layer, height, dbh, age, species, meteo_daily );

							/* daily modifier */
							modifiers ( c, layer, height, dbh, age, species, meteo_daily, meteo_annual );

							if ( ! g_settings->PSN_mod )
							{
								/**********************************************************************/
								/** FvCB + Jarvis **/

								/* note: following Piao et al., 2010 */
								/* (Maint Resp)->(Growth Resp = (GPP - Maint Resp) * eff_grperc)->(NPP) */

								/* maintenance respiration */
								maintenance_respiration ( c, layer, height, dbh, age, species, meteo_daily );

								/* canopy water fluxes */
								canopy_evapotranspiration ( c, layer, height, dbh, age, species, meteo_daily );

								/* canopy carbon assimilation ( Farquhar Von Caemmerer and Berry - FvCB - approach ) */
								photosynthesis_FvCB     ( c, height, dbh, age, species, meteo_daily, meteo_annual );

								/**********************************************************************/
							}
							else if ( g_settings->PSN_mod == 1 )
							{
								/**********************************************************************/
								/** LUE + Jarvis **/

								/* canopy water fluxes */
								canopy_evapotranspiration ( c, layer, height, dbh, age, species, meteo_daily );

								/* canopy carbon assimilation ( Monteith - LUE - approach ) */
								photosynthesis_LUE      ( c, layer, height, dbh, age, species, meteo_annual );

								/* maintenance respiration */
								maintenance_respiration ( c, layer, height, dbh, age, species, meteo_daily );

								/**********************************************************************/
							}
							else
							{
								//test
								/**********************************************************************/
								/** FvCB + Ball Berry / Medlyn */

								/* note: following Piao et al., 2010 */
								/* (Maint Resp)->(Growth Resp = (GPP - Maint Resp) * eff_grperc)->(NPP) */

								/* maintenance respiration */
								maintenance_respiration  ( c, layer, height, dbh, age, species, meteo_daily );

								/* canopy carbon assimilation ( Farquhar Von Caemmerer and Berry - FvCB - approach ) */
								photosynthesis_FvCB_BB    ( c, height, dbh, age, species, meteo_daily, meteo_annual );

								/* canopy water fluxes */
								canopy_evapotranspiration ( c, layer, height, dbh, age, species, meteo_daily );

								/**********************************************************************/
							}

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
							growth_respiration                 ( c, layer, height, dbh, age, species );

							/* autotrophic respiration */
							autotrophic_respiration            ( c, layer, height, dbh, age, species, meteo_daily );

							/* carbon fluxes */
							carbon_fluxes                      ( c, height, dbh, age, species );

							/* C productivity */
							carbon_productivity                ( c, height, dbh, age, species );



							if ( c->doy == ( IS_LEAP_YEAR ( c->years[year].year ) ? 366 : 365) )
							{

#if 1    //set to 1 . set to 0 only for testing purpose


                                        //if ( g_settings->management == MANAGEMENT_VAR || g_settings->management == MANAGEMENT_VAR1)
                                        if ( g_settings->management == MANAGEMENT_VAR)
                                        {

                                            // compute last year of available stand density data
                                            // (in the stand.txt file each layer or class has to have
                                            // the same number of stand density data)

                                             row = g_dataset->rows_count ;

                                             year_dens_fin = g_dataset->rows[row-1].year_stand;


                                                if (c->years[year].year > year_dens_fin)
                                                    {


                                    /* Mortality based on tree Age (LPJ) */
									age_mortality        ( c, height, dbh, age, species );


                                    /* Mortality based on stochasticity */
									stochastic_mortality ( c, height, dbh, age, species );

                                                        }

                                                    }
                                                    else
                                                    { // with management ON or OFF always compute mortality

									/* Mortality based on tree Age (LPJ) */
									age_mortality        ( c, height, dbh, age, species );


									/* Mortality based on stochasticity */
									//stochastic_mortality ( c, height, dbh, age, species );

									//Saplings mortality
							        saplings_mortality (c, height, dbh, age, species);

							      }


#else   // this option is wrongly formulated: to be excluded. use only for testing

								/* ISIMIP: exclude age mortality function when management is "var"
								//  and year < year start management */
								if ( c->years[year].year > g_settings->year_start_management && g_settings->management != MANAGEMENT_VAR )
								{
									/* Mortality based on tree Age (LPJ) */
									age_mortality        ( c, height, dbh, age, species );


									/* Mortality based on stochasticity */
									stochastic_mortality ( c, height, dbh, age, species );
								}
#endif
							}

							/* allocate daily carbon */


							// carbon_allocation       ( c, a, s, day, month, year ); // old subroutine

							//note: this is basically the new function in version v.5.5
							carbon_allocation_new   ( c, a, s, day, month, year );

							/* allocate daily nitrogen */
							nitrogen_allocation     ( c, s, day, month, year );

							/* MORTALITY */
							/* Mortality based on growth efficiency */
							/* note: when it happens the overall class is removed */

							if ( ! growth_efficiency_mortality ( c, height, dbh, age, species ) )
							{

								/* turnover */
								// turnover ( c, a, s, day, month, year ); // old function

								//note: this is basically the new function in version v.5.5
								sapwood_turnover ( c, a, s, day, month, year );
								livewood_turnover ( c, a, s, day, month, year );

								/* carbon use efficiency */
								carbon_use_efficiency ( c, height, dbh, age, species, day, month, year );

								/* biomass production efficiency */
								biomass_production_efficiency ( c, height, dbh, age, species, day, month, year );

								/* water use efficiency */
								water_use_efficiency  ( c, height, dbh, age, species, day, month, year );

                                // update canopy cover projection considering the
                                // number of living trees after natural mortality

                                  if ( c->doy == ( IS_LEAP_YEAR ( c->years[year].year ) ? 366 : 365) )
                                    {

                                     canopy_cover    ( c, height, dbh, age, species );

                                    }

								/* update Leaf Area Index */
								daily_lai             ( c, a, s );

								/* litter fluxes and pools */
								littering             ( c, s );

								/* tree level dendrometry */
								dendrometry_old       ( c, layer, height, dbh, age, species, year );

#if 0
    /*************************************** NATURAL REGENERATION PROCESS **************************************************/


                        if ( g_settings->Natural_regeneration) {

                            if ( c->doy == ( IS_LEAP_YEAR ( c->years[year].year ) ? 366 : 365) )

								{
                                   // Seeds production function
                                   Fruit_to_seeds_function ( c, a, s, day, month, year);
                                   //printf("Carbon fruit = %f\n", s->value[MAX_FRUIT_C]);
                                   //printf("Seed = %ld\n", s->counter[N_SEED]);

                                   } //end of leap year

                                     //Compute spring thermic sum before germination
                                     Thermic_sum_spring (c, meteo_daily, day, month, year);
                                     //printf("Thermicsum = %f\n", meteo_daily->spring_thermic_sum);

                                     if (year && meteo_daily->spring_thermic_sum >= s->value[GDD_SEED])
//if ( meteo_daily->spring_thermic_sum >= s->value[GDD_SEED])
                                        {
                                         //Seeds germination
                                         germination (c, meteo_daily, s, day, month, year);
                                         //printf("Seedlings = \t%d\n", s->counter[SEEDLINGS]);
                                          // ddalmo
                                           //printf("Sono in gertmination");
                                        }
                                        else
                                        {

                                        c->Seedlings_Number = 0;

                                    }

                                   //Compute winter soil mean temperature before germination
                                   //Soil_winter_temperature (c, meteo_daily, day, month, year);
                                   //printf("Winter soil temp = %f\n", meteo_daily->winter_soil);

                                   //Compute par that reach the soil in summer
                                   Seedling_soil_par (c, meteo_daily, day, month, year);

                                   //Compute air temperature in summer
                                   Seedling_temp (c, meteo_daily, day, month, year);

                                   // Establishment
                                   establishment (c, meteo_daily, s, day, month, year);

                                   //Mortality Establishment
                                   //seedlings_mortality (c, meteo_daily, s, day, month, year);
                                   // ddalmo test to see when we let the new stratum tu appear
                                  // if ( year==1 && c->doy == ( IS_LEAP_YEAR ( c->years[year].year ) ? 366 : 365) ) // questo funziona
                                    if ( year && c->doy == ( IS_LEAP_YEAR ( c->years[year].year ) ? 366 : 365) )  // dal secondo anno
				    {
                                  printf("in recr 4 \n") ;
                                  //if (!day && !month && year)
                                  //{
                                    if( !recruitment (c, s, day, month, year) )
                                     {
                                        logger_error(g_debug_log, "unable to add new regeneration class! (exit)\n");
                                        exit(1);
                                    }
                                    }
                                // }



                                  } /*end natural regeneration settings (on || off)*/


#endif

#//------------------------------------------------------------------------------------------------------------------------------------------------

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

                                     //printf("Carbon fruit = %f\n", s->value[MAX_FRUIT_C]);
                                     //printf("Seed = %ld\n", s->counter[N_SEED]);
                                     //printf("tank =   %ld\n", s->counter[TANK_SEEDS]);
                                     //printf("Seedlings = \t%d\n", s->counter[SEEDLINGS]);
                                     //printf("Seedlings pool = %d\n", s->counter[SEEDLINGS_POOL]);
                                     //printf("Seedlings surv =  \t%d\n", s->counter[SEEDLINGS_SURV]);

                                } /*end of year*/




                        /***************************************************************************************************/

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
								/* 6 */  //fixme if ( ! check_tree_class_nitrogen_mass_balance  ( c, layer, height, dbh, age, species ) ) return 0;
							}
							else
							{
								//FIXME here model should remove class just after have checked that the balances are closed
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


    /*************************************** NATURAL REGENERATION PROCESS **************************************************/


        if ( g_settings->Regeneration) {

            if ( c->doy == ( IS_LEAP_YEAR ( c->years[year].year ) ? 366 : 365) )

                {
                  /** Compute seeds production function **/
                  Fruit_to_seeds_function_npp ( c, a, s, day, month, year);
                  //printf("Carbon fruit = %f\n", s->value[MAX_FRUIT_C]);
                  //printf("Seed = %ld\n", s->counter[N_SEED]);

                    } //end of leap year

                      /** Compute spring thermic sum before germination **/
                      Thermic_sum_spring (c, meteo_daily, day, month, year);
                      //printf("Thermicsum = %f\n", meteo_daily->spring_thermic_sum);

                        if (year && meteo_daily->spring_thermic_sum >= s->value[GDD_SEED])

                         {
                           /** Compute seeds germination **/
                           germination (c, meteo_daily, s, day, month, year);
                           //printf("Seedlings = \t%ld\n", s->counter[SEEDLINGS]);

                            }
                            else
                            {
                             s->counter[SEEDLINGS_SURV] = 0;
                             //c->Seedlings_Number = 0;

                            }

                             /** Compute winter soil mean temperature before germination **/
                             //Soil_winter_temperature (c, meteo_daily, day, month, year);
                             //printf("Winter soil temp = %f\n", meteo_daily->winter_soil);

                              //Compute par that reach the soil in summer
                              Seedling_soil_par (c, meteo_daily, day, month, year);

                              //Compute air temperature in summer
                              Seedling_temp (c, meteo_daily, day, month, year);

                              /** Establishment appearing **/
                              establishment (c, meteo_daily, s, day, month, year);
                              //printf("Seedlings survived = %ld\n", s->counter[SEEDLINGS_SURV]);

                              //Mortality Establishment
                              //saplings_mortality (c, height, dbh, age, species);

                            if ( year ==1 && c->doy == ( IS_LEAP_YEAR ( c->years[year].year ) ? 366 : 365) )

                                    //if ( year == 1 && c->doy == 1 ) test Sapo
                                    {

                                     //if (!day && !month && year) test Sapo
                                     //{
                                     // Recruitment //
                                    if( !recruitment (c, s, day, month, year) )
                                    {
                                     logger_error(g_debug_log, "unable to add new regeneration class! (exit)\n");
                                      exit(1);
                                     }
                                //printf("Name species = %s\n", s->name);

                                //Mortality of saplings
                                //saplings_mortality (c, height, dbh, age, species);
                              }
                        // }


                    } /*end natural regeneration settings (on || off)*/


	/* ok */
	return 1;
}
/*****************************************************************************************************************/
/*****************************************************************************************************************/
/*****************************************************************************************************************/


