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
#include "regeneration.h"
#include "decomposition.h"
#include "littering.h"
#include "CN-balance.h"

extern logger_t* g_debug_log;
//extern soil_settings_t* g_soil_settings;
extern settings_t* g_settings;
extern dataset_t* g_dataset;

static const char sz_management[] = "TCN";


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

            //if ( (c->dos == 1) & ((MANAGEMENT_VAR == g_settings->management) || (MANAGEMENT_ON == g_settings->management)))
	  if ( (c->dos == 1))
	{
	   c->PRINT_MAN_HEADER = 0  ;    // serve per eventualmente stampare l header dell output file 

	   EOY_print_output_class_level_management(c, 0, 0, 0, 0, year);
	}

	c->PRINT_MAN_HEADER = 1 ; 


               // ddalmo at some point: let's perform the if statement here
               //if (MANAGEMENT_VAR == g_settings->management || MANAGEMENT_ON == g_settings->management)
               // {
		    /* management blocks */
	     	forest_management ( c, day, month, year );
             // }

	}

    // i.e. this condition might occurr when harvest happens, without replanting 
	if ( c->n_trees ==0 ) goto end_tree;


	/****************************************************************************/

	/* annual forest structure (only the year after the first) */
	if ( ( ! day && ! month && year ) && ( ! c->harvesting ) )
	{

  //printf(" IN TREE DAILY ENTERING ANNUAL FOREST STRUCTURE \n");

		annual_forest_structure ( c, year );


	}

	/* daily forest structure*/
	daily_forest_structure ( c,  meteo_daily, year);

	/* print forest cell data */
	print_daily_forest_data ( c );

	/* prephenology */
	prephenology ( c, meteo_daily, day, month, year );


	/***********************************************************************************************/

	/* sort heights in descending order */
#ifndef USE_NEW_OUTPUT
	qsort ( c->heights, c->heights_count, sizeof (height_t), sort_by_heights_desc );
#endif

 //printf(" c->tree_layers_count             = %d \n ", c->tree_layers_count);


	/* loop on each cell layers starting from highest to lower */
	for ( layer = c->tree_layers_count -1 ; layer >= 0; --layer )
	{
		/* assign shortcut */
		l = &m->cells[cell].tree_layers[layer];

		//printf(" IN TREE_DAILY   layer              = %d  \n",  layer );


		logger(g_debug_log,"*****************************************************************************\n"
				"                                layer = %d                                 \n"
				"*****************************************************************************\n",layer);

		l->layer_z = layer;

		/* loop on each heights starting from highest to lower */
		for ( height = 0; height <  c->heights_count; ++height )
		{

			//printf(" BEGINNING CYCLE height              = %d  \n", height);
		    //printf(" c->heights_count              = %d \n ", c->heights_count);
			//printf(" l->layer_z             = %d \n ", l->layer_z);
            //printf("in TREE DAILY  c->heights[height] %g,\n",c->heights[height].value);

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

				//printf(" TREE HEIGHT CLASS MATCH CELL LAYER \n");

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

							//printf(" IN TREE MODEL SPECIES  %s!!!\n", s->name);
							//printf(" IN TREE MODEL management type for the species  %c !!!\n", sz_management[s->management]);

							logger(g_debug_log,"*****************************************************************************\n"
									"*                              species = %s                         *\n"
									"*****************************************************************************\n", s->name);
							/**********************************/

                          // printf("in TREE DAILY  HEIGHT %g,\n",c->heights[height].value);
						//	printf("in TREE DAILY  AGE %d,\n",c->heights[height].dbhs[dbh].ages[age].value);
						   //printf("in TREE DAILY  DBH  %g,\n",c->heights[height].dbhs[dbh].value);
						  //  printf("in TREE DAILY  NTREE  %d,\n",s->counter[N_TREE]);

                             // printf("in TREE DAILY  CROWN_AREA_PROJ  %g,\n", s->value[CROWN_AREA_PROJ]);
                             


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

							//modifiers ( c, layer, height, dbh, age, species, meteo_daily, meteo_annual );
                                                        // 5p606
                                                         modifiers ( c, layer, height, dbh, age, species, meteo_daily, meteo_annual, month, day );

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

							/* allocate daily carbon */

							// carbon_allocation       ( c, a, s, day, month, year ); // old subroutine

							//note: this is basically the new function in version v.5.5
							carbon_allocation_new   ( c, a, s, day, month, year );

							/* allocate daily nitrogen */
							nitrogen_allocation     ( c, s, day, month, year );

							/* MORTALITY */
							/* Mortality based on growth efficiency */
							/* note: when it happens the overall class is removed */

                            //printf("BEFORE mortaliy GEff  height              = %d  \n", height);
							// printf("BEFORE mortaliy GEff  m->cells[cell].heights_count               = %d \n ", m->cells[cell].heights_count);
								 

							if ( ! growth_efficiency_mortality ( c, height, dbh, age, species ,year) )
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

								/* update Leaf Area Index */
								daily_lai             ( c, a, s );

								  //printf("IN TREE MODEL BEFORE DENDRO OLD \n");
                    
								/* tree level dendrometry */
								dendrometry_old       ( c, layer, height, dbh, age, species, year );

								
								/** IF END OF YEAR **/

								/* last day of the year */
								if ( c->doy == ( IS_LEAP_YEAR ( c->years[year].year ) ? 366 : 365) )
								{

                                     /* annual GROSS volume increment  */  // before mortality happen
									annual_gross_tree_increment ( c, height, dbh, age, species, year );


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
											age_mortality        ( c, height, dbh, age, species ,year);


							          		 /* Mortality based on stochasticity */


											stochastic_mortality ( c, height, dbh, age, species ,year);


											if (s->counter[DEAD_TREE] > 0 )
											{                             
												//new var to add 
												c->dead_tree_to_print     =       s->counter[DEAD_TREE];  
												c->dead_stem_b_to_print   =       s->counter[DEAD_TREE] * s->value[TREE_STEM_C]  ; 
												c->dead_branch_b_to_print  =      s->counter[DEAD_TREE] * s->value[TREE_BRANCH_C]  ; 
												c->dead_croot_b_to_print   =      s->counter[DEAD_TREE] * s->value[TREE_CROOT_C] ; 

												c->thinned_tree_to_print   =       0    ;     
												c->hwp_to_print            =       0.;
												c->thinned_branch_to_print =      0.  ; /*  stem volume removed (m3/ha/yr) */
												c->thinned_stem_to_print   =   0. ;
												c->thinned_stem2_to_print  =    0.;

												EOY_print_output_class_level_management(c, height, dbh, age, species, year);
											}

                                         }

                                    }
                                        else
                                    {         // with management ON or OFF always compute mortality
	
											/* Mortality based on tree Age (LPJ) */
											age_mortality        ( c, height, dbh, age, species ,year);


											/* Mortality based on stochasticity */
											stochastic_mortality ( c, height, dbh, age, species,year );

                                   			 // PRINT MORTALITY DATA 


											 if (s->counter[DEAD_TREE] > 0 )
			   								{  

												//new var to add 
												c->dead_tree_to_print     =       s->counter[DEAD_TREE];  
												c->dead_stem_b_to_print   =       s->counter[DEAD_TREE] * s->value[TREE_STEM_C]  ; 
												c->dead_branch_b_to_print  =      s->counter[DEAD_TREE] * s->value[TREE_BRANCH_C]  ; 
												c->dead_croot_b_to_print   =      s->counter[DEAD_TREE] * s->value[TREE_CROOT_C] ; 

												c->thinned_tree_to_print   =       0    ;     
												c->hwp_to_print            =       0.;
												c->thinned_branch_to_print =      0.  ; /*  stem volume removed (m3/ha/yr) */
												c->thinned_stem_to_print   =   0. ;
												c->thinned_stem2_to_print  =    0.;

												EOY_print_output_class_level_management(c, height, dbh, age, species, year);
		 									}
							      	}

                                    // update 

									s->value[STAND_BASAL_AREA]    = s->value[BASAL_AREA]    * s->counter[N_TREE];
	                                s->value[STAND_BASAL_AREA_m2] = s->value[BASAL_AREA_m2] * s->counter[N_TREE];


									// update canopy cover projection 
                                     
                                    canopy_cover    ( c, height, dbh, age, species );
                                    
									/* update Leaf Area Index */
									daily_lai             ( c, a, s );

									/* above ground-below ground stocks */
									abg_bgb_biomass ( c, height, dbh, age, species );

									/* annual branch and bark fraction */
									tree_branch_and_bark ( c, height, dbh, age, species );

									/* annual volume, MAI and CAI */
									annual_tree_increment ( c, height, dbh, age, species, year );
								}

								//printf("TREE DAILZ s->counter[N_TREE]           = %d \n",s->counter[N_TREE] );
                                //printf("TREE DAILZ s->value[sTEM_C]           = %g \n",s->value[STEM_C] );

                               // printf("IN TREE MODEL c->daily_leaf_to_litrC           = %g tC/cell/day\n",c->daily_leaf_to_litrC );
                               // printf("IN TREE MODEL BEFORE LITTERING \n");

								/* litter fluxes and pools */
								    littering             ( c, s );

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

                                // GROWTH EFFICIENCY MORTALITY 
								 
								//FIXME here model should remove class just after have checked that the balances are closed
								//so model has to include c fluxes that go out to litter and cwd

                                //printf("in mortaliy GEff  height              = %d  \n", height);
								 //printf("in mortaliy GEff  m->cells[cell].heights_count               = %d \n ", m->cells[cell].heights_count);
								//printf(" in mortaliy GEff c->tree_layers_count             = %d \n ", c->tree_layers_count);

                                 c->GREFFMORT_HAPPENS = 1 ;
								 //printf(" in mortaliy GEff c->GREFFMORT_HAPPENS            = %d \n ", c->GREFFMORT_HAPPENS);


								//printf("in mortaliy GEff  height              = %d  \n", height);
								 //printf("in mortaliy GEff  m->cells[cell].heights_count               = %d \n ", m->cells[cell].heights_count);
                                 // annual_forest_structure ( c, year );
                                 
								  
								if ( height >= m->cells[cell].heights_count ) 
								   { 
									printf("GO TO HEIGHT END  \n ") ;
									goto height_end; 
									}
								if ( dbh >= m->cells[cell].heights[height].dbhs_count ) goto dbh_end;
								if ( age >= m->cells[cell].heights[height].dbhs[dbh].ages_count ) goto age_end;
								//if ( species >= m->cells[cell].heights[height].dbhs[dbh].ages[age].species_count ) goto species_end;
							    
								 height = height -1 ;  // per poter procedere alle classi di altezze piu' basse 
								 
								 //c->cell_heights_count --;
								//l->layer_n_height_class =  l->layer_n_height_class -1 ;
								//l->layer_height_class_counter --;
								

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
		//printf(" FINITO GIRO SULLE ALTEZZE  \n");
		 //printf(" height              = %d  \n", height);
								 

		height_end:
		logger(g_debug_log, "****************END OF HEIGHT CLASS***************\n");
	}
	// printf(" IN TREE DAILY ENDE  c->tree_layers_count %d\n",c->tree_layers_count);
	logger(g_debug_log, "****************END OF LAYER CLASS***************\n");
	/* ok */

    //Compute par that reach the soil in summer
    Seedling_soil_par (c, meteo_daily, month, year);


    

    // Regeneration: check if conditions for regeneration sussist
	// 5.7 we consider regerenration as prescribed and with MAN = OFF 

    if ( c->doy == ( IS_LEAP_YEAR ( c->years[year].year ) ? 366 : 365) )

       	{

		 //printf(" END OF THE YEAR CALCULATING SEEDLINGD PAR! \n");


         double seedlings_par_cum ;
		 double seedlings_par_ave ;
		 double seedl_par_threshold;

		 seedl_par_threshold = 4.0 ;  // mmol etetc (intanto lo settiamo internamente)

         seedlings_par_cum = c->years[year].seedling_par ;

		 seedlings_par_ave = c->years[year].seedling_par/90. ;   // NOTE: in first approximation 90 days
         

         //printf(" seedlings_par_ave= %g \n",  seedlings_par_ave);

         if (seedlings_par_ave >= seedl_par_threshold) 
		 {
           c->seedl_reg = 1; 
         //printf(" CONDITION FOR REGENERATION IS ON!  c->seedl_reg= %d \n",   c->seedl_reg);


		 } else {

            c->seedl_reg = 0; 
		  	
		 }
		   //printf(" CONDITION FOR REGENERATION:  c->seedl_reg= %d \n",   c->seedl_reg);
		}
  
    

	end_tree:  
	return 1;
}

/*****************************************************************************************************************/
/*****************************************************************************************************************/
/*****************************************************************************************************************/


