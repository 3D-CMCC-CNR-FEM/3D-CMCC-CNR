/*management.c*/

/* includes */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <math.h>
#include <assert.h>
#include "new_forest_tree_class.h"
#include "g-function.h"
#include "management.h"
#include "constants.h"
#include "settings.h"
#include "common.h"
#include "logger.h"
#include "remove_tree_class.h"
#include "biomass.h"
#include "allometry.h"
#include "structure.h"
#include "initialization.h"

extern management_t* g_management;
extern settings_t* g_settings;
extern logger_t* g_debug_log;
extern dataset_t* g_dataset;

static int harvesting (cell_t *const c, const int height, const int dbh, const int age, const int species)
{
	int nat_man;   /* natural or managed mortality 0 = natural; 1 = managed */
	species_t *s;

	nat_man = 1;


	s = &c->heights[height].dbhs[dbh].ages[age].species[species];

	/* at the moment it considers a complete harvesting for all classes (if considered) */
	logger(g_debug_log, "\n** Management options: Harvesting ** \n");

	/* add harvested trees */
	s->counter[THINNED_TREE]  += s->counter[N_TREE];
	s->counter[HARVESTING_HAPPENS] = 1;

	/* update C and N biomass */
	tree_biomass_remove ( c, height, dbh, age, species, s->counter[N_TREE], nat_man );

        // update cell level C and N pools
    #if 0
	/*** update cell level carbon pools (tC/cell) ***/
	c->leaf_carbon              -= (s->value[C_LEAF_TO_LITR]   * 1e6 / g_settings->sizeCell);
	c->froot_carbon             -= (s->value[C_FROOT_TO_LITR]  * 1e6 / g_settings->sizeCell);
	c->stem_carbon              -= (s->value[C_STEM_TO_CWD]    * 1e6 / g_settings->sizeCell);
	c->branch_carbon            -= (s->value[C_BRANCH_TO_CWD]  * 1e6 / g_settings->sizeCell);
	c->croot_carbon             -= (s->value[C_CROOT_TO_CWD]   * 1e6 / g_settings->sizeCell);
	c->reserve_carbon           -= (s->value[C_RESERVE_TO_CWD] * 1e6 / g_settings->sizeCell);
	c->fruit_carbon             -= (s->value[C_FRUIT_TO_CWD]   * 1e6 / g_settings->sizeCell);

	/* check */
	CHECK_CONDITION ( c->leaf_carbon,    < , ZERO );
	CHECK_CONDITION ( c->froot_carbon,   < , ZERO );
	CHECK_CONDITION ( c->stem_carbon,    < , ZERO );
	CHECK_CONDITION ( c->branch_carbon,  < , ZERO );
	CHECK_CONDITION ( c->croot_carbon,   < , ZERO );
	CHECK_CONDITION ( c->reserve_carbon, < , ZERO );
	CHECK_CONDITION ( c->fruit_carbon,   < , ZERO );
   #endif
        //FIXME : update the cell level N pools too!

	/* remove completely all trees */

	return tree_class_remove (c, height, dbh, age, species );
}


int forest_management (cell_t *const c, const int day, const int month, const int year)
{

	//int layer;
	int height;
	int dbh;
	int age;
	int species;
        int THINNING_INDEX =0;

        //v5.6
        double harvested_carbon;
        double harvested_carbon_cum ;
        double harvested_volume ;
        double harvested_volume_cum ;

        int thinned_plants ;
        int row = 0;
        int year_dens_fin = 0;  // only used if MANAGEMENT == VAR

	height_t *h;
	dbh_t *d;
	age_t *a;
	species_t *s;


	// sort by above or below ?
	// ALESSIOR TO ALESSIOC FIXME
	// ALESSIOR commented on 10/07/2017
#ifndef USE_NEW_OUTPUT
	qsort ( c->heights, c->heights_count, sizeof (height_t),

			(THINNING_REGIME_ABOVE == g_settings->thinning_regime) ? sort_by_heights_asc : sort_by_heights_desc
	);
#endif

	/* loop on each heights starting from highest to lower */  // NOTE: this works with man = ON, monolayer and multilayer, with one specie for each layer.
                                                                   // in the latter case a list of replanting setting for each specie should be
                                                                   // provided (and the code should be changed
                                                                   // as currently only one replanted-specie info can be provided)
                                                                   // and with man = VAR && regeneration = OFF

        if (!g_settings->regeneration)
        {

         for ( height = c->heights_count -1 ; height >= 0; --height )
         {

		/* assign shortcut */
		h = &c->heights[height];

		/* loop on each dbh starting from highest to lower */
		for ( dbh = h->dbhs_count - 1; dbh >= 0; --dbh )
		{
			/* assign shortcut */
			d = &h->dbhs[dbh];

			/* loop on each age class */
			for ( age = d->ages_count - 1 ; age >= 0 ; --age )
			{
				/* assign shortcut */
				a = &d->ages[age];

				/* loop on each species class */
				for ( species = 0; species < a->species_count; ++species )
				{
					int flag = 0;

					/* assign shortcut */
					s = &a->species[species];

					/* initilIze management */
					s->counter[THINNING_HAPPENS] = 0;
					c->harvesting                = 0;

					if ( MANAGEMENT_ON == g_settings->management )
					{
						/* check at the beginning of simulation */
						if( ! year )
						{
							CHECK_CONDITION ( c->years[year].year, >, g_settings->year_start_management );
							CHECK_CONDITION ( (g_settings->year_start_management - g_settings->year_start), >, s->value[THINNING] );
						}

						if ( ( c->years[year].year == g_settings->year_start_management) ||

							(( c->years[year].year >= g_settings->year_start_management )
							&& ( s->value[THINNING] == s->counter[YEARS_THINNING] )) )
						{
							flag = 1;
						}
					}
					else if ( (MANAGEMENT_VAR == g_settings->management) || (MANAGEMENT_VAR1 == g_settings->management) )
					{
					        // compute last year of available stand density data
					        // it works with one layer/specie only

                            row = g_dataset->rows_count ;

                            year_dens_fin = g_dataset->rows[row-1].year_stand;

                            // if (c->years[year].year > year_dens_fin)
                            //{
						/* management forced by stand data */  // if data are available
						// if ( year )
						if ( year && c->years[year].year <= year_dens_fin)
						{
                         
							s->counter[THINNING_HAPPENS] = 1;

							prescribed_thinning ( c, height, dbh, age, species, c->years[year].year );
						}

                        // Management prescribed by external file: only with one class or combined with regeneration option

						if ( g_management && g_management->thinning_years_count )
						{

							int i;
							for ( i = 0; i < g_management->thinning_years_count; i++ )
							{
								if ( c->years[year].year == g_management->thinning_years[i] )
								{
									flag = 1;
                                    //s->value[THINNING_INDEX]= i;
                                    THINNING_INDEX= i;          // to set e.g. intensity of the thinning from external file
									break;
								}
							}
						}
					}

					/* thinning */
					if ( flag )
					{
						logger(g_debug_log,"**FOREST MANAGEMENT**\n");
						logger(g_debug_log,"**THINNING**\n");

						s->counter[THINNING_HAPPENS] = 1;

						  
						thinning ( c, height, dbh, age, species, year ,THINNING_INDEX);

						/* reset counter */
						s->counter[YEARS_THINNING] = 0;
					}

					/* increment counter */
					++s->counter[YEARS_THINNING];

					/*************************************** HARVESTING AND REPLANTING *************************************/

					flag = 0;
					if ( MANAGEMENT_ON == g_settings->management )
					{
						/* check */
						CHECK_CONDITION( s->counter[YEARS_THINNING], >, s->value[ROTATION] );

						/* check for harvesting */
						if ( a->value == s->value[ROTATION] )
						{
							flag = 1;
						}
					}
					else if ( MANAGEMENT_VAR == g_settings->management )
					{
						if ( g_management && g_management->harvesting_years_count )
						{
							int i;
							for ( i = 0; i < g_management->harvesting_years_count; i++ )
							{
								if ( c->years[year].year == g_management->harvesting_years[i] )
								{
									flag = 1;
									break;
								}
							}
						}
					}
					if ( flag )
					{
						int rsi;               /* replanted species index */

						logger(g_debug_log,"**FOREST MANAGEMENT**\n");
						logger(g_debug_log,"**HARVESTING**\n");

                                          	// august 2021, Note:
                                          	// the following part of code has been written - very likely - with the idea to provide replanted settings for
                                          	// different species, and then trough the loop identify rsi in a such way that the replanted specie
                                          	// corresponds to the current specie in the layer. However, this does not allow the change in specie after clear cut

#if 0
						/* get replanted_species_index */
						for ( rsi = 0; rsi < g_settings->replanted_count; rsi++ )
						{
							if ( ! string_compare_i(c->heights[height].dbhs[dbh].ages[age].species[species].name
														, g_settings->replanted[rsi].species) )
							{
								/* index found */
								break;
							}
						}


						/* species found ? */
						CHECK_CONDITION(rsi, ==, g_settings->replanted_count );
#endif
                                                thinned_plants =       s->counter[N_TREE] ; //save number of harvested trees

						/* remove tree class */
						//if (  ! harvesting ( c, height, dbh, age, species, rsi ) ) //ddamo: no need to pass rsi
                                                if (  ! harvesting ( c, height, dbh, age, species ) )
						{
							logger_error(g_debug_log, "unable to harvesting! (exit)\n");
							exit(1);
						}

						/* note: RESET c->dos */
						c->dos = 0;

                                                // august 2021

                                                // with 'harvest' we remove the class, with replanting we create a new class in the same
                                                // virtual space (as result the value of e.g. s->value[C_HWP] is simply set to 0 once the replanted class
                                                // is added. For this reason We save hence some key variables belonging to the old removed layer
                                                // and we save it in the new replanted layer object.
                                                // no tests have been yet performed With multilayer or multispecie

                                                harvested_carbon =      s->value[C_HWP]     ;      /* woody biomass removed (tC/ha/yr) */
                                                harvested_carbon_cum =  s->value[CUM_C_HWP] ;
                                                harvested_volume =      s->value[VOLUME_HWP]     ; /*  stem volume removed (m3/ha/yr) */
                                                harvested_volume_cum =  s->value[CUM_VOLUME_HWP] ;


                                               // v5.6
                                            	// As currently set, we need to assume that the clear cut is performed in one class only, and
                                            	// only one new class is added (thus in the setting file, only one replanting dataset is provided)
                                            	// The replanted specie can be different from the removed old layer.
                                            	// rsi is simply g_settings->replanted_count -1

                                            	rsi =  g_settings->replanted_count -1 ;

						/* check that all mandatory variables are filled */
						CHECK_CONDITION (g_settings->replanted[rsi].n_tree, <, ZERO);
						CHECK_CONDITION (g_settings->replanted[rsi].height, <, 1.3);
						CHECK_CONDITION (g_settings->replanted[rsi].avdbh,  <, ZERO);
						CHECK_CONDITION (g_settings->replanted[rsi].age,    <, ZERO);

						/* re-planting tree class */
						if( g_settings->replanted[rsi].n_tree )
						  {
							if ( ! add_tree_class_for_replanting( c , day, month, year, rsi ) )
							{
								logger_error(g_debug_log, "unable to add new replanted class! (exit)\n");
								exit(1);
							}

							// indexes
							h = &c->heights[height];
							d = &h->dbhs[dbh];
							a = &d->ages[age];
							s = &a->species[species];

							/* reset years_for_thinning */
							s->counter[YEARS_THINNING] = 1;

                                                  	// save in the new class object the following variables wich refer to the
                                                  	// removed old class.

                                                  	s->value[C_HWP]          =  harvested_carbon      ;
                                                  	s->value[CUM_C_HWP]      =  harvested_carbon_cum  ;
                                                  	s->value[VOLUME_HWP]     =  harvested_volume      ;
                                                  	s->value[CUM_VOLUME_HWP] =  harvested_volume_cum  ;
                                                  	s->counter[THINNED_TREE] =  thinned_plants        ;

                                                }

                                                s->counter[THINNING_HAPPENS] = 1;
						 c->harvesting                = 1;

			                }

				}
			}
		}
	 }

        }
        else  // SHELTERWOOD CASE WITH PRESCRIBED REGENERATION : only in combination with MAN = VAR
        {

       // 5p6 regeneration+MAN=VAr option works only on a monoclass forest
       // Management is performed in the dominant layer only.

         height = c->heights_count -1 ;

		/* assign shortcut */
		h = &c->heights[height];

		/* loop on each dbh starting from highest to lower */
		for ( dbh = h->dbhs_count - 1; dbh >= 0; --dbh )
		{
			/* assign shortcut */
			d = &h->dbhs[dbh];

			/* loop on each age class */
			for ( age = d->ages_count - 1 ; age >= 0 ; --age )
			{
				/* assign shortcut */
				a = &d->ages[age];

				/* loop on each species class */
				for ( species = 0; species < a->species_count; ++species )
				{
					int flag = 0;

					/* assign shortcut */
					s = &a->species[species];

					/* initilIze management */
					s->counter[THINNING_HAPPENS] = 0;
					c->harvesting                = 0;

					//MANAGEMENT = VAR when REGENERATION = ON

				        /* management forced by stand data */
					if ( year )
					 {
						 s->counter[THINNING_HAPPENS] = 1; // FIXME the removed trees as prescribed from stand_file.txt
						// are then saved as HWP, to be consistent with the option MAN=VAR & REG=OFF
						// however the observed density in the stand could be due to mortality and not thinning
                        //printf(" SONO QUA MANAGEMENT.c s->counter[THINNING_HAPPENS]  **************************    = %d \n", s->counter[THINNING_HAPPENS]);

						prescribed_thinning ( c, height, dbh, age, species, c->years[year].year );
					  }

                                        // management prescribed by external file
					if ( g_management && g_management->thinning_years_count )
					 {

					   int i;
						for ( i = 0; i < g_management->thinning_years_count; i++ )
						 {
						   if ( c->years[year].year == g_management->thinning_years[i] )
						     {
							flag = 1;
                                                        //s->value[THINNING_INDEX]= i;
                                                        THINNING_INDEX= i;

							break;
						     }
						 }
					}


					/* thinning */
					if ( flag )
					{
						logger(g_debug_log,"**FOREST MANAGEMENT**\n");
						logger(g_debug_log,"**THINNING**\n");

						s->counter[THINNING_HAPPENS] = 1;

						thinning ( c, height, dbh, age, species, year ,THINNING_INDEX);

						/* reset counter */
						s->counter[YEARS_THINNING] = 0;
					}

					/* increment counter */
					++s->counter[YEARS_THINNING];

					/*************************************** HARVESTING *************************************/

					flag = 0;

					if ( g_management && g_management->harvesting_years_count )
					{
						int i;
						for ( i = 0; i < g_management->harvesting_years_count; i++ )
						{
							if ( c->years[year].year == g_management->harvesting_years[i] )
							{
								flag = 1;
								break;
							}
						}
					}

					if ( flag )
					{
						logger(g_debug_log,"**FOREST MANAGEMENT**\n");
						logger(g_debug_log,"**HARVESTING AFTER ADDED REGENERATION LAYER**\n");


                                                thinned_plants =       s->counter[N_TREE] ; //save number of harvested trees

						/* remove tree class */

                                                if (  ! harvesting ( c, height, dbh, age, species ) )
						{
							logger_error(g_debug_log, "unable to harvesting! (exit)\n");
							exit(1);
						}

						 /* note: RESET c->dos */
						 c->dos = 0;    //comment: however, when there is a simulated regeneration, is not stricly
                                                               //speaking starting from scratch


                                                 // save data of the removed layer


                                                 harvested_carbon =      s->value[C_HWP]     ;      /* woody biomass removed (tC/ha/yr) */
                                                 harvested_carbon_cum =  s->value[CUM_C_HWP] ;
                                                 harvested_volume =      s->value[VOLUME_HWP]     ; /*  stem volume removed (m3/ha/yr) */
                                                 harvested_volume_cum =  s->value[CUM_VOLUME_HWP] ;

                                                 annual_forest_structure ( c, year );

                                                 height = c->heights_count -  1;         //a questo punto dovrei avere solo un layer, quello di regen che diventa0
	                                         dbh = c->heights[height].dbhs_count - 1;
	                                         age = c->heights[height].dbhs[dbh].ages_count - 1;
	                                         species = c->heights[height].dbhs[dbh].ages[age].species_count - 1;


                                                 // indexes of the remaining layer after harves, i.e. the regeneration layer

					                         h = &c->heights[height];
						                     d = &h->dbhs[dbh];
						                     a = &d->ages[age];
					                         s = &a->species[species];

                                                 // save in the class object of the regeneration layer
                                                 // the data of the removed old class.

                                                 s->value[C_HWP]          =  harvested_carbon      ;
                                                 s->value[CUM_C_HWP]      =  harvested_carbon_cum  ;
                                                 s->value[VOLUME_HWP]     =  harvested_volume      ;
                                                 s->value[CUM_VOLUME_HWP] =  harvested_volume_cum  ;
                                                 s->counter[THINNED_TREE] =  thinned_plants        ;

                                                 s->counter[THINNING_HAPPENS] = 1;
						                         c->harvesting                = 1;
					}

                                        /*************************************** PRESCRIBED REGENERATION *************************************/

                                        //v. 5p6

                                        //if( g_settings->regeneration && (MANAGEMENT_VAR == g_settings->management))
                                        //{

                                        flag = 0;
                                        if ( g_management->regeneration_years_count )
						{
							int i;
							for ( i = 0; i < g_management->regeneration_years_count; i++ )
							{
								if ( c->years[year].year == g_management->regeneration_years[i] )
								{
									flag = 1;
									break;
								}
							}
						}


                                          if ( flag )   //add layer for regeneration
                                          {

                                          /* check that all mandatory variables are filled */

						CHECK_CONDITION (g_settings->regeneration_n_tree, <, ZERO);
						CHECK_CONDITION (g_settings->regeneration_height, <, 1.3);
						CHECK_CONDITION (g_settings->regeneration_avdbh,  <, ZERO);
						CHECK_CONDITION (g_settings->regeneration_age,    <, ZERO);

                                          /* mimic natural regeneration with a re-planting tree class */

	                                        if ( ! add_tree_class_for_replanting_reg( c , day, month, year) )
						{
							logger_error(g_debug_log, "unable to add new replanted class! (exit)\n");
							exit(1);
						}

						// indexes
						h = &c->heights[height];  // height should be = 0
						d = &h->dbhs[dbh];
						a = &d->ages[age];
						s = &a->species[species];

                                                /* reset years_for_thinning */
						s->counter[YEARS_THINNING] =   1;

 						c->harvesting              =   1;    //  in this way the annual_structure function is not
                                                                                    // called in the tree_model.c

                                          }
                                        //}

				}
			}
		}

        }
	return 0;

}

/*****************************************************************************************************************************************/

void thinning (cell_t *const c, const int height, const int dbh, const int age, const int species, const int year, const int THINNING_INDEX)
{
	int trees_to_remove = 0;
	int nat_man;                        /* natural or managed mortality 0 = natural; 1 = managed */
        int thinning_intensity_prescribed;

        int index= 0 ;


	species_t *s;

	nat_man = 1;

	s = &c->heights[height].dbhs[dbh].ages[age].species[species];

	//printf("** IN THINNING  **\n");
					
	/* thinning function based on basal area */

	//TODO
	if (THINNING_REGIME_ABOVE == s->value[THINNING_REGIME] )
	{

	}
	else
	{

	}

	logger(g_debug_log, "** Management options: Thinning ** \n");

	/* BAU MANAGEMENT */

        // currently management_type is not set, hence it is automatically equal to 0

	if ( 0 == g_settings->management_type )
	{
               // 5p6: added the possibility to prescribe thinning intensity for each thinning event via the management file.
               if ( g_management && g_management->thinning_intensity_count )
		{
                                      index =	 THINNING_INDEX;

                                      thinning_intensity_prescribed= g_management->thinning_intensity[index];


                /* compute integer number of trees to remove */
		trees_to_remove = ROUND((thinning_intensity_prescribed / 100. ) * s->counter[N_TREE]);

               }
               else
               {
                // thinning intensity constant for each thinning event as set in the specie-file
		/* compute integer number of trees to remove */
		trees_to_remove = ROUND((s->value[THINNING_INTENSITY] / 100. ) * s->counter[N_TREE]);
               }
	}
	else
	{
		/* LOCAL MANAGEMENT */
		if ( ! IS_INVALID_VALUE(s->value[MINSTOCKGROW]) && ! IS_INVALID_VALUE(s->value[VOLUME]) )
		{
			int perc;
			double v;

			perc = 0;
			v = (s->value[VOLUME] - s->value[MINSTOCKGROW]) * 100. / s->value[MINSTOCKGROW];

			if ( v >= 80 )
			{
				perc = (int)s->value[THINNING_80];
			}
			else if ( (v >= 60) && (v < 80) )
			{
				perc = (int)s->value[THINNING_60_80];
			}
			else if ( (v >= 40) && (v < 60) )
			{
				perc = (int)s->value[THINNING_40_60];
			}
			else if ( (v >= 20) && (v < 40) )
			{
				perc = (int)s->value[THINNING_40_20];
			}
			else
			{
				perc = (int)s->value[THINNING_0];
			}

			if ( IS_INVALID_VALUE(perc) || ! perc )
			{
				trees_to_remove = 0;
			}
			else
			{
				trees_to_remove = s->counter[N_TREE] * perc / 100.;
			}
		}
	}

	/* added thinned trees */
	s->counter[THINNED_TREE] += trees_to_remove;

	if ( trees_to_remove < s->counter[N_TREE] )
	{
		/* update C and N biomass */   //FIXME: no update of the stand data is performed in tree_biomass_remove!! The update is done at
		                               // the end of the day, after GPP and NPP are computed.
		tree_biomass_remove ( c, height, dbh, age, species, trees_to_remove, nat_man );

		/* remove trees */
		s->counter[N_TREE] -= trees_to_remove;


	}
	else   // comment: this should not happen! anyhow we should set a minimum we have to leave in the stand.
	{
		/* update C and N biomass */
		tree_biomass_remove ( c, height, dbh, age, species, s->counter[N_TREE], nat_man );

		/* remove completely all trees */
		tree_class_remove   (c, height, dbh, age, species );
	}

	/* check */
	CHECK_CONDITION(s->counter[N_TREE], <, ZERO );

	/*********************************************************************************************************************************************************************/
	/* update stand trees */
	c->n_trees          -= trees_to_remove;
	c->annual_dead_tree += trees_to_remove;

}
/*****************************************************************************************************************************************/

void prescribed_thinning (cell_t *const c, const int height, const int dbh, const int age, const int species, const int year)
{
	int row;
	int nat_man;   /* natural or managed mortality 0 = natural; 1 = managed */

	nat_man = 1;

	              // comment: note that it could be a simple prescription of observed density rather then thinning for HPW removal.
	              // this means, that the removed trees should not be counted as HPW. (and nat_man should be 0! FIXME)

	assert(g_dataset);

    printf("** IN PRESCRIBED_THINNING **\n");

	for ( row = 0; row < g_dataset->rows_count; ++row )
	{
		if ( year == g_dataset->rows[row].year_stand )
		{
			if ( ! string_compare_i(g_dataset->rows[row].species, c->heights[height].dbhs[dbh].ages[age].species[species].name) )
			{
				if ( g_dataset->rows[row].n != c->heights[height].dbhs[dbh].ages[age].species[species].counter[N_TREE] )
				{
					int tree_remove;

					CHECK_CONDITION(c->heights[height].dbhs[dbh].ages[age].species[species].counter[N_TREE], <, g_dataset->rows[row].n);

					/* compute number of tree to remove */
					tree_remove = c->heights[height].dbhs[dbh].ages[age].species[species].counter[N_TREE] - g_dataset->rows[row].n;
                     
					//printf(" SONO QUA MANAGEMENT.c tree_remove]  = %d \n", tree_remove);
                    //printf(" SONO QUA MANAGEMENT.c g_dataset->rows[row].n  = %d \n", g_dataset->rows[row].n);

					/* added thinned trees */
					c->heights[height].dbhs[dbh].ages[age].species[species].counter[THINNED_TREE] += tree_remove;

					logger(g_debug_log, "\n** Management options: Prescribed Thinning **\n");

					tree_biomass_remove(c, height, dbh, age, species, tree_remove, nat_man);

					c->heights[height].dbhs[dbh].ages[age].species[species].counter[N_TREE] = g_dataset->rows[row].n;

                    /* update stand trees */
	                c->n_trees           -= tree_remove;
	                c->annual_dead_tree += tree_remove;

					// added on 7 September 2017
					// with 'var1' model gets also changes in dbh and height
					// on september 2024 it has been removed the option VAR1 as prescribing DBH 
					// would imply correcting the allocation of carbon from and to stem.
					if ( MANAGEMENT_VAR1 == g_settings->management )
					{
						c->heights[height].value = g_dataset->rows[row].height;
						c->heights[height].dbhs[dbh].value = g_dataset->rows[row].avdbh;
						c->heights[height].dbhs[dbh].ages[age].value = g_dataset->rows[row].age;

						/* initialize power function */
						allometry_power_function                ( c );

						/* initialize carbon pool fraction */
						carbon_pool_fraction                    ( c );

						/* initialize forest structure */
						annual_forest_structure                 ( c, year );

						/* initialize class carbon pools */
						initialization_forest_class_C           ( c, height, dbh, age, species );

						/* initialize cell carbon pools */
						initialization_forest_cell_C            ( c, height, dbh, age, species );

						/* initialize class nitrogen pools */
						initialization_forest_class_N           ( c, height, dbh, age, species );

						/* initialize cell nitrogen pools */
						initialization_forest_cell_N            ( c, height, dbh, age, species );
					}
				}
			}
		}
	}
}


void management_free(management_t* p)
{
	if ( p )
	{
		if ( p->harvesting_years_count )
		{
			free(p->harvesting_years);
		}

		if  ( p->thinning_years_count )
		{
			free(p->thinning_years);
		}

                if  ( p->thinning_intensity_count )
		{
			free(p->thinning_intensity);
		}


		free(p);
	}
}

#if 0  //the new function is reported below
management_t* management_load(const char* const filename)
{
#define BUFFER_SIZE 512
	char* token;
	char* p;
	char buffer[BUFFER_SIZE]; // should be enough
	int** p_years;
	int* p_years_count;
	int thinning_flag = 1;
	FILE *f;
	management_t* management;

	const char sz_delimiters[] = " ,\t\r\n";
	const char sz_harvesting[] = "harvesting";
	const char sz_thinning[] = "thinning";

        //const char sz_thinning_int[] = "thinning_int";   //

	assert(filename);

	f = fopen(filename, "r");
	if ( ! f )
	{
		logger_error(g_debug_log, "file not found\n");
		return NULL;
	}

	management = malloc(sizeof*management);
	if ( ! management )
	{
		logger_error(g_debug_log, "out of memory\n");
		fclose(f);
		return NULL;
	}
	memset(management, 0, sizeof(management_t));

	while ( fgets(buffer, BUFFER_SIZE, f) )
	{
		// skip empty lines
		if ( '\0' == buffer[0] )
		{
			continue;
		}

		token = string_tokenizer(buffer, sz_delimiters, &p);
		if ( ! token ) continue;

		if ( ! string_compare_i(token, sz_harvesting) )
		{
			// parse harvesting
			thinning_flag = 0;
		}
		else if ( ! string_compare_i(token, sz_thinning) )
		{
			// parse thinning
			thinning_flag = 1;
		}
		else
		{
			// do nothing
			continue;
		}

		// set pointers
		if ( thinning_flag )
		{
			p_years = &management->thinning_years;
			p_years_count = &management->thinning_years_count;
		}
		else
		{
			p_years = &management->harvesting_years;
			p_years_count = &management->harvesting_years_count;
		}

		if ( *p_years_count )
		{
			logger_error(g_debug_log, "too many rows of %s specified\n"
						, thinning_flag ? sz_thinning : sz_harvesting
			);
			management_free(management);
			fclose(f);
			return NULL;
		}

		while ( ( token = string_tokenizer(NULL, sz_delimiters, &p ) ) )
		{
			int err;
			int* int_no_leak;
			int year = convert_string_to_int(token, &err);
			if ( err )
			{
				logger_error(g_debug_log, "unable to convert year '%s' for %s\n"
							, token
							, thinning_flag ? sz_thinning : sz_harvesting
				);
				management_free(management);
				fclose(f);
				return NULL;
			}

			int_no_leak = realloc(*p_years, (*p_years_count+1)*sizeof*int_no_leak);
			if ( ! int_no_leak )
			{
				logger_error(g_debug_log, "out of memory\n");
				management_free(management);
				fclose(f);
				return NULL;

			}
			*p_years = int_no_leak;
			(*p_years)[*p_years_count] = year;
			++*p_years_count;
		}
	}
	fclose(f);

	if ( ! management->harvesting_years_count && ! management->thinning_years_count )
	{
		logger_error(g_debug_log, "file is empty ?\n");
		management_free(management);
		management = NULL;
	}

	return management;
#undef BUFFER_SIZE
}
#endif

// New load function to load management external information of harvest,thinning and thinning intensity

management_t* management_load_presc(const char* const filename)
{
#define BUFFER_SIZE 512
	char* token;
	char* p;
	char buffer[BUFFER_SIZE]; // should be enough
	int** p_years;
	int* p_years_count;
        int** p_years_int;        // for each year the intensity of removal is provided
        int* p_years_int_count;
        int** p_years_reg;        // the year when regeneration has to happened is provided
        int* p_years_reg_count;

	int thinning_flag     = 0;
        int harvesting_flag   = 0;
        int thinning_int_flag = 0;
        int regeneration_flag = 0;

        FILE *f;
	management_t* management;

	const char sz_delimiters[] = " ,\t\r\n";
	const char sz_harvesting[] = "harvesting";
	const char sz_thinning[]   = "thinning";
        const char sz_thinning_int[] = "thinning_int";
        const char sz_regeneration[] = "regeneration";   //august 2021 :add also years where a new regeneration layer can be added (shelterwood cut)

	assert(filename);

	f = fopen(filename, "r");
	if ( ! f )
	{
		logger_error(g_debug_log, "file not found\n");
		return NULL;
	}

	management = malloc(sizeof*management);
	if ( ! management )
	{
		logger_error(g_debug_log, "out of memory\n");
		fclose(f);
		return NULL;
	}
	memset(management, 0, sizeof(management_t));

	while ( fgets(buffer, BUFFER_SIZE, f) )
	{
		// skip empty lines
		if ( '\0' == buffer[0] )
		{
			continue;
		}

		token = string_tokenizer(buffer, sz_delimiters, &p);
		if ( ! token ) continue;

		if ( ! string_compare_i(token, sz_harvesting) )
		{
			// parse harvesting
                        harvesting_flag = 1;
			thinning_flag = 0;
                        thinning_int_flag = 0;
                        regeneration_flag = 0;
		}
		else if ( ! string_compare_i(token, sz_thinning) )
		{
			// parse thinning
                        harvesting_flag = 0;
			thinning_flag = 1;
                        thinning_int_flag = 0;
                        regeneration_flag = 0;
		}
                else if ( ! string_compare_i(token, sz_thinning_int) )  //5p6
		{
			// parse thinning_intensity
                        harvesting_flag = 0;
			thinning_int_flag = 1;
                        thinning_flag = 0;
                        regeneration_flag = 0;
		}
                else if ( ! string_compare_i(token, sz_regeneration) )  //5p6
		{
			// parse regeneration
                        harvesting_flag = 0;
			thinning_int_flag = 0;
                        thinning_flag = 0;
                        regeneration_flag = 1;
		}
		else
		{
			// do nothing
			continue;
		}

                //if (!thinning_int_flag)
                //{
		// set pointers
		if ( thinning_flag )
		{
			p_years = &management->thinning_years;
			p_years_count = &management->thinning_years_count;
		}
		else if ( harvesting_flag )
		{
			p_years = &management->harvesting_years;
			p_years_count = &management->harvesting_years_count;
		}
                else if ( thinning_int_flag  )
		{
			p_years = &management->thinning_intensity;
			p_years_count = &management->thinning_intensity_count;
		}
                else if ( regeneration_flag )
		{
			p_years = &management->regeneration_years;
			p_years_count = &management->regeneration_years_count;
		}
                else
		{
			// do nothing
			continue;
		}

		if ( *p_years_count )
		{
			logger_error(g_debug_log, "too many rows of %s specified\n"
						, thinning_flag ? sz_thinning : sz_harvesting
			);
			management_free(management);
			fclose(f);
			return NULL;
		}

		while ( ( token = string_tokenizer(NULL, sz_delimiters, &p ) ) )
		{
			int err;
			int* int_no_leak;
			int year = convert_string_to_int(token, &err);
			if ( err )
			{
				logger_error(g_debug_log, "unable to convert year '%s' for %s\n"
							, token
							, thinning_flag ? sz_thinning : sz_harvesting
				);
				management_free(management);
				fclose(f);
				return NULL;
			}

			int_no_leak = realloc(*p_years, (*p_years_count+1)*sizeof*int_no_leak);
			if ( ! int_no_leak )
			{
				logger_error(g_debug_log, "out of memory\n");
				management_free(management);
				fclose(f);
				return NULL;

			}
			*p_years = int_no_leak;
			(*p_years)[*p_years_count] = year;
			++*p_years_count;
		}







	}

	fclose(f);

        //DATA INPUT CHECK!

	if ( ! management->harvesting_years_count && ! management->thinning_years_count)
	{
		logger_error(g_debug_log, "file is empty? Or missing data?(e.g. only thinning intensity/regeneration years data have been provided)\n");
		management_free(management);
		management = NULL;
	}

        if ( ! management->thinning_years_count && management->thinning_intensity_count)
	{
		logger_error(g_debug_log, "Missing prescribed years for thinning while thinning intensity data is provided\n");
		management_free(management);
		management = NULL;
	}

        if ( g_settings->regeneration && ! management->regeneration_years_count)
	{
		logger_error(g_debug_log, "Missing prescribed years for regeneration while regeneration= ON\n");
		management_free(management);
		management = NULL;
	}


	return management;
#undef BUFFER_SIZE
}


