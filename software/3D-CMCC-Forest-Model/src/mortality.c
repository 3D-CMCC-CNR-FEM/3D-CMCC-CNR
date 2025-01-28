/*mortality.c*/

/* includes */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "matrix.h"
#include "constants.h"
#include "common.h"
#include "g-function.h"
#include "settings.h"
#include "logger.h"
#include "mortality.h"
#include "tree_model.h"
#include "remove_tree_class.h"
#include "biomass.h"
#include "canopy_cover.h"
#include "allometry.h"
#include "littering.h"

extern settings_t* g_settings;
extern logger_t* g_debug_log;

static const char sz_management[] = "TCN";


/******************************************************************************************************************/
/******************************************************************************************************************/

void self_thinning_mortality_new(cell_t *const c, const int layer, const int year , const int  tree_remove_st)
{
	
	// ddalmo oct.2024 
	
	int height;
	int dbh;
	int age;
	int species;
	int deadtree;
	int livetree;
	int nat_man;   /* natural or managed mortality 0 = natural; 1 = managed */

	int tree_to_rm2 = 0 ;

	height_t *h;
	dbh_t *d;
	age_t *a;
	species_t *s;

	

	nat_man = 0;

	// l->tree_remove_crowded 

    //printf("** IN SELF THINNING MORTALITY **\n");

	/* "First, large plants suppress small plants.
	 * The result is a “hierarchy of dominance and suppression” in which the smaller plants
	 * are at an accumulating disadvantage and finally die".
	 * Westoby, 1984, Advances in Ecological Research */

	/* NOTE: it mainly follows rationale of TREEDYN 3 "Crowding Competition mortality function",
	 * Bossel, 1996, Ecological Modelling (eq. 29) */

	/* This is a sort of self-thinning rules not based on biomass but based on canopy cover */

	/* the model makes die trees of the lower height class for that layer because
	it passes through the function sort_by_height_desc the height classes starting from the lowest */

	logger(g_debug_log, "\n\n*****SELF THINNING MORTALITY for layer %d*****\n", layer);

	//printf("in SELFTHINNING heights_count             = %d  \n",  c->heights_count);
	 //	 printf("  IN SELF THINNING tree_remove_st %d!!!\n", tree_remove_st);
		//printf("in SELFTHINNING BEFORE REMOVING CLASS heights_count             = %d  \n",  c->heights_count);



	for ( height = c->heights_count - 1; height >= 0; --height )
	{
		/* first of all sort by descending height */
		/* model makes die lower height-class in each layer first */

#ifndef USE_NEW_OUTPUT
			qsort (c->heights, c->heights_count, sizeof (height_t), sort_by_heights_desc);
#endif

		h = &c->heights[height];

			//printf(" in self thin c->heights[height] %f!!!\n",c->heights[height].value);
		// ddalmo : added, as the selfthinning bc of crown competition is whitin the same layer...

		if( layer == c->heights[height].height_z )
			{


			for ( dbh = h->dbhs_count - 1; dbh >= 0; --dbh)
			{
				d = &c->heights[height].dbhs[dbh];

				// Skip selt-thinning tree removal is large trees

				if (d >= self_thinning_treshold)  goto next_dbh;


				for ( age = 0; age < d->ages_count ; ++age )
				{
					a = &c->heights[height].dbhs[dbh].ages[age];

					for ( species = 0; species < a->species_count; ++species )
					{
						s = &c->heights[height].dbhs[dbh].ages[age].species[species];

						// note: for how we work, for the lowest height in the layer we will be having only one class.
						// 
						logger(g_debug_log, "MORTALITY BASED ON HIGH CANOPY COVER height %f species %s dbh %f !!!\n", h->value, s->name, d->value);

						/* initialize variables */
						livetree = s->counter[N_TREE];
						//deadtree = 0;
						// we know that usually starting from the last height class, there is only one DBH and one age...
					
						// for forest navigator test: we use a relaxation factor. so that we remove smoothly the trees
						// after some test, with 10% we have an exponential decrease of living trees

						tree_to_rm2 = round(tree_remove_st*0.1)  ;
						deadtree = tree_to_rm2;

						//deadtree = tree_remove_st;
					
					
						livetree -= deadtree;
		
						// check if this is larger than the number of trees in the class
						// if so, remove the entire class and go on.
						// 
						//
						// directly remove the plants, and we should call then 
						
						// livetree = max(0, livetree) ;

						if ( livetree > 0 )
						{

                            /* update at cell level */
							c->daily_dead_tree += deadtree;

							/* update layer trees */
							c->tree_layers[layer].layer_n_trees -=deadtree;

							/* update layer density */
							c->tree_layers[layer].layer_density = c->tree_layers[layer].layer_n_trees / g_settings->sizeCell;

							/* update layer cover proj */
							c->tree_layers[layer].layer_cover_proj  -= s->value[CANOPY_COVER_PROJ];


							/* recompute dbhdc eff */
							dbhdc_function         (c, layer, height, dbh, age, species, year);

							/* recompute crown area */
							crown_allometry        ( c, height, dbh, age, species );

							/* recompute canopy cover with */
							s->value[CANOPY_COVER_PROJ] = s->value[CROWN_AREA_PROJ] * livetree / g_settings->sizeCell;

							/* check for recompued canopy cover */
							c->tree_layers[layer].layer_cover_proj += s->value[CANOPY_COVER_PROJ];


							//printf(" SONO THIN PRIMA DI TREE BIOMASS REMOVE \n");

							//FIXME
							/* remove dead C and N biomass */
							tree_biomass_remove ( c, height, dbh, age, species, deadtree, nat_man );

							/* update live and dead tree (do not move above) */
							s->counter[DEAD_TREE] += deadtree;
							s->counter[N_TREE]     = livetree;
							c->n_trees            -= deadtree;

							deadtree = 0;
							livetree = 0;

							/* check */
							CHECK_CONDITION( s->value[CANOPY_COVER_PROJ] ,  > , s->value[CANOPY_COVER_PROJ] + eps );
							goto height_end; 

						} else {   // removing more trees than the ones in the class


                     

							deadtree = s->counter[N_TREE];
							// remove the ENTIRE CLASS; 
							//the issue is, go on, 

							/* update at cell level */
							c->daily_dead_tree +=  s->counter[N_TREE];

							/* update layer trees */
							c->tree_layers[layer].layer_n_trees -= s->counter[N_TREE];


							/* update layer density */
							c->tree_layers[layer].layer_density = c->tree_layers[layer].layer_n_trees / g_settings->sizeCell;

							/* update layer cover proj */
							c->tree_layers[layer].layer_cover_proj  -= s->value[CANOPY_COVER_PROJ];



											/* reset to zero n_trees */
							s->counter[DEAD_TREE] += s->counter[N_TREE];
							s->counter[N_TREE]    = 0;
							c->n_trees -= deadtree;

							// jan 2025 save mortality data / save info for printing
							// if entire class is remove, call here the output print function
							// save info for printing (auxiliary variables)

							//new var to add 
							c->dead_tree_to_print     =       deadtree ;  
							c->dead_stem_b_to_print   =       deadtree * s->value[TREE_STEM_C]  ; 
							c->dead_branch_b_to_print  =      deadtree * s->value[TREE_BRANCH_C]  ; 
							c->dead_croot_b_to_print   =      deadtree * s->value[TREE_CROOT_C] ; 

							c->thinned_tree_to_print   =       0    ;     
							c->hwp_to_print            =       0.   ;
							c->thinned_branch_to_print =      0.  ; /*  branch volume removed (m3/ha/yr) */
							c->thinned_stem_to_print   =      0. ;
							c->thinned_stem2_to_print  =      0.;

							EOY_print_output_class_level_mortality(c, height, dbh, age, species, year);
		

							/* remove dead C and N biomass */
							tree_biomass_remove ( c, height, dbh, age, species, deadtree , nat_man);

						
							/* litter fluxes and pools */
							littering             ( c, s );

							// ad this point we should have the printing of the tree biomass removed and everything

							/* call remove_tree_class */
							if ( ! tree_class_remove(c, height, dbh, age, species) )
							{
								logger_error(g_debug_log, "unable to remove tree class");
								exit(1);
							}

		
							deadtree = 0;
							livetree = 0;

							/* check */
							CHECK_CONDITION( s->value[CANOPY_COVER_PROJ] ,  > , s->value[CANOPY_COVER_PROJ] + eps );
							
							c->GREFFMORT_HAPPENS = 1 ;

							goto height_end; 

						}	
				}
			}
			next_dbh :
		}

		} // loop in the same layer for all the heigth
	}


    height_end:

	 							//printf("in mortaliy GEff  height              = %d  \n", height);
								// printf("in mortaliy GEff  m->cells[cell].heights_count               = %d \n ", m->cells[cell].heights_count);
								//printf(" in SELFTHINNING  c->tree_layers_count             = %d \n ", c->tree_layers_count);

                                 //c->GREFFMORT_HAPPENS = 1 ;
								 //printf(" in mortaliy GEff c->GREFFMORT_HAPPENS            = %d \n ", c->GREFFMORT_HAPPENS);


								//printf("in mortaliy GEff  height              = %d  \n", height);


	/* reset values for layer (they are recomputed below) */
	c->tree_layers[layer].layer_n_height_class = 0;
	c->tree_layers[layer].layer_n_trees        = 0;
	c->tree_layers[layer].layer_density        = 0;

	/* REcompute numbers of height classes, tree number and density after mortality within each layer */
	logger(g_debug_log, "REcompute numbers of height classes, tree number and density after mortality within each layer\n\n");


	for ( height = c->heights_count -1; height >= 0 ; --height )
	{
		if( layer == c->heights[height].height_z )
		{
			/* recompute number of height classes */
			c->tree_layers[layer].layer_n_height_class += 1;
		}
		for ( dbh = c->heights[height].dbhs_count - 1; dbh >= 0; --dbh)
		{
			for ( age = 0; age < c->heights[height].dbhs[dbh].ages_count ; ++age )
			{
				for ( species = 0; species < c->heights[height].dbhs[dbh].ages[age].species_count; ++species )
				{
					s = &c->heights[height].dbhs[dbh].ages[age].species[species];

					if( layer == c->heights[height].height_z )
					{
						/* recompute number of trees for each layer */
						c->tree_layers[layer].layer_n_trees += s->counter[N_TREE];
					}
				}
			}
		}
	}

	/*recompute density for each layer */
	c->tree_layers[layer].layer_density = c->tree_layers[layer].layer_n_trees / g_settings->sizeCell;

	logger(g_debug_log, "-layer            = %d \n", layer);
	logger(g_debug_log, "-height class(es) = %d layer\n", c->tree_layers[layer].layer_n_height_class);
	logger(g_debug_log, "-number of trees  = %d layer\n", c->tree_layers[layer].layer_n_trees);
	logger(g_debug_log, "-density          = %f layer\n", c->tree_layers[layer].layer_density);
	logger(g_debug_log, "-Dead tree(s)     = %d trees\n", deadtree);

     //printf("c->daily_dead_tree SELF THINNING %d \n",c->daily_dead_tree );
	/* reset dead tree */
	deadtree = 0;


}

void self_thinning_mortality(cell_t *const c, const int layer, const int year)
{
	int height;
	int dbh;
	int age;
	int species;
	int deadtree;
	int livetree;
	int nat_man;   /* natural or managed mortality 0 = natural; 1 = managed */

	height_t *h;
	dbh_t *d;
	age_t *a;
	species_t *s;

	nat_man = 0;


    //printf("** IN SELF THINNING MORTALITY **\n");

	/* "First, large plants suppress small plants.
	 * The result is a “hierarchy of dominance and suppression” in which the smaller plants
	 * are at an accumulating disadvantage and finally die".
	 * Westoby, 1984, Advances in Ecological Research */

	/* NOTE: it mainly follows rationale of TREEDYN 3 "Crowding Competition mortality function",
	 * Bossel, 1996, Ecological Modelling (eq. 29) */

	/* This is a sort of self-thinning rules not based on biomass but based on canopy cover */

	/* the model makes die trees of the lower height class for that layer because
	it passes through the function sort_by_height_desc the height classes starting from the lowest */

	logger(g_debug_log, "\n\n*****SELF THINNING MORTALITY for layer %d*****\n", layer);

	for ( height = c->heights_count - 1; height >= 0; --height )
	{
		/* first of all sort by descending height */
		/* model makes die lower height-class in each layer first */
#ifndef USE_NEW_OUTPUT
		qsort (c->heights, c->heights_count, sizeof (height_t), sort_by_heights_desc);
#endif

		h = &c->heights[height];

		//	printf(" in self thin c->heights[height] %f!!!\n",c->heights[height].value);

		for ( dbh = h->dbhs_count - 1; dbh >= 0; --dbh)
		{
			d = &c->heights[height].dbhs[dbh];

			for ( age = 0; age < d->ages_count ; ++age )
			{
				a = &c->heights[height].dbhs[dbh].ages[age];

				for ( species = 0; species < a->species_count; ++species )
				{
					s = &c->heights[height].dbhs[dbh].ages[age].species[species];

					logger(g_debug_log, "MORTALITY BASED ON HIGH CANOPY COVER height %f species %s dbh %f !!!\n", h->value, s->name, d->value);

					/* initialize variables */
					livetree = s->counter[N_TREE];
					deadtree = 0;

					//printf(" SELF THINNING 1 species %s!!!\n", s->name);

					/* mortality */
					//FIXME FOR MULTICLASS IN THE SAME LAYER POURPOSES
					//					while ( c->tree_layers[layer].layer_cover_proj >= g_settings->max_layer_cover )
					//					{
					while (s->value[DBHDC_EFF] <  s->value[DBHDCMIN] )    //ddalmo may23: test include if layer_cc_proj > max_layer_cc_proj
					{

                        //printf(" SELF THINNING 2 species %s!!!\n", s->name);

						//printf("SONO IN SELFTHINNING \n");

                      /**FIXME**/

						/* remove one single tree at each run */
						++deadtree;
						--livetree;

						 //printf("deadtree %d!!!\n", deadtree);


						/* check */
						CHECK_CONDITION(livetree, ==, 0);

						/* impose DBHDC_EFF = DBHDCMIN */
						//s->value[DBHDC_EFF] = s->value[DBHDCMIN];

						/* update at cell level */
						++c->daily_dead_tree;

						/* update layer trees */
						--c->tree_layers[layer].layer_n_trees;

						/* update layer density */
						c->tree_layers[layer].layer_density = c->tree_layers[layer].layer_n_trees / g_settings->sizeCell;

						/* update layer cover proj */
						c->tree_layers[layer].layer_cover_proj  -= s->value[CANOPY_COVER_PROJ];

						if ( livetree > 0 )
						{
							/* recompute dbhdc eff */
							dbhdc_function         (c, layer, height, dbh, age, species, year);

							/* recompute crown area */
							crown_allometry        ( c, height, dbh, age, species );

							/* recompute canopy cover with */
							s->value[CANOPY_COVER_PROJ] = s->value[CROWN_AREA_PROJ] * livetree / g_settings->sizeCell;

							/* check for recompued canopy cover */
							c->tree_layers[layer].layer_cover_proj += s->value[CANOPY_COVER_PROJ];
						}
					}

					//printf(" SONO THIN PRIMA DI TREE BIOMASS REMOVE \n");

					//FIXME
					/* remove dead C and N biomass */
					tree_biomass_remove ( c, height, dbh, age, species, deadtree, nat_man );

					/* update live and dead tree (do not move above) */
					s->counter[DEAD_TREE] += deadtree;
					s->counter[N_TREE]     = livetree;
					c->n_trees            -= deadtree;

					deadtree = 0;
					livetree = 0;

					/* check */
					CHECK_CONDITION( s->value[CANOPY_COVER_PROJ] ,  > , s->value[CANOPY_COVER_PROJ] + eps );

					     // ddalmo
			            // printf(" in self_thinning \n");
				}
			}
		}
	
	}

	/* reset values for layer (they are recomputed below) */
	c->tree_layers[layer].layer_n_height_class = 0;
	c->tree_layers[layer].layer_n_trees        = 0;
	c->tree_layers[layer].layer_density        = 0;

	/* REcompute numbers of height classes, tree number and density after mortality within each layer */
	logger(g_debug_log, "REcompute numbers of height classes, tree number and density after mortality within each layer\n\n");


	for ( height = c->heights_count -1; height >= 0 ; --height )
	{
		if( layer == c->heights[height].height_z )
		{
			/* recompute number of height classes */
			c->tree_layers[layer].layer_n_height_class += 1;
		}
		for ( dbh = c->heights[height].dbhs_count - 1; dbh >= 0; --dbh)
		{
			for ( age = 0; age < c->heights[height].dbhs[dbh].ages_count ; ++age )
			{
				for ( species = 0; species < c->heights[height].dbhs[dbh].ages[age].species_count; ++species )
				{
					s = &c->heights[height].dbhs[dbh].ages[age].species[species];

					if( layer == c->heights[height].height_z )
					{
						/* recompute number of trees for each layer */
						c->tree_layers[layer].layer_n_trees += s->counter[N_TREE];
					}
				}
			}
		}
	}

	/*recompute density for each layer */
	c->tree_layers[layer].layer_density = c->tree_layers[layer].layer_n_trees / g_settings->sizeCell;

	logger(g_debug_log, "-layer            = %d \n", layer);
	logger(g_debug_log, "-height class(es) = %d layer\n", c->tree_layers[layer].layer_n_height_class);
	logger(g_debug_log, "-number of trees  = %d layer\n", c->tree_layers[layer].layer_n_trees);
	logger(g_debug_log, "-density          = %f layer\n", c->tree_layers[layer].layer_density);
	logger(g_debug_log, "-Dead tree(s)     = %d trees\n", deadtree);

	/* reset dead tree */
	deadtree = 0;


}

/******************************************************************************************************************/
/******************************************************************************************************************/

int growth_efficiency_mortality ( cell_t *const c, const int height, const int dbh, const int age, const int species , const int year)
{
	/* this function superimpose mortality for all trees in class when reserves
	 * go under zero assuming that reserve pool hasn't be refilled during the day
	 * and make the tree class die */

	int nat_man;   /* natural or managed mortality 0 = natural; 1 = managed */
    int deadtree;

	species_t* s;

	nat_man = 0;
   

	s = &c->heights[height].dbhs[dbh].ages[age].species[species];

	if( s->value[RESERVE_C] < 0 )
	{
		//printf("** IN GROWTH EFFICIENCY MORTALITY **\n");
		printf("growth efficiency mortality species %s!!!\n", s->name);
  
        deadtree = s->counter[N_TREE];
       
		/* reset to zero n_trees */
		s->counter[DEAD_TREE] += deadtree;
		s->counter[N_TREE]    = 0;
		c->n_trees -= deadtree;

        // jan 2025 save mortality data / save info for printing
        // if entire class is remove, call here the output print function
        // save info for printing (auxiliary variables)

		//new var to add 
		c->dead_tree_to_print     =       deadtree ;  
		c->dead_stem_b_to_print   =       deadtree * s->value[TREE_STEM_C]  ; 
		c->dead_branch_b_to_print  =      deadtree * s->value[TREE_BRANCH_C]  ; 
		c->dead_croot_b_to_print   =      deadtree * s->value[TREE_CROOT_C] ; 

		c->thinned_tree_to_print   =       0    ;     
		c->hwp_to_print            =       0.;
		c->thinned_branch_to_print =      0.  ; /*  stem volume removed (m3/ha/yr) */
		c->thinned_stem_to_print   =   0. ;
		c->thinned_stem2_to_print  =    0.;

		EOY_print_output_class_level_mortality(c, height, dbh, age, species, year);


		/* remove dead C and N biomass */
		tree_biomass_remove ( c, height, dbh, age, species,deadtree, nat_man);

		/*** update cell level carbon fluxes (gC/m2/day) ***/
		// c->daily_leaf_carbon        -= (s->value[C_LEAF_TO_LITR]   * 1e6 / g_settings->sizeCell);
		//c->daily_froot_carbon       -= (s->value[C_FROOT_TO_LITR]  * 1e6 / g_settings->sizeCell);
		//c->daily_stem_carbon        -= (s->value[C_STEM_TO_CWD]    * 1e6 / g_settings->sizeCell);
		//c->daily_croot_carbon       -= (s->value[C_CROOT_TO_CWD]   * 1e6 / g_settings->sizeCell);
		//c->daily_branch_carbon      -= (s->value[C_BRANCH_TO_CWD]  * 1e6 / g_settings->sizeCell);
		////c->daily_reserve_carbon     -= (s->value[C_RESERVE_TO_CWD] * 1e6 / g_settings->sizeCell);
		//c->daily_reserve_carbon     -= (s->value[C_RESERVE_TO_LITR] * 1e6 / g_settings->sizeCell);    
		//c->daily_fruit_carbon       -= (s->value[C_FRUIT_TO_CWD]   * 1e6 / g_settings->sizeCell);

        //FIXME : to update the cell level N pools too!

		/* litter fluxes and pools */
		littering             ( c, s );

		/* call remove_tree_class */
		if ( ! tree_class_remove(c, height, dbh, age, species) )
		{
			logger_error(g_debug_log, "unable to remove tree class");
			exit(1);
		}
		return 1;
	}
	return 0;
}

/******************************************************************************************************************/
/******************************************************************************************************************/

// This function is not called/used

int annual_growth_efficiency_mortality ( cell_t *const c, const int height, const int dbh, const int age, const int species )
{
	/* this function superimpose mortality for all trees in class when reserves
	 * go under zero assuming that reserve pool hasn't be refilled during the year
	 * and make the tree class die */

	int nat_man;   /* natural or managed mortality 0 = natural; 1 = managed */
	int deadtree;

	species_t* s;
   
	nat_man = 0;

	s = &c->heights[height].dbhs[dbh].ages[age].species[species];

	if( s->value[RESERVE_C] < 0 )
	{
		deadtree = s->counter[N_TREE]; 
		/* remove dead C and N biomass */
		tree_biomass_remove ( c, height, dbh, age, species, deadtree, nat_man );

		/* call remove_tree_class */
		if ( ! tree_class_remove(c, height, dbh, age, species) )
		{
			logger_error(g_debug_log, "unable to remove tree class");
			exit(1);
		}

		/* reset to zero n_trees */
		s->counter[DEAD_TREE] += deadtree;
		s->counter[N_TREE]     = 0;
		c->n_trees -= deadtree;

		return 1;
	}
	return 0;
}

/******************************************************************************************************************/
/******************************************************************************************************************/

/* Age mortality function from LPJ-GUESS */
void age_mortality ( cell_t *const c, const int height, const int dbh, const int age, const int species, const int year )
{
	int livetree       = 0;
	int deadtree       = 0;
	double mort_factor = 3;     /* age mortality factor 3 as in Smith et al., (LPJ-Guess) */
	double mort_log    = 0.001; /* age mortality log 0.001 as in Smith et al., (LPJ-Guess) */
	double age_mort;
	int nat_man;                /* natural or managed mortality 0 = natural; 1 = managed */

	age_t *a;
	species_t *s;

	nat_man = 0;

	a = &c->heights[height].dbhs[dbh].ages[age];

	//printf("** IN AGE RELATED MORTALITY **\n");

	// ALESSIOR
	//if ( ! a ) return;
	s = &a->species[species];

	/* Age probability function */
	age_mort = ( - ( mort_factor * log ( mort_log ) ) / (s->value[MAXAGE])) * pow (((double)a->value /s->value[MAXAGE]), 2.);

   // printf("ddalmo age_mort  = %f \n",age_mort);

	livetree = s->counter[N_TREE];

	//printf("ddalmo LIVETREES  = %d \n",livetree);

	deadtree = (int)(livetree * age_mort);

	logger(g_debug_log, "dead trees = %d\n", deadtree);

    if ( ( deadtree ) > 1 )  //FIXME set >=1
	{
		logger(g_debug_log, "**MORTALITY based on Tree Age (LPJ)**\n");

		if ( livetree > deadtree)
		{
			/* update C and N biomass */
			tree_biomass_remove ( c, height, dbh, age, species, deadtree, nat_man );

            /* update at class level */
		    s->counter[DEAD_TREE] += deadtree;
		    s->counter[N_TREE]    -= deadtree;

		    /* check */
		    CHECK_CONDITION(s->counter[N_TREE],    <=, 0);
		    CHECK_CONDITION(s->counter[DEAD_TREE], <, 0);
			 
		}
		else
		{
            deadtree =    s->counter[N_TREE];

			// jan 2025 save mortality data / save info for printing
        	// if entire class is remove, call here the output print function
        	// save info for printing (auxiliary variables)

			//new var to add 
			c->dead_tree_to_print     =       deadtree ;  
			c->dead_stem_b_to_print   =       deadtree * s->value[TREE_STEM_C]  ; 
			c->dead_branch_b_to_print  =      deadtree * s->value[TREE_BRANCH_C]  ; 
			c->dead_croot_b_to_print   =      deadtree * s->value[TREE_CROOT_C] ; 

			c->thinned_tree_to_print   =       0    ;     
			c->hwp_to_print            =       0.;
			c->thinned_branch_to_print =      0.  ; /*  stem volume removed (m3/ha/yr) */
			c->thinned_stem_to_print   =   0. ;
			c->thinned_stem2_to_print  =    0.;

			EOY_print_output_class_level_mortality(c, height, dbh, age, species, year);
		

			/* update C and N biomass */
			tree_biomass_remove ( c, height, dbh, age, species, deadtree, nat_man );

			/* litter fluxes and pools */
			littering             ( c, s );

			if ( ! tree_class_remove(c, height, dbh, age, species) )
			{
				logger_error(g_debug_log, "unable to remove tree class");
				exit(1);
			}
		}

		

		/* update at cell level */
		c->daily_dead_tree   += deadtree;
		c->monthly_dead_tree += deadtree;
		c->annual_dead_tree  += deadtree;
		c->n_trees           -= deadtree;
        
		/* check */
		CHECK_CONDITION(c->daily_dead_tree  , <, 0);
		CHECK_CONDITION(c->monthly_dead_tree, <, 0);
		CHECK_CONDITION(c->annual_dead_tree , <, 0);

	}
	//printf("c->daily_dead_tree 2  %d \n",c->daily_dead_tree );
}

/******************************************************************************************************************/
/******************************************************************************************************************/

void stochastic_mortality ( cell_t *const c, const int height, const int dbh, const int age, const int species, const int year )
{
	int livetree         = 0;
	int deadtree         = 0;
	double age_mort_fact = 0.01;   /* ntrees years-1 */

	int nat_man;   /* natural or managed mortality 0 = natural; 1 = managed */

	/* stochastic mortality factor which considers that (see BIOME-BGC) a stochastic mortality happens
	 * BUT for reasons differet (e.g. pests) to age mortality factor */

	age_t *a;
	species_t *s;

	nat_man = 0;

	a = &c->heights[height].dbhs[dbh].ages[age];
	// ALESSIOR
	//if ( ! a ) return;
	s = &a->species[species];

	livetree = s->counter[N_TREE];

	/* stochastic probability function */
	deadtree = (int)( livetree * age_mort_fact );
	//printf("livetree %d \n",livetree );

	//printf("s->counter[DEAD_TREE] %d \n",s->counter[DEAD_TREE]);

	logger(g_debug_log, "dead trees = %d\n", deadtree);

	if ( ( deadtree ) >= 1 )  // //FIXME set >1
	{
		logger(g_debug_log, "**MORTALITY based on stochasticity **\n");

		if ( livetree > deadtree)
		{

			/* remove C and N biomass */
			tree_biomass_remove ( c, height, dbh, age, species, deadtree, nat_man );

			
			/* update at class level */
			s->counter[DEAD_TREE] += deadtree;
			s->counter[N_TREE]    -= deadtree;

			//printf("s->counter[N_TREE]  %d \n",s->counter[N_TREE] );
			//printf("c->daily_dead_tree 1 %d \n",c->daily_dead_tree );
			/* check */
			CHECK_CONDITION(s->counter[N_TREE],    <=, 0);
			CHECK_CONDITION(s->counter[DEAD_TREE], <, 0);

		}
		else
		{
			
			deadtree =   s->counter[N_TREE];

			// jan 2025 save mortality data / save info for printing
       		// if entire class is remove, call here the output print function
        	// save info for printing (auxiliary variables)

			//new var to add 
			c->dead_tree_to_print     =       deadtree ;  
			c->dead_stem_b_to_print   =       deadtree * s->value[TREE_STEM_C]  ; 
			c->dead_branch_b_to_print  =      deadtree * s->value[TREE_BRANCH_C]  ; 
			c->dead_croot_b_to_print   =      deadtree * s->value[TREE_CROOT_C] ; 

			c->thinned_tree_to_print   =       0    ;     
			c->hwp_to_print            =       0.;
			c->thinned_branch_to_print =      0.  ; /*  stem volume removed (m3/ha/yr) */
			c->thinned_stem_to_print   =   0. ;
			c->thinned_stem2_to_print  =    0.;

			EOY_print_output_class_level_mortality(c, height, dbh, age, species, year);
		

			/* update C and N biomass */
			tree_biomass_remove ( c, height, dbh, age, species, deadtree, nat_man );

		    /* litter fluxes and pools */
		    littering             ( c, s );

			if ( ! tree_class_remove(c, height, dbh, age, species) )
			{
				logger_error(g_debug_log, "unable to remove tree class");
				exit(1);
			}
		}

		/* update at cell level */
		c->daily_dead_tree   += deadtree;
		c->monthly_dead_tree += deadtree;
		c->annual_dead_tree  += deadtree;
		c->n_trees           -= deadtree;

		/* check */
		CHECK_CONDITION(c->daily_dead_tree  , <, 0);
		CHECK_CONDITION(c->monthly_dead_tree, <, 0);
		CHECK_CONDITION(c->annual_dead_tree , <, 0);

	}
	//printf("c->daily_dead_tree 2  %d \n",c->daily_dead_tree );
}

/******************************************************************************************************************/
/******************************************************************************************************************/

void self_pruning( cell_t *const c, const int height, const int dbh, const int age, const int species, const double old, const double current )
{
	double self_pruning_ratio;

	//ddalmo october 2024 

	species_t *s;
	s = &c->heights[height].dbhs[dbh].ages[age].species[species];

	logger(g_debug_log, "\n\n*****SELF PRUNING*****\n");

	//printf(" SELF PRUNING  species %s!!!\n", s->name);


	//printf(" ** IN SELF_PRUNING ** \n");

	/* reduce proportionally branch biomass to the crown area reduction */

	/* compute FRACTION in crown area reduction for self-pruning */
	self_pruning_ratio = current / old;

	/* check for precision control */
	if ( self_pruning_ratio > 1. ) self_pruning_ratio = 1.;

	if (s->counter[N_TREE] < 1 )   // only for testing
	{

	/*** branch self-pruning ***/
	/** carbon **/

	/* remove biomass from branch */
	// values do not have to sum up but computed ex novo
   // printf(" BEFORE PRUNING  s->value[BRANCH_C]  %f!!!\n",s->value[BRANCH_C] );

   	s->value[BRANCH_C]    -= ( s->value[BRANCH_C] * ( 1. - self_pruning_ratio ) );

	// printf(" AFTER PRUNING  s->value[BRANCH_C]  %f!!!\n",s->value[BRANCH_C] );

     s->value[BRANCH_SAPWOOD_C] -= ( s->value[BRANCH_SAPWOOD_C]     * ( 1. - self_pruning_ratio ) );
   	s->value[BRANCH_LIVEWOOD_C] -=( s->value[BRANCH_LIVEWOOD_C]     * ( 1. - self_pruning_ratio ) );
    s->value[BRANCH_HEARTWOOD_C] -=( s->value[BRANCH_HEARTWOOD_C]     * ( 1. - self_pruning_ratio ) );
	 s->value[BRANCH_DEADWOOD_C] -=( s->value[BRANCH_DEADWOOD_C]     * ( 1. - self_pruning_ratio ) );

     s->value[TOT_SAPWOOD_C]  -= ( s->value[BRANCH_SAPWOOD_C]     * ( 1. - self_pruning_ratio ) );
	   s->value[TOT_LIVEWOOD_C] -=( s->value[BRANCH_LIVEWOOD_C]     * ( 1. - self_pruning_ratio ) );
    s->value[TOT_HEARTWOOD_C] -=( s->value[BRANCH_HEARTWOOD_C]     * ( 1. - self_pruning_ratio ) );
	 s->value[TOT_DEADWOOD_C] -=( s->value[BRANCH_DEADWOOD_C]     * ( 1. - self_pruning_ratio ) );


    

	//s->value[BRANCH_C_TO_REMOVE]             += ( s->value[BRANCH_C] * ( 1. - self_pruning_ratio ) );

	/* adding to BRANCH C transfer pools */
	s->value[C_BRANCH_TO_RESERVE]            = ( s->value[BRANCH_C] * ( 1. - self_pruning_ratio ) )  * C_FRAC_TO_RETRANSL;
	s->value[C_BRANCH_TO_CWD]                += ( s->value[BRANCH_C] * ( 1. - self_pruning_ratio ) )   * ( 1. - C_FRAC_TO_RETRANSL );

	/* adding BRANCH SAPWOOD and HEARTWOOD to CWD */
	// TODO to keep it simple we assume the retranslocation occurrs from heartwood too
	// the retranslocation should strictly occurr from the livewood pool alone

	s->value[C_BRANCH_SAPWOOD_TO_CWD]       += ( s->value[BRANCH_SAPWOOD_C]   * ( 1. - self_pruning_ratio ) ) * ( 1. - C_FRAC_TO_RETRANSL );
	s->value[C_BRANCH_HEARTWOOD_TO_CWD]     += ( s->value[BRANCH_HEARTWOOD_C] * ( 1. - self_pruning_ratio ) ) * ( 1. - C_FRAC_TO_RETRANSL );

	s->value[C_BRANCH_LIVEWOOD_TO_CWD]      += ( s->value[BRANCH_LIVEWOOD_C]   * ( 1. - self_pruning_ratio ) ) * ( 1. - C_FRAC_TO_RETRANSL );
	s->value[C_BRANCH_DEADWOOD_TO_CWD]      += ( s->value[BRANCH_DEADWOOD_C] * ( 1. - self_pruning_ratio ) ) * ( 1. - C_FRAC_TO_RETRANSL );

	/* adding to reserve and CWD pool */
	//s->value[C_TO_RESERVE]                   += s->value[C_BRANCH_TO_RESERVE];
	s->value[C_TO_CWD]                       +=  ( s->value[BRANCH_C] * ( 1. - self_pruning_ratio ) )   * ( 1. - C_FRAC_TO_RETRANSL );
    s->value[RESERVE_C]   += s->value[C_BRANCH_TO_RESERVE];
	/** nitrogen **/
	// TODO include Live-dead-hart-sap wood pool too

	/* remove biomass from branch */
	s->value[BRANCH_N_TO_REMOVE]             = ( s->value[BRANCH_N]           * ( 1. - self_pruning_ratio ) );

	/* adding to BRANCH_N C transfer pools */
	s->value[N_BRANCH_TO_RESERVE]            = ( s->value[BRANCH_N]           * ( 1. - self_pruning_ratio ) )  * N_FRAC_TO_RETRANSL;
	s->value[N_BRANCH_TO_CWD]                += s->value[BRANCH_N_TO_REMOVE]   * ( 1. - N_FRAC_TO_RETRANSL );

	s->value[N_TO_RESERVE]                   += s->value[N_BRANCH_TO_RESERVE];
	s->value[N_TO_CWD]                       += s->value[N_BRANCH_TO_CWD];

	/*** coarse root self-pruning ***/
	/** carbon **/

	/* remove biomass from coarse root */
	// s->value[CROOT_C_TO_REMOVE]             += ( s->value[CROOT_C]             * ( 1. - self_pruning_ratio ) );
   s->value[CROOT_C]     -= ( s->value[CROOT_C]             * ( 1. - self_pruning_ratio ) );
   s->value[CROOT_SAPWOOD_C] -= ( s->value[CROOT_SAPWOOD_C]     * ( 1. - self_pruning_ratio ) );
   s->value[CROOT_LIVEWOOD_C] -=( s->value[CROOT_LIVEWOOD_C]     * ( 1. - self_pruning_ratio ) );
    s->value[CROOT_HEARTWOOD_C] -=( s->value[CROOT_HEARTWOOD_C]     * ( 1. - self_pruning_ratio ) );
	 s->value[CROOT_DEADWOOD_C] -=( s->value[CROOT_DEADWOOD_C]     * ( 1. - self_pruning_ratio ) );

    s->value[TOT_SAPWOOD_C] -= ( s->value[CROOT_SAPWOOD_C]     * ( 1. - self_pruning_ratio ) );
   s->value[TOT_LIVEWOOD_C] -=( s->value[CROOT_LIVEWOOD_C]     * ( 1. - self_pruning_ratio ) );
    s->value[TOT_HEARTWOOD_C] -=( s->value[CROOT_HEARTWOOD_C]     * ( 1. - self_pruning_ratio ) );
	 s->value[TOT_DEADWOOD_C] -=( s->value[CROOT_DEADWOOD_C]     * ( 1. - self_pruning_ratio ) );


	/* adding to CROOT_C transfer pools */
	s->value[C_CROOT_TO_RESERVE]            = ( s->value[CROOT_C] * ( 1. - self_pruning_ratio ) )     * C_FRAC_TO_RETRANSL;
	s->value[C_CROOT_TO_CWD]                += ( s->value[CROOT_C]* ( 1. - self_pruning_ratio ) )   * ( 1. - C_FRAC_TO_RETRANSL );

	/* adding CROOT_C SAPWOOD and HEARTWOOD to CWD */
	s->value[C_CROOT_SAPWOOD_TO_CWD]        += ( s->value[CROOT_SAPWOOD_C]     * ( 1. - self_pruning_ratio ) )* ( 1. - C_FRAC_TO_RETRANSL );
	s->value[C_CROOT_HEARTWOOD_TO_CWD]      += ( s->value[CROOT_HEARTWOOD_C]   * ( 1. - self_pruning_ratio ) )* ( 1. - C_FRAC_TO_RETRANSL );

	s->value[C_CROOT_LIVEWOOD_TO_CWD]      += ( s->value[CROOT_LIVEWOOD_C]   * ( 1. - self_pruning_ratio ) ) * ( 1. - C_FRAC_TO_RETRANSL );
	s->value[C_CROOT_DEADWOOD_TO_CWD]      += ( s->value[CROOT_DEADWOOD_C] * ( 1. - self_pruning_ratio ) ) * ( 1. - C_FRAC_TO_RETRANSL );

	/* adding to reserve and CWD pool */
	//s->value[C_TO_RESERVE]                  += s->value[C_CROOT_TO_RESERVE];
	s->value[C_TO_CWD]                      += ( s->value[CROOT_C]* ( 1. - self_pruning_ratio ) )   * ( 1. - C_FRAC_TO_RETRANSL );

    	s->value[RESERVE_C]   +=  s->value[C_CROOT_TO_RESERVE];

	/** nitrogen **/
        // TODO include Live-dead-hart-sap wood pool too

	/* remove biomass from coarse root */
	s->value[CROOT_N_TO_REMOVE]             = ( s->value[CROOT_N]             * ( 1. - self_pruning_ratio ) );

	/* adding to CROOT transfer pools */
	s->value[N_CROOT_TO_RESERVE]            = s->value[CROOT_N_TO_REMOVE]     * N_FRAC_TO_RETRANSL;
	s->value[N_CROOT_TO_CWD]                += s->value[CROOT_N_TO_REMOVE]     * ( 1. - N_FRAC_TO_RETRANSL );

	s->value[N_TO_RESERVE]                  += s->value[N_CROOT_TO_RESERVE];
	s->value[N_TO_CWD]                      += s->value[N_CROOT_TO_CWD];


	// update tree level biomass 


	/*** update class level carbon mass pools ***/
	 /*** update cell level carbon fluxes (gC/m2/day)***/
      	
	//c->daily_leaf_carbon        += (s->value[C_TO_LEAF]    * 1e6 / g_settings->sizeCell);
	//c->daily_froot_carbon       += (s->value[C_TO_FROOT]   * 1e6 / g_settings->sizeCell);
	//c->daily_stem_carbon        += (s->value[C_TO_STEM]    * 1e6 / g_settings->sizeCell);
	c->daily_croot_carbon       += (s->value[C_TO_CROOT]   * 1e6 / g_settings->sizeCell);
	c->daily_branch_carbon      += (s->value[C_TO_BRANCH]  * 1e6 / g_settings->sizeCell);
	c->daily_reserve_carbon     += (s->value[C_TO_RESERVE] * 1e6 / g_settings->sizeCell);
	//c->daily_fruit_carbon       += (s->value[C_TO_FRUIT]   * 1e6 / g_settings->sizeCell);

	/*** update cell level carbon pools (gC/m2)***/
	//c->leaf_carbon              += (s->value[C_TO_LEAF]    * 1e6 / g_settings->sizeCell);
	//c->froot_carbon             += (s->value[C_TO_FROOT]   * 1e6 / g_settings->sizeCell);
	//c->stem_carbon              += (s->value[C_TO_STEM]    * 1e6 / g_settings->sizeCell);
	 //printf(" BEFORE PRUNING  c->branch_carbon  %f!!!\n",c->branch_carbon );

	c->branch_carbon            -= (( s->value[BRANCH_C] * ( 1. - self_pruning_ratio ) )  * 1e6 / g_settings->sizeCell);
	//printf(" AFTER PRUNING  c->branch_carbon  %f!!!\n",c->branch_carbon );
	//printf(" BEFORE PRUNING  c->croot_carbon  %f!!!\n",c->croot_carbon );

	c->croot_carbon             -= (( s->value[CROOT_C] * ( 1. - self_pruning_ratio ) )  * 1e6 / g_settings->sizeCell);
	//printf(" AFTER PRUNING  c->croot_carbon  %f!!!\n",c->croot_carbon );
	//printf(" BEFORE PRUNING  c->reserve_carbon  %f!!!\n",c->reserve_carbon );
     //printf(" BEFORE PRUNING  flyx to reserve_carbon  %f!!!\n",((( s->value[C_CROOT_TO_RESERVE]+ s->value[C_BRANCH_TO_RESERVE])) * 1e6 / g_settings->sizeCell) );
	c->reserve_carbon           += ((( s->value[C_CROOT_TO_RESERVE]+ s->value[C_BRANCH_TO_RESERVE])) * 1e6 / g_settings->sizeCell);
	//printf(" AFTER PRUNING  c->reserve_carbon  %f!!!\n",c->reserve_carbon );
	//c->fruit_carbon             += (s->value[C_TO_FRUIT]   * 1e6 / g_settings->sizeCell);



	
					/* compute tree average C pools */
					//s->value[TREE_LEAF_C]                = (s->value[LEAF_C]             / (double)s->counter[N_TREE]);
					//s->value[TREE_STEM_C]                = (s->value[STEM_C]             / (double)s->counter[N_TREE]);
					//s->value[TREE_FROOT_C]               = (s->value[FROOT_C]            / (double)s->counter[N_TREE]);
					s->value[TREE_CROOT_C]               = (s->value[CROOT_C]            / (double)s->counter[N_TREE]);
					s->value[TREE_RESERVE_C]             = (s->value[RESERVE_C]          / (double)s->counter[N_TREE]);
					s->value[TREE_BRANCH_C]              = (s->value[BRANCH_C]           / (double)s->counter[N_TREE]);
					//s->value[TREE_FRUIT_C]               = (s->value[FRUIT_C]            / (double)s->counter[N_TREE]);

				//	s->value[TREE_STEM_SAPWOOD_C]        = (s->value[STEM_SAPWOOD_C]     / (double)s->counter[N_TREE]);
				//	s->value[TREE_STEM_HEARTWOOD_C]      = (s->value[STEM_HEARTWOOD_C]   / (double)s->counter[N_TREE]);
					s->value[TREE_CROOT_SAPWOOD_C]       = (s->value[CROOT_SAPWOOD_C]    / (double)s->counter[N_TREE]);
					s->value[TREE_CROOT_HEARTWOOD_C]     = (s->value[CROOT_HEARTWOOD_C]  / (double)s->counter[N_TREE]);
					s->value[TREE_BRANCH_SAPWOOD_C]      = (s->value[BRANCH_SAPWOOD_C]   / (double)s->counter[N_TREE]);
					s->value[TREE_BRANCH_HEARTWOOD_C]    = (s->value[BRANCH_HEARTWOOD_C] / (double)s->counter[N_TREE]);
					s->value[TREE_SAPWOOD_C]             = (s->value[TOT_SAPWOOD_C]      / (double)s->counter[N_TREE]);
					s->value[TREE_HEARTWOOD_C]           = (s->value[TOT_HEARTWOOD_C]    / (double)s->counter[N_TREE]);

				//	s->value[TREE_STEM_LIVEWOOD_C]       = (s->value[STEM_LIVEWOOD_C]    / (double)s->counter[N_TREE]);
				//	s->value[TREE_STEM_DEADWOOD_C]       = (s->value[STEM_DEADWOOD_C]    / (double)s->counter[N_TREE]);
					s->value[TREE_CROOT_LIVEWOOD_C]      = (s->value[CROOT_LIVEWOOD_C]   / (double)s->counter[N_TREE]);
					s->value[TREE_CROOT_DEADWOOD_C]      = (s->value[CROOT_DEADWOOD_C]   / (double)s->counter[N_TREE]);
					s->value[TREE_BRANCH_LIVEWOOD_C]     = (s->value[BRANCH_LIVEWOOD_C]  / (double)s->counter[N_TREE]);
					s->value[TREE_BRANCH_DEADWOOD_C]     = (s->value[BRANCH_DEADWOOD_C]  / (double)s->counter[N_TREE]);
					s->value[TREE_TOT_LIVEWOOD_C]        = (s->value[TOT_LIVEWOOD_C]     / (double)s->counter[N_TREE]);
					s->value[TREE_TOT_DEADWOOD_C]        = (s->value[TOT_DEADWOOD_C]     / (double)s->counter[N_TREE]);

	} // only for testing
}

void self_pruning_old ( cell_t *const c, const int height, const int dbh, const int age, const int species, const double old, const double current )
{
	double self_pruning_ratio;

	species_t *s;
	s = &c->heights[height].dbhs[dbh].ages[age].species[species];

	logger(g_debug_log, "\n\n*****SELF PRUNING*****\n");

	//printf(" SELF PRUNING  species %s!!!\n", s->name);


	//printf(" ** IN SELF_PRUNING ** \n");

	/* reduce proportionally branch biomass to the crown area reduction */

	/* compute FRACTION in crown area reduction for self-pruning */
	self_pruning_ratio = current / old;

	/* check for precision control */
	if ( self_pruning_ratio > 1. ) self_pruning_ratio = 1.;

	/*** branch self-pruning ***/
	/** carbon **/

	/* remove biomass from branch */
	// values do not have to sum up but computed ex novo

	s->value[BRANCH_C_TO_REMOVE]             += ( s->value[BRANCH_C] * ( 1. - self_pruning_ratio ) );

	/* adding to BRANCH C transfer pools */
	s->value[C_BRANCH_TO_RESERVE]            = ( s->value[BRANCH_C] * ( 1. - self_pruning_ratio ) )  * C_FRAC_TO_RETRANSL;
	s->value[C_BRANCH_TO_CWD]                += ( s->value[BRANCH_C] * ( 1. - self_pruning_ratio ) )   * ( 1. - C_FRAC_TO_RETRANSL );

	/* adding BRANCH SAPWOOD and HEARTWOOD to CWD */
	// TODO to keep it simple we assume the retranslocation occurrs from heartwood too
	// the retranslocation should strictly occurr from the livewood pool alone

	s->value[C_BRANCH_SAPWOOD_TO_CWD]       += ( s->value[BRANCH_SAPWOOD_C]   * ( 1. - self_pruning_ratio ) ) * ( 1. - C_FRAC_TO_RETRANSL );
	s->value[C_BRANCH_HEARTWOOD_TO_CWD]     += ( s->value[BRANCH_HEARTWOOD_C] * ( 1. - self_pruning_ratio ) ) * ( 1. - C_FRAC_TO_RETRANSL );

	s->value[C_BRANCH_LIVEWOOD_TO_CWD]      += ( s->value[BRANCH_LIVEWOOD_C]   * ( 1. - self_pruning_ratio ) ) * ( 1. - C_FRAC_TO_RETRANSL );
	s->value[C_BRANCH_DEADWOOD_TO_CWD]      += ( s->value[BRANCH_DEADWOOD_C] * ( 1. - self_pruning_ratio ) ) * ( 1. - C_FRAC_TO_RETRANSL );

	/* adding to reserve and CWD pool */
	s->value[C_TO_RESERVE]                   += s->value[C_BRANCH_TO_RESERVE];
	s->value[C_TO_CWD]                       +=  ( s->value[BRANCH_C] * ( 1. - self_pruning_ratio ) )   * ( 1. - C_FRAC_TO_RETRANSL );

	/** nitrogen **/
	// TODO include Live-dead-hart-sap wood pool too

	/* remove biomass from branch */
	s->value[BRANCH_N_TO_REMOVE]             = ( s->value[BRANCH_N]           * ( 1. - self_pruning_ratio ) );

	/* adding to BRANCH_N C transfer pools */
	s->value[N_BRANCH_TO_RESERVE]            = ( s->value[BRANCH_N]           * ( 1. - self_pruning_ratio ) )  * N_FRAC_TO_RETRANSL;
	s->value[N_BRANCH_TO_CWD]                += s->value[BRANCH_N_TO_REMOVE]   * ( 1. - N_FRAC_TO_RETRANSL );

	s->value[N_TO_RESERVE]                   += s->value[N_BRANCH_TO_RESERVE];
	s->value[N_TO_CWD]                       += s->value[N_BRANCH_TO_CWD];

	/*** coarse root self-pruning ***/
	/** carbon **/

	/* remove biomass from coarse root */
	s->value[CROOT_C_TO_REMOVE]             += ( s->value[CROOT_C]             * ( 1. - self_pruning_ratio ) );

	/* adding to CROOT_C transfer pools */
	s->value[C_CROOT_TO_RESERVE]            = ( s->value[CROOT_C] * ( 1. - self_pruning_ratio ) )     * C_FRAC_TO_RETRANSL;
	s->value[C_CROOT_TO_CWD]                += ( s->value[CROOT_C]* ( 1. - self_pruning_ratio ) )   * ( 1. - C_FRAC_TO_RETRANSL );

	/* adding CROOT_C SAPWOOD and HEARTWOOD to CWD */
	s->value[C_CROOT_SAPWOOD_TO_CWD]        += ( s->value[CROOT_SAPWOOD_C]     * ( 1. - self_pruning_ratio ) )* ( 1. - C_FRAC_TO_RETRANSL );
	s->value[C_CROOT_HEARTWOOD_TO_CWD]      += ( s->value[CROOT_HEARTWOOD_C]   * ( 1. - self_pruning_ratio ) )* ( 1. - C_FRAC_TO_RETRANSL );

	s->value[C_CROOT_LIVEWOOD_TO_CWD]      += ( s->value[CROOT_LIVEWOOD_C]   * ( 1. - self_pruning_ratio ) ) * ( 1. - C_FRAC_TO_RETRANSL );
	s->value[C_CROOT_DEADWOOD_TO_CWD]      += ( s->value[CROOT_DEADWOOD_C] * ( 1. - self_pruning_ratio ) ) * ( 1. - C_FRAC_TO_RETRANSL );

	/* adding to reserve and CWD pool */
	s->value[C_TO_RESERVE]                  += s->value[C_CROOT_TO_RESERVE];
	s->value[C_TO_CWD]                      += ( s->value[CROOT_C]* ( 1. - self_pruning_ratio ) )   * ( 1. - C_FRAC_TO_RETRANSL );

	/** nitrogen **/
        // TODO include Live-dead-hart-sap wood pool too

	/* remove biomass from coarse root */
	s->value[CROOT_N_TO_REMOVE]             = ( s->value[CROOT_N]             * ( 1. - self_pruning_ratio ) );

	/* adding to CROOT transfer pools */
	s->value[N_CROOT_TO_RESERVE]            = s->value[CROOT_N_TO_REMOVE]     * N_FRAC_TO_RETRANSL;
	s->value[N_CROOT_TO_CWD]                += s->value[CROOT_N_TO_REMOVE]     * ( 1. - N_FRAC_TO_RETRANSL );

	s->value[N_TO_RESERVE]                  += s->value[N_CROOT_TO_RESERVE];
	s->value[N_TO_CWD]                      += s->value[N_CROOT_TO_CWD];


}

/******************************************************************************************************************/
/******************************************************************************************************************/


/* Self-thinnig mortality function from 3PG */
/*-----------------------------------------------------------------------------*/
/*                                                                             */
/*                     Self-thinnig mortality (3-PG)                           */
/*This function determines the number of stems to remove to ensure the         */
/*self-thinning rule is satisfied. It applies the Newton-Rhapson method        */
/*to solve for N to an accuracy of 1 stem or less. To change this,             */
/*change the value of "accuracy".                                              */
/*This was the old mortality function:                                         */
/*getMortality = oldN - 1000 * (wSx1000 * oldN / oldW / 1000) ^ (1 / thinPower)*/
/*which has been superceded by the following ..                                */
/*                                                                             */
/*-----------------------------------------------------------------------------*/
/* NOT USED */
//void self_thinning_mortality_3PG (species_t *const s, int years)
//{
//
//	int i;
//	double fN, dfN;
//	double dN, n, x1, x2;
//
//
//	s->value[WS_MAX] = s->value[WSX1000] * pow ( ( 1000 / (double) s->counter[N_TREE] ), s->value[THINPOWER] );
//
//	//ALESSIOC ??????????????
//	//modifified version for 1Km^2 spatial resolution
//	s->value[WS_MAX] = s->value[WSX1000];
//
//	if ( s->value[AV_STEM_MASS_KgDM] > s->value[WS_MAX])
//	{
//		logger(g_log, "*Self-thinnig mortality (3-PG)*\n");
//		logger(g_log, "Average Stem Mass > WSMax\n");
//		logger(g_log, "WS MAX = %f kgC/tree\n", s->value[WS_MAX]);
//		logger(g_log, "Average Stem Mass = %f kgC stem/tree\n", s->value[AV_STEM_MASS_KgDM]);
//		logger(g_log, "Tree Number before Mortality Function = %d\n", s->counter[N_TREE]);
//
//		n = (double)s->counter[N_TREE] / 1000;
//		logger(g_log, "n = %f\n", n);
//
//		x1 = 1000 * s->value[MS] * s->value[STEM_C] / (double)s->counter[N_TREE];
//		logger(g_log, "x1 = %f\n", x1);
//		i = 0;
//		while ( 1 )
//		{
//			i = i + 1;
//			logger(g_log, "i = %d\n", i);
//			x2 = s->value[WSX1000] * pow (n, (1 - s->value[THINPOWER]));
//			logger(g_log, "X2 = %f\n", x2);
//			fN = x2 - x1 * n - (1 - s->value[MS]) * s->value[STEM_C];
//			logger(g_log, "fN = %f\n", fN);
//			dfN = (1 - s->value[THINPOWER]) * x2 / n - x1;
//			logger(g_log, "dfN = %f\n", dfN);
//			dN = -fN / dfN;
//			logger(g_log, "dN = %f\n", dN);
//			n = n + dN;
//			logger(g_log, "n = %f\n", n);
//			if ((fabs(dN) <= eps) || (i >= 5))
//				break;
//		}
//
//		s->counter[DEAD_STEMS] = (int)(s->counter[N_TREE] - 1000 * n);
//		logger(g_log, "Dead Tree In Mortality Function = %d trees \n", s->counter[DEAD_STEMS]);
//
//		//SERGIO CONTROL: if DEAD_STEMS < 0 set it to its minimum plausible value; that is 0
//		if (s->counter[DEAD_STEMS] < 0)
//		{
//			s->counter[DEAD_STEMS]	 = 0;
//		}
//		//control
//		if (s->counter[DEAD_STEMS] > s->counter[N_TREE])
//		{
//			logger(g_log, "ERROR Number of Dead Trees > N Trees\n");
//			logger(g_log, "Dead Trees = %d\n", s->counter[DEAD_STEMS]);
//			logger(g_log, "Live Trees = %d\n", s->counter[N_TREE]);
//		}
//		else
//		{
//			s->counter[N_TREE] = s->counter[N_TREE] - s->counter[DEAD_STEMS];
//			logger(g_log, "Number of Trees  after mortality = %d trees\n", s->counter[N_TREE]);
//			s->value[BIOMASS_FOLIAGE_tDM] = s->value[BIOMASS_FOLIAGE_tDM] - s->value[MF] * s->counter[DEAD_STEMS] * (s->value[BIOMASS_FOLIAGE_tDM] / s->counter[N_TREE]);
//			s->value[BIOMASS_ROOTS_TOT_tDM] = s->value[BIOMASS_ROOTS_TOT_tDM] - s->value[MR] * s->counter[DEAD_STEMS] * (s->value[BIOMASS_ROOTS_TOT_tDM] / s->counter[N_TREE]);
//			s->value[BIOMASS_STEM_tDM] = s->value[BIOMASS_STEM_tDM] - s->value[MS] * s->counter[DEAD_STEMS] * (s->value[BIOMASS_STEM_tDM] / s->counter[N_TREE]);
//			logger(g_log, "Wf after dead = %f tDM/ha\n", s->value[BIOMASS_FOLIAGE_tDM]);
//			logger(g_log, "Wr after dead = %f tDM/ha\n", s->value[BIOMASS_ROOTS_TOT_tDM]);
//			logger(g_log, "Ws after dead = %f tDM/ha\n", s->value[BIOMASS_STEM_tDM]);
//		}
//
//
//
//		//----------------Number of trees after mortality---------------------
//
//
//
//
//		//--------------------------------------------------------------------
//
//		//deselected algorithm for 1Km^2 spatial resolution
//		//s->value[WS_MAX] = s->value[WSX1000] * pow((1000 / (double)s->counter[N_TREE]), s->value[THINPOWER]);
//
//		//modifified version for 1Km^2 spatial resolution
//		s->value[WS_MAX] = s->value[WSX1000];
//
//		logger(g_log, "wSmax = %f KgDM/tree\n", s->value[WS_MAX]);
//		s->value[AV_STEM_MASS_KgDM] = s->value[BIOMASS_STEM_tDM] * 1000.0 / (double)s->counter[N_TREE];
//		logger(g_log, "AvStemMass after dead = %f Kg/tree\n", s->value[AV_STEM_MASS_KgDM]);
//	}
//	else
//	{
//		logger(g_log, "NO MORTALITY based SELF-THINNING RULE\n");
//		logger(g_log, "Average Stem Mass < WSMax\n");
//	}
//
//	logger(g_log, "**********************************\n");
//
//
//
//}

void pruning_daily(matrix_t* m, int cell_index, pruning_t* p)
{
	double self_pruning_ratio;
	cell_t* c;
	species_t *s;

	assert(m);
	assert(p);

	c = &m->cells[cell_index];
	CHECK_CONDITION(c->heights_count, !=, 1);
	CHECK_CONDITION(c->heights[0].dbhs_count, !=, 1);
	CHECK_CONDITION(c->heights[0].dbhs[0].ages_count, !=, 1);
	CHECK_CONDITION(c->heights[0].dbhs[0].ages[0].species_count, !=, 1);

	s = &c->heights[0].dbhs[0].ages[0].species[0];

	logger(g_debug_log, "\n\n*****PRUNING DAILY (%04d-%02d-%02d)*****\n"
								, p->year
								, p->month
								, p->day
	);

	/* reduce proportionally branch biomass to the crown area reduction */

	/* compute FRACTION in crown area reduction for self-pruning */
	self_pruning_ratio = p->intensity / 100;

	/* check for precision control */
	if ( self_pruning_ratio > 1. ) self_pruning_ratio = 1.;

	/*** branch self-pruning ***/
	/** carbon **/

	/* remove biomass from branch */
	s->value[BRANCH_C_TO_REMOVE]             += ( s->value[BRANCH_C]           * ( 1. - self_pruning_ratio ) );

	/* adding to BRANCH C transfer pools */
	s->value[C_BRANCH_TO_RESERVE]            += s->value[BRANCH_C_TO_REMOVE]   * C_FRAC_TO_RETRANSL;
	s->value[C_BRANCH_TO_CWD]                += s->value[BRANCH_C_TO_REMOVE]   * ( 1. - C_FRAC_TO_RETRANSL );

	/* adding BRANCH SAPWOOD and HEARTWOOD to CWD */
	s->value[C_BRANCH_SAPWOOD_TO_CWD]       += ( s->value[BRANCH_SAPWOOD_C]   * ( 1. - self_pruning_ratio ) );
	s->value[C_BRANCH_HEARTWOOD_TO_CWD]     += ( s->value[BRANCH_HEARTWOOD_C] * ( 1. - self_pruning_ratio ) );

	/* adding to reserve and CWD pool */
	s->value[C_TO_RESERVE]                   += s->value[C_BRANCH_TO_RESERVE];
	s->value[C_TO_CWD]                       += s->value[C_BRANCH_TO_CWD];

	/** nitrogen **/

	/* remove biomass from branch */
	s->value[BRANCH_N_TO_REMOVE]             += ( s->value[BRANCH_N]           * ( 1. - self_pruning_ratio ) );

	/* adding to BRANCH_N C transfer pools */
	s->value[N_BRANCH_TO_RESERVE]            += s->value[BRANCH_N_TO_REMOVE]   * N_FRAC_TO_RETRANSL;
	s->value[N_BRANCH_TO_CWD]                += s->value[BRANCH_N_TO_REMOVE]   * ( 1. - N_FRAC_TO_RETRANSL );

	s->value[N_TO_RESERVE]                   += s->value[N_BRANCH_TO_RESERVE];
	s->value[N_TO_CWD]                       += s->value[N_BRANCH_TO_CWD];

	/*** coarse root self-pruning ***/
	/** carbon **/

	/* remove biomass from coarse root */
	s->value[CROOT_C_TO_REMOVE]             += ( s->value[CROOT_C]             * ( 1. - self_pruning_ratio ) );

	/* adding to CROOT_C transfer pools */
	s->value[C_CROOT_TO_RESERVE]            += s->value[CROOT_C_TO_REMOVE]     * C_FRAC_TO_RETRANSL;
	s->value[C_CROOT_TO_CWD]                += s->value[CROOT_C_TO_REMOVE]     * ( 1. - C_FRAC_TO_RETRANSL );;

	/* adding CROOT_C SAPWOOD and HEARTWOOD to CWD */
	s->value[C_CROOT_SAPWOOD_TO_CWD]        += ( s->value[CROOT_SAPWOOD_C]     * ( 1. - self_pruning_ratio ) );
	s->value[C_CROOT_HEARTWOOD_TO_CWD]      += ( s->value[CROOT_HEARTWOOD_C]   * ( 1. - self_pruning_ratio ) );

	/* adding to reserve and CWD pool */
	s->value[C_TO_RESERVE]                  += s->value[C_CROOT_TO_RESERVE];
	s->value[C_TO_CWD]                      += s->value[C_CROOT_TO_CWD];

	/** nitrogen **/

	/* remove biomass from coarse root */
	s->value[CROOT_N_TO_REMOVE]             += ( s->value[CROOT_N]             * ( 1. - self_pruning_ratio ) );

	/* adding to CROOT transfer pools */
	s->value[N_CROOT_TO_RESERVE]            += s->value[CROOT_N_TO_REMOVE]     * N_FRAC_TO_RETRANSL;
	s->value[N_CROOT_TO_CWD]                += s->value[CROOT_N_TO_REMOVE]     * ( 1. - N_FRAC_TO_RETRANSL );

	s->value[N_TO_RESERVE]                  += s->value[N_CROOT_TO_RESERVE];
	s->value[N_TO_CWD]                      += s->value[N_CROOT_TO_CWD];
}
