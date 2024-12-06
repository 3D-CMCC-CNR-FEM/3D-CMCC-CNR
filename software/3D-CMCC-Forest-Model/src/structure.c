/* structure.c */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include "matrix.h"
#include "constants.h"
#include "settings.h"
#include "logger.h"
#include "g-function.h"
#include "mortality.h"
#include "structure.h"
#include "g-function.h"
#include "allometry.h"
#include "canopy_cover.h"

extern settings_t* g_settings;
extern logger_t* g_debug_log;
extern dataset_t* g_dataset;

//extern int MonthLength [];
//extern int MonthLength_Leap [];

/* */
int alloc_struct(void** t, int* count, int* avail, unsigned int size) {
	void *no_leak;

	if ( *avail ) {
		--*avail;
		++*count;
	} else {
		no_leak = realloc(*t, ++*count*size);
		if ( ! no_leak ) {
			--*count;
			return 0;
		}
		*t = no_leak;
	}
	return 1;
}

int layer_add(cell_t* const c)
{
	int ret;
	static tree_layer_t t_layer = { 0 };

	assert(c);

	ret = alloc_struct((void **)&c->tree_layers, &c->tree_layers_count, &c->t_layers_avail, sizeof(tree_layer_t));
	if ( ret )
	{
		c->tree_layers[c->tree_layers_count-1] = t_layer;
	}

	return ret;
}
/*************************************************************************************************************************/
int annual_forest_structure(cell_t* const c, const int year)
{
	int layer;
	int height;
	int dbh;
	int age;
	int species;
	int light_tol;
	int zeta_count = 0;
        int row;
        int year_dens_fin = 0;  // only used if MANAGEMENT == VAR or VAR1
 
    int tree_remove_st = 0 ;   // tree to remove because of self thinning 
         //ddalmo
       assert(g_dataset);

	species_t *s;

	/* this function compute annually:
	 * -the number of forest layers comparing the tree height values of all tree height classes
	 * -layer density
	 * -DBHDC_EFF based on overall layer cover
	 * -crown diameter
	 * -crown area
	 * -class cover
	 * -layer cover
	 * -cell cover
	 * */


	logger(g_debug_log, "\n***ANNUAL FOREST STRUCTURE***\n");

	// TODO
	// ALESSIOR
	// check if following code is useful!

    #if 0	
	assert( ! c->tree_layers_count );

	for ( height = 0; height < c->heights_count; ++height )
	{
		assert( ! c->heights[height].height_z );
		exit(1);
	}
	#endif 

    //printf(" ANNUAL STRUCTURE height_count  %d\n",c->heights_count);

	/***************************************************************************************************************/

	for ( height = 0; height < c->heights_count ; ++height )
	{
		for ( dbh = 0; dbh < c->heights[height].dbhs_count; ++dbh )
		{
			for ( age = 0; age < c->heights[height].dbhs[dbh].ages_count ; ++age )
			{
				for ( species = 0; species < c->heights[height].dbhs[dbh].ages[age].species_count; ++species )
				{

					light_tol = (int)c->heights[height].dbhs[dbh].ages[age].species[species].value[LIGHT_TOL];

					switch ( light_tol )
					{
					case 1:
						/* very shade tolerant */
						c->heights[height].dbhs[dbh].ages[age].species[species].value[MAX_LAYER_COVER] = 1.1;
						break;
					case 2:
						/* shade tolerant */
						c->heights[height].dbhs[dbh].ages[age].species[species].value[MAX_LAYER_COVER] = 1.0;
						break;
					case 3:
						/* shade intolerant */
						c->heights[height].dbhs[dbh].ages[age].species[species].value[MAX_LAYER_COVER] = 0.9;
						break;
					case 4:
						/* very shade intolerant */
						c->heights[height].dbhs[dbh].ages[age].species[species].value[MAX_LAYER_COVER] = 0.8;
						break;
					}
				}
			}
		}
	}

	/*********************************************************************************************************/
   // printf(" c->dos           = %d \n ", c->dos);
  
	/** compute number of annual layers **/

	/* note: it starts from the lowest height values up to the highest */
    //printf(" ANNUAL STRUCTURE 2 c->tree_layers_count %d\n",c->tree_layers_count);
    //printf(" ANNUAL STRUCTURE 2 c->heights_count %d\n",c->heights_count);
    // printf(" ANNUAL STRUCTURE 2 c->n_trees  %d\n",c->n_trees );

	if ( c->tree_layers_count )   // eg. in case of replanting of multiple species, we already have a layer
	{ 
     // reset and compute below, this happens when we add multiple classes during replanting
   reset_annual_layer_variables ( c );


	}	
  //printf(" ANNUAL STRUCTURE  c->tree_layers_count %d\n",c->tree_layers_count);
  // in questo modo ho messo a zero anche il tree_layers_count, ma mi serv iva per mettere a zero tutte le altre variabili
  // at layer level. TOBE FIXED and TESTED as
  ///* reset values for layer (they are recomputed below) */
  //	c->tree_layers[layer].layer_n_height_class = 0;
  //	c->tree_layers[layer].layer_n_trees        = 0;
  //	c->tree_layers[layer].layer_density        = 0;


    if ( !c->tree_layers_count )   // eg. in case of replanting of multiple species, we already have a layer
	{ 
	/* add 1 layer by default */
	if ( ! layer_add(c) ) return 0;
	}
	//printf(" ANNUAL STRUCTURE 2b c->tree_layers_count %d\n",c->tree_layers_count);

	logger(g_debug_log, "*compute height_z*\n");

	/* note: it must starts from the lowest tree height class */
	qsort(c->heights, c->heights_count, sizeof(height_t), sort_by_heights_asc);

	logger(g_debug_log, "*compute height_z*\n");

	/* compute zeta counter */

   // printf(" ANNUAL STRUCTURE 2c c->heights_count %d\n",c->heights_count);

	//printf(" ANNUAL STRUCTURE 2d  c->heights[0].value  %f\n",c->heights[0].value);
//printf(" ANNUAL STRUCTURE 2d  c->heights[1].value  %f\n",c->heights[1].value);
//if (c->heights_count>2)
//{
	//printf(" ANNUAL STRUCTURE 2d  c->heights[2].value  %f\n",c->heights[2].value);
//}
	if (c->heights_count > 1)
	{
		for ( height = 0; height < c->heights_count-1; ++height )
		{
			logger(g_debug_log, "*value %f*\n\n", c->heights[height].value);
           // printf(" ANNUAL STRUCTURE 2d  c->heights[height+1].value  %f\n",c->heights[height+1].value);
			// ALESSIOR TO ALESSIOC...this give error
			// on +1 YOU MUST remove -1 from count!
			if ( (c->heights[height+1].value - c->heights[height].value) > g_settings->tree_layer_limit )
			{
				++zeta_count;
              //     printf(" ANNUAL STRUCTURE 4 zeta_count  %d\n",zeta_count);
				/* compute layer number and alloc memory for each one */
				if ( ! layer_add(c) ) return 0;
			}
		}
	}
	   // printf(" ANNUAL STRUCTURE 3 c->tree_layers_count %d\n",c->tree_layers_count);
       //printf(" ANNUAL STRUCTURE 3 c->heights_count %d\n",c->heights_count);
       //  printf(" ANNUAL STRUCTURE 3 zeta_count  %d\n",zeta_count);
	// printf(" in structure c->tree_layers_count             = %d \n ", c->tree_layers_count);

	logger(g_debug_log, "*zeta_count %d*\n\n", zeta_count);
	logger(g_debug_log, "*c->t_layers_count %d*\n\n", c->tree_layers_count);

	// check if a regereneration layer already exist 5p7

    //if (c->tree_layers_count >1)   // more than one layer 
	//{
     // in realtaà potrei avere che la simulazione parte con il solo layer di rigenerazione 
    // decidere con Elia, come mettiamo il limite per far apparire il layer di rigenerazione

    // TEST!!!!!!!! 
	//printf(" in ANNUAL_FOREST_STRUCTURE: ETA' PIU BASSA c->cell_age_min) !! %d\n",c->cell_age_min);

    if ( (c->cell_age_min) < reg_threshold)
    {
	 c->seedl_layer=1 ;  // regeneration layer already exist or we do not let the possibility to have 
	                     // the layer, we might want to consider the smaller DBH class treshold.
     //printf("STRUCTURE 1 c->seedl_layer) !! %d\n",c->seedl_layer);
		} else {
		c->seedl_layer=0; 

	}

    //printf(" in ANNUAL_FOREST_STRUCTURE:  2 c->seedl_layer) !! %d\n",c->seedl_layer);
 
 	/*****************************************************************************************/

	/** assign zeta for each height class **/

	for ( height = c->heights_count -1 ; height >= 0; --height )
	{

		//printf("SONO IN ASSIGN ZETA \n");
		// ALESSIOR...on height == 0 you can't go to previous values (height-1)
		//if ( (c->heights[height].value - c->heights[height-1].value) > g_settings->tree_layer_limit )
		if ( height && (c->heights[height].value - c->heights[height-1].value) > g_settings->tree_layer_limit )
		{
			c->heights[height].height_z = zeta_count;
			--zeta_count;
		}
		else
		{
			c->heights[height].height_z = zeta_count;
		}
	}

	/* check */
	CHECK_CONDITION(c->tree_layers_count, <, 1);
	CHECK_CONDITION(c->tree_layers_count, >, c->heights_count);

	/*************************************************************************************/
  
   //printf(" BEOFRE COMPUTING N heigths in layer c->tree_layers[layer].layer_n_height_class %d\n", c->tree_layers[0].layer_n_height_class);

	/** compute numbers of height classes within each layer **/

	logger(g_debug_log, "*compute numbers of height classes within each layer*\n\n");

	for ( layer = c->tree_layers_count - 1; layer >= 0 ; --layer )
	{
		for ( height = c->heights_count -1; height >= 0 ; --height )
		{
			if( layer == c->heights[height].height_z )
			{
				++c->tree_layers[layer].layer_n_height_class;
			}
		}
		logger(g_debug_log, "-layer %d height class(es) = %d\n", layer, c->tree_layers[layer].layer_n_height_class);

		/* check */
		CHECK_CONDITION(c->tree_layers[layer].layer_n_height_class, <, ZERO);
	}
	logger(g_debug_log, "**************************************\n\n");

    // printf(" AFTER  COMPUTING N heigths in layer c->tree_layers[0].layer_n_height_class %d\n", c->tree_layers[0].layer_n_height_class);
    //printf(" AFTER  COMPUTING N heigths in layer c->tree_layers[1.layer_n_height_class %d\n", c->tree_layers[1].layer_n_height_class);


	/*************************************************************************************/

	/** compute numbers of trees within each layer **/



	logger(g_debug_log, "*compute numbers of trees within each layer*\n\n");

	for ( layer = c->tree_layers_count - 1; layer >= 0; --layer )
	{
		for ( height = 0; height < c->heights_count ; ++height )
		{
			//printf("in annual structure compute numbers of trees within each layer  c->heights[height] %g,\n",c->heights[height].value);
			if( layer == c->heights[height].height_z )
			{
				for ( dbh = 0; dbh < c->heights[height].dbhs_count; ++dbh )
				{


					for ( age = 0; age < c->heights[height].dbhs[dbh].ages_count ; ++age )
					{

						for ( species = 0; species < c->heights[height].dbhs[dbh].ages[age].species_count; ++species )
						{
							//if( layer == c->heights[height].height_z ) //ddalmo: there is already this condition
							//{
								s = &c->heights[height].dbhs[dbh].ages[age].species[species];

                               //  printf("in annual structure compute numbers of trees within each layer  s->counter[N_TREE] %d,\n",s->counter[N_TREE]);
								c->tree_layers[layer].layer_n_trees += s->counter[N_TREE];
							//}
						}
					}
				}
			}
		}
		logger(g_debug_log, "-layer %d number of trees = %d\n", layer, c->tree_layers[layer].layer_n_trees);
	}
	logger(g_debug_log, "**************************************\n\n");
     
	 // printf(" AFTER  COMPUTING N heigths in layer c->tree_layers[0].llayer_n_trees%d\n", c->tree_layers[0].layer_n_trees);
   // printf(" AFTER  COMPUTING N heigths in layer c->tree_layers[1].layer_n_trees %d\n", c->tree_layers[1].layer_n_trees);


	/*************************************************************************************/

	/** compute density within each layer **/

	logger(g_debug_log, "*compute density within each layer*\n\n");

	for (layer = c->tree_layers_count - 1; layer >= 0; layer --)
	{
		c->tree_layers[layer].layer_density = c->tree_layers[layer].layer_n_trees / g_settings->sizeCell;

		logger(g_debug_log, "-layer %d density = %f layer\n", layer, c->tree_layers[layer].layer_density);
	}
	logger(g_debug_log, "**************************************\n\n");

	/*************************************************************************************/

	/** compute effective dbh/crown diameter ratio within each class based on layer density (class level) **/

	logger(g_debug_log, "compute effective dbh/crown diameter ratio within each class based on layer density (class level)\n\n");

	for ( layer = c->tree_layers_count - 1; layer >= 0; --layer )
	{
		logger(g_debug_log, "----------------------------------\n");

#ifndef USE_NEW_OUTPUT
		qsort (c->heights, c->heights_count, sizeof (height_t), sort_by_heights_desc);
#endif

		for ( height = 0; height < c->heights_count ; ++height )
		{
			if( layer == c->heights[height].height_z )
			{
				for ( dbh = 0; dbh < c->heights[height].dbhs_count; ++dbh )
				{
					for ( age = 0; age < c->heights[height].dbhs[dbh].ages_count ; ++age )
					{
						for ( species = 0; species < c->heights[height].dbhs[dbh].ages[age].species_count; ++species )
						{
							dbhdc_function (c, layer, height, dbh, age, species, year);
						}
					}
				}
			}
		}
	}
	logger(g_debug_log, "**************************************\n\n");
	
	//printf("SONO ANNUAL STRUCTURE DOPO DBHDC FUNCTION ----------------------\n");
	/*************************************************************************************/

	/** Compute Crown and Canopy allometry **/

	logger(g_debug_log, "*Compute Crown and Canopy allometry*\n");

	for ( layer = c->tree_layers_count - 1; layer >= 0; --layer )
	{
		for ( height = 0; height < c->heights_count ; ++height )
		{
			if( layer == c->heights[height].height_z )
			{
				for ( dbh = 0; dbh < c->heights[height].dbhs_count; ++dbh )
				{
					for ( age = 0; age < c->heights[height].dbhs[dbh].ages_count ; ++age )
					{
						for ( species = 0; species < c->heights[height].dbhs[dbh].ages[age].species_count; ++species )
						{
							crown_allometry ( c, height, dbh, age, species );
							canopy_cover    ( c, height, dbh, age, species );

						}
					}
				}
			}
		}
	}
	logger(g_debug_log, "**************************************\n\n");

	/*************************************************************************************/

	/** compute layer canopy cover within each layer (layer level) **/

	logger(g_debug_log, "*compute layer canopy cover within each layer (layer level)*\n\n");

        // note; this variable is however not used. Yet it should be considered in the canopy radiative subroutine (e.g. competition when overlapping classes)

	for (layer = c->tree_layers_count - 1; layer >= 0; --layer)
	{
		for ( height = 0; height < c->heights_count ; ++height )
		{
		//	printf(" annual LOOP height %d!!!\n",height); 
		//	printf(" annual LOOP c->heights[height] %f!!!\n",c->heights[height].value);
			if( layer == c->heights[height].height_z )
			{
				for ( dbh = 0; dbh < c->heights[height].dbhs_count; ++dbh )
				{
					for ( age = 0; age < c->heights[height].dbhs[dbh].ages_count ; ++age )
					{
						for ( species = 0; species < c->heights[height].dbhs[dbh].ages[age].species_count; ++species )
						{
							s = &c->heights[height].dbhs[dbh].ages[age].species[species];

							c->tree_layers[layer].layer_cover_proj += s->value[CANOPY_COVER_PROJ];
							logger(g_debug_log, "layer %d cover_proj         = %f\n", layer, c->tree_layers[layer].layer_cover_proj);

						}
					}
				}
			}
		}
	}
	logger(g_debug_log, "**************************************\n");

	/*************************************************************************************/

	/** check if layer cover exceeds maximum layer cover **/


	logger(g_debug_log, "*check if layer cover exceeds maximum layer cover*\n\n");  // ddalmo: where it is checked?
                                                                                        // the DBHDC_eff is already computed fixing the maximum layer cover of the SPECIES!!

//printf(" annual structure   c->tree_layers_count %d!!!\n", c->tree_layers_count);
//printf(" annual structure  c->heights_count %d!!!\n",c->heights_count);
//printf(" annual structure c->heights[height] %d!!!\n",c->heights[height]);

	for (layer = c->tree_layers_count - 1; layer >= 0; --layer)
	{

		/* note: 04 Oct 2016 */
		/* call of this function is due to the assumption that:
			-overall layer canopy cover cannot exceeds its maximum   //DDALMO: where is this condition checked????
			-DBHDC_EFF cannot be < DBHDCMIN (otherwise self-thinning mortality is called)
		 */

		/* start to reduce DBHDC_EFF from the lowest height class */
#ifndef USE_NEW_OUTPUT
		qsort (c->heights, c->heights_count, sizeof (height_t), sort_by_heights_asc);
#endif

		for ( height = 0; height < c->heights_count ; ++height )
		{
		//	printf(" annual structure DDALMO c->heights[height] %f!!!\n",c->heights[height].value);
            // parte da height = 0 che è il più basso all interno dello stand
			if( layer == c->heights[height].height_z )
			{
				//	printf(" SONO QUA!! layer  !!! %d\n",layer);
				for ( dbh = 0; dbh < c->heights[height].dbhs_count; ++dbh )
				{
					for ( age = 0; age < c->heights[height].dbhs[dbh].ages_count ; ++age )
					{
						for ( species = 0; species < c->heights[height].dbhs[dbh].ages[age].species_count; ++species )
						{
							s = &c->heights[height].dbhs[dbh].ages[age].species[species];

							/*************** self-thinning ****************/

                                                        // compute last year of available stand density data
                                                        // note: in case of multi-class forest, for each class
                                                        // the same number of stand observations has to be provided
                                                        row = g_dataset->rows_count ;

                                                        year_dens_fin = g_dataset->rows[row-1].year_stand;


                                                        if ( c->years[year].year > year_dens_fin )
							  {
							//	printf(" annual structure self thin 3 species %s!!!\n", s->name);
							//	printf(" annual structure self thin  3 s->value[DBHDC_EFF] %f!!!\n",s->value[DBHDC_EFF]);

								if ( s->value[DBHDC_EFF] <= s->value[DBHDCMIN] )   //ddalmo may23: test include if layer_cc_proj > max_layer_cc_proj
								{
							//		printf(" annual structure self thin  LAYER %d!!!\n",layer); 
							//	printf(" ENTRO IN SELF THINNING MORTALITY !!!\n");

								// i transfer as well the number of trees to be remove 
                                    tree_remove_st = s->value[tree_remove_crowded] ; 
							//		 printf(" ENTRO IN SELF THINNING tree_remove_st %d!!!\n", tree_remove_st);
									self_thinning_mortality_new( c, layer, year , tree_remove_st);
									//self_thinning_mortality ( c, layer, year );

								}
							  }

							/* note: special case for ISIMIP, avoid self thinning when management is 'var' */
							/*if ( ( ( c->years[year].year >= g_settings->year_start_management ) && ( MANAGEMENT_VAR == g_settings->management ) )
									|| ( MANAGEMENT_OFF == g_settings->management ) )
							{

								if ( s->value[DBHDC_EFF] <= s->value[DBHDCMIN] )
								{
									self_thinning_mortality ( c, layer, year );
								}
							} */


						}
					}
				}
			}
		}
	}
	logger(g_debug_log, "**************************************\n");
	logger(g_debug_log, "**************************************\n");
	/*************************************************************************************/

	/** compute total number of trees **/

	//printf(" compute total number of trees !!!\n");


	logger(g_debug_log, "*compute total number of trees and total basal area *\n\n");

	c->n_trees = 0;

	c->basal_area = 0.;

	for ( height = 0; height < c->heights_count ; ++height )
	{

		for ( dbh = 0; dbh < c->heights[height].dbhs_count; ++dbh )
		{

			for ( age = 0; age < c->heights[height].dbhs[dbh].ages_count ; ++age )
			{
				for ( species = 0; species < c->heights[height].dbhs[dbh].ages[age].species_count; ++species )
				{
					s = &c->heights[height].dbhs[dbh].ages[age].species[species];

					c->n_trees += s->counter[N_TREE];
					c->basal_area += s->value[STAND_BASAL_AREA_m2] ; // ( ( pow( ( d->value / 2.), 2. ) ) * Pi)* 0.0001* s->counter[N_TREE];
				}
			}
		}
	}
    //printf("  c->n_trees      = %d  \n",c->n_trees);
	logger(g_debug_log,"n_trees = %d \n", c->n_trees);
	logger(g_debug_log, "**************************************\n");
	logger(g_debug_log, "**************************************\n");
	/*************************************************************************************/

	/** compute overall cell cover (cell level) **/

	logger(g_debug_log, "*daily overall cell cover (cell level))*\n");

	for ( layer = c->tree_layers_count - 1; layer >= 0; --layer )
	{
		/* note: overall cell cover can't exceed its area */
		c->cell_cover_proj += c->tree_layers[layer].layer_cover_proj;

		if ( c->cell_cover_proj > 1 + eps )
		{
			c->cell_cover_proj = 1;
		}
	}

	/*************************************************************************************/

	/** compute overall bare soil cover **/

	c->bare_soil_cover = 1. - c->cell_cover_proj;

	logger(g_debug_log, "-Number of trees cell level    = %d trees/cell\n", c->n_trees);
	logger(g_debug_log, "-Canopy cover at cell level    = %f %%\n", c->cell_cover_proj * 100.);
	logger(g_debug_log, "-Bare soil cover at cell level = %f %%\n", c->bare_soil_cover * 100.);
	logger(g_debug_log, "**************************************\n");
	logger(g_debug_log, "**************************************\n");

	return 1;
}
/*************************************************************************************************************************/
int monthly_forest_structure (cell_t* const c)
{
	return 1;
}
/*************************************************************************************************************************/

#define TEST 1 //test 1 = old 0 = new
int daily_forest_structure ( cell_t *const c, const meteo_daily_t *const meteo_daily , const int year)
{
	int layer;
	int height;
	int dbh;
	int age;
	int species;

	int layer_height_class_counter;

	tree_layer_t *l;
	height_t *h;
	species_t *s;

	/* This function computes at class level and at daily scale :
	 * -Canopy Cover Projected (based on daily LAI)
	 * -Canopy Cover Exposed (based on daily LAI)
	 * (being LAI variables at daily scale Canopy Cover must be updated)
	 */


	logger(g_debug_log, "\n***DAILY FOREST STRUCTURE***\n");


    if (c->GREFFMORT_HAPPENS ) 
	{

    //printf(" IN DAILY STRUCTURE  WHEN MORTALITy GROWTH HAPPENS    \n");

	 int layer;
	int height;
	int dbh;
	int age;
	int species;
	int light_tol;
	int zeta_count = 0;
        int row;
        int year_dens_fin = 0;  // only used if MANAGEMENT == VAR or VAR1

         //ddalmo
       assert(g_dataset);

	species_t *s;
	age_t *a;
	//height_t *h;
	//dbh_t *d;

	/* this function compute annually:
	 * -the number of forest layers comparing the tree height values of all tree height classes
	 * -layer density
	 * -DBHDC_EFF based on overall layer cover
	 * -crown diameter
	 * -crown area
	 * -class cover
	 * -layer cover
	 * -cell cover
	 * */


	logger(g_debug_log, "\n***ANNUAL FOREST STRUCTURE***\n");


	 reset_annual_layer_variables (c);

   //l->layer_n_height_class =  l->layer_n_height_class -1 ;
								//l->layer_height_class_counter --;

    //	annual_forest_structure ( c, year );

   printf(" c->doy           = %d \n ", c->doy);
    
	/***************************************************************************************************************/

	for ( height = 0; height < c->heights_count ; ++height )
	{
		for ( dbh = 0; dbh < c->heights[height].dbhs_count; ++dbh )
		{
			for ( age = 0; age < c->heights[height].dbhs[dbh].ages_count ; ++age )
			{
				for ( species = 0; species < c->heights[height].dbhs[dbh].ages[age].species_count; ++species )
				{

					light_tol = (int)c->heights[height].dbhs[dbh].ages[age].species[species].value[LIGHT_TOL];

					switch ( light_tol )
					{
					case 1:
						/* very shade tolerant */
						c->heights[height].dbhs[dbh].ages[age].species[species].value[MAX_LAYER_COVER] = 1.1;
						break;
					case 2:
						/* shade tolerant */
						c->heights[height].dbhs[dbh].ages[age].species[species].value[MAX_LAYER_COVER] = 1.0;
						break;
					case 3:
						/* shade intolerant */
						c->heights[height].dbhs[dbh].ages[age].species[species].value[MAX_LAYER_COVER] = 0.9;
						break;
					case 4:
						/* very shade intolerant */
						c->heights[height].dbhs[dbh].ages[age].species[species].value[MAX_LAYER_COVER] = 0.8;
						break;
					}
				}
			}
		}
	}

	/*********************************************************************************************************/
   //  printf(" c->dos           = %d \n ", c->dos);
 
	/** compute number of annual layers **/

	/* note: it starts from the lowest height values up to the highest */

	/* add 1 layer by default */
	if ( ! layer_add(c) ) return 0;

	logger(g_debug_log, "*compute height_z*\n");

	/* note: it must starts from the lowest tree height class */
	qsort(c->heights, c->heights_count, sizeof(height_t), sort_by_heights_asc);

	logger(g_debug_log, "*compute height_z*\n");

	/* compute zeta counter */

	if (c->heights_count > 1)
	{
		for ( height = 0; height < c->heights_count-1; ++height )
		{
			logger(g_debug_log, "*value %f*\n\n", c->heights[height].value);

			// ALESSIOR TO ALESSIOC...this give error
			// on +1 YOU MUST remove -1 from count!
			if ( (c->heights[height+1].value - c->heights[height].value) > g_settings->tree_layer_limit )
			{
				++zeta_count;

				/* compute layer number and alloc memory for each one */
				if ( ! layer_add(c) ) return 0;
			}
		}
	}
	
	//printf(" in structure c->tree_layers_count             = %d \n ", c->tree_layers_count);

	logger(g_debug_log, "*zeta_count %d*\n\n", zeta_count);
	logger(g_debug_log, "*c->t_layers_count %d*\n\n", c->tree_layers_count);


 	/*****************************************************************************************/

	/** assign zeta for each height class **/

	for ( height = c->heights_count -1 ; height >= 0; --height )
	{
		// ALESSIOR...on height == 0 you can't go to previous values (height-1)
		//if ( (c->heights[height].value - c->heights[height-1].value) > g_settings->tree_layer_limit )
		if ( height && (c->heights[height].value - c->heights[height-1].value) > g_settings->tree_layer_limit )
		{
			c->heights[height].height_z = zeta_count;
			--zeta_count;
		}
		else
		{
			c->heights[height].height_z = zeta_count;
		}
	}

	/* check */
	CHECK_CONDITION(c->tree_layers_count, <, 1);
	CHECK_CONDITION(c->tree_layers_count, >, c->heights_count);

	/*************************************************************************************/

	/** compute numbers of height classes within each layer **/

	logger(g_debug_log, "*compute numbers of height classes within each layer*\n\n");

	for ( layer = c->tree_layers_count - 1; layer >= 0 ; --layer )
	{
		for ( height = c->heights_count -1; height >= 0 ; --height )
		{
			if( layer == c->heights[height].height_z )
			{
				++c->tree_layers[layer].layer_n_height_class;
			}
		}
		logger(g_debug_log, "-layer %d height class(es) = %d\n", layer, c->tree_layers[layer].layer_n_height_class);

		/* check */
		CHECK_CONDITION(c->tree_layers[layer].layer_n_height_class, <, ZERO);
	}
	logger(g_debug_log, "**************************************\n\n");

	/*************************************************************************************/

	/** compute numbers of trees within each layer **/



	logger(g_debug_log, "*compute numbers of trees within each layer*\n\n");

	for ( layer = c->tree_layers_count - 1; layer >= 0; --layer )
	{
		for ( height = 0; height < c->heights_count ; ++height )
		{
			//printf("in annual structure compute numbers of trees within each layer  c->heights[height] %g,\n",c->heights[height].value);
			if( layer == c->heights[height].height_z )
			{
				for ( dbh = 0; dbh < c->heights[height].dbhs_count; ++dbh )
				{


					for ( age = 0; age < c->heights[height].dbhs[dbh].ages_count ; ++age )
					{

						for ( species = 0; species < c->heights[height].dbhs[dbh].ages[age].species_count; ++species )
						{
							//if( layer == c->heights[height].height_z ) //ddalmo: there is already this condition
							//{
								s = &c->heights[height].dbhs[dbh].ages[age].species[species];
                                // printf(" IN DAILY STRUCTURE  species %s!!!\n", s->name);
                                //printf("in annual structure compute numbers of trees within each layer  s->counter[N_TREE] %d,\n",s->counter[N_TREE]);
								c->tree_layers[layer].layer_n_trees += s->counter[N_TREE];
							//}
						}
					}
				}
			}
		}
		logger(g_debug_log, "-layer %d number of trees = %d\n", layer, c->tree_layers[layer].layer_n_trees);
	}
	logger(g_debug_log, "**************************************\n\n");

	/*************************************************************************************/

	/** compute density within each layer **/

	logger(g_debug_log, "*compute density within each layer*\n\n");

	for (layer = c->tree_layers_count - 1; layer >= 0; layer --)
	{
		c->tree_layers[layer].layer_density = c->tree_layers[layer].layer_n_trees / g_settings->sizeCell;
       // printf("in annual structure compute numbers of trees within each layer c->tree_layers[layer].layer_density  %f,\n",c->tree_layers[layer].layer_density );
		logger(g_debug_log, "-layer %d density = %f layer\n", layer, c->tree_layers[layer].layer_density);
	}
	logger(g_debug_log, "**************************************\n\n");

	/*************************************************************************************/
	/*************************************************************************************/
    /** Compute Crown and Canopy allometry **/
#if 0
	logger(g_debug_log, "*Compute Crown and Canopy allometry*\n");

	for ( layer = c->tree_layers_count - 1; layer >= 0; --layer )
	{
		for ( height = 0; height < c->heights_count ; ++height )
		{
			if( layer == c->heights[height].height_z )
			{
				for ( dbh = 0; dbh < c->heights[height].dbhs_count; ++dbh )
				{
					for ( age = 0; age < c->heights[height].dbhs[dbh].ages_count ; ++age )
					{
						for ( species = 0; species < c->heights[height].dbhs[dbh].ages[age].species_count; ++species )
						{
							crown_allometry ( c, height, dbh, age, species );
							canopy_cover    ( c, height, dbh, age, species );
							s = &c->heights[height].dbhs[dbh].ages[age].species[species];
							printf(" IN DAILY STRUCTURE  species %s!!!\n", s->name);
                            printf(" IN DAILY RAD s->value[CANOPY_COVER_PROJ] %f!!!\n",s->value[CANOPY_COVER_PROJ]);


						}
					}
				}
			}
		}
	}
	logger(g_debug_log, "**************************************\n\n");

#endif

	/*************************************************************************************/

	/*************************************************************************************/
    int eta_temp ;

	/** compute total number of trees **/ 
	// and daily canopy cover projection (which is equal to canopy_cover_projection)

	logger(g_debug_log, "*compute total number of trees*\n\n");

	c->n_trees = 0;

	c->basal_area = 0.;

	for ( height = 0; height < c->heights_count ; ++height )
	{
	//	  h = &c->heights[height];
		for ( dbh = 0; dbh < c->heights[height].dbhs_count; ++dbh )
		{
		//	d = &h->dbhs[dbh];
			for ( age = 0; age < c->heights[height].dbhs[dbh].ages_count ; ++age )
			{
				for ( species = 0; species < c->heights[height].dbhs[dbh].ages[age].species_count; ++species )
				{
					s = &c->heights[height].dbhs[dbh].ages[age].species[species];
					eta_temp = c->heights[height].dbhs[dbh].ages[age].value;

				//	printf("eta_temp  %d\n", eta_temp);
					
					//s = &m->cells[cell].heights[height].dbhs[dbh].ages[age].species[species];
          
		           s = &c->heights[height].dbhs[dbh].ages[age].species[species];
				   a = &c->heights[height].dbhs[dbh].ages[age]; 
                    
					canopy_cover    ( c, height, dbh, age, species );
                    //printf("  XXXX s->value[SLA_AVG1]   = %g\n", s->value[SLA_AVG1]  );
					// printf("  XXXX s->value[SLA_AVG0]   = %g\n", s->value[SLA_AVG0]  );
					//printf("CANOPY_COVER_PROJ = %f\n", s->value[CANOPY_COVER_PROJ]);
                  
				    if (c->doy >1 )
					{
                    // NB i can call daily_lai only if it is not the first day of the year
					// i call daily_lay only if class mortality occurr during the year, and not
					// because of e.g. self thinning 
                    daily_lai       ( c, a,  s );
					}  
                    
					s->value[DAILY_CANOPY_COVER_PROJ] = s->value[CANOPY_COVER_PROJ];

					c->n_trees += s->counter[N_TREE];
						c->basal_area +=  s->value[STAND_BASAL_AREA_m2] ; // ( ( pow( ( d->value / 2.), 2. ) ) * Pi)* 0.0001* s->counter[N_TREE];
				}
			}
		}
	}

	/** compute layer canopy cover within each layer (layer level) **/

	logger(g_debug_log, "*compute layer canopy cover within each layer (layer level)*\n\n");

        // note; this variable is however not used. Yet it should be considered in the canopy radiative subroutine (e.g. competition when overlapping classes)

	for (layer = c->tree_layers_count - 1; layer >= 0; --layer)
	{
		for ( height = 0; height < c->heights_count ; ++height )
		{

			if( layer == c->heights[height].height_z )
			{
				for ( dbh = 0; dbh < c->heights[height].dbhs_count; ++dbh )
				{
					for ( age = 0; age < c->heights[height].dbhs[dbh].ages_count ; ++age )
					{
						for ( species = 0; species < c->heights[height].dbhs[dbh].ages[age].species_count; ++species )
						{
							s = &c->heights[height].dbhs[dbh].ages[age].species[species];
							//  printf(" IN DAILY RAD2 c->tree_layers[layer].daily_layer_cover_proj  %f !!!\n", c->tree_layers[layer].daily_layer_cover_proj  );

                           // printf(" IN DAILY RAD2 s->value[CANOPY_COVER_PROJ] %f!!!\n",s->value[CANOPY_COVER_PROJ]);
               //printf(" IN DAILY RAD2 layer  %d!!!\n", layer );

							c->tree_layers[layer].layer_cover_proj += s->value[CANOPY_COVER_PROJ];
							logger(g_debug_log, "layer %d cover_proj         = %f\n", layer, c->tree_layers[layer].layer_cover_proj);
                 //           printf(" IN DAILY RAD3 c->tree_layers[layer].daily_layer_cover_proj  %f !!!\n", c->tree_layers[layer].layer_cover_proj );

						}
					}
				}
			}
		}
		
	}
	logger(g_debug_log, "**************************************\n");

	  //printf(" IN DAILY RAD4 c->tree_layers[0].daily_layer_cover_proj  %f !!!\n", c->tree_layers[0].layer_cover_proj );

	  //printf(" IN DAILY RAD4 c->tree_layers[1].daily_layer_cover_proj  %f !!!\n", c->tree_layers[1].layer_cover_proj );



	//printf(" IN DAILY STRUCTURE BECAUSE OF GREFF MORTALITY N TOT TREES %d \n", c->n_trees);

	logger(g_debug_log, "**************************************\n");
	logger(g_debug_log, "**************************************\n");
	/*************************************************************************************/

	/** compute overall cell cover (cell level) **/

	logger(g_debug_log, "*daily overall cell cover (cell level))*\n");

	for ( layer = c->tree_layers_count - 1; layer >= 0; --layer )
	{
		/* note: overall cell cover can't exceed its area */
		c->cell_cover_proj += c->tree_layers[layer].layer_cover_proj;

		if ( c->cell_cover_proj > 1 + eps )
		{
			c->cell_cover_proj = 1;
		}
	}

	/*************************************************************************************/

	/** compute overall bare soil cover **/

	c->bare_soil_cover = 1. - c->cell_cover_proj;

	logger(g_debug_log, "-Number of trees cell level    = %d trees/cell\n", c->n_trees);
	logger(g_debug_log, "-Canopy cover at cell level    = %f %%\n", c->cell_cover_proj * 100.);
	logger(g_debug_log, "-Bare soil cover at cell level = %f %%\n", c->bare_soil_cover * 100.);
	logger(g_debug_log, "**************************************\n");
	logger(g_debug_log, "**************************************\n");

    
	c->GREFFMORT_HAPPENS =0;

	//return 1;

   




}

//printf(" ******************** COMPUTE VERTICAL COMPETITION **********************************\n");


	/*******************************COMPUTE VERTICAL COMPETITION*******************************/
#if TEST //test layer level
	/** compute average layer tree height **/

	logger(g_debug_log, "\n***VERTICAL COMPUTATION***\n");

	for (layer = c->tree_layers_count - 1; layer >= 0; --layer)
	{
		l = &c->tree_layers[layer];

		layer_height_class_counter = 0;

		for ( height = 0; height < c->heights_count ; ++height )
		{
			if ( layer == c->heights[height].height_z )
			{
				h = &c->heights[height];
				logger(g_debug_log, "-Layer %d tree_height           = %f\n", layer, h->value);

				/* cumulate class height values */
				l->layer_avg_tree_height += h->value;

				/* increment counter */
				++layer_height_class_counter;
			}
		}

		/* compute mean */
		l->layer_avg_tree_height /= l->layer_n_height_class;

		logger(g_debug_log, "-Layer %d layer_avg_tree_height = %f\n", layer, l->layer_avg_tree_height);
		logger(g_debug_log, "-Layer %d layer_n_height_class  = %d\n", layer, l->layer_n_height_class);
	}
	logger(g_debug_log, "**************************************\n");
	logger(g_debug_log, "**************************************\n");

	/*****************************************************************************/

	/** compute TREE HEIGHT VERTICAL MODIFIER **/

	logger(g_debug_log, "\n*TREE HEIGHT VERTICAL MODIFIER*\n");

	//note: do not add to layer > 0 layer >= 0 !!!
	for (layer = c->tree_layers_count - 1; layer >= 0; --layer)
	{
		l = &c->tree_layers[layer];

		logger(g_debug_log, "layers count = %d\n", c->tree_layers_count);

		/* only if more than one layer is present */
		if ( layer > 0 )
		{
			/* following Wallace et al., 1991, Cannell and Grace 1993, Can. J. For. Res. Vol. 23 pag. 1976 */
			/* note: it considers differences in average values of tree height among all over the layers */
			/* note: it currently works among layers and not within the same class or layer */
			/* fixme: it could works also within the same layer, to be fixed once... */

			logger(g_debug_log, "-Layer %d minus layer %d = %f\n", layer, layer - 1, c->tree_layers[layer].layer_avg_tree_height - c->tree_layers[layer-1].layer_avg_tree_height);

			/* upper layer */
			l->layer_tree_height_modifier = 0.5 * ( 1. + pow ( 2. , ( - c->tree_layers[layer-1].layer_avg_tree_height / c->tree_layers[layer].layer_avg_tree_height ) ) -
					( pow ( 2. , ( - c->tree_layers[layer].layer_avg_tree_height / c->tree_layers[layer-1].layer_avg_tree_height ) ) ) );
			logger(g_debug_log, "-Layer %d light vertical modifier = %f\n", layer, l->layer_tree_height_modifier);

			/* lower layer */
			c->tree_layers[layer-1].layer_tree_height_modifier = 1. - l->layer_tree_height_modifier;
			logger(g_debug_log, "-Layer %d light vertical modifier = %f\n", layer - 1, c->tree_layers[layer-1].layer_tree_height_modifier);

		}
		if( c->tree_layers_count == 1 )
		{
			/* no vertical competition */
			l->layer_tree_height_modifier = 1.;
			logger(g_debug_log, "-Only one layer %d light vertical modifier = %f\n", layer, l->layer_tree_height_modifier);
		}
	}
#else //test tree height class level 22 May 2016

	/* following Wallace et al., 1991, Cannell and Grace 1993, Can. J. For. Res. Vol. 23 pag. 1976 */
	/* note: it considers differences in average values of tree height among all over the layers */
	/* note: it currently works among layers and not within the same class or layer */


	/* sort heights in descending order */
	qsort ( c->heights, c->heights_count, sizeof (height_t), sort_by_heights_desc );

	// ALESSIOR TO ALESSIOC...this give error
	// on +1 YOU MUST remove -1 from count!
	for ( height = 0; height < c->heights_count-1; ++height )
	{
		for ( dbh = 0; dbh < c->heights[height].dbhs_count; ++dbh )
		{
			for ( age = 0; age < c->heights[height].dbhs[dbh].ages_count ; ++age )
			{
				for ( species = 0; species < c->heights[height].dbhs[dbh].ages[age].species_count; ++species )
				{
					if ( c->heights_count > 1 && c->heights[height+1].dbhs[dbh].ages[age].species[species].value[LAI_PROJ] )
					{
						/* upper height class */
						c->heights[height].dbhs[dbh].ages[age].species[species].value[F_LIGHT_VERT] = 0.5 * ( 1. + pow ( 2. , ( - c->heights[height+1].value / c->heights[height].value ) ) -
								( pow ( 2. , ( - c->heights[height].value / c->heights[height+1].value ) ) ) );

						/* lower height class */
						c->heights[height+1].dbhs[dbh].ages[age].species[species].value[F_LIGHT_VERT] = 1 - c->heights[height].dbhs[dbh].ages[age].species[species].value[F_LIGHT_VERT];

						logger(g_debug_log, "-%s light vertical modifier = %f\n", c->heights[height].dbhs[dbh].ages[age].species[species].name, c->heights[height].dbhs[dbh].ages[age].species[species].value[F_LIGHT_VERT]);
						logger(g_debug_log, "-%s light vertical modifier = %f\n", c->heights[height+1].dbhs[dbh].ages[age].species[species].name, c->heights[height+1].dbhs[dbh].ages[age].species[species].value[F_LIGHT_VERT]);
					}
					else
					{
						c->heights[height].dbhs[dbh].ages[age].species[species].value[F_LIGHT_VERT] = 1.;
						logger(g_debug_log, "-Only one Height light vertical modifier = %f\n", c->heights[height].dbhs[dbh].ages[age].species[species].value[F_LIGHT_VERT]);
					}
				}
			}
		}
	}
#endif
	logger(g_debug_log, "**************************************\n");
	logger(g_debug_log, "**************************************\n");

	/******************************COMPUTE HORIZONTAL COMPETITION******************************/
//printf(" ******************** COMPUTE HORIZONTAL COMPETITION **********************************\n");


	/** overall cell level within each layer **/

	logger(g_debug_log, "\n***HORIZONTAL COMPUTATION***\n");

	for (layer = c->tree_layers_count - 1; layer >= 0; --layer)
	{
		l = &c->tree_layers[layer];

		for ( height = 0; height < c->heights_count ; ++height )
		{
			if ( layer == c->heights[height].height_z )
			{
				h = &c->heights[height];

				for ( dbh = 0; dbh < c->heights[height].dbhs_count; ++dbh )
				{
					for ( age = 0; age < c->heights[height].dbhs[dbh].ages_count ; ++age )
					{
						for ( species = 0; species < c->heights[height].dbhs[dbh].ages[age].species_count; ++species )
						{
							s = &c->heights[height].dbhs[dbh].ages[age].species[species];

							/* check for Veg period otherwise don't consider the class */
							if ( s->value[LAI_PROJ] )
							{
								/*******************************************************************************************/
								/* compute daily canopy projected cover */
								s->value[DAILY_CANOPY_COVER_PROJ] = s->value[CANOPY_COVER_PROJ];
  //                                printf(" IN DAILY STRUCTURE COMPETITION s->value[CANOPY_COVER_PROJ]  %f !!!\n", s->value[CANOPY_COVER_PROJ] );


								/* canopy cannot absorbs more than 100% of incoming flux (e.g. rain) */
								if( s->value[DAILY_CANOPY_COVER_PROJ] > 1. ) s->value[DAILY_CANOPY_COVER_PROJ] = 1.;

								/*******************************************************************************************/

								/* compute daily canopy exposed cover */
								//s->value[DAILY_CANOPY_COVER_EXP] = s->value[CANOPY_COVER_EXP];

								/* canopy cannot absorb more than 100% of incoming flux (e.g. light) */
								//if( s->value[DAILY_CANOPY_COVER_EXP] > 1. ) s->value[DAILY_CANOPY_COVER_EXP] = 1.;

								/*******************************************************************************************/
#if TEST
								/* integrate with layer TREE HEIGHT MODIFIER */
								//if ( layer == c->heights[height].height_z ) s->value[DAILY_CANOPY_COVER_EXP] *= l->layer_tree_height_modifier ;
#endif
								/*******************************************************************************************/

								/* canopy cannot absorb more than 100% of incoming flux (e.g. light) */
								//if( s->value[DAILY_CANOPY_COVER_EXP] > 1. ) s->value[DAILY_CANOPY_COVER_EXP] = 1.;

								/*******************************************************************************************/

								/* sum all over canopy cover projected for each layer */
								c->tree_layers[layer].daily_layer_cover_proj += s->value[DAILY_CANOPY_COVER_PROJ];
    //                             printf(" IN DAILY STRUCTURE COMPETITION c->tree_layers[layer].daily_layer_cover_proj %f !!!\n", c->tree_layers[layer].daily_layer_cover_proj);


								/* sum all over canopy cover exposed for each layer */
								//c->tree_layers[layer].daily_layer_cover_exp  += s->value[DAILY_CANOPY_COVER_EXP];
							}

							logger(g_debug_log, "-Species %s DAILY_CANOPY_COVER_PROJ = %f\n", s->name, s->value[DAILY_CANOPY_COVER_PROJ]);
							//logger(g_debug_log, "-Species %s DAILY_CANOPY_COVER_EXP  = %f\n", s->name, s->value[DAILY_CANOPY_COVER_EXP]);

							logger(g_debug_log, "-Layer %d daily_layer_cover_proj = %f\n", layer, l->daily_layer_cover_proj);
							//logger(g_debug_log, "-Layer %d daily_layer_cover_exp  = %f\n", layer, l->daily_layer_cover_exp);
						}
					}
				}
			}
		}
	}
	logger(g_debug_log, "**************************************\n");
	logger(g_debug_log, "**************************************\n");

//	/** compute fraction of class cover (exposed and projected) **/
//
//	for (layer = c->tree_layers_count - 1; layer >= 0; --layer)
//	{
//		l = &c->tree_layers[layer];
//
//		for ( height = 0; height < c->heights_count ; ++height )
//		{
//			if ( layer == c->heights[height].height_z )
//			{
//				h = &c->heights[height];
//
//				for ( dbh = 0; dbh < c->heights[height].dbhs_count; ++dbh )
//				{
//					for ( age = 0; age < c->heights[height].dbhs[dbh].ages_count ; ++age )
//					{
//						for ( species = 0; species < c->heights[height].dbhs[dbh].ages[age].species_count; ++species )
//						{
//							if ( s->value[LAI_PROJ] )
//							{
//								int i;
//								int odd;
//								double diff_cover;
//								double max_cover;
//								double min_cover;
//								double temp_diff_cover;
//								double temp_can_cover;
//								double cum_temp_cover = 0.;
//								int daylength = ROUND(meteo_daily->daylength);
//
//								s = &c->heights[height].dbhs[dbh].ages[age].species[species];
//
//								/*****************************************************************************************/
//
//
//								logger(g_debug_log, "\n..integrating all over the day\n");
//
//								//fixme IT MUST BE INTEGRATED ALL OVER THE LAYER!!!!!!!!!!
//
//								/** now integrating all over the daylength **/
//								max_cover       = s->value[DAILY_CANOPY_COVER_EXP];
//								min_cover       = s->value[DAILY_CANOPY_COVER_PROJ];
//								diff_cover      = max_cover - min_cover;
//								temp_diff_cover = diff_cover / ( meteo_daily->daylength / 2. );
//
//								/* check if odd */
//								if(daylength %2 != 0) odd = 0; /* ok is a odd number */
//								else odd = 1;                  /* in NOT an odd number */
//
//								/* morning to nadir */
//								for ( i = 0; i < daylength ; ++i)
//								{
//									if ( ! odd )
//									{
//										if ( i == 0 )
//										{
//											temp_can_cover  = max_cover;
//										}
//										else if ( i <= ( daylength  / 2. ) )
//										{
//											/* morning to nadir */
//											temp_can_cover -= ( temp_diff_cover );
//										}
//										else
//										{
//											/* from nadir to evening */
//											temp_can_cover += ( temp_diff_cover );
//										}
//									}
//									else
//									{
//										if ( i == 0 )
//										{
//											temp_can_cover  = max_cover;
//										}
//										else if ( i <= ( daylength  / 2. ) -1 )
//										{
//											/* morning to nadir */
//											temp_can_cover -= ( temp_diff_cover );
//										}
//										else if (i <= ( daylength  / 2. ))
//										{
//											/* nothing to do */
//										}
//										else
//										{
//											/* from nadir to evening */
//											temp_can_cover += ( temp_diff_cover );
//										}
//									}
//									/* cumulate */
//									cum_temp_cover += temp_can_cover;
//								}
//
//								/** compute integrating all over the day length corresponding canopy intercepting cover **/
//								s->value[DAILY_CANOPY_COVER_EXP] = cum_temp_cover / daylength;
//							}
//
//							logger(g_debug_log, "-Species %s DAILY_CANOPY_COVER_PROJ = %f\n", s->name, s->value[DAILY_CANOPY_COVER_PROJ]);
//							logger(g_debug_log, "-Species %s DAILY_CANOPY_COVER_EXP  = %f\n", s->name, s->value[DAILY_CANOPY_COVER_EXP]);
//
//							/* check */
//							CHECK_CONDITION ( s->value[DAILY_CANOPY_COVER_EXP], > , 1 )
//						}
//					}
//				}
//			}
//		}
//	}

	/******************************COUPLING HORIZONTAL+VERTICAL******************************/
#if TEST
#else

	/** coupling vertical and horizontal modifiers **/

	logger(g_debug_log, "\n..coupling vertical modifier with canopy cover exp\n");

	for ( height = 0; height < c->heights_count ; ++height )
	{
		for ( dbh = 0; dbh < c->heights[height].dbhs_count; ++dbh )
		{
			for ( age = 0; age < c->heights[height].dbhs[dbh].ages_count ; ++age )
			{
				for ( species = 0; species < c->heights[height].dbhs[dbh].ages[age].species_count; ++species )
				{
					s = &c->heights[height].dbhs[dbh].ages[age].species[species];

					s->value[DAILY_CANOPY_COVER_EXP] *=s->value[F_LIGHT_VERT];
					logger(g_debug_log, "-Species %s test_daily_conopy_cover_exp  = %f\n", s->name, s->value[DAILY_CANOPY_COVER_EXP]);
				}
			}
		}
	}
	logger(g_debug_log, "**************************************\n");
	logger(g_debug_log, "**************************************\n");
#endif
	return 1;
}


