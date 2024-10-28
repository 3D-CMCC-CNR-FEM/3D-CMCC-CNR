/*new_forest_tree_class.c*/

/* includes */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "new_forest_tree_class.h"
#include "matrix.h"
#include "management.h"
#include "constants.h"
#include "settings.h"
#include "common.h"
#include "logger.h"
#include "initialization.h"
#include "structure.h"
#include "biomass.h"
#include "allometry.h"
#include "g-function.h"
#include "print.h"

extern settings_t* g_settings;
extern logger_t* g_debug_log;

static int fill_cell_for_replanting(cell_t *const c, const int species_index)
{
	char* p;
	height_t* h;
	dbh_t* d;
	age_t* a;

	static height_t height   = { 0 };
	static dbh_t dbh         = { 0 };
	static age_t age         = { 0 };
	static species_t species = { 0 };

	assert(c);

	/* alloc memory for heights */
	if ( ! alloc_struct((void **)&c->heights, &c->heights_count, &c->heights_avail, sizeof(height_t)) )
	{
		return 0;
	}
	c->heights[c->heights_count-1]       = height;
	c->heights[c->heights_count-1].value = g_settings->replanted[species_index].height;
	h = &c->heights[c->heights_count-1];

	/* alloc memory for dbhs */
	if ( ! alloc_struct((void **)&h->dbhs, &h->dbhs_count, &h->dbhs_avail, sizeof(dbh_t)) )
	{
		return 0;
	}
	h->dbhs[h->dbhs_count-1]       = dbh;
	h->dbhs[h->dbhs_count-1].value = g_settings->replanted[species_index].avdbh;
	d = &h->dbhs[h->dbhs_count-1];

	/* alloc memory for ages */
	if ( ! alloc_struct((void **)&d->ages, &d->ages_count, &d->ages_avail, sizeof(age_t)) )
	{
		return 0;
	}
	d->ages[d->ages_count-1]       = age;
	d->ages[d->ages_count-1].value = (int)g_settings->replanted[species_index].age;
	a = &d->ages[d->ages_count-1];

	/* alloc memory for species */
	if ( ! alloc_struct((void **)&a->species, &a->species_count, &a->species_avail, sizeof(species_t)) )
	{
		return 0;
	}

	p = string_copy(g_settings->replanted[species_index].species);
	if ( ! p ) return 0;

	a->species[a->species_count-1]                  = species;
	a->species[a->species_count-1].management       = T;
	a->species[a->species_count-1].name             = p;
	a->species[a->species_count-1].counter[N_TREE]  = (int)g_settings->replanted[species_index].n_tree;
	a->species[a->species_count-1].counter[N_STUMP] = 0;
	a->species[a->species_count-1].value[LAI_PROJ]  = g_settings->replanted[species_index].lai;

	return 1;
}

int add_tree_class_for_replanting (cell_t *const c, const int day, const int month, const int year, const int rsi)
{
	int height;
	int dbh;
	int age;
	int species;
	int day_temp;
	int month_temp;
	int DaysInMonth [] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };

	logger(g_debug_log, "**ADD NEW TREE CLASS (REPLANTATION)**\n");
      
	/* it is used only with "human" regeneration */
	logger(g_debug_log, "Human management\n");

	if ( ! fill_cell_for_replanting ( c, rsi ) ) return 0;  

	height = c->heights_count - 1;
	dbh = c->heights[height].dbhs_count - 1;
	age = c->heights[height].dbhs[dbh].ages_count - 1;
	species = c->heights[height].dbhs[dbh].ages[age].species_count - 1;

	  
	species_t *s;
	s = &c->heights[height].dbhs[dbh].ages[age].species[species];


	/* fill with species values from parameterization file */
	if ( ! fill_species_from_file ( &c->heights[height].dbhs[dbh].ages[age].species[species]) )
	{
		return 0;
	}   

	/* check for veg days */
	for (month_temp = 0; month_temp < 12; ++month_temp)
	{
		/* for handling leap years */
		int days_per_month;

		days_per_month = DaysInMonth[month_temp];
		if ( (FEBRUARY == month_temp) && IS_LEAP_YEAR(c->years[year].year) )
		{
			++days_per_month;
		}

		for ( day_temp = 0; day_temp < days_per_month; ++day_temp )
		{
			/* compute annually the days for the growing season before any other process */
			Veg_Days ( c , day_temp, month_temp, year );
		}
	}
	/* initialize new power function */
	allometry_power_function           ( c );

	/* initialize new carbon pool fraction */
	carbon_pool_fraction               ( c );

	/* initialize new forest structure */
	initialization_forest_structure    (c , day, month, year);
  
	/* initialize new forest class pools */
	initialization_forest_class_C      ( c, height, dbh, age, species );
       
	/* initialize new nitrogen pools */
	initialization_forest_class_N      ( c, height, dbh, age, species );

	/* initialize new litter pools */
	//initialization_forest_class_litter ( c, height, dbh, age, species );

        //ddalmo august 2021 
        /* update forest cell pool */
	initialization_forest_cell_C ( c, height, dbh, age, species ); 

	 c->annual_Ctree_repl +=  (( (s->value[STEM_C]  * 1e6 )/ g_settings->sizeCell)
                            + ((s->value[CROOT_C]  * 1e6 )/ g_settings->sizeCell)
                            + ((s->value[FROOT_C]  * 1e6 )/ g_settings->sizeCell)
							+ ((s->value[BRANCH_C]  * 1e6 )/ g_settings->sizeCell)
							+ ((s->value[FRUIT_C]  * 1e6 )/ g_settings->sizeCell)
                            + ((s->value[LEAF_C]  * 1e6 )/ g_settings->sizeCell)
		                    + ((s->value[RESERVE_C]  * 1e6 )/ g_settings->sizeCell)
	                          )    ; 
							  

	/* print new forest class dataset */
	print_new_daily_forest_class_data  ( c, height, dbh, age, species );

	return 1;
}
/************************************************************************************************************/
//  5p6 August 2021
//  it is included the possibility to mimic the regeneration via replanting process (i.e. adding a new class of saplings).
//  the data of the new class are taken from the regeneration 'settings'.
//  Currently it is advisable to add the regeneration class to an existing layer alone, as the multilayer option has not been tested yet (FIXME)
//  Currently it is assumed that regeneration setting is used in combination with management (e.g.shelterwood)
//  hence the information about when to add the new young (saplings -DBH age) class is provided in the management.txt file
//  The only crucial variable to be provided is the density of the saplings.

static int fill_cell_for_replanting_reg(cell_t *const c )
{
	char* p;
	height_t* h;
	dbh_t* d;
	age_t* a;

	static height_t height   = { 0 };
	static dbh_t dbh         = { 0 };
	static age_t age         = { 0 };
	static species_t species = { 0 };

	assert(c);

	/* alloc memory for heights */
	if ( ! alloc_struct((void **)&c->heights, &c->heights_count, &c->heights_avail, sizeof(height_t)) )
	{
		return 0;
	}
	c->heights[c->heights_count-1]       = height;
	c->heights[c->heights_count-1].value =  g_settings->regeneration_height; //g_settings->replanted[species_index].height;
	h = &c->heights[c->heights_count-1];

	/* alloc memory for dbhs */
	if ( ! alloc_struct((void **)&h->dbhs, &h->dbhs_count, &h->dbhs_avail, sizeof(dbh_t)) )
	{
		return 0;
	}
	h->dbhs[h->dbhs_count-1]       = dbh;
	h->dbhs[h->dbhs_count-1].value =  g_settings->regeneration_avdbh; //g_settings->replanted[species_index].avdbh;
	d = &h->dbhs[h->dbhs_count-1];

	/* alloc memory for ages */
	if ( ! alloc_struct((void **)&d->ages, &d->ages_count, &d->ages_avail, sizeof(age_t)) )
	{
		return 0;
	}
	d->ages[d->ages_count-1]       = age;
	d->ages[d->ages_count-1].value =  (int)g_settings->regeneration_age; //(int)g_settings->replanted[species_index].age;
	a = &d->ages[d->ages_count-1];

	/* alloc memory for species */
	if ( ! alloc_struct((void **)&a->species, &a->species_count, &a->species_avail, sizeof(species_t)) )
	{
		return 0;
	}

	p =string_copy(g_settings->regeneration_species); //string_copy(g_settings->replanted[species_index].species);
	if ( ! p ) return 0;

	a->species[a->species_count-1]                  = species;
	a->species[a->species_count-1].management       = T;
	a->species[a->species_count-1].name             = p;
	a->species[a->species_count-1].counter[N_TREE]  =(int)g_settings->regeneration_n_tree; //(int)g_settings->replanted[species_index].n_tree;
	a->species[a->species_count-1].counter[N_STUMP] = 0;
	a->species[a->species_count-1].value[LAI_PROJ]  = g_settings->regeneration_lai; //g_settings->replanted[species_index].lai;

	return 1;
}

int add_tree_class_for_replanting_reg (cell_t *const c, const int day, const int month, const int year)
{
	int height;
	int dbh;
	int age;
	int species;
	int day_temp;
	int month_temp;
	int DaysInMonth [] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };

	logger(g_debug_log, "**ADD NEW TREE CLASS (REPLANTATION AS GENERATION)**\n");
      
	/* it is used only mimicing natural regeneration by mean of a replanting function */
	logger(g_debug_log, "Replanting as regeneration\n");

	if ( ! fill_cell_for_replanting_reg ( c ) ) return 0;  

	height = c->heights_count - 1;
	dbh = c->heights[height].dbhs_count - 1;
	age = c->heights[height].dbhs[dbh].ages_count - 1;
	species = c->heights[height].dbhs[dbh].ages[age].species_count - 1;
  

  species_t *s;
	s = &c->heights[height].dbhs[dbh].ages[age].species[species];

        // comment: at this point the new class simply has been added and has 
        // an height index = c->heights_count - 1

	/* fill with species values from parameterization file */
	if ( ! fill_species_from_file ( &c->heights[height].dbhs[dbh].ages[age].species[species]) )
	{
		return 0;
	}   

	/* check for veg days */
	for (month_temp = 0; month_temp < 12; ++month_temp)
	{
		/* for handling leap years */
		int days_per_month;

		days_per_month = DaysInMonth[month_temp];
		if ( (FEBRUARY == month_temp) && IS_LEAP_YEAR(c->years[year].year) )
		{
			++days_per_month;
		}

		for ( day_temp = 0; day_temp < days_per_month; ++day_temp )
		{
			/* compute annually the days for the growing season before any other process */
			Veg_Days ( c , day_temp, month_temp, year );
		}
	}

	/* initialize new power function */
	allometry_power_function           ( c );

	/* initialize new carbon pool fraction */
	carbon_pool_fraction               ( c );

	/* initialize new forest structure */ 
	initialization_forest_structure    (c , day, month, year);

        // comment: in the initialization_forest_structure the new added layer has height index = 0 and the 
        // dominant layers index >= 1 ordered according to descending height

        height= 0;  // so to initialize the new layer only! Which is the new added 'regeneration-layer'

	/* initialize new forest class pools */
	initialization_forest_class_C      ( c, height, dbh, age, species );

#if 0
	printf(" ADD REG LAYER c->value[STEM_C]  %f!!!\n", (s->value[STEM_C]  * 1e6 )/ g_settings->sizeCell);
    printf(" ADD REG LAYER c->value[CROOT_C]  %f!!!\n", (s->value[CROOT_C]  * 1e6 )/ g_settings->sizeCell);
	printf(" ADD REG LAYER c->value[BRANCH_C]  %f!!!\n", (s->value[BRANCH_C]  * 1e6 )/ g_settings->sizeCell);
     printf(" ADD REG LAYER c->value[RESERVE_C]  %f!!!\n", (s->value[RESERVE_C]  * 1e6 )/ g_settings->sizeCell);
    printf(" ADD REG LAYER c->value[LEAF_C]  %f!!!\n", (s->value[LEAF_C]  * 1e6 )/ g_settings->sizeCell);
    printf(" ADD REG LAYER c->value[FROOT_C]  %f!!!\n", (s->value[FROOT_C]  * 1e6 )/ g_settings->sizeCell);
printf(" ADD REG LAYER c->value[FRUIT_C]  %f!!!\n", (s->value[FRUIT_C]  * 1e6 )/ g_settings->sizeCell);
#endif
	/* initialize new nitrogen pools */
	initialization_forest_class_N      ( c, height, dbh, age, species );

	/* initialize new litter pools */
	//initialization_forest_class_litter ( c, height, dbh, age, species );

        //ddalmo august 2021 
        /* update forest cell pool*/
		// adding the new class pool
	initialization_forest_cell_C ( c, height, dbh, age, species ); 

   #if 0
	printf(" ADD REG LAYER  c->leaf_carbon   %f!!!\n",c->leaf_carbon  );
			printf("ADD REG LAYER    c->croot_carbon   %f!!!\n",c->croot_carbon   );
			printf(" ADD REG LAYER   c->froot_carbon   %f!!!\n",c->froot_carbon   );
			printf("ADD REG LAYER   c->stem_carbon   %f!!!\n",c->stem_carbon  );
			printf(" ADD REG LAYER   c->branch_carbon   %f!!!\n",c->branch_carbon  );
			printf(" ADD REG LAYER    c->reserve_carbon   %f!!!\n",c->reserve_carbon  );
			printf(" ADD REG LAYER    c->fruit_carbon   %f!!!\n",c->fruit_carbon  );
#endif
    c->annual_Ctree_repl +=  (( (s->value[STEM_C]  * 1e6 )/ g_settings->sizeCell)
                            + ((s->value[CROOT_C]  * 1e6 )/ g_settings->sizeCell)
                            + ((s->value[FROOT_C]  * 1e6 )/ g_settings->sizeCell)
							+ ((s->value[BRANCH_C]  * 1e6 )/ g_settings->sizeCell)
							+ ((s->value[FRUIT_C]  * 1e6 )/ g_settings->sizeCell)
                            + ((s->value[LEAF_C]  * 1e6 )/ g_settings->sizeCell)
		                    + ((s->value[RESERVE_C]  * 1e6 )/ g_settings->sizeCell)
	                          )    ; 
	/* print new forest class dataset */
	print_new_daily_forest_class_data  ( c, height, dbh, age, species );
      
	return 1;
}




/************************************************************************************************************/

static int fill_cell_for_regeneration( cell_t *const c )
{
	char* p;

	height_t* h;
	dbh_t* d;
	age_t* a;

	static height_t height = { 0 };
	static dbh_t dbh = { 0 };
	static age_t age = { 0 };
	static species_t species = { 0 };

	assert(c);

	logger(g_debug_log, "\n**Fill cell from regeneration**\n");

	/* alloc memory for heights */
	if ( ! alloc_struct((void **)&c->heights, &c->heights_count, &c->heights_avail, sizeof(height_t)) )
	{
		return 0;
	}
	c->heights[c->heights_count-1] = height;
	c->heights[c->heights_count-1].value = g_settings->regeneration_height;
	h = &c->heights[c->heights_count-1];

	/* alloc memory for dbhs */
	if ( ! alloc_struct((void **)&h->dbhs, &h->dbhs_count, &h->dbhs_avail, sizeof(dbh_t)) )
	{
		return 0;
	}
	h->dbhs[h->dbhs_count-1] = dbh;
	h->dbhs[h->dbhs_count-1].value = g_settings->regeneration_avdbh;
	d = &h->dbhs[h->dbhs_count-1];

	/* alloc memory for ages */
	if ( ! alloc_struct((void **)&d->ages, &d->ages_count, &d->ages_avail, sizeof(age_t)) )
	{
		return 0;
	}
	d->ages[d->ages_count-1] = age;
	d->ages[d->ages_count-1].value = (int)g_settings->regeneration_age;
	a = &d->ages[d->ages_count-1];

	/* alloc memory for species */
	if ( ! alloc_struct((void **)&a->species, &a->species_count, &a->species_avail, sizeof(species_t)) )
	{
		return 0;
	}

	p = string_copy(g_settings->regeneration_species);
	if ( ! p ) return 0;

	a->species[a->species_count-1]                  = species;
	a->species[a->species_count-1].management       = T;
	a->species[a->species_count-1].name             = p;
	a->species[a->species_count-1].counter[N_TREE]  = (int)g_settings->regeneration_n_tree;
	a->species[a->species_count-1].counter[N_STUMP] = 0;
	a->species[a->species_count-1].value[LAI_PROJ]  = g_settings->regeneration_lai;
	logger(g_debug_log, "-height = %f\n", c->heights[c->heights_count-1].value);
	logger(g_debug_log, "-dbh    = %f\n", h->dbhs[h->dbhs_count-1].value);
	logger(g_debug_log, "-age    = %d\n", d->ages[d->ages_count-1].value);
	logger(g_debug_log, "-name   = %s\n", a->species[a->species_count-1].name);
	logger(g_debug_log, "-N-tree = %d\n", a->species[a->species_count-1].counter[N_TREE]);
	logger(g_debug_log, "-lai    = %f\n", a->species[a->species_count-1].value[LAI_PROJ]);


	return 1;
}


int add_tree_class_for_regeneration ( cell_t *const c )
{
	logger(g_debug_log, "**ADD NEW TREE CLASS (REGENERATION)**\n");

	/* it is used only with natural regeneration */
	if ( ! fill_cell_for_regeneration( c ) ) return 0;

	/* fill with species values from parameterization file */
	if ( ! fill_species_from_file(&c->heights[c->heights_count-1].dbhs[0].ages[0].species[0]) )
	{
		return 0;
	}

	/* initialize power function */
	allometry_power_function           ( c );

	/* initialize carbon pool fraction */
	carbon_pool_fraction               ( c );

	/* initialize new forest class pools */
	initialization_forest_class_C      ( c, c->heights_count-1, 0, 0, 0 );

	/* initialize nitrogen pools */
	initialization_forest_class_N      ( c, c->heights_count-1, 0, 0, 0 );

	/* print new forest class dataset */
	print_new_daily_forest_class_data  ( c, c->heights_count-1, 0, 0, 0 );

	return 1;
}


