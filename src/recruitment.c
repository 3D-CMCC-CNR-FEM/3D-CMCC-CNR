/* includes */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "matrix.h"
#include "constants.h"
#include "settings.h"
#include "structure.h"
#include "logger.h"
#include "common.h"
#include "initialization.h"
#include "biomass.h"
#include "allometry.h"
#include "g-function.h"
#include "print.h"
#include "C-fruit-partitioning.h"
#include "regeneration.h"
#include "recruitment.h"
#include "management.h"
#include "remove_tree_class.h"

extern management_t* g_management;
extern settings_t* g_settings;
extern logger_t* g_debug_log;
extern dataset_t* g_dataset;

#if 0
/**/
void Establishment_LPJ (cell_t *const c, const int layer, const int height, const int dbh, const int age, const int species,const meteo_daily_t *const meteo_daily)
{

	//double FProCov;   //LPJ Foliage Projective  Cover for Seed Establishment
	double EstabRate;   //Seed establishment rate Under Dominant Canopy
	static int Nsapling;

	species_t *s;
	s = &c->heights[height].dbhs[dbh].ages[age].species[species];

	logger(g_debug_log, "\n--LPJ ESTABLISHMENT--\n");

	logger(g_debug_log, "Intrinsic Germinability Rate = %f %% \n", s->value[GERMCAPACITY] * 100);

	logger(g_debug_log, "Annual Number of seeds using LPJ  = %d seeds/ha/year\n", s->counter[N_SEED]);

	//da ricontrollare per vedere quale Lai prende di quale layer e semmai non utilizzare il FProCov

	/*
	   LightAbsorb = 1 -  (exp(- s->value[K] * s->value[LAI]));

	   logger(g_debug_log, "Light Absorb in Dominant Canopy = %f \n", LightAbsorb);
	   logger(g_debug_log, "Canopy Cover in  Dominant layer with DBHDC function = %f \n", canopy_cover_dominant);
	   FProCov = canopy_cover_dominant * LightAbsorb;
	   logger(g_debug_log, "Foliage Projective Cover = %f \n", FProCov);


	   logger(g_debug_log, "LPJ Fractional Projective Cover FPC = %f \n", FProCov);

	   EstabRate = s->value[GERMCAPACITY] * (1 - exp((-5) * (1 - FProCov))) * (1 - FProCov);
	   logger(g_debug_log, "Seed Establishment Rate from LPJ = %f saplings/m^2 \n", EstabRate);
	 */

	EstabRate = (s->value[GERMCAPACITY] * (1 - exp((-5) * (1 - meteo_daily->par))) * (1 - meteo_daily->par)) / s->counter[N_TREE];
	logger(g_debug_log, "Light absorb for establishment = %f \n", meteo_daily->par);
	logger(g_debug_log, "Seed Establishment Rate from LPJ = %f saplings/m^2 \n", EstabRate);

	s->counter[N_TREE_SAP] = (int)(s->counter[N_SEED] * EstabRate);
	logger(g_debug_log, "Annual Number of Saplings per hectare using LPJ = %d Saplings/year hectare\n", Nsapling);
	logger(g_debug_log, "Annual Number of Saplings using LPJ = %f Saplings/year m^2\n", (double) Nsapling / g_settings->sizeCell );
	logger(g_debug_log, "Percentage of seeds survived using LPJ = %f %% seeds/year hectare\n", ((double)Nsapling * 100)/(double)s->counter[N_SEED] );
}

#endif // 0

 /*********************************** SAPONARO 04/2023 RECRUITMENT MODULE*********************************************/

#if 0      //***********This is just a test to allocation dynamic memory (DO NOT USE) *********************//

int recruitment (cell_t *const c, const int day, const int month, const int year) {

   int Seedling_surv;
   int* Seedling_ptr;

   float S_height;
   float* Seedling_h_ptr;

   float S_dbh;
   float* Seedling_d_prt;

   int S_age;
   int* Seedling_a_ptr;

   int S_species;
   int* Seedling_s_ptr;


    /*alloc memory for heights */
    Seedling_ptr = (int*) malloc(Seedling_surv * sizeof(int));

	if (Seedling_ptr == NULL) {
        printf("Memory not allocated.\n");
        exit(0);
    }
    else {

    // Memory has been successfully allocated
    printf("Memory successfully allocated using malloc.\n");
    }

	/* alloc memory for heights */
	Seedling_h_ptr = (float*) malloc(S_height * sizeof(float));

	if (Seedling_h_ptr == NULL) {
        printf("Memory not allocated.\n");
        exit(0);
    }
    else {

    // Memory has been successfully allocated
    printf("Memory successfully allocated using malloc.\n");
    }

	/* alloc memory for dbhs */
	Seedling_d_prt = (float*) malloc(S_dbh * sizeof(float));

	if (Seedling_h_ptr == NULL) {
        printf("Memory not allocated.\n");
        exit(0);
    }
    else {

    // Memory has been successfully allocated
    printf("Memory successfully allocated using malloc.\n");
	}

	/* alloc memory for ages */
	Seedling_a_ptr = (int*) malloc(S_age * sizeof(int));

	if (Seedling_a_ptr == NULL) {
        printf("Memory not allocated.\n");
        exit(0);
    }
    else {

    // Memory has been successfully allocated
    printf("Memory successfully allocated using malloc.\n");
	}

	/* alloc memory for species */
    Seedling_s_ptr = (int*) malloc(S_species * sizeof(int));

	if (Seedling_s_ptr == NULL) {
        printf("Memory not allocated.\n");
        exit(0);
    }
    else {

    // Memory has been successfully allocated
    printf("Memory successfully allocated using malloc.\n");

 }
 return 0;
}

#endif // 0

#if 1     //THIS FUNCTION WORKS

static int fill_cell_for_regeneration (cell_t *const c, species_t *const s) {


	char* p;      //Creates a variable p that point a string
	height_t* h;  //Creates a pointer for height
	dbh_t* d;     //Same for dbh
	age_t* a;     //Same for age


	static height_t height   = { 0 }; //Initializes a variable named height of typedef struct height_t and sets its value to { 0 }
	static dbh_t dbh         = { 0 }; //Same for dbh
	static age_t age         = { 0 }; //Same for age
	static species_t species = { 0 }; //Same for species

	//Is a preprocessor macro (<assert.h>) that checks the truthiness of the expression (c)
	assert(c);

	//Allocation memory for heights using the structure.c module function
	if ( ! alloc_struct((void **)&c->heights, &c->heights_count, &c->heights_avail, sizeof(height_t)) )
    {
        return 0;
     }

	c->heights[c->heights_count-1]       = height;
	c->heights[c->heights_count-1].value = (double)1.3;
	h = &c->heights[c->heights_count-1];

	//Allocation memory for dbhs using the structure.c module function
	if ( ! alloc_struct((void **)&h->dbhs, &h->dbhs_count, &h->dbhs_avail, sizeof(dbh_t)) )
	{
         return 0;
	 }

	h->dbhs[h->dbhs_count-1]       = dbh;
	h->dbhs[h->dbhs_count-1].value = (double)1.0;
	d = &h->dbhs[h->dbhs_count-1];

	//Allocation memory for ages using the structure.c module function
	if ( ! alloc_struct((void **)&d->ages, &d->ages_count, &d->ages_avail, sizeof(age_t)) )
	{
         return 0;
     }

	d->ages[d->ages_count-1]       = age;
	d->ages[d->ages_count-1].value = (int)1;
	a = &d->ages[d->ages_count-1];

	//Allocation memory for species using the structure.c module function
	if ( ! alloc_struct((void **)&a->species, &a->species_count, &a->species_avail, sizeof(species_t)) )
	{
         return 0;
     }

    //Species name is got by the structure
	p = string_copy(s->name);
	if ( ! p )
	{
	 return 0;
	//printf("Name species = %s\n", p);
	 }

	//All of these pointer and vector below refers to species typedef struct
	a->species[a->species_count-1]                  = species;
	a->species[a->species_count-1].management       = T;
	a->species[a->species_count-1].name             = p;
    a->species[a->species_count-1].counter[N_TREE]  = s->counter[SEEDLINGS_SURV]; //The limit seems to be 700.000 seedlings
	a->species[a->species_count-1].counter[N_STUMP] = 0;
	a->species[a->species_count-1].value[LAI_PROJ]  = (double)1.5;

/*
        printf("altezza = %f\n", c->heights[c->heights_count-1].value);
        printf("dbh =      %f\n", h->dbhs[h->dbhs_count-1].value);
        printf("age =      \t%d\n", d->ages[d->ages_count-1].value);
        printf("Number of seedlings =  \t%d\n", a->species[a->species_count-1].counter[N_TREE]);
        printf("LAI =   %f\n", a->species[a->species_count-1].value[LAI_PROJ]);
*/


	return 1;
}

int recruitment (cell_t *const c, species_t *s, const int day, const int month, const int year)
{

    int height;  //Used to assign a position in a vector of heights
	int dbh;     // Same for dbhs
	int age;     // Same for ages
	int species; // Same for species
	int day_temp;  // Used to calculates vegetative days
	int month_temp; // Same above
	int DaysInMonth [] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 }; //Vector of days with its value (days per month)
    int recr = 1;



	height_t *h; //Used to assign pointer to a shortcut height
	dbh_t *d;    //Same for dbh
	age_t *a;    //Same for age
    //species_t *s; //Same for species


        if ( ! fill_cell_for_regeneration ( c, s ) )
	{
       return 0;
	 }


	 height = c->heights_count - 1;   //Assigned the value of the last index of the array
	 dbh = c->heights[height].dbhs_count - 1;
	 age = c->heights[height].dbhs[dbh].ages_count - 1;
	 species = c->heights[height].dbhs[dbh].ages[age].species_count - 1;

	 h = &c->heights[height]; //Here you're accessing the structure nested (typedef) to read the data
     d = &h->dbhs[dbh];
     a = &d->ages[age];
     s = &a->species[species];

       //printf("Name species = %s\n", s->name);
       //printf("altezza = %f\n", h);
       //printf("dbh =      %f\n", h->dbhs[h->dbhs_count-1].value);
       //printf("age =      \t%d\n", d->ages[d->ages_count-1].value);
       //printf("Number of seedlings =  \t%d\n", a->species[a->species_count-1].counter[N_TREE]);
       //printf("LAI =   %f\n", a->species[a->species_count-1].value[LAI_PROJ]);



	/* fill with species values from parameterization file */
	if ( ! fill_species_from_file ( &c->heights[height].dbhs[dbh].ages[age].species[species]) )
	 {
      return 0;
     }

	// check for veg days
	for (month_temp = 0; month_temp < 12; ++month_temp)
	{
		// for handling leap years
		int days_per_month;

		days_per_month = DaysInMonth[month_temp];
		if ( (FEBRUARY == month_temp) && IS_LEAP_YEAR(c->years[year].year) )
		{
			++days_per_month;
		}

		for ( day_temp = 0; day_temp < days_per_month; ++day_temp )
		{
			// compute annually the days for the growing season before any other process
			Veg_Days ( c , day_temp, month_temp, year );
		}

	} //end for

	// initialize new power function
	allometry_power_function           ( c );

	// initialize new carbon pool fraction
	carbon_pool_fraction               ( c );

    //Flag in structure.c
    c->recr = 1;

	// initialize new forest structure
	initialization_forest_structure    (c , day, month, year);
	//annual_forest_structure ( c, year );

	//Need to put zero height for the new layer
	height = 0;

	// initialize new forest class pools
	initialization_forest_class_C      ( c, height, dbh, age, species );

	// initialize new nitrogen pools
	initialization_forest_class_N      ( c, height, dbh, age, species );

	// initialize new litter pools
	initialization_forest_class_litter ( c, height, dbh, age, species );

    //ddalmo august 2021 update forest cell pool
	initialization_forest_cell_C ( c, height, dbh, age, species );

	// print new forest class dataset
	print_new_daily_forest_class_data  ( c, height, dbh, age, species );


	return 1;
}

#endif // 1

#if 0
                 //*********** REAL SEEDLINGS TO SAPLINGS RECRUITMENT FUNCTION *********************//

static int recruitment_memory (cell_t *const c)
{

    char* p;
	height_t* h;
	dbh_t* d;
	age_t* a;
	species_t* s;

	static height_t height   = { 0 };
	static dbh_t dbh         = { 0 };
	static age_t age         = { 0 };
	static species_t species = { 0 };

    assert(c);

    // Alloc memory for number of Seedlings survived
    s->counter[SEEDLINGS_TO_SAPLINGS] = malloc(sizeof(s->counter[SEEDLINGS_TO_SAPLINGS]));

	if (s->counter[SEEDLINGS_TO_SAPLINGS] == NULL) {
        printf("Memory not allocated.\n");
        exit(0);
    } else {

    // Memory has been successfully allocated
    printf("Memory successfully allocated using malloc.\n");
    }

	// Alloc memory for heights
	height_t* height = malloc(sizeof(height_t));

	if (height == NULL) {
        printf("Memory not allocated.\n");
        exit(0);
    } else {

    // Memory has been successfully allocated
    printf("Memory successfully allocated using malloc.\n");
    }

	// Alloc memory for dbhs
	dbh_t* dbh = malloc(sizeof(dbh_t));

	if (dbh == NULL) {
        printf("Memory not allocated.\n");
        exit(0);
    } else {

    // Memory has been successfully allocated
    printf("Memory successfully allocated using malloc.\n");
	}

	// Alloc memory for ages
	age_t* age = malloc(sizeof(age_t));


	if (age == NULL) {
        printf("Memory not allocated.\n");
        exit(0);
    } else {

    // Memory has been successfully allocated
    printf("Memory successfully allocated using malloc.\n");
	}

	// Alloc memory for species
    species_t* species = malloc(sizeof(species_t));

	if (species == NULL) {
        printf("Memory not allocated.\n");
        exit(0);
    } else {

    // Memory has been successfully allocated
    printf("Memory successfully allocated using malloc.\n");
    }

/*
    // alloc memory for species_name
    name_s = (char*) malloc(species_t * sizeof(char));

	if (name_s == NULL) {
        printf("Memory not allocated.\n");
        exit(0);
    } else {

    // Memory has been successfully allocated
    printf("Memory successfully allocated using malloc.\n");
 } free(name_s);
*/
 return 0;
}
#endif // 0

#if 0

 int recruitment (cell_t *const c, const int day, const int month, const int year) {

    int Seedling_surv;
    float height = 10.; //cm
	float dbh = 1.; //cm
	int age = 2;
	int species;
	int day_temp;
	int month_temp;
	int DaysInMonth [] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };


 // Check for veg days
	for (month_temp = 0; month_temp < 12; ++month_temp)
	{
		// for handling leap years
		int days_per_month;

		days_per_month = DaysInMonth[month_temp];
		if ( (FEBRUARY == month_temp) && IS_LEAP_YEAR(c->years[year].year) )
		{
			++days_per_month;
		}

		for ( day_temp = 0; day_temp < days_per_month; ++day_temp )
		{
			// compute annually the days for the growing season before any other process
			Veg_Days ( c , day_temp, month_temp, year );
		}
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
	//print_new_daily_forest_class_data  ( c, c->heights_count-1, 0, 0, 0 );

	return 1;
}

#endif // 1

#if 0   //TESTARE
static int recruitment_memory (cell_t *const c)
{

   height_t* allocate_height() {
    height_t* height = malloc(sizeof(height_t));
    if (height == NULL) {
        // Errore di allocazione della memoria
         return 0;
    }
    return height;
 }

   dbh_t* allocate_dbh() {
    dbh_t* dbh = malloc(sizeof(dbh_t));
    if (dbh == NULL) {
        // Errore di allocazione della memoria
         return 0;
    }
    return dbh;
 }

   age_t* allocate_age() {
    age_t* age = malloc(sizeof(age_t));
    if (age == NULL) {
        // Errore di allocazione della memoria
         return 0;
    }
    return age;
 }

   species_t* allocate_species() {
    species_t* species = malloc(sizeof(species_t));
    if (species == NULL) {
        // Errore di allocazione della memoria
         return 0;
    }
    return species;
 }
}


int recruitment (cell_t *const c, const int day, const int month, const int year)
{

	int height = 1.2;
	int dbh = 1.;
	int age = 1;
	int species;
/*	int day_temp;
	int month_temp;
	int DaysInMonth [] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
*/

	if (! recruitment_memory (c)) return 0;

	height = c->heights_count - 1;
	dbh = c->heights[height].dbhs_count - 1;
	age = c->heights[height].dbhs[dbh].ages_count - 1;
	species = c->heights[height].dbhs[dbh].ages[age].species_count - 1;

	if ( ! fill_species_from_file(&c->heights[c->heights_count-1].dbhs[0].ages[0].species[0]) )
	{
		return 0;
	}
/*
	// check for veg days
	for (month_temp = 0; month_temp < 12; ++month_temp)
	{
		// for handling leap years
		int days_per_month;

		days_per_month = DaysInMonth[month_temp];
		if ( (FEBRUARY == month_temp) && IS_LEAP_YEAR(c->years[year].year) )
		{
			++days_per_month;
		}

		for ( day_temp = 0; day_temp < days_per_month; ++day_temp )
		{
			// compute annually the days for the growing season before any other process
			Veg_Days ( c , day_temp, month_temp, year );
		}
	}

/*
	// initialize new power function
	allometry_power_function           ( c );

	// initialize new carbon pool fraction
	carbon_pool_fraction               ( c );

	// initialize new forest structure
	initialization_forest_structure    (c , day, month, year);

    // comment: in the initialization_forest_structure the new added layer has height index = 0 and the
    // dominant layers index >= 1 ordered according to descending height

    height = 0;  // so to initialize the new layer only! Which is the new added 'regeneration-layer'

	// initialize new forest class pools
	initialization_forest_class_C      ( c, height, dbh, age, species );

	// initialize new nitrogen pools
	initialization_forest_class_N      ( c, height, dbh, age, species );

	// initialize new litter pools
	initialization_forest_class_litter ( c, height, dbh, age, species );

        //ddalmo august 2021
        // update forest cell pool
	initialization_forest_cell_C ( c, height, dbh, age, species );

	// print new forest class dataset
	print_new_daily_forest_class_data  ( c, height, dbh, age, species );
*/
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
#endif // 1
