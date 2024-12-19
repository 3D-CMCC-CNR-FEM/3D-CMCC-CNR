/*
 * regeneration.c
 *
 *  Created on: 14/nov/2012
 *      Author: alessio
 */

/* includes */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include "matrix.h"
#include "constants.h"
#include "logger.h"
#include "common.h"
#include "settings.h"
#include "regeneration.h"
#include "new_forest_tree_class.h"

extern logger_t* g_debug_log;
extern settings_t* g_settings;

#if 0

void regeneration (cell_t *const c, const int height, const int dbh, const int age, const int species)
{
	double fruit_gDM;
	int seeds_number;
	int saplings_number;

	species_t *s;

        // TODO UNDER DEVELOPMENT 
	// ALESSIOR
	//if ( ! c->heights[height].dbhs[dbh].ages ) return;
	s = &c->heights[height].dbhs[dbh].ages[age].species[species];

	logger(g_debug_log, "\n**REGENERATION**\n");

	/* compute number of seeds */
	/*convert fruit pool from tC to gDM */
	fruit_gDM = s->value[FRUIT_C] * 1000000 * GC_GDM;
	logger(g_debug_log, "fruit_gDM = %g gDM\n", fruit_gDM);

	seeds_number = (int)(fruit_gDM / s->value[WEIGHTSEED]);
	logger(g_debug_log, "fruit biomass= %g tC\n", s->value[FRUIT_C]);
	logger(g_debug_log, "number of seeds = %d\n", seeds_number);

	/* reset annually fruit pool */
	s->value[FRUIT_C] = 0.;

	/* compute number of saplings based on germination capacity */
	saplings_number = (int)(seeds_number * s->value[GERMCAPACITY]);
	logger(g_debug_log, "number of saplings = %d\n", saplings_number);

	logger(g_debug_log, "name = %s\n", s->name);

	/* assign values */
	//FIXME currently prescribed, 
	g_settings->regeneration_n_tree = saplings_number;

	/* it gets name of species that produces seeds */
	//FIXME
	strncpy(g_settings->regeneration_species, ((const char*)s->name), SETTINGS_REGENERATION_SPECIES_MAX_SIZE-1);

	/* replanting tree class */
	if( g_settings->regeneration_n_tree )
	{
		if ( ! add_tree_class_for_regeneration( c ) )
		{
			logger_error(g_debug_log, "unable to add new regeneration class! (exit)\n");
			exit(1);
		}
	}
}

#endif // 0

// SAPONARO 12/2022
        /*******************************************REGENERATION MODULE*******************************************************/

#if 1                                 /*******************USE THIS********************/

int germination (cell_t *const c, const meteo_daily_t *const meteo_daily, species_t *const s, const int day, const int month, const int year) {


    int N_seed_vitality = 0;
    long int Seedlings_Number = 0;



            //Perform number of seeds vitality
            N_seed_vitality = s->counter[N_SEED] * EMPTY_SEEDS;
            //printf("Seed vitality = %d\n", N_seed_vitality);

            //Perform number of seedlings
            Seedlings_Number = N_seed_vitality * s->value[GERMCAPACITY];

            //Assigned seedling number to a species vector
            s->counter[SEEDLINGS] = Seedlings_Number;
            //printf("Seedlings = \t%d\n", s->counter[SEEDLINGS]);


    return 0;
}

#endif // 0

#if 0

int germination (cell_t *const c, const meteo_daily_t *const meteo_daily, species_t *const s, const int day, const int month, const int year) {


 int Seedlings_Number = 0.;          //Number of seedlings per species
 int Seedlings_pool = 0.;            //Seedlings pool

 //Check whether the seedlings are zero in the array
// s->counter[SEEDLINGS] = Seedlings_Number;


  meteo_t *met;
  met = c->years[year].m;

   // Winter condition (Jan-March)
   if( c->doy >= 1 && c->doy <= 120) {

     float SumTerm = meteo_daily->thermic_sum;
     //printf("Thermicsum = %f\n", SumTerm);
     float Soil_T = c->years[year].yearly_mean.tsoil;
     //printf("Soiltemp = %f\n", Soil_T);
     float swp = c->psi;
     //printf("swp: %f\n", swp);


	    // FIXME replace asw with Soil Water Potential (see E.Falleri et al.,2004)
        if (SumTerm >= s->value[GDD_SEED] && Soil_T > 5. && swp >= -0.2) {

                //Germination capacity is the average between 0 and -0.2 MPa (see table 2 in E.Falleri et al.,2004)
    			Seedlings_Number = (s->counter[TANK_SEEDS] * 0.62);
                //Empty Tank seeds when seedlings appear
                s->counter[TANK_SEEDS] = 0.;
    			//Assigned seedling number to a species vector
                s->counter[SEEDLINGS] = Seedlings_Number;
                //Seedlings appeared over years going into the pool
                s->counter[SEEDLINGS_POOL] += s->counter[SEEDLINGS];
    			}

        else if (SumTerm >= s->value[GDD_SEED] && Soil_T > 5. && swp < -0.2 && swp >= -0.6) {
    			Seedlings_Number = (s->counter[TANK_SEEDS] * 0.48);

    			s->counter[TANK_SEEDS] = 0.;

                s->counter[SEEDLINGS] = Seedlings_Number;

                s->counter[SEEDLINGS_POOL] += s->counter[SEEDLINGS];
    			}

        else if (SumTerm >= s->value[GDD_SEED] && Soil_T > 5. && swp < -0.6 && swp >= -1.2) {
    			Seedlings_Number = (s->counter[TANK_SEEDS] * 0.20);

                s->counter[TANK_SEEDS] = 0.;

                s->counter[SEEDLINGS] = Seedlings_Number;

                s->counter[SEEDLINGS_POOL] += s->counter[SEEDLINGS];
                }

         else {
                Seedlings_Number = 0.;
                printf("Germination failed due to environmental condition.\n");
              }
            }

  if ( c->doy == ( IS_LEAP_YEAR( c->years[year].year ) ? 366 : 365 )) {

  printf("Seeds = %ld\n", s->counter[N_SEED]);
  printf("Tankseeds = %ld\n", s->counter[TANK_SEEDS]);
  printf("Seedlings = %d\n", s->counter[SEEDLINGS]);
  printf("Seedlings_pool = %d\n", s->counter[SEEDLINGS_POOL]);
  }
 return 0;
}

#endif // 0



#if 0 // ******* Test those hypotesis FIXME, DO NOT USE ******** //

int germination (cell_t *const c, const meteo_daily_t *const meteo_daily, species_t *const s, const int year) {


  int Seedlings_Number = 0.;          //Number of seedlings per species
  int Seedlings_pool = 0.;            //Pool of seedlings

  meteo_t *met;

  //Variable inizialized in matrix.h
  Seedlings_Number = s->counter[SEEDLINGS];
  Seedlings_pool = s->counter[SEEDLINGS_POOL];

  //Start the seedling germination after the first year of simulation
  if ( !years && ( c->doy == ( IS_LEAP_YEAR( c->years[year].year ) ? 366 : 365))) {

    //Temperature, light and water condition to perform the seeds germination
    if ( meteo_daily->thermic_sum >= s->value[GDD_SEED] &&  meteo_daily->sw_downward_W >= 0. && c->asw >= 0. && c->psi < min_psi) {

      //The model gets the seeds produced the year before
      if (*yos_count > 1) {

        yos[*yos_count-1].m[month].d[day].s->counter[N_SEED] = yos[*yos_count-2].m[month].d[day].s->counter[N_SEED];

         //Perform seeds germination based on germcapacity
         Seedlings_Number = (s->counter[TANK_SEEDS] * s->value[GERMCAPACITY]);

         //Assigned seedling number to a species vector
         s->counter[SEEDLINGS] = Seedlings_Number;
           }
         }
         else {

         s->counter[SEEDLINGS] = 0;

         }
         //Condition if in year of simulation the Seedlingpool is filled hence the seeds do not accumulate again
         do {

         //Accumulation of seedling pool
         s->counter[SEEDLINGS_POOL] += s->counter[SEEDLINGS];

         } while (s->counter[TANK_SEEDS] = 0);

    }

 return 0;
}

#endif

// ********* FIXME, DO NOT USE******** //
#if 0
int establishment(cell_t *const c, const meteo_daily_t *const meteo_daily, species_t *const s, const int day, const int month, const int year) {

 int Saplings_Number = 0.;            //Number of saplings per species per year
 int Saplings_pool = 0.;              //Saplings' pool per year


 meteo_t *met;
 met = c->years[year].m;


//To be more precisely put the months (in this case is summer).
 if (c->doy >= 150 && c->doy <= 245) {

 //Condition in which seedlings can grow to become saplings.
 //Note: Radiation that reach the soil surface is in mol/m2*day, in Muffler L. et al., 2021 is in (80)micromol/m2*sec.
 //It is important to convert the Muffler's number in mol/m2*day as in the model calculation (in this case is 7 mol/m2*day).
  if (meteo_daily->par >= 7 && meteo_daily->rh_f >= 65) {

   //Number of saplings estabilished
   Saplings_Number = s->counter[SEEDLINGS_POOL] * 0.70; //Value observed in MUffler L. er al., 2021

   s->counter[SAPLINGS] = Saplings_Number;

        do {

   //Accumulation of saplings over year
   s->counter[SAPLINGS_POOL] += s->counter[SAPLINGS];
   //printf("Saplings = %d\n", s->counter[SAPLINGS]);

       } while (s->counter[SEEDLINGS_POOL] = 0.);

     } else {

    s->counter[SAPLINGS] = 0;

   }
     // printf("parsoil = %f\n", meteo_daily->par);
     // printf("doy = %d\n", c->doy);
     // printf("rh = %f\n", meteo_daily->rh_f);


  } if ( c->doy == ( IS_LEAP_YEAR( c->years[year].year ) ? 366 : 365 )) {

     // printf("Seeds = %d\n", s->counter[N_SEED]);
     // printf("Tankseeds = %d\n", s->counter[TANK_SEEDS]);
     // printf("Seedlings = %d\n", s->counter[SEEDLINGS]);
     // printf("Seedlings_pool = %d\n", s->counter[SEEDLINGS_POOL]);
     // printf("Saplings = %d\n", s->counter[SAPLINGS]);
     // printf("Saplings_pool = %d\n", s->counter[SAPLINGS_POOL]);

  }

  //printf("parsoil = %f\n", meteo_daily->par);
  //printf("doy = %d\n", c->doy);
     //printf("rh = %f\n", meteo_daily->rh_f);


  return 0;
}


#endif // 0

#if 1        /**********************************USE THIS****************************/

int establishment (cell_t *const c, const meteo_daily_t *const meteo_daily, species_t *const s, const int day, const int month, const int year) {


    long int Seedlings_surv = 0; //Seedlings number that survive after first year

       //Note: Radiation that reachs the soil surface in Muffler L. et al., 2021 is in (max surv)80 micromol/m2*sec.
       //The model calculates the radiation in mol/m2*day (in this case is 7 mol/m2*day).

        //Optimal condition for light and temperature (case1)
		if (meteo_daily->seedling_par >= s->value[SURV_PAR] && meteo_daily->seedling_temp <= s->value[SURV_TEMP]) {

			Seedlings_surv = (s->counter[SEEDLINGS] * 0.01); //L.Muffler et al.2021 0.70

            s->counter[SEEDLINGS_SURV] = Seedlings_surv;

             }

        //Optimal condition for light but not for temperature (case2)
        else if (meteo_daily->seedling_par >= s->value[SURV_PAR] && meteo_daily->seedling_temp > s->value[SURV_TEMP]) {

	        Seedlings_surv = (s->counter[SEEDLINGS] * 0.01); //L.Muffler et al.2021 0.55
            s->counter[SEEDLINGS_SURV] = Seedlings_surv;

             }

        //Optimal condition for temperature but not for light (case3)
        else if (meteo_daily->seedling_par < s->value[SURV_PAR] && meteo_daily->seedling_temp < s->value[SURV_TEMP]) {

			Seedlings_surv = (s->counter[SEEDLINGS] * 0.01); //L.Muffler et al.2021 0.30
            s->counter[SEEDLINGS_SURV] = Seedlings_surv;

             }

		//Non-Optimal condition for Light and Temperature (case 4)
        else if (meteo_daily->seedling_par < s->value[SURV_PAR] && meteo_daily->seedling_temp > s->value[SURV_TEMP]) {

			Seedlings_surv = (s->counter[SEEDLINGS] * 0.01); //L.Muffler et al.2021 0.20
            s->counter[SEEDLINGS_SURV] = Seedlings_surv;


       } else {

	    Seedlings_surv = 0;
	    //printf("All seedlings died.\n");

    }

    //printf("Rad = %f\n", meteo_daily->par);
    //printf("Temperature = \t%f\n", meteo_daily->tavg);
/*
    if ( c->doy == ( IS_LEAP_YEAR ( c->years[year].year ) ? 366 : 365) ) {

      //Accumulate seedlings until become saplings
      s->counter[SEEDLINGS_POOL] += s->counter[SEEDLINGS_SURV];

    }
*/
    //printf("Soil par seedlings = \t%f\n", meteo_daily->seedling_par);
    //printf(" Temp seedlings = \t%f\n", meteo_daily->seedling_temp);



 return 0;
}
#endif
