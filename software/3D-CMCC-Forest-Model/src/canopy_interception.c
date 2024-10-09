/*
 * canopy_interception.c
 *
 *  Created on: 19/ago/2016
 *      Author: alessio-cmcc
 */


/* includes */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "common.h"
#include "canopy_evapotranspiration.h"
#include "constants.h"
#include "logger.h"
#include "Penman_Monteith.h"
#include "heat_fluxes.h"


void canopy_interception(cell_t *const c, const int layer, const int height, const int dbh, const int age, const int species, meteo_daily_t *const meteo_daily)
{

	tree_layer_t *l;
	species_t *s;

	l = &c->tree_layers[layer];
	s = &c->heights[height].dbhs[dbh].ages[age].species[species];

	//printf(" c->canopy_int_cell_height_class_counter %d\n",c->canopy_int_cell_height_class_counter);

	/* it computes canopy rain interception with a slightly different rationale used for canopy radiation absorption */
	/* note: it is still used cause rainfall is considered always vertical from zenith angle */

	/* CANOPY INTERCEPTION */

	/* first height class in the cell is processed */
	if( ! l->canopy_int_layer_height_class_counter && ! c->canopy_int_cell_height_class_counter )
	{
		/* reset temporary values when the first height class in layer is processed */
		c->temp_int_rain = 0.;
		c->temp_int_snow = 0.;

		/* assign meteo variables to cell variables */
		/* assign incoming rain */
		//c->rain = meteo_daily->rain;
		/* assign incoming snow */
		//c->snow = meteo_daily->snow;
	}

	/*****************************************************************************************************************/
	/* counters */
	l->canopy_int_layer_height_class_counter ++;
	c->canopy_int_cell_height_class_counter++;

	/*************************************************************************/
	/* shared functions among all class/layers */

	//note: not accumulation of rain or snow happens if canopy is wet or snowed */
	if( ( meteo_daily->prcp > 0.) && ( s->value[LAI_PROJ] > 0. ) && ( ! s->value[CANOPY_WATER] ) && ( ! s->value[CANOPY_SNOW] ) )
	{
		/* for rain */
		if( meteo_daily->rain )
		{
			/* compute maximum water storage (interception) (mm/m2 area covered) */

			double Int_max_rain;    /* maximum intercepted rain (mm/m2) */
#if 0
			/* (1) compute rain interception for dry canopy (Lawrence et al., 2006) */
			s->value[CANOPY_INT]      = s->value[INT_COEFF] * meteo_daily->rain * ( 1. - exp(-0.5 * s->value[LAI_PROJ])) * s->value[DAILY_CANOPY_COVER_PROJ];

			/* (2) compute rain interception for dry canopy (Tatarinov et al., 2007) */
			Int_max_rain              = 0.3 * s->value[ALL_LAI_PROJ];
			s->value[CANOPY_INT_RAIN] = MIN( Int_max_rain , meteo_daily->rain );

#else
			//june 2017
			/* following Jiao et al., 2016, Water Eq. [8] pg. 9 */
			{
				Int_max_rain = 0.284 + 0.092 * s->value[LAI_PROJ] * ( 1. - exp ( -0.231 * meteo_daily->rain ) );
				s->value[CANOPY_INT_RAIN] = MIN( Int_max_rain , meteo_daily->rain );
			}

#endif

			// 5p606
			// check condition when we have several classes within the same layer and the total canopy cover projection is higher than 1
			
			if ( ( s->value[CANOPY_INT_RAIN] >  meteo_daily->rain )  &&  (c->tree_layers[layer].daily_layer_cover_proj > 1)) ;
			 {
			 
			  s->value[CANOPY_INT_RAIN] =0.;   // we set to 0 the 'lower' (as z) class in the layer).
			 
			 }
                         

			/* update pool */
			s->value[CANOPY_WATER]    = s->value[CANOPY_INT_RAIN];

			/* make zero snow */
			s->value[CANOPY_INT_SNOW] = 0.;

			CHECK_CONDITION( s->value[CANOPY_INT_RAIN], > , meteo_daily->rain );
		}
		/* for snow */
		else
		{
			/* following Dewire (PhD thesis) and Pomeroy et al., 1998., Hedstrom & Pomeroy, 1998 */

			double Int_max_snow;      /* maximum intercepted snow (mm/m2)*/

			Int_max_snow              = 4.4 * s->value[LAI_PROJ];

			//test
#if 0
			s->value[CANOPY_INT_SNOW] = 0.;
#else
			//s->value[CANOPY_INT_SNOW] = s->value[CANOPY_SNOW] + 0.7 * ( Int_max_snow - s->value[CANOPY_SNOW] ) * (1. - exp( - ( meteo_daily->snow / Int_max_snow ) ) ) * s->value[DAILY_CANOPY_COVER_PROJ];
			
			//5p6 we multiply later for the cc_proj
			s->value[CANOPY_INT_SNOW] = s->value[CANOPY_SNOW] + 0.7 * ( Int_max_snow - s->value[CANOPY_SNOW] ) * (1. - exp( - ( meteo_daily->snow / Int_max_snow ) ) ) ;
#endif

 
 
                         //5p606
			// check condition when we have several classes within the same layer and the total canopy cover projection is higher than 1
			
			if ( ( s->value[CANOPY_INT_SNOW] >  meteo_daily->snow )  &&  (c->tree_layers[layer].daily_layer_cover_proj > 1)) ;
			 {
			 
			  s->value[CANOPY_INT_SNOW] =0.;   // we set to 0 the 'lower' (as z) class in the layer).
			 
			 }
                        
                         

			/* update pool */
			//fixme for now assuming no snow accumulation
			s->value[CANOPY_SNOW]     = s->value[CANOPY_INT_SNOW];

			/* make zero rain */
			s->value[CANOPY_INT_RAIN] = 0.;

			CHECK_CONDITION(s->value[CANOPY_INT_SNOW], > , meteo_daily->snow);
		}
		 // 5p6 convert stored canopy snow and water as value per m2 of grid cell
                 // note that canopy water and snow can not accumulate
                 s->value[CANOPY_SNOW]  = s->value[CANOPY_SNOW]  * s->value[DAILY_CANOPY_COVER_PROJ];
                 s->value[CANOPY_WATER] = s->value[CANOPY_WATER] * s->value[DAILY_CANOPY_COVER_PROJ];
	}

	/**********************************************************************************************************/

	/* update temporary rain and snow */

	// 5p6: please note that CANOPY_INT_RAIN/SNOW refer to m2 of canopy, hence, in order to compute the average amount of rain
	// reaching the lower layer, it should account for the CC_PROJ
	
	// TODO the overlapping of classes within the same layer!
	// this is also important in order to compute ASW at cell level.

	/** TODO **/ //check and correct if whithin the same laver overlapping of classe occurrs

	c->temp_int_rain += s->value[CANOPY_INT_RAIN]* s->value[DAILY_CANOPY_COVER_PROJ];
	c->temp_int_snow += s->value[CANOPY_INT_SNOW]* s->value[DAILY_CANOPY_COVER_PROJ];

	CHECK_CONDITION( c->temp_int_rain, > , meteo_daily->rain );  // controll on the value at cell/layer-level
        CHECK_CONDITION( c->temp_int_snow, > , meteo_daily->snow );


	//c->temp_int_rain += s->value[CANOPY_INT_RAIN];
	//c->temp_int_snow += s->value[CANOPY_INT_SNOW];

	/**********************************************************************************************************/

	/* when matches the last height class in the layer is processed */
	if ( l->layer_n_height_class == l->canopy_int_layer_height_class_counter )
	{
		/* compute interceptable rain for lower layers */
		//c->daily_canopy_rain_int += s->value[CANOPY_INT_RAIN];
		c->daily_canopy_rain_int += s->value[CANOPY_INT_RAIN]* s->value[DAILY_CANOPY_COVER_PROJ];
		meteo_daily->rain        -= c->temp_int_rain;

                CHECK_CONDITION( meteo_daily->rain , <, 0. );

		/* compute interceptable snow for lower layers */
		//c->daily_canopy_snow_int += s->value[CANOPY_INT_SNOW];
		c->daily_canopy_snow_int += s->value[CANOPY_INT_SNOW]* s->value[DAILY_CANOPY_COVER_PROJ];
		meteo_daily->snow        -= c->temp_int_snow;

		/* reset temporary values when the last height class in layer is processed */
		c->temp_int_rain = 0.;
		c->temp_int_snow = 0.;

		/* reset counter */
		l->canopy_int_layer_height_class_counter = 0;
	}

	/**********************************************************************************************************/
	/* when it matches the last height class in the cell is processed */
	if ( c->heights_count == c->canopy_int_cell_height_class_counter )
	{
		/* last height class in cell processed */

		/* reset counter */
		c->canopy_int_cell_height_class_counter = 0;
	}
	/*****************************************************************************************************************/

	/* cumulate */
	s->value[MONTHLY_CANOPY_INT] += (s->value[CANOPY_INT_RAIN] + s->value[CANOPY_INT_SNOW]);
	s->value[YEARLY_CANOPY_INT]  += (s->value[CANOPY_INT_RAIN] + s->value[CANOPY_INT_SNOW]);

	// 5p6 convert intercepted canopy snow and water as value per m2 of grid cell
        s->value[CANOPY_INT_SNOW] = s->value[CANOPY_INT_SNOW] * s->value[DAILY_CANOPY_COVER_PROJ];
        s->value[CANOPY_INT_RAIN] = s->value[CANOPY_INT_RAIN] * s->value[DAILY_CANOPY_COVER_PROJ];
}
