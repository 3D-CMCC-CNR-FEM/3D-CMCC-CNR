/*
 * canopy_cover.c
 *
 *  Created on: 02/mag/2017
 *      Author: alessio-cmcc
 */

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

extern logger_t* g_debug_log;
extern settings_t* g_settings;


void dbhdc_function ( cell_t *const c, const int layer, const int height, const int dbh, const int age, const int species, const int year)
{
	double temp_crown_area     = 0.;
	double temp_crown_radius   = 0.;
	double temp_crown_diameter = 0.;
	double previous_dbhdc_eff  = 0.;
	double max_dbhdc_incr      = 0.1;      /* fraction of maximum dbhdc increment */ //note: this was 0.001 in v.5.4
	//double max_dbhdc_decr      = 0.001;  /* fraction of maximum dbhdc decrement */


	dbh_t *d;
	species_t *s;
	d = &c->heights[height].dbhs[dbh];
	s = &c->heights[height].dbhs[dbh].ages[age].species[species];

	logger(g_debug_log,"\n*DBHDC FUNCTION for %s for %d*\n", s->name, year);

	/************************************************************************************************************************/
	/* note: 04 Oct 2016 */
	/* compute potential maximum and minimum density for DBHDC function */
	/* new DBHDC function */
	/* this function in mainly based on the assumptions that trees tend to occupy */
	/* all space they can, if they cannot then fixed values constrain their crown */
	/* see also: Lhotka and Loewenstein, 1997; Lhotka and Loewenstein, 2008 */

	/* previous dbhdc eff */
	previous_dbhdc_eff  = s->value[DBHDC_EFF];
	logger(g_debug_log,"-DBHDC (old)         = %f\n", s->value[DBHDC_EFF]);


	temp_crown_area     = ( s->value[MAX_LAYER_COVER] * g_settings->sizeCell ) / ( c->tree_layers[layer].layer_density * g_settings->sizeCell );
	logger(g_debug_log,"-temp_crown_area     = %f\n", temp_crown_area);



	temp_crown_radius   = sqrt(temp_crown_area / Pi);
	logger(g_debug_log,"-temp_crown_radius   = %f\n", temp_crown_radius);

	temp_crown_diameter = temp_crown_radius * 2.;
	logger(g_debug_log,"-temp_crown_diameter = %f\n", temp_crown_diameter);

	s->value[DBHDC_EFF] = temp_crown_diameter / d->value;
	logger(g_debug_log,"-DBHDC (new)         = %f\n", s->value[DBHDC_EFF]);


	/* check if current dbhdc_eff grows too much (case when there's thinning) */
	/* this is checked to avoid unrealistic crown area increment */

	/* note: max_dbhdc_incr corresponds to an arbitrary increment of n value */
	/* note: not used in the first year of simulation */
	if ( ( s->counter[YOS] ) && ( s->value[DBHDC_EFF] > ( previous_dbhdc_eff + ( previous_dbhdc_eff * max_dbhdc_incr ) ) ) )
	{
		s->value[DBHDC_EFF] = previous_dbhdc_eff + ( previous_dbhdc_eff * max_dbhdc_incr );
	}


	/***************************************************************************************************/
	//note test: 18 June 2018
	//data obtained by:	Ritter and Nothdurft et al., 2018 forests
	//dbhdcmax decreases as dbh increases

#if 0

	if ( s->value[PHENOLOGY] == 0.1 || s->value[PHENOLOGY] == 0.2 )
	{
		s->value[DBHDCMAX] = 0.9667 * pow ( d->value , -0.287 );
	}
	else
	{
		s->value[DBHDCMAX] = 0.8543 * pow ( d->value , -0.254 );
	}

#endif


       // 5p606 
       // calibrate equation on literature data
       // the equation defines how fast the canopy can horizontally develop when DBH increase (it is suppose than when the
       // the density is low), and it defines in turn the chance
       // for the regeneration layer to growh

#if 1

	if ( s->value[PHENOLOGY] == 0.1 || s->value[PHENOLOGY] == 0.2  )
	{
	
            //  s->value[DBHDCMAX] = 0.5639 * pow ( d->value , -0.285 );   //  5p605
              
                s->value[DBHDCMAX] = 0.4991 * pow ( d->value , -0.202 );   // Castanea Sativa (upper limit) 5p606
                 
            //     if (d->value < 10.)
	    //    {
	    //	s->value[DBHDCMAX] = 0.4991 * pow ( 10. ,  -0.202);  // 5p607 
            //}
	
	    if (s->value[DBHDCMAX] < s->value[DBHDCMIN] )   s->value[DBHDCMAX] = s->value[DBHDCMIN] ;
        }
	else
	{

	 s->value[DBHDCMAX] = 0.55 * pow ( d->value , -0.309 );  // Ritter's equation for conifers 5p606

	    if (s->value[DBHDCMAX] < s->value[DBHDCMIN] )   s->value[DBHDCMAX] = s->value[DBHDCMIN] ;

	}

#endif
	/**************************************************************************************************/

	/* check */
	if ( s->value[DBHDC_EFF] > s->value[DBHDCMAX] )
	{
		logger(g_debug_log,"-DBHDC effective (%f) > DBHDCMAX (%f) \n", s->value[DBHDC_EFF] , s->value[DBHDCMAX]);
		s->value[DBHDC_EFF] = s->value[DBHDCMAX];

	}


	/************************************************************************************************************************/

	/* if dbhdc decreases than with the same proportion also branch and coarse root fractions decrease */

	if ( s->counter[YOS] && ( previous_dbhdc_eff > s->value[DBHDC_EFF] ) )
	{
		/******** self pruning ********/

		self_pruning ( c, height, dbh, age, species, previous_dbhdc_eff, s->value[DBHDC_EFF] );

	}

	logger(g_debug_log,"-DBHDC effective     = %f\n", s->value[DBHDC_EFF]);

	/* check */
	CHECK_CONDITION (s->value[DBHDC_EFF], < , ZERO);

}

void canopy_cover (cell_t *const c, const int height, const int dbh, const int age, const int species)
{

	species_t *s;
	s = &c->heights[height].dbhs[dbh].ages[age].species[species];

	/*** Canopy cover ***/

	logger(g_debug_log,"\n*CANOPY COVER for %s *\n", s->name);

	/* Canopy and Soil Projected Cover using DBH-DC (at zenith angle)  it is never higher than max_cover_proj */
	s->value[CANOPY_COVER_PROJ] = s->value[CROWN_AREA_PROJ] * s->counter[N_TREE] / g_settings->sizeCell;
	logger(g_debug_log, "-Canopy Projected Cover   = %f %%\n", s->value[CANOPY_COVER_PROJ] * 100.);

	/* check */
	CHECK_CONDITION( s->value[CANOPY_COVER_PROJ] ,  > , s->value[CANOPY_COVER_PROJ] + eps );

	/* Canopy cover able to absorb light (integrated all over all viewing angles) */
	//s->value[CANOPY_COVER_EXP] = (s->value[CROWN_AREA_EXP] * s->counter[N_TREE]) / g_settings->sizeCell ;
	//logger(g_debug_log, "-Canopy Exposed Cover     = %f %%\n", s->value[CANOPY_COVER_EXP]  * 100.);

	//printf("CANOPY_COVER_EXP = %f\n", s->value[CANOPY_COVER_EXP]);

	//if (s->value[CANOPY_COVER_EXP] > 1)s->value[CANOPY_COVER_EXP] = 1.;


}
