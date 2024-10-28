/*
 * canopy_radiation_sw_band.c
 *
 *  Created on: 05 lug 2016
 *      Author: alessio
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "common.h"
#include "constants.h"
#include "logger.h"
#include "soil_settings.h"
#include "topo.h"
#include "canopy_radiation_sw_band.h"

extern logger_t* g_debug_log;

#define TEST 0

void canopy_sw_band_abs_trans_refl_radiation(cell_t *const c, const int layer,const int height, const int dbh, const int age, const int species, const meteo_daily_t *const meteo_daily,
		double Light_abs_frac, double Light_abs_frac_sun, double Light_abs_frac_shade, double Light_refl_par_frac, double Light_refl_par_frac_sun, double Light_refl_par_frac_shade,
		double Light_refl_sw_frac, double Light_refl_sw_frac_sun, double Light_refl_sw_frac_shade)
{
	species_t *s;
	s = &c->heights[height].dbhs[dbh].ages[age].species[species];

	// note: This function works at layer/class level computing absorbed, transmitted and reflected PAR, Short Wave
	// and PPFD through different height classes/layers considering at square meter takes into account coverage,
	// it means that a square meter grid cell * represents overall grid cell (see Duursma and Makela, 2007) 
	// The scheme is similar to CLM(ED) Norman scheme.
	// the radiation onto the understory layer is a spatial (canopy-projecton) weighted value 

	/* it follows a little bit different rationale compared to BIOME-BGC approach
	 * in BIOME-BGC:
	 * apar = par * (1 - (exp(- K * LAI)));
	 * apar_sun = par * (1 - (exp(- K * LAI_SUN)));
	 * apar_shade = apar- apar_sun;
	 *
	 * in 3D-CMCC FEM:
	 * apar_sun = par * (1 - (exp(- K * LAI_SUN)));
	 * par_transm_sun  = par - apar_sun;
	 * apar_shade = par_transm_sun * (1 - (exp(- K * LAI_SHADE)));
	 * apar = apar_sun + apar_shade;
	 *
	 * then it consider that an amount of sunlit leaf are not completely outside the canopy
	 * but there's an exponential decay of absorption also for sunlit leaves	 *
	 */

	/***********************************************************************************************/
    double scal_cover = 1. ;
double temp_scale = 1. ;
	// when multispecies in the same layer overlapp, they might give a higher absorpion than what is entering
    // this might happen when some species are particular shade tollerant
    //printf("ddalmo c->tree_layers[layer].daily_layer_cover_proj   = %f \n",c->tree_layers[layer].daily_layer_cover_proj );
	//printf("ddalmo s->value[DAILY_CANOPY_COVER_PROJ]  = %f \n",s->value[DAILY_CANOPY_COVER_PROJ]);

   if (c->tree_layers[layer].daily_layer_cover_proj > 1.)

   {
	scal_cover = 1/c->tree_layers[layer].daily_layer_cover_proj  ;
   }   else  {

   scal_cover = 1. ;
   }
   //printf(" IN CANOPY RAD  layer  %d!!!\n", layer );
   //printf(" IN CANOPY RAD SCALE COVER %f!!!\n", scal_cover );
   //printf(" IN CANOPY RAD c->tree_layers[layer].daily_layer_cover_proj  %f!!!\n", c->tree_layers[layer].daily_layer_cover_proj  );
 //printf(" IN CANOPY RAD2 c->tree_layers[layer].daily_layer_cover_proj  %f !!!\n", c->tree_layers[layer].layer_cover_proj );

	/* light reflection, absorption and transmission */
	logger(g_debug_log,"\n*Light reflection, absorption and transmission*\n");

	/***********************************************************************************************/
	/* PAR computation (molPAR/m2 covered/day) */
	logger(g_debug_log,"\n-PAR-\n");

	/** available par **/
	s->value[PAR]             = meteo_daily->par  * s->value[DAILY_CANOPY_COVER_PROJ]*scal_cover;
    temp_scale=s->value[DAILY_CANOPY_COVER_PROJ]*scal_cover;

   //printf("in canopy radiation  meteo_daily->par              = %g  \n", meteo_daily->par) ;
   //printf("in canopy radiation  s->value[PAR]              = %g  \n", s->value[PAR]) ;
     // printf("in canopy radiation  s->value[DAILY_CANOPY_COVER_PROJ]            = %g  \n", s->value[DAILY_CANOPY_COVER_PROJ]) ;
//printf("in canopy radiation TEMP SCALE          = %g  \n", temp_scale ) ;


	/** sun leaves **/
	s->value[PAR_REFL_SUN]     = s->value[PAR]     * Light_refl_par_frac_sun;
	s->value[PAR_SUN]          = s->value[PAR]     - s->value[PAR_REFL_SUN];
	s->value[APAR_SUN]         = s->value[PAR_SUN] * Light_abs_frac_sun;
	s->value[TRANSM_PAR_SUN]   = s->value[PAR_SUN] - s->value[APAR_SUN];

	/* check PAR balance for sun leaves */
	CHECK_CONDITION ( fabs ( ( s->value[PAR_SUN] - s->value[TRANSM_PAR_SUN] ) - s->value[APAR_SUN] ), >, eps );
	CHECK_CONDITION ( fabs ( ( s->value[PAR]     - s->value[TRANSM_PAR_SUN] ) - ( s->value[PAR_REFL_SUN] + s->value[APAR_SUN] ) ), >, eps );

	/** shaded leaves **/
	s->value[PAR_REFL_SHADE]   = s->value[TRANSM_PAR_SUN] * Light_refl_par_frac_shade;
	s->value[PAR_SHADE]        = s->value[TRANSM_PAR_SUN] - s->value[PAR_REFL_SHADE];
	s->value[APAR_SHADE]       = s->value[PAR_SHADE]      * Light_abs_frac_shade;
	s->value[TRANSM_PAR_SHADE] = s->value[PAR_SHADE]      - s->value[APAR_SHADE];

	/* check PAR balance for shaded leaves */
	CHECK_CONDITION ( fabs ( ( s->value[PAR_SHADE]      - s->value[TRANSM_PAR_SHADE] ) - s->value[APAR_SHADE] ), >, eps );
	CHECK_CONDITION ( fabs ( ( s->value[TRANSM_PAR_SUN] - s->value[TRANSM_PAR_SHADE] ) - ( s->value[PAR_REFL_SHADE] + s->value[APAR_SHADE] ) ), >, eps );

	/** overall canopy **/

	s->value[APAR]             = s->value[APAR_SUN]     + s->value[APAR_SHADE];
	s->value[PAR_REFL]         = s->value[PAR_REFL_SUN] + s->value[PAR_REFL_SHADE];
	s->value[TRANSM_PAR]       = s->value[TRANSM_PAR_SHADE];

    //printf("in canopy radiation  s->value[APAR]  = %g  \n", s->value[APAR] ) ;

//printf("in canopy radiation s->value[PAR_REFL]  = %g  \n", s->value[PAR_REFL]) ;

//printf("in canopy radiation s->value[TRANSM_PAR]  = %g  \n", s->value[TRANSM_PAR] ) ;

	logger(g_debug_log,"APAR = %f par/m2\n", s->value[APAR]);

	/* compute fAPAR */
	if ( s->value[PAR] )
	{
		s->value[fAPAR]        = s->value[APAR] / s->value[PAR];
	}
	else
	{
		s->value[fAPAR]        = 0.;
	}

	s->value[YEARLY_APAR]     += s->value[APAR];

	/* check */
	CHECK_CONDITION(s->value[APAR],       <, ZERO );
	CHECK_CONDITION(s->value[TRANSM_PAR], <, ZERO );
	CHECK_CONDITION(s->value[APAR] + s->value[TRANSM_PAR], <, ZERO );
	CHECK_CONDITION(fabs((s->value[APAR] + s->value[TRANSM_PAR] + s->value[PAR_REFL] )-s->value[PAR]), >, eps);


	/***********************************************************************************************/
	/* Short Wave computation (W/m2 covered) */
	logger(g_debug_log,"\n-Short Wave-\n");

	/** available Short Wave **/
	s->value[SW_RAD]              = meteo_daily->sw_downward_W * s->value[DAILY_CANOPY_COVER_PROJ]*scal_cover;

	/** sun leaves **/
	s->value[SW_RAD_REFL_SUN]     = s->value[SW_RAD]     * Light_refl_sw_frac_sun;
	s->value[SW_RAD_SUN]          = s->value[SW_RAD]     - s->value[SW_RAD_REFL_SUN];
	s->value[SW_RAD_ABS_SUN]      = s->value[SW_RAD_SUN] * Light_abs_frac_sun;
	s->value[SW_RAD_TRANSM_SUN]   = s->value[SW_RAD_SUN] - s->value[SW_RAD_ABS_SUN];

	/* check Short Wave balance for sun leaves */
	CHECK_CONDITION ( fabs ( ( s->value[SW_RAD_SUN] - s->value[SW_RAD_TRANSM_SUN] ) - s->value[SW_RAD_ABS_SUN] ), >, eps );
	CHECK_CONDITION ( fabs ( ( s->value[SW_RAD]     - s->value[SW_RAD_TRANSM_SUN] ) - ( s->value[SW_RAD_REFL_SUN] + s->value[SW_RAD_ABS_SUN] ) ), >, eps );

	/** shaded leaves **/
	s->value[SW_RAD_REFL_SHADE]   = s->value[SW_RAD_TRANSM_SUN] * Light_refl_sw_frac_shade;
	s->value[SW_RAD_SHADE]        = s->value[SW_RAD_TRANSM_SUN] - s->value[SW_RAD_REFL_SHADE];
	s->value[SW_RAD_ABS_SHADE]    = s->value[SW_RAD_SHADE]      * Light_abs_frac_shade;
	s->value[SW_RAD_TRANSM_SHADE] = s->value[SW_RAD_SHADE]      - s->value[SW_RAD_ABS_SHADE];

	/* check Short Wave balance for shaded leaves */
	CHECK_CONDITION ( fabs ( ( s->value[SW_RAD_SHADE]      - s->value[SW_RAD_TRANSM_SHADE] ) - s->value[SW_RAD_ABS_SHADE] ), >, eps );
	CHECK_CONDITION ( fabs ( ( s->value[SW_RAD_TRANSM_SUN] - s->value[SW_RAD_TRANSM_SHADE] ) - ( s->value[SW_RAD_REFL_SHADE] + s->value[SW_RAD_ABS_SHADE] ) ), >, eps );

	/** overall canopy **/
	s->value[SW_RAD_ABS]          = s->value[SW_RAD_ABS_SUN]  + s->value[SW_RAD_ABS_SHADE];
	s->value[SW_RAD_REFL]         = s->value[SW_RAD_REFL_SUN] + s->value[SW_RAD_REFL_SHADE];
	s->value[SW_RAD_TRANSM]       = s->value[SW_RAD_TRANSM_SHADE];

	logger(g_debug_log,"SW_RAD_ABS = %f W/m2\n", s->value[SW_RAD_ABS]);

//    printf("SW_RAD_ABS  species %s!!!\n", s->name);
//	printf("ddalmo SW_RAD  = %f W/m2\n", s->value[SW_RAD]);
//	printf("ddalmo meteo_daily->sw_downward_W  = %f W/m2\n", meteo_daily->sw_downward_W);
//	printf("ddalmo DAILY_CANOPY_COVER_PROJ  = %f \n", s->value[DAILY_CANOPY_COVER_PROJ]);
//	printf("ddalmo SW_RAD_ABS  = %f W/m2\n", s->value[SW_RAD_ABS]);
//	printf("ddalmo s->value[SW_RAD_TRANSM]   = %f W/m2\n", s->value[SW_RAD_TRANSM] );
//	printf("ddalmo s->value[SW_RAD_REFL]   = %f W/m2\n", s->value[SW_RAD_REFL] );

	/* check */
	CHECK_CONDITION(s->value[SW_RAD_ABS],    <, ZERO );
	CHECK_CONDITION(s->value[SW_RAD_TRANSM], <, ZERO);
	CHECK_CONDITION(s->value[SW_RAD_ABS] + s->value[SW_RAD_TRANSM], <, ZERO );
	CHECK_CONDITION(fabs((s->value[SW_RAD_ABS] + s->value[SW_RAD_TRANSM] + s->value[SW_RAD_REFL] ) - s->value[SW_RAD]), >, eps);

	/***********************************************************************************************/
	//test 25 June 2017: currently use and modify on february 2022
	
#if 1
	// Net Radiation computation (W/m2 covered) computed with 3-PG method //
	// see Landsberg Terrestrial ecology Chapter 2: Weather and Energy Balance //
	
	// Note: this is a NET RADIATION, hence reflection is already included (no need to compute it)
	// This is a very rough approximation and has to be FIXME
		
	logger(g_debug_log,"\n-Net Radiation-\n");

	/** available Net Radiation **/
	s->value[NET_RAD]              = meteo_daily->Net_rad_threePG * s->value[DAILY_CANOPY_COVER_PROJ]*scal_cover;

        // to avoid negative Transpiration values (e.g. at night transpiration is assumed 0)
        // yet, net radiation can be negative and lead to e.g. dew formation)
        
        if ( s->value[NET_RAD] <= 0. )
	{
				s->value[NET_RAD] = 0.;
	}
			
	/** sun leaves **/
	s->value[NET_RAD_REFL_SUN]     = 0.  ; // s->value[NET_RAD] * Light_refl_sw_frac_sun;
	s->value[NET_RAD_SUN]          = s->value[NET_RAD]     - s->value[NET_RAD_REFL_SUN];
	s->value[NET_RAD_ABS_SUN]      = s->value[NET_RAD_SUN] * Light_abs_frac_sun;
	s->value[NET_RAD_TRANSM_SUN]   = s->value[NET_RAD_SUN] - s->value[NET_RAD_ABS_SUN];

	/* check Net Radiation balance for sun leaves */
	CHECK_CONDITION ( fabs ( ( s->value[NET_RAD_SUN] - s->value[NET_RAD_TRANSM_SUN] ) - s->value[NET_RAD_ABS_SUN] ), >, eps );
	CHECK_CONDITION ( fabs ( ( s->value[NET_RAD]     - s->value[NET_RAD_TRANSM_SUN] ) - ( s->value[NET_RAD_REFL_SUN] + s->value[NET_RAD_ABS_SUN] ) ), >, eps );

	/** shaded leaves **/
	s->value[NET_RAD_REFL_SHADE]   = 0. ; // s->value[NET_RAD_TRANSM_SUN] * Light_refl_sw_frac_shade;
	s->value[NET_RAD_SHADE]        = s->value[NET_RAD_TRANSM_SUN] - s->value[NET_RAD_REFL_SHADE];
	s->value[NET_RAD_ABS_SHADE]    = s->value[NET_RAD_SHADE]      * Light_abs_frac_shade;
	s->value[NET_RAD_TRANSM_SHADE] = s->value[NET_RAD_SHADE]      - s->value[NET_RAD_ABS_SHADE];

	/* check Net Radiation balance for shaded leaves */
	CHECK_CONDITION ( fabs ( ( s->value[NET_RAD_SHADE]      - s->value[NET_RAD_TRANSM_SHADE] ) - s->value[NET_RAD_ABS_SHADE] ), >, eps );
	CHECK_CONDITION ( fabs ( ( s->value[NET_RAD_TRANSM_SUN] - s->value[NET_RAD_TRANSM_SHADE] ) - ( s->value[NET_RAD_REFL_SHADE] + s->value[NET_RAD_ABS_SHADE] ) ), >, eps );

	/** overall canopy **/
	s->value[NET_RAD_ABS]          = s->value[NET_RAD_ABS_SUN]  + s->value[NET_RAD_ABS_SHADE];
	s->value[NET_RAD_REFL]         = s->value[NET_RAD_REFL_SUN] + s->value[NET_RAD_REFL_SHADE];
	s->value[NET_RAD_TRANSM]       = s->value[NET_RAD_TRANSM_SHADE];

	logger(g_debug_log,"NET_RAD_ABS = %f W/m2\n", s->value[NET_RAD_ABS]);

	/* check */
	
	CHECK_CONDITION(s->value[NET_RAD_ABS],    <, ZERO );   // NET_RAD forced to 0 if negative to avoid Transpiration values <0 
	CHECK_CONDITION(s->value[NET_RAD_TRANSM], <, ZERO);
	CHECK_CONDITION(s->value[NET_RAD_ABS] + s->value[NET_RAD_TRANSM], <, ZERO );
	CHECK_CONDITION(fabs((s->value[NET_RAD_ABS] + s->value[NET_RAD_TRANSM] + s->value[NET_RAD_REFL] )-s->value[NET_RAD]), >, eps);
	
#endif

	/* it follows rationale of BIOME-BGC to obtain m2 instead m2/m2 */
	//fixme then recompute transmitted fraction!!!!!!!!!!!!
	/*
	if(s->value[PPFD_ABS_SHADE] < 0.0)
	{
		s->value[PPFD_ABS_SHADE]  = 0.0;
	}
	if(s->value[LAI_SUN] > 0.0 && s->value[LAI_SHADE] > 0.0)
	{
		s->value[PPFD_ABS_SUN] /= s->value[LAI_SUN];
		s->value[PPFD_ABS_SHADE] /= s->value[LAI_SHADE];
	}
	else
	{
		s->value[PPFD_ABS_SUN] = s->value[PPFD_ABS_SHADE] = 0.0;

	}
	s->value[PPFD_ABS] = s->value[PPFD_ABS_SUN] + s->value[PPFD_ABS_SHADE];
	 */
}

void canopy_radiation_sw_band(cell_t *const c, const int layer, const int height, const int dbh, const int age, const int species,  meteo_daily_t *const meteo_daily)
{
	double Light_abs_frac, Light_abs_frac_sun, Light_abs_frac_shade;       /* (ratio) fraction of PAR and Short Wave radiation absorbed */
	double Light_trasm_frac, Light_trasm_frac_sun, Light_trasm_frac_shade; /* (ratio) fraction of PAR and Short Wave radiation transmitted */

	double Light_refl_par_frac;                                            /* (ratio) fraction of PAR reflected from canopy */
	double Light_refl_par_frac_sun;                                        /* (ratio) fraction of PAR reflected from sun leaves */
	double Light_refl_par_frac_shade;                                      /* (ratio) fraction of PAR reflected from shaded leaves */
	double Light_refl_sw_frac;                                             /* (ratio) fraction of Short Wave radiation reflected from the canopy */
	double Light_refl_sw_frac_sun;                                         /* (ratio) fraction of Short Wave radiation reflected from sun leaves */
	double Light_refl_sw_frac_shade;                                       /* (ratio) fraction of Short Wave radiation reflected from shaded leaves */
	double Light_refl_net_rad_frac;                                        /* (ratio) fraction of Net radiation reflected from the canopy */
	double Light_refl_net_rad_frac_sun;                                    /* (ratio) fraction of Net radiation reflected from sun leaves */
	double Light_refl_net_rad_frac_shade;                                  /* (ratio) fraction of Net radiation reflected from shaded leaves */
	double k;
	//double k_eff;

	tree_layer_t *l;
	species_t *s;

	l = &c->tree_layers[layer];
	s = &c->heights[height].dbhs[dbh].ages[age].species[species];

	//following Ritchie et al., 1998 and Hydi et al.,2016
	//double actual_albedo;

	logger(g_debug_log, "\n**SHORT WAVE BAND RADIATION**\n");

	/***********************************************************************************************************/

	k = s->value[K];

	/***********************************************************************************************************/
	/* SHORT WAVE RADIATION FRACTIONS */
	/* compute fractions of light intercepted, transmitted and reflected from the canopy */

	/** fraction of light transmitted through the canopy **/
	/* note: we currently use approach for homogeneous canopies that improves representation when canopy is not closed: */
	/* see method from Cannel and Grace, Can. J. For: Res. Vol. 23, 1993 [Eq.8]  */
	/* see method from Duursma and Makela, Tree Phys: Vol. 27, 859-870, 2007 [Eq.4] */
    
        //5p606
#if 0
	//FIXME this is a very bad assumption (I know!!) to avoid that under low density gpp is too low
	Light_trasm_frac       = exp ( - k * s->value[LAI_PROJ]) * s->value[DAILY_CANOPY_COVER_PROJ];
	Light_trasm_frac_sun   = exp ( - k * s->value[LAI_SUN_PROJ]) * s->value[DAILY_CANOPY_COVER_PROJ];
	Light_trasm_frac_shade = exp ( - k * s->value[LAI_SHADE_PROJ]) * s->value[DAILY_CANOPY_COVER_PROJ] ;
#else
	Light_trasm_frac       = exp ( - k * s->value[LAI_PROJ]);
	Light_trasm_frac_sun   = exp ( - k * s->value[LAI_SUN_PROJ]);
	Light_trasm_frac_shade = exp ( - k * s->value[LAI_SHADE_PROJ]);
#endif

	/** fraction of light absorbed by the canopy **/
	Light_abs_frac         = 1. - Light_trasm_frac;
	Light_abs_frac_sun     = 1. - Light_trasm_frac_sun;
	Light_abs_frac_shade   = 1. - Light_trasm_frac_shade;

	/** fraction of light reflected by the canopy **/
	/* for Short Wave radiation and PAR */
	/* following BIOME albedo for PAR is 1/3 of albedo. the absorbed PAR is */
	/* calculated similarly to sw except that albedo is 1/3 for PAR because less */
	/* PAR is reflected than sw_radiation (Jones 1992) */

	Light_refl_sw_frac             = s->value[ALBEDO];
	Light_refl_sw_frac_sun         = s->value[ALBEDO]         /* * Light_abs_frac_sun */;
	Light_refl_sw_frac_shade       = s->value[ALBEDO]         /* * Light_abs_frac_shade */;

	Light_refl_net_rad_frac        = s->value[ALBEDO];
	Light_refl_net_rad_frac_sun    = s->value[ALBEDO]         /* * Light_abs_frac_sun */;
	Light_refl_net_rad_frac_shade  = s->value[ALBEDO]         /* * Light_abs_frac_shade */;

	Light_refl_par_frac            = (s->value[ALBEDO] / 3. );
	Light_refl_par_frac_sun        = (s->value[ALBEDO] / 3. ) /* * Light_abs_frac_sun */;
	Light_refl_par_frac_shade      = (s->value[ALBEDO] / 3. ) /* * Light_abs_frac_shade */;

	//fixme set that if gapcover is bigger then 0.5 albedo should be considered also in dominated layer!!!!
	//fixme following MAESPA (Duursma et al.,) and from Campbell & Norman (2000, p. 259) dominated layers should have just shaded leaves

	/* RADIATION */
	/*****************************************************************************************************************/
	/* first height class in the cell is processed */
	if( ! l->layer_height_class_counter && ! c->cell_height_class_counter )
	{

	//	 printf(" ************************ first height class in the cell is processed  \n");
		/* reset temporary values when the first height class in layer is processed */
		c->temp_apar         = 0.;
		c->temp_par_refl     = 0.;
		c->temp_sw_rad_abs   = 0.;
		c->temp_sw_rad_refl  = 0.;
		c->temp_net_rad_abs  = 0.;
		c->temp_net_rad_refl = 0.;
	}
	/*****************************************************************************************************************/

	/* shared functions among all class/layers */
	/* counters */
	++l->layer_height_class_counter;
	++c->cell_height_class_counter;

	/*************************************************************************/
	/* compute reflected, absorbed and transmitted Par, Short Wave radiation and PPFD class level */
	canopy_sw_band_abs_trans_refl_radiation (c, layer, height, dbh, age, species, meteo_daily, Light_abs_frac, Light_abs_frac_sun, Light_abs_frac_shade,
			Light_refl_par_frac, Light_refl_par_frac_sun, Light_refl_par_frac_shade, Light_refl_sw_frac, Light_refl_sw_frac_sun, Light_refl_sw_frac_shade);

	/*************************************************************************/
	/* temporary absorbed and reflected values */

	/* update temporary absorbed and reflected PAR for lower layer */
	c->temp_apar       += s->value[APAR];
	c->apar            += s->value[APAR];
	//	 printf(" s->value[APAR]             = %g \n ", s->value[APAR]);
	// printf(" c->apar            = %g \n ",c->apar);

	logger(g_debug_log,"cum apar = %f\n", c->apar);

	c->temp_par_refl   += s->value[PAR_REFL];
	c->par_refl        += s->value[PAR_REFL];
	// printf(" s->value[PAR_REFL]             = %g \n ", s->value[PAR_REFL]);
	// printf(" c->par_refl            = %g \n ",c->par_refl );
	logger(g_debug_log,"cum par_refl = %f\n", c->par_refl);

	/* update temporary absorbed and transmitted Short Wave radiation lower layer */
	c->temp_sw_rad_abs  += s->value[SW_RAD_ABS];
	c->sw_rad_abs       += s->value[SW_RAD_ABS];
	logger(g_debug_log,"cum sw_rad_abs = %f\n", c->sw_rad_abs);

	c->temp_sw_rad_refl += s->value[SW_RAD_REFL];
	c->sw_rad_refl      += s->value[SW_RAD_REFL];
	logger(g_debug_log,"cum sw_rad_refl = %f\n", c->sw_rad_refl);

	//	printf("ddalmo c->temp_sw_rad_abs  = %f W/m2\n", c->temp_sw_rad_abs);
	//		printf("ddalmo c->sw_rad_refl   = %f W/m2\n",c->sw_rad_refl );

	/* update temporary absorbed and transmitted net radiation lower layer */
	c->temp_net_rad_abs += s->value[NET_RAD_ABS];
	c->net_rad_abs      += s->value[NET_RAD_ABS];
	logger(g_debug_log,"cum net_rad_abs = %f\n", c->net_rad_abs);

	c->temp_net_rad_refl += s->value[NET_RAD_REFL];
	c->net_rad_refl      += s->value[NET_RAD_REFL];
	logger(g_debug_log,"cum net_rad_refl = %f\n", c->net_rad_refl);

	/*****************************************************************************************************************/
	/* when it matches the last height class in the layer is processed */

	// printf("  l->layer_n_height_class           = %d \n ", l->layer_n_height_class);
	//  printf(" l->layer_height_class_counter        = %d \n ",l->layer_height_class_counter );
 //printf(" c->cell_height_class_counter       = %d \n ",c->cell_height_class_counter );

	if ( l->layer_n_height_class == l->layer_height_class_counter )
	{
		logger(g_debug_log,"\n************************************\n");
		logger(g_debug_log,"last height class in layer processed\n");
		logger(g_debug_log,"update radiation values for lower layer\n");


   // printf(" ********last height class in layer processed ************    \n " );

		/* compute values for lower layer when last height class in layer is processed */
		/* compute par for lower layer */
	//	printf(" BEFORE END LAYER PAR TO THE LOWER LAYER meteo_daily->par  = %f \n", meteo_daily->par );

		meteo_daily->par           -= (c->temp_apar + c->temp_par_refl);

	//	printf(" END LAYER PAR TO THE LOWER LAYER meteo_daily->par  = %f \n", meteo_daily->par );
      //  printf(" END LAYER PAR TO THE LOWER LAYER c->temp_apar  = %f \n", c->temp_apar );
      //  printf(" END LAYER PAR TO THE LOWER LAYER c->temp_par_refl = %f \n", c->temp_par_refl);

        //printf(" SONO ULTIMO ALTEZZA NEL LAYER \n") ;
		/* compute Short Wave radiation for lower layesr */
		meteo_daily->sw_downward_W -= (c->temp_sw_rad_abs + c->temp_sw_rad_refl);

		/* compute Net radiation for lower layesr */
		meteo_daily->Net_rad_threePG -= (c->temp_net_rad_abs + c->temp_net_rad_refl);

		logger(g_debug_log, "Radiation for lower layers\n");
		logger(g_debug_log, "PAR                   = %f molPAR/m2/day\n", meteo_daily->par);
		logger(g_debug_log, "Short Wave            = %f W/m2\n", meteo_daily->sw_downward_W);
		logger(g_debug_log, "Short Net_rad_threePG = %f W/m2\n", meteo_daily->Net_rad_threePG);

		/* reset temporary values when the last height class in layer is processed */
		c->temp_apar         = 0.;
		c->temp_par_refl     = 0.;
		c->temp_sw_rad_abs   = 0.;
		c->temp_sw_rad_refl  = 0.;
		c->temp_net_rad_abs  = 0.;
		c->temp_net_rad_refl = 0.;

		/* reset counter */
		l->layer_height_class_counter = 0;
	}

    
	/*************************************************************************/
	/* when it matches the last height class in the cell is processed */
	//fixme sometimes it doesn't go in caused by the a jump in "cell_height_class_counter"
	//as it is now is used just for print data but it should be fixed
	if ( c->heights_count == c->cell_height_class_counter )
	{
		logger(g_debug_log,"\n************************************\n");
		logger(g_debug_log, "last height class in cell processed\n");
		logger(g_debug_log, "Radiation for soil\n");
		logger(g_debug_log, "PAR                   = %f molPAR/m2/day\n", meteo_daily->par);
		logger(g_debug_log, "Short Wave            = %f W/m2\n", meteo_daily->sw_downward_W);
		logger(g_debug_log, "Net_rad_threePG       = %f W/m2\n", meteo_daily->Net_rad_threePG);
		logger(g_debug_log,"\n***********************************\n");
	}
	/*****************************************************************************************************************/
}
