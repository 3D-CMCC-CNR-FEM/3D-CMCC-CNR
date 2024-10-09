/*
 * canopy_evapotranspiration.c
 *
 *  Created on: 23/mar/2016
 *      Author: alessio
 */


/* includes */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "common.h"
#include "settings.h"
#include "canopy_evapotranspiration.h"
#include "canopy_interception.h"
#include "constants.h"
#include "logger.h"
#include "Penman_Monteith.h"
#include "heat_fluxes.h"

extern logger_t* g_debug_log;
extern settings_t* g_settings;

#define TEST 0
/*************************************************************************************************************************/

void gs_Jarvis (species_t *const s, const double gl_x, const double g_corr)
{
	double m_final;
	double m_final_sun;
	double m_final_shade;

	/* following Jarvis 1997 + Frank et al., 2013 + Hidy et al., 2016 GMD */

	/* leaf stomatal conductance: first generate multipliers, then apply them to maximum stomatal conductance */
	/* apply all multipliers to the maximum stomatal conductance */
	/* differently from BIOME we use F_T that takes into account not only minimum temperature effects */
	/* differently from BIOME we use also F_AGE */

	/* original BIOME version (Jarvis method) */
	m_final       = s->value[F_LIGHT]       * s->value[F_SW] * s->value[F_T] * s->value[F_VPD] * s->value[F_AGE];
	m_final_sun   = s->value[F_LIGHT_SUN]   * s->value[F_SW] * s->value[F_T] * s->value[F_VPD] * s->value[F_AGE];
	m_final_shade = s->value[F_LIGHT_SHADE] * s->value[F_SW] * s->value[F_T] * s->value[F_VPD] * s->value[F_AGE];

	/* check */
	//if (m_final       < eps) m_final       = eps;
	//if (m_final_sun   < eps) m_final_sun   = eps;
	//if (m_final_shade < eps) m_final_shade = eps;
	
	//  Set m_corr to 0 when m_final is == eps 
        //  this avoid having gpp values larger than 0 in winter 
        // e.g. in boreal forests
#if 1
	
        if (m_final       < eps) {m_final       = eps;
                                  s->value[m_corr]= 0. ;
                                  }
	if (m_final_sun   < eps) {m_final_sun       = eps;
                                  s->value[m_corr_sun]= 0.;
                                  }
	if (m_final_shade < eps) {m_final_shade       = eps;
                                  s->value[m_corr_shade] =0.;
                                  }
        
#else 
         
         s->value[m_corr]= 1.;
         s->value[m_corr_sun]= 1.;
         s->value[m_corr_shade] =1.;
          
#endif
         
	s->value[STOMATAL_CONDUCTANCE]       = gl_x * m_final       * g_corr;
	s->value[STOMATAL_SUN_CONDUCTANCE]   = gl_x * m_final_sun   * g_corr;
	s->value[STOMATAL_SHADE_CONDUCTANCE] = gl_x * m_final_shade * g_corr;
}

/*************************************************************************************************************************/

void gs_Ball_Berry (species_t *const s, const double gl_x, const double g_corr)
{


}

/*************************************************************************************************************************/

void canopy_evapotranspiration(cell_t *const c, const int layer, const int height, const int dbh, const int age, const int species, meteo_daily_t *const meteo_daily)
{
	double g_corr;
	double gl_bl;
	double gl_x;                                                          /* maximum stomatal conductance */
	//double gl_s, gl_s_sun, gl_s_shade;
	double gl_c;
	//double m_final, m_final_sun, m_final_shade;
	double gl_e_wv;
	//double gl_t_wv, gl_t_wv_sun, gl_t_wv_shade;
	double gl_sh;
	double gc_e_wv;
	double gc_sh;
	double net_rad;
	double rv;
	double rh;                                                             /* Leaf-Canopy resistance to sensible heat */
	double leaf_transp;                                                    /* (kg/m2/s) leaf transpiration flux */
	double evap_daylength_sec;                                             /* (sec) daylength time need for canopy rain evaporation */
	double transp_daylength_sec;                                           /* (sec) daylength time need for canopy transpiration */
	double subl_melt_daylength_sec;                                        /* (sec) daylength time need for canopy snow sublimation or melting */
	double evapo;                                                          /* (kg/m2/s) rain evaporation or leaf transpiration flux */
	double melt;                                                           /* (kg/m2/s) snow sublimation flux */
	double subl;                                                           /* (kg/m2/s) snow melt flux */
	double subl_melt;                                                      /* (kg/m2/s) snow sublimation or melt flux */


	species_t *s;
	s = &c->heights[height].dbhs[dbh].ages[age].species[species];


	/* it mainly follows rationale and algorithms of BIOME-BGC v.4.2 */
	/* it computes canopy interception, evaporation and transpiration */
	/********************************************************************************************************/

	logger(g_debug_log, "\n**CANOPY EVAPO-TRANSPIRATION**\n");

	/********************************************************************************************************/

	/* call canopy interception */
	canopy_interception( c, layer, height, dbh, age, species, meteo_daily );

	/********************************************************************************************************/
	/* temperature and pressure correction factor for conductances */
	g_corr = pow( ( meteo_daily->tday + TempAbs ) / 293.15, 1.75) * 101300. / meteo_daily->air_pressure;

	/* calculate leaf- and canopy-level conductances to water vapor and
		sensible heat fluxes */

	/* leaf aerodynamic-boundary-layer conductance */
	gl_bl = s->value[BLCOND] * g_corr;

#if 0
	/**************************************************************************/
	/* aerodynamic boundary layer conductance corrected for wind speed */
	/* following Shi et al., 2008 Journal of Biophysical Research */
	if ( meteo_daily->windspeed != NO_DATA)
	{
		double zero = 20;    /* height measurement of wind speed (m) */ //fixme add to settings.txt data
		double rough;        /* roughness length */
		double plane;        /* zero place displacement */
		height_t *h;
		h = &c->heights[height];

		/* compute roughness length (m) */
		rough = 0.032 * h->value;

		/* compute zero plane displacement (m) */
		plane = 0.88 * h->value;

		/* aerodynamic boundary layer conductance */
		gl_bl = pow ( KARM , 2.) / pow ( log ( ( zero - plane) / rough ) , 2. );

		//gl_bl *= g_corr;
		//printf("wind speed = %f gl_bl corrected for wind speed and corrected for g_corr %f\n", meteo_daily->windspeed, gl_bl);

	}
	/**************************************************************************/
#endif

	/* canopy boundary layer conductance */
	s->value[CANOPY_BLCOND] = gl_bl * s->value[LAI_PROJ];

	/* leaf cuticular conductance corrected for temperature and air pressure */
	gl_c = s->value[CUTCOND] * g_corr;

	if ( g_settings->PSN_mod != 2 )
	{
		/** correction for [CO2] stomatal conductance **/
		/* correct maximum stomatal conductance for CO2 concentration */
		gl_x = (s->value[F_CO2_TR] / 0.9116) * s->value[MAXCOND];

		/** Jarvis effective stomatal conductance */
		gs_Jarvis (s, gl_x, g_corr);
	}
	else
	{
		/** Ball-Woodrow-Berry effective stomatal conductance */
		//gs_Ball_Berry (s, gl_x, g_corr);
	}

	/* calculate leaf-and canopy-level conductances to water vapor and
		sensible heat fluxes, to be used in Penman-Monteith calculations of
		canopy evaporation and canopy transpiration */

	/* Leaf conductance to evaporated water vapor, per unit projected LAI */
	gl_e_wv       = gl_bl;

	/** leaf conductance **/
	/* Leaf conductance to transpired water vapor, per unit projected
		LAI.  This formula is derived from stomatal and cuticular conductances
		in parallel with each other, and both in series with leaf boundary
		layer conductance. */
	s->value[LEAF_CONDUCTANCE]       = (gl_bl * (s->value[STOMATAL_CONDUCTANCE]       + gl_c)) / (gl_bl + s->value[STOMATAL_CONDUCTANCE]       + gl_c);
	s->value[LEAF_SUN_CONDUCTANCE]   = (gl_bl * (s->value[STOMATAL_SUN_CONDUCTANCE]   + gl_c)) / (gl_bl + s->value[STOMATAL_SUN_CONDUCTANCE]   + gl_c);
	s->value[LEAF_SHADE_CONDUCTANCE] = (gl_bl * (s->value[STOMATAL_SHADE_CONDUCTANCE] + gl_c)) / (gl_bl + s->value[STOMATAL_SHADE_CONDUCTANCE] + gl_c);

	/* Leaf conductance to sensible heat, per unit all-sided LAI */
	gl_sh         = gl_bl;

	/* Canopy conductance to evaporated water vapor */
	gc_e_wv       = gl_e_wv * s->value[LAI_PROJ];

	/* Canopy conductance to sensible heat */
	gc_sh         = gl_sh   * s->value[LAI_PROJ];

	/* Canopy evaporation, if any water was intercepted */
	/* Calculate Penman-Monteith evaporation, given the canopy conductances to
	evaporated water and sensible heat. Calculate the time required to
	evaporate all the canopy water at the daily average conditions, and
	subtract that time from the daylength to get the effective daylength for
	transpiration. */

	/* note: as in Campbell and Norman, Environmental Biophysics,
	model should use Absorbed Radiation (instead just incident) but
	composed of Short and Long wave radiation */
	
	// comment. TODO consider dew formation (when net-radiation is negative)

	if( s->value[LAI_PROJ] > 0. )
	{
		/********************************************************************************************************************************/
		/* if canopy has water */
		if( s->value[CANOPY_WATER] > 0. )
		{
			/* convert conductance to resistance */
			rv = 1. / gc_e_wv;
			rh = 1. / gc_sh;

			/* radiation */
			/* note: Biome-bgc used SW radiation */
#if 0
			net_rad = s->value[SW_RAD_ABS];
#else
			net_rad = s->value[NET_RAD_ABS];
#endif

			/* call Penman-Monteith function, it returns Potential evaporation in kg/m2/s for evaporation and W/m2 for latent heat */
			evapo = Penman_Monteith ( meteo_daily, rv, rh, net_rad );

			/* check for negative values */
			if( evapo < 0. ) evapo  = 0.;    //FIXME dew formation

			s->value[CANOPY_EVAPO] = evapo;

			/* calculate the time required to evaporate all the canopy water */
			evap_daylength_sec     = s->value[CANOPY_WATER] / s->value[CANOPY_EVAPO];

			/* day not long enough to evap. all int. water */
			if( evap_daylength_sec > meteo_daily->daylength_sec )
			{
				/* NO Canopy transpiration (Canopy completely Wet) */

				++ s->counter[DAYS_WITH_CANOPY_WET];

				/* adjust day length for transpiration */
				transp_daylength_sec             = 0.;
				s->value[CANOPY_FRAC_DAY_TRANSP] = 0.;

				/* no time left for transpiration */
				s->value[CANOPY_TRANSP] = 0.;

				/* day length limits canopy evaporation */
				evapo *= meteo_daily->daylength_sec;
				s->value[CANOPY_EVAPO] = evapo;

				s->value[CANOPY_EVAPO_TRANSP] = s->value[CANOPY_EVAPO] + s->value[CANOPY_TRANSP];

				/* remove evaporated water from canopy water pool */
				s->value[CANOPY_WATER] -= s->value[CANOPY_EVAPO];

			}
			/* all intercepted water evaporated and there's time for transpiration */
			else
			{
				/* all intercepted water evaporated */
				s->counter[DAYS_WITH_CANOPY_WET] = 0;

				/* all canopy water evaporates */
				evapo = s->value[CANOPY_WATER];

				s->value[CANOPY_EVAPO]  = evapo;

				/* reset canopy water from canopy water pool */
				s->value[CANOPY_WATER] -= s->value[CANOPY_EVAPO];

				/* Canopy transpiration (Partial Canopy Wet) */

				/* adjust day length for transpiration */
				transp_daylength_sec   = meteo_daily->daylength_sec - evap_daylength_sec;

				s->value[CANOPY_FRAC_DAY_TRANSP] = transp_daylength_sec / meteo_daily->daylength_sec;

				/* Leaf-Canopy resistance to sensible heat */
				rh = 1. / gl_sh;

				/************************************************************************************/
				/* calculate transpiration using adjusted day length */
				/** LAI SUN **/

				/* resistance for sun canopy fraction */
				rv = 1. / s->value[LEAF_SUN_CONDUCTANCE];

				/* note: Net Rad is Short wave flux in Biome-bgc*/
				/* convert radiation to stomatal scale */
				// FIXME: compute correctly the net_radiation
#if 0
				net_rad = ( s->value[SW_RAD_ABS_SUN]  / s->value[LAI_SUN_PROJ] ) * s->value[F_LIGHT_SUN_MAKELA];
#else
				net_rad = ( s->value[NET_RAD_ABS_SUN] / s->value[LAI_SUN_PROJ] ) * s->value[F_LIGHT_SUN_MAKELA];
#endif
				/* call Penman-Monteith function, returns e in kg/m2/s for transpiration and W/m2 for latent heat */
				//fixme use correct net radiation
				leaf_transp                 = Penman_Monteith ( meteo_daily, rv, rh, net_rad );
				s->value[CANOPY_TRANSP_SUN] = leaf_transp *  ( transp_daylength_sec * s->value[LAI_SUN_PROJ] ); //* s->value[DAILY_CANOPY_COVER_PROJ] );
                //s->value[CANOPY_TRANSP_SUN] = leaf_transp *  ( transp_daylength_sec * s->value[LAI_SUN_PROJ] )* s->value[DAILY_CANOPY_COVER_PROJ] ;
				
				/** LAI SHADE **/

				/* resistance for shaded canopy fraction */
				rv = 1. / s->value[LEAF_SHADE_CONDUCTANCE];

				/* note: Net Rad is Short wave flux in Biome-bgc*/
				/* convert radiation to stomatal scale */
				// FIXME: compute correctly the net_radiation
#if 0
				net_rad = ( s->value[SW_RAD_ABS_SHADE]  / s->value[LAI_SHADE_PROJ] ) * s->value[F_LIGHT_SHADE_MAKELA];
#else
				net_rad = ( s->value[NET_RAD_ABS_SHADE] / s->value[LAI_SHADE_PROJ] ) * s->value[F_LIGHT_SHADE_MAKELA];
#endif

				/* call Penman-Monteith function, returns e in kg/m2/s for transpiration and W/m2 for latent heat*/
				//fixme use correct net radiation
				leaf_transp                   = Penman_Monteith ( meteo_daily, rv, rh, net_rad );
				s->value[CANOPY_TRANSP_SHADE] = leaf_transp * ( transp_daylength_sec * s->value[LAI_SHADE_PROJ] ); //* s->value[DAILY_CANOPY_COVER_PROJ] );
                //s->value[CANOPY_TRANSP_SHADE] = leaf_transp * ( transp_daylength_sec * s->value[LAI_SHADE_PROJ] )* s->value[DAILY_CANOPY_COVER_PROJ] ;

				/************************************************************************************/
				
				
				s->value[CANOPY_TRANSP_SUN] = s->value[CANOPY_TRANSP_SUN] *s->value[m_corr_sun] ;
				s->value[CANOPY_TRANSP_SHADE] = s->value[CANOPY_TRANSP_SHADE] *s->value[m_corr_shade] ;
				
				/* overall canopy */
				s->value[CANOPY_TRANSP] = s->value[CANOPY_TRANSP_SUN] + s->value[CANOPY_TRANSP_SHADE];

				/* check for negative values */
				if( s->value[CANOPY_TRANSP] < 0. ) s->value[CANOPY_TRANSP]  = 0.;

				/* canopy evapotranspiration */
				s->value[CANOPY_EVAPO_TRANSP] = s->value[CANOPY_EVAPO] + s->value[CANOPY_TRANSP];
			}

			/* note special case for deciduous: assuming that all canopy water evaporates during last day of leaffall */
			if ( ( s->value[PHENOLOGY] == 0.1 || s->value[PHENOLOGY] == 0.2 ) &&
					( s->counter[DAYS_LEAFFALL] == s->counter[LEAF_FALL_COUNTER] ) )
			{
				/* last day of leaf fall */
				s->value[CANOPY_EVAPO]       += s->value[CANOPY_WATER];
				s->value[CANOPY_WATER]        = 0.;
				s->value[CANOPY_EVAPO_TRANSP] = s->value[CANOPY_EVAPO] + s->value[CANOPY_TRANSP];
			}

			/* update canopy water pool */
			c->canopy_water_stored           += ( s->value[CANOPY_INT_RAIN] - s->value[CANOPY_EVAPO] );
			/********************************************************************************************************************************/
		}
		/* if canopy has snow */
		else if ( s->value[CANOPY_SNOW] > 0. )
		{
#if 1
			/********************************************************************************************************************************/
			/* CANOPY SNOW SUBLIMATION OR MELTING (Canopy Snow) */

			rv = 1. / gc_e_wv;
			rh = 1. / gc_sh;

			/* radiation */
			/* convert from W/m2 --> KJ/m2/sec */
#if 0                  // Biome-bgc used SW Radiation
			net_rad = s->value[SW_RAD_ABS]  * SNOW_ABS * 0.001 ;
#else
			net_rad = s->value[NET_RAD_ABS] * SNOW_ABS * 0.001 ;
#endif

			/* snow sublimation or melting */
			subl = net_rad / meteo_daily->lh_sub;
			melt = net_rad / meteo_daily->lh_fus;

			/* check for negative values */
			if( subl < 0. ) subl  = 0.;
			if( melt < 0. ) melt  = 0.;

			if ( meteo_daily->tavg < 0. ) subl_melt = subl;
			else                          subl_melt = melt;

			//FIXME
			s->value[CANOPY_EVAPO] = subl_melt;

			/* calculate the time required to sublimate all the canopy snow */
			subl_melt_daylength_sec     = s->value[CANOPY_SNOW] / s->value[CANOPY_EVAPO];

			/* day not long enough to sublimation or melt all int. snow */
			if( subl_melt_daylength_sec > meteo_daily->daylength_sec )
			{
				/* day not long enough to sublimate or melt all snow intercepted */

				/* NO CANOPY TRANSPIRATION (Canopy completely Snowed)!!! */

				++s->counter[DAYS_WITH_CANOPY_SNOW];

				/* adjust day length for transpiration */
				transp_daylength_sec             = 0.;
				s->value[CANOPY_FRAC_DAY_TRANSP] = 0.;

				/* no time left for transpiration */
				s->value[CANOPY_TRANSP] = 0.;

				/* day length limits canopy sublimation */
				//TODO remove in latent heat computation the sublimation since it is computed still here!!!!
				subl_melt             *= meteo_daily->daylength_sec;
				s->value[CANOPY_EVAPO] = subl_melt;

				s->value[CANOPY_EVAPO_TRANSP] = s->value[CANOPY_EVAPO] + s->value[CANOPY_TRANSP];

				/* remove evaporated water from canopy water pool */
				s->value[CANOPY_SNOW] -= s->value[CANOPY_EVAPO];

			}
			/* all intercepted water evaporated and there's time for transpiration */
			else
			{
				/* all intercepted snow evaporated */
				s->counter[DAYS_WITH_CANOPY_SNOW] = 0;

				/* all canopy snow sublimates or melts */
				subl_melt = s->value[CANOPY_SNOW];

				s->value[CANOPY_EVAPO]  = subl_melt;

				/* reset canopy water from canopy water pool */
				s->value[CANOPY_SNOW]  -= s->value[CANOPY_EVAPO];

				/* adjust day length for transpiration */
				transp_daylength_sec      = meteo_daily->daylength_sec - subl_melt_daylength_sec;

				s->value[CANOPY_FRAC_DAY_TRANSP] = subl_melt_daylength_sec / meteo_daily->daylength_sec;

				/* Leaf-Canopy resistance to sensible heat */
				rh = 1. / gl_sh;

				/************************************************************************************/
				/* calculate transpiration using adjusted day length */

				/** LAI SUN **/

				/* resistance for sun canopy fraction */
				rv = 1. / s->value[LEAF_SUN_CONDUCTANCE];

				/* note: Net Rad is Short wave flux in Biome-bgc*/
				/* convert radiation to stomatal scale */
				// FIXME: compute correctly the net_radiation
#if 0
				net_rad = ( s->value[SW_RAD_ABS_SUN]  / s->value[LAI_SUN_PROJ] ) * s->value[F_LIGHT_SUN_MAKELA];
#else
				net_rad = ( s->value[NET_RAD_ABS_SUN] / s->value[LAI_SUN_PROJ] ) * s->value[F_LIGHT_SUN_MAKELA];
#endif

				/* call Penman-Monteith function, returns e in kg/m2/s for transpiration and W/m2 for latent heat */
				//fixme use correct net radiation
				leaf_transp                 = Penman_Monteith ( meteo_daily, rv, rh, net_rad );
				s->value[CANOPY_TRANSP_SUN] = leaf_transp * ( transp_daylength_sec * s->value[LAI_SUN_PROJ] ); //* s->value[DAILY_CANOPY_COVER_PROJ] );
               //s->value[CANOPY_TRANSP_SUN] = leaf_transp * ( transp_daylength_sec * s->value[LAI_SUN_PROJ] )* s->value[DAILY_CANOPY_COVER_PROJ] ;

				/** LAI SHADE **/

				/* resistance for shaded canopy fraction */
				rv = 1. / s->value[LEAF_SHADE_CONDUCTANCE];

				/* note: Net Rad is Short wave flux in Biome-bgc*/
				/* convert radiation to stomatal scale */
				// FIXME: compute correctly the net_radiation
#if 0
				net_rad = ( s->value[SW_RAD_ABS_SHADE]  / s->value[LAI_SHADE_PROJ] ) * s->value[F_LIGHT_SHADE_MAKELA];
#else
				net_rad = ( s->value[NET_RAD_ABS_SHADE] / s->value[LAI_SHADE_PROJ] ) * s->value[F_LIGHT_SHADE_MAKELA];
#endif
				/* call Penman-Monteith function, returns e in kg/m2/s for transpiration and W/m2 for latent heat*/
				//fixme use correct net radiation
				leaf_transp                   = Penman_Monteith ( meteo_daily, rv, rh, net_rad );
				s->value[CANOPY_TRANSP_SHADE] = leaf_transp * ( transp_daylength_sec * s->value[LAI_SHADE_PROJ] ); //* s->value[DAILY_CANOPY_COVER_PROJ] );
                //s->value[CANOPY_TRANSP_SHADE] = leaf_transp * ( transp_daylength_sec * s->value[LAI_SHADE_PROJ])* s->value[DAILY_CANOPY_COVER_PROJ];

				/************************************************************************************/
				
				
				s->value[CANOPY_TRANSP_SUN] = s->value[CANOPY_TRANSP_SUN] *s->value[m_corr_sun] ;
				s->value[CANOPY_TRANSP_SHADE] = s->value[CANOPY_TRANSP_SHADE] *s->value[m_corr_shade] ;
				
				/* overall canopy */
				s->value[CANOPY_TRANSP] = s->value[CANOPY_TRANSP_SUN] + s->value[CANOPY_TRANSP_SHADE];

				/* check for negative values */
				if( s->value[CANOPY_TRANSP] < 0. ) s->value[CANOPY_TRANSP]  = 0.;

				/* canopy sublimation or melting + transpiration */
				s->value[CANOPY_EVAPO_TRANSP] = s->value[CANOPY_EVAPO] + s->value[CANOPY_TRANSP];
			}

			/* note special case for deciduous: assuming that all canopy snow sublimates or melts during last day of leaffall */
			if ( ( s->value[PHENOLOGY] == 0.1 || s->value[PHENOLOGY] == 0.2 ) &&
					( s->counter[DAYS_LEAFFALL] == s->counter[LEAF_FALL_COUNTER] ) )
			{
				/* last day of leaf fall */
				s->value[CANOPY_EVAPO]       += s->value[CANOPY_SNOW];
				s->value[CANOPY_SNOW]         = 0.;
				s->value[CANOPY_EVAPO_TRANSP] = s->value[CANOPY_EVAPO] + s->value[CANOPY_TRANSP];
			}

			/* update canopy snow pool */
			c->canopy_snow_stored                += ( s->value[CANOPY_INT_SNOW] - s->value[CANOPY_EVAPO] );
			/********************************************************************************************************************************/
#endif
		}
		/** if canopy is dry **/
		else
		{
			/********************************************************************************************************************************/
			/* no canopy evaporation or sublimation occurs */
			evapo = 0.;
			subl  = 0.;
			melt  = 0.;

			s->value[CANOPY_EVAPO] = evapo + subl + melt;

			/* Canopy transpiration (Canopy Dry) */

			/* all day transpiration */
			transp_daylength_sec             = meteo_daily->daylength_sec;

			s->value[CANOPY_FRAC_DAY_TRANSP] = 1.;

			/* Leaf-Canopy resistance to sensible heat */
			rh = 1. / gl_sh;

			/************************************************************************************/
			/* calculate transpiration using adjusted day length */
			/* compute only transpiration */

			/** LAI SUN **/
			/* Transpiration for LAI sun */

			/* resistance for sun canopy fraction */
			rv = 1. / s->value[LEAF_SUN_CONDUCTANCE];

			/* note: Net Rad is Short wave flux in Biome-bgc*/
			/* convert radiation to stomatal scale */
			// FIXME: compute correctly the net_radiation
#if 0
			net_rad = ( s->value[SW_RAD_ABS_SUN] / s->value[LAI_SUN_PROJ] ) * s->value[F_LIGHT_SUN_MAKELA];
#else
			net_rad = ( s->value[NET_RAD_ABS_SUN] / s->value[LAI_SUN_PROJ] ) * s->value[F_LIGHT_SUN_MAKELA];
#endif

			/* call Penman-Monteith function, returns e in kg/m2/s for transpiration and W/m2 for latent heat*/
			//fixme use correct net radiation
			leaf_transp                 = Penman_Monteith (meteo_daily, rv, rh, net_rad);
			s->value[CANOPY_TRANSP_SUN] = leaf_transp * ( transp_daylength_sec * s->value[LAI_SUN_PROJ] ); //* s->value[DAILY_CANOPY_COVER_PROJ] );
          // s->value[CANOPY_TRANSP_SUN] = leaf_transp * ( transp_daylength_sec * s->value[LAI_SUN_PROJ] )* s->value[DAILY_CANOPY_COVER_PROJ] ;
			/** LAI SHADE **/
			/* Transpiration for LAI shade */

			/* resistance for shaded canopy fraction */
			rv = 1. / s->value[LEAF_SHADE_CONDUCTANCE];

			/* note: Net Rad is Short wave flux in Biome-bgc*/
			/* convert radiation to stomatal scale */
			// FIXME: compute correctly the net_radiation
#if 0
			net_rad = ( s->value[SW_RAD_ABS_SHADE] / s->value[LAI_SHADE_PROJ] ) * s->value[F_LIGHT_SHADE_MAKELA];
#else
			net_rad = ( s->value[NET_RAD_ABS_SHADE] / s->value[LAI_SHADE_PROJ] ) * s->value[F_LIGHT_SHADE_MAKELA];
#endif

			/* call Penman-Monteith function, returns e in kg/m2/s for transpiration and W/m2 for latent heat*/
			leaf_transp                   = Penman_Monteith (meteo_daily, rv, rh, net_rad);
			s->value[CANOPY_TRANSP_SHADE] = leaf_transp * ( transp_daylength_sec * s->value[LAI_SHADE_PROJ] ); //* s->value[DAILY_CANOPY_COVER_PROJ] );
  //s->value[CANOPY_TRANSP_SHADE] = leaf_transp * ( transp_daylength_sec * s->value[LAI_SHADE_PROJ] )* s->value[DAILY_CANOPY_COVER_PROJ];
			/************************************************************************************/
			
			
		        s->value[CANOPY_TRANSP_SUN] = s->value[CANOPY_TRANSP_SUN] *s->value[m_corr_sun] ;
		        s->value[CANOPY_TRANSP_SHADE] = s->value[CANOPY_TRANSP_SHADE] *s->value[m_corr_shade] ;
				
			/* overall canopy */
			s->value[CANOPY_TRANSP] = s->value[CANOPY_TRANSP_SUN] + s->value[CANOPY_TRANSP_SHADE];

			/* check for negative values */
			if(s->value[CANOPY_TRANSP] < 0.) s->value[CANOPY_TRANSP] = 0.;

			/* canopy evapotranspiration */
			s->value[CANOPY_EVAPO_TRANSP] = s->value[CANOPY_EVAPO] + s->value[CANOPY_TRANSP];
			/********************************************************************************************************************************/
		}
	}
	else
	{
		s->value[CANOPY_INT_RAIN]     = 0.;
		s->value[CANOPY_WATER]        = 0.;
		s->value[CANOPY_INT_SNOW]     = 0.;
		s->value[CANOPY_SNOW]         = 0.;
		s->value[CANOPY_TRANSP]       = 0.;
		s->value[CANOPY_TRANSP_SHADE] = 0.;
		s->value[CANOPY_EVAPO]        = 0.;
		s->value[CANOPY_EVAPO_TRANSP] = 0.;
	}

	/**********************************************************************************/
	/**********************************************************************************/
	/* control if canopy transpiration exceeds available soil water */
	/* in case it happens firstly down-regulate transpiration to soil water modifier
	 * if it is not enough than superimpose to zero */
	// TODO correct for the multiclass case. or set a proper priority/competition rule.
	// NOTE in this case, all the cohort/class uptake from the same layer.
	 
	if ( s->value[CANOPY_TRANSP] > c->asw )
	{
		printf("ATTENTION CANOPY TRANSPIRATION EXCEEDS ASW!!!!!!\n");

		s->value[CANOPY_TRANSP] *= s->value[F_SW];

		if ( s->value[CANOPY_TRANSP] > c->asw )
		{
			s->value[CANOPY_TRANSP] = 0.;
			printf("CANOPY_TRANSP = 0 !!!!!!\n");
		}
		/* check */
		CHECK_CONDITION ( s->value[CANOPY_TRANSP] , > , c->asw );
	}
	/**********************************************************************************/
	/**********************************************************************************/

	/* compute latent heat fluxes for canopy */
	Canopy_latent_heat_fluxes   (s, meteo_daily);

	/* compute sensible heat fluxes for canopy */
	Canopy_sensible_heat_fluxes (c, layer, height, dbh, age, species, meteo_daily);

	s->value[MONTHLY_CANOPY_TRANSP]       += s->value[CANOPY_TRANSP];
	s->value[MONTHLY_CANOPY_TRANSP_SUN]   += s->value[CANOPY_TRANSP_SUN];
	s->value[MONTHLY_CANOPY_TRANSP_SHADE] += s->value[CANOPY_TRANSP_SHADE];
	s->value[MONTHLY_CANOPY_EVAPO_TRANSP] += s->value[CANOPY_EVAPO_TRANSP];
	s->value[YEARLY_CANOPY_TRANSP]        += s->value[CANOPY_TRANSP];
	s->value[YEARLY_CANOPY_TRANSP_SUN]    += s->value[CANOPY_TRANSP_SUN];
	s->value[YEARLY_CANOPY_TRANSP_SHADE]  += s->value[CANOPY_TRANSP_SHADE];
	s->value[YEARLY_CANOPY_EVAPO_TRANSP]  += s->value[CANOPY_EVAPO_TRANSP];

	c->daily_canopy_evapo                 += s->value[CANOPY_EVAPO];
	c->daily_canopy_transp                += s->value[CANOPY_TRANSP];
	c->daily_canopy_et                    += s->value[CANOPY_EVAPO_TRANSP];

}

