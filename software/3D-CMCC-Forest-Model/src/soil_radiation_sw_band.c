/*
 * soil_radiation_sw_band.c
 *
 *  Created on: 10/set/2016
 *      Author: alessio-cmcc
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "math.h"
#include "common.h"
#include "soil_radiation_sw_band.h"
#include "constants.h"
#include "logger.h"

extern logger_t* g_debug_log;

void soil_radiation_sw_band ( cell_t *const c, meteo_daily_t *meteo_daily )
{
        // FIXME compute net radiation correctly
        // so to consider effectively the incoming radiation and the amount
        // that has to be reflected
        
	//fixme move soil and snow albedo into soil.txt file
	const double soil_albedo = 0.15;                                                      /* (ratio) soil albedo without snow (see MAESPA model) */
	const double snow_albedo = 0.65;                                                      /* (ratio) snow albedo as an average between freshly 
	                                                                                          fallen snow and 0.4 for melting snow */
	double Light_refl_rad_soil_frac = 0.;                                                 /* (ratio) fraction of Short Wave radiation reflected 
	                                                                                         from the soil */
	double Light_refl_rad_snow_frac = 0.;                                                 /* (ratio) fraction of Short Wave radiation reflected 
	                                                                                         from the snow */

	/* check parameters */
	assert ( c );

	/* radiation for soil */
	logger (g_debug_log, "**SOIL RADIATION**\n");

	/* fraction of light reflected */

	if ( !c->snow_pack )
	{
		/* fraction of light reflected by the soil */
		Light_refl_rad_soil_frac = soil_albedo;
	}
	else
	{
		/* fraction of light reflected by the snow */
		Light_refl_rad_snow_frac = snow_albedo;
	}
	logger(g_debug_log, "LightReflec_soil = %g %%\n", Light_refl_rad_soil_frac * 100);
	logger(g_debug_log, "LightReflec_snow = %g %%\n", Light_refl_rad_snow_frac * 100);
	logger(g_debug_log, "******************************************************\n");

	/* compute values for soil layer when last height class in cell is processed */
	// the input radiation-values take into account the fraction absorbed by vegetation
	
	/* Par Short Wave radiation and PPFD for the soil */

	logger(g_debug_log,"*incoming light for soil*\n");
	logger(g_debug_log,"incoming PAR for soil         = %g molPAR/m^2/day\n", meteo_daily->par);
	logger(g_debug_log,"incoming Short Wave radiation = %g W/m2\n", meteo_daily->sw_downward_W);
	logger(g_debug_log,"incoming       Net radiation  = %g W/m2\n", meteo_daily->Net_rad_threePG);

	// accumulate par at the soil and compute average value for the summer period )

	
	
         // remove reflected part

	if ( ! c->snow_pack )
	{
		c->par_refl_soil     = meteo_daily->par             * Light_refl_rad_soil_frac;
		c->sw_rad_refl_soil  = meteo_daily->sw_downward_W   * Light_refl_rad_soil_frac;
		c->net_rad_refl_soil = meteo_daily->Net_rad_threePG * Light_refl_rad_soil_frac;
		logger(g_debug_log,"reflected light from the soil\n");
		logger(g_debug_log,"par_refl_soil         = %g molPAR/m^2/day\n", c->par_refl_soil);
		logger(g_debug_log,"sw_rad_for_soil_refl  = %g W/m2\n", c->sw_rad_refl_soil);
		logger(g_debug_log,"net_rad_for_soil_refl = %g W/m2\n", c->net_rad_refl_soil);
	}
	else
	{
		c->par_refl_snow     = meteo_daily->par             * Light_refl_rad_snow_frac;
		c->sw_rad_refl_snow  = meteo_daily->sw_downward_W   * Light_refl_rad_snow_frac;
		c->net_rad_refl_snow = meteo_daily->Net_rad_threePG * Light_refl_rad_snow_frac;
		logger(g_debug_log,"reflected light from the snow\n");
		logger(g_debug_log,"par_refl_snow         = %g molPAR/m^2/day\n", c->par_refl_snow);
		logger(g_debug_log,"sw_rad_for_snow_refl  = %g W/m2\n", c->sw_rad_refl_snow);
		logger(g_debug_log,"net_rad_for_snow_refl = %g W/m2\n", c->net_rad_refl_snow);
	}



	if ( ! c->snow_pack )
	{
		/* compute absorbed part from soil */
		c->apar_soil        = meteo_daily->par             - c->par_refl_soil;
		c->sw_rad_abs_soil  = meteo_daily->sw_downward_W   - c->sw_rad_refl_soil;
		c->net_rad_abs_soil = meteo_daily->Net_rad_threePG - c->net_rad_refl_soil;
		logger(g_debug_log, "absorbed soil radiation *\n");
	        logger(g_debug_log, "soil aPAR                            = %g molPAR/m^2/day\n", c->apar_soil );
	        logger(g_debug_log, "soil absorbed Short Wave radiation   = %g W/m2\n", c->sw_rad_abs_soil);
	        logger(g_debug_log, "soil absorbed Net radiation for soil = %g W/m2\n", c->net_rad_abs_soil);
	}
	else
	{
		/* compute absorbed part from snow */
		c->apar_snow        = meteo_daily->par             - c->par_refl_snow;
		c->sw_rad_abs_snow  = meteo_daily->sw_downward_W   - c->sw_rad_refl_snow;
		c->net_rad_abs_snow = meteo_daily->Net_rad_threePG - c->net_rad_refl_snow;
	}
	

}
