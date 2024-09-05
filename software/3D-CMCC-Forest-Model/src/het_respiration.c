/*
 * het_respiration.c
 *
 *  Created on: 15 mar 2017
 *      Author: alessio
 */

/* includes */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "aut_respiration.h"
#include "het_respiration.h"
#include "common.h"
#include "constants.h"
#include "settings.h"
#include "logger.h"

//extern settings_t* g_settings;
//extern logger_t* g_debug_log;

void heterotrophic_respiration ( cell_t *const c, const meteo_daily_t *const meteo_daily )
{

}

void litter_heterotrophic_respiration ( cell_t *const c, const meteo_daily_t *const meteo_daily )
{
	/* summarize litter heterotrophic respiration */
	c->daily_litr_het_resp = c->daily_litr1_het_resp + c->daily_litr2_het_resp + c->daily_litr4_het_resp;

	/* summarize carbon litter to soil */
	c->daily_litr_to_soilC = c->daily_litr1C_to_soil1C + c->daily_litr2C_to_soil2C + c->daily_litr4C_to_soil3C;

	/* summarize nitrogen litter to soil */
	c->daily_litr_to_soilN = c->daily_litr1N_to_soil1N + c->daily_litr2N_to_soil2N + c->daily_litr4N_to_soil3N;

	/*******************************************************************************************************************************/

	/* total heterotrophic respiraton */

	/* daily */
	c->daily_het_resp        += c->daily_litr_het_resp;

	/* monthly */
	c->monthly_het_resp      += c->daily_litr_het_resp;

	/* annual */
	c->annual_het_resp       += c->daily_litr_het_resp;

	/*******************************************************************************************************************************/
}

void soil_heterotrophic_respiration ( cell_t *const c, const meteo_daily_t *const meteo_daily )
{
	/* summarize soil heterotrophic respiration */
	c->daily_soil_het_resp = c->daily_soil1_het_resp + c->daily_soil2_het_resp + c->daily_soil3_het_resp + c->daily_soil4_het_resp;

	/*******************************************************************************************************************************/
	
	/* total heterotrophic respiraton */

	/* daily */
	c->daily_het_resp        += c->daily_soil_het_resp;

	/* monthly */
	c->monthly_het_resp      += c->daily_soil_het_resp;

	/* annual */
	c->annual_het_resp       += c->daily_soil_het_resp;

	/*******************************************************************************************************************************/
}

