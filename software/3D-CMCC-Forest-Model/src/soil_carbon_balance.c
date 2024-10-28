/*
 * soil_carbon_balance.c
 *
 *  Created on: 17 nov 2017
 *      Author: alessio
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "matrix.h"
#include "common.h"
#include "constants.h"
#include "logger.h"
#include "soil_model.h"

extern logger_t* g_debug_log;

void soil_carbon_balance (cell_t *const c)
{

	logger(g_debug_log, "\n**SOIL CARBON BALANCE**\n");

	/*********************************************************************************************************/

	/* update carbon fluxes balance of fast soil microbial recycling pool */
	c->daily_to_soil1C          = c->daily_litr1C_to_soil1C
			- ( c-> daily_soil1_het_resp + c->daily_soil1C_to_soil2C );

	/* update carbon mass of fast soil microbial pool */
	c->soil1C                  += c->daily_to_soil1C;
	/*********************************************************************************************************/

	/* update carbon fluxes balance of medium soil microbial recycling pool */
	c->daily_to_soil2C          = ( c->daily_litr2C_to_soil2C + c->daily_soil1C_to_soil2C )
			- ( c-> daily_soil2_het_resp + c->daily_soil2C_to_soil3C );

	/* update carbon mass of medium soil microbial pool */
	c->soil2C                  += c->daily_to_soil2C;
	/*********************************************************************************************************/

	/* update carbon fluxes balance of slow soil microbial recycling pool */
	c->daily_to_soil3C          = ( c->daily_litr4C_to_soil3C + c->daily_soil2C_to_soil3C )
			- ( c-> daily_soil3_het_resp + c->daily_soil3C_to_soil4C );

	/* update carbon mass of slow soil microbial pool */
	c->soil3C                  += c->daily_to_soil3C;
	/*********************************************************************************************************/

	/* update carbon fluxes balance of recalcitrant SOM pool */
	c->daily_to_soil4C          = ( c->daily_soil3C_to_soil4C )
			- ( c-> daily_soil4_het_resp );

	/* update carbon mass of recalcitrant SOM pool */
	c->soil4C                  += c->daily_to_soil4C;
	/*********************************************************************************************************/

	c->daily_to_soilC           = c->daily_litr1C_to_soil1C +
			c->daily_litr2C_to_soil2C +
			c->daily_litr4C_to_soil3C ;

	/* total mass */
	c->soilC                    = c->soil1C + c->soil2C + c->soil3C + c->soil4C;

	 //printf("LITR TO SOIL           = %g tC/cell/day\n", c->daily_litr1C_to_soil1C + c->daily_litr2C_to_soil2C + c->daily_litr4C_to_soil3C );


}
