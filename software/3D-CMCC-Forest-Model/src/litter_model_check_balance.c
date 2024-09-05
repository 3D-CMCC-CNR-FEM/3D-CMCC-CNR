/*
 * litter_model_check_balance.c
 *
 *  Created on: 03 apr 2017
 *      Author: alessio
 */

/* includes */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "common.h"
#include "check_balance.h"
#include "constants.h"
#include "settings.h"
#include "logger.h"

extern logger_t* g_debug_log;

int check_litter_carbon_flux_balance(cell_t *const c)
{
	double in;
	double out;
	double store;
	double balance;

	logger(g_debug_log, "\n*********CHECK LITTER CARBON BALANCE************\n");

	/* check complete litter and cwd level carbon flux balance */

	/* sum of carbon sources */
	in      = c->daily_leaf_to_litrC + c->daily_froot_to_litrC + c->daily_cwd_to_litrC + c->daily_reserve_to_litrC;

	/* sum of carbon sinks */
	out     = c->daily_litr_het_resp + c->daily_litr_to_soilC;

	/* sum of current carbon storage */
	store   = c->daily_to_litr1C +
			c->daily_to_litr2C   +
			c->daily_to_litr3C   +
			c->daily_to_litr4C;

	balance = in - out - store;

	if ( ( fabs( balance ) > eps ) && ( c->dos > 1) )
	{
		error_log("\nLITTER CARBON FLUX BALANCE\n");

		error_log("DOY                       = %d\n", c->doy);
		error_log("\nin\n");
		error_log("c->daily_leaf_to_litrC    = %f gC/m2/day\n", c->daily_leaf_to_litrC);
		error_log("c->daily_froot_to_litrC   = %f gC/m2/day\n", c->daily_froot_to_litrC);
		error_log("c->daily_cwd_to_litrC     = %f gC/m2/day\n", c->daily_cwd_to_litrC);
		error_log("\nout\n");
		error_log("c->daily_litr_to_soilC    = %f gC/m2/day\n", c->daily_litr_to_soilC);
		error_log("c->daily_litr_het_resp    = %f gC/m2/day\n", c->daily_litr_het_resp);
		error_log("\nstore\n");
		error_log("c->daily_to_litr1C        = %f gC/m2/day\n", c->daily_to_litr1C);
		error_log("c->daily_to_litr2C        = %f gC/m2/day\n", c->daily_to_litr2C);
		error_log("c->daily_to_litr3C        = %f gC/m2/day\n", c->daily_to_litr3C);
		error_log("c->daily_to_litr4C        = %f gC/m2/day\n", c->daily_to_litr4C);
		error_log("\nbalance\n");
		error_log("carbon in                 = %f gC/m2/day\n", in);
		error_log("carbon out                = %f gC/m2/day\n", out);
		error_log("carbon store              = %f gC/m2/day\n", store);
		error_log("carbon_balance            = %f gC/m2/day\n", balance);
		error_log("...FATAL ERROR in 'Litter_ model_daily' carbon balance (exit)\n");
		CHECK_CONDITION(fabs( balance ), > , eps);

		return 0;
	}
	else
	{
		logger(g_debug_log, "...ok 'Litter_ model_daily' carbon balance\n");
	}
	/* ok */
	return 1;
}

int check_litter_carbon_mass_balance(cell_t *const c)
{
	/* DAILY CHECK ON LITTER LEVEL CARBON MASS BALANCE */
	/* check complete litter level carbon mass balance */

	/* sum of sources */
	c->litr_carbon_in    = c->daily_leaf_to_litrC + c->daily_froot_to_litrC + c->daily_cwd_to_litrC + c->daily_reserve_to_litrC;

	/* sum of sinks */
	c->litr_carbon_out   = c->daily_litr_het_resp + c->daily_litr_to_soilC;

	/* sum of current storage */
	c->litr_carbon_store =
			c->litr1C +
			c->litr2C +
			c->litr3C +
			c->litr4C;

	/* check carbon pool balance */
	c->litr_carbon_balance = c->litr_carbon_in - c->litr_carbon_out - ( c->litr_carbon_store - c->litr_carbon_old_store );

	/* check for carbon mass balance closure */
	if ( ( fabs( c->litr_carbon_balance ) > eps ) && ( c->dos > 1 ) )
	{
		error_log("\nLITTER LEVEL CARBON MASS BALANCE\n");

		error_log("DOS                    = %d\n", c->dos);
		error_log("\nin\n");
		error_log("in                     = %f gC/m2/day\n", c->litr_carbon_in);
		error_log("daily_leaf_to_litrC    = %f gC/m2\n",     c->daily_leaf_to_litrC);
		error_log("daily_froot_to_litrC   = %f gC/m2\n",     c->daily_froot_to_litrC);
		error_log("daily_cwd_to_litrC     = %f gC/m2/day\n", c->daily_cwd_to_litrC);
		error_log("\nout\n");
		error_log("out                    = %f gC/m2/day\n", c->litr_carbon_out);
		error_log("daily_het_resp         = %f gC/m2/day\n", c->daily_litr_het_resp);
		error_log("daily_to_soilC         = %f gC/m2/day\n", c->daily_to_soilC);
		error_log("\nstore\n");
		error_log("litr1C                 = %f gC/m2/day\n", c->litr1C);
		error_log("litr2C                 = %f gC/m2/day\n", c->litr2C);
		error_log("litr3C                 = %f gC/m2/day\n", c->litr3C);
		error_log("litr4C                 = %f gC/m2/day\n", c->litr4C);
		error_log("\nbalance\n");
		error_log("carbon in              = %f gC/m2/day\n", c->litr_carbon_in);
		error_log("carbon out             = %f gC/m2/day\n", c->litr_carbon_out);
		error_log("delta carbon store     = %f gC/m2/day\n", c->litr_carbon_store - c->litr_carbon_old_store);
		error_log("carbon_balance         = %f gC/m2/day\n", c->litr_carbon_balance);
		error_log("...FATAL ERROR in 'Litter_model_daily' carbon mass balance (gC/m2/day) (exit)\n");
		CHECK_CONDITION(fabs( c->litr_carbon_balance ), > , eps);

		return 0;
	}
	else
	{
		c->litr_carbon_old_store = c->litr_carbon_store;
		logger(g_debug_log, "...ok in 'Litter_model_daily' carbon mass balance (gC/m2/day)\n");
	}
	/* ok */
	return 1;

}

/******************************************************************************************/
/******************************************************************************************/

int check_litter_nitrogen_flux_balance(cell_t *const c)
{
	double in;
	double out;
	double store;
	double balance;

	logger(g_debug_log, "\n*********CHECK LITTER NITROGEN BALANCE************\n");

	/* check complete litter level nitrogen flux balance */

	/* sum of nitrogen sources */
	in      = c->daily_leaf_to_litrN + c->daily_froot_to_litrN + c->daily_cwd_to_litrN /* fixme include Ndepostiion daily NOT here but in the soil balance */;

	/* sum of nitrogen sinks */
	out     = c->daily_litr_to_soilN /* fixme include N plants demand */;

	/* sum of current nitrogen storage */
	store   = c->daily_to_litr1N +
			c->daily_to_litr2N   +
			c->daily_to_litr3N   +
			c->daily_to_litr4N;

	balance = in - out - store;

	if ( ( fabs( balance ) > eps ) && ( c->dos > 1) )
	{
		error_log("\nLITTER NITROGEN FLUX BALANCE\n");

		error_log("DOY                       = %d\n", c->doy);
		error_log("\nin\n");
		error_log("c->daily_leaf_to_litrN    = %f gN/m2/day\n", c->daily_leaf_to_litrN);
		error_log("c->daily_froot_to_litrN   = %f gN/m2/day\n", c->daily_froot_to_litrN);
		error_log("c->daily_cwd_to_litrN     = %f gN/m2/day\n", c->daily_cwd_to_litrN);
		error_log("\nout\n");
		error_log("c->daily_litr_to_soilN    = %f gN/m2/day\n", c->daily_litr_to_soilN);
		error_log("\nstore\n");
		error_log("c->daily_to_litr1N        = %f gN/m2/day\n", c->daily_to_litr1N);
		error_log("c->daily_to_litr2N        = %f gN/m2/day\n", c->daily_to_litr2N);
		error_log("c->daily_to_litr3N        = %f gN/m2/day\n", c->daily_to_litr3N);
		error_log("c->daily_to_litr4N        = %f gN/m2/day\n", c->daily_to_litr4N);
		error_log("\nbalance\n");
		error_log("nitrogen in               = %f gN/m2/day\n", in);
		error_log("nitrogen out              = %f gN/m2/day\n", out);
		error_log("nitrogen store            = %f gN/m2/day\n", store);
		error_log("nitrogen_balance          = %f gN/m2/day\n", balance);
		error_log("...FATAL ERROR in 'Litter_ model_daily' nitrogen balance (exit)\n");
		CHECK_CONDITION(fabs( balance ), > , eps);

		return 0;
	}
	else
	{
		logger(g_debug_log, "...ok 'Litter_ model_daily' nitrogen balance\n");
	}
	/* ok */
	return 1;
}

int check_litter_nitrogen_mass_balance(cell_t *const c)
{
	/* DAILY CHECK ON LITTER LEVEL NITROGEN MASS BALANCE */
	/* check complete litter level nitrogen mass balance */

	/* sum of sources */
	c->litr_nitrogen_in    = c->daily_leaf_to_litrN + c->daily_froot_to_litrN + c->daily_cwd_to_litrN /* fixme include Ndepostiion daily */;

	//fixme
	/* sum of sinks */
	c->litr_nitrogen_out   = c->daily_litr_to_soilN /* fixme include N plants demand */;

	//todo include cwd nitrogen pool
	/* sum of current storage */
	c->litr_nitrogen_store = c->litr1N +
			c->litr2N  +
			c->litr3N  +
			c->litr4N ;

	/* check nitrogen pool balance */
	c->litr_nitrogen_balance = c->litr_nitrogen_in - c->litr_nitrogen_out - ( c->litr_nitrogen_store - c->litr_nitrogen_old_store );

	/* check for nitrogen mass balance closure */
	if ( ( fabs( c->litr_nitrogen_balance ) > eps ) && ( c->dos > 1 ) )
	{
		error_log("\nLITTER LEVEL NITROGEN MASS BALANCE\n");

		error_log("DOS                       = %d\n", c->dos);
		error_log("\nin                      = %f gN/m2/day\n", c->litr_nitrogen_in);
		error_log("daily_leaf_to_litrN       = %f gN/m2\n",     c->daily_leaf_to_litrN);
		error_log("c->daily_froot_to_litrN   = %f gN/m2/day\n", c->daily_froot_to_litrN);
		error_log("c->daily_cwd_to_litrN     = %f gN/m2/day\n", c->daily_cwd_to_litrN);
		error_log("\nout                     = %f gN/m2/day\n", c->litr_nitrogen_out);
		error_log("daily_litr_to_soilN       = %f gN/m2/day\n", c->daily_litr_to_soilN);
		error_log("\nstore                   = %f gN/m2/day\n", c->litr_nitrogen_store);
		error_log("leaf_litr1N               = %f gN/m2/day\n", c->litr1N);
		error_log("leaf_litr2N               = %f gN/m2/day\n", c->litr2N);
		error_log("leaf_litr3N               = %f gN/m2/day\n", c->litr3N);
		error_log("leaf_litr4N               = %f gN/m2/day\n", c->litr4N);
		error_log("\nbalance                 = %f gN/m2/day\n", c->litr_nitrogen_balance);
		error_log("nitrogen in               = %f gN/m2/day\n", c->litr_nitrogen_in);
		error_log("nitrogen out              = %f gN/m2/day\n", c->litr_nitrogen_out);
		error_log("delta nitrogen store      = %f gN/m2/day\n", c->litr_nitrogen_store - c->litr_nitrogen_old_store);
		error_log("nitrogen_balance          = %f gN/m2/day\n", c->litr_nitrogen_balance);
		error_log("...FATAL ERROR in 'Litter_model_daily' nitrogen mass balance (exit)\n");
		CHECK_CONDITION(fabs( c->litr_nitrogen_balance ), > , eps);

		return 0;
	}
	else
	{
		c->litr_nitrogen_old_store = c->litr_nitrogen_store;
		logger(g_debug_log, "...ok in 'Litter_model_daily' nitrogen mass balance (gN/m2/day)\n");
	}
	/* ok */
	return 1;
}
