/*
 * cell_model_check_balance.c
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

/* note: when model checks fluxes there's no need to use global variables */
/* note: when model checks for stocks and it uses at least one stock for balance or it
 *  simulates more than one cell it needs to use global variables */

int check_cell_radiation_flux_balance (cell_t *const c, const meteo_daily_t *const meteo_daily)
{
	double in;
	double out;
	double store;
	double balance;

	logger(g_debug_log,"\n*********CHECK CELL RADIATIVE BALANCE************\n");
	/* DAILY CHECK ON RADIATIVE BALANCE */

	/* PAR RADIATIVE BALANCE */

	/* sum of sources */
	in = meteo_daily->incoming_par;

	/* sum of sinks */
	out = c->par_refl + c->par_refl_soil + c->par_refl_snow;

	/* sum of current storage */
	store = c->apar + c->apar_soil + c->apar_snow ;

	balance = in - out -store;

	logger(g_debug_log, "\nCELL RADIATIVE BALANCE (PAR)\n");

	if ( fabs( balance ) > eps )
	{
		error_log("DOY = %d\n", c->doy);
		error_log("\nin\n");
		error_log("incoming par     = %f molPAR/m2/day\n", meteo_daily->incoming_par);
		error_log("\nout\n");
		error_log("c->par_refl      = %f molPAR/m2/day\n", c->par_refl);
		error_log("c->par_refl_soil = %f molPAR/m2/day\n", c->par_refl_soil);
		error_log("c->par_refl_snow = %f molPAR/m2/day\n", c->par_refl_snow);
		error_log("\nstore\n");
		error_log("c->apar          = %f molPAR/m2/day\n", c->apar);
		error_log("c->par_for_soil  = %f molPAR/m2/day\n", c->apar_soil);
		error_log("c->par_for_snow  = %f molPAR/m2/day\n", c->apar_snow);
		error_log("\npar in         = %f molPAR/m2/day\n", in);
		error_log("par out          = %f molPAR/m2/day\n", out);
		error_log("par store        = %f molPAR/m2/day\n", store);
		error_log("par balance      = %f molPAR/m2/day\n",balance);
		error_log("...FATAL ERROR IN PAR radiative balance (exit)\n");
		CHECK_CONDITION(fabs( balance ), > , eps);

		return 0;
	}
	else
	{
		logger(g_debug_log, "...ok PAR radiative balance\n");
	}

	/*******************************************************************************************************************************************/

	/* SHORT WAVE RADIATION RADIATIVE BALANCE */

	/* sum of sources */
	in = meteo_daily->incoming_sw_downward_W;

	/* sum of sinks */
	out = c->sw_rad_refl + c->sw_rad_refl_soil + c->sw_rad_refl_snow;

	/* sum of current storage */
	store = c->sw_rad_abs + c->sw_rad_abs_soil + c->sw_rad_abs_snow;

	balance = in - out -store;

	logger(g_debug_log, "\nCELL RADIATIVE BALANCE (Short Wave Radiation)\n");

	if ( fabs( balance ) > eps )
	{
		error_log("DOY = %d\n", c->doy);
		error_log("\nin\n");
		error_log("incoming radiation      = %f W/m2\n", meteo_daily->incoming_sw_downward_W);
		error_log("\nout\n");
		error_log("c->sw_rad_refl          = %f W/m2\n", c->sw_rad_refl);
		error_log("c->sw_rad_for_soil_refl = %f W/m2\n", c->sw_rad_refl_soil);
		error_log("c->sw_rad_for_soil_refl = %f W/m2\n", c->sw_rad_refl_snow);
		error_log("\nstore\n");
		error_log("c->sw_rad_abs           = %f W/m2\n", c->sw_rad_abs);
		error_log("c->sw_rad_for_soil  = %f W/m2\n", c->sw_rad_abs_soil);
		error_log("\nradiation in          = %f W/m2\n", in);
		error_log("radiation out           = %f W/m2\n", out);
		error_log("net radiation store     = %f W/m2\n", store);
		error_log("radiation balance       = %f W/m2\n", balance);
		error_log("...FATAL ERROR IN Short Wave Radiation radiative balance (exit)\n");
		CHECK_CONDITION(fabs( balance ), > , eps);

		return 0;
	}
	else
	{
		logger(g_debug_log, "...ok Short Wave Radiation radiative balance\n");
	}

	/* ok */
	return 1;
}

int check_cell_carbon_flux_balance(cell_t *const c)
{
	double in;
	double out;
	double store;
	double balance;

	logger(g_debug_log, "\n*********CHECK CELL CARBON BALANCE************\n");

	/* check complete cell level carbon flux balance */

	/* sum of carbon sources */
	in      = c->daily_gpp;

	/* sum of carbon sinks */
	out     = c->daily_aut_resp + c->daily_het_resp;

	/* sum of current carbon storage */
	store   = - c->daily_nee;

	balance = in - out - store;

	logger(g_debug_log, "\nCELL CARBON FLUX BALANCE\n");

	if ( fabs( balance ) > eps )
	{
		error_log("DOY = %d\n", c->doy);
		error_log("\nin\n");
		error_log("daily_gpp             = %f gC/m2/day\n", c->daily_gpp);
		error_log("\nout\n");
		error_log("daily_tot_aut_resp    = %f gC/m2/day\n", c->daily_aut_resp);
		error_log("daily_tot_het_resp    = %f gC/m2/day\n", c->daily_het_resp);
		error_log("\nstore\n");
		error_log("daily_nee             = %f gC/m2/day\n", c->daily_nee);
		error_log("\ncarbon in           = %f gC/m2/day\n", in);
		error_log("carbon out            = %f gC/m2/day\n", out);
		error_log("carbon store          = %f gC/m2/day\n", store);
		error_log("carbon_balance        = %f gC/m2/day\n", balance);
		error_log("...FATAL ERROR in 'Cell_model_daily' carbon flux balance (first) (exit)\n");
		CHECK_CONDITION(fabs( balance ), > , eps);

		return 0;
	}
	else
	{
		logger(g_debug_log, "...ok 'Cell_model_daily' carbon balance\n");
	}

	/***************************************************************************************************************/
//
//	/* sum of carbon sources */
//	in = c->daily_gpp;
//
//	/* sum of carbon sinks */
//	out = c->daily_maint_resp + c->daily_growth_resp  + c->daily_het_resp;
//
//	/* sum of current carbon storage */
//	store = c->daily_leaf_carbon   + c->daily_stem_carbon +
//			c->daily_froot_carbon  + c->daily_croot_carbon +
//			c->daily_branch_carbon + c->daily_reserve_carbon +
//			c->daily_fruit_carbon  + c->daily_to_litrC +
//			c->daily_to_soilC;
//
//	balance = in - out - store;
//
//	if ( fabs( balance ) > eps )
//	{
//		logger(g_debug_log, "\nCELL CARBON FLUX BALANCE\n");
//
//		error_log("DOY = %d\n", c->doy);
//		error_log("\nin\n");
//		error_log("c->daily_gpp                = %f gC/m2/day\n", c->daily_gpp);
//		error_log("\nout\n");
//		error_log("c->daily_tot_aut_resp       = %f gC/m2/day\n", c->daily_maint_resp + c->daily_growth_resp);
//		error_log("c->daily_tot_het_resp       = %f gC/m2/day\n", c->daily_het_resp);
//		error_log("\nstore\n");
//		error_log("c->daily_leaf_carbon        = %f gC/m2/day\n", c->daily_leaf_carbon);
//		error_log("c->daily_stem_carbon        = %f gC/m2/day\n", c->daily_stem_carbon);
//		error_log("c->daily_fine_root_carbon   = %f gC/m2/day\n", c->daily_froot_carbon);
//		error_log("c->daily_coarse_root_carbon = %f gC/m2/day\n", c->daily_croot_carbon);
//		error_log("c->daily_branch_carbon      = %f gC/m2/day\n", c->daily_branch_carbon);
//		error_log("c->daily_reserve_carbon     = %f gC/m2/day\n", c->daily_reserve_carbon);
//		error_log("c->daily_litrC              = %f gC/m2/day\n", c->daily_to_litrC);
//		error_log("c->daily_soilC              = %f gC/m2/day\n", c->daily_to_soilC);
//		error_log("\ncarbon in                 = %f gC/m2/day\n", in);
//		error_log("carbon out                  = %f gC/m2/day\n", out);
//		error_log("carbon store                = %f gC/m2/day\n", store);
//		error_log("carbon_balance              = %f gC/m2/day\n",balance);
//		error_log("...FATAL ERROR in 'Cell_model_daily' carbon flux balance (second) (exit)\n");
//		CHECK_CONDITION(fabs( balance ), > , eps);
//
//		return 0;
//	}
//	else
//	{
//		logger(g_debug_log, "...ok carbon balance\n");
//	}
	/* ok */
	return 1;
}

int check_cell_carbon_mass_balance(cell_t *const c)
{
	static double leaf_carbon ;
	static double froot_carbon   ;
	static double croot_carbon   ;
	static double stem_carbon    ;
	static double branch_carbon  ;
	static double reserve_carbon ;
	static double fruit_carbon   ;
	static double litrC          ;
	static double soilC          ;

	/* DAILY CHECK ON CELL LEVEL CARBON MASS BALANCE */
	/* check complete cell level carbon mass balance */

	/* sum of sources */
	c->cell_carbon_in    = c->daily_gpp;

	/* sum of sinks */
	c->cell_carbon_out   = c->daily_aut_resp + c->daily_het_resp;

	/* sum of current storage */
	c->cell_carbon_store = c->leaf_carbon +
			c->froot_carbon   +
			c->croot_carbon   +
			c->stem_carbon    +
			c->branch_carbon  +
			c->reserve_carbon +
			c->fruit_carbon   +
			c->litrC          +
			c->soilC          ;

	/* check carbon pool balance */
	c->cell_carbon_balance = c->cell_carbon_in - c->cell_carbon_out - ( c->cell_carbon_store - c->cell_carbon_old_store );

	/* check for carbon mass balance closure */
	if ( ( fabs( c->cell_carbon_balance ) > 1e3 ) && ( c->dos > 1 ) )
	{
		logger(g_debug_log, "\nCELL LEVEL CARBON MASS BALANCE\n");

		error_log("DOS = %d\n", c->dos);
		error_log("\nin\n");
		error_log("daily_gpp_carbon   = %f gC/m2/day\n",     c->daily_gpp);
		error_log("\nout\n");
		error_log("daily_aut_resp     = %f gC/m2/day\n", c->daily_aut_resp);
		error_log("daily_het_resp     = %f gC/m2/day\n", c->daily_het_resp);
		error_log("\nstore\n");
		error_log("leaf_carbon        = %f gC/m2\n", c->leaf_carbon    - leaf_carbon);
		error_log("froot_carbon       = %f gC/m2\n", c->froot_carbon   - froot_carbon);
		error_log("croot_carbon       = %f gC/m2\n", c->croot_carbon   - croot_carbon);
		error_log("stem_carbon        = %f gC/m2\n", c->stem_carbon    - stem_carbon);
		error_log("branch_carbon      = %f gC/m2\n", c->branch_carbon  - branch_carbon);
		error_log("reserve_carbon     = %f gC/m2\n", c->reserve_carbon - reserve_carbon);
		error_log("fruit_carbon       = %f gC/m2\n", c->fruit_carbon   - fruit_carbon);
		error_log("litr_carbon        = %f gC/m2\n", c->litrC          - litrC);
		error_log("soil_carbon        = %f gC/m2\n", c->soilC          - soilC);
		error_log("\ncarbon in        = %f gC/m2\n", c->cell_carbon_in);
		error_log("carbon out         = %f gC/m2\n", c->cell_carbon_out);
		error_log("delta carbon store = %f gC/m2\n", c->cell_carbon_store - c->cell_carbon_old_store);
		error_log("carbon_balance     = %f gC/m2\n", c->cell_carbon_balance);
		error_log("...FATAL ERROR in 'Cell_model_daily' carbon mass balance (gC/m2/day) (exit)\n");
		CHECK_CONDITION(fabs( c->cell_carbon_balance ), > , eps);

		return 0;
	}
	else
	{
		 leaf_carbon    = c->leaf_carbon;
		 froot_carbon   = c->froot_carbon;
		 croot_carbon   = c->croot_carbon;
		 stem_carbon    = c->stem_carbon;
		 branch_carbon  = c->branch_carbon;
		 reserve_carbon = c->reserve_carbon;
		 fruit_carbon   = c->fruit_carbon;
		 litrC          = c->litrC;
		 soilC          = c->soilC;
		c->cell_carbon_old_store = c->cell_carbon_store;
		logger(g_debug_log, "...ok in 'Cell_model_daily' carbon mass balance (gC/m2/day)\n");
	}
	/* ok */
	return 1;
}

int check_cell_nitrogen_flux_balance(cell_t *const c, const meteo_daily_t *const meteo_daily)
{
	//todo as for carbon
	double in;
	double out;
	double store;
	double balance;

	logger(g_debug_log, "\n*********CHECK CELL NITROGEN BALANCE************\n");

	/* check complete cell level nitrogen flux balance */

	/* sum of sources */
	in    = meteo_daily->Ndeposition + c->Nfix;

	/* sum of sinks */
	out   = c->Nleach + c->Nvol;

	/* sum of current nitrogen storage */
	store = c->daily_to_litrN + c->daily_to_soilN;

	balance = in - out - store;

	if ( fabs( balance ) > eps )
	{
		logger(g_debug_log, "\nCELL NITROGEN FLUX BALANCE\n");

		error_log("DOS = %d\n", c->dos);
		error_log("\nin\n");
		error_log("meteo_daily->Ndep = %f gN/m2\n",     meteo_daily->Ndeposition);
		error_log("c->Nfix           = %f gN/m2\n",     c->Nfix);
		error_log("\nout\n");
		error_log("Nleach            = %f gN/m2/day\n", c->Nleach);
		error_log("Nvol              = %f gN/m2/day\n", c->Nvol);
		error_log("\nstore\n");
		error_log("c->daily_to_litrN = %f gC/m2/day\n", c->daily_to_litrN);
		error_log("c->daily_to_soilN = %f gC/m2/day\n", c->daily_to_soilN);
		error_log("\ncarbon in       = %f gC/m2/day\n", in);
		error_log("carbon out        = %f gC/m2/day\n", out);
		error_log("carbon store      = %f gC/m2/day\n", store);
		error_log("carbon_balance    = %f gC/m2/day\n", balance);
		error_log("...FATAL ERROR in 'Cell_model_daily' nitrogen flux balance (first) (exit)\n");
		CHECK_CONDITION(fabs( balance ), > , eps);

		return 0;
	}
	else
	{
		logger(g_debug_log, "...ok 'Cell_model_daily' nitrogen balance\n");
	}
//
//	/***************************************************************************************************************/
//
//	/* sum of carbon sources */
//	in = c->daily_gpp;
//
//	/* sum of carbon sinks */
//	out = c->daily_maint_resp + c->daily_growth_resp;
//
//	/* sum of current carbon storage */
//	store = c->daily_leaf_carbon   + c->daily_stem_carbon +
//			c->daily_froot_carbon  + c->daily_croot_carbon +
//			c->daily_branch_carbon + c->daily_reserve_carbon +
//			c->daily_fruit_carbon  + c->daily_litrC +
//			c->daily_soilC         + c->daily_cwdC;
//
//	balance = in - out - store;
//
//	logger(g_debug_log, "\nCELL CARBON FLUX BALANCE\n");
//
//	if ( fabs( balance ) > eps )
//	{
//		error_log("DOY = %d\n", c->doy);
//		error_log("\nin\n");
//		error_log("c->daily_gpp = %f gC/m2/day\n", c->daily_gpp);
//		error_log("\nout\n");
//		error_log("c->daily_tot_aut_resp = %f gC/m2/day\n",c->daily_maint_resp + c->daily_growth_resp);
//		error_log("\nstore\n");
//		error_log("c->daily_leaf_carbon = %f gC/m2/day\n",        c->daily_leaf_carbon);
//		error_log("c->daily_stem_carbon = %f gC/m2/day\n",        c->daily_stem_carbon);
//		error_log("c->daily_fine_root_carbon = %f gC/m2/day\n",   c->daily_froot_carbon);
//		error_log("c->daily_coarse_root_carbon = %f gC/m2/day\n", c->daily_croot_carbon);
//		error_log("c->daily_branch_carbon = %f gC/m2/day\n",      c->daily_branch_carbon);
//		error_log("c->daily_reserve_carbon = %f gC/m2/day\n",     c->daily_reserve_carbon);
//		error_log("c->daily_litrC = %f gC/m2/day\n",              c->daily_litrC);
//		error_log("c->daily_cwdC = %f gC/m2/day\n",               c->daily_cwdC);
//		error_log("c->daily_soilC = %f gC/m2/day\n",              c->daily_soilC);
//		error_log("\ncarbon in = %f gC/m2/day\n", in);
//		error_log("carbon out = %f gC/m2/day\n", out);
//		error_log("carbon store = %f gC/m2/day\n", store);
//		error_log("carbon_balance = %f gC/m2/day\n",balance);
//		error_log("...FATAL ERROR in 'Cell_model_daily' carbon flux balance (second) (exit)\n");
//		CHECK_CONDITION(fabs( balance ), > , eps);
//
//		return 0;
//	}
//	else
//	{
//		logger(g_debug_log, "...ok carbon balance\n");
//	}
	/* ok */
	return 1;
}

int check_cell_nitrogen_mass_balance(cell_t *const c, const meteo_daily_t *const meteo_daily)
{
	static double leaf;
	static double froot;
	static double croot;
	static double branch;
	static double reserve;
	static double stem;
	static double fruit;
	static double litr;
	static double cwd;
	static double soil;

	/* DAILY CHECK ON CELL LEVEL NITROGEN MASS BALANCE */
	/* check complete cell level nitrogen mass balance */

	/* sum of sources */
	c->cell_nitrogen_in    = meteo_daily->Ndeposition + c->Nfix;

	/* sum of sinks */
	c->cell_nitrogen_out   = c->Nleach + c->Nvol;

	/* sum of current storage */
	c->cell_nitrogen_store = c->leaf_nitrogen +
			c->froot_nitrogen   +
			c->croot_nitrogen   +
			c->stem_nitrogen    +
			c->branch_nitrogen  +
			c->reserve_nitrogen +
			c->fruit_nitrogen   +
			c->litrN            +
			c->soilN            ;

	/* check nitrogen pool balance */
	c->cell_nitrogen_balance = c->cell_nitrogen_in - c->cell_nitrogen_out - ( c->cell_nitrogen_store - c->cell_nitrogen_old_store );

	logger(g_debug_log, "\nCELL LEVEL NITROGEN MASS BALANCE\n");

	/* check for nitrogen mass balance closure */
	if ( ( fabs( c->cell_nitrogen_balance ) > eps ) && ( c->dos > 1 ) )
	{
		error_log("DOS = %d\n", c->dos);
		error_log("\nin\n");
		error_log("Ndep = %f gN/m2\n",                     meteo_daily->Ndeposition);
		error_log("c->Nfix = %f gN/m2\n",                  c->Nfix);
		error_log("\nout\n");
		error_log("Nleach = %f gN/m2/day\n",               c->Nleach);
		error_log("Nvol = %f gN/m2/day\n",                 c->Nvol);
		error_log("\nstore\n");
		error_log("leaf_nitrogen = %f gN/m2/day\n",        leaf - c->leaf_nitrogen);
		error_log("froot_nitrogen = %f gN/m2/day\n",       froot - c->froot_nitrogen);
		error_log("croot_nitrogen = %f gN/m2/day\n",       croot - c->croot_nitrogen);
		error_log("stem_nitrogen = %f gN/m2/day\n",        stem - c->stem_nitrogen);
		error_log("branch_nitrogen = %f gN/m2/day\n",      branch - c->branch_nitrogen);
		error_log("reserve_nitrogen = %f gN/m2/day\n",     reserve - c->reserve_nitrogen);
		error_log("fruit_nitrogen = %f gN/m2/day\n",       fruit - c->fruit_nitrogen);
		error_log("litr_nitrogen = %f gN/m2/day\n",        litr - c->litrN);
		error_log("soil_nitrogen = %f gN/m2/day\n",        soil - c->soilN);
		error_log("\nnitrogen in = %f gN/m2/day\n",        c->cell_nitrogen_in);
		error_log("nitrogen out = %f gN/m2/day\n",         c->cell_nitrogen_out);
		error_log("delta nitrogen store = %f gN/m2/day\n", c->cell_nitrogen_store - c->cell_nitrogen_old_store);
		error_log("nitrogen_balance = %f gN/m2/day\n",     c->cell_nitrogen_balance);
		error_log("...FATAL ERROR in 'Cell_model_daily' nitrogen mass balance (gN/m2/day) (exit)\n");
		CHECK_CONDITION(fabs( c->cell_nitrogen_balance ), > , eps);

		return 0;
	}
	else
	{
		leaf    = c->leaf_nitrogen;
		froot   = c->froot_nitrogen;
		croot   = c->croot_nitrogen;
		branch  = c->branch_nitrogen;
		reserve = c->reserve_nitrogen;
		stem    = c->stem_nitrogen;
		fruit   = c->fruit_nitrogen;
		litr    = c->litrC;
		soil    = c->soilC;

		c->cell_nitrogen_old_store = c->cell_nitrogen_store;
		logger(g_debug_log, "...ok in 'Cell_model_daily' nitrogen mass balance (gC/m2/day)\n");
	}
	/* ok */
	return 1;
}

int check_cell_water_flux_balance(cell_t *const c, const meteo_daily_t *const meteo_daily)
{

	logger(g_debug_log, "\n*********CHECK CELL SOIL WATER BALANCE************\n");
	/* DAILY CHECK ON SOIL POOL-ATMOSPHERE WATER BALANCE */
	/* it takes into account soil-atmosphere fluxes */

	/* sum of sources (rain + snow) */
	c->cell_water_in      = meteo_daily->prcp;

	/* sum of sinks */
	c->cell_water_out     = c->daily_canopy_transp + c->daily_canopy_evapo + c->daily_soil_evapo + c->daily_snow_subl + c->daily_out_flow;

	/* sum of current storage in cell pools */
	c->cell_water_store   = ( c->asw + c->canopy_water_stored + c->canopy_snow_stored + c->snow_pack );

	/* check soil pool water balance */
	c->cell_water_balance = c->cell_water_in - c->cell_water_out - ( c->cell_water_store - c->cell_water_old_store );

	logger(g_debug_log, "\nCELL SOIL POOL WATER BALANCE\n");

	/* check for soil water pool water balance */
	if ( ( fabs( c->cell_water_balance ) > eps ) && ( c->dos > 1 ) )
	{
		error_log("DOY = %d\n", c->dos);
		error_log("\nin\n");
		error_log("meteo_daily->prcp      = %f\n", meteo_daily->prcp);
		error_log("\nout\n");
		error_log("c->daily_canopy_transp = %f\n", c->daily_canopy_transp);
		error_log("c->daily_canopy_evapo  = %f\n", c->daily_canopy_evapo);
		error_log("c->soil_evaporation    = %f\n", c->daily_soil_evapo);
		error_log("c->daily_snow_subl     = %f\n", c->daily_snow_subl);
		error_log("c->out_flow            = %f\n", c->daily_out_flow);
		error_log("\nstore\n");
		error_log("c->canopy_water_stored = %f\n", c->canopy_water_stored);
		error_log("c->canopy_snow_stored  = %f\n", c->canopy_snow_stored);
		error_log("c->asw                 = %f\n", c->asw);
		error_log("c->snow_pack           = %f\n", c->snow_pack);
		error_log("soil water in          = %f\n", c->cell_water_in);
		error_log("soil water out         = %f\n", c->cell_water_out  );
		error_log("soil water store       = %f\n", c->cell_water_store);
		error_log("delta soil water store = %f\n", c->cell_water_store - c->cell_water_old_store);
		error_log("soil water balance     = %f\n", c->cell_water_balance);
		error_log("...FATAL ERROR IN 'Cell_model_daily' soil water balance (exit)\n");
		CHECK_CONDITION(fabs( c->cell_water_balance ), > , eps);

		return 0;
	}
	else
	{
		c->cell_water_old_store = c->cell_water_store;
		logger(g_debug_log, "...ok 'Cell_model_daily' soil water balance\n");
	}

	logger(g_debug_log, "\n*********CHECK CELL SNOW WATER BALANCE************\n");
	/* DAILY CHECK ON SNOW POOL-SOIL POOL-ATMOSPHERE WATER BALANCE */
	/* it takes into account snow-soil-atmosphere fluxes */

	/* sum of sources (snow) */
	c->cell_snow_in      = meteo_daily->snow;

	/* sum of sinks */
	c->cell_snow_out     = c->daily_snow_melt + c->daily_snow_subl;

	/* sum of current storage in snow pack */
	c->cell_snow_store   = c->snow_pack;

	/* check snow pool water balance */
	c->cell_snow_balance = c->cell_snow_in - c->cell_snow_out - ( c->cell_snow_store - c->cell_snow_old_store);

	logger(g_debug_log, "\nCELL SNOW POOL WATER BALANCE\n");

	/* check for soil water pool water balance */
	if ( ( fabs( c->cell_snow_balance ) > eps ) && ( c->dos > 1 ) )
	{
		error_log("DOY = %d\n", c->doy);
		error_log("\nin\n");
		error_log("meteo_daily->snow  = %f\n", meteo_daily->snow);
		error_log("\nout\n");
		error_log("c->daily_snow_subl = %f\n", c->daily_snow_subl);
		error_log("c->daily_snow_melt = %f\n", c->daily_snow_melt);
		error_log("\nstore (as a difference between old and current)\n");
		error_log("delta c->asw       = %f\n", c->snow_pack);
		error_log("soil water in      = %f\n", c->cell_snow_in);
		error_log("soil water out     = %f\n", c->cell_snow_out);
		error_log("soil water store   = %f\n", c->cell_snow_store);
		error_log("soil water balance = %f\n", c->cell_snow_balance);
		error_log("...FATAL ERROR IN 'Cell_model_daily' snow water balance (exit)\n");
		CHECK_CONDITION(fabs( c->cell_snow_balance ), > , eps);

		return 0;
	}
	else
	{
		c->cell_snow_old_store = c->cell_snow_store;
		logger(g_debug_log, "...ok 'Cell_model_daily' snow water balance\n");
	}
	/* ok */
	return 1;
}
