/*
 * decomposition.c
 *
 *  Created on: 13 mar 2017
 *      Author: alessio
 */


/* includes */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "common.h"
#include "constants.h"
#include "settings.h"
#include "logger.h"
#include "littering.h"

extern settings_t* g_settings;
//extern soil_settings_t *g_soil_settings;

#define BIOME 0 /* 0 is off and uses CENTURY approach */

double decomposition ( cell_t *const c, const meteo_daily_t *const meteo_daily, int flag )
{
	double temp_scalar;           /* soil temperature scalar */
	double water_scalar;          /* soil water scalar */
	double soil_scalar;           /* soil physic scalar (only for soil layer and when using century method) */
	double rate_scalar;           /* final rate scalar as the product of the temperature and water scalars */
	double min_psi, max_psi;
	double TsoilK;
	double E0;                    /* activation-energy-type parameter of Lloyd and Taylor [1994] (K-1) */
	double T0;                    /* the lower temperature limit for the soil respiration (K) */
	double Tbase;

	/************************************************************************************************/

	/*** decomposition ***/

#if BIOME

	E0    = 308.56;
	T0    = 227.13;
	Tbase = 71.02;   /* for Biome it is 71.02 but for Lloyd and Taylor it was 56.02 */

	/* calculate the rate constant scalar for soil temperature,
		assuming that the base rate constants are assigned for non-moisture
		limiting conditions at 25 C. The function used here is taken from
		Lloyd, J., and J.A. Taylor, 1994. On the temperature dependence of
		soil respiration. Functional Ecology, 8:315-323.
		This equation is a modification of their eqn. 11, changing the base
		temperature from 10 C to 25 C, since most of the microcosm studies
		used to get the base decomp rates were controlled at 25 C.
		note: this implies to change value value from 56.02 to 71.02 K */


	if ( meteo_daily->tsoil < -10. )
	{
		/* no decomposition processes for tsoil < -10.0 C */
		temp_scalar = 0.0;
	}
	else
	{
		TsoilK      = meteo_daily->tsoil + TempAbs;
		temp_scalar = exp ( E0 * ( ( 1. / Tbase ) - ( 1. / ( TsoilK - T0 ) ) ) );
	}
	//fixme bug in original BIOME-BGC CODE
	if ( temp_scalar > 1 ) temp_scalar = 1.;

#else

	/* litter */
	/* following Parton et al., 1987 and Epron et al., 2001 effects of temperature on litter and soil decomposition rate */
	if ( flag == 0 )
	{
		temp_scalar = pow((( 45. - meteo_daily->tsoil ) / 20. ) , 0.2 ) * exp (0.076 * ( 1. - pow((( 45. - meteo_daily->tsoil ) / 20. ) , 4.9 ) ) );
	}
	/* soil */
	else
	{
		temp_scalar = pow((( 45. - meteo_daily->tsoil ) / 10. ) , 0.2 ) * exp (0.076 * ( 1. - pow((( 45. - meteo_daily->tsoil ) / 10. ) , 2.63 ) ) );
	}
	if ( temp_scalar > 1 ) temp_scalar = 1.;

#endif

#if BIOME

	/* calculate the rate constant scalar for soil water content.
		Uses the log relationship with water potential given in
		Andren, O., and K. Paustian, 1987. Barley straw decomposition in the field:
		a comparison of models. Ecology, 68(5):1190-1200.
		and supported by data in
		Orchard, V.A., and F.J. Cook, 1983. Relationship between soil respiration
		and soil moisture. Soil Biol. Biochem., 15(4):447-453.
	 */

	/* set the maximum and minimum values for water potential limits (MPa) */
	min_psi = -10.;
	max_psi = c->psi_sat;

	if (c->psi < min_psi)
	{
		/* no respiration below the minimum soil water potential */
		water_scalar = 0.;
	}
	else if ( c->psi > max_psi )
	{
		/* this shouldn't ever happen, but just in case... */
		water_scalar = 1.;
	}
	else
	{
		water_scalar = log ( min_psi / c->psi ) / log ( min_psi / max_psi );
	}

	soil_scalar = 1.;

#else

	/* following Parton et al., 1987 and Epron et al., 2001 effects of temperature on litter and soil decomposition rate */
	water_scalar = 1. / ( 1. + 30. * exp ( -8.5 * ( c->asw / c->max_asw_fc ) ) );
	if ( water_scalar > 1 ) water_scalar = 1.;


	/* soil physic scalar for soil pool as in Parton et al., 1987 and Epron et al., 2001 */
	/*litter */
	if ( flag == 0 )
	{
		/* no soil physic scalar for litter pool */
		soil_scalar = 1.;
	}
	else
	{
		soil_scalar = ( 1. - 0.75 * ( ( c->clay_perc / 100. ) + ( c->silt_perc / 100. ) ) );
	}
	if ( soil_scalar > 1 ) soil_scalar = 1.;

#endif

	/* calculate the final rate scalar as the product of the temperature and water scalars */
	rate_scalar     = temp_scalar * water_scalar * soil_scalar;

	/* check */
	CHECK_CONDITION (temp_scalar , > , 1);
	CHECK_CONDITION (water_scalar, > , 1);
	CHECK_CONDITION (soil_scalar , > , 1);
	CHECK_CONDITION (rate_scalar , > , 1);

	return (rate_scalar);
}


void litter_decomposition (cell_t *const c, const meteo_daily_t *const meteo_daily)
{
	double rate_scalar;
	double cn_cwd = 0.;
	double cn_cwd2 = 0.;
	double cn_cwd3 = 0.;
	double cn_cwd4 = 0.;
	double cn_litr1 = 0.;
	double cn_litr2 = 0.;
	double cn_litr4 = 0.;

	double deadwood_fragm_rate;   /* deadwood physical fragmentation */
	double litt_decomp_rate1;     /* labile litter decomposition rate */
	double litt_decomp_rate2;     /* cellulose litter decomposition rate */
	double litt_decomp_rate4;     /* lignin litter decomposition rate */
	double pot_litr1C_loss;       /* potential labile litter loss */
	double pot_litr2C_loss;       /* potential unshielded litter loss */
	double pot_litr4C_loss;       /* potential lignin litter loss */

	double pot_min_litr1_soil1;   /* potential mineralization rate */
	double pot_min_litr2_soil2;   /* potential mineralization rate */
	double pot_min_litr4_soil3;   /* potential mineralization rate */
	double pot_min_soil1_soil2;   /* potential mineralization rate */
	double pot_min_soil2_soil3;   /* potential mineralization rate */
	double pot_min_soil3_soil4;   /* potential mineralization rate */
	double pot_min_soil4;         /* potential mineralization rate */
	double pot_immobilization;
	double mineralized;
	double Nlimit;


	/******************************************************************************************************************/

	/* when spinup is on */
	/* recompute deadwood + litter fractions */
	if ( g_settings->spinup ) spinup_littering ( c );

	/******************************************************************************************************************/

	/* calculate the final rate scalar as the product of the temperature and water scalars */
	rate_scalar     = decomposition ( c, meteo_daily, 0 );

	/* calculate C:N ratios */
	if ( c->cwd_N  > 0. ) cn_cwd       = c->cwd_C  / c->cwd_N;
	if ( c->cwd_2N > 0. ) cn_cwd2      = c->cwd_2C / c->cwd_2N;
	if ( c->cwd_3N > 0. ) cn_cwd3      = c->cwd_3C / c->cwd_3N;
	if ( c->cwd_4N > 0. ) cn_cwd4      = c->cwd_4C / c->cwd_4N;
	if ( c->litr1N      > 0. ) cn_litr1     = c->litr1C      / c->litr1N;
	if ( c->litr2N      > 0. ) cn_litr2     = c->litr2C      / c->litr2N;
	if ( c->litr4N      > 0. ) cn_litr4     = c->litr4C      / c->litr4N;

	/******************************************************************************************************************/

	/* calculate the flux from deadwood to litter lignin and cellulose compartments, due to physical fragmentation */
	deadwood_fragm_rate              = KFRAG_BASE     * rate_scalar;

	/* check */
	CHECK_CONDITION ( ( c->cwd_2C + c->cwd_3C + c->cwd_4C ) - c->cwd_C, > , eps);

	/* coarse woody debris carbon to carbon litter poool */
	c->daily_cwd_to_litrC       = c->cwd_C  * deadwood_fragm_rate;
	c->daily_cwd_to_litr2C      = c->cwd_2C * deadwood_fragm_rate;
	c->daily_cwd_to_litr3C      = c->cwd_3C * deadwood_fragm_rate;
	c->daily_cwd_to_litr4C      = c->cwd_4C * deadwood_fragm_rate;

	/* check */
	CHECK_CONDITION ( c->daily_cwd_to_litr2C + c->daily_cwd_to_litr3C + c->daily_cwd_to_litr4C, > , c->daily_cwd_to_litrC + eps);

	/* coarse woody debris nitrogen to nitrogen litter pool */
	if ( cn_cwd  > 0. ) c->daily_cwd_to_litrN  = c->daily_cwd_to_litrC  / cn_cwd;
	if ( cn_cwd2 > 0. ) c->daily_cwd_to_litr2N = c->daily_cwd_to_litr2C / cn_cwd2;
	if ( cn_cwd3 > 0. ) c->daily_cwd_to_litr3N = c->daily_cwd_to_litr3C / cn_cwd3;
	if ( cn_cwd4 > 0. ) c->daily_cwd_to_litr4N = c->daily_cwd_to_litr4C / cn_cwd4;

	/******************************************************************************************************************/

	/* litter decomposition rate */
	litt_decomp_rate1     = KL1_BASE * rate_scalar;
	litt_decomp_rate2     = KL2_BASE * rate_scalar;
	litt_decomp_rate4     = KL4_BASE * rate_scalar;

	/* compute litter potential carbon loss */

	/* 1. labile litter to fast microbial recycling pool */
	if ( c->litr1C > 0. )
	{
		pot_litr1C_loss      = c->litr1C * litt_decomp_rate1;
	}

	/* 2. cellulose litter to medium microbial recycling pool */
	if ( c->litr2C > 0. )
	{
		pot_litr2C_loss      = c->litr2C * litt_decomp_rate2;
	}

	/* 3. lignin litter to slow microbial recycling pool */
	if ( c->litr4C > 0. )
	{
		pot_litr4C_loss      = c->litr4C * litt_decomp_rate4;
	}

	/* calculate the non-nitrogen limited fluxes between litter and
		soil compartments. These will be ammended for N limitation if it turns
		out the potential gross immobilization is greater than potential gross
		mineralization. */

	/* 1 labile litter fluxes */
	if ( c->litr1C > 0. )
	{
		/* carbon */
		c->daily_litr1_het_resp   = pot_litr1C_loss * RFL1S1;
		c->daily_litr1C_to_soil1C = pot_litr1C_loss * ( 1. - RFL1S1 );
		/* nitrogen */
		if ( c->litr1N > 0. ) c->daily_litr1N_to_soil1N = pot_litr1C_loss / cn_litr1;
		else                  c->daily_litr1N_to_soil1N = 0.;
	}

	/* cellulose litter fluxes */
	if ( c->litr2C > 0. )
	{
		/* carbon */
		c->daily_litr2_het_resp   = pot_litr2C_loss * RFL2S2;
		c->daily_litr2C_to_soil2C = pot_litr2C_loss * ( 1. - RFL2S2 );
		/* nitrogen */
		if ( c->litr2N > 0. ) c->daily_litr2N_to_soil2N = pot_litr2C_loss / cn_litr2;
		else                  c->daily_litr2N_to_soil2N = 0.;
	}

	/* release of shielded cellulose litter, tied to the decay rate of lignin litter */
	if ( c->litr3C > 0. )
	{
		//fixme fixme fixme fixme miss competition (see biome)
		c->daily_litr3C_to_litr2C = litt_decomp_rate4 * c->litr3C;
	}

	/* lignin litter fluxes */
	if ( c->litr4C > 0. )
	{
		/* carbon */
		c->daily_litr4_het_resp   = pot_litr4C_loss * RFL4S3;
		c->daily_litr4C_to_soil3C = pot_litr4C_loss * ( 1. - RFL4S3 );
	
		/* nitrogen */
		if ( c->litr4N > 0. ) c->daily_litr4N_to_soil3N = pot_litr4C_loss / cn_litr4;
		else                  c->daily_litr4N_to_soil3N = 0.;
	}
}

void soil_decomposition (cell_t *const c, const meteo_daily_t *const meteo_daily)
{
	double rate_scalar;
	double soil_decomp_rate1;     /* labile soil decomposition rate */
	double soil_decomp_rate2;     /* cellulose unshielded soil decomposition rate */
	double soil_decomp_rate3;     /* cellulose schielded soil decomposition rate */
	double soil_decomp_rate4;     /* lignin soil decomposition rate */
	double pot_soil1C_loss;       /* potential labile soil loss */
	double pot_soil2C_loss;       /* potential unshielded soil loss */
	double pot_soil3C_loss;       /* potential shielded soil loss */
	double pot_soil4C_loss;       /* potential lignin soil loss */

	/* calculate the final rate scalar as the product of the temperature and water scalars */
	rate_scalar     = decomposition ( c, meteo_daily, 1 );

	/* soil decomposition rate */
	soil_decomp_rate1     = KS1_BASE * rate_scalar;
	soil_decomp_rate2     = KS2_BASE * rate_scalar;
	soil_decomp_rate3     = KS3_BASE * rate_scalar;
	soil_decomp_rate4     = KS4_BASE * rate_scalar;

	/* compute soil potential carbon loss */

	/* 1. fast microbial recycling pool to medium microbial recycling pool */
	if ( c->soil1C > 0.)
	{
		pot_soil1C_loss = soil_decomp_rate1 * c->soil1C;
		//pmnf_s1s2       = ( pot_soil1C_loss * (1. - RFS1S2 - ( SOIL2_CN / SOIL1_CN ) ) ) / SOIL2_CN;
	}

	/* 2. medium microbial recycling pool to slow microbial recycling pool */
	if ( c->soil2C > 0. )
	{
		pot_soil2C_loss = soil_decomp_rate2 * c->soil2C;
		//pmnf_s2s3       = (pot_soil2C_loss * ( 1. - RFS2S3 - ( SOIL3_CN / SOIL2_CN ) ) ) / SOIL3_CN;
	}

	/* 3. slow microbial recycling pool to recalcitrant SOM pool */
	if ( c->soil3C > 0. )
	{
		pot_soil3C_loss = soil_decomp_rate3 * c->soil3C;
		//pmnf_s3s4       = (pot_soil3C_loss * ( 1. - RFS3S4 - ( SOIL4_CN / SOIL3_CN ) ) ) / SOIL4_CN;
	}

	/* 4. mineralization of recalcitrant SOM */
	if ( c->soil4C > 0. )
	{
		pot_soil4C_loss = soil_decomp_rate4 * c->soil4C;
		//pmnf_s4 = -pot_soil4C_loss / SOIL4_CN;
	}

	/***********************************************/

	/* fast microbial recycling pool */
	if ( c->soil1C > 0. )
	{
		/* carbon */
		c->daily_soil1_het_resp   = pot_soil1C_loss * RFS1S2;
		c->daily_soil1C_to_soil2C = pot_soil1C_loss * ( 1. - RFS1S2 );
		/* nitrogen */
		c->daily_soil1N_to_soil2N = pot_soil1C_loss / SOIL1_CN;
	}

	/* medium microbial recycling pool */
	if ( c->soil2C > 0. )
	{
		/* carbon */
		c->daily_soil2_het_resp   = pot_soil2C_loss * RFS2S3;
		c->daily_soil2C_to_soil3C = pot_soil2C_loss * ( 1. - RFS2S3 );
		/* nitrogen */
		c->daily_soil2N_to_soil3N = pot_soil2C_loss / SOIL2_CN;
	}

	/* medium microbial recycling pool */
	if ( c->soil3C > 0. )
	{
		/* carbon */
		c->daily_soil2_het_resp   = pot_soil2C_loss * RFS2S3;
		c->daily_soil3C_to_soil4C = pot_soil2C_loss * ( 1. - RFS2S3 );
		/* nitrogen */
		c->daily_soil3N_to_soil4N = pot_soil3C_loss / SOIL3_CN;
	}

	/* recalcitrant SOM pool (rf = 1.0, always mineralizing) */
	if ( c->soil4C > 0. )
	{
		/* carbon */
		c->daily_soil4_het_resp     = pot_soil4C_loss;
		/* nitrogen */
		c->daily_soil4N_to_soilMinN = pot_soil4C_loss / SOIL4_CN;
	}
}

