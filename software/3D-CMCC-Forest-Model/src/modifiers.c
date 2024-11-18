/*modifiers.c*/

/* includes */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "soil_settings.h"
#include "modifiers.h"
#include "constants.h"
#include "common.h"
#include "settings.h"
#include "logger.h"

/* externs */
extern settings_t* g_settings;
extern logger_t* g_debug_log;
//extern soil_settings_t *g_soil_settings;

void modifiers(cell_t *const c, const int layer, const int height, const int dbh, const int age, const int species, const meteo_daily_t *const meteo_daily,
		const meteo_annual_t *const meteo_annual , const int month , const int day)     
		
{
	double RelAge;

	/* constants and variables for Veroustraete's CO2 modifier computation */
	/* CO2 dependence */
	double KmCO2;	                       /* (ppmv CO2) affinity coefficients temperature dependent according to Arrhenius relationship */
	const double Ea1   = 59400.0;                /* (J mol-1) Activation energy for CO2 fixation (KmCO2 temp dependence)  Badger and Collatz 1977 */
	const double A1    = 2.419 * pow(10,13);     /* (ppmv) Arrhenius constant for KmCO2 tø dependence in ppm for t>=15  */
	const double Ea2   = 109600.0;	           /* (J mol-1) Activation energy for CO2 fixation for t<15 C */
	const double A2    = 1.976 * pow(10,22);     /* (ppmv) Arrhenius constant for KmCO2 tø dependence in ppm for t<15 C */
	/* O2 dependence */
	double KO2;	                           /* (% O2) Inhibition constant for 02 */
	const double EaKO2 = 13913.5;                /* (J mol-1) Activation energy for O2 inhibition */
	const double AKO2  = 8240;                   /* Arrhenius constantfor KO2 to dependence */

	/* CO2/O2 dependence */
	double tau;	                           /* CO2/O2 specifity ratio */
	const double Eatau = -42896.9;               /* (J mol-1) Activation energy for CO2/O2 specificity */
	const double Atau  = 7.87 * pow(10,-5);      /* (dimensionless) Arrhenius constant */

	double tleaf;                          /* (°C) daily leaf temperature */
	double tleaf_K;                        /* (K) daily leaf temperature */
	double v1, v2;

	/* constants for Wang et al., 2016 Nature Plants CO2 modifiers computation */
	const double beta   = 240.;   /* (dimensionless) definition ? */
	const double ni     = 0.89;   /* (dimensionless) vicosity of water relative to its value at 25 °C */
	const double ci     = 0.41;   /* (dimensionless) carbon cost unit for the maintenance of electron transport capacity */
	const double q10Kc  = 2.1;    /* (DIM) Q_10 for Kc */
	const double q10Ko  = 1.2;    /* (DIM) Q_10 for Ko */

	/* Badger and Collatz 1977 */
	const double Kc25   = 404;    /* (ubar or umol mol-1) Michaelis-Menten const carboxylase, 25 deg C  Badger and Collatz value*/
	const double Ea_Kc  = 59400;  /* (J mol-1) Activation energy for carboxylase */
	const double Ko25   = 248000; /* (ubar or umol mol-1) Michaelis-Menten const oxygenase, 25 deg C 248 Badger and Collatz, 278.4 Bernacchi et al., 2001 */
	const double Ea_Ko  = 36000;  /* (J mol-1) Activation energy for oxygenase */

//	/* Bernacchi et al., 2001 */
//	const double Kc25          = 404.9;  /* (ubar or umol mol-1) Michaelis-Menten const carboxylase, 25 deg C  Badger and Collatz value*/
//	const double Ea_Kc         = 79430;  /* (kJ mol-1) Activation energy for carboxylase */
//	const double Ko25          = 278400; /* (ubar or umol mol-1) Michaelis-Menten const oxygenase, 25 deg C 248 Badger and Collatz, 278.4 Bernacchi et al., 2001 */
//	const double Ea_Ko         = 36380;  /* (J mol-1) Activation energy for oxygenase */

	/* variables for Wang et al., 2016 Nature Plants CO2 modifiers computation */
	double Kc;           /* (ppm) effective Michaelis-Menten coefficienct for Rubisco */
	double gamma_star;   /* (ppm) Photorespiratory compensation point */
	double O2;           /* (Pa) atmospheric [O2] */
	double Ko;           /* (ppm) michaelis-menten constant for oxygenase reaction */
	double m0;
	double vpd;          /* (Pa) vpd */
	double Ca;           /* (ppmv) CO2 atmospheric concentration */
	double Ca_ref;       /* (ppmv) CO2 atmospheric concentration */

	double Av_rel, Aj_rel;
 

	static int modifier;
	static int test_assimilation;
            
        double tau_acl      ; 
        double S_acl_MAX    ; 
        double XO           ;  
        
        double St_acl ;   
        double  X_acl ; 

	age_t *a;
	species_t *s;

	a = &c->heights[height].dbhs[dbh].ages[age];
	s = &c->heights[height].dbhs[dbh].ages[age].species[species];

	logger(g_debug_log, "\n**DAILY_MODIFIERS**\n");

	Ca     = meteo_annual->co2Conc;
	Ca_ref = g_settings->co2Conc;

	tleaf    = meteo_daily->tday;
	tleaf_K  = meteo_daily->tday + TempAbs;

	/* convert vpd from hPa to Pa */
	vpd = meteo_daily->vpd / 100.;

	/* calculate atmospheric O2 in ppm, assumes 20.9% O2 by volume */
	O2  = O2CONC * 10000;

	/************************** ASSIMILATION CO2 MODIFIERS *************************/

	/***************************** VEROUSTRAETE'S VERSION **************************/
	/*******************************************************************************/

	/* CO2 MODIFIER FOR ASSIMILATION  */
	/* fertilization effect with rising CO2 from: Veroustraete 1994,
	 * Veroustraete et al., 2002, Remote Sensing of Environment */

	/* compute effective Arrhenius coefficienct for Rubisco */
	/* Dependence of KmCO2 on temperature data from Badger and Collatz (1976) */

	if ( tleaf >= 15 )
	{
		KmCO2 = A1 * exp( - Ea1 / ( Rgas * tleaf_K ) );
	}
	else
	{
		KmCO2 = A2 * exp ( - Ea2 / ( Rgas * tleaf_K ) );
	}

	/* Dependence of KO2 on temperature data */
	KO2 = AKO2 * exp ( - EaKO2 / ( Rgas * tleaf_K ) );

	/* dependence of assimilation rate on atmospheric carbon dioxyde concentration and competition by O2 */
	tau = Atau * exp ( - Eatau / ( Rgas * tleaf_K) );

	v1  = ( Ca - ( O2CONC / ( 2. * tau ) ) ) / ( Ca_ref - ( O2CONC / ( 2. * tau ) ) );
	v2  = ( KmCO2 * ( 1 + ( O2CONC / KO2 ) ) + Ca_ref ) / ( KmCO2 * ( 1. + ( O2CONC / KO2 ) ) + Ca );

	/* CO2 assimilation modifier */
	s->value[F_CO2_VER] = v1 * v2;

	/*******************************************************************************/

#if 0
	/*************************** MICHAELIS MENTEN KINETICS *************************/
	/*******************************************************************************/

	/* correct kinetic constants for temperature, and do unit conversions */
	Ko  = Ko25 * pow ( q10Ko , ( tleaf - 25. ) / 10. );

	/* compute effective Michaelis-Menten coefficient for Rubisco as in Collatz et al., (1991) see von Caemmerer 2000 "Biochemical model of leaf photosynthesis" */

	if ( tleaf > 15. )
	{
		Kc  = Kc25  * pow ( q10Kc , ( tleaf - 25. ) / 10. );
	}
	else
	{
		Kc  = Kc25  * pow ( 1.8 * q10Kc , ( tleaf - 15. ) / 10.) / q10Kc;
	}

	/*******************************************************************************/
#else

	/****************************** ARRHENIUS KINETICS *****************************/
	/*******************************************************************************/
	/* the enzyme kinetics built into this model are based on Medlyn (1999) */

	/* correct oxygenase kinetic constants for temperature */
	Ko  = Ko25 * exp ( Ea_Ko * ( tleaf - 25. ) / ( 298. * Rgas * ( tleaf + TempAbs ) ) );

	/* correct carboxylase kinetic constants for temperature */
	Kc  = Kc25 * exp ( Ea_Kc * ( tleaf - 25. ) / ( 298. * Rgas * ( tleaf + TempAbs ) ) );

#endif
	/*******************************************************************************/

	/* convert ubar umol --> Pa */
	Ko *= 0.1;

	/* convert ubar umol --> Pa */
	Kc *= 0.1;

	/* fixme  conversions */

	/* convert Pa --> ppm */
	Ko /= ( meteo_daily->air_pressure / 1e6 );

	/* convert Pa --> ppm */
	Kc /= ( meteo_daily->air_pressure / 1e6 );

	/*******************************************************************************/

	/*******************************************************************************/

	/* calculate gamma (ppm) CO2 compensation point due to photorespiration, in the absence of respiration */

#if 1
	/* it assumes Vomax/Vcmax = 0.21; Badger & Andrews (1974) */
	/* 0.5 because with 1 mol of oxygenations assumed to release 0.5 molCO2 by glycine decarboxilation (Farquhar and Busch, 2017) */

	gamma_star = 0.5 * 0.21 * O2 * Kc / Ko;

#else

	if ( tleaf < -1. )
	{
		/* note:  Maespa */
		gamma_star = 36.9 + 1.88*(-26.0) + 0.036*(-26.0)*(-26.0);
	}
	else
	{
		/* note: dePury and Farquhar 1997 method */
		gamma_star = 36.9 + 1.88 * ( tleaf - 25. ) + 0.036 * pow( ( tleaf - 25. ) , 2. );
	}

	/* note: Bernacchi et al., 2001 method (in umol/mol)*/
	gamma_star = 42.75 * exp ( 37830 * ( tleaf - 25.) / ( Rgas * ( tleaf + TempAbs ) * ( 25. + TempAbs ) ) );

	/* convert from umol --> Pa -->  ppm */
	gamma_star *=  0.1;
	gamma_star /= ( meteo_daily->air_pressure / 1e6 );

#endif

	/****************************** 'FRANKS ET AL' VERSION *************************/
	/*******************************************************************************/
	/* CO2 MODIFIER FOR ASSIMILATION  */
	/* see Franks et al., (2013) NewPhytologist Eq 5 */

	/* solve for Av and Aj using the quadratic equation, substitution for Ci
		from A = g(Ca-Ci) into the equations from Farquhar, von Caemmerer and Berry (1980)
		and Bernacchi et al (2003) (for values 4.5 and 10.5) :

		       Vmax (Ci - gamma)
		Av =  -------------------   -   Rd
		      Ci + Kc (1 + O2/Ko)


		         J (Ci - gamma)
		Aj  =  -------------------  -   Rd
	           4.5 Ci + 10.5 gamma
	 */

	/*******************************************************************************/
	/* note: all variables are in ppm */

	/* compute FCO2 modifier (Aj_rel) */

	/* note: original Franks et al., (2013) */
	//Aj_rel = ( ( Ca - gamma_star ) * ( Ca_ref + 2 * gamma_star ) ) / ( ( Ca + 2. * gamma_star ) * ( Ca_ref - gamma_star ) );

	/* note Following parameterization from Bernacchi et al., (2001) and applied to Franks et al., (2013) */
	Aj_rel = ( ( Ca - gamma_star ) * ( ( 4.5 * Ca_ref ) + ( 10.5 * gamma_star ) ) ) / ( ( ( 4.5 * Ca )  + ( 10.5 * gamma_star ) ) * ( Ca_ref - gamma_star ) );

	/*******************************************************************************/

	/* compute F CO2 modifier (Av_rel) */

	Av_rel = ( ( Ca - gamma_star ) * ( Ca_ref + Kc * ( 1. + ( O2 / Ko ) ) ) ) / ( ( Ca + Kc * ( 1. + ( O2 / Ko ) ) ) * ( Ca_ref - gamma_star ) );

	/*******************************************************************************/

	test_assimilation = 0;

	switch ( test_assimilation )
	{
	case 0:
		/* estimate A as the minimum of (Av,Aj) */

		s->value[F_CO2_FRANKS] = MIN ( Av_rel, Aj_rel );
		break;
	case 1:
		/* estimate A as Av */

		s->value[F_CO2_FRANKS] = Av_rel;
		break;
	case 2:
		/* estimate A as Aj */

		s->value[F_CO2_FRANKS] = Aj_rel;
		break;
	}

	/******************************************************************************************/

	modifier = 1; /* 0 for Veroustraete; 1 for Franks */

	/* selection for CO2 modifier to be used */
	switch ( modifier )
	{
	case 0:
		s->value[F_CO2] = s->value[F_CO2_VER];
		break;
	case 1:

		s->value[F_CO2] = s->value[F_CO2_FRANKS];
		break;
	}

	/****************************************************************************************************************/

	/************************** TRANSPIRATION CO2 MODIFIERS ***************************/

	/* CO2 MODIFIER FOR TRANSPIRATION  */
	/* limitation effects on maximum stomatal conductance from:
	 * Frank et al., 2013 New Phytologist
	 * Hidy et al., 2016 Geosc. Model Dev.
	 */

	s->value[F_CO2_TR] = 39.43 * pow ( Ca , -0.64 );

	/****************************************************************************************************************/

	/* LIGHT MODIFIER */
	/* (Following Makela et al. , 2008, Peltioniemi et al. 2012) */
	s->value[F_LIGHT_MAKELA]        = 1. / ( ( s->value[GAMMA_LIGHT] * s->value[APAR])       + 1. );
	s->value[F_LIGHT_SUN_MAKELA]    = 1. / ( ( s->value[GAMMA_LIGHT] * s->value[APAR_SUN])   + 1. );
	s->value[F_LIGHT_SHADE_MAKELA]  = 1. / ( ( s->value[GAMMA_LIGHT] * s->value[APAR_SHADE]) + 1. );

	/* LIGHT MODIFIER (BIOME METHOD) */
	/* following Biome-BGC */
	/* photosynthetic photon flux density conductance control */

	/* for overall leaves */
	s->value[F_LIGHT]       = s->value[APAR]        / ( PARD50 + s->value[APAR] );
	logger(g_debug_log, "fLIGHT       = %f \n", s->value[F_LIGHT]);

	/* for sun leaves */
	s->value[F_LIGHT_SUN]   =  s->value[APAR_SUN]   / ( PARD50 + s->value[APAR_SUN] );
	logger(g_debug_log, "fLIGHT_SUN   = %f \n", s->value[F_LIGHT_SUN]);

	/* for shaded leaves */
	s->value[F_LIGHT_SHADE] =  s->value[APAR_SHADE] / ( PARD50 + s->value[APAR_SHADE] );
	logger(g_debug_log, "fLIGHT_SHADE = %f \n", s->value[F_LIGHT_SHADE]);

	/********************************************************************************************/

	/* TEMPERATURE MODIFIER */
	if ( meteo_daily->tday == NO_DATA )
	{
		if ( ( meteo_daily->tavg <= s->value[GROWTHTMIN]) || (meteo_daily->tavg >= s->value[GROWTHTMAX] ) )
		{
			s->value[F_T] = 0.;
		}
		else
		{
			s->value[F_T] = ( ( meteo_daily->tavg - s->value[GROWTHTMIN]) / (s->value[GROWTHTOPT] - s->value[GROWTHTMIN] ) ) *
					pow ( ( ( s->value[GROWTHTMAX] - meteo_daily->tavg) / (s->value[GROWTHTMAX] - s->value[GROWTHTOPT] ) ),
							( ( s->value[GROWTHTMAX] - s->value[GROWTHTOPT]) / (s->value[GROWTHTOPT] - s->value[GROWTHTMIN] ) ) );
		}
	}
	else
	{
		if ( ( meteo_daily->tday <= s->value[GROWTHTMIN]) || (meteo_daily->tday >= s->value[GROWTHTMAX] ) )
		{
			s->value[F_T] = 0.;
		}
		else
		{
			s->value[F_T] = ( ( meteo_daily->tday - s->value[GROWTHTMIN]) / (s->value[GROWTHTOPT] - s->value[GROWTHTMIN] ) ) *
					pow ( ( ( s->value[GROWTHTMAX] - meteo_daily->tday) / (s->value[GROWTHTMAX] - s->value[GROWTHTOPT] ) ),
							( ( s->value[GROWTHTMAX] - s->value[GROWTHTOPT]) / (s->value[GROWTHTOPT] - s->value[GROWTHTMIN] ) ) );
		}
	}
	logger(g_debug_log, "fT = %f\n", s->value[F_T]);

	/* check */
	CHECK_CONDITION(s->value[F_T], >, 1);
	CHECK_CONDITION(s->value[F_T], <, ZERO);

	/********************************************************************************************/

	/* VPD MODIFIER */
	//logger(g_debug_log, "--RH = %f %%\n", met[month].rh);
	//logger(g_debug_log, "--vpd = %f mbar", vpd);

	//The input VPD data is in KPa
	//if the VPD is in KPa
	//Convert to mbar
	//1 Kpa = 10 mbar
	//1 hPa = 1 mbar
	//s->value[F_VPD] = exp (- s->value[COEFFCOND] * vpd) * 10);
	//convert also COEFFCOND multiply it for
	s->value[F_VPD] = exp ( - s->value[COEFFCOND] * meteo_daily->vpd );
	logger(g_debug_log, "fVPD = %f\n", s->value[F_VPD]);

	/* check */
	CHECK_CONDITION(s->value[F_VPD], >, 1);
	CHECK_CONDITION(s->value[F_VPD], <, ZERO);

	/********************************************************************************************/

	//test following biome-bgc it doesn't seems to work properly here (too many higher values for gpp and le
	/* vapor pressure deficit multiplier, vpd in Pa */
#if 0
	if (meteo_daily->vpd < vpd_open)    /* no vpd effect */
		s->value[F_VPD] = 1.;
	else if (meteo_daily->vpd > vpd_close)   /* full vpd effect */
		s->value[F_VPD] = 0.;
	else                   /* partial vpd effect */
		s->value[F_VPD] = (vpd_close - meteo_daily->vpd) / (vpd_close - vpd_open);
#endif

	/* AGE MODIFIER */
	/* note: in some ways it reflects the Hydraulic Limitation Hypothesis (HLH) of Yoder et al., (1994)
	 * that asserts that increasing tree height and branch length (and then age), xylem path and length
	 * and thus hydraulic resistance increases (see also Grote et al., 2016) */

	if ( a->value != 0 )
	{
		/* for TIMBER */
		/* AGE FOR TIMBER IS THE EFFECTIVE AGE */
		RelAge = (double)a->value / s->value[MAXAGE];
		s->value[F_AGE] = ( 1. / ( 1. + pow ( ( RelAge / (double)s->value[RAGE] ), (double)s->value[NAGE] ) ) );
		logger(g_debug_log, "fAge = %f\n", s->value[F_AGE]);
	}
	else
	{
		s->value[F_AGE] = 1.;
		logger(g_debug_log, "no data for age F_AGE = 1\n");
	}

	/* check */
	CHECK_CONDITION(s->value[F_AGE], >, 1.);
	CHECK_CONDITION(s->value[F_AGE], <, ZERO);
	
	/********************************************************************************************/
	
	// Acclimatation factor to high temperature in spring, after dormance period
	// ref. Mäkelä et al.2004, 2008 ; Luo et al 2023 (New Phytologist)
	// Gennaretti et al.2017 (applied to T in Vcmax)

	
	//5p7
	if ( g_settings->Cold_accl )
	{ 
	tau_acl =   s->value[TAU]     ; 
        S_acl_MAX = s->value[Smax]    ; 
        XO  =       s->value[X0]      ;
        
        //CURRENTLY NOT USED
        //#if 0
	
	    //s->value[F_ACCL] = 1. ;  
	
         //#else 

	 if (c->doy == 1 )  
	
	 {

	 s->value[F_ACCL] = 0. ;
	 s->value[old_X_acl] = meteo_daily->tmin ;
	
	 }
	 
	 if ( c->doy > 1   && month < JULY ) 
	 {

         X_acl = s->value[old_X_acl]  + (1/tau_acl) *(meteo_daily->tmin - s->value[old_X_acl] ) ;
         St_acl  = X_acl-XO ;
         St_acl  = MAX(St_acl,0.)    ;
            
         s->value[F_ACCL] = St_acl/S_acl_MAX;
         s->value[F_ACCL] = MIN(s->value[F_ACCL],1.) ;
         
         s->value[old_X_acl]  = X_acl;
   
 	 } else {
	 //old_St
	 s->value[F_ACCL] = 1. ;    // at 1 january put the value to 1
	 }

  //#endif 
	

    }  //TODO end flag cold acclimatation

	/********************************************************************************************/

        //TODO add flag for N modifier
        // currently not used
	
//#if 1
	
	// s->value[F_NUTR] = 1. ;  
	
//#else

if ( g_settings->N_reg )
	{

	/* SOIL NUTRIENT MODIFIER */
	s->value[F_NUTR] = 1. - ( 1. - c->fn0 ) * pow ( ( 1. - c->fr), c->fnn );
	logger(g_debug_log, "fNutr = %f\n", s->value[F_NUTR]);

	/* check */
	CHECK_CONDITION( s->value[F_NUTR], >, 1 );
	CHECK_CONDITION( s->value[F_NUTR], <, ZERO );	

 //#endif

	} else {

	s->value[F_NUTR] = 1. ; 

	}

	//test 25 nov 2016
#if 0
	/* SOIL NUTRIENT MODIFIER */

	/* compute fertility rate based on N available and N demand */
	s->values[F_NUTR] = c->soilN / s->value[NPP_gN_DEMAND];

	/* check */
	if ( s->values[F_NUTR] > 1.) s->values[F_NUTR] = 1.;
	logger(g_debug_log, "fNutr = %f\n", s->value[F_NUTR]);

	/* check */
	CHECK_CONDITION( s->value[F_NUTR] , > , 1);
	CHECK_CONDITION( s->value[F_NUTR] , < , ZERO);
#endif
	/********************************************************************************************/
	/* (MPa) water potential of soil and leaves */
	/* SOIL MATRIC POTENTIAL */

	/* convert kg/m2 or mm --> m3/m2 --> m3/m3 */
	/* 100 mm H20 m^-2 = 100 kg H20 m^-2 */
	/* calculate the soil pressure-volume coefficients from texture data */
	/* Uses the multivariate regressions from Cosby et al., 1984 */
	/* volumetric water content */


	// c->vwc = c->asw / c->max_asw_fc; //this was the old (bugged) equation

	/* volumetric water content */
	c->vwc = c-> asw/((c->soil_depth/100.)*1000.);

	/* soil matric potential */
	
	c->psi = c->psi_sat * pow((c->vwc/c->vwc_sat), c->soil_b);


	/* no water stress */
	
	if (c->psi > s->value[SWPOPEN])
	{
		s->value[F_PSI] = 1.;
	}
	/* full water stress */
	
	else if (c->psi <= s->value[SWPCLOSE])
	{
		/* forced to  0.3 to avoid zero values */
		/* see: Clark et al., 2011 for JULES model impose 0.2 */
		s->value[F_PSI] = WATER_STRESS_LIMIT ;
	}
	/* partial water stress */
	else
	{
		s->value[F_PSI] = ( s->value[SWPCLOSE] - c->psi ) / ( s->value[SWPCLOSE] - s->value[SWPOPEN] );

		/* for consistency with complete stress values */
		if(s->value[F_PSI] < WATER_STRESS_LIMIT) s->value[F_PSI] = WATER_STRESS_LIMIT;
	}


	s->value[F_SW] = s->value[F_PSI];
	logger(g_debug_log, "fSW = %f\n", s->value[F_PSI]);

	/* check */
	CHECK_CONDITION(s->value[F_SW], >, 1.);
	CHECK_CONDITION(s->value[F_SW], <, ZERO);

	/********************************************************************************************/

	/* PHYSIOLOGICAL MODIFIER */
	s->value[PHYS_MOD] = MIN (s->value[F_VPD], (s->value[F_SW] * s->value[F_AGE]));

	/* check */
	CHECK_CONDITION(s->value[PHYS_MOD], >, 1.);
	CHECK_CONDITION(s->value[PHYS_MOD], <, ZERO);

	/********************************************************************************************/

	/*SOIL DROUGHT MODIFIER*/
	//(see Duursma et al., 2008)rev_Angelo
	/*
		//to put in species.txt
		//numbers are not real just used to compile!!!!!!!!
		double leaf_res = 1; //leaf specific plant hydraulic resistance
		double min_leaf_pot = 1; //minimum leaf water potential


		//to put in site.txt ?????
		double soil_res = 1; //soil hydraulic resistance
		double psi0 = 2; //dry soil water potential in MPa
		double soil_coeff = 1; //empirical soil coefficient

		//soil average dimension particle
		//value are averaged from limits in site.txt
		double clay_dim = 0.001; //clay avg dimension of particle
		double silt_dim =  0.026;//silt avg dimension of particle
		double sand_dim =  1.025;//sand avg dimension of particle

	    double bulk_pot; //bulk soil water potential
	    double asw_vol; //available soil water in volume
	    double eq;
		double eq1;
		double eq2;
		double sat_soil_water_cont; //saturated soil water content (m^3 m^-3)
		double soil_avg_dim; //soil mean particle diameter in mm
		double sigma_g; //geometric standard deviation in particle size distribution (mm)
		double pentry_temp; //soil entry water potential (MPa)
		double pentry; //correction for bulk density effects
		double bsl; //coefficient in soil water release curve (-)
		double soil_water_pot_sat; //soil water potential at saturation
		double sat_conduct; //saturated conductivity
		double specific_soil_cond; //specific soil hydraulic conductance
		double leaf_specific_soil_cond;

		//compute soil hydraulic characteristics from soil granulometry
		//from model Hydrall
		eq1 = (g_soil_settings->values[SOIL_CLAY_PERC] * log(clay_dim)) + (g_soil_settings->values[SOIL_silt_perc * log(silt_dim)) + (g_soil_settings->values[SOIL_sand_perc * log(sand_dim));
		logger(g_debug_log, "eq1 = %f\n", eq1);

		//soil mean particle diameter in mm
		soil_avg_dim = exp(eq1);
		logger(g_debug_log, "soil_avg_dim = %f\n", soil_avg_dim);

	    eq2 = sqrt ((pow ((g_soil_settings->values[SOIL_CLAY_PERC] * log(clay_dim)),2)) + (pow ((g_soil_settings->values[SOIL_sand_perc * log(sand_dim)),2)) + (pow ((g_soil_settings->values[SOIL_silt_perc * log(silt_dim)),2)));
	    logger(g_debug_log, "eq2 = %f\n", eq2);

	    //geometric standard deviation in particle size distribution (mm)
	    sigma_g = exp(eq2);
	    logger(g_debug_log, "sigma_g = %f\n", sigma_g);

	    //soil entry water potential (MPa)
	    pentry_temp = -0.5 / sqrt(soil_avg_dim)/1000;
	    logger(g_debug_log, "pentry_temp = %f\n", pentry_temp);
	    //correction for bulk density effects with dens = 1.49 g/cm^3
	    pentry = pentry_temp * pow ((c->bulk_density / 1.3), (0.67 * bsl));
	    logger(g_debug_log, "pentry = %f\n", pentry);

	    bsl = -2 * (pentry * 1000) + 0.2 * sigma_g;
	    logger(g_debug_log, "bsl = %f\n", bsl);

	    //saturated soil water content
	    sat_soil_water_cont= 1.0 - (c->bulk_density/2.56);
	    logger(g_debug_log, "soil water content at saturation = %f\n", sat_soil_water_cont);

	    eq = pentry * pow ((sat_soil_water_cont / c->soil_moist_ratio), bsl);
	    logger(g_debug_log, "eq = %f\n", eq);

	    //compute bulk soil water potential
	    //for psi see Magnani xls
	    bulk_pot = Maximum (eq, min_leaf_pot);
	    logger(g_debug_log, "bulk soil water potential = %f\n", bulk_pot);

	    //compute leaf-specific soil hydraulic conductance
		leaf_specific_soil_cond = sat_conduct * pow ((soil_water_pot_sat / bulk_pot), (2 + (3 / soil_coeff)));
		logger(g_debug_log, "leaf-specific soil hydraulic conductance = %f\n", leaf_specific_soil_cond);


		s->value[F_DROUGHT] = (leaf_res * (bulk_pot - min_leaf_pot)) / (- min_leaf_pot * ((leaf_res + soil_res) * bulk_pot));
		logger(g_debug_log, "F_DROUGHT = %f\n", s->value[F_DROUGHT]);
	 */

	logger(g_debug_log, "-------------------\n");
}
