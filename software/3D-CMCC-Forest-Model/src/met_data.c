/* met_data.c */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "math.h"
#include "constants.h"
#include "soil_settings.h"
#include "topo.h"
#include "settings.h"
#include "logger.h"
#include "matrix.h"
#include "met_data.h"
#include "common.h"

#define WEIGHTED_DAYS	10
#define AVERAGED_DAYS	10

//extern soil_settings_t *g_soil_settings;
extern topo_t *g_topo;
extern logger_t* g_debug_log;

void Sat_vapour_pressure(cell_t *const c, const int day, const int month, const int year)
{
	double e0max;                                    /* saturation vapour pressure at the maximum air temperature (KPa) */
	double e0min;                                    /* saturation vapour pressure at the minimum air temperature (KPa) */
	const double A = 17.27;
	const double Tstroke = 36;
	double TmaxK;
	double TminK;

	meteo_t *met;
	met = (meteo_t*) c->years[year].m;

	TmaxK = met[month].d[day].tmax + TempAbs;
	TminK = met[month].d[day].tmin + TempAbs;

	/* following Allen et al., 1998 */
	/* compute saturation vapour pressure at the maximum and minimum air temperature (KPa) */
	e0max = ESTAR * exp((A*met[month].d[day].tmax)/(TmaxK - Tstroke));
	e0min = ESTAR * exp((A*met[month].d[day].tmin)/(TminK - Tstroke));

	/* compute weighted mean saturation vapour pressure at the air temperature (KPa)*/
	met[month].d[day].es = ((e0max*met[month].d[day].ni) + (e0min*(1.0-met[month].d[day].ni)));

	/* compute actual vapour pressure derived from relative humidity data (KPa) */
	met[month].d[day].ea = (met[month].d[day].rh_f/100.)*met[month].d[day].es;

	if (met[month].d[day].ea < 0.)
	{
		met[month].d[day].ea = 0.;
	}
}

static double get_daily_potential_radiation(const double latitude, const double longitude, const int d_, const double t_) {
	double localstandardtime;
	double localapparentsolartime;
	double tthet;
	double signedLAS;
	double omega;
	double decl_rad;
	double lat_rad;
	double theta_rad;
	double Rpot;
	double Rpot_h;

	localstandardtime = t_;
	tthet = 2.*Pi*(d_-1.) / 365.;

	localapparentsolartime = localstandardtime;
	signedLAS = 12.-localapparentsolartime;
	signedLAS = fabs(signedLAS);

	omega = -15.*signedLAS;
	decl_rad = 0.006918-0.399912*cos(tthet)+0.070257*sin(tthet)-0.006758*cos(2*tthet)+0.000907*sin(2*tthet)-0.002697*cos(3*tthet)+0.00148*sin(3*tthet);
	lat_rad = latitude*Pi/180.;

	theta_rad = acos(sin(decl_rad)*sin(lat_rad)+cos(decl_rad)*cos(lat_rad)*cos(omega*Pi/180.));

	/* Compute Potential short wave downward radiation [W m-2] (short wave clear sky radiation) */
	Rpot = Q0*(1.00011+0.034221*cos(tthet)+0.00128*sin(tthet)+0.000719*cos(2*tthet)+0.000077*sin(2*tthet));
	Rpot_h = Rpot*cos(theta_rad);

	return (Rpot_h > 0) ? Rpot_h : 0.0;
}

static double compute_potential_rad(const double lat, const double lon, const int day) {
#define ROWS_PER_DAY	(60*24)
	int i;
	double rpot;

	/* potential radiation is computed for each minute */
	rpot = 0.;
	for ( i = 0; i < ROWS_PER_DAY; ++i ) {
		double doy;
		double hrs;
		doy = ((double)(day*ROWS_PER_DAY+i) / ROWS_PER_DAY) + 1;
		hrs = (double)((day*ROWS_PER_DAY+i) % ROWS_PER_DAY) / 60.;
		rpot += get_daily_potential_radiation(lat, lon, doy, hrs);
	}
	rpot /= ROWS_PER_DAY;
	return rpot;

#undef ROWS_PER_DAY
}

void compute_atm_lw_downward_W(cell_t *const c, const int day, const int month, const int year) {

	/*calculates longwave radiation from Vapour Pressure and Air Temp. and Radiation
	-------------------------------------------------------------------------------
	implemented by Sönke Zaehle (szaehle@bgc-jena.mpg.de)
	MPI for Biogeochemistry, Jena Germany

	DISCLAIMER:
	This software comes with absolutely no warrenty whatsoever wrt to its performance
	or scientific correctness. Any error associated with results produced from
	this routine is not in the responsibility of the author of this routine
	-------------------------------------------------------------------------------
	Input:
	 -Tair: Air temperature [°C]
	 -sw_downward_W: Shortwave downward radiation [W m-2]
	 -Rg_pot: Potential shortwave downward radiation [W m-2]
	 -VP:  Vapour Pressure [Pa]
	Output:
	 -atm_lw_downward_W: Longwave downward radiation [W m-2]
	-------------------------------------------------------------------------------
	Calculation of Longwave follows the JSBACH algorithm
	Downward long wave radiation flux "R_d" [W/m^2] is according to [1],[2] computed by
	(1) R_d = r_cloud * epsA * SBC_W * T^4,
	where "TairK" is the air temperature [K] (in the bottom layer?), SBC_W is the Stefan-Boltzmann constant (in W),
	"epsA" is the emissivity of the cloudless atmosphere given by
	(2) epsA = epsA0 * (e_A/T)^(1/7).
	Here "e_vap" is the vapor pressure in [Pa] from above, and epsA0 = 0.64 (units: [(K/Pa)^(1/7)]) a constant (see [2]).
	The factor "r_cloud" in (1) corrects for clouds and is given by
	(3) r_cloud = 1+0.22*n^2,
	where "n" is the cloudcover.
	[1] W. Knorr, "Satellite remote sensing and modelling of the global CO2 exchange of land vegetation", Examensarbeit 49,
	(Max Planck Institut fuer Meteorologie, Hamburg, 1998).
	[2] W. Brutsaert, "Evaporation into the Atmosphere", (Reidel, Dordrecht, 1982), pp. 299.
	Calculation of Vapour Pressure follows Monteith & Unsworth 2008, page 11f [3] */

	double fpar;
	//double cloud_cover;
	//double r_cloud;
	double vp;
	//double epsA;                                /* emissivity of the cloudless atmosphere */
	//double lw_downward_W;

	/* input met data */
	double sw_in;
	double ta;
	double vpd;                        /* deficit vapour pressure at the air temperature in hPa */
	double sw_pot_in;
	double esat;                        /* saturation vapour pressure at the air temperature (KPa)*/
	double TairK;

	meteo_t *met;
	met = (meteo_t*) c->years[year].m;

	sw_in = met[month].d[day].sw_downward_W;
	ta = met[month].d[day].tavg;
	vpd = met[month].d[day].vpd;
	sw_pot_in = met[month].d[day].sw_pot_downward_W;
	esat = met[month].d[day].es;
	TairK = ta + TempAbs;


	fpar = 0.;
	if ( ! (0. == sw_pot_in) && ! IS_INVALID_VALUE(sw_in) )
	{
		fpar = sw_in / sw_pot_in;
		if ( fpar < 0. )
		{
			fpar = 0.;
		}
	}

	/* Cloud cover and cloud correction factor Eq. (3) */
	met[month].d[day].cloud_cover_frac = 1. - (fpar - 0.5) / 0.4;
	if ( met[month].d[day].cloud_cover_frac > 1. )
	{
		met[month].d[day].cloud_cover_frac = 1.;
	}
	if ( met[month].d[day].cloud_cover_frac < 0. )
	{
		met[month].d[day].cloud_cover_frac = 0.;
	}
	met[month].d[day].cloud_cover_frac_corr = 1. + 0.22 * pow( met[month].d[day].cloud_cover_frac, 2. );

	/* Saturation and actual Vapour pressure [3], and associated emissivity Eq. [2] */
	/* convert esat kPa --> Pa */
	esat *= 1000.0;

	/* convert vpd hPa --> Pa */
	vpd *= 100.0;

	//fixme it is .ea
	vp = esat - vpd;
	if ( vp < 0.0 )
	{
		vp = 3.3546e-004;
	}
	met[month].d[day].emis_atm_clear_sky = 0.64 * pow( vp / TairK, 0.14285714 );

	/* compute atmopsheric emissivity based on cloud cover corrected */
	met[month].d[day].emis_atm = met[month].d[day].cloud_cover_frac_corr * met[month].d[day].emis_atm_clear_sky;

	/* Atmospheric Longwave radiation flux downward Eq. [1] */
	met[month].d[day].atm_lw_downward_W = met[month].d[day].emis_atm * SBC_W * pow(TairK, 4);
	if ( (met[month].d[day].atm_lw_downward_W < 10.) || (met[month].d[day].atm_lw_downward_W > 1000.) )
	{
		met[month].d[day].atm_lw_downward_W = INVALID_VALUE;
	}
}

void Radiation (cell_t *const c, const int day, const int month, const int year)
{
	double a = 107.;                                           /* (W/m)  empirical constants for long wave radiation computation */
	double b = 0.2;                                            /* (unit less) empirical constants for long wave radiation computation */
	double TmaxK, TminK;

	meteo_t *met;
	met = (meteo_t*) c->years[year].m;

	TmaxK = met[month].d[day].tmax + TempAbs;
	TminK = met[month].d[day].tmin + TempAbs;

	/***************************************************************************************************************************************/

	/* SHORT WAVE RADIATION */
	/* INCOMING SHORT WAVE RADIATION */

	/* compute downward Short-Wave radiation */
	met[month].d[day].sw_downward_MJ         = met[month].d[day].solar_rad;

	/* convert into W/m2 */
	met[month].d[day].incoming_sw_downward_W = met[month].d[day].sw_downward_MJ * MJ_TO_W;
	met[month].d[day].sw_downward_W          = met[month].d[day].incoming_sw_downward_W ;

	met[month].d[day].sw_pot_downward_W      = compute_potential_rad(c->lat, c->lon, day);

	/* convert incoming Short-Wave flux in PAR from MJ/m2/day to molPAR/m2/day (Biome-BGC method) */
	met[month].d[day].incoming_par           = met[month].d[day].sw_downward_MJ * RAD2PAR * EPAR;
	met[month].d[day].par                    = met[month].d[day].incoming_par;

	/***************************************************************************************************************************************/

	/* LONG WAVE RADIATION */
	/* INCOMING LONG WAVE RADIATION */

	/* compute Long-Wave radiation flux downward (following JSBACH model) */
	compute_atm_lw_downward_W(c, day, month, year);

	/* NET LONG WAVE RADIATION */

	//fixme to avoid crash in model for negative "ea" values use different calculation of long_wave_radiation following Prentice (IT HAS TO BE SOLVED ANYWAY)
	if(met[month].d[day].ea > 0.0)
	{
		/* following Allen et al., 1998 */
		/* it represents the outgoing part of long wave radiation */
		met[month].d[day].lw_net_MJ = SBC_MJ * ( pow ( ( ( TmaxK + TminK ) / 2.) , 4 ) ) * ( 0.34 - 0.14 * ( sqrt(met[month].d[day].ea ) ) ) * met[month].d[day].cloud_cover_frac;


		/* convert into W/m2 */
		met[month].d[day].lw_net_W = met[month].d[day].lw_net_MJ * MJ_TO_W;

		/*****************************************************************************************/
	}
	else
	{
		//todo check it Prentice says "net upward long-wave flux"
		/* following Prentice et al., 1993 */
		met[month].d[day].lw_net_W = ( b + ( 1. - b ) * met[month].d[day].ni ) * ( a - met[month].d[day].tavg );

		/* convert into MJ/m^2 day */
		met[month].d[day].lw_net_MJ = met[month].d[day].lw_net_W * W_TO_MJ;

		/*****************************************************************************************/
	}

	/***************************************************************************************************************************************/

	/* NET RADIATION */

	// net radiation based on 3-PG method // see Landsberg Terrestrial ecology Chapter 2: Weather and Energy Balance //
	// Note: Jan 2022  we set the long radiation terms QA as an average value, 
	// for both coniferous and broadleaved forests, -60 instead of -90 (QA) 
	// This has of course to be FIXME  
	
	met[month].d[day].Net_rad_threePG = QA + QB * ( met[month].d[day].solar_rad * pow ( 10., 6 ) / ( met[month].d[day].daylength * 3600. ) );

	/***************************************************************************************************************************************/


}

void Check_prcp(cell_t *c, const int day, const int month, const int year)
{
	meteo_t *met;
	met = c->years[year].m;
	if(met[month].d[day].prcp > 0.)
	{
		if (met[month].d[day].tavg > 0.)
		{
			met[month].d[day].rain = met[month].d[day].prcp;
			met[month].d[day].snow = 0.;
		}
		else
		{
			met[month].d[day].snow = met[month].d[day].prcp;
			met[month].d[day].rain = 0.;
		}
	}
	else
	{
		met[month].d[day].rain = 0.;
		met[month].d[day].snow = 0.;
	}

	/* following Lagergren et al., 2006 */
	/* impose zero value for tsoil in case of snow presence */
	//todo move to soil when snow melt and subl will go to soil
	//if(c->snow_pack != 0) met[month].d[day].tsoil = 0.0;
}

void Day_Length(cell_t *c, const int day, const int month, const int year)
{

	double ampl;  //seasonal variation in Day Length from 12 h

	//double adjust_latitude;

	/* BIOME-BGC version */
	/* Running-Coughlan 1988, Ecological Modelling */

	meteo_t *met;
	met = c->years[year].m;

	/* compute doy for DayLength function */
	if ( ! day && month == JANUARY )
	{
		c->doy_daylength = 0;
	}
	++c->doy_daylength;


	//4/apr/2016
	//test following Schwalm & Ek 2004 instead of only geographical latitude adjusted latitude is used
	// for every 125m in altitude 1° in latitude is added
	//	adjust_latitude = g_topo->values[TOPO_ELEV] / 125.0;
	//	ampl = (exp (7.42 + (0.045 * (g_soil_settings->values[SOIL_LAT]+adjust_latitude)))) / 3600;
	//	met[month].d[day].daylength = ampl * (sin ((doy - 79) * 0.01721)) + 12;
	//logger(g_debug_log, "with altitude = %f\n", met[month].d[day].daylength);

	ampl = ( exp ( 7.42 + ( 0.045 * c->lat ) ) ) / 3600.;

	/* compute daylength in hours */
	met[month].d[day].daylength     = ampl * ( sin ( ( c->doy_daylength - 79. ) * 0.01721 ) ) + 12.;

	/* compute daylength in seconds */
	met[month].d[day].daylength_sec = met[month].d[day].daylength * 3600.;


	/* compute fraction of daytime */
	met[month].d[day].ni = met[month].d[day].daylength / 24.;

}

/* following Running et al., 1987 */
void Daily_avg_temperature(meteo_t *met, const int day, const int month)
{
	if ( NO_DATA == met[month].d[day].tavg )
	{
		if ( (NO_DATA == met[month].d[day].tmax) && (NO_DATA == met[month].d[day].tmin) )
		{
			logger(g_debug_log, "NO DATA FOR TEMPERATURE!\n");
		} else {
			met[month].d[day].tavg =  (0.606 * met[month].d[day].tmax) + (0.394 * met[month].d[day].tmin);
		}
	}
}

void Psychrometric(meteo_t *met, const int day, const int month)
{
	/* compute psychrometric (KPa/°C) constant as in Allen et al., 1998 */
	met[month].d[day].psych = ( ( CP / 1e6 ) * ( met[month].d[day].air_pressure / 1e3 ) ) / ( MWratio * ( met[month].d[day].lh_vap / 1e6 ) );
}

void Daylight_avg_temperature(meteo_t *const met, const int day, const int month)
{
	/* BIOME-BGC version */
	/* Running-Coughlan 1988, Ecological Modelling */

	if ( met[month].d[day].tmax != NO_DATA && met[month].d[day].tmin != NO_DATA )
	{
		met[month].d[day].tday = 0.45 * ( met[month].d[day].tmax - met[month].d[day].tavg ) + met[month].d[day].tavg;

		/* but see also Running et al., (1987) as in CenW Kirschbaum model */
		//met[month].d[day].tday = ( 0.606 * met[month].d[day].tmax ) + ( 0.394 * met[month].d[day].tmin );
	}
	else
	{
		met[month].d[day].tday = NO_DATA;
		logger(g_debug_log, "NO TMAX and TMIN can't compute TDAY!!! \n");
	}
}

void Nightime_avg_temperature(meteo_t *const met, const int day, const int month)
{
	/* BIOME-BGC version */
	/* Running-Coughlan 1988, Ecological Modelling */

	if (met[month].d[day].tday != NO_DATA )
	{
		met[month].d[day].tnight = ( met[month].d[day].tday + met[month].d[day].tmin ) / 2. ;
	}
	else
	{
		met[month].d[day].tnight = NO_DATA;
		logger(g_debug_log, "NO TMAX and TMIN can't compute TNIGHT!!! \n");
	}
}

void Thermic_sum (cell_t * c, meteo_t *met, const int day, const int month, const int year)
{

	if ( !day && !month )
	{
		met[month].d[day].thermic_sum = 0.;

		if(met[month].d[day].tavg > GDD_BASIS)
		{
			met[month].d[day].thermic_sum = met[month].d[day].tavg - GDD_BASIS;
			c->previous_thermic_sum = met[month].d[day].thermic_sum;
		}
		else
		{
			met[month].d[day].thermic_sum = 0.;
			c->previous_thermic_sum          = 0.;
		}
		if (met[month].d[day].tavg == NO_DATA)
			logger(g_debug_log, "tavg NO_DATA!!\n");
	}
	else
	{
		if(met[month].d[day].tavg > GDD_BASIS)
		{
			met[month].d[day].thermic_sum = c->previous_thermic_sum + (met[month].d[day].tavg - GDD_BASIS);
			c->previous_thermic_sum = met[month].d[day].thermic_sum;
		}
		else
		{
			met[month].d[day].thermic_sum = c->previous_thermic_sum;
		}
		if (met[month].d[day].tavg == NO_DATA)
			logger(g_debug_log, "tavg NO_DATA!!\n");
	}
}

void Air_pressure(meteo_t *met, const int day, const int month)
{
	double t1, t2;

	/* compute air pressure */
	/* BIOME-BGC version */

	/* daily atmospheric pressure (Pa) as a function of elevation (m) */
	/* From the discussion on atmospheric statics in:
	Iribane, J.V., and W.L. Godson, 1981. Atmospheric Thermodynamics, 2nd
		Edition. D. Reidel Publishing Company, Dordrecht, The Netherlands.
		(p. 168)*/

	t1 = 1. - ( LR_STD * g_topo->values[TOPO_ELEV] ) / T_STD;
	t2 = G_STD / ( LR_STD * ( Rgas / AIRMASS ) );
	met[month].d[day].air_pressure = P_STD * pow ( t1, t2 );
}


void Air_density (meteo_t *met, const int day, const int month) {

	/* compute density of air (in kg/m3) */
	/* following Solantie R., 2004, Boreal Environmental Research, 9: 319-333, the model uses tday if available */

	if(met[month].d[day].tday == NO_DATA)
	{
		met[month].d[day].rho_air = 1.292 - (0.00428 * met[month].d[day].tavg);
	}
	else
	{
		met[month].d[day].rho_air = 1.292 - (0.00428 * met[month].d[day].tday);
	}
}

void Latent_heat(meteo_t *met, const int day, const int month)
{
	/*BIOME-BGC APPROACH*/
	/*compute latent heat of vaporization (J/Kg)*/

	met[month].d[day].lh_vap      = 2.5023e6 - 2430.54 * met[month].d[day].tavg;
	met[month].d[day].lh_vap_soil = 2.5023e6 - 2430.54 * met[month].d[day].tsoil;
	/*latent heat of fusion (KJ/Kg)*/
	met[month].d[day].lh_fus      = 335.;
	/*latent heat of sublimation (KJ/Kg)*/
	met[month].d[day].lh_sub      = 2845.;
}
void Soil_temperature(const cell_t *const c, int day, int month, int year)
{
	int i;
	int day_avg = WEIGHTED_DAYS;
	int current_day = day;
	int current_month = month;
	int current_year_index = year;
	double weighted_avg;
	extern int days_per_month[];


	meteo_t *met;
	met = (meteo_t*) c->years[year].m;

	/* compute Tsoil only if Tsoil from met data is not available */

	if ( -9999 == met[month].d[day].ts_f )
	{

		/* following BIOME-bgc 4.2 */
		/* for this version, an 10-day running weighted average of daily
	average temperature is used as the soil temperature at 10 cm.
	For days 1-10, a 1-10 day running weighted average is used instead.
	The tail of the running average is weighted linearly from 1 to 11.
	There are no corrections for snowpack or vegetation cover.
		 */

		//TODO CHECK SOIL TEMPERATURE CORRECTION FROM BIOME
		/* soil temperature correction using difference from annual average tair */
		/*file bgc.c*/
		/* original biome_bgc version */
		/* *
		tdiff = tair_avg - metv.tsoil;
		if (ws.snoww)
		{
			metv.tsoil += 0.83 * tdiff;
		}
		else
		{
			metv.tsoil += 0.2 * tdiff;
		}
		 */


		assert(c);

		i = 0;
		weighted_avg = 0.;
		do
		{
			i += day_avg;
			weighted_avg += ( c->years[year].m[month].d[day].tavg * day_avg );

			if ( --day < 0 )
			{
				if ( --month < 0 )
				{
					if ( --year < 0 )
					{
						break;
					}
					month = 11; // zero based index
				}
				day = days_per_month[month];
				if ( IS_LEAP_YEAR(c->years[year].year) && (1 == month) )
				{
					++day;
				}
				--day; // zero based index
			}
			--day_avg;
		} while ( day_avg > 0 );
		c->years[current_year_index].m[current_month].d[current_day].tsoil = weighted_avg / i;
	}
	else
	{
		c->years[current_year_index].m[current_month].d[current_day].tsoil = met[month].d[day].ts_f;
	}

}

void Weighted_average_temperature(const cell_t *const c, const e_weighted_average_var var, int day, int month, int year)
{
	int i;
	int day_avg = WEIGHTED_DAYS;
	int current_day = day;
	int current_month = month;
	int current_year_index = year;
	double weighted_avg;
	extern int days_per_month[];

	assert(((var >= 0) && (var < WEIGHTED_MEAN_COUNT)) && c);

	i = 0;
	weighted_avg = 0.;
	do
	{
		double v;

		i += day_avg;

		switch ( var ) {
		case WEIGHTED_MEAN_TAVG:
			v = c->years[year].m[month].d[day].tavg;
			break;

		case WEIGHTED_MEAN_TDAY:
			v = c->years[year].m[month].d[day].tday;
			break;

		case WEIGHTED_MEAN_TNIGHT:
			v = c->years[year].m[month].d[day].tnight;
			break;

		case WEIGHTED_MEAN_TSOIL:
			v = c->years[year].m[month].d[day].tsoil;
			break;
		}

		weighted_avg += (v*day_avg);

		if ( --day < 0 ) {
			if ( --month < 0 ) {
				if ( --year < 0 ) {
					break;
				}
				month = 11; // zero based index
			}
			day = days_per_month[month];
			if ( IS_LEAP_YEAR(c->years[year].year) && (1 == month) ) {
				++day;
			}
			--day; // zero based index
		}
		--day_avg;

	} while ( day_avg > 0 );

	switch ( var ) {
	case WEIGHTED_MEAN_TAVG:
		c->years[current_year_index].m[current_month].d[current_day].ten_day_weighted_avg_tavg = weighted_avg / i;
		break;

	case WEIGHTED_MEAN_TDAY:
		c->years[current_year_index].m[current_month].d[current_day].ten_day_weighted_avg_tday = weighted_avg / i;
		break;

	case WEIGHTED_MEAN_TNIGHT:
		c->years[current_year_index].m[current_month].d[current_day].ten_day_weighted_avg_tnight = weighted_avg / i;
		break;

	case WEIGHTED_MEAN_TSOIL:
		c->years[current_year_index].m[current_month].d[current_day].ten_day_weighted_avg_tsoil = weighted_avg / i;
		break;
	}
}

void Averaged_temperature(const cell_t *const c, const e_averaged_var var, int day, int month, int year)
{
	int i;
	int day_avg = AVERAGED_DAYS;
	int current_day = day;
	int current_month = month;
	int current_year_index = year;
	double averaged;
	extern int days_per_month[];

	assert(((var >= 0) && (var < AVERAGED_COUNT)) && c);

	i = 0;
	averaged = 0.;

	do
	{
		double v;

		++i;

		switch ( var ) {
		case AVERAGED_TAVG:
			v = c->years[year].m[month].d[day].tavg;
			break;

		case AVERAGED_TDAY:
			v = c->years[year].m[month].d[day].tday;
			break;

		case AVERAGED_TNIGHT:
			v = c->years[year].m[month].d[day].tnight;
			break;

		case AVERAGED_TSOIL:
			v = c->years[year].m[month].d[day].tsoil;
			break;
		}

		averaged += v;

		if ( --day < 0 ) {
			if ( --month < 0 ) {
				if ( --year < 0 ) {
					break;
				}
				month = 11; // zero based index
			}
			day = days_per_month[month];
			if ( IS_LEAP_YEAR(c->years[year].year) && (1 == month) ) {
				++day;
			}
			--day; // zero based index
		}
		--day_avg;

	} while ( day_avg > 0 );

	switch ( var ) {
	case AVERAGED_TAVG:
		c->years[current_year_index].m[current_month].d[current_day].ten_day_avg_tavg = averaged / i;
		break;

	case AVERAGED_TDAY:
		c->years[current_year_index].m[current_month].d[current_day].ten_day_avg_tday = averaged / i;
		break;

	case AVERAGED_TNIGHT:
		c->years[current_year_index].m[current_month].d[current_day].ten_day_avg_tnight = averaged / i;
		break;

	case AVERAGED_TSOIL:
		c->years[current_year_index].m[current_month].d[current_day].ten_day_avg_tsoil = averaged / i;
		break;
	}
}

void Dew_temperature (meteo_t *const met, const int day, const int month)
{
	/* dew point temperature based on Allen et al., 1998; Bosen, 1958; Murray, 1967 */
	met[month].d[day].tdew = (116.91 + 237.3 * log(met[month].d[day].ea))/(16.78 - log(met[month].d[day].ea));
}

void Daily_Ndeposition (const cell_t *const c, int day, int month, int year)
{
	int doy;

	if ( IS_LEAP_YEAR(c->years[year].year)) doy = 366;
	else doy = 365;

	c->years[year].m[month].d[day].Ndeposition = c->years[year].Ndep / doy;

}


/********************* CLIMATIC VARIABLE FOR REGENERATION ******************/

// 5p7 SAPONARO (con integrazioni/modifiche ddalmo)


#if 1

#if 0

// this void function needs to tbe included in the daily met structure

void Thermic_sum_spring (const cell_t * c, const meteo_daily_t *const meteo_daily, const int day, const int month, const int year) {

   double spring_thermic_sum = 0.;

   // Spring condition (March-May)
   if ( day == 30 && month == MAY) {

   //Adding this chain of vector can export the results in other modules
   c->years[year].m[month].d[day].spring_thermic_sum = c->years[year].m[month].d[day].thermic_sum;

   } else {

   c->years[year].m[month].d[day].thermic_sum = 0.;
   c->years[year].m[month].d[day].spring_thermic_sum = 0.;

   }

   return spring_thermic_sum;
   //printf("Thermicsum = %f\n", meteo_daily->spring_thermic_sum);
}

void Soil_winter_temperature (const cell_t *const c, const meteo_daily_t *const meteo_daily, int day, int month, int year) {


   double winter_soil = 0.;

    //Winter condition (January-March)
    if (month >= JANUARY && month <= MARCH) {

    c->years[year].m[month].d[day].winter_soil = c->years[year].m[month].d[day].tsoil;

    } else {

    winter_soil = 0.;

    }
  //printf("Winter soil temp = %f\n", c->years[year].m[month].d[day].winter_soil);
}


void Seedling_temp (const cell_t *const c, const meteo_daily_t *const meteo_daily , int day, int month, int year) {

  double seedling_temp = 0.;

  //Summer condition (June-August)
  if (month >= JUNE && month <= AUGUST) {

   c->years[year].m[month].d[day].seedling_temp = c->years[year].m[month].d[day].tavg;

   } else {

  seedling_temp = 0.;

  }
  //printf(" Temp seedlings = \t%f\n", meteo_daily->seedling_temp);

}

#endif 

	void Seedling_soil_par (const cell_t *const c, const meteo_daily_t *const meteo_daily, int month, int year) {

  	//double seedling_par = 0.;

  	//Summer condition (June-August)
  	if (month >= JUNE && month <= AUGUST) {

  	//c->years[year].m[month].d[day].seedling_par = meteo_daily->par;

  	c->years[year].seedling_par += meteo_daily->par;

  	} 


}

#endif 