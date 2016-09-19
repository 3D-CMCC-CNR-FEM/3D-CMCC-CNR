/* yos.c */
#include "yos.h"
#include <stdio.h>
#include <stdlib.h>
#include "logger.h"
#include "netcdf.h"
#include "common.h"
#include "constants.h"
#include "settings.h"
#include "g-function.h"
#include <string.h>
#include <assert.h>
#include <math.h>

#define RG_RANGE_MIN         0
#define RG_RANGE_MAX        50
#define TA_RANGE_MIN       -20
#define TA_RANGE_MAX        50
#define TMAX_RANGE_MIN     -20
#define TMAX_RANGE_MAX      50
#define TMIN_RANGE_MIN     -20
#define TMIN_RANGE_MAX      50
#define VPD_RANGE_MIN        0
#define VPD_RANGE_MAX      150
#define PRECIP_RANGE_MIN     0
#define PRECIP_RANGE_MAX   250
#define RH_RANGE_MIN         0
#define RH_RANGE_MAX       100

/* do not change this order */
enum {
	YEAR = 0
	, MONTH
	, DAY
	, RG_F
	, TA_F
	, TMAX
	, TMIN
	, VPD_F
	, TS_F
	, PRECIP
	, SWC
	, NDVI_LAI
	, ET
	, WS_F
	, RH_F

	, MET_COLUMNS_COUNT
};


/*	ALESSIOR: please note that leap years are handled correctly
	so do not change 28 (february) to 29! */
static int days_per_month [] = {
		31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31
};

static const char *sz_month_names[YOS_MONTHS_COUNT] = {
		"JANUARY", "FEBRUARY", "MARCH", "APRIL", "MAY", "JUNE", "JULY"
		, "AUGUST", "SEPTEMBER", "OCTOBER", "NOVEMBER", "DECEMBER"
};

extern logger_t* g_log;
extern settings_t* g_settings;
extern const char sz_err_out_of_memory[];
extern char *g_sz_program_path;
extern char *g_sz_input_path;
extern char *g_sz_co2_conc_file;

/* do not change this order */
static const char *sz_met_columns[MET_COLUMNS_COUNT+2] = {
		"Year"
		, "Month"
		, "n_days"
		, "Rg_f"
		, "Ta_f"
		, "Tmax"
		, "Tmin"
		, "VPD_f"
		, "Ts_f"
		, "Precip"
		, "SWC"
		, "LAI"
		, "ET"
		, "WS_f"
		, "RH_f"

		/* hack */
		, "LAT"
		, "LON"
};

static void timestamp_split(const double value, int *const YYYY, int *const MM, int *const DD) {
	*YYYY = (int)value/10000;
	*MM = ((int)value - *YYYY*10000)/100;
	*DD = (int)value - (*YYYY*10000) - (*MM*100);
}

static double get_co2_conc(const int year, int*const err) {
	char buf[256];
	int _year;
	int flag;
	double co2_conc;
	FILE *f;

	assert(err);

	flag = 0;
	co2_conc = 0.;
	f = NULL;

	*err = 0;

	if ( ! g_sz_co2_conc_file ) { *err = 1; goto quit; }

	if ( g_sz_input_path ) {
		int len = strlen(g_sz_input_path);
		int _flag = (('/' == g_sz_input_path[len-1]) || ('\\' == g_sz_input_path[len-1]));
		sprintf(buf, "%s%s%s", g_sz_input_path, _flag ? "" : FOLDER_DELIMITER, g_sz_co2_conc_file);
		f = fopen(buf, "r");
	} else {
		f = fopen(g_sz_co2_conc_file, "r");
	}
	
	if ( ! f )  { *err = 1; goto quit; }
	while ( fgets(buf, 256, f) ) {
		if ( 2 == sscanf(buf, "%d\t%lf", &_year, &co2_conc) ) {
			if ( year == _year ) {
				flag = 1;
				break;
			}
		}
	}

	if ( ! flag ) *err = 1;

quit:
	if ( f )fclose(f);

	return co2_conc;
}

static void yos_clear(yos_t *const yos) {
	if ( yos ) {
		int i;
		int y;
		yos->year = 0;
		yos->co2_conc = INVALID_VALUE;
		for ( i = 0; i < YOS_MONTHS_COUNT; ++i ) {
			for ( y = 0; y < YOS_DAYS_COUNT; ++y ) {
				yos->m[i].d[y].n_days = INVALID_VALUE;
				yos->m[i].d[y].solar_rad = INVALID_VALUE;
				yos->m[i].d[y].tavg = INVALID_VALUE;
				yos->m[i].d[y].tmax = INVALID_VALUE;
				yos->m[i].d[y].tmin = INVALID_VALUE;
				yos->m[i].d[y].tday = INVALID_VALUE;
				yos->m[i].d[y].tnight = INVALID_VALUE;
				yos->m[i].d[y].vpd = INVALID_VALUE;
				yos->m[i].d[y].ts_f = INVALID_VALUE;
				yos->m[i].d[y].prcp = INVALID_VALUE;
				yos->m[i].d[y].swc = INVALID_VALUE;
				yos->m[i].d[y].ndvi_lai = INVALID_VALUE;
				yos->m[i].d[y].daylength = INVALID_VALUE;
				yos->m[i].d[y].thermic_sum = INVALID_VALUE;
				yos->m[i].d[y].rho_air = INVALID_VALUE;
				yos->m[i].d[y].tsoil = INVALID_VALUE;
				yos->m[i].d[y].et = INVALID_VALUE;
				yos->m[i].d[y].windspeed = INVALID_VALUE;
				yos->m[i].d[y].rh_f = INVALID_VALUE;
				yos->m[i].d[y].lh_vap = INVALID_VALUE;
				yos->m[i].d[y].lh_vap_soil = INVALID_VALUE;
				yos->m[i].d[y].lh_fus = INVALID_VALUE;
				yos->m[i].d[y].lh_sub = INVALID_VALUE;
				yos->m[i].d[y].air_pressure = INVALID_VALUE;
			}
		}
	}
}

static void compute_vpd(double *const values, const int rows_count, const int columns_count) {
#define VALUE_AT(r,c)	((r)+((c)*rows_count))

	int i;
	double value;

	for ( i = 0; i < rows_count; i++ ) {
		double ta = values[VALUE_AT(i, TA_F)];
		double rh = values[VALUE_AT(i, RH_F)];
		value = INVALID_VALUE;

		if ( ! IS_INVALID_VALUE(ta) && ! IS_INVALID_VALUE(rh) ) {
			/* 6.1076 is for hPa */
			value = 6.1076 * exp(17.26938818 * ta / (237.3 + ta));
			value *= (1 - rh / 100.0);
			/* convert NaN to invalid value */
			if ( value != value ) {
				value = INVALID_VALUE;
			}

			/* check range */
			if ( ! IS_INVALID_VALUE(value) ) {
				if ( value < VPD_RANGE_MIN || value > VPD_RANGE_MAX ) {
					value = INVALID_VALUE;
				}
			}
		}
		values[VALUE_AT(i, VPD_F)] = value;
	}

#undef VALUE_AT
}

static void compute_rh(double *const values, const int rows_count, const int columns_count) {
#define VALUE_AT(r,c)	((r)+((c)*rows_count))

	int i;
	double value;
	//double svp, vp;
	double rel_hum;

	double e0max, e0min;
	double es;
	double ea;

	for ( i = 0; i < rows_count; i++ ) {
		double ta = values[VALUE_AT(i, TA_F)];
		double tmax = values[VALUE_AT(i, TMAX)];
		double tmin = values[VALUE_AT(i, TMIN)];
		double vpd = values[VALUE_AT(i, VPD_F)];
		value = INVALID_VALUE;

		if ( ! IS_INVALID_VALUE(ta)
				&& ! IS_INVALID_VALUE(tmax)
				&& ! IS_INVALID_VALUE(tmin)
				&& ! IS_INVALID_VALUE(vpd))
		{

			//			/* see Zhang et al., 2008 Ecological Modelling */
			//			/* saturation vapour pressure at the air temperature T (mbar-hPa) */
			//			svp = 6.1076 * exp((17.26 * ta) / (237.3 + ta));
			//			printf("ta = %g\n", ta);
			//			printf("svp = %f\n", svp);
			//
			//			/* compute vapour pressure */
			//			vp = svp - vpd;
			//			printf("vp = %g\n", vp);
			//			printf("vpd = %g\n", vpd);
			//
			//			/* compute relative humidity */
			//			rel_hum = (ea/es)*100.0;
			//			value = rel_hum;
			//			printf("rel_hum = %g\n", rel_hum);


			/* compute saturation vapour pressure at the maximum and minimum air temperature (hPa) */
			e0max = 6.1076 * exp((17.27*tmax)/(tmax+237.3));
			e0min = 6.1076 * exp((17.27*tmin)/(tmin+237.3));
			//printf("e0max = %g\n", e0max);
			//printf("e0min = %g\n", e0min);

			/* compute weighted mean saturation vapour pressure at the air temperature (hPa) */
			//todo ((e0max*c->ni) + (e0min*(1.0-c->ni)));
			es = (e0max + e0min)/2.0;
			//printf("es = %g\n", es);

			/* compute actual vapour pressure (hPa) */
			ea = es - vpd;
			//printf("ea = %g\n", ea);
			//printf("vpd = %g\n", vpd);

			rel_hum = (ea/es)*100.0;
			value = rel_hum;
			//printf("rel_hum = %g %%\n", rel_hum);

			CHECK_CONDITION(rel_hum, < RH_RANGE_MIN);
			CHECK_CONDITION(rel_hum, > RH_RANGE_MAX);

			/* convert NaN to invalid value */
			if ( value != value ) {
				value = INVALID_VALUE;
			}

			/* check range */
			if ( ! IS_INVALID_VALUE(value) ) {
				if ( value < RH_RANGE_MIN || value > RH_RANGE_MAX ) {
					value = INVALID_VALUE;
				}
			}
		}
		values[VALUE_AT(i, RH_F)] = value;
	}
#undef VALUE_AT
}

static int yos_from_arr(double *const values, const int rows_count, const int columns_count, yos_t** p_yos, int *const yos_count) {
#define VALUE_AT(r,c)	((r)+((c)*rows_count))
	yos_t *yos_no_leak;
	yos_t *yos;
	int row;
	int year;
	int month;
	int day;
	int current_year;

	static double previous_solar_rad,
	previous_tavg,
	previous_tmax,
	previous_tmin,
	previous_vpd,
	previous_ts_f,
	previous_prcp,
	previous_swc,
	previous_ndvi_lai;

	assert(p_yos && yos_count);

	year = 0;
	current_year = -1;
	yos = *p_yos;

	/* check imported values */
	{
		int row;
		int column;
		int *flag;
		/* remove date, month and year from columns count */
		flag = malloc((columns_count-3)*sizeof*flag);
		if ( ! flag ) {
			logger(g_log, sz_err_out_of_memory);
			free(yos);
			return 0;
		}

		/* skip day, month and year */
		for ( column = 0; column < columns_count-3; ++column ) {
			flag[column] = 0;
			for ( row = 0; row < rows_count; ++row ) {
				if ( IS_INVALID_VALUE(values[VALUE_AT(row, column+3)]) ) {
					++flag[column];
				}
			}
			if ( rows_count == flag[column]  ) {
				flag[column] = 1; /* column is entirely invalid! */
			} else {
				flag[column] = 0;
			}
		}

		/* check for missing vars */
		if ( flag[VPD_F-3] && flag[RH_F-3] ) {
			logger(g_log, "VPD and RH columns are missing!\n");
			free(flag);
			free(yos);
			return 0;
		}

		if ( flag[TA_F-3] && flag[TMIN-3] && flag[TMAX-3] ) {
			logger(g_log, "TA, TMIN and TMAX columns are missing!\n");
			free(flag);
			free(yos);
			return 0;
		}

		/* compute ta */
		if ( flag[TA_F-3] ) {
			for ( row = 0; row < rows_count; ++row ) {
				if ( ! IS_INVALID_VALUE(values[VALUE_AT(row, TMAX)])
					&& ! IS_INVALID_VALUE(values[VALUE_AT(row, TMIN)]) ) {
						values[VALUE_AT(row, TA_F)] = (0.606 * values[VALUE_AT(row, TMAX)]) + (0.394 * values[VALUE_AT(row, TMIN)]);
				}
			}
		}

		/* ALESSIOR: fix rh 'cause can be imported out of range */
		if ( ! flag[RH_F-3] ) {
			for ( row = 0; row < rows_count; ++row ) {
				if ( ! IS_INVALID_VALUE(values[VALUE_AT(row, RH_F)]) ) {
					if ( values[VALUE_AT(row, RH_F)] < RH_RANGE_MIN ) values[VALUE_AT(row, RH_F)] = RH_RANGE_MIN;
					else if ( values[VALUE_AT(row, RH_F)] > RH_RANGE_MAX ) values[VALUE_AT(row, RH_F)] = RH_RANGE_MAX;
				}
			}
		}

		/* compute vpd ? or rh ?*/
		/* please note that we must have TA, so computing is done after computing TA (if needed) */
		if ( flag[VPD_F-3] ) {
			compute_vpd(values, rows_count, columns_count);
		} else if ( flag[RH_F-3] ) {
			compute_rh(values, rows_count, columns_count);
		}

		free(flag);
	}

	for ( row = 0; row < rows_count; ++row ) {
		year = (int)values[VALUE_AT(row, YEAR)];
		if ( ! year ) {
			puts("year cannot be zero!\n");
			logger(g_log, "year cannot be zero!\n");
			free(yos);
			return 0;
		}

		month = (int)values[VALUE_AT(row, MONTH)];
		if ( month < 1 || month > YOS_MONTHS_COUNT ) {
			printf("bad month for year %d\n\n", year);
			logger(g_log, "bad month for year %d\n\n", year);
			free(yos);
			return 0;

		}
		--month;

		day = (int)values[VALUE_AT(row, DAY)];
		if ( (day <= 0) || day > days_per_month[month] + (((1 == month) && IS_LEAP_YEAR(year)) ? 1 :0 ) ) {
			printf("bad day for %s %d\n\n", sz_month_names[month], year);
			logger(g_log, "bad day for %s %d\n\n", sz_month_names[month], year);
			free(yos);
			return 0;
		}
		--day;

		if ( current_year != year ) {
			yos_no_leak = realloc(yos, (*yos_count+1)*sizeof*yos_no_leak);
			if ( ! yos_no_leak )
			{
				puts(sz_err_out_of_memory);
				free(yos);
				return 0;
			}
			yos = yos_no_leak;
			yos_clear(&yos[*yos_count]);

			yos[*yos_count].year = year;
			++*yos_count;

			current_year = year;
		}

		yos[*yos_count-1].m[month].d[day].n_days = day+1;
		if (yos[*yos_count-1].m[month].d[day].n_days > YOS_DAYS_COUNT)
		{
			logger(g_log, "ERROR IN N_DAYS DATA!!\n");
			free(yos);
			return 0;
		}

		/* case RG_F: //Rg_f - solar_rad -daily average solar radiation */
		yos[*yos_count-1].m[month].d[day].solar_rad = values[VALUE_AT(row,RG_F)];
		if ( IS_INVALID_VALUE (yos[*yos_count-1].m[month].d[day].solar_rad) && (!((day == 0) && (1 == *yos_count) && (month == 0))))
		{

			//the model gets the value of the day before
			//Log ("* SOLAR RAD -NO DATA in year %d month %s, day %d!!!!\n", yos[*yos_count-1].year, MonthName[month], day);
			//logger(g_log, "Getting previous day values.. !!\n");
			yos[*yos_count-1].m[month].d[day].solar_rad = previous_solar_rad;
			//logger(g_log, "..value of the previous day = %f\n", yos[*yos_count-1].m[month].d[day].solar_rad);
			if ( IS_INVALID_VALUE (yos[*yos_count-1].m[month].d[day].solar_rad))
			{
				//Log ("********* SOLAR RAD -NO DATA- in previous day!!!!\n" );

				//the model gets the value of the year before
				if ( *yos_count > 1 ) {
					yos[*yos_count-1].m[month].d[day].solar_rad = yos[*yos_count-2].m[month].d[day].solar_rad;
					if (IS_INVALID_VALUE (yos[*yos_count-1].m[month].d[day].solar_rad))
					{
						//Log ("********* SOLAR RAD -NO DATA- in previous year!!!!\n" );
						yos[*yos_count-1].m[month].d[day].solar_rad = NO_DATA;
					}
				} else {
					yos[*yos_count-1].m[month].d[day].solar_rad = yos[*yos_count-2].m[month].d[day-1].solar_rad;
					/*
					if (IS_INVALID_VALUE (yos[*yos_count-1].m[month].d[day-1].solar_rad))
					{
						Log ("********* SOLAR RAD -NO DATA- in previous day!!!!\n" );
						exit(1);
					}
					 */
				}
			}
		}
		/* check if values are outside ranges */
		else if(yos[*yos_count-1].m[month].d[day].solar_rad < RG_RANGE_MIN || yos[*yos_count-1].m[month].d[day].solar_rad > RG_RANGE_MAX)
		{
			logger(g_log, "BAD DATA FOR RG = %f in day = %d month = %d year = %d\n", yos[*yos_count-1].m[month].d[day].solar_rad, day+1, month+1, year);
			exit(1);
		}
		else
		{
			previous_solar_rad = yos[*yos_count-1].m[month].d[day].solar_rad;
		}

		/* case TA_F: //Ta_f -  temperature average */
		yos[*yos_count-1].m[month].d[day].tavg = values[VALUE_AT(row,TA_F)];
		if (IS_INVALID_VALUE (yos[*yos_count-1].m[month].d[day].tavg) && (!((day == 0) && (*yos_count == 1) && (month == 0))))
		{
			//the model gets the value of the day before
			//Log ("* TAVG -NO DATA in year %d month %s, day %d!!!!\n", yos[*yos_count-1].year, MonthName[month], day);

			yos[*yos_count-1].m[month].d[day].tavg = previous_tavg;
			if ( IS_INVALID_VALUE (yos[*yos_count-1].m[month].d[day].tavg))
			{
				//the model gets the value of the year before
				yos[*yos_count-1].m[month].d[day].tavg = yos[*yos_count-2].m[month].d[day].tavg;

				if (IS_INVALID_VALUE (yos[*yos_count-1].m[month].d[day].tavg))
				{
					//Log ("********* TAVG -NO DATA- in previous year!!!!\n" );
					yos[*yos_count-1].m[month].d[day].tavg = NO_DATA;
				}
			}
		}
		/* check if values are outside ranges */
		else if(yos[*yos_count-1].m[month].d[day].tavg < TA_RANGE_MIN || yos[*yos_count-1].m[month].d[day].tavg > TA_RANGE_MAX)
		{
			logger(g_log, "BAD DATA FOR Tavg = %f in day = %d month = %d year = %d\n", yos[*yos_count-1].m[month].d[day].tavg, day+1, month+1, year);
			exit(1);
		}
		else
		{
			previous_tavg = yos[*yos_count-1].m[month].d[day].tavg;
		}

		/* case TMAX: //TMAX -  maximum temperature */
		yos[*yos_count-1].m[month].d[day].tmax = values[VALUE_AT(row,TMAX)];
		if ( IS_INVALID_VALUE (yos[*yos_count-1].m[month].d[day].tmax) && (!((day == 0) && (*yos_count == 1)&& (month == 0))))
		{
			//the model gets the value of the day before
			//Log ("* TMAX -NO DATA in year %d month %s, day %d!!!!\n", yos[*yos_count-1].year, MonthName[month], day);
			//logger(g_log, "Getting previous day values.. !!\n");
			yos[*yos_count-1].m[month].d[day].tmax = yos[*yos_count-1].m[month].d[day].tavg;
			//logger(g_log, "..using tavg = %f\n", yos[*yos_count-1].m[month].d[day].tavg);
			if ( IS_INVALID_VALUE (yos[*yos_count-1].m[month].d[day].tmax))
			{
				//the model gets the value of the year before
				yos[*yos_count-1].m[month].d[day].tmax = yos[*yos_count-2].m[month].d[day].tmax;

				if (IS_INVALID_VALUE (yos[*yos_count-1].m[month].d[day].tmax))
				{
					//Log ("********* TMAX -NO DATA- in previous year!!!!\n" );
					yos[*yos_count-1].m[month].d[day].tmax = NO_DATA;
				}
			}
		}
		/* check if values are outside ranges */
		else if(yos[*yos_count-1].m[month].d[day].tmax < TMAX_RANGE_MIN || yos[*yos_count-1].m[month].d[day].tmax > TMAX_RANGE_MAX)
		{
			logger(g_log, "BAD DATA FOR Tmax = %f in day = %d month = %d year = %d\n", yos[*yos_count-1].m[month].d[day].tmax, day+1, month+1, year);
			exit(1);
		}
		else
		{
			previous_tmax = yos[*yos_count-1].m[month].d[day].tmax;
		}

		//case TMIN: //TMIN -  minimum temperature
		yos[*yos_count-1].m[month].d[day].tmin = values[VALUE_AT(row,TMIN)];
		if ( IS_INVALID_VALUE (yos[*yos_count-1].m[month].d[day].tmin) && (!((day == 0) && (*yos_count == 1)&& (month == 0))))
		{
			//the model gets the value of the day before
			//Log ("* TMIN -NO DATA in year %d month %s, day %d!!!!\n", yos[*yos_count-1].year, MonthName[month], day);
			//logger(g_log, "Getting previous day values.. !!\n");
			yos[*yos_count-1].m[month].d[day].tmin = yos[*yos_count-1].m[month].d[day].tavg;
			//logger(g_log, "..using tavg = %f\n", yos[*yos_count-1].m[month].d[day].tavg);
			if ( IS_INVALID_VALUE (yos[*yos_count-1].m[month].d[day].tmin))
			{
				//the model gets the value of the year before
				yos[*yos_count-1].m[month].d[day].tmin = yos[*yos_count-2].m[month].d[day].tmin;

				if (IS_INVALID_VALUE (yos[*yos_count-1].m[month].d[day].tmin))
				{
					//Log ("********* TMIN -NO DATA- in previous year!!!!\n" );
					yos[*yos_count-1].m[month].d[day].tmin = NO_DATA;
				}
			}
		}
		/* check if values are outside ranges */
		else if(yos[*yos_count-1].m[month].d[day].tmin < TMIN_RANGE_MIN || yos[*yos_count-1].m[month].d[day].tmin > TMIN_RANGE_MAX)
		{
			logger(g_log, "BAD DATA FOR Tmin = %f in day = %d month = %d year = %d\n", yos[*yos_count-1].m[month].d[day].tmin, day+1, month+1, year);
			exit(1);
		}
		else
		{
			previous_tmin = yos[*yos_count-1].m[month].d[day].tmin;
		}

		/* case VPD_F: //RH_f - RH */
		yos[*yos_count-1].m[month].d[day].vpd = values[VALUE_AT(row,VPD_F)];
		if ( IS_INVALID_VALUE (yos[*yos_count-1].m[month].d[day].vpd) && (!((day == 0) && (*yos_count == 1)&& (month == 0))))
		{
			//the model gets the value of the day before
			//Log ("* VPD -NO DATA in year %d month %s, day %d!!!!\n", yos[*yos_count-1].year, MonthName[month], day);
			//logger(g_log, "Getting previous day values.. !!\n");
			yos[*yos_count-1].m[month].d[day].vpd = previous_vpd;
			//logger(g_log, "..value of the previous day = %f\n", yos[*yos_count-1].m[month].d[day].vpd);
			if ( IS_INVALID_VALUE (yos[*yos_count-1].m[month].d[day].vpd))
			{
				//Log ("********* VPD -NO DATA- in previous year!!!!\n" );

				//the model gets the value of the year before
				yos[*yos_count-1].m[month].d[day].vpd = yos[*yos_count-2].m[month].d[day].vpd;

				if (IS_INVALID_VALUE (yos[*yos_count-1].m[month].d[day].vpd))
				{
					//Log ("********* VPD -NO DATA- in previous year!!!!\n" );
					yos[*yos_count-1].m[month].d[day].vpd = NO_DATA;
				}

			}
		}
		/* check if values are outside ranges */
		else if(yos[*yos_count-1].m[month].d[day].vpd < VPD_RANGE_MIN || yos[*yos_count-1].m[month].d[day].vpd > VPD_RANGE_MAX)
		{
			logger(g_log, "BAD DATA FOR vpd = %f in day = %d month = %d year = %d\n", yos[*yos_count-1].m[month].d[day].vpd, day+1, month+1, year);
			exit(1);
		}
		else
		{
			previous_vpd = yos[*yos_count-1].m[month].d[day].vpd;
		}
		//logger(g_log, "%d-%s-vpd = %f\n",yos[*yos_count-1].m[month].d[day].n_days, MonthName[month], yos[*yos_count-1].m[month].d[day].vpd);

		/* case TS_F: // ts_f   Soil temperature */
		yos[*yos_count-1].m[month].d[day].ts_f = values[VALUE_AT(row,TS_F)];
		if ( IS_INVALID_VALUE (yos[*yos_count-1].m[month].d[day].ts_f) && (!((day == 0) && (*yos_count == 1)&& (month == 0))))
		{
			//the model gets the value of the day before
			//Log ("* TS_F -NO DATA in year %s month %s!!!!\n", year, MonthName[month] );
			yos[*yos_count-1].m[month].d[day].ts_f = previous_ts_f;
			/*
			if ( IS_INVALID_VALUE (yos[*yos_count-1].m[month].d[day].ts_f))
			{
				//Log ("********* TS_F -NO DATA- in previous year!!!!\n" );

				//the model gets the value of the year before
				yos[*yos_count-1].m[month].d[day].ts_f = yos[*yos_count-1].m[month].d[day].ts_f;

				if (IS_INVALID_VALUE (yos[*yos_count-1].m[month].d[day].ts_f))
				{
					//Log ("********* TS_F -NO DATA- in previous year!!!!\n" );
					yos[*yos_count-1].m[month].d[day].ts_f = NO_DATA;
				}
			}
			 */
		}
		else
		{
			previous_ts_f = yos[*yos_count-1].m[month].d[day].ts_f;
		}

		/* case PRECIP:  //Precip - rain */
		yos[*yos_count-1].m[month].d[day].prcp = values[VALUE_AT(row,PRECIP)];
		if ( IS_INVALID_VALUE (yos[*yos_count-1].m[month].d[day].prcp) && (!((day == 0) && (*yos_count == 1)&& (month == 0))))
		{
			//the model gets the value of the day before
			//Log ("* PRECIPITATION -NO DATA in year %d month %s, day %d!!!!\n", yos[*yos_count-1].year, MonthName[month], day+1);
			//logger(g_log, "Getting previous day values.. !!\n");
			yos[*yos_count-1].m[month].d[day].prcp = previous_prcp;
			//logger(g_log, "..value of the previous day = %f\n", yos[*yos_count-1].m[month].d[day].prcp);
			if ( IS_INVALID_VALUE (yos[*yos_count-1].m[month].d[day].prcp))
			{
				logger(g_log, "********* PRECIPITATION -NO DATA- in previous year!!!!\n" );

				//the model gets the value of the year before
				yos[*yos_count-1].m[month].d[day].prcp = yos[*yos_count-2].m[month].d[day].prcp;

				if (IS_INVALID_VALUE (yos[*yos_count-1].m[month].d[day].prcp))
				{
					logger(g_log, "********* RAIN -NO DATA- in previous year!!!!\n" );
					yos[*yos_count-1].m[month].d[day].prcp = NO_DATA;
				}

			}
			//logger(g_log, "precipitation of previous year = %f mm\n", yos[*yos_count-1].m[month].rain);
		}
		/* check if values are outside ranges */
		else if(yos[*yos_count-1].m[month].d[day].prcp < PRECIP_RANGE_MIN || yos[*yos_count-1].m[month].d[day].prcp > PRECIP_RANGE_MAX)
		{
			logger(g_log, "BAD DATA FOR prcp = %f in day = %d month = %d year = %d\n", yos[*yos_count-1].m[month].d[day].prcp, day+1, month+1, year);
			exit(1);
		}
		else
		{
			previous_prcp = yos[*yos_count-1].m[month].d[day].prcp;
		}

		/* case SWC: //Soil Water Content (%) */
		yos[*yos_count-1].m[month].d[day].swc = values[VALUE_AT(row,SWC)];
		/*if ( error_flag )
		{
			printf("unable to convert value \"%s\" at column %d for %s day %d\n", token2, column+1, MonthName[month], day);
			logger(g_log, "unable to convert value \"%s\" at column %d for %s day %d\n", token2, column+1, MonthName[month], day);
			free(yos);
			fclose(f);
			return 0;
		}*/
		if ( IS_INVALID_VALUE (yos[*yos_count-1].m[month].d[day].swc) && (!((day == 0) && (*yos_count == 1)&& (month == 0))))
		{
			//the model gets the value of the day before
			//Log ("********* SWC -NO DATA in year %s month %s!!!!\n", year, MonthName[month] );
			//logger(g_log, "Getting previous years values !!\n");
			yos[*yos_count-1].m[month].d[day].swc = previous_swc;
			/*
			if ( IS_INVALID_VALUE (yos[*yos_count-1].m[month].d[day].swc))
			{
				//Log ("* SWC -NO DATA- in previous year!!!!\n" );
				//the model gets the value of the year before
				yos[*yos_count-1].m[month].d[day].swc = yos[*yos_count-1].m[month].d[day].swc;

				if (IS_INVALID_VALUE (yos[*yos_count-1].m[month].d[day].swc))
				{
					//Log ("********* RAIN -NO DATA- in previous year!!!!\n" );
					yos[*yos_count-1].m[month].d[day].swc = NO_DATA;
				}
			}
			 */
		}
		else
		{
			previous_swc = yos[*yos_count-1].m[month].d[day].swc;
		}

		//logger(g_log, "%d-%s-swc= %f\n",yos[*yos_count-1].m[month].d[day].n_days, MonthName[month], yos[*yos_count-1].m[month].d[day].swc);
		//break;
		//case NDVI_LAI: //Get LAI in spatial version
		//todo to remove
		if ( 's' == g_settings->spatial )
		{
			yos[*yos_count-1].m[month].d[day].ndvi_lai = values[VALUE_AT(row,NDVI_LAI)];
			//if is not the first year the model get the previous year value
			if (*yos_count > 1)
			{

				//control in lai data if is an invalid value
				if ( IS_INVALID_VALUE (yos[*yos_count-1].m[month].d[day].ndvi_lai) && (!((day == 0) && (month == 0))))
				{
					//the model gets the value of the day before
					//Log ("********* LAI -NO DATA in year %d month %s, day %d!!!!\n", yos[*yos_count-1].year, MonthName[month], day+1 );
					//logger(g_log, "Getting previous years values !!\n");
					yos[*yos_count-1].m[month].d[day].ndvi_lai = previous_ndvi_lai;
					if ( IS_INVALID_VALUE (yos[*yos_count-1].m[month].d[day].ndvi_lai))
					{
						//Log ("* LAI -NO DATA- in previous year!!!!\n" );
						yos[*yos_count-1].m[month].d[day].ndvi_lai = NO_DATA;
					}
				}
				else
				{
					previous_ndvi_lai = yos[*yos_count-1].m[month].d[day].ndvi_lai;
				}
			}
		}

		//case ET: //ET for stand alone RothC (todo remove)
		yos[*yos_count-1].m[month].d[day].et = values[VALUE_AT(row,ET)];
		if ( IS_INVALID_VALUE (yos[*yos_count-1].m[month].d[day].et) && (!((day == 0) && (*yos_count == 1)&& (month == 0))))
		{
			//the model gets the value of the day before
			//Log ("* ET -NO DATA in year %s month %s, day %d!!!!\n", year, MonthName[month], day);
		}
		//logger(g_log, "%d-%s-tavg = %f\n",yos[*yos_count-1].m[month].d[day].n_days, MonthName[month], yos[*yos_count-1].m[month].d[day].tavg);

		//case WS_F: //windspeed
		yos[*yos_count-1].m[month].d[day].windspeed = values[VALUE_AT(row,WS_F)];
		if ( IS_INVALID_VALUE (yos[*yos_count-1].m[month].d[day].windspeed) && (!((day == 0) && (*yos_count == 1)&& (month == 0))))
		{
			//the model gets the value of the day before
			//Log ("* windspeed -NO DATA in year %s month %s, day %d!!!!\n", yos[*yos_count-1].year, MonthName[month], day);
		}
		//logger(g_log, "%d-%s-tavg = %f\n",yos[*yos_count-1].m[month].d[day].n_days, MonthName[month], yos[*yos_count-1].m[month].d[day].tavg);

		// RH_f
		yos[*yos_count-1].m[month].d[day].rh_f = values[VALUE_AT(row,RH_F)];
		if ( IS_INVALID_VALUE (yos[*yos_count-1].m[month].d[day].rh_f) && (!((day == 0) && (*yos_count == 1)&& (month == 0))))
		{
			// ALESSIOR
			// TODO: ask ALESSIOC
		}
		// ALESSIOR:
		// fixed inside yos_from_arr, so compute vpd is safe!
		/*
		else if(yos[*yos_count-1].m[month].d[day].rh_f < RH_RANGE_MIN || yos[*yos_count-1].m[month].d[day].rh_f > RH_RANGE_MAX)
		{
			logger(g_log, "BAD DATA FOR rh = %g in day = %d month = %d year = %d\n", yos[*yos_count-1].m[month].d[day].rh_f, day+1, month+1, year);

			//fixme
			if(yos[*yos_count-1].m[month].d[day].rh_f < RH_RANGE_MIN)
			{
				logger(g_log, "WARNING!!!!! forced RH to 0%%\a\n");
				yos[*yos_count-1].m[month].d[day].rh_f = 0.0;
			}
			if(yos[*yos_count-1].m[month].d[day].rh_f > RH_RANGE_MAX)
			{
				logger(g_log, "WARNING!!!!! forced RH to 100%%\a\n");
				yos[*yos_count-1].m[month].d[day].rh_f = 100.0;
			}

			//fixme
			//exit(1);
		}
		*/
	}
	*p_yos = yos;
	return 1;
#undef VALUE_AT
}

static int import_nc(const char* const filename, yos_t** pyos, int* const yos_count) {
#define COLUMN_AT(c)	((c)*dims_size[ROWS_DIM])
#define VALUE_AT(r,c)	((r)+(COLUMN_AT((c))))
	int i;
	int y;
	int z;
	int id_file;
	int ret;

	int dims_count;	/* dimensions */
	int vars_count;
	int atts_count;	/* attributes */
	int unl_count;	/* unlimited dimensions */
	char name[NC_MAX_NAME+1];
	nc_type type;
	size_t size;
	int *i_values = NULL; /* required */
	double *values = NULL; /* required */
	int columns[MET_COLUMNS_COUNT];

	/* DO NOT CHANGE THIS ORDER */
	enum {
		ROWS_DIM,
		X_DIM,
		Y_DIM,

		DIMS_COUNT
	};

	int dims_size[DIMS_COUNT];
	const char *sz_dims[DIMS_COUNT] = { "row", "x", "y" }; /* DO NOT CHANGE THIS ORDER...please see top */

	/* */
	ret = nc_open(filename, NC_NOWRITE, &id_file);
	if ( ret != NC_NOERR ) goto quit;

	ret = nc_inq(id_file, &dims_count, &vars_count, &atts_count, &unl_count);
	if ( ret != NC_NOERR ) goto quit;

	if ( ! dims_count || ! vars_count ) {
		printf("bad nc file! %d dimensions and %d vars\n\n", dims_count, vars_count);
		nc_close(id_file);
		return 0;
	}

	/* check if vars count are at least MET_COLUMNS+2 */
	if ( vars_count < MET_COLUMNS_COUNT+2 ) {
		printf("bad nc file! Vars count should be %d at least.\n\n", vars_count);
		nc_close(id_file);
		return 0;
	}

	/* reset */
	for ( i = 0; i < MET_COLUMNS_COUNT; ++i ) {
		columns[i] = 0;
	}
	for ( i = 0; i < DIMS_COUNT; ++i ) {
		dims_size[i] = -1;
	}

	/* get dimensions */
	for ( i = 0; i < dims_count; ++i ) {
		ret = nc_inq_dim(id_file, i, name, &size);
		if ( ret != NC_NOERR ) goto quit;
		for ( y = 0; y < DIMS_COUNT; ++y ) {
			if ( ! string_compare_i(sz_dims[y], name) ) {
				if ( dims_size[y] != -1 ) {
					printf("dimension %s already found!\n", sz_dims[y]);
					nc_close(id_file);
					return 0;
				}
				dims_size[y] = size;
				break;
			}
		}
	}

	/* check if we have all dimensions */
	for ( i = 0; i < DIMS_COUNT; ++i ) {
		if ( -1 == dims_size[i] ) {
			printf("dimension %s not found!\n", sz_dims[i]);
			nc_close(id_file);
			return 0;
		}
	}

	/* check x and y */
	if( (dims_size[X_DIM] != 1) || (dims_size[Y_DIM] != 1) ) {
		puts("x and y must be 1!");
		nc_close(id_file);
		return 0;
	}

	/* alloc memory for int values */
	i_values = malloc(dims_size[ROWS_DIM]*sizeof*values);
	if ( ! i_values ) {
		puts(sz_err_out_of_memory);
		nc_close(id_file);
		return 0;
	}

	/* alloc memory for double values */
	/* please note that we alloc double for year,month,day too */
	values = malloc(dims_size[ROWS_DIM]*MET_COLUMNS_COUNT*sizeof*values);
	if ( ! values ) {
		puts(sz_err_out_of_memory);
		free(i_values);
		nc_close(id_file);
		return 0;
	}

	/* set all double values to -9999 */
	for ( i = 0; i < dims_size[ROWS_DIM]*MET_COLUMNS_COUNT; ++i ) {
		values[i] = INVALID_VALUE;
	}

	/* get vars */
	for ( i = 0; i < vars_count; ++i ) {
		ret = nc_inq_var(id_file, i, name, &type, NULL, NULL, NULL);
		if ( ret != NC_NOERR ) goto quit;
		/* check if we need that var */
		for ( y = 0; y  < MET_COLUMNS_COUNT; ++y ) {
			if ( ! string_compare_i(name, sz_met_columns[y]) ) {
				/* check if we've already get that var */
				if ( columns[y] ) {
					logger(g_log, "column %s already imported!", name);
					free(values);
					free(i_values);
					nc_close(id_file);
					return 0;
				} else {
					columns[y] = 1;
				}
				/* get values */
				if ( NC_DOUBLE == type ) {
					ret = nc_get_var_double(id_file, i, &values[COLUMN_AT(y)]);
					if ( ret != NC_NOERR ) goto quit;
				} else if ( NC_INT == type ) {
					ret = nc_get_var_int(id_file, i, i_values);
					if ( ret != NC_NOERR ) goto quit;
					for ( z = 0; z < dims_size[ROWS_DIM]; ++z ) {
						values[VALUE_AT(z, y)] = (double)i_values[z];
					}
				} else {
					/* type format not supported! */
					logger(g_log, "type format for %s column not supported", name);
					free(values);
					free(i_values);
					nc_close(id_file);
					return 0;
				}

				break;
			}
		}
	}
	free(i_values);
	i_values = NULL;
#if 0
	/* check if we've all needed vars */
	for ( i = 0; i < MET_COLUMNS_COUNT; ++i ) {
		if ( ((VPD_F == i) && (-1 == columns[RH_F])) || ((RH_F == i) && (-1 == columns[VPD_F])) ) {
			logger(g_log, "met columns %s and %s are missing!\n\n", sz_met_columns[VPD_F], sz_met_columns[RH_F]);
			free(values);
			free(i_values);
			nc_close(id_file);
			return 0;
		} else if ( -1 == columns[i] ) {
			logger(g_log, "met column %s not found.\n\n", sz_met_columns[i]);
			free(values);
			free(i_values);
			nc_close(id_file);
			return 0;
		}
	}

	/* compute vpd ?*/
	if ( -1 == columns[VPD_F] ) {
		compute_vpd(values, dims_size[ROWS_DIM], MET_COLUMNS_COUNT);
	}
#endif

	/* save file */
#if 0
	{
		FILE *f;
		int row;
		f = fopen("debug_file", "w");
		if ( ! f ) {
			puts("unable to create output file!");
			nc_close(id_file);
			free(values);
			free(i_values);
			return 0;
		}
		/* write header */
		fputs("Year\tMonth\tn_days\tRg_f\tTa_f\tTmax\tTmin\tVPD_f\tTs_f\tPrecip\tSWC\tLAI\tET\tWS_F\n", f);
		for ( row = 0; row < dims_size[ROWS_DIM]; ++row ) {
			fprintf(f, "%d\t%d\t%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
					/*
						(int)values[VALUE_AT(0,0,row,YEAR)]
						, (int)values[VALUE_AT(0,0,row,MONTH)]
						, (int)values[VALUE_AT(0,0,row,DAY)]
						, values[VALUE_AT(0,0,row,RG_F)]
						, values[VALUE_AT(0,0,row,TA_F)]
						, values[VALUE_AT(0,0,row,TMAX)]
						, values[VALUE_AT(0,0,row,TMIN)]
						, values[VALUE_AT(0,0,row,VPD_F)]
						, values[VALUE_AT(0,0,row,TS_F)]
						, values[VALUE_AT(0,0,row,PRECIP)]
						, values[VALUE_AT(0,0,row,NDVI_LAI)]
						, values[VALUE_AT(0,0,row,ET)]
						, values[VALUE_AT(0,0,row,WS)]
					 */
					(int)values[VALUE_AT(row,YEAR)]
					, (int)values[VALUE_AT(row,MONTH)]
					, (int)values[VALUE_AT(row,DAY)]
					, values[VALUE_AT(row,RG_F)]
					, values[VALUE_AT(row,TA_F)]
					, values[VALUE_AT(row,TMAX)]
					, values[VALUE_AT(row,TMIN)]
					, values[VALUE_AT(row,VPD_F)]
					, values[VALUE_AT(row,TS_F)]
					, values[VALUE_AT(row,PRECIP)]
					, values[VALUE_AT(row,SWC)]
					, values[VALUE_AT(row,NDVI_LAI)]
					, values[VALUE_AT(row,ET)]
					/* ALESSIOC */
					, values[VALUE_AT(row,WS_F)]
			);
		}
		fclose(f);
	}
#endif
	if ( ! yos_from_arr(values, dims_size[ROWS_DIM], MET_COLUMNS_COUNT, pyos, yos_count) ) {
		free(values);
		return 0;
	}
	/* hack */
	ret = 0;

	quit:
	nc_close(id_file);
	free(values);
	free(i_values);
	if ( ret ) {
		puts(nc_strerror(ret));
		ret = 1;
	}
	return ! ret;

#undef COLUMN_AT
#undef VALUE_AT
}


static int import_lst(const char *const filename, yos_t** p_yos, int *const yos_count, const int x_cell, const int y_cell) {
#define VARS_COUNT		((MET_COLUMNS_COUNT)-3)	/* we remove first 3 columns: year, month and day */
#define COLUMN_AT(c)	((c)*rows_count)
#define VALUE_AT(r,c)	((r)+(COLUMN_AT(c)))

	/* TODO: fix sz_path and sz_nc_filename size...can't be same of buffer */
	char buffer[256];
	char sz_path[256];
	char sz_nc_filename[256];
	char *p;
	char *p2;
	int i;
	int z;
	FILE *f;
	int rows_count;
	int vars[VARS_COUNT];
	float *f_values;
	double *values;
	enum {
		X_DIM,
		Y_DIM,
		TIME_DIM,
		HEIGHT_DIM,

		DIMS_COUNT
	};
	int dims_id[DIMS_COUNT];
	int dims_size[DIMS_COUNT];
	const char *sz_lat = "lat";
	const char *sz_lon = "lon";
	const char *sz_time = "time";
	const char *sz_dims[DIMS_COUNT] = { "x", "y", "time", "height_2m" }; /* DO NOT CHANGE THIS ORDER...please see above */
	const char *sz_vars[VARS_COUNT] = { "RADS"
										, "T_2M"
										, "TMAX_2M"
										, "TMIN_2M"
										, "VPD"
										, "TSOIL"
										, "TOT_PREC"
										, "SWC"
										, "LAI"
										, "ET"
										, "WS_F"
										, "RH"
	};

	int y;
	int id_file;
	int ret;
	int dims_count;	/* dimensions */
	int vars_count;
	int atts_count;	/* attributes */
	int unl_count;	/* unlimited dimensions */
	char name[NC_MAX_NAME+1];
	nc_type type;
	size_t size;
	int n_dims;
	int ids[NC_MAX_VAR_DIMS];
	int flag;
	int date_imported;
	float lat;
	float lon;

	/* init */
	rows_count = 0;
	values = NULL;
	f_values = NULL;
	date_imported = 0;
	for ( i = 0; i < VARS_COUNT; i++ ) {
		vars[i] = 0;
	}

	/* open lst file */
	f = fopen(filename, "r");
	if ( ! f ) {
		logger(g_log, "unable to open met data file !\n");
		return 0;
	}

	/* get filename path */
	sz_path[0] = '\0';
	if ( filename[1] != ':' ) {
		strncpy(sz_path, g_sz_program_path, 256);
	}
	/* check for buffer overflow */
	strcat(sz_path, filename);
	p = (strrchr(sz_path, '\\'));
	if ( p ) ++p;
	p2 = (strrchr(sz_path, '/'));
	if ( p2 ) ++p2;
	if  ( p2 > p ) p = p2;
	*p = '\0';

	/* parse lst file */
	while ( fgets(buffer, 256, f) ) {
		/* remove \r\n and skip blank line */
		for ( i = 0; buffer[i]; ++i ) {
			if ( '\r' == buffer[i] || '\n' == buffer[i] ) {
				buffer[i] = '\0';
				break;
			}
		}
		if ( ! buffer[0] ) {
			continue;
		}

		/* TODO: check for buffer overflow */
		sprintf(sz_nc_filename, "%s%s", sz_path, buffer);

		/* try to open nc file */
		ret = nc_open(sz_nc_filename, NC_NOWRITE, &id_file);
		if ( ret != NC_NOERR ) goto quit;

		ret = nc_inq(id_file, &dims_count, &vars_count, &atts_count, &unl_count);
		if ( ret != NC_NOERR ) goto quit;

		if ( ! dims_count || ! vars_count ) {
			logger(g_log, "bad nc file! %d dimensions and %d vars\n\n", dims_count, vars_count);
			goto quit_no_nc_err;
		}

		/* reset */
		for ( i = 0; i < DIMS_COUNT; ++i ) {
			dims_id[i] = -1;
			dims_size[i] = -1;
		}

		/* get dimensions */
		for ( i = 0; i < dims_count; ++i ) {
			ret = nc_inq_dim(id_file, i, name, &size);
			if ( ret != NC_NOERR ) goto quit;
			for ( y = 0; y < DIMS_COUNT; ++y ) {
				if ( ! string_compare_i(sz_dims[y], name) ) {
					if ( dims_size[y] != -1 ) {
						logger(g_log, "dimension %s already found!\n", sz_dims[y]);
						goto quit_no_nc_err;
					}
					dims_size[y] = size;
					dims_id[y] = i;
					break;
				}
			}
		}

		/* check if we have all dimensions */
		for ( i = 0; i < DIMS_COUNT; ++i ) {
			/* height_2m can be missing */
			if ( (-1 == dims_size[i]) && (i != HEIGHT_DIM) ) {
				logger(g_log, "dimension %s not found!\n", sz_dims[i]);
				goto quit_no_nc_err;
			}
		}

		/* check if x_cell is >= x_dim */
		if ( x_cell >= dims_size[X_DIM] ) {
			printf("x_cell >= x_dim: %d,%d\n", x_cell, dims_size[X_DIM]);
			goto quit_no_nc_err;
		}

		/*
		if ( ! x_cells_count ) {
			x_cells_count = dims_size[X_DIM];
			y_cells_count = dims_size[Y_DIM];
		}
		*/

		/* check if y_cell is >= y_dim */
		if ( y_cell >= dims_size[Y_DIM] ) {
			logger(g_log, "y_cell >= y_dim: %d,%d\n", y_cell, dims_size[Y_DIM]);
			goto quit_no_nc_err;
		}

		/* if we have height, it cannot be > 1 */
		if ( dims_size[HEIGHT_DIM] > 1 ) {
			logger(g_log, "height_2m cannot be > 1 : (%d)\n", dims_size[HEIGHT_DIM]);
			goto quit_no_nc_err;
		}

		/* check rows_count */
		if ( ! rows_count ) {
			rows_count = dims_size[TIME_DIM];
			/* alloc memory for double values */
			values = malloc(rows_count*MET_COLUMNS_COUNT*sizeof*values);
			if ( ! values ) {
				puts(sz_err_out_of_memory);
				nc_close(id_file);
				fclose(f);
				return 0;
			}
			/* set all double values to -9999 */
			for ( i = 0; i < rows_count*MET_COLUMNS_COUNT; ++i ) {
				values[i] = INVALID_VALUE;
			}
			/* alloc memory for float values */
			f_values = malloc(rows_count*sizeof*f_values);
			if ( ! f_values ) {
				puts(sz_err_out_of_memory);
				nc_close(id_file);
				fclose(f);
				free(values);
				return 0;
			}
			/* set all float values to -9999..not really needed */
			for ( i = 0; i < rows_count; ++i ) {
				f_values[i] = INVALID_VALUE;
			}
		} else {
			if ( rows_count != dims_size[TIME_DIM] ) {
				logger(g_log, "rows count inside %s should be %d not %d\n\n", buffer, rows_count, dims_size[TIME_DIM]);
				goto quit_no_nc_err;
			}
		}

		/* get var */
		flag = 0;
		for ( i = 0; i < vars_count; ++i ) {
			size_t start[DIMS_COUNT];
			size_t count[DIMS_COUNT];

			ret = nc_inq_var(id_file, i, name, &type, &n_dims, ids, NULL);
			if ( ret != NC_NOERR ) goto quit;
			if ( ! string_compare_i(name, sz_lat) ) {
				/* n_dims can be only 2 and ids only x and y */
				if ( 2 != n_dims ) {
					logger(g_log, "bad %s dimension size. It should be 2 not %d\n", sz_lat, n_dims);
					goto quit_no_nc_err;
				}
				/* who cames first ? x or y ? */
				if ( dims_id[X_DIM] == ids[0] ) {
					/* x first */
					start[0] = x_cell;
					start[1] = y_cell;
				} else {
					/* y first */
					start[0] = y_cell;
					start[1] = x_cell;
				}
				count[0] = count[1] = 1;
				ret = nc_get_vara_float(id_file, i, start, count, &lat);
				if ( ret != NC_NOERR ) goto quit;
			} else if ( ! string_compare_i(name, sz_lon) ) {
				/* n_dims can be only 2 and ids only x and y */
				if ( 2 != n_dims ) {
					logger(g_log, "bad %s dimension size. It should be 2 not %d\n", sz_lon, n_dims);
					goto quit_no_nc_err;
				}
				/* who cames first ? x or y ? */
				if ( dims_id[X_DIM] == ids[0] ) {
					/* x first */
					start[0] = x_cell;
					start[1] = y_cell;
				} else {
					/* y first */
					start[0] = y_cell;
					start[1] = x_cell;
				}
				count[0] = count[1] = 1;
				ret = nc_get_vara_float(id_file, i, start, count, &lon);
				if ( ret != NC_NOERR ) goto quit;
			} else if ( ! string_compare_i(name, sz_time) ) {
				if ( ! date_imported ) {
					if ( type != NC_DOUBLE ) {
						logger(g_log, "type format in %s for time column not supported\n\n", buffer);
						goto quit_no_nc_err;
					}
					ret = nc_get_var_double(id_file, i, &values[COLUMN_AT(YEAR)]);
					if ( ret != NC_NOERR ) goto quit;

					/* adjust time to YYYY,MM,DD */
					for ( z = 0; z < rows_count; ++z ) {
						int YYYY;
						int MM;
						int DD;
						timestamp_split(values[VALUE_AT(z, YEAR)], &YYYY, &MM, &DD);
						values[VALUE_AT(z, YEAR)] = YYYY;
						values[VALUE_AT(z, MONTH)] = MM;
						values[VALUE_AT(z, DAY)] = DD;
					}
					date_imported = 1;
				} else {
					/* ALESSIOR...add a datetime compare for each variable ? */
				}
			} else {
				for ( y = 0; y  < VARS_COUNT; ++y ) {
					if ( ! string_compare_i(name, sz_vars[y]) ) {
						/* check if we already have imported that var */
						if ( vars[y] ) {
							logger(g_log, "var %s already imported\n", sz_vars[y]);
							goto quit_no_nc_err;
						}
						/* ALESSIOR...do a clean up here...too much code! */
						if ( dims_id[X_DIM] == ids[0] ) {
							start[0] = x_cell;
							count[0] = 1;
						} else if ( dims_id[Y_DIM] == ids[0] ) {
							start[0] = y_cell;
							count[0] = 1;
						} else if ( dims_id[TIME_DIM] == ids[0] ) {
							start[0] = 0;
							count[0] = dims_size[TIME_DIM];
						} else {
							start[ids[0]] = 0;
							count[ids[0]] = dims_size[HEIGHT_DIM];
						}
						if ( dims_id[X_DIM] == ids[1] ) {
							start[1] = x_cell;
							count[1] = 1;
						} else if ( dims_id[Y_DIM] == ids[1] ) {
							start[1] = y_cell;
							count[1] = 1;
						} else if ( dims_id[TIME_DIM] == ids[1] ) {
							start[1] = 0;
							count[1] = dims_size[TIME_DIM];
						} else {
							start[1] = 0;
							count[1] = dims_size[HEIGHT_DIM];
						}
						if ( dims_id[X_DIM] == ids[2] ) {
							start[2] = x_cell;
							count[2] = 1;
						} else if ( dims_id[Y_DIM] == ids[2] ) {
							start[2] = y_cell;
							count[2] = 1;
						} else if ( dims_id[TIME_DIM] == ids[2] ) {
							start[2] = 0;
							count[2] = dims_size[TIME_DIM];
						} else {
							start[2] = 0;
							count[2] = dims_size[HEIGHT_DIM];
						}
						if ( dims_id[HEIGHT_DIM] != -1 ) {
							if ( dims_id[X_DIM] == ids[3] ) {
								start[3] = x_cell;
								count[3] = 1;
							} else if ( dims_id[Y_DIM] == ids[3] ) {
								start[3] = y_cell;
								count[3] = 1;
							} else if ( dims_id[TIME_DIM] == ids[3] ) {
								start[3] = 0;
								count[3] = dims_size[TIME_DIM];
							} else {
								start[3] = 0;
								count[3] = dims_size[HEIGHT_DIM];
							}
						}
						/* get values */
						if ( NC_FLOAT == type ) {
							ret = nc_get_vara_float(id_file, i, start, count, f_values);
							if ( ret != NC_NOERR ) goto quit;
							for ( z = 0; z < rows_count; ++z ) {
								values[VALUE_AT(z, y + 3)] = f_values[z];
							}
						} else if ( NC_DOUBLE == type ) {
							ret = nc_get_vara_double(id_file, i, start, count, &values[COLUMN_AT(y + 3)]);
							if ( ret != NC_NOERR ) goto quit;
						} else {
							/* type format not supported! */
							logger(g_log, "type format in %s for %s column not supported\n\n", buffer, sz_vars[y]);
							goto quit_no_nc_err;
						}
						vars[y] = 1;
						flag = 1;
						break;
					}
				}
			}
		}
		if ( ! flag ) {
			logger(g_log, "var not found inside %s\n\n", buffer);
			goto quit_no_nc_err;
		}
		nc_close(id_file);
	}
	fclose(f);
	free(f_values);
	f_values = NULL;

#if 0
	/* check for TA, TMIN && TMAX */
	if ( ! vars[TA_F-3] ) {
		if ( ! vars[VPD_F-3] && ! vars[VPD_F-3] ) {
			logger(g_log, "VPD and RH columns are missing!\n\n");
			free(values);
			return 0;
		}
	}

	/* check for missing vars */
	for ( i = 0; i < VARS_COUNT; ++i ) {
		switch ( i ) {
			case VPD_F-3:
			case RH_F-3:
			if ( ! vars[VPD_F-3] && ! vars[RH_F-3] ) {
				logger(g_log, "VPD and RH columns are missing!\n\n");
				free(values);
				return 0;
			}
			break;

			case TA_F-3:
			if ( ! vars[i] && ! vars[TMIN-3] && ! vars[TMAX-3] ) {
				logger(g_log, "TA, TMIN and TMAX columns are missing!\n\n");
				free(values);
				return 0;
			}
			break;

			case TMIN-3:
			case TMAX-3:
			if ( ! vars[i] && ! vars[TA_F-3]) {
				logger(g_log, "%s is missing!\n\n", sz_vars[i]);
				free(values);
				return 0;
			}
			break;

			default:
				if ( ! vars[i] ) {
					logger(g_log, "met columns %s is missing!\n\n", sz_vars[i]);
					free(values);
					return 0;
				}
		}
	}

	/* check if RH is valid ( not all -9999 ) */
	if ( vars[RH_F-3] ) {
		int current_row = 0; /* used as rows count for valid RH */
		for ( i = 0; i < rows_count; ++i ) {
			if ( ! IS_INVALID_VALUE(values[VALUE_AT(i, RH_F)]) ) {
				++current_row;
			}
		}
		if ( ! current_row ) {
			vars[RH_F-3] = 0;
		}
	}

	/* check if VPD is valid ( not all -9999 ) */
	if ( vars[VPD_F] ) {
		int current_row = 0; /* used as rows count for valid VPD */
		for ( i = 0; i < rows_count; ++i ) {
			if ( ! IS_INVALID_VALUE(values[VALUE_AT(i, VPD_F)]) ) {
				++current_row;
			}
		}
		if ( ! current_row ) {
			vars[VPD_F-3] = 0;
		}
	}

	if ( ! vars[RH_F-3] && ! vars[VPD_F-3] ) {
		logger(g_log, "rh and vpd not found!");
		free(values);
		return 0;
	}

	/* compute ta ? */
	if ( vars[TA_F-3] ) {
		int missings_count = 0;
		for ( i = 0; i < rows_count; ++i ) {
			if ( IS_INVALID_VALUE(values[VALUE_AT(i, TA_F)]) ) {
				++missings_count;
			}
		}
		if ( missings_count == rows_count ) {
			vars[TA_F-3] = 0;
		}
	}

	/* compute ta ! */
	if ( ! vars[TA_F-3] ) {
		for ( i = 0; i < rows_count; ++i ) {
			if ( ! IS_INVALID_VALUE(values[VALUE_AT(i, TMAX)])
					&& ! IS_INVALID_VALUE(values[VALUE_AT(i, TMIN)]) ) {
				values[VALUE_AT(i, TA_F)] = (0.606 * values[VALUE_AT(i, TMAX)]) + (0.394 * values[VALUE_AT(i, TMIN)]);
			}
		}
	}

	/* compute vpd ? or rh ?*/
	/* please note that we must have TA, so computing is done after computing TA (if needed) */
	if ( ! vars[VPD_F-3] ) {
		compute_vpd(values, rows_count, MET_COLUMNS_COUNT);
	} else if ( ! vars[RH_F-3] ) {
		compute_rh(values, rows_count, MET_COLUMNS_COUNT);
	}
#endif

#ifdef _WIN32
#ifdef _DEBUG
	{
		FILE *f;
		int row;
		char buffer[64];
		sprintf(buffer, "debug_file_%g_%g_%d_%d.txt", lat, lon, x_cell, y_cell);
		f = fopen(buffer, "w");
		if ( ! f ) {
			puts("unable to create output file!");
			free(values);
			return 0;
		}
		/* write header */
		fputs("LAT,LON,DATE,ET,LAI,RADS,SWC,TMAX,TMIN,TOT_PREC,TSOIL,VPD,WS_f\n", f);
		for ( row = 0; row < dims_size[TIME_DIM]; ++row ) {
			fprintf(f, "%g,%g,%02d/%02d/%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n"
					, lat
					, lon
					,(int)values[VALUE_AT(row,DAY)]
					, (int)values[VALUE_AT(row,MONTH)]
					, (int)values[VALUE_AT(row,YEAR)]
					, values[VALUE_AT(row,ET)]
					, values[VALUE_AT(row,NDVI_LAI)]
					, values[VALUE_AT(row,RG_F)]
					, values[VALUE_AT(row,SWC)]
					, values[VALUE_AT(row,TMAX)]
					, values[VALUE_AT(row,TMIN)]
					, values[VALUE_AT(row,PRECIP)]
					, values[VALUE_AT(row,TS_F)]
					, values[VALUE_AT(row,VPD_F)]
					, values[VALUE_AT(row,WS_F)]

			);
		}
		fclose(f);
	}
#endif
#endif

	i = yos_from_arr(values, rows_count, MET_COLUMNS_COUNT, p_yos, yos_count);
	free(values);
	return i;

quit:
	logger(g_log, nc_strerror(ret));
quit_no_nc_err:
	nc_close(id_file);
	if ( f_values ) free(f_values);
	if ( values ) free(values);
	return 0;

#undef VALUE_AT
#undef COLUMN_AT
#undef VARS_COUNT
}

static int import_txt(const char *const filename, yos_t** p_yos, int *const yos_count) {
#define BUFFER_SIZE	1024
#define VALUE_AT(r,c)	((r)+((c)*rows_count))
	int i = 0,
			column = 0,
			year = 0,
			error_flag = 0,
			columns[MET_COLUMNS_COUNT];

	int no_year_column;
	const char met_delimiter[] = " ,\t\r\n";

	int current_row;

	char *token2 = NULL,
			*p,
			buffer[BUFFER_SIZE] = { 0 };
	int rows_count;
	double *values;

	FILE *f;

	assert(p_yos);

	// get rows count
	rows_count = file_get_rows_count(filename);
	if ( rows_count <= 0 ) {
		logger(g_log, "unable to open met data file !\n");
		return 0;
	}

	// remove header
	--rows_count;

	// alloc memory for values
	values = malloc(rows_count*MET_COLUMNS_COUNT*sizeof*values);
	if ( ! values ) {
		logger(g_log, sz_err_out_of_memory);
		return 0;
	}

	// clear values with INVALID_VALUES
	for ( i = 0; i < rows_count*MET_COLUMNS_COUNT; ++i ) {
		values[i] = INVALID_VALUE;
	}

	// open file
	f = fopen(filename, "r");
	if ( !f )
	{
		logger(g_log, "unable to open met data file !\n");
		free(values);
		return 0;
	}

	// get header
	if ( ! fgets(buffer, BUFFER_SIZE, f) )
	{
		logger(g_log, "empty met data file ?\n");
		free(values);
		fclose(f);
		return 0;
	}

	// reset
	year = 0;
	for ( i = 0; i < MET_COLUMNS_COUNT; ++i )
	{
		columns[i] = -1;
	}

	// parse header
	for ( column = 0, token2 = string_tokenizer(buffer, met_delimiter, &p); token2; token2 = string_tokenizer(NULL, met_delimiter, &p), column++ )
	{
		for ( i = 0; i < MET_COLUMNS_COUNT; ++i )
		{
			if ( ! string_compare_i(token2, sz_met_columns[i]) )
			{
				if  ( -1 != columns[i] )
				{
					printf("met column %s already assigned.\n\n", token2);
					logger(g_log, "met column %s already assigned.\n\n", token2);
					free(values);
					fclose(f);
					return 0;
				}
				columns[i] = column;
				break;
			}
		}
	}

	// check if each columns was assigned
	no_year_column = 0;
	for ( i = 0; i < MET_COLUMNS_COUNT; ++i )
	{
		if ( -1 == columns[i] )
		{
			/* get year from filename */
			if ( YEAR == i )
			{
				char *p = strrchr(filename, '_');
				if ( p )
				{
					int y;

					++p;
					for ( y = 0; p[y]; ++y );
					if ( y > 4 )
					{
						p[4] = '\0';
						year = convert_string_to_int(p, &error_flag);
						if ( ! error_flag )
						{
							no_year_column = 1;
							continue;
						}
					}
				}
			}
			else if ( ((VPD_F == i) && (-1 != columns[RH_F])) || ((RH_F == i) && (-1 != columns[VPD_F])) ) {
				continue;
			}
			printf("met column %s not found.\n\n", sz_met_columns[i]);
			free(values);
			fclose(f);
			return 0;
		}
	}

	if ( (! no_year_column && (YEAR != columns[YEAR])) || (MONTH-no_year_column != columns[MONTH]) || (DAY-no_year_column != columns[DAY]) )
	{
		puts("date must be on first columns!\n\n");
		free(values);
		fclose(f);
		return 0;
	}

	// parse rows
	current_row = 0;
	while ( fgets(buffer, BUFFER_SIZE, f) )
	{
		// remove carriage return and newline ( if any )
		for ( i = 0; buffer[i]; i++ )
		{
			if ( ('\n' == buffer[i]) || ('\r' == buffer[i]) )
			{
				buffer[i] = '\0';
				break;
			}
		}

		// skip empty lines
		if ( '\0' == buffer[0] )
		{
			continue;
		}
		if ( ++current_row > rows_count ) {
			puts("too many rows found!");
			logger(g_log, "too many rows found!");
			free(values);
			fclose(f);
			return 0;
		}

		// parse string
		for ( column = no_year_column, token2 = string_tokenizer(buffer, met_delimiter, &p); token2; token2 = string_tokenizer(NULL, met_delimiter, &p), ++column )
		{
			if ( column >= MET_COLUMNS_COUNT ) {
				printf("too many columns at row %d\n", current_row+1);
				free(values);
				fclose(f);
				return 0;
			}

			for ( i = 0; i < MET_COLUMNS_COUNT; ++i ) {
				if ( column == columns[i] ) break;
			}

			values[VALUE_AT(current_row-1, i)] = convert_string_to_float(token2, &error_flag);

			if ( error_flag )
			{
				printf("unable to convert value \"%s\" for %s column\n", token2, sz_met_columns[i+no_year_column]);
				logger(g_log, "unable to convert value \"%s\" for %s column\n", token2, sz_met_columns[i+no_year_column]);
				free(values);
				fclose(f);
				return 0;
			}
		}
	}

	fclose(f);

	if ( rows_count != current_row ) {
		printf("rows count should be %d not %d\n", rows_count, current_row);
		logger(g_log, "rows count should be %d not %d\n", rows_count, current_row);
		free(values);
		return 0;
	}

	/* fix year */
	if ( no_year_column ) {
		for ( i = 0; i < rows_count; ++i ) {
			values[VALUE_AT(i, YEAR)] = year;
		}
	}

#ifdef _WIN32
#ifdef _DEBUG
	{
		FILE *f;
		int i;
		int row;
		f = fopen("debug_import_file_txt.csv", "w");
		if ( ! f ) {
			puts("unable to create output file!");
			free(values);
			return 0;
		}
		/* write header */
		fputs("Year,Month,n_days,", f);
		for ( i = 0; i < MET_COLUMNS_COUNT; ++i ) {
			fprintf(f, "%s", sz_met_columns[i]);
			if ( i < MET_COLUMNS_COUNT-1 ) {
				fputs(",", f);
			}
		}
		fputs("\n", f);
			
		for ( row = 0; row < rows_count; ++row ) {
			for ( i = 0; i < MET_COLUMNS_COUNT; ++i ) {
				fprintf(f, "%g", values[VALUE_AT(row,i)]);
				if ( i < MET_COLUMNS_COUNT-1 ) {
					fputs(",", f);
				}
			}
			fputs("\n", f);
		}
		fclose(f);
	}
#endif
#endif

#if 0
	/* check if RH is valid ( not all -9999 ) */
	if ( columns[RH_F] != -1 ) {
		current_row = 0; /* used as rows count for valid RH */
		for ( i = 0; i < rows_count; ++i ) {
			if ( ! IS_INVALID_VALUE(values[VALUE_AT(i, RH_F)]) ) {
				++current_row;
			}
		}
		if ( ! current_row ) {
			columns[RH_F] = -1;
		}
	}

	/* check if VPD is valid ( not all -9999 ) */
	if ( columns[VPD_F] != -1 ) {
		current_row = 0; /* used as rows count for valid VPD */
		for ( i = 0; i < rows_count; ++i ) {
			if ( ! IS_INVALID_VALUE(values[VALUE_AT(i, VPD_F)]) ) {
				++current_row;
			}
		}
		if ( ! current_row ) {
			columns[VPD_F] = -1;
		}
	}

	if ( (-1 == columns[RH_F]) 
			&& (-1 == columns[VPD_F]) ) {
		logger(g_log, "rh and vpd not found!");
		free(values);
		return 0;
	}

	/* compute vpd ?*/
	if ( -1 == columns[VPD_F] ) {
		compute_vpd(values, rows_count, MET_COLUMNS_COUNT);
	} else if ( -1 == columns[RH_F] ) {
		compute_rh(values, rows_count, MET_COLUMNS_COUNT);
	}
#endif
	if ( ! yos_from_arr(values, rows_count, MET_COLUMNS_COUNT, p_yos, yos_count) ) {
		free(values);
		return 0;
	}

	free(values);

	/* ok */
	return 1;
#undef BUFFER_SIZE
#undef VALUE_AT
}

/* file is the comma separated files list!!! not a single file, initially yos_count is equal to 0 */
yos_t* yos_import(const char *const file, int *const yos_count, const int x, const int y) {
	int i;
	char *token;
	char *p;
	char *p2;
	char *filename;
	char *temp;
	yos_t *yos;
	const char comma_delimiter[] = ",\r\n";

	assert(file && yos_count);

	temp = string_copy(file);
	if ( ! temp ) {
		logger(g_log, sz_err_out_of_memory);
		return NULL;
	}

	yos = NULL;
	*yos_count = 0;
	for ( token = string_tokenizer(temp, comma_delimiter, &p); token; token = string_tokenizer(NULL, comma_delimiter, &p) ) {
		i = strlen(token);
		if ( ! i ) {
			continue;
		}

		// Get only filename.txt from X:\folder\filename.txt
		// or /this/is/the/complete/path/of/filename.txt
		filename = (strrchr(token, '\\'));
		if ( filename ) ++filename;
		p2 = (strrchr(token, '/'));
		if ( p2 ) ++p2;
		if  ( p2 > filename ) {
			filename = p2;
		}
		if ( ! filename ) {
			filename = token;
		}
		i = 0; /* flag for netcdf file */
		p2 = strrchr(filename, '.');
		if ( p2 ) {
			++p2;
			if ( ! string_compare_i(p2, "nc") || ! string_compare_i(p2, "nc4") ) {
				i = 1; /* is a netcdf file */
			}
		}
		if ( i ) {
			i = import_nc(token, &yos, yos_count);
		} else if ( ! string_compare_i(p2, "lst") ) {
			i = import_lst(token, &yos, yos_count, x, y);
		} else {
			i = import_txt(token, &yos, yos_count);
		}
		if ( ! i ) {
			free(temp);
			free(yos);
			return NULL;
		}
	}
	free(temp);

	/* import co2 conc */
	if ( ! string_compare_i(g_settings->CO2_mod, "on") ) {
		if ( ! string_compare_i(g_settings->CO2_fixed, "off") ) {
			int err;

			if ( ! g_sz_co2_conc_file ) {
				logger(g_log, "co2 concentration file not specified!");
				free(yos);
				return NULL;
			}

			for ( i = 0; i < *yos_count; ++i ) {
				yos[i].co2_conc = get_co2_conc(yos[i].year, &err);
				if ( err ) {
					logger(g_log, "unable to get co2 concentration for year %d\n", yos[i].year);
					free(yos);
					return NULL;
				}
			}
		} else {
			for ( i = 0; i < *yos_count; ++i ) {
				yos[i].co2_conc = g_settings->co2Conc;
			}
		}
	}

	if ( *yos_count > 1 ) {
		qsort(yos, *yos_count, sizeof*yos, sort_by_years);
	}
#if 0
	/* save imported yos for debugging purposes */
	{
		char buffer[64];
		int i;
		int month;
		int z;
		FILE *f;

		for ( i = 0; i < *yos_count; ++i ) {
			sprintf(buffer, "debug_met_daily_%d.csv", yos[i].year);
			f = fopen(buffer, "w");
			if ( ! f ) {
				printf("unable to create %s\n", buffer);
				free(yos);
				return NULL;
			}
			fputs("year,month,day,co2_conc,n_days,solar_rad,tavg,tmax,tmin,tday,tnight,vpd,ts_f"
					",rain,swc,ndvi_lai,daylength,thermic_sum"
					",rho_air,tsoil,et,windspeed\n", f);

			for ( month = 0; month < 12; ++month ) {
				for ( z = 0; z < 31; ++z ) {
					if ( z == yos[i].m[month].d[z].n_days )
					{
						break;
					}
					fprintf(f, "%d,%d,%d,%g,%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n"
							, yos[i].year
							, month+1
							, z+1
							, yos[i].co2_conc
							, yos[i].m[month].d[z].n_days
							, yos[i].m[month].d[z].solar_rad
							, yos[i].m[month].d[z].tavg
							, yos[i].m[month].d[z].tmax
							, yos[i].m[month].d[z].tmin
							, yos[i].m[month].d[z].tday
							, yos[i].m[month].d[z].tnight
							, yos[i].m[month].d[z].vpd
							, yos[i].m[month].d[z].ts_f
							, yos[i].m[month].d[z].rain
							, yos[i].m[month].d[z].swc
							, yos[i].m[month].d[z].ndvi_lai
							, yos[i].m[month].d[z].daylength
							, yos[i].m[month].d[z].thermic_sum
							, yos[i].m[month].d[z].rho_air
							, yos[i].m[month].d[z].tsoil
							, yos[i].m[month].d[z].et
							//ALESSIOC
							, yos[i].m[month].d[z].windspeed);
				}
			}
			fclose(f);
		}
	}
#endif

	return yos;
}
