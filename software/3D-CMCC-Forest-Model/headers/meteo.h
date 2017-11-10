/* meteo.h */
#ifndef METEO_H_
#define METEO_H_

#define METEO_MONTHS_COUNT		12
#define METEO_DAYS_COUNT		31
#define METEO_HOURS_COUNT		24
#define METEO_HALFHOURS_COUNT	2

#define METEO_COMMON_MEMBERS																		\
		int n_days;																					\
		double solar_rad;						   /* (MJ/m2/day) incoming short wave radiation */ \
		double tavg;					           /* (deg C) daily average air temperature */ \
		double tdew;					           /* (deg C) daily dew average daily air temperature */ \
		double vpd;						           /* (hPa-mbar) daily vapour pressure dificit */ \
		double rh_f;					           /* (%) daily relative humidity */ \
		double ts_f;					           /* (deg C) daily soil temperature from data file */ \
		double prcp;					           /* (mm/m2) daily precipitation */ \
		double rain;					           /* (mm/m2) daily rain */ \
		double snow;					           /* (cm/m2) daily snow */ \
		double swc;						           /* (% vol) daily soil wtaer content from data file */ \
		double rho_air;					           /* (kg/m3) daily air density */ \
		double tsoil;					           /* (deg C) daily soil temperature computed by the model */ \
		double et;						           /* (mm) daily evapo-transpiration (not used)*/ \
		double windspeed;				           /* (m/sec) daily wind speed */ \
		double lh_vap;					           /* (J/Kg) daily latent heat of vaporization trees */ \
		double lh_vap_soil;				           /* (J/Kg) daily latent heat of vaporization soil */ \
		double lh_fus;					           /* (KJ/Kg) daily latent heat of fusion snow */ \
		double lh_sub;					           /* (KJ/Kg) daily latent heat of sublimation snow */ \
		double air_pressure;			           /* (Pa) daily air pressure */ \
		double es;						           /* (KPa) daily weighted mean saturation vapour pressure at the air temperature */ \
		double ea;						           /* (KPa) daily actual vapour pressure derived from relative humidity data */ \
		double psych;						       /* (KPa/°C) daily psychrometric constant */ \
		double sw_pot_downward_W;		           /* (W/m2) potential shortwave downward radiation */ \
		double sw_downward_MJ;			           /* (MJ/m2/day) daily downward short wave radiation */ \
		double incoming_sw_downward_W;	           /* (W/m2) incoming daily downward short wave radiation */ \
		double sw_downward_W;                      /* (W/m2) daily downward short wave radiation */ \
		double lw_downward_MJ;                     /* (MJ/m2/day) daily downward long wave radiation */ \
		double atm_lw_downward_W;                  /* (W/m2) daily downward atmopsheric long wave radiation */ \
		double lw_net_MJ;                          /* (MJ/m2/day) daily net long wave radiation */ \
		double lw_net_W;                           /* (W/m2) daily net long wave radiation */ \
		double incoming_par;                       /* (molPAR/m2/day) incoming PAR */ \
		double par;                                /* (molPAR/m2/day) PAR */ \
		double incoming_ppfd;                      /* (umol/m2/sec) incoming PPFD */ \
		double ppfd;                               /* (umol/m2/sec) PPFD */ \
		double emis_atm_clear_sky;                 /* (ratio) atmospheric emissivity with clear sky */ \
		double emis_atm;                           /* (ratio) atmospheric emissivity corrected for cloud cover */ \
		double cloud_cover_frac;                   /* (ratio) daily cloud cover fraction */ \
		double cloud_cover_frac_corr;              /* (ratio) daily cloud cover fraction */ 

typedef struct {
	METEO_COMMON_MEMBERS
	double tmax;                                /* (deg C) daily maximum air temperature */
	double tmin;                                /* (deg C) daily minimum air temperature */
	double tday;                                /* (deg C) daylight  daily air temperature */
	double tnight;                              /* (deg C) nightime average daily air temperature */
	double daylength;                           /* (hours) day length */
	double daylength_sec;						/* (seconds) day length */
	double ni;                                  /* (ratio) fraction of day length */
	double thermic_sum;                         /* (deg C) daily thermic sum */
	double ndvi_lai;                            /* (m2/m2) daily LAI from data file */
	double ten_day_avg_tavg;                    /* (°C) ten day avg tavg */
	double ten_day_avg_tsoil;                   /* (°C) ten day avg tsoil */
	double ten_day_avg_tday;                    /* (°C) ten day avg tday */
	double ten_day_avg_tnight;                  /* (°C) ten day avg tnight */
	double ten_day_weighted_avg_tavg;           /* (°C) ten day weighted avg tavg */
	double ten_day_weighted_avg_tsoil;          /* (°C) ten day weighted avg tsoil */
	double ten_day_weighted_avg_tday;           /* (°C) ten day weighted avg tday */
	double ten_day_weighted_avg_tnight;         /* (°C) ten day weighted avg tnight */
	double Ndeposition;                         /* (gN/m2/day) nitrogen deposition scaled from annual */
	double Net_rad_threePG;
} meteo_daily_t;

typedef struct {
	meteo_daily_t d[METEO_DAYS_COUNT];
} meteo_d_t;

/* hourly and halfhourly common struct */

typedef struct {
	METEO_COMMON_MEMBERS
	double tmax;                                /* (deg C) daily maximum air temperature */
	double tmin;                                /* (deg C) daily minimum air temperature */
} meteo_h_hh_t;

/* hourly struct */

typedef struct {
	meteo_h_hh_t h[METEO_HOURS_COUNT];
} meteo_day_h_t;

typedef struct {
	meteo_day_h_t d[METEO_DAYS_COUNT];
} meteo_h_t;

/* halfhourly struct */

typedef struct {
	meteo_h_hh_t hh[METEO_HALFHOURS_COUNT];
} meteo_hour_t;

typedef struct {
	meteo_hour_t h[METEO_HOURS_COUNT];
} meteo_day_hh_t;

typedef struct {
	meteo_day_hh_t d[METEO_DAYS_COUNT];
} meteo_hh_t;

typedef struct {
	double solar_rad;
	double tavg;
	double tmax;
	double tmin;
	double tday;
	double tnight;
	double vpd;
	double prcp;
	double tsoil;
	double rh_f;
	double incoming_par;
	double par;
	double incoming_ppfd;
	double ppfd;
} meteo_mean_t;

typedef struct {
	int year;
	double co2Conc;	                             /* (ppmv) annual atmospheric CO2 concentration */
	double Ndep;	                             /* (kgN/m2/year) annual nitrogen deposition */
	
	//void* m;
	meteo_hh_t* halfhourly;
	meteo_h_t* hourly;
	meteo_d_t daily[METEO_MONTHS_COUNT];
	meteo_mean_t monthly[METEO_MONTHS_COUNT];
	meteo_mean_t yearly;
} meteo_annual_t;

meteo_annual_t* meteo_annual_import(const char *const file, int *const yos_count, const int x, const int y);
void meteo_annual_free(meteo_annual_t* p, const int count);

#endif /* METEO_H_ */
