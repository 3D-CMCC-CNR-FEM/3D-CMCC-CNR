/* settings.h */
#ifndef SETTINGS_H_
#define SETTINGS_H_

enum { CO2_TRANS_OFF, CO2_TRANS_ON, CO2_TRANS_VAR };
enum { MANAGEMENT_OFF, MANAGEMENT_ON, MANAGEMENT_VAR, MANAGEMENT_VAR1 };
enum { IRRIGATION_OFF, IRRIGATION_ON, IRRIGATION_VAR };// carlo_irrigation
enum { THINNING_REGIME_ABOVE, THINNING_REGIME_BELOW };

#define SETTINGS_SITENAME_MAX_SIZE					(255+1)
#define SETTINGS_REPLANTED_SPECIES_MAX_SIZE			(32+1)
#define SETTINGS_REGENERATION_SPECIES_MAX_SIZE		(32+1)

typedef struct
{
	char species[SETTINGS_REPLANTED_SPECIES_MAX_SIZE];
	char management;           /*this is a mandatory variables that each class has to have AND DETERMINES HABITUS OF PLANT*/
	double n_tree;
	double age;
	double avdbh;
	double lai;
	double height;

} replanted_vars_t;

typedef struct
{
	int day;
	int month;
	int year;

	double intensity;
	double height; /* height after pruning */

} pruning_t;

typedef struct
{
	int day;
	int month;
	int year;

	double amount; 

} irrigation_t;

typedef struct {
	char sitename[SETTINGS_SITENAME_MAX_SIZE];
	char version;                        /* fixme to remove */
	char spatial;	                     /* spatial or unspatial must be 's' or 'u' (not longer used) fixme to remove */
	char time;		                     /* time scale simulation must be 'm' (monthòy),'d' (daily), 'h' (hourly), 's' (semihourly) */
	char spinup;
	int spinup_years;
	char screen_output;                  /* screen output must be 'on' or 'off' */
	char debug_output;                   /* debug output must be 'on' or 'off' */
	char daily_output;                   /* daily output must be 'on' or 'off' */
	char monthly_output;                 /* monthly output must be 'on' or 'off' */
	char yearly_output;                  /* annual output must be 'on' or 'off' */
	char yearly_man_output;              /* annual management output must be 'on' or 'off' */
	char soil_output;                    /* soil output must be 'on' or 'off' */
	char netcdf_output;                  /* netcdf output must be 'on' or 'off' */
	int year_start;                      /* starting year of simulation */
	int year_end;                        /* ending year of simulation */
	int year_restart;                    /* year to restart simulation */
	int PSN_mod;                         /* photosynthesis (0 = FvCB; 1 = LUE) */
	char CO2_trans;                      /* CO2 transient must be 'on' or 'off' or 'var' */
	int year_start_co2_fixed;            /* for CO2_trans = var, year at which CO2 is fixed */
	char Ndep_fixed;                     /* fixed Nitrogen deposition must be 'on' or 'off' */
	double Tbase_resp;                   /* T base temperature for respiration (°C) */
	char Photo_accl;
	char Resp_accl;                      /* acclimation of repiration must be 'on' or 'off' */
	char Cold_accl;                      /* acclimation of GPP from winter to spring must be 'on' or 'off' */    //5p7
	char N_reg;                          /* downregulation of the GPP via soil fertility acclimation of GPP from winter to spring must be 'on' or 'off' */    //5p7
	char regeneration;                   /* regeneration must be 'on' or 'off' */
	char management;                     /* management must be 'on' of 'off' */
	char management_type;                /* optional: 0 ( default ) for bau, 1 for local */
	int year_start_management;           /* start year management (if management is on) */
	char Prog_Aut_Resp;                  /* Prognostic autotrophic respiration */
	char thinning_regime;                /* thinning regime above = A or below = B */
	int  thinning_branch_rem;                /* percentage of branches left in the stand as CWD after thinning harvest*/
	char regeneration_species[SETTINGS_REGENERATION_SPECIES_MAX_SIZE]; /* species name of regeneration species */
	
	// ALESSIOR: use e_management from matrix.h not char
	char regeneration_management;        /*this is a mandatory variables that each class has to have AND DETERMINES HABITUS OF PLANT*/

	/* DO NOT MODIFY BELOW, PLEASE INSERT STUFF IN SPACE ABOVE */

	double sizeCell;                     /* side length of the cell (meters) */

	double Fixed_Aut_Resp_rate;          /* It is the "Y"values for fixed autotrophic respiration */

	double co2Conc;                      /* reference CO2 concentration (ppmv) */
	double co2_incr;                     /* annual CO2 increment fixme to remove */

	double init_frac_maxasw;             /* minimum fraction of available soil water at the beginning of simulation */

	/* layer settings */
	double tree_layer_limit;             /* differences among tree heights before to consider different canopy layers (meters) */
	double number_of_soil_layer;         /* number of soil layers to simulate */
	double max_layer_cover;              /* maximum layer cover */

	/* replanted input (Management options)*/
	replanted_vars_t* replanted;
	int replanted_count;

	/* regeneration input (Regeneration options)*/
	double regeneration_n_tree;
	double regeneration_age;             /* age of regeneration (1 year) */
	double regeneration_avdbh;           /* dbh for regeneration (cm) */
	double regeneration_lai;             /* regeneration LAI (m2/m2) (mandatory for evergreen) */
	double regeneration_height;          /* regeneration height (meter) */

	/* pruning */
	pruning_t* prunings;
	int prunings_count;

	/* irrigation */
	irrigation_t* irrigations;
	int irrigations_count;

} settings_t;

settings_t* settings_import(const char *const filename);
void settings_free(settings_t* s);

#endif /* SETTINGS_H_ */
