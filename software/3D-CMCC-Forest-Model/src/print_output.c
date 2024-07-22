/* cumulative_balance.c */
#ifdef USE_NEW_OUTPUT

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "print_output.h"
#include "common.h"
#include "settings.h"
#include "soil_settings.h"
#include "logger.h"
#include "g-function.h"
#include "nc.h"

extern settings_t* g_settings;
extern logger_t* g_daily_log;
extern logger_t* g_monthly_log;
extern logger_t* g_annual_log;
extern logger_t* g_daily_soil_log;
extern logger_t* g_monthly_soil_log;
extern logger_t* g_annual_soil_log;
extern char *g_sz_dataset_file;
extern char *g_sz_soil_file;
extern char *g_sz_input_met_file;
extern char *g_sz_settings_file;
extern char *g_sz_topo_file;
extern char *g_sz_co2_conc_file;

extern const char sz_launched[];

extern int MonthLength [];
extern int MonthLength_Leap [];

static const char sz_management[] = "TCN";

static const char* get_filename(const char *const s)
{
	const char *p;
	const char *p2;

	p = NULL;

	if ( s ) {
		p = strrchr(s, '/');
		if ( p ) ++p;
		p2 = strrchr(s, '\\');
		if ( p2 ) ++p2;
		if ( p2 > p ) p = p2;
		if ( ! p ) p = s;
	}

	return p;
}

void print_model_paths(logger_t *const _log)
{
	//assert(_log);

	logger(_log, "\n\n#site: %s\n", g_settings->sitename);
	if ( g_sz_dataset_file )
		logger(_log, "#--input files--\n");
	logger(_log, "#input file = %s\n", get_filename(g_sz_dataset_file));
	logger(_log, "#soil file = %s\n", get_filename(g_sz_soil_file));
	logger(_log, "#topo file = %s\n", get_filename(g_sz_topo_file));
	logger(_log, "#met file = %s\n", get_filename(g_sz_input_met_file));
	logger(_log, "#settings file = %s\n", get_filename(g_sz_settings_file));
	if ( g_settings->CO2_trans )
		logger(_log, "#CO2 file = %s\n", get_filename(g_sz_co2_conc_file));
}

void print_model_settings(logger_t*const log)
{
	//assert(log);

	logger(log, "#--model settings--\n");
	logger(log, "#CO2_mod = %s\n", g_settings->CO2_mod ? "on" : "off");
	logger(log, "#CO2 trans = %s\n", (CO2_TRANS_VAR == g_settings->CO2_trans) ? "var" : (CO2_TRANS_ON == g_settings->CO2_trans) ? "on" : "off");
	if ( CO2_TRANS_OFF == g_settings->CO2_trans )
	{
		logger(log, "#fixed co2 concentration = %g ppmv\n", g_settings->co2Conc);
	}
	else if ( CO2_TRANS_VAR == g_settings->CO2_trans )
	{
		logger(log, "#year %d at which co2 concentration is fixed at value = %g ppmv\n", g_settings->year_start_co2_fixed, g_settings->co2Conc);
	}
	logger(log, "#Resp accl = %s\n", g_settings->Resp_accl ? "on" : "off");
	logger(log, "#regeneration = %s\n", g_settings->regeneration ? "on" : "off");
	logger(log, "#Management = %s\n", (MANAGEMENT_VAR == g_settings->management) ? "var" : (MANAGEMENT_ON == g_settings->management) ? "on" : "off");
	if ( g_settings->management )
	{
		logger(log, "#Year Start Management = %d\n", g_settings->year_start_management);
	}
	if ( g_settings->year_restart != -1 )
	{
		logger(log, "#Year restart = %d\n", g_settings->year_restart);
	}
	else
	{
		logger(log, "#Year restart = off\n");
	}
}

static int check_height_index(const int index, const height_t* const heights, const int heights_count)
{
	int i;

	assert(heights);

	for ( i = 0; i < heights_count; ++i )
	{
		if ( index == heights[i].index )
		{
			return 1;
		}
	}
	return 0;
}

static int check_dbh_index(const int index, const dbh_t* const dbhs, const int dbhs_count)
{
	int i;

	assert(dbhs);

	for ( i = 0; i < dbhs_count; ++i )
	{
		if ( index == dbhs[i].index )
		{
			return 1;
		}
	}
	return 0;
}

static int check_age_index(const int index, const age_t* const ages, const int ages_count)
{
	int i;

	assert(ages);

	for ( i = 0; i < ages_count; ++i )
	{
		if ( index == ages[i].index )
		{
			return 1;
		}
	}
	return 0;
}

static int check_species_index(const int index, const species_t* const species, const int species_count)
{
	int i;

	assert(species);

	for ( i = 0; i < species_count; ++i )
	{
		if ( index == species[i].index )
		{
			return 1;
		}
	}
	return 0;
}

void EOD_print_output_cell_level(cell_t *const c, const int day, const int month, const int year, const int years_of_simulation )
{
	int layer;
	int height;
	int dbh;
	int age;
	int species;

	static int years_counter;

	species_t *s;

	/* return if daily logging is off*/
	if ( ! g_daily_log ) return;

	/* heading */
	if ( ! day && ! month && ! year )
	{
		logger(g_daily_log, "X,Y,YEAR,MONTH,DAY");

		for ( layer = c->tree_layers_count - 1; layer >= 0; --layer )
		{
			/* heading for layers */
			logger(g_daily_log, ",LAYER");

			for ( height = 0; height < c->heights_count+c->heights_avail; ++height )
			{
				if ( check_height_index(height, c->heights, c->heights_count) )
				{
					if( layer == c->heights[height].height_z )
					{							
						/* heading for heights value */
						logger(g_daily_log,",HEIGHT");

						for ( dbh = 0; dbh < c->heights[height].dbhs_count+c->heights[height].dbhs_avail; ++dbh )
						{							
							logger(g_daily_log,",DBH");

							if ( check_dbh_index(dbh, c->heights[height].dbhs, c->heights[height].dbhs_count) )
							{
								for ( age = 0; age < c->heights[height].dbhs[dbh].ages_count+c->heights[height].dbhs[dbh].ages_avail; ++age )
								{								
									/* heading for ages */
									logger(g_daily_log,",AGE");

									if ( check_age_index(age, c->heights[height].dbhs[dbh].ages, c->heights[height].dbhs[dbh].ages_count) )
									{									
										for ( species = 0; species < c->heights[height].dbhs[dbh].ages[age].species_count+c->heights[height].dbhs[dbh].ages[age].species_avail; ++species )
										{										
											if ( check_species_index(species, c->heights[height].dbhs[dbh].ages[age].species, c->heights[height].dbhs[dbh].ages[age].species_count)  )
											{
												/* heading for species name */
												logger(g_daily_log,",SPECIES");

												/* heading for management */
												logger(g_daily_log, ",MANAGEMENT");

												logger(g_daily_log,
														",GPP"
														",RG"
														",RM"
														",RA"
														",NPP"
														",CUE"
														",LAI_PROJ"
														",PEAK-LAI_PROJ"
														",LAI_EXP"
														",D-CC_P"
														",DBHDC"
														",CROWN_AREA_PROJ"
														",CROWN_AREA_EXP"
														",PAR"
														",APAR"
														",fAPAR"
														",Ntree"
														",VEG_D"
														",C_INT"
														",C_WAT"
														",C_EVA"
														",C_TRA"
														",C_ET"
														",C_LE"
														",WUE"
														",WRes"
														",WS"
														",WSsap"
														",WSL"
														",WSD"
														",WL"
														",WFR"
														",WCR"
														",WCRsap"
														",WCRL"
														",WCRD"
														",WBB"
														",WBBsap"
														",WBBL"
														",WBBD"
														",WFru"
														",dWRes"
														",dWS"
														",dWL"
														",dWFR"
														",dWCR"
														",dWBB"
														",dFRUIT"
														",SAR"
														",LAR"
														",FRAR"
														",CRAR"
														",BBAR"
														",FCO2"
														",FCO2_TR"
														",FLIGHT"
														",FAGE"
														",FT"
														",FVPD"
														",FN"
														",FSW"
														",LITR_C"
														",CWD_C"
														",SOIL_C"
												);
											}
										}
									}
								}
							}
						}
					}
				}
			}
			/*
							if ( c->heights[height].dbhs[dbh].ages[age].species_count > 1 ) logger(g_daily_log,",*");
						}
						if ( c->heights[height].dbhs[dbh].ages_count > 1 ) logger(g_daily_log,",**");
					}
					if ( c->heights[height].dbhs_count > 1 ) logger(g_daily_log,",***");
				}
				if ( c->tree_layers[layer].layer_n_height_class > 1 ) logger(g_daily_log,",****");
			}
			if ( c->tree_layers_count > 1 ) logger(g_daily_log,",*****");
			 */
		}
		/************************************************************************/
		/* heading variables only at cell level */
		logger(g_daily_log,",gpp,npp,ar,et,le,soil_evapo,snow_pack,asw,iWue,litrC,cwdC,soilC,litrN,soilN,Tsoil,Daylength\n");
	}
	/*****************************************************************************************************/

	/* values */
	logger(g_daily_log, "%d,%d,%d,%d,%d", c->x, c->y, c->years[year].year, month + 1, day + 1);

	// ALESSIOR commented on 10/07/2017
	//qsort(c->heights, c->heights_count, sizeof(height_t), sort_by_heights_index_asc);

	/* print class level values */
	for ( layer = c->tree_layers_count - 1; layer >= 0; --layer )
	{
		int i;

		/* print layer */
		logger(g_daily_log,",%d", layer);

		for ( height = 0; height < c->heights_count+c->heights_avail; ++height )
		{
			/* print height */
			if ( check_height_index(height, c->heights, c->heights_count) )
			{
				if( layer == c->heights[height].height_z )
				{
					logger(g_daily_log,",%g", c->heights[height].value);

					// ALESSIOR commented on 10/07/2017
					//qsort(c->heights[height].dbhs, c->heights[height].dbhs_count, sizeof(dbh_t), sort_by_dbhs_index_asc);

					for ( dbh = 0; dbh < c->heights[height].dbhs_count+c->heights[height].dbhs_avail; ++dbh )
					{
						if ( check_dbh_index(dbh, c->heights[height].dbhs, c->heights[height].dbhs_count) )
						{
							/* print dbh */
							logger(g_daily_log, ",%g", c->heights[height].dbhs[dbh].value);

							// ALESSIOR commented on 10/07/2017
							//qsort(c->heights[height].dbhs[dbh].ages, c->heights[height].dbhs[dbh].ages_count, sizeof(age_t), sort_by_ages_index_asc);

							for ( age = 0; age < c->heights[height].dbhs[dbh].ages_count+c->heights[height].dbhs[dbh].ages_avail; ++age )
							{
								if ( check_age_index(age, c->heights[height].dbhs[dbh].ages, c->heights[height].dbhs[dbh].ages_count) )
								{
									/* print age */
									logger(g_daily_log,",%d", c->heights[height].dbhs[dbh].ages[age].value);

									// ALESSIOR commented on 10/07/2017
									//qsort(c->heights[height].dbhs[dbh].ages[age].species, c->heights[height].dbhs[dbh].ages[age].species_count, sizeof(species_t), sort_by_species_index_asc);

									for ( species = 0; species < c->heights[height].dbhs[dbh].ages[age].species_count+c->heights[height].dbhs[dbh].ages[age].species_avail; ++species )
									{
										if ( check_species_index(species, c->heights[height].dbhs[dbh].ages[age].species, c->heights[height].dbhs[dbh].ages[age].species_count)  )
										{
											s = &c->heights[height].dbhs[dbh].ages[age].species[species];

											/* print species name */
											logger(g_daily_log,",%s", c->heights[height].dbhs[dbh].ages[age].species[species].name);

											/* print management */
											logger(g_daily_log,",%c", sz_management[c->heights[height].dbhs[dbh].ages[age].species[species].management]);

											/* print variables at layer-class level */
											logger(g_daily_log,",%6.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,"
													"%d,%d,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f"
													",%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f"
													",%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f"
													",%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f",
													s->value[GPP],
													s->value[TOTAL_GROWTH_RESP],
													s->value[TOTAL_MAINT_RESP],
													s->value[TOTAL_AUT_RESP],
													s->value[NPP],
													s->value[CUE],
													s->value[LAI_PROJ],
													s->value[PEAK_LAI_PROJ],
													s->value[LAI_EXP],
													s->value[DAILY_CANOPY_COVER_PROJ],
													s->value[DBHDC_EFF],
													s->value[CROWN_AREA_PROJ],
													s->value[CROWN_AREA_EXP],
													s->value[PAR],
													s->value[APAR],
													s->value[fAPAR],
													s->counter[N_TREE],
													s->counter[VEG_DAYS],
													s->value[CANOPY_INT_RAIN],
													s->value[CANOPY_WATER],
													s->value[CANOPY_EVAPO],
													s->value[CANOPY_TRANSP],
													s->value[CANOPY_EVAPO_TRANSP],
													s->value[CANOPY_LATENT_HEAT],
													s->value[WUE],
													s->value[RESERVE_C],
													s->value[STEM_C],
													s->value[STEM_SAPWOOD_C],
													s->value[STEM_LIVEWOOD_C],
													s->value[STEM_DEADWOOD_C],
													s->value[LEAF_C],
													s->value[FROOT_C],
													s->value[CROOT_C],
													s->value[CROOT_SAPWOOD_C],
													s->value[CROOT_LIVEWOOD_C],
													s->value[CROOT_DEADWOOD_C],
													s->value[BRANCH_C],
													s->value[BRANCH_SAPWOOD_C],
													s->value[BRANCH_LIVEWOOD_C],
													s->value[BRANCH_DEADWOOD_C],
													s->value[FRUIT_C],
													s->value[C_TO_RESERVE],
													s->value[C_TO_STEM],
													s->value[C_TO_LEAF],
													s->value[C_TO_FROOT],
													s->value[C_TO_CROOT],
													s->value[C_TO_BRANCH],
													s->value[C_TO_FRUIT],
													s->value[STEM_AUT_RESP],
													s->value[LEAF_AUT_RESP],
													s->value[FROOT_AUT_RESP],
													s->value[CROOT_AUT_RESP],
													s->value[BRANCH_AUT_RESP],
													s->value[F_CO2],
													s->value[F_CO2_TR],
													s->value[F_LIGHT_MAKELA],
													s->value[F_AGE],
													s->value[F_T],
													s->value[F_VPD],
													s->value[F_NUTR],
													s->value[F_SW],
													s->value[LITR_C],
													s->value[CWD_C],
													s->value[SOIL_C]
											);
										}
										else
										{
											for ( i = 0; i < 66; ++i )
											{
												logger(g_daily_log, ",%d", -9999);
											}
										}
									}
								}
								else
								{
									for ( i = 0; i < 67; ++i )
									{
										logger(g_daily_log, ",%d", -9999);
									}
								}
							}
						}
						else
						{
							for ( i = 0; i < 68; ++i )
							{
								logger(g_daily_log, ",%d", -9999);
							}
						}
					}
				}
			}
			else
			{
				for ( i = 0; i < 69; ++i )
				{
					logger(g_daily_log,",%d", -9999);
				}
			}
		}
	}
	/************************************************************************/
	/* printing variables only at cell level */

	logger(g_daily_log, ",%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f\n",
			c->daily_gpp,
			c->daily_npp,
			c->daily_aut_resp,
			c->daily_et,
			c->daily_lh_flux,
			c->daily_soil_evapo,
			c->snow_pack,
			c->asw,
			c->daily_iwue,
			c->litrC,
			c->cwdC,
			c->soilC,
			c->litrN,
			c->soilN,
			c->years[year].m[month].d[day].tsoil,
			c->years[year].m[month].d[day].daylength
	);
	/************************************************************************/

	if ( c->doy == ( IS_LEAP_YEAR( c->years[year].year ) ? 366 : 365 ) )
	{
		++years_counter;

		if ( years_counter ==  years_of_simulation)
		{
			logger(g_daily_log, sz_launched, netcdf_get_version(), datetime_current());
			print_model_paths(g_daily_log);
			//const char* p;
			//p = file_get_name_only(g_daily_log->filename);
			logger(g_daily_log, "#output file = %s\n", g_daily_log->filename);
			print_model_settings(g_daily_log);
		}
	}
}

void EOM_print_output_cell_level(cell_t *const c, const int month, const int year, const int years_of_simulation )
{
	int layer;
	int height;
	int dbh;
	int age;
	int species;

	static int years_counter;

	species_t *s;

	/* return if monthly logging is off*/
	if ( ! g_monthly_log ) return;

	/* heading */
	if ( !month && !year )
	{
		logger(g_monthly_log, "X,Y,YEAR,MONTH");

		for ( layer = c->tree_layers_count - 1; layer >= 0; --layer )
		{
			/* heading for layers */
			logger(g_monthly_log,",LAYER");

			for ( height = 0; height < c->heights_count+c->heights_avail; ++height )
			{
				if ( check_height_index(height, c->heights, c->heights_count)  )
				{
					if ( layer == c->heights[height].height_z )
					{
						/* heading for heights value */
						logger(g_monthly_log,",HEIGHT");

						for ( dbh = 0; dbh < c->heights[height].dbhs_count+c->heights[height].dbhs_avail; ++dbh )
						{
							/* heading for dbhs value */
							logger(g_monthly_log, ",DBH");

							if ( check_dbh_index(dbh, c->heights[height].dbhs, c->heights[height].dbhs_count)  )
							{
								for ( age = 0; age < c->heights[height].dbhs[dbh].ages_count+c->heights[height].dbhs[dbh].ages_avail; ++age )
								{
									/* heading for ages */
									logger(g_monthly_log,",AGE");

									if ( check_age_index(age, c->heights[height].dbhs[dbh].ages, c->heights[height].dbhs[dbh].ages_count) )
									{
										for ( species = 0; species < c->heights[height].dbhs[dbh].ages[age].species_count+c->heights[height].dbhs[dbh].ages[age].species_avail; ++species )
										{
											if ( check_species_index(species, c->heights[height].dbhs[dbh].ages[age].species, c->heights[height].dbhs[dbh].ages[age].species_count)  )
											{
												/* heading for species name */
												logger(g_monthly_log,",SPECIES");

												/* heading for management */
												logger(g_monthly_log,",MANAGEMENT");

												logger(g_monthly_log,
														",GPP"
														",RA"
														",NPP"
														",CUE"
														",CTRANSP"
														",CET"
														",CLE"
														",CC"
														",DBHDC"
														",HD_EFF"
														",HDMAX"
														",HDMIN"
														",N_TREE"
														",WUE"
														",WRes"
														",WS"
														",WSL"
														",WSD"
														",PWL"
														",PWFR"
														",WCR"
														",WCRL"
														",WCRD"
														",WBB"
														",WBBL"
														",WBBD"
												);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		/************************************************************************/
		/* heading variables at cell level only if there's more than one layer */

		if( c->heights_count > 1 )
		{
			logger(g_monthly_log,",***,gpp,npp,ar");
		}
		/* heading variables at cell level also if there's more than one layer */
		else
		{
			logger(g_monthly_log,",*****");
		}
		/* heading variables only at cell level */
		logger(g_monthly_log,",et,le,asw,iWue\n");
	}

	/************************************************************************/
	/* values */
	logger ( g_monthly_log, "%d,%d,%d,%d", c->x, c->y, c->years[year].year, month +1 );

	/* print class level values */
	for ( layer = c->tree_layers_count - 1; layer >= 0; --layer )
	{
		int i;

		/* print layer */
		logger(g_monthly_log,",%d", layer);

		for ( height = 0; height < c->heights_count+c->heights_avail; ++height )
		{
			if ( check_height_index(height, c->heights, c->heights_count) )
			{
				if( layer == c->heights[height].height_z )
				{
					/* print height */
					logger(g_monthly_log,",%g", c->heights[height].value);

					for ( dbh = 0; dbh < c->heights[height].dbhs_count+c->heights[height].dbhs_avail; ++dbh )
					{
						/* print dbh */
						if ( check_dbh_index(dbh, c->heights[height].dbhs, c->heights[height].dbhs_count) )
						{
							logger(g_monthly_log,",%g", c->heights[height].dbhs[dbh].value);

							for ( age = 0; age < c->heights[height].dbhs[dbh].ages_count+c->heights[height].dbhs[dbh].ages_avail; ++age )
							{
								if ( check_age_index(age, c->heights[height].dbhs[dbh].ages, c->heights[height].dbhs[dbh].ages_count) )
								{
									/* print age */
									logger(g_monthly_log,",%d", c->heights[height].dbhs[dbh].ages[age].value);

									for ( species = 0; species < c->heights[height].dbhs[dbh].ages[age].species_count+c->heights[height].dbhs[dbh].ages[age].species_avail; ++species )
									{
										if ( check_species_index(species, c->heights[height].dbhs[dbh].ages[age].species, c->heights[height].dbhs[dbh].ages[age].species_count) )
										{
											s  = &c->heights[height].dbhs[dbh].ages[age].species[species];

											/* print species name */
											logger(g_monthly_log,",%s", c->heights[height].dbhs[dbh].ages[age].species[species].name);

											/* print management */
											logger(g_monthly_log,",%c", sz_management[c->heights[height].dbhs[dbh].ages[age].species[species].management]);

											/* print variables at layer-class level */
											logger(g_monthly_log,",%6.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%d,%3.4f,%3.4f,%3.4f,%3.4f"
													",%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f",
													s->value[MONTHLY_GPP],
													s->value[MONTHLY_TOTAL_AUT_RESP],
													s->value[MONTHLY_NPP],
													s->value[MONTHLY_CUE],
													s->value[MONTHLY_CANOPY_TRANSP],
													s->value[MONTHLY_CANOPY_EVAPO_TRANSP],
													s->value[MONTHLY_CANOPY_LATENT_HEAT],
													s->value[CANOPY_COVER_PROJ],
													s->value[DBHDC_EFF],
													s->value[HD_EFF],
													s->value[HD_MAX],
													s->value[HD_MIN],
													s->counter[N_TREE],
													s->value[MONTHLY_WUE],
													s->value[RESERVE_C],
													s->value[STEM_C],
													s->value[STEM_LIVEWOOD_C],
													s->value[STEM_DEADWOOD_C],
													s->value[MAX_LEAF_C],
													s->value[MAX_FROOT_C],
													s->value[CROOT_C],
													s->value[CROOT_LIVEWOOD_C],
													s->value[CROOT_DEADWOOD_C],
													s->value[BRANCH_C],
													s->value[BRANCH_LIVEWOOD_C],
													s->value[BRANCH_DEADWOOD_C]);
										}
										else
										{
											for ( i = 0; i < 28; ++i )
											{
												logger(g_monthly_log, ",%d", -9999);
											}
										}
									}

								}
								else
								{
									for ( i = 0; i < 29; ++i )
									{
										logger(g_monthly_log, ",%d", -9999);
									}
								}
							}
						}
						else
						{
							for ( i = 0; i < 30; ++i )
							{
								logger(g_monthly_log, ",%d", -9999);
							}
						}
					}
				}
			}
			else
			{
				for ( i = 0; i < 31; ++i )
				{
					logger(g_monthly_log, ",%d", -9999);
				}
			}
		}
	}
	/************************************************************************/
	/* printing variables at cell level only if there's more than one layer */

	if( c->heights_count > 1 )
	{
		logger(g_monthly_log, ",***,%3.4f,%3.4f,%3.4f",
				c->monthly_gpp,
				c->monthly_npp,
				c->monthly_aut_resp);

	}
	/* printing variables at cell level also if there's more than one layer */
	else
	{
		logger(g_monthly_log,",*****");
	}
	/* printing variables only at cell level */
	logger(g_monthly_log, ",%3.2f,%3.2f,%3.2f,%3.2f\n",
			c->monthly_et,
			c->monthly_lh_flux,
			c->asw,
			c->monthly_iwue);

	/************************************************************************/


	if ( ( IS_LEAP_YEAR( c->years[year].year ) ? (MonthLength_Leap[11]) : (MonthLength[11] )) == c->doy )
	{
		++years_counter;

		if ( years_counter ==  years_of_simulation)
		{
			logger(g_monthly_log, sz_launched, netcdf_get_version(), datetime_current());
			print_model_paths(g_monthly_log);
			//const char* p;
			//p = file_get_name_only(g_monthly_log->filename);
			logger(g_monthly_log, "#output file = %s\n", g_monthly_log->filename);
			print_model_settings(g_monthly_log);
		}
	}
}

void EOY_print_output_cell_level(cell_t *const c, const int year, const int years_of_simulation )
{
	int layer;
	int height;
	int dbh;
	int age;
	int species;

	static int years_counter;

	species_t *s;

	/* return if annual logging is off*/
	if ( ! g_annual_log ) return;



	/* test */
	//note it can be used only if no other classes are added!!
	/*
	if ( !year)
	{
		c->initial_tree_layers_count = c->tree_layers_count;
		c->initial_heights_count = c->heights_count;

		for ( height = c->heights_count - 1; height >= 0 ; --height )
		{
			c->heights[height].initial_dbhs_count = c->heights[height].dbhs_count;
			for ( dbh = c->heights[height].dbhs_count - 1; dbh >= 0; --dbh )
			{
				c->heights[height].dbhs[dbh].initial_ages_count = c->heights[height].dbhs[dbh].ages_count;

				for ( age = 0; age < c->heights[height].dbhs[dbh].ages_count ; ++age )
				{
					c->heights[height].dbhs[dbh].ages[age].initial_species_count = c->heights[height].dbhs[dbh].ages[age].species_count;
				}
			}
		}
	}
	 */

	/* heading */
	if ( ! year )
	{
		logger(g_annual_log, "X,Y,YEAR");

		for ( layer = c->tree_layers_count - 1; layer >= 0; --layer )
		{
			/* heading for layers */
			logger(g_annual_log,",LAYER");

			for ( height = 0; height < c->heights_count+c->heights_avail; ++height )
			{
				if ( check_height_index(height, c->heights, c->heights_count) )
				{
					if( layer == c->heights[height].height_z )
					{
						/* heading for heights value */
						logger(g_annual_log,",HEIGHT");

						for ( dbh = 0; dbh < c->heights[height].dbhs_count+c->heights[height].dbhs_avail; ++dbh )
						{
							/* heading for dbhs value */
							logger(g_annual_log, ",DBH");

							if ( check_dbh_index(dbh, c->heights[height].dbhs, c->heights[height].dbhs_count) )
							{
								for ( age = 0; age < c->heights[height].dbhs[dbh].ages_count+c->heights[height].dbhs[dbh].ages_avail; ++age )
								{
									/* heading for ages */
									logger(g_annual_log,",AGE");

									if ( check_age_index(age, c->heights[height].dbhs[dbh].ages, c->heights[height].dbhs[dbh].ages_count)  )
									{
										for ( species = 0; species < c->heights[height].dbhs[dbh].ages[age].species_count+c->heights[height].dbhs[dbh].ages[age].species_avail; ++species )
										{
											if ( check_species_index(species, c->heights[height].dbhs[dbh].ages[age].species, c->heights[height].dbhs[dbh].ages[age].species_count)  )
											{

												/* heading for species name */
												logger(g_annual_log,",SPECIES");

												/* heading for management name */
												logger(g_annual_log,",MANAGEMENT");

												logger(g_annual_log,
														",GPP"
														",GR"
														",MR"
														",RA"
														",NPP"
														",CUE"
														",Y(perc)"
														",PeakLAI"
														",MaxLAI"
														",SAPWOOD-AREA"
														",CC-Proj"
														",DBHDC"
														",CROWN_DIAMETER"
														",CROWN_HEIGHT"
														",CROWN_AREA_PROJ"
														",APAR"
														",LiveTree"
														",DeadTree"
														",ThinnedTree"
														",VEG_D"
														",FIRST_VEG_DAY"
														",CTRANSP"
														",CINT"
														",CLE"
														",WUE"
														",MIN_RESERVE_C"
														",RESERVE_C"
														",STEM_C"
														",STEM_SAP_C"
														",STEM_HEA_C"
														",STEM_LIVE_C"
														",STEM_DEAD_C"
														",MAX_LEAF_C"
														",MAX_FROOT_C"
														",CROOT_C"
														",CROOT_LIVE_C"
														",CROOT_DEAD_C"
														",BRANCH_C"
														",BRANCH_LIVE_C"
														",BRANCH_DEAD_C"
														",FRUIT_C"
														",STANDING_WOOD"
														",DELTA_WOOD"
														",CUM_DELTA_WOOD"
														",BASAL_AREA"
														",TREE_CAI"
														",TREE_MAI"
														",CAI"
														",MAI"
														",VOLUME"
														",TREE_VOLUME"
														",DELTA_TREE_VOL(perc)"
														",DELTA_AGB"
														",DELTA_BGB"
														",AGB"
														",BGB"
														",BGB:AGB"
														",DELTA_TREE_AGB"
														",DELTA_TREE_BGB"
														",C_HWP"
														",VOLUME_HWP"
														",STEM_RA"
														",LEAF_RA"
														",FROOT_RA"
														",CROOT_RA"
														",BRANCH_RA");
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		/************************************************************************/
		/* heading cell variables */
		logger(g_annual_log,",gpp,npp,ar,et,le,soil-evapo,asw,iWue,vol,cum_vol,run_off");
		/************************************************************************/
		/* heading meteo variables */
		logger(g_annual_log,",solar_rad,tavg,tmax,tmin,tday,tnight,vpd,prcp,tsoil,rh,[CO2]");
		/************************************************************************/
		/* end heading */
		logger(g_annual_log,"\n");
		/************************************************************************/

	}

	/*****************************************************************************************************/

	/* values */
	logger(g_annual_log, "%d,%d,%d", c->x, c->y, c->years[year].year);

	/* print class level values */
	if ( c->heights_count )
	{
		for ( layer = c->tree_layers_count - 1; layer >= 0; --layer )
		{
			int i;

			/* print layer */
			logger(g_annual_log,",%d", layer);

			// ALESSIOR commented on 10/07/2017
			//qsort(c->heights, c->heights_count, sizeof(height_t), sort_by_heights_desc);

			for ( height = 0; height <c->heights_count+c->heights_avail; ++height )
			{
				if ( check_height_index(height, c->heights, c->heights_count) )
				{
					if ( layer == c->heights[height].height_z )
					{
						/* print height */
						logger(g_annual_log,",%g", c->heights[height].value);

						for ( dbh = 0; dbh < c->heights[height].dbhs_count+c->heights[height].dbhs_avail; ++dbh )
						{
							if ( check_dbh_index(dbh, c->heights[height].dbhs, c->heights[height].dbhs_count) )
							{
								/* print dbh */
								logger(g_annual_log,",%g", c->heights[height].dbhs[dbh].value);

								for ( age = 0; age < c->heights[height].dbhs[dbh].ages_count+c->heights[height].dbhs[dbh].ages_avail; ++age )
								{
									if ( check_age_index(age, c->heights[height].dbhs[dbh].ages, c->heights[height].dbhs[dbh].ages_count) )
									{
										/* print age */
										logger(g_annual_log,",%d", c->heights[height].dbhs[dbh].ages[age].value);

										for ( species = 0; species < c->heights[height].dbhs[dbh].ages[age].species_count+c->heights[height].dbhs[dbh].ages[age].species_avail; ++species )
										{
											if ( check_species_index(species, c->heights[height].dbhs[dbh].ages[age].species, c->heights[height].dbhs[dbh].ages[age].species_count) )
											{
												s  = &c->heights[height].dbhs[dbh].ages[age].species[species];

												/* print species name */
												logger(g_annual_log,",%s", c->heights[height].dbhs[dbh].ages[age].species[species].name);
                                  
												/* print management */
												logger(g_annual_log,",%c", sz_management[c->heights[height].dbhs[dbh].ages[age].species[species].management]);

												/* print variables at layer-class level */
												logger(g_annual_log,",%6.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%d,%d,%d,%d,%d,%3.4f"
														",%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f"
														",%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f",
														s->value[YEARLY_GPP],
														s->value[YEARLY_TOTAL_GROWTH_RESP],
														s->value[YEARLY_TOTAL_MAINT_RESP],
														s->value[YEARLY_TOTAL_AUT_RESP],
														s->value[YEARLY_NPP],
														s->value[YEARLY_CUE],
														s->value[YEARLY_TOTAL_AUT_RESP]/s->value[YEARLY_GPP]*100.,
														s->value[PEAK_LAI_PROJ],
														s->value[MAX_LAI_PROJ],
														s->value[SAPWOOD_AREA],
														s->value[CANOPY_COVER_PROJ],
														s->value[DBHDC_EFF],
														s->value[CROWN_DIAMETER],
														s->value[CROWN_HEIGHT],
														s->value[CROWN_AREA_PROJ],
														s->value[YEARLY_APAR],
														s->counter[N_TREE],
														s->counter[DEAD_TREE],
														s->counter[THINNED_TREE],
														s->counter[YEARLY_VEG_DAYS],
														s->counter[FIRST_VEG_DAYS],
														s->value[YEARLY_CANOPY_TRANSP],
														s->value[YEARLY_CANOPY_INT],
														s->value[YEARLY_CANOPY_LATENT_HEAT],
														s->value[YEARLY_WUE],
														s->value[MIN_RESERVE_C],
														s->value[RESERVE_C],
														s->value[STEM_C],
														s->value[STEM_SAPWOOD_C],
														s->value[STEM_HEARTWOOD_C],
														s->value[STEM_LIVEWOOD_C],
														s->value[STEM_DEADWOOD_C],
														s->value[MAX_LEAF_C],
														s->value[MAX_FROOT_C],
														s->value[CROOT_C],
														s->value[CROOT_LIVEWOOD_C],
														s->value[CROOT_DEADWOOD_C],
														s->value[BRANCH_C],
														s->value[BRANCH_LIVEWOOD_C],
														s->value[BRANCH_DEADWOOD_C],
														s->value[FRUIT_C],
														s->value[STANDING_WOOD],
														s->value[YEARLY_C_TO_WOOD],
														s->value[CUM_YEARLY_C_TO_WOOD],
														s->value[STAND_BASAL_AREA_m2],
														s->value[TREE_CAI],
														s->value[TREE_MAI],
														s->value[CAI],
														s->value[MAI],
														s->value[VOLUME],
														s->value[TREE_VOLUME],
														(s->value[TREE_CAI]/s->value[TREE_VOLUME])*100.,
														s->value[DELTA_AGB],
														s->value[DELTA_BGB],
														s->value[AGB],
														s->value[BGB],
														s->value[BGB]/s->value[AGB],
														s->value[DELTA_TREE_AGB],
														s->value[DELTA_TREE_BGB],
														s->value[C_HWP],
														s->value[VOLUME_HWP],
														s->value[YEARLY_STEM_AUT_RESP],
														s->value[YEARLY_LEAF_AUT_RESP],
														s->value[YEARLY_FROOT_AUT_RESP],
														s->value[YEARLY_CROOT_AUT_RESP],
														s->value[YEARLY_BRANCH_AUT_RESP]
												);
											}
											else
											{
												for ( i = 0; i < 68; ++i )
												{
													logger(g_annual_log, ",%d", -9999);
												}
											}
										}
									}
									else
									{
										for ( i = 0; i < 70; ++i )
										{
											logger(g_annual_log, ",%d", -9999);
										}
									}
								}
							}
							else
							{
								for ( i = 0; i < 71; ++i )
								{
									logger(g_annual_log, ",%d", -9999);
								}
							}
						}
					}
				}
				else
				{
					for ( i = 0; i < 72; ++i )
					{
						logger(g_annual_log, ",%d", -9999);
					}
				}
			}
		}
	}
	else
	{
		//ALESSIOC TO ALLESSIOR PRINT EMPTY SPACES WHEN N_TREE = 0
	}
	/************************************************************************/
	/* printing variables at cell level only if there's more than one layer */
	logger(g_annual_log, ",%3.4f,%3.4f,%3.4f,%3.4f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f",
			c->annual_gpp,
			c->annual_npp,
			c->annual_aut_resp,
			c->annual_et,
			c->annual_lh_flux,
			c->annual_soil_evapo,
			c->asw,
			c->annual_iwue,
			c->volume,
			c->cum_volume,
			c->annual_out_flow);
	/************************************************************************/
	/* print meteo variables at cell level */
	logger(g_annual_log, ",%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f",
			c->years[year].yearly_mean.solar_rad     ,
			c->years[year].yearly_mean.tavg          ,
			c->years[year].yearly_mean.tmax          ,
			c->years[year].yearly_mean.tmin          ,
			c->years[year].yearly_mean.tday          ,
			c->years[year].yearly_mean.tnight        ,
			c->years[year].yearly_mean.vpd           ,
			c->years[year].yearly_mean.prcp          ,
			c->years[year].yearly_mean.tsoil         ,
			c->years[year].yearly_mean.rh_f          ,
			c->years[year].co2Conc);
	/************************************************************************/
	/* end print */
	logger(g_annual_log,"\n");
	/************************************************************************/

	++years_counter;

	if ( years_counter ==  years_of_simulation)
	{
		g_annual_log->std_output = 1;
		logger(g_annual_log, sz_launched, netcdf_get_version(), datetime_current());
		print_model_paths(g_annual_log);
		//const char* p;
		//p = file_get_name_only(g_annual_log->filename);
		logger(g_annual_log, "#output file = %s\n", g_annual_log->filename);
		print_model_settings(g_annual_log);
	}
}


/*************************************************************************************************************************************/
/*************************************************************************************************************************************/

void EOD_print_output_soil_cell_level(cell_t *const c, const int day, const int month, const int year, const int years_of_simulation )
{
	static int years_counter;

	/* return if monthly logging is off*/
	if ( ! g_daily_soil_log ) return;

	/* heading */
	if ( !day && !month && !year )
	{
		logger(g_daily_soil_log, "X,Y,YEAR,MONTH,DAY");

		/************************************************************************/
		/* heading variables at cell level only if there's more than one layer */

		/* carbon pools */
		logger(g_daily_soil_log,
				",***,"
				"leaf_litr1C,"
				"leaf_litr2C,"
				"leaf_litr3C,"
				"leaf_litr4C,"
				"froot_litr1C,"
				"froot_litr2C,"
				"froot_litr3C,"
				"froot_litr4C,"
				"litrC,"
				"litr1C,"
				"litr2C,"
				"litr3C,"
				"litr4C,"
				"soilC,"
				"soil1C,"
				"soil2C,"
				"soil3C,"
				"soil4C");

		/* nitrogen pools */
		logger(g_daily_soil_log,
				",***,"
				"leaf_litr1N,"
				"leaf_litr2N,"
				"leaf_litr3N,"
				"leaf_litr4N,"
				"froot_litr1N,"
				"froot_litr2N,"
				"froot_litr3N,"
				"froot_litr4N,"
				"litrN,"
				"litr1N,"
				"litr2N,"
				"litr3N,"
				"litr4N,"
				"soilN,"
				"soil1N,"
				"soil2N,"
				"soil3N,"
				"soil4N\n");
	}

	/* values */
	logger(g_daily_soil_log, "%d,%d,%d,%d,%d", c->x, c->y, c->years[year].year, month + 1, day + 1);

	/************************************************************************/

	/* carbon pools */
	logger(g_daily_soil_log, ",***,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f",
			c->leaf_litr1C,
			c->leaf_litr2C,
			c->leaf_litr3C,
			c->leaf_litr4C,
			c->froot_litr1C,
			c->froot_litr2C,
			c->froot_litr3C,
			c->froot_litr4C,
			c->litrC,
			c->litr1C,
			c->litr2C,
			c->litr3C,
			c->litr4C,
			c->soilC,
			c->soil1C,
			c->soil2C,
			c->soil3C,
			c->soil4C);

	/* nitrogen pools */
	logger(g_daily_soil_log, ",***,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f\n",
			c->leaf_litr1N,
			c->leaf_litr2N,
			c->leaf_litr3N,
			c->leaf_litr4N,
			c->froot_litr1N,
			c->froot_litr2N,
			c->froot_litr3N,
			c->froot_litr4N,
			c->litrN,
			c->litr1N,
			c->litr2N,
			c->litr3N,
			c->litr4N,
			c->soilN,
			c->soil1N,
			c->soil2N,
			c->soil3N,
			c->soil4N);

	if ( c->doy == ( IS_LEAP_YEAR( c->years[year].year ) ? 366 : 365 ) )
	{
		++years_counter;

		if ( years_counter ==  years_of_simulation)
		{
			logger(g_daily_soil_log, sz_launched, netcdf_get_version(), datetime_current());
			print_model_paths(g_daily_soil_log);
			//const char* p;
			//p = file_get_name_only(g_daily_soil_log->filename);
			logger(g_daily_soil_log, "#output file = %s\n", g_daily_soil_log->filename);
			print_model_settings(g_daily_soil_log);
		}
	}
}

void EOM_print_output_soil_cell_level(cell_t *const c, const int month, const int year, const int years_of_simulation )
{
	static int years_counter;

	/* return if monthly logging is off*/
	if ( ! g_monthly_soil_log ) return;

	/* heading */
	if ( !month && !year )
	{
		logger(g_monthly_soil_log, "X,Y,YEAR,MONTH");

		/************************************************************************/
		/* heading variables at cell level only if there's more than one layer */

		/* carbon pools */
		logger(g_monthly_soil_log,
				",***,"
				"leaf_litr1C,"
				"leaf_litr2C,"
				"leaf_litr3C,"
				"leaf_litr4C,"
				"froot_litr1C,"
				"froot_litr2C,"
				"froot_litr3C,"
				"froot_litr4C,"
				"litrC,"
				"litr1C,"
				"litr2C,"
				"litr3C,"
				"litr4C,"
				"soilC,"
				"soil1C,"
				"soil2C,"
				"soil3C,"
				"soil4C");

		/* nitrogen pools */
		logger(g_monthly_soil_log,
				",***,"
				"leaf_litr1N,"
				"leaf_litr2N,"
				"leaf_litr3N,"
				"leaf_litr4N,"
				"froot_litr1N,"
				"froot_litr2N,"
				"froot_litr3N,"
				"froot_litr4N,"
				"litrN,"
				"litr1N,"
				"litr2N,"
				"litr3N,"
				"litr4N,"
				"soilN,"
				"soil1N,"
				"soil2N,"
				"soil3N,"
				"soil4N\n");
	}

	/* values */
	logger (g_monthly_soil_log, "%d,%d,%d,%d", c->years[year].year, month +1 );

	/************************************************************************/

	/* carbon pools */
	logger(g_monthly_soil_log, ",***,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f",
			c->leaf_litr1C,
			c->leaf_litr2C,
			c->leaf_litr3C,
			c->leaf_litr4C,
			c->froot_litr1C,
			c->froot_litr2C,
			c->froot_litr3C,
			c->froot_litr4C,
			c->litrC,
			c->litr1C,
			c->litr2C,
			c->litr3C,
			c->litr4C,
			c->soilC,
			c->soil1C,
			c->soil2C,
			c->soil3C,
			c->soil4C);

	/* nitrogen pools */
	logger(g_monthly_soil_log, ",***,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f\n",
			c->leaf_litr1N,
			c->leaf_litr2N,
			c->leaf_litr3N,
			c->leaf_litr4N,
			c->froot_litr1N,
			c->froot_litr2N,
			c->froot_litr3N,
			c->froot_litr4N,
			c->litrN,
			c->litr1N,
			c->litr2N,
			c->litr3N,
			c->litr4N,
			c->soilN,
			c->soil1N,
			c->soil2N,
			c->soil3N,
			c->soil4N);

	if ( ( IS_LEAP_YEAR( c->years[year].year ) ? (MonthLength_Leap[11]) : (MonthLength[11] )) == c->doy )
	{
		++years_counter;

		if ( years_counter ==  years_of_simulation)
		{
			logger(g_monthly_soil_log, sz_launched, netcdf_get_version(), datetime_current());
			print_model_paths(g_monthly_soil_log);
			//const char* p;
			//p = file_get_name_only(g_monthly_soil_log->filename);
			logger(g_monthly_soil_log, "#output file = %s\n", g_monthly_soil_log->filename);
			print_model_settings(g_monthly_soil_log);
		}
	}
}

void EOY_print_output_soil_cell_level(cell_t *const c, const int year, const int years_of_simulation )
{
	static int years_counter;

	/* return if annual logging is off*/
	if ( ! g_annual_soil_log ) return;

	/* heading */
	if ( ! year )
	{
		logger(g_annual_soil_log, "X,Y,YEAR");

		/************************************************************************/
		/* heading variables at cell level only if there's more than one layer */

		/* carbon pools */
		logger(g_annual_soil_log,
				",***,"
				"leaf_litr1C,"
				"leaf_litr2C,"
				"leaf_litr3C,"
				"leaf_litr4C,"
				"froot_litr1C,"
				"froot_litr2C,"
				"froot_litr3C,"
				"froot_litr4C,"
				"litrC,"
				"litr1C,"
				"litr2C,"
				"litr3C,"
				"litr4C,"
				"soilC,"
				"soil1C,"
				"soil2C,"
				"soil3C,"
				"soil4C");

		/* nitrogen pools */
		logger(g_annual_soil_log,
				",***,"
				"leaf_litr1N,"
				"leaf_litr2N,"
				"leaf_litr3N,"
				"leaf_litr4N,"
				"froot_litr1N,"
				"froot_litr2N,"
				"froot_litr3N,"
				"froot_litr4N,"
				"litrN,"
				"litr1N,"
				"litr2N,"
				"litr3N,"
				"litr4N,"
				"soilN,"
				"soil1N,"
				"soil2N,"
				"soil3N,"
				"soil4N\n");
	}
	/*****************************************************************************************************/

	/* values */
	logger(g_annual_soil_log, "%d,%d,%d", c->x, c->y, c->years[year].year);

	/************************************************************************/

	/* carbon pools */
	logger(g_annual_soil_log, ",***,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f",
			c->leaf_litr1C,
			c->leaf_litr2C,
			c->leaf_litr3C,
			c->leaf_litr4C,
			c->froot_litr1C,
			c->froot_litr2C,
			c->froot_litr3C,
			c->froot_litr4C,
			c->litrC,
			c->litr1C,
			c->litr2C,
			c->litr3C,
			c->litr4C,
			c->soilC,
			c->soil1C,
			c->soil2C,
			c->soil3C,
			c->soil4C);
	/* nitrogen pools */
	logger(g_annual_soil_log, ",***,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f\n",
			c->leaf_litr1N,
			c->leaf_litr2N,
			c->leaf_litr3N,
			c->leaf_litr4N,
			c->froot_litr1N,
			c->froot_litr2N,
			c->froot_litr3N,
			c->froot_litr4N,
			c->litrN,
			c->litr1N,
			c->litr2N,
			c->litr3N,
			c->litr4N,
			c->soilN,
			c->soil1N,
			c->soil2N,
			c->soil3N,
			c->soil4N);

	/************************************************************************/

	++years_counter;

	if ( years_counter == years_of_simulation )
	{
		g_annual_soil_log->std_output = 1;
		logger(g_annual_soil_log, sz_launched, netcdf_get_version(), datetime_current());
		print_model_paths(g_annual_soil_log);
		//const char* p;
		//p = file_get_name_only(g_annual_soil_log->filename);
		logger(g_annual_log, "#output file = %s\n", g_annual_soil_log->filename);
		print_model_settings(g_annual_soil_log);
	}
}
#else    //USE_NEW_OUTPUT //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "constants.h"
#include "print_output.h"
#include "common.h"
#include "settings.h"
#include "soil_settings.h"
#include "logger.h"
#include "g-function.h"
#include "nc.h"

extern settings_t* g_settings;
extern logger_t* g_daily_log;
extern logger_t* g_monthly_log;
extern logger_t* g_annual_log;
extern logger_t* g_daily_soil_log;
extern logger_t* g_monthly_soil_log;
extern logger_t* g_annual_soil_log;
extern char *g_sz_dataset_file;
extern char *g_sz_soil_file;
extern char *g_sz_input_met_file;
extern char *g_sz_settings_file;
extern char *g_sz_topo_file;
extern char *g_sz_co2_conc_file;

extern const char sz_launched[];

//extern int MonthLength [];
//extern int MonthLength_Leap [];

static const char sz_management[] = "TCN";

//extern const int CO2_MODIFIER;


static const char* get_filename(const char *const s)
{
	const char *p;
	const char *p2;

	p = NULL;

	if ( s ) {
		p = strrchr(s, '/');
		if ( p ) ++p;
		p2 = strrchr(s, '\\');
		if ( p2 ) ++p2;
		if ( p2 > p ) p = p2;
		if ( ! p ) p = s;
	}

	return p;
}

void print_model_paths(logger_t *const _log)
{
	//assert(_log);

	logger(_log, "\n\n#site: %s\n", g_settings->sitename);
	if ( g_sz_dataset_file )
		logger(_log, "#--input files--\n");
	logger(_log, "#input file = %s\n", get_filename(g_sz_dataset_file));
	logger(_log, "#soil file = %s\n", get_filename(g_sz_soil_file));
	logger(_log, "#topo file = %s\n", get_filename(g_sz_topo_file));
	logger(_log, "#met file = %s\n", get_filename(g_sz_input_met_file));
	logger(_log, "#settings file = %s\n", get_filename(g_sz_settings_file));
	if ( g_settings->CO2_trans )
		logger(_log, "#CO2 file = %s\n", get_filename(g_sz_co2_conc_file));
}

void print_model_settings(logger_t*const log)
{
	//assert(log);

	logger(log, "#--model settings--\n");

	if ( !g_settings->PSN_mod )
	{
		logger(log, "#PSN_mod = Faquhar von Caemmerer and Berry (FvCB)\n");
	}
	else
	{
		logger(log, "#PSN_mod = Monteith (LUE)\n");

//		if ( !c->CO2_modifier )
//		{
//			logger(log, "#CO2 modifier = Wang et al' method\n");
//		}
//		else
//		{
//			logger(log, "#CO2 modifier = Veroustraete et al' method\n");
//		}
	}
	logger(log, "#CO2 trans = %s\n", (CO2_TRANS_VAR == g_settings->CO2_trans) ? "var" : (CO2_TRANS_ON == g_settings->CO2_trans) ? "on" : "off");

	if ( CO2_TRANS_OFF == g_settings->CO2_trans )
	{
		logger(log, "#fixed co2 concentration = %g ppmv\n", g_settings->co2Conc);
	}
	else if ( CO2_TRANS_VAR == g_settings->CO2_trans )
	{
		logger(log, "#year %d at which co2 concentration is fixed at value = %g ppmv\n", g_settings->year_start_co2_fixed, g_settings->co2Conc);
	}

	logger(log, "#Photo accl = %s\n", g_settings->Photo_accl ? "on" : "off");
	logger(log, "#Resp  accl = %s\n", g_settings->Resp_accl ? "on" : "off");
	logger(log, "#regeneration = %s\n", g_settings->regeneration ? "on" : "off");
	logger(log, "#Management = %s\n", (MANAGEMENT_VAR == g_settings->management) ? "var" : (MANAGEMENT_ON == g_settings->management) ? "on" : "off");
	if ( g_settings->management )
	{
		logger(log, "#Year Start Management = %d\n", g_settings->year_start_management);
	}
	if ( g_settings->year_restart != -1 )
	{
		logger(log, "#Year restart = %d\n", g_settings->year_restart);
	}
	else
	{
		logger(log, "#Year restart = off\n");
	}
}

void EOD_print_output_cell_level(cell_t *const c, const int day, const int month, const int year, const int years_of_simulation )
{
	int layer;
	int height;
	int dbh;
	int age;
	int species;

	static int print_header = 0;

	species_t *s;

	/* return if daily logging is off*/
	if ( ! g_daily_log ) return;

	/* heading */
	if ( ! print_header )
	{
		print_header = 1;

		logger(g_daily_log, "X,Y,YEAR,MONTH,DAY");

		for ( layer = c->tree_layers_count - 1; layer >= 0; --layer )
		{
			for ( height = c->heights_count - 1; height >= 0 ; --height )
			{
				if( layer == c->heights[height].height_z )
				{
					/* heading for layers */
					logger(g_daily_log, ",LAYER");

					/* heading for heights value */
					logger(g_daily_log,",HEIGHT");

					for ( dbh = c->heights[height].dbhs_count - 1; dbh >= 0; --dbh )
					{
						logger(g_daily_log,",DBH");

						for ( age = 0; age < c->heights[height].dbhs[dbh].ages_count ; ++age )
						{
							/* heading for ages */
							logger(g_daily_log,",AGE");

							for ( species = 0; species < c->heights[height].dbhs[dbh].ages[age].species_count; ++species )
							{
								/* heading for species name */
								logger(g_daily_log,",SPECIES");

								/* heading for management */
								logger(g_daily_log, ",MANAGEMENT");

								logger(g_daily_log,
										",GPP"
										",Av_TOT"
										",Aj_TOT"
										",A_TOT"
										",RG"
										",RM"
										",RA"
										",NPP"
										",BP"
										",CUE"
										",BPE"
										",LAI_PROJ"
										",PEAK-LAI_PROJ"
										",LAI_EXP"
										",D-CC_P"
										",DBHDC"
										",CROWN_AREA_PROJ"
										",PAR"
										",APAR"
										",fAPAR"
										",NTREE"
										",VEG_D"
										",INT"
										",WAT"
										",EVA"
										",TRA"
										",ET"
										",LE"
										",WUE"
										",RESERVE_C"
										",STEM_C"
										",STEMSAP_C"
										",STEMLIVE_C"
										",STEMDEAD_C"
										",LEAF_C"
										",FROOT_C"
										",CROOT_C"
										",CROOTSAP_C"
										",CROOTLIVE_C"
										",CROOTDEAD_C"
										",BRANCH_C"
										",BRANCHSAP_C"
										",BRANCHLIVE_C"
										",BRANCHDEAD_C"
										",FRUIT_C"
										",TOT_SAPWOOD_C"
										",DELTA_RESERVE_C"
										",DELTA_STEM_C"
										",DELTA_LEAF_C"
										",DELTA_FROOT_C"
										",DELTA_CROOT_C"
										",DELTA_BRANCH_C"
										",DELTA_FRUIT_C"
										",RESERVE_N"
										",STEM_N"
										",STEMLIVE_N"
										",STEMDEAD_N"
										",LEAF_N"
										",FROOT_N"
										",CROOT_N"
										",CROOTLIVE_N"
										",CROOTDEAD_N"
										",BRANCH_N"
										",BRANCHLIVE_N"
										",BRANCHDEAD_N"
										",FRUIT_N"
										",DELTA_RESERVE_N"
										",DELTA_STEM_N"
										",DELTA_LEAF_N"
										",DELTA_FROOT_N"
										",DELTA_CROOT_N"
										",DELTA_BRANCH_N"
										",DELTA_FRUIT_N"
										",STEM_AR"
										",LEAF_AR"
										",FROOT_AR"
										",CROOT_AR"
										",BRANCH_AR"
										",F_CO2"
										",F_CO2_VER"
										",F_CO2_FRA"
										",FCO2_TR"
										",FLIGHT"
										",FAGE"
										",FT"
										",FVPD"
										",FN"
										",FSW"
										",LITR_C"
										",CWD_C"
								);

							}
							if ( c->heights[height].dbhs[dbh].ages[age].species_count > 1 ) {
								logger(g_daily_log,",*");
							}
						}
						if ( c->heights[height].dbhs[dbh].ages_count > 1 ) {
							logger(g_daily_log, ",**");
						}
					}
					if ( c->heights[height].dbhs_count > 1 ) {
						logger(g_daily_log,",***");
					}
				}
				if ( c->tree_layers[layer].layer_n_height_class > 1 ) {
					logger(g_daily_log,",****");
				}
			}
			if ( c->tree_layers_count > 1 ) {
				logger(g_daily_log,",*****");
			}
			/*
							if ( c->heights[height].dbhs[dbh].ages[age].species_count > 1 ) logger(g_daily_log,",*");
						}
						if ( c->heights[height].dbhs[dbh].ages_count > 1 ) logger(g_daily_log,",**");
					}
					if ( c->heights[height].dbhs_count > 1 ) logger(g_daily_log,",***");
				}
				if ( c->tree_layers[layer].layer_n_height_class > 1 ) logger(g_daily_log,",****");
			}
			if ( c->tree_layers_count > 1 ) logger(g_daily_log,",*****");
			 */
		}
		/************************************************************************/
                if (!g_settings->regeneration)
                {  
		/* heading variables only at cell level */
		logger(g_daily_log,",gpp,npp,ar,hr,rsoil,reco,nee,nep,et,le,soil_evapo,snow_pack,asw,moist_ratio,iWue,"
				"litrC,litr1C,litr2C,litr3C,litr4C,deadwoodC,deadwood2C,deadwood3C,deadwood4C,soilC,soil1C,soil2C,soil3C,soil4C,"
				"litrN,litr1N,litr2N,litr3N,litr4N,deadwoodN,deadwood2N,deadwood3N,deadwood4N,soilN,soil1N,soil2N,soil3N,soil4N,"
				"Tsoil,Daylength");

                }
               	/* end heading */
		logger(g_daily_log,"\n");
	}
	/*****************************************************************************************************/



	/* values */
	logger(g_daily_log, "%d,%d,%d,%d,%d", c->x, c->y, c->years[year].year, month + 1, day + 1);

	/* print class level values */
	for ( layer = c->tree_layers_count - 1; layer >= 0; --layer )
	{
		qsort(c->heights, c->heights_count, sizeof(height_t), sort_by_heights_desc);

		for ( height = 0; height < c->heights_count; ++height )
		{
			if( layer == c->heights[height].height_z )
			{
				/* print layer */
				logger(g_daily_log,",%d", layer);

				/* print height */
				logger(g_daily_log,",%g", c->heights[height].value);

				for ( dbh = c->heights[height].dbhs_count - 1; dbh >= 0; --dbh )
				{
					/* print dbh */
					logger(g_daily_log,",%g", c->heights[height].dbhs[dbh].value);

					for ( age = 0; age < c->heights[height].dbhs[dbh].ages_count ; ++age )
					{
						/* print age */
						logger(g_daily_log,",%d", c->heights[height].dbhs[dbh].ages[age].value);

						for ( species = 0; species < c->heights[height].dbhs[dbh].ages[age].species_count; ++species )
						{
							s  = &c->heights[height].dbhs[dbh].ages[age].species[species];


							/* print species name */
							logger(g_daily_log,",%s", c->heights[height].dbhs[dbh].ages[age].species[species].name);

							/* print management */
							logger(g_daily_log,",%c", sz_management[c->heights[height].dbhs[dbh].ages[age].species[species].management]);

							/* print variables at layer-class level */

							logger(g_daily_log,",%6.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%d,"
									"%d,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f"
									",%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f"
									",%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f",
									s->value[GPP],
									s->value[Av_TOT],
									s->value[Aj_TOT],
									s->value[A_TOT],
									s->value[TOTAL_GROWTH_RESP],
									s->value[TOTAL_MAINT_RESP],
									s->value[TOTAL_AUT_RESP],
									s->value[NPP],
									s->value[BP],
									s->value[CUE],
									s->value[BPE],
									s->value[LAI_PROJ],
									s->value[PEAK_LAI_PROJ],
									s->value[LAI_EXP],
									s->value[DAILY_CANOPY_COVER_PROJ],
									s->value[DBHDC_EFF],
									s->value[CROWN_AREA_PROJ],
									s->value[PAR],
									s->value[APAR],
									s->value[fAPAR],
									s->counter[N_TREE],
									s->counter[VEG_DAYS],
									s->value[CANOPY_INT_RAIN],
									s->value[CANOPY_WATER],
									s->value[CANOPY_EVAPO],
									s->value[CANOPY_TRANSP],
									s->value[CANOPY_EVAPO_TRANSP],
									s->value[CANOPY_LATENT_HEAT],
									s->value[WUE],
									s->value[RESERVE_C],
									s->value[STEM_C],
									s->value[STEM_SAPWOOD_C],
									s->value[STEM_LIVEWOOD_C],
									s->value[STEM_DEADWOOD_C],
									s->value[LEAF_C],
									s->value[FROOT_C],
									s->value[CROOT_C],
									s->value[CROOT_SAPWOOD_C],
									s->value[CROOT_LIVEWOOD_C],
									s->value[CROOT_DEADWOOD_C],
									s->value[BRANCH_C],
									s->value[BRANCH_SAPWOOD_C],
									s->value[BRANCH_LIVEWOOD_C],
									s->value[BRANCH_DEADWOOD_C],
									s->value[FRUIT_C],
									s->value[TOT_SAPWOOD_C],
									s->value[C_TO_RESERVE],
									s->value[C_TO_STEM],
									s->value[C_TO_LEAF],
									s->value[C_TO_FROOT],
									s->value[C_TO_CROOT],
									s->value[C_TO_BRANCH],
									s->value[C_TO_FRUIT],
									s->value[RESERVE_N],
									s->value[STEM_N],
									s->value[STEM_LIVEWOOD_N],
									s->value[STEM_DEADWOOD_N],
									s->value[LEAF_N],
									s->value[FROOT_N],
									s->value[CROOT_N],
									s->value[CROOT_LIVEWOOD_N],
									s->value[CROOT_DEADWOOD_N],
									s->value[BRANCH_N],
									s->value[BRANCH_LIVEWOOD_N],
									s->value[BRANCH_DEADWOOD_N],
									s->value[FRUIT_N],
									s->value[N_TO_RESERVE],
									s->value[N_TO_STEM],
									s->value[N_TO_LEAF],
									s->value[N_TO_FROOT],
									s->value[N_TO_CROOT],
									s->value[N_TO_BRANCH],
									s->value[N_TO_FRUIT],
									s->value[STEM_AUT_RESP],
									s->value[LEAF_AUT_RESP],
									s->value[FROOT_AUT_RESP],
									s->value[CROOT_AUT_RESP],
									s->value[BRANCH_AUT_RESP],
									s->value[F_CO2],
									s->value[F_CO2_VER],
									s->value[F_CO2_FRANKS],
									s->value[F_CO2_TR],
									s->value[F_LIGHT_MAKELA],
									s->value[F_AGE],
									s->value[F_T],
									s->value[F_VPD],
									s->value[F_NUTR],
									s->value[F_SW],
									s->value[LITR_C],
									s->value[CWD_C]);
                                                                        //s->value[psi]);
						}

						if ( c->heights[height].dbhs[dbh].ages[age].species_count > 1 ) {
							logger(g_daily_log,",*");
						}
					}
					if ( c->heights[height].dbhs[dbh].ages_count > 1 ) {
						logger(g_daily_log,",**");
					}
				}
				if ( c->heights[height].dbhs_count > 1 ) {
					logger(g_daily_log,",***");
				}
			}
			if ( c->tree_layers[layer].layer_n_height_class > 1 ) {
				logger(g_daily_log,",****");
			}
		}
		if ( c->tree_layers_count > 1 ) {
			logger(g_daily_log,",*****");
		}
	}
	/************************************************************************/
        if (!g_settings->regeneration)
        { 
	/* printing variables only at cell level */

	logger(g_daily_log, ",%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,"
			"%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f",
			c->daily_gpp,
			c->daily_npp,
			c->daily_aut_resp,
			c->daily_het_resp,
			c->daily_soil_resp,
			c->daily_r_eco,
			c->daily_nee,
			c->daily_nep,
			c->daily_et,
			c->daily_lh_flux,
			c->daily_soil_evapo,
			c->snow_pack,
			c->asw,
			c->soil_moist_ratio,
			c->daily_iwue,
			c->litrC,
			c->litr1C,
			c->litr2C,
			c->litr3C,
			c->litr4C,
			c->cwd_C,
			c->cwd_2C,
			c->cwd_3C,
			c->cwd_4C,
			c->soilC,
			c->soil1C,
			c->soil2C,
			c->soil3C,
			c->soil4C,
			c->litrN,
			c->litr1N,
			c->litr2N,
			c->litr3N,
			c->litr4N,
			c->cwd_N,
			c->cwd_2N,
			c->cwd_3N,
			c->cwd_4N,
			c->soilN,
			c->soil1N,
			c->soil2N,
			c->soil3N,
			c->soil4N,
			c->years[year].m[month].d[day].tsoil,
			c->years[year].m[month].d[day].daylength
	);
        }
	/************************************************************************/
        /* end heading */
	logger(g_daily_log,"\n"); 
}

void EOD_print_output_cell_level_ddalmo(cell_t *const c, const int day, const int month, const int year, const int years_of_simulation )
{
	int layer;
	int height;
	int dbh;
	int age;
	int species;

	static int print_header = 0;

	species_t *s;

	/* return if daily logging is off*/
	if ( ! g_daily_log ) return;

	/* heading */
	if ( ! print_header )
	{
		print_header = 1;

		//logger(g_daily_log, "X,Y,YEAR,MONTH,DAY");

		for ( layer = c->tree_layers_count - 1; layer >= 0; --layer )
		{
                      	logger(g_daily_log, "X,Y,YEAR,MONTH,DAY");
  
			for ( height = c->heights_count - 1; height >= 0 ; --height )
			{
				if( layer == c->heights[height].height_z )
				{
					/* heading for layers */
					logger(g_daily_log, ",LAYER");

					/* heading for heights value */
					logger(g_daily_log,",HEIGHT");

					for ( dbh = c->heights[height].dbhs_count - 1; dbh >= 0; --dbh )
					{
						logger(g_daily_log,",DBH");

						for ( age = 0; age < c->heights[height].dbhs[dbh].ages_count ; ++age )
						{
							/* heading for ages */
							logger(g_daily_log,",AGE");

							for ( species = 0; species < c->heights[height].dbhs[dbh].ages[age].species_count; ++species )
							{
								/* heading for species name */
								logger(g_daily_log,",SPECIES");

								/* heading for management */
								logger(g_daily_log, ",MANAGEMENT");

								logger(g_daily_log,
										",GPP"
										",Av_TOT"
										",Aj_TOT"
										",A_TOT"
										",RG"
										",RM"
										",RA"
										",NPP"
										",BP"
										",CUE"
										",BPE"
										",LAI_PROJ"
										",PEAK-LAI_PROJ"
										",LAI_EXP"
										",D-CC_P"
										",DBHDC"
										",CROWN_AREA_PROJ"
										",PAR"
										",APAR"
										",fAPAR"
										",NTREE"
										",VEG_D"
										",INT"
										",WAT"
										",EVA"
										",TRA"
										",ET"
										",LE"
										",WUE"
										",RESERVE_C"
										",STEM_C"
										",STEMSAP_C"
										",STEMLIVE_C"
										",STEMDEAD_C"
										",LEAF_C"
										",FROOT_C"
										",CROOT_C"
										",CROOTSAP_C"
										",CROOTLIVE_C"
										",CROOTDEAD_C"
										",BRANCH_C"
										",BRANCHSAP_C"
										",BRANCHLIVE_C"
										",BRANCHDEAD_C"
										",FRUIT_C"
										",TOT_SAPWOOD_C"
										",DELTA_RESERVE_C"
										",DELTA_STEM_C"
										",DELTA_LEAF_C"
										",DELTA_FROOT_C"
										",DELTA_CROOT_C"
										",DELTA_BRANCH_C"
										",DELTA_FRUIT_C"
										",RESERVE_N"
										",STEM_N"
										",STEMLIVE_N"
										",STEMDEAD_N"
										",LEAF_N"
										",FROOT_N"
										",CROOT_N"
										",CROOTLIVE_N"
										",CROOTDEAD_N"
										",BRANCH_N"
										",BRANCHLIVE_N"
										",BRANCHDEAD_N"
										",FRUIT_N"
										",DELTA_RESERVE_N"
										",DELTA_STEM_N"
										",DELTA_LEAF_N"
										",DELTA_FROOT_N"
										",DELTA_CROOT_N"
										",DELTA_BRANCH_N"
										",DELTA_FRUIT_N"
										",STEM_AR"
										",LEAF_AR"
										",FROOT_AR"
										",CROOT_AR"
										",BRANCH_AR"
										",F_CO2"
										",F_CO2_VER"
										",F_CO2_FRA"
										",FCO2_TR"
										",FLIGHT"
										",FAGE"
										",FT"
										",FT_ACCL"
										",FVPD"
										",FN"
										",FSW"
										",LITR_C"
										",CWD_C"
								);

							}
							if ( c->heights[height].dbhs[dbh].ages[age].species_count > 1 ) {
								logger(g_daily_log,",*");
							}
						}
						if ( c->heights[height].dbhs[dbh].ages_count > 1 ) {
							logger(g_daily_log, ",**");
						}
					}
					if ( c->heights[height].dbhs_count > 1 ) {
						logger(g_daily_log,",***");
					}
				}
				if ( c->tree_layers[layer].layer_n_height_class > 1 ) {
					logger(g_daily_log,",****");
				}
			}
			//if ( c->tree_layers_count > 1 ) {
		        //		logger(g_daily_log,",*****");
			//}
			/*
							if ( c->heights[height].dbhs[dbh].ages[age].species_count > 1 ) logger(g_daily_log,",*");
						}
						if ( c->heights[height].dbhs[dbh].ages_count > 1 ) logger(g_daily_log,",**");
					}
					if ( c->heights[height].dbhs_count > 1 ) logger(g_daily_log,",***");
				}
				if ( c->tree_layers[layer].layer_n_height_class > 1 ) logger(g_daily_log,",****");
			}
			if ( c->tree_layers_count > 1 ) logger(g_daily_log,",*****");
			 */
		//}

		  /************************************************************************/
                  if (!g_settings->regeneration)
                  {  
		  /* heading variables only at cell level */
		  logger(g_daily_log,",gpp,npp,ar,hr,rsoil,reco,nee,nep,et,le,soil_evapo,snow_pack,asw,moist_ratio,iWue,"
				"litrC,litr1C,litr2C,litr3C,litr4C,deadwoodC,deadwood2C,deadwood3C,deadwood4C,soilC,soil1C,soil2C,soil3C,soil4C,"
				"litrN,litr1N,litr2N,litr3N,litr4N,deadwoodN,deadwood2N,deadwood3N,deadwood4N,soilN,soil1N,soil2N,soil3N,soil4N,"
				"Tsoil,Daylength");

                  }
               	 /* end heading */
		 logger(g_daily_log,"\n");
                } //loop on layer   
	}
	/*****************************************************************************************************/



	/* values */
	//logger(g_daily_log, "%d,%d,%d,%d,%d", c->x, c->y, c->years[year].year, month + 1, day + 1);

	/* print class level values */
	for ( layer = c->tree_layers_count - 1; layer >= 0; --layer )
	{
                /* values */
        	logger(g_daily_log, "%d,%d,%d,%d,%d", c->x, c->y, c->years[year].year, month + 1, day + 1);

		qsort(c->heights, c->heights_count, sizeof(height_t), sort_by_heights_desc);

		for ( height = 0; height < c->heights_count; ++height )
		{
			if( layer == c->heights[height].height_z )
			{
				/* print layer */
				logger(g_daily_log,",%d", layer);

				/* print height */
				logger(g_daily_log,",%g", c->heights[height].value);

				for ( dbh = c->heights[height].dbhs_count - 1; dbh >= 0; --dbh )
				{
					/* print dbh */
					logger(g_daily_log,",%g", c->heights[height].dbhs[dbh].value);

					for ( age = 0; age < c->heights[height].dbhs[dbh].ages_count ; ++age )
					{
						/* print age */
						logger(g_daily_log,",%d", c->heights[height].dbhs[dbh].ages[age].value);

						for ( species = 0; species < c->heights[height].dbhs[dbh].ages[age].species_count; ++species )
						{
							s  = &c->heights[height].dbhs[dbh].ages[age].species[species];


							/* print species name */
							logger(g_daily_log,",%s", c->heights[height].dbhs[dbh].ages[age].species[species].name);

							/* print management */
							logger(g_daily_log,",%c", sz_management[c->heights[height].dbhs[dbh].ages[age].species[species].management]);

							/* print variables at layer-class level */
					
							logger(g_daily_log,",%6.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%d,"
									"%d,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f"
									",%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f"
									",%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f",
									s->value[GPP],
									s->value[Av_TOT],
									s->value[Aj_TOT],
									s->value[A_TOT],
									s->value[TOTAL_GROWTH_RESP],
									s->value[TOTAL_MAINT_RESP],
									s->value[TOTAL_AUT_RESP],
									s->value[NPP],
									s->value[BP],
									s->value[CUE],
									s->value[BPE],
									s->value[LAI_PROJ],
									s->value[PEAK_LAI_PROJ],
									s->value[LAI_EXP],
									s->value[DAILY_CANOPY_COVER_PROJ],
									s->value[DBHDC_EFF],
									s->value[CROWN_AREA_PROJ],
									s->value[PAR],
									s->value[APAR],
									s->value[fAPAR],
									s->counter[N_TREE],
									s->counter[VEG_DAYS],
									s->value[CANOPY_INT_RAIN],
									s->value[CANOPY_WATER],
									s->value[CANOPY_EVAPO],
									s->value[CANOPY_TRANSP],
									s->value[CANOPY_EVAPO_TRANSP],
									s->value[CANOPY_LATENT_HEAT],
									s->value[WUE],
									s->value[RESERVE_C],
									s->value[STEM_C],
									s->value[STEM_SAPWOOD_C],
									s->value[STEM_LIVEWOOD_C],
									s->value[STEM_DEADWOOD_C],
									s->value[LEAF_C],
									s->value[FROOT_C],
									s->value[CROOT_C],
									s->value[CROOT_SAPWOOD_C],
									s->value[CROOT_LIVEWOOD_C],
									s->value[CROOT_DEADWOOD_C],
									s->value[BRANCH_C],
									s->value[BRANCH_SAPWOOD_C],
									s->value[BRANCH_LIVEWOOD_C],
									s->value[BRANCH_DEADWOOD_C],
									s->value[FRUIT_C],
									s->value[TOT_SAPWOOD_C],
									s->value[C_TO_RESERVE],
									s->value[C_TO_STEM],
									s->value[C_TO_LEAF],
									s->value[C_TO_FROOT],
									s->value[C_TO_CROOT],
									s->value[C_TO_BRANCH],
									s->value[C_TO_FRUIT],
									s->value[RESERVE_N],
									s->value[STEM_N],
									s->value[STEM_LIVEWOOD_N],
									s->value[STEM_DEADWOOD_N],
									s->value[LEAF_N],
									s->value[FROOT_N],
									s->value[CROOT_N],
									s->value[CROOT_LIVEWOOD_N],
									s->value[CROOT_DEADWOOD_N],
									s->value[BRANCH_N],
									s->value[BRANCH_LIVEWOOD_N],
									s->value[BRANCH_DEADWOOD_N],
									s->value[FRUIT_N],
									s->value[N_TO_RESERVE],
									s->value[N_TO_STEM],
									s->value[N_TO_LEAF],
									s->value[N_TO_FROOT],
									s->value[N_TO_CROOT],
									s->value[N_TO_BRANCH],
									s->value[N_TO_FRUIT],
									s->value[STEM_AUT_RESP],
									s->value[LEAF_AUT_RESP],
									s->value[FROOT_AUT_RESP],
									s->value[CROOT_AUT_RESP],
									s->value[BRANCH_AUT_RESP],
									s->value[F_CO2],
									s->value[F_CO2_VER],
									s->value[F_CO2_FRANKS],
									s->value[F_CO2_TR],
									s->value[F_LIGHT_MAKELA],
									s->value[F_AGE],
									s->value[F_T],
									s->value[F_ACCL],
									s->value[F_VPD],
									s->value[F_NUTR],
									s->value[F_SW],
									s->value[LITR_C],
									s->value[CWD_C]);
                                                                        //s->value[psi]);
						}

						if ( c->heights[height].dbhs[dbh].ages[age].species_count > 1 ) {
							logger(g_daily_log,",*");
						}
					}
					if ( c->heights[height].dbhs[dbh].ages_count > 1 ) {
						logger(g_daily_log,",**");
					}
				}
				if ( c->heights[height].dbhs_count > 1 ) {
					logger(g_daily_log,",***");
				}
			}
			if ( c->tree_layers[layer].layer_n_height_class > 1 ) {
				logger(g_daily_log,",****");
			}
		}
		//if ( c->tree_layers_count > 1 ) {
		//	logger(g_daily_log,",*****");
		//}
	 //}
	  /************************************************************************/
          if (!g_settings->regeneration)
          { 
	  /* printing variables only at cell level */

	  logger(g_daily_log, ",%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,"
			"%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f",
			c->daily_gpp,
			c->daily_npp,
			c->daily_aut_resp,
			c->daily_het_resp,
			c->daily_soil_resp,
			c->daily_r_eco,
			c->daily_nee,
			c->daily_nep,
			c->daily_et,
			c->daily_lh_flux,
			c->daily_soil_evapo,
			c->snow_pack,
			c->asw,
			c->soil_moist_ratio,
			c->daily_iwue,
			c->litrC,
			c->litr1C,
			c->litr2C,
			c->litr3C,
			c->litr4C,
			c->cwd_C,
			c->cwd_2C,
			c->cwd_3C,
			c->cwd_4C,
			c->soilC,
			c->soil1C,
			c->soil2C,
			c->soil3C,
			c->soil4C,
			c->litrN,
			c->litr1N,
			c->litr2N,
			c->litr3N,
			c->litr4N,
			c->cwd_N,
			c->cwd_2N,
			c->cwd_3N,
			c->cwd_4N,
			c->soilN,
			c->soil1N,
			c->soil2N,
			c->soil3N,
			c->soil4N,
			c->years[year].m[month].d[day].tsoil,
			c->years[year].m[month].d[day].daylength
	   );
           }
	   /************************************************************************/
           /* end heading */
	   logger(g_daily_log,"\n"); 
        }   // loop on layer     
}

// currently used in the 5p6 
void EOD_print_output_cell_level_mc(cell_t *const c, const int day, const int month, const int year, const int years_of_simulation )
{
	int layer;
	int height;
	int dbh;
	int age;
	int species;

	static int print_header = 0;
	
	int print_cell = 0;   // print cell level value   

	species_t *s;

	/* return if daily logging is off*/
	if ( ! g_daily_log ) return;

	/* heading */
	if ( ! print_header )
	{
		print_header = 1;

	
                      	logger(g_daily_log, "X,Y,YEAR,MONTH,DAY");
  
			
					/* heading for layers */
					logger(g_daily_log, ",LAYER");

					/* heading for heights value */
					logger(g_daily_log,",HEIGHT");

					
						logger(g_daily_log,",DBH");

						
							/* heading for ages */
							logger(g_daily_log,",AGE");

								/* heading for species name */
								logger(g_daily_log,",SPECIES");

								/* heading for management */
								logger(g_daily_log, ",MANAGEMENT");

								logger(g_daily_log,
										",GPP"
										",Av_TOT"
										",Aj_TOT"
										",A_TOT"
										",RG"
										",RM"
										",RA"
										",NPP"
										",BP"
										",CUE"
										",BPE"
										",LAI_PROJ"
										",PEAK-LAI_PROJ"
										",LAI_EXP"
										",D-CC_P"
										",DBHDC"
										",CROWN_AREA_PROJ"
										",PAR"
										",APAR"
										",fAPAR"
										",NTREE"
										",VEG_D"
										",INT"
										",WAT"
										",EVA"
										",TRA"
										",ET"
										",LE"
										",WUE"
										",RESERVE_C"
										",STEM_C"
										",STEMSAP_C"
										",STEMLIVE_C"
										",STEMDEAD_C"
										",LEAF_C"
										",FROOT_C"
										",CROOT_C"
										",CROOTSAP_C"
										",CROOTLIVE_C"
										",CROOTDEAD_C"
										",BRANCH_C"
										",BRANCHSAP_C"
										",BRANCHLIVE_C"
										",BRANCHDEAD_C"
										",FRUIT_C"
										",TOT_SAPWOOD_C"
										",DELTA_RESERVE_C"
										",DELTA_STEM_C"
										",DELTA_LEAF_C"
										",DELTA_FROOT_C"
										",DELTA_CROOT_C"
										",DELTA_BRANCH_C"
										",DELTA_FRUIT_C"
										",RESERVE_N"
										",STEM_N"
										",STEMLIVE_N"
										",STEMDEAD_N"
										",LEAF_N"
										",FROOT_N"
										",CROOT_N"
										",CROOTLIVE_N"
										",CROOTDEAD_N"
										",BRANCH_N"
										",BRANCHLIVE_N"
										",BRANCHDEAD_N"
										",FRUIT_N"
										",DELTA_RESERVE_N"
										",DELTA_STEM_N"
										",DELTA_LEAF_N"
										",DELTA_FROOT_N"
										",DELTA_CROOT_N"
										",DELTA_BRANCH_N"
										",DELTA_FRUIT_N"
										",STEM_AR"
										",LEAF_AR"
										",FROOT_AR"
										",CROOT_AR"
										",BRANCH_AR"
										",F_CO2"
										",F_CO2_VER"
										",F_CO2_FRA"
										",FCO2_TR"
										",FLIGHT"
										",FAGE"
										",FT"
										",FT_ACCL"
										",FVPD"
										",FN"
										",FSW"
										",LITR_C"
										",CWD_C"
								);
								
		  /************************************************************************/
                 // if (!g_settings->regeneration)
                  if (!print_cell)
                  {  
		  /* heading variables only at cell level */
		  logger(g_daily_log,",gpp,npp,ar,hr,rsoil,reco,nee,nep,et,le,soil_evapo,snow_pack,asw,psi,moist_ratio,iWue,"
				"litrC,litr1C,litr2C,litr3C,litr4C,deadwoodC,deadwood2C,deadwood3C,deadwood4C,soilC,soil1C,soil2C,soil3C,soil4C,"
				"litrN,litr1N,litr2N,litr3N,litr4N,deadwoodN,deadwood2N,deadwood3N,deadwood4N,soilN,soil1N,soil2N,soil3N,soil4N,"
				"Tsoil,Daylength");

                  }
               	 /* end heading */
		 logger(g_daily_log,"\n");
													
					
	}
							
	/*****************************************************************************************************/

	/* values */
	//logger(g_daily_log, "%d,%d,%d,%d,%d", c->x, c->y, c->years[year].year, month + 1, day + 1);

	/* print class level values */
	for ( layer = c->tree_layers_count - 1; layer >= 0; --layer )
	{
               
		qsort(c->heights, c->heights_count, sizeof(height_t), sort_by_heights_desc);

		for ( height = 0; height < c->heights_count; ++height )
		{
			if( layer == c->heights[height].height_z )
			{
				

				for ( dbh = c->heights[height].dbhs_count - 1; dbh >= 0; --dbh )
				{
					
					for ( age = 0; age < c->heights[height].dbhs[dbh].ages_count ; ++age )
					{
						
						for ( species = 0; species < c->heights[height].dbhs[dbh].ages[age].species_count; ++species )
						{
							s  = &c->heights[height].dbhs[dbh].ages[age].species[species];

                                                        /* values */
           				        	logger(g_daily_log, "%d,%d,%d,%d,%d", c->x, c->y, c->years[year].year, month + 1, day + 1);

                                                        /* print layer */
							logger(g_daily_log,",%d", layer);

							/* print height */
							logger(g_daily_log,",%g", c->heights[height].value); 
							
							/* print dbh */
					                logger(g_daily_log,",%g", c->heights[height].dbhs[dbh].value);

                                                       /* print age */
						        logger(g_daily_log,",%d", c->heights[height].dbhs[dbh].ages[age].value);

							/* print species name */
							logger(g_daily_log,",%s", c->heights[height].dbhs[dbh].ages[age].species[species].name);

							/* print management */
							logger(g_daily_log,",%c", sz_management[c->heights[height].dbhs[dbh].ages[age].species[species].management]);

							/* print variables at layer-class level */

							logger(g_daily_log,",%6.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%d,"
									"%d,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f"
									",%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f"
									",%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f",
									s->value[GPP],
									s->value[Av_TOT],
									s->value[Aj_TOT],
									s->value[A_TOT],
									s->value[TOTAL_GROWTH_RESP],
									s->value[TOTAL_MAINT_RESP],
									s->value[TOTAL_AUT_RESP],
									s->value[NPP],
									s->value[BP],
									s->value[CUE],
									s->value[BPE],
									s->value[LAI_PROJ],
									s->value[PEAK_LAI_PROJ],
									s->value[LAI_EXP],
									s->value[DAILY_CANOPY_COVER_PROJ],
									s->value[DBHDC_EFF],
									s->value[CROWN_AREA_PROJ],
									s->value[PAR],
									s->value[APAR],
									s->value[fAPAR],
									s->counter[N_TREE],
									s->counter[VEG_DAYS],
									s->value[CANOPY_INT_RAIN],
									s->value[CANOPY_WATER],
									s->value[CANOPY_EVAPO],
									s->value[CANOPY_TRANSP],
									s->value[CANOPY_EVAPO_TRANSP],
									s->value[CANOPY_LATENT_HEAT],
									s->value[WUE],
									s->value[RESERVE_C],
									s->value[STEM_C],
									s->value[STEM_SAPWOOD_C],
									s->value[STEM_LIVEWOOD_C],
									s->value[STEM_DEADWOOD_C],
									s->value[LEAF_C],
									s->value[FROOT_C],
									s->value[CROOT_C],
									s->value[CROOT_SAPWOOD_C],
									s->value[CROOT_LIVEWOOD_C],
									s->value[CROOT_DEADWOOD_C],
									s->value[BRANCH_C],
									s->value[BRANCH_SAPWOOD_C],
									s->value[BRANCH_LIVEWOOD_C],
									s->value[BRANCH_DEADWOOD_C],
									s->value[FRUIT_C],
									s->value[TOT_SAPWOOD_C],
									s->value[C_TO_RESERVE],
									s->value[C_TO_STEM],
									s->value[C_TO_LEAF],
									s->value[C_TO_FROOT],
									s->value[C_TO_CROOT],
									s->value[C_TO_BRANCH],
									s->value[C_TO_FRUIT],
									s->value[RESERVE_N],
									s->value[STEM_N],
									s->value[STEM_LIVEWOOD_N],
									s->value[STEM_DEADWOOD_N],
									s->value[LEAF_N],
									s->value[FROOT_N],
									s->value[CROOT_N],
									s->value[CROOT_LIVEWOOD_N],
									s->value[CROOT_DEADWOOD_N],
									s->value[BRANCH_N],
									s->value[BRANCH_LIVEWOOD_N],
									s->value[BRANCH_DEADWOOD_N],
									s->value[FRUIT_N],
									s->value[N_TO_RESERVE],
									s->value[N_TO_STEM],
									s->value[N_TO_LEAF],
									s->value[N_TO_FROOT],
									s->value[N_TO_CROOT],
									s->value[N_TO_BRANCH],
									s->value[N_TO_FRUIT],
									s->value[STEM_AUT_RESP],
									s->value[LEAF_AUT_RESP],
									s->value[FROOT_AUT_RESP],
									s->value[CROOT_AUT_RESP],
									s->value[BRANCH_AUT_RESP],
									s->value[F_CO2],
									s->value[F_CO2_VER],
									s->value[F_CO2_FRANKS],
									s->value[F_CO2_TR],
									s->value[F_LIGHT_MAKELA],
									s->value[F_AGE],
									s->value[F_T],
									s->value[F_ACCL], 
									s->value[F_VPD],
									s->value[F_NUTR],
									s->value[F_SW],
									s->value[LITR_C],
									s->value[CWD_C]);
                                                                  
                                                                /************************************************************************/
        
          if (!print_cell)
          { 
	  /* printing variables only at cell level */

	  logger(g_daily_log, ",%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,"
			"%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f",
			c->daily_gpp,
			c->daily_npp,
			c->daily_aut_resp,
			c->daily_het_resp,
			c->daily_soil_resp,
			c->daily_r_eco,
			c->daily_nee,
			c->daily_nep,
			c->daily_et,
			c->daily_lh_flux,
			c->daily_soil_evapo,
			c->snow_pack,
			c->asw,
			c->psi,   
			c->soil_moist_ratio,
			c->daily_iwue,
			c->litrC,
			c->litr1C,
			c->litr2C,
			c->litr3C,
			c->litr4C,
			c->cwd_C,
			c->cwd_2C,
			c->cwd_3C,
			c->cwd_4C,
			c->soilC,
			c->soil1C,
			c->soil2C,
			c->soil3C,
			c->soil4C,
			c->litrN,
			c->litr1N,
			c->litr2N,
			c->litr3N,
			c->litr4N,
			c->cwd_N,
			c->cwd_2N,
			c->cwd_3N,
			c->cwd_4N,
			c->soilN,
			c->soil1N,
			c->soil2N,
			c->soil3N,
			c->soil4N,
			c->years[year].m[month].d[day].tsoil,
			c->years[year].m[month].d[day].daylength
	   );
           }
	   /************************************************************************/
           /* end heading */
	   logger(g_daily_log,"\n");   
                                                                  
                                                                  
                                                                  
						}

					}
					
				}
				
			}
			
		}
		
	  
        }   // loop on layer     
}
void EOD_cell_msg(void)
{
	if ( g_daily_log )
	{
		g_daily_log->std_output = 1;
		logger(g_daily_log, sz_launched, netcdf_get_version(), datetime_current());
		print_model_paths(g_daily_log);
		//const char* p;
		//p = file_get_name_only(g_daily_log->filename);
		logger(g_daily_log, "#output file = %s\n", g_daily_log->filename);
		print_model_settings(g_daily_log);
	}
}

void EOM_print_output_cell_level(cell_t *const c, const int month, const int year, const int years_of_simulation )
{
	int layer;
	int height;
	int dbh;
	int age;
	int species;

	static int print_header = 0;

	species_t *s;

	/* return if monthly logging is off*/
	if ( ! g_monthly_log ) return;

	/* heading */
	if ( ! print_header )
	{
		print_header = 1;

		logger(g_monthly_log, "X,Y,YEAR,MONTH");

		for ( layer = c->tree_layers_count - 1; layer >= 0; --layer )
		{
			for ( height = c->heights_count - 1; height >= 0 ; --height )
			{
				if( layer == c->heights[height].height_z )
				{
					/* heading for layers */
					logger(g_monthly_log,",LAYER");

					/* heading for heights value */
					logger(g_monthly_log,",HEIGHT");

					for ( dbh = c->heights[height].dbhs_count - 1; dbh >= 0; --dbh )
					{
						/* heading for dbhs value */
						logger(g_monthly_log, ",DBH");

						for ( age = 0; age < c->heights[height].dbhs[dbh].ages_count ; ++age )
						{
							/* heading for ages */
							logger(g_monthly_log,",AGE");

							for ( species = 0; species < c->heights[height].dbhs[dbh].ages[age].species_count; ++species )
							{
								/* heading for species name */
								logger(g_monthly_log,",SPECIES");

								/* heading for management */
								logger(g_monthly_log,",MANAGEMENT");

								logger(g_monthly_log,
										",GPP"
										",NET_ASS"
										",RA"
										",NPP"
										",CUE"
										",CTRANSP"
										",CET"
										",CLE"
										",LAI"
										",CC"
										",DBHDC"
										",HD_EFF"
										",HDMAX"
										",HDMIN"
										",N_TREE"
										",WUE"
										",WRes"
										",WS"
										",WSL"
										",WSD"
										",PWL"
										",PWFR"
										",WCR"
										",WCRL"
										",WCRD"
										",WBB"
										",WBBL"
										",WBBD"
								);

							}
							if ( c->heights[height].dbhs[dbh].ages[age].species_count > 1 ) {
								logger(g_monthly_log,",*");
							}
						}
						if ( c->heights[height].dbhs[dbh].ages_count > 1 ) {
							logger(g_monthly_log,",**");
						}
					}
					if ( c->heights[height].dbhs_count > 1 ) {
						logger(g_monthly_log,",***");
					}
				}
				if ( c->tree_layers[layer].layer_n_height_class > 1 ) {
					logger(g_monthly_log,",****");
				}
			}
			if ( c->tree_layers_count > 1 ) {
				logger(g_monthly_log,",*****");
			}
		}
		/************************************************************************/
                if (!g_settings->regeneration)
                { 
		/* heading variables at cell level only if there's more than one layer */

		//if( c->heights_count > 1 )
		//{
			logger(g_monthly_log,",***,gpp,npp,ar");
		//}
		/* heading variables at cell level also if there's more than one layer */
		//else
		//{
		//	logger(g_monthly_log,",*****");
		//}
		/* heading variables only at cell level */
		logger(g_monthly_log,",et,le,asw,iWue");
                }
                /* end heading */  
                logger(g_monthly_log,"\n");
                
	}

	/************************************************************************/
	/* values */
	logger ( g_monthly_log, "%d,%d,%d,%d", c->x, c->y, c->years[year].year, month +1 );

	/* print class level values */
	for ( layer = c->tree_layers_count - 1; layer >= 0; --layer )
	{
		for ( height = c->heights_count - 1; height >= 0 ; --height )
		{
			if( layer == c->heights[height].height_z )
			{
				/* print layer */
				logger(g_monthly_log,",%d", layer);

				/* print height */
				logger(g_monthly_log,",%g", c->heights[height].value);

				for ( dbh = c->heights[height].dbhs_count - 1; dbh >= 0; --dbh )
				{
					/* print age */
					logger(g_monthly_log,",%g", c->heights[height].dbhs[dbh].value);

					for ( age = 0; age < c->heights[height].dbhs[dbh].ages_count ; ++age )
					{
						/* print age */
						logger(g_monthly_log,",%d", c->heights[height].dbhs[dbh].ages[age].value);

						for ( species = 0; species < c->heights[height].dbhs[dbh].ages[age].species_count; ++species )
						{
							s  = &c->heights[height].dbhs[dbh].ages[age].species[species];

							/* print species name */
							logger(g_monthly_log,",%s", c->heights[height].dbhs[dbh].ages[age].species[species].name);

							/* print management */
							logger(g_monthly_log,",%c", sz_management[c->heights[height].dbhs[dbh].ages[age].species[species].management]);

							/* print variables at layer-class level */
							logger(g_monthly_log,",%6.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%d,%3.4f,%3.4f,%3.4f,%3.4f"
									",%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f",
									s->value[MONTHLY_GPP],
									s->value[MONTHLY_ASSIMILATION],
									s->value[MONTHLY_TOTAL_AUT_RESP],
									s->value[MONTHLY_NPP],
									s->value[MONTHLY_CUE],
									s->value[MONTHLY_CANOPY_TRANSP],
									s->value[MONTHLY_CANOPY_EVAPO_TRANSP],
									s->value[MONTHLY_CANOPY_LATENT_HEAT],
									s->value[MONTHLY_LAI_PROJ]/31.,   //5p6 
									s->value[CANOPY_COVER_PROJ],
									s->value[DBHDC_EFF],
									s->value[HD_EFF],
									s->value[HD_MAX],
									s->value[HD_MIN],
									s->counter[N_TREE],
									s->value[MONTHLY_WUE],
									s->value[RESERVE_C],
									s->value[STEM_C],
									s->value[STEM_LIVEWOOD_C],
									s->value[STEM_DEADWOOD_C],
									s->value[MAX_LEAF_C],
									s->value[MAX_FROOT_C],
									s->value[CROOT_C],
									s->value[CROOT_LIVEWOOD_C],
									s->value[CROOT_DEADWOOD_C],
									s->value[BRANCH_C],
									s->value[BRANCH_LIVEWOOD_C],
									s->value[BRANCH_DEADWOOD_C]);
						}
						if ( c->heights[height].dbhs[dbh].ages[age].species_count > 1 ) {
							logger(g_monthly_log,",*");
						}
					}
					if ( c->heights[height].dbhs[dbh].ages_count > 1 ) {
						logger(g_monthly_log,",**");
					}
				}
				if ( c->heights[height].dbhs_count > 1 ) {
					logger(g_monthly_log,",***");
				}
			}
			if ( c->tree_layers[layer].layer_n_height_class > 1 ) {
				logger(g_monthly_log,",****");
			}
		}
		if ( c->tree_layers_count > 1 ) {
			logger(g_monthly_log,",*****");
		}
	}
	/************************************************************************/
         if (!g_settings->regeneration)
                 { 
	/* printing variables at cell level only if there's more than one layer */

	//if( c->heights_count > 1 )
	//{
		logger(g_monthly_log, ",***,%3.4f,%3.4f,%3.4f",
				c->monthly_gpp,
				c->monthly_npp,
				c->monthly_aut_resp);

	//}
	/* printing variables at cell level also if there's more than one layer */
	//else
	//{
	//	logger(g_monthly_log,",*****");
	//}
	/* printing variables only at cell level */
	logger(g_monthly_log, ",%3.2f,%3.2f,%3.2f,%3.2f",
			c->monthly_et,
			c->monthly_lh_flux,
			c->asw,
			c->monthly_iwue);
        }
        /* end print */
	logger(g_monthly_log,"\n");
}


void EOM_print_output_cell_level_ddalmo(cell_t *const c, const int month, const int year, const int years_of_simulation )
{
	int layer;
	int height;
	int dbh;
	int age;
	int species;

	static int print_header = 0;

	species_t *s;

	/* return if monthly logging is off*/
	if ( ! g_monthly_log ) return;

	/* heading */
	if ( ! print_header )
	{
		print_header = 1;

		//logger(g_monthly_log, "X,Y,YEAR,MONTH");

		for ( layer = c->tree_layers_count - 1; layer >= 0; --layer )
		{
                        logger(g_monthly_log, "X,Y,YEAR,MONTH");

			for ( height = c->heights_count - 1; height >= 0 ; --height )
			{
				if( layer == c->heights[height].height_z )
				{
					/* heading for layers */
					logger(g_monthly_log,",LAYER");

					/* heading for heights value */
					logger(g_monthly_log,",HEIGHT");

					for ( dbh = c->heights[height].dbhs_count - 1; dbh >= 0; --dbh )
					{
						/* heading for dbhs value */
						logger(g_monthly_log, ",DBH");

						for ( age = 0; age < c->heights[height].dbhs[dbh].ages_count ; ++age )
						{
							/* heading for ages */
							logger(g_monthly_log,",AGE");

							for ( species = 0; species < c->heights[height].dbhs[dbh].ages[age].species_count; ++species )
							{
								/* heading for species name */
								logger(g_monthly_log,",SPECIES");

								/* heading for management */
								logger(g_monthly_log,",MANAGEMENT");

								logger(g_monthly_log,
										",GPP"
										",NET_ASS"
										",RA"
										",NPP"
										",CUE"
										",CTRANSP"
										",CET"
										",CLE"
										",LAI"
										",CC"
										",DBHDC"
										",HD_EFF"
										",HDMAX"
										",HDMIN"
										",N_TREE"
										",WUE"
										",WRes"
										",WS"
										",WSL"
										",WSD"
										",PWL"
										",PWFR"
										",WCR"
										",WCRL"
										",WCRD"
										",WBB"
										",WBBL"
										",WBBD"
								);

							}
							if ( c->heights[height].dbhs[dbh].ages[age].species_count > 1 ) {
								logger(g_monthly_log,",*");
							}
						}
						if ( c->heights[height].dbhs[dbh].ages_count > 1 ) {
							logger(g_monthly_log,",**");
						}
					}
					if ( c->heights[height].dbhs_count > 1 ) {
						logger(g_monthly_log,",***");
					}
				}
				if ( c->tree_layers[layer].layer_n_height_class > 1 ) {
					logger(g_monthly_log,",****");
				}
			}
			//if ( c->tree_layers_count > 1 ) {
			//	logger(g_monthly_log,",*****");
			//}
		//}
		/************************************************************************/
                if (!g_settings->regeneration)
                { 
		/* heading variables at cell level only if there's more than one layer */

		//if( c->heights_count > 1 )
		//{
			logger(g_monthly_log,",***,gpp,npp,ar");
		//}
		/* heading variables at cell level also if there's more than one layer */
		//else
		//{
		//	logger(g_monthly_log,",*****");
		//}
		/* heading variables only at cell level */
		logger(g_monthly_log,",et,le,asw,iWue");
                }
                /* end heading */  
                logger(g_monthly_log,"\n");
               } //loop layer
                
	}

	/************************************************************************/
	/* values */
	//logger ( g_monthly_log, "%d,%d,%d,%d", c->x, c->y, c->years[year].year, month +1 );

	/* print class level values */
	for ( layer = c->tree_layers_count - 1; layer >= 0; --layer )
	{
                /* values */
         	logger ( g_monthly_log, "%d,%d,%d,%d", c->x, c->y, c->years[year].year, month +1 );
 
		for ( height = c->heights_count - 1; height >= 0 ; --height )
		{
			if( layer == c->heights[height].height_z )
			{
				/* print layer */
				logger(g_monthly_log,",%d", layer);

				/* print height */
				logger(g_monthly_log,",%g", c->heights[height].value);

				for ( dbh = c->heights[height].dbhs_count - 1; dbh >= 0; --dbh )
				{
					/* print age */
					logger(g_monthly_log,",%g", c->heights[height].dbhs[dbh].value);

					for ( age = 0; age < c->heights[height].dbhs[dbh].ages_count ; ++age )
					{
						/* print age */
						logger(g_monthly_log,",%d", c->heights[height].dbhs[dbh].ages[age].value);

						for ( species = 0; species < c->heights[height].dbhs[dbh].ages[age].species_count; ++species )
						{
							s  = &c->heights[height].dbhs[dbh].ages[age].species[species];

							/* print species name */
							logger(g_monthly_log,",%s", c->heights[height].dbhs[dbh].ages[age].species[species].name);

							/* print management */
							logger(g_monthly_log,",%c", sz_management[c->heights[height].dbhs[dbh].ages[age].species[species].management]);

							/* print variables at layer-class level */
							logger(g_monthly_log,",%6.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%d,%3.4f,%3.4f,%3.4f,%3.4f"
									",%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f",
									s->value[MONTHLY_GPP],
									s->value[MONTHLY_ASSIMILATION],
									s->value[MONTHLY_TOTAL_AUT_RESP],
									s->value[MONTHLY_NPP],
									s->value[MONTHLY_CUE],
									s->value[MONTHLY_CANOPY_TRANSP],
									s->value[MONTHLY_CANOPY_EVAPO_TRANSP],
									s->value[MONTHLY_CANOPY_LATENT_HEAT],
								       s->value[MONTHLY_LAI_PROJ]/31.,   //5p6 divided by 31, otherwise we have value higher than the PEAK_LAI
									s->value[CANOPY_COVER_PROJ],
									s->value[DBHDC_EFF],
									s->value[HD_EFF],
									s->value[HD_MAX],
									s->value[HD_MIN],
									s->counter[N_TREE],
									s->value[MONTHLY_WUE],
									s->value[RESERVE_C],
									s->value[STEM_C],
									s->value[STEM_LIVEWOOD_C],
									s->value[STEM_DEADWOOD_C],
									s->value[MAX_LEAF_C],
									s->value[MAX_FROOT_C],
									s->value[CROOT_C],
									s->value[CROOT_LIVEWOOD_C],
									s->value[CROOT_DEADWOOD_C],
									s->value[BRANCH_C],
									s->value[BRANCH_LIVEWOOD_C],
									s->value[BRANCH_DEADWOOD_C]);
						}
						if ( c->heights[height].dbhs[dbh].ages[age].species_count > 1 ) {
							logger(g_monthly_log,",*");
						}
					}
					if ( c->heights[height].dbhs[dbh].ages_count > 1 ) {
						logger(g_monthly_log,",**");
					}
				}
				if ( c->heights[height].dbhs_count > 1 ) {
					logger(g_monthly_log,",***");
				}
			}
			if ( c->tree_layers[layer].layer_n_height_class > 1 ) {
				logger(g_monthly_log,",****");
			}
		}
		//if ( c->tree_layers_count > 1 ) {
		//	logger(g_monthly_log,",*****");
		//}
	//}
	/************************************************************************/
         if (!g_settings->regeneration)
                 { 
	/* printing variables at cell level only if there's more than one layer */

	//if( c->heights_count > 1 )
	//{
		logger(g_monthly_log, ",***,%3.4f,%3.4f,%3.4f",
				c->monthly_gpp,
				c->monthly_npp,
				c->monthly_aut_resp);

	//}
	/* printing variables at cell level also if there's more than one layer */
	//else
	//{
	//	logger(g_monthly_log,",*****");
	//}
	/* printing variables only at cell level */
	logger(g_monthly_log, ",%3.2f,%3.2f,%3.2f,%3.2f",
			c->monthly_et,
			c->monthly_lh_flux,
			c->asw,
			c->monthly_iwue);
        }
        /* end print */
	logger(g_monthly_log,"\n");
       } //loop layer
}

void EOM_print_output_cell_level_mc(cell_t *const c, const int month, const int year, const int years_of_simulation )
{
	int layer;
	int height;
	int dbh;
	int age;
	int species;

	static int print_header = 0;
        int print_cell = 0;   // print cell level value   

	species_t *s;

	/* return if monthly logging is off*/
	if ( ! g_monthly_log ) return;

	/* heading */
	if ( ! print_header )
	{
		print_header = 1;

		logger(g_monthly_log, "X,Y,YEAR,MONTH");
		/* heading for layers */
								logger(g_monthly_log,",LAYER");

								/* heading for heights value */
								logger(g_monthly_log,",HEIGHT");
								
								/* heading for dbhs value */
								logger(g_monthly_log, ",DBH");
					
							       /* heading for ages */
							        logger(g_monthly_log,",AGE");
								/* heading for species name */
								logger(g_monthly_log,",SPECIES");

								/* heading for management */
								logger(g_monthly_log,",MANAGEMENT");

								logger(g_monthly_log,
										",GPP"
										",NET_ASS"
										",RA"
										",NPP"
										",CUE"
										",CTRANSP"
										",CET"
										",CLE"
										",LAI"
										",CC"
										",DBHDC"
										",HD_EFF"
										",HDMAX"
										",HDMIN"
										",N_TREE"
										",WUE"
										",WRes"
										",WS"
										",WSL"
										",WSD"
										",PWL"
										",PWFR"
										",WCR"
										",WCRL"
										",WCRD"
										",WBB"
										",WBBL"
										",WBBD"
								);
								if (!print_cell)
                						{ 

							        logger(g_monthly_log,",gpp,npp,ar");
									
								logger(g_monthly_log,",et,le,asw,iWue");
               						 }
                						/* end heading */  
                						logger(g_monthly_log,"\n");

       }
							

	/************************************************************************/
	/* values */
	//logger ( g_monthly_log, "%d,%d,%d,%d", c->x, c->y, c->years[year].year, month +1 );
	
	
	//if ( c->heights_count )    //FIXME only used if we want to print info at cell level when there are more than one class only
	//{
	
	/* print class level values */
        for ( layer = c->tree_layers_count - 1; layer >= 0; --layer )
	{
                 
        // qsort(c->heights, c->heights_count, sizeof(height_t), sort_by_heights_desc);  // lo mettiamo?

		for ( height = c->heights_count - 1; height >= 0 ; --height )
		{
			if( layer == c->heights[height].height_z )
			{
				
				for ( dbh = c->heights[height].dbhs_count - 1; dbh >= 0; --dbh )
				{
				
					for ( age = 0; age < c->heights[height].dbhs[dbh].ages_count ; ++age )
					{
						
						for ( species = 0; species < c->heights[height].dbhs[dbh].ages[age].species_count; ++species )
						{
							s  = &c->heights[height].dbhs[dbh].ages[age].species[species];
							
												
							  /* values */
         	                                         logger ( g_monthly_log, "%d,%d,%d,%d", c->x, c->y, c->years[year].year, month +1 );
         	
         	                                         /* print layer */
				                         logger(g_monthly_log,",%d", layer);

				                         /* print height */
				                         logger(g_monthly_log,",%g", c->heights[height].value);
				
					                 /* print dbh */
					                 logger(g_monthly_log,",%g", c->heights[height].dbhs[dbh].value);
					
					                 /* print age */
						         logger(g_monthly_log,",%d", c->heights[height].dbhs[dbh].ages[age].value);

							/* print species name */
							logger(g_monthly_log,",%s", c->heights[height].dbhs[dbh].ages[age].species[species].name);

							/* print management */
							logger(g_monthly_log,",%c", sz_management[c->heights[height].dbhs[dbh].ages[age].species[species].management]);

							/* print variables at layer-class level */
							logger(g_monthly_log,",%6.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%d,%3.4f,%3.4f,%3.4f,%3.4f"
									",%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f",
									s->value[MONTHLY_GPP],
									s->value[MONTHLY_ASSIMILATION],
									s->value[MONTHLY_TOTAL_AUT_RESP],
									s->value[MONTHLY_NPP],
									s->value[MONTHLY_CUE],
									s->value[MONTHLY_CANOPY_TRANSP],
									s->value[MONTHLY_CANOPY_EVAPO_TRANSP],
									s->value[MONTHLY_CANOPY_LATENT_HEAT],
									s->value[MONTHLY_LAI_PROJ]/31.,   //divided by 31, otherwise we have value higher than the PEAK_LAI
									s->value[CANOPY_COVER_PROJ],
									s->value[DBHDC_EFF],
									s->value[HD_EFF],
									s->value[HD_MAX],
									s->value[HD_MIN],
									s->counter[N_TREE],
									s->value[MONTHLY_WUE],
									s->value[RESERVE_C],
									s->value[STEM_C],
									s->value[STEM_LIVEWOOD_C],
									s->value[STEM_DEADWOOD_C],
									s->value[MAX_LEAF_C],
									s->value[MAX_FROOT_C],
									s->value[CROOT_C],
									s->value[CROOT_LIVEWOOD_C],
									s->value[CROOT_DEADWOOD_C],
									s->value[BRANCH_C],
									s->value[BRANCH_LIVEWOOD_C],
									s->value[BRANCH_DEADWOOD_C]);
									/************************************************************************/

         if (!print_cell)
                 { 
               /* printing variables only at cell level */
		logger(g_monthly_log, ",%3.4f,%3.4f,%3.4f",
				c->monthly_gpp,
				c->monthly_npp,
				c->monthly_aut_resp);

	
	
	logger(g_monthly_log, ",%3.2f,%3.2f,%3.2f,%3.2f",
			c->monthly_et,
			c->monthly_lh_flux,
			c->asw,
			c->monthly_iwue);
        }
        /* end print */
	logger(g_monthly_log,"\n");	
									
									
						}
						
					}
					
				}
				
			}
			
		}
		

       } //loop layer
       
   // } // check on height_count   
}
void EOM_cell_msg(void)
{
	if ( g_monthly_log )
	{
		g_monthly_log->std_output = 1;
		logger(g_monthly_log, sz_launched, netcdf_get_version(), datetime_current());
		print_model_paths(g_monthly_log);
		//const char* p;
		//p = file_get_name_only(g_monthly_log->filename);
		logger(g_monthly_log, "#output file = %s\n", g_monthly_log->filename);
		print_model_settings(g_monthly_log);
	}
}
#if 0
void EOY_print_output_cell_level(cell_t *const c, const int year, const int years_of_simulation )
{
	int layer;
	int height;
	int dbh;
	int age;
	int species;

	static int print_header = 0;

	species_t *s;

	/* return if annual logging is off*/
	if ( ! g_annual_log ) return;

	/* test */
	//note it can be used only if no other classes are added!!
	/*
	if ( !year)
	{
		c->initial_tree_layers_count = c->tree_layers_count;
		c->initial_heights_count = c->heights_count;

		for ( height = c->heights_count - 1; height >= 0 ; --height )
		{
			c->heights[height].initial_dbhs_count = c->heights[height].dbhs_count;
			for ( dbh = c->heights[height].dbhs_count - 1; dbh >= 0; --dbh )
			{
				c->heights[height].dbhs[dbh].initial_ages_count = c->heights[height].dbhs[dbh].ages_count;

				for ( age = 0; age < c->heights[height].dbhs[dbh].ages_count ; ++age )
				{
					c->heights[height].dbhs[dbh].ages[age].initial_species_count = c->heights[height].dbhs[dbh].ages[age].species_count;
				}
			}
		}
	}
	 */

	/* heading */
	if ( ! print_header )
	{
		print_header = 1;

		logger(g_annual_log, "X,Y,YEAR");

		for ( layer = c->tree_layers_count - 1; layer >= 0; --layer )

		{
			for ( height = c->heights_count - 1; height >= 0 ; --height )
			{
				if( layer == c->heights[height].height_z )
				{
					/* heading for layers */
					logger(g_annual_log,",LAYER");

					/* heading for heights value */
					logger(g_annual_log,",HEIGHT");

					for ( dbh = c->heights[height].dbhs_count - 1; dbh >= 0; --dbh )
					{
						/* heading for dbhs value */
						logger(g_annual_log, ",DBH");

						for ( age = 0; age < c->heights[height].dbhs[dbh].ages_count ; ++age )
						{
							/* heading for ages */
							logger(g_annual_log,",AGE");

							 for ( species = 0; species < c->heights[height].dbhs[dbh].ages[age].species_count; ++species )
                           
							{
								/* heading for species name */
								logger(g_annual_log,",SPECIES");
                                                             
								/* heading for management name */
								logger(g_annual_log,",MANAGEMENT");

								logger(g_annual_log,
										",GPP"
										",GPP_SUN:GPP"
										",GPP_SHADE:GPP"
										",Av_SUN:A_SUN"
										",Aj_SUN:A_SUN"
										",Av_SHADE:A_SHADE"
										",Aj_SHADE:A_SHADE"
										",Av_TOT:A_TOT"
										",Aj_TOT:A_TOT"
										",GR"
										",MR"
										",RA"
										",NPP"
										",BP"
										",reser_as_diff"
										",ResAlloc"
										",ResDeple"
										",ResUsage"
										",BP/NPP"
										",ResAlloc/NPP"
										",ResAlloc/BP"
										",ResDeple/NPP"
										",ResDeple/BP"
										",ResUsage/NPP"
										",ResUsage/BP"
										",CUE"
										",BPE"
										",diffCUE-BPE"
										",Y(perc)"
										",MAX_NSC_CONC"
										",MIN_NSC_CONC"
										",PeakLAI"
										",MaxLAI"
										",SLA"
										",SAPWOOD-AREA"
										",HEARTWOOD-AREA"
										",CC-Proj"
										",DBHDC"
										",CROWN_DIAMETER"
										",CROWN_HEIGHT"
										",CROWN_AREA_PROJ"
										",APAR"
										",LIVETREE"
										",DEADTREE"
										",THINNEDTREE"
										",VEG_D"
										",FIRST_VEG_DAY"
										",CTRANSP"
										",CINT"
										",CLE"
										",WUE"
										",MAX_ANN_RESERVE_C"
                                                                                ",MIN_ANN_RESERVE_C"
										",MIN_RESERVE_C"
										",RESERVE_C"
										",STEM_C"
										",STEMSAP_C"
										",STEMHEART_C"
										",STEMSAP_PERC"
										",STEMLIVE_C"
										",STEMDEAD_C"
										",STEMLIVE_PERC"
										",MAX_LEAF_C"
										",MAX_FROOT_C"
										",CROOT_C"
										",CROOTLIVE_C"
										",CROOTDEAD_C"
										",CROOTLIVE_PERC"
										",BRANCH_C"
										",BRANCHLIVE_C"
										",BRANCHDEAD_C"
										",BRANCHLIVE_PERC"
										",FRUIT_C"
										",MAX_FRUIT_C"
										",TOT_SAPWOOD_C"
										",TOT_HEARTWOOD_C"
										",RESERVE_N"
										",STEM_N"
										",STEMLIVE_N"
										",STEMDEAD_N"
										",CROOT_N"
										",CROOTLIVE_N"
										",CROOTDEAD_N"
										",BRANCH_N"
										",BRANCHLIVE_N"
										",BRANCHDEAD_N"
										",FRUIT_N"
										",STANDING_WOOD"
										",DELTA_WOOD"
										",CUM_DELTA_WOOD"
										",BASAL_AREA"
										",TREE_CAI"
										",TREE_MAI"
										",CAI"
										",MAI"
										",VOLUME"
										",TREE_VOLUME"
										",DELTA_TREE_VOL(perc)"
										",DELTA_AGB"
										",DELTA_BGB"
										",AGB"
										",BGB"
										",BGB:AGB"
										",DELTA_TREE_AGB"
										",DELTA_TREE_BGB"
										",C_HWP"
										",VOLUME_HWP"
										",STEM_RA"
										",LEAF_RA"
										",FROOT_RA"
										",CROOT_RA"
										",BRANCH_RA");
							}
							if ( c->heights[height].dbhs[dbh].ages[age].species_count > 1 ) {
								logger(g_annual_log,",*");
							}
						}
						if ( c->heights[height].dbhs[dbh].ages_count > 1 ) {
							logger(g_annual_log,",**");
						}
					}
					if ( c->heights[height].dbhs_count > 1 ) {
						logger(g_annual_log,",***");
					}
				}
				if ( c->tree_layers[layer].layer_n_height_class > 1 ) {
					logger(g_annual_log,",****");
				}
			}
			if ( c->tree_layers_count > 1 ) {
				logger(g_annual_log,",*****");
			}
		}
		/************************************************************************/
		/* heading cell variables */
		logger(g_annual_log,",gpp,npp,ar,hr,rsoil,rsoilCO2,reco,nee,nep,et,le,soil-evapo,avg_asw,iWue,vol,cum_vol,run_off,"   // put avg_asw
				"litrC,litr1C,litr2C,litr3C,litr4C,cwdC,cwd2C,cwd3C,cwd4C,soilC,soil1C,soil2C,soil3C,soil4C,"
				"litrN,litr1N,litr2N,litr3N,litr4N,cwdN,cwd2N,cwd3N,cwd4N,soilN,soil1N,soil2N,soil3N,soil4N,soil_depth"); 
		/************************************************************************/
		/* heading meteo variables */
		logger(g_annual_log,",solar_rad,tavg,tmax,tmin,tday,tnight,vpd,prcp,tsoil,rh,[CO2]"); // removed avg_asw
		/************************************************************************/
		/* end heading */
		logger(g_annual_log,"\n");
		/************************************************************************/

	}

	/*****************************************************************************************************/


	/* values */
	logger(g_annual_log, "%d,%d,%d", c->x, c->y, c->years[year].year);

	/* print class level values */
	if ( c->heights_count )
	{
		for ( layer = c->tree_layers_count - 1; layer >= 0; --layer )
		{
			qsort(c->heights, c->heights_count, sizeof(height_t), sort_by_heights_desc);

			for ( height = 0; height < c->heights_count; ++height )
			{
				if( layer == c->heights[height].height_z )
				{
					/* print layer */
					logger(g_annual_log,",%d", layer);

					/* print height */
					logger(g_annual_log,",%g", c->heights[height].value);

					for ( dbh = c->heights[height].dbhs_count - 1; dbh >= 0; --dbh )
					{
						/* print dbh */
						logger(g_annual_log,",%g", c->heights[height].dbhs[dbh].value);

						// ALESSIOR TO ALESSIOC
						// SHOULD THIS BE CHANGED TO
						// start from c->heights[height].dbhs[dbh].ages_count-1 ?
						for ( age = 0; age < c->heights[height].dbhs[dbh].ages_count ; ++age )
						{
							/* print age */
							logger(g_annual_log,",%d", c->heights[height].dbhs[dbh].ages[age].value);

							for ( species = 0; species < c->heights[height].dbhs[dbh].ages[age].species_count; ++species )
							{
								s  = &c->heights[height].dbhs[dbh].ages[age].species[species];

								/* print species name */
								logger(g_annual_log,",%s", c->heights[height].dbhs[dbh].ages[age].species[species].name);

								/* print management */
								logger(g_annual_log,",%c", sz_management[c->heights[height].dbhs[dbh].ages[age].species[species].management]);

								/* print variables at layer-class level */
								logger(g_annual_log,",%6.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%d,%d,%d,%d,%d,%3.4f"
										",%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f"
										",%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f"
										",%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f",
										s->value[YEARLY_GPP],
										s->value[YEARLY_GPP_SUN]   / s->value[YEARLY_GPP],
										s->value[YEARLY_GPP_SHADE] / s->value[YEARLY_GPP],
										s->value[YEARLY_Av_SUN]    / s->value[YEARLY_A_SUN],
										s->value[YEARLY_Aj_SUN]    / s->value[YEARLY_A_SUN],
										s->value[YEARLY_Av_SHADE]  / s->value[YEARLY_A_SHADE],
										s->value[YEARLY_Aj_SHADE]  / s->value[YEARLY_A_SHADE],
										s->value[YEARLY_Av_TOT]    / s->value[YEARLY_A_TOT],
										s->value[YEARLY_Aj_TOT]    / s->value[YEARLY_A_TOT],
										s->value[YEARLY_TOTAL_GROWTH_RESP],
										s->value[YEARLY_TOTAL_MAINT_RESP],
										s->value[YEARLY_TOTAL_AUT_RESP],
										s->value[YEARLY_NPP],
										s->value[YEARLY_BP],
										(s->value[YEARLY_NPP] - s->value[YEARLY_BP]),
										s->value[YEARLY_RESERVE_ALLOC],
										s->value[YEARLY_RESERVE_DEPLE],
										s->value[YEARLY_RESERVE_USAGE],
										(s->value[YEARLY_BP] / s->value[YEARLY_NPP]),
										(s->value[YEARLY_RESERVE_ALLOC] / s->value[YEARLY_NPP]),
										(s->value[YEARLY_RESERVE_ALLOC] / s->value[YEARLY_BP]),
										(s->value[YEARLY_RESERVE_DEPLE] / s->value[YEARLY_NPP]),
										(s->value[YEARLY_RESERVE_DEPLE] / s->value[YEARLY_BP]),
										(s->value[YEARLY_RESERVE_USAGE] / s->value[YEARLY_NPP]),
										(s->value[YEARLY_RESERVE_USAGE] / s->value[YEARLY_BP]),
										s->value[YEARLY_CUE],
										s->value[YEARLY_BPE],
										(s->value[YEARLY_CUE] - s->value[YEARLY_BPE]),
										s->value[YEARLY_TOTAL_AUT_RESP] / s->value[YEARLY_GPP] * 100.,
										s->value[MAX_RESERVE_C_CONC],
										s->value[MIN_RESERVE_C_CONC],
										s->value[PEAK_LAI_PROJ],
										s->value[MAX_LAI_PROJ],
										s->value[SLA_PROJ],
										s->value[SAPWOOD_AREA],
										s->value[HEARTWOOD_AREA],
										s->value[CANOPY_COVER_PROJ],
										s->value[DBHDC_EFF],
										s->value[CROWN_DIAMETER],
										s->value[CROWN_HEIGHT],
										s->value[CROWN_AREA_PROJ],
										s->value[YEARLY_APAR],
										s->counter[N_TREE],
										s->counter[DEAD_TREE],
										s->counter[THINNED_TREE],
										s->counter[YEARLY_VEG_DAYS],
										s->counter[FIRST_VEG_DAYS],
										s->value[YEARLY_CANOPY_TRANSP],
										s->value[YEARLY_CANOPY_INT],
										s->value[YEARLY_CANOPY_LATENT_HEAT],
										s->value[YEARLY_WUE],
										s->value[MAX_ANN_RESERVE_C],
                                                                               s->value[MIN_ANN_RESERVE_C],  //5p6
										s->value[MIN_RESERVE_C],
										s->value[RESERVE_C],
										s->value[STEM_C],
										s->value[STEM_SAPWOOD_C],
										s->value[STEM_HEARTWOOD_C],
										(s->value[STEM_SAPWOOD_C] * 100. ) / s->value[STEM_C],
										s->value[STEM_LIVEWOOD_C],
										s->value[STEM_DEADWOOD_C],
										(s->value[STEM_LIVEWOOD_C] * 100. ) / s->value[STEM_C],
										s->value[MAX_LEAF_C],
										s->value[MAX_FROOT_C],
										s->value[CROOT_C],
										s->value[CROOT_LIVEWOOD_C],
										s->value[CROOT_DEADWOOD_C],
										(s->value[CROOT_LIVEWOOD_C] * 100. ) / s->value[CROOT_C],
										s->value[BRANCH_C],
										s->value[BRANCH_LIVEWOOD_C],
										s->value[BRANCH_DEADWOOD_C],
										(s->value[BRANCH_LIVEWOOD_C] * 100. ) / s->value[BRANCH_C],
										s->value[FRUIT_C],
										s->value[MAX_FRUIT_C],
										s->value[TOT_SAPWOOD_C],
										s->value[TOT_HEARTWOOD_C],
										s->value[RESERVE_N],
										s->value[STEM_N],
										s->value[STEM_LIVEWOOD_N],
										s->value[STEM_DEADWOOD_N],
										s->value[CROOT_N],
										s->value[CROOT_LIVEWOOD_N],
										s->value[CROOT_DEADWOOD_N],
										s->value[BRANCH_N],
										s->value[BRANCH_LIVEWOOD_N],
										s->value[BRANCH_DEADWOOD_N],
										s->value[FRUIT_N],
										s->value[STANDING_WOOD],
										s->value[YEARLY_C_TO_WOOD],
										s->value[CUM_YEARLY_C_TO_WOOD],
										s->value[STAND_BASAL_AREA_m2],
										s->value[TREE_CAI],
										s->value[TREE_MAI],
										s->value[CAI],
										s->value[MAI],
										s->value[VOLUME],
										s->value[TREE_VOLUME],
										(s->value[TREE_CAI]*100.)/s->value[TREE_VOLUME],
										s->value[DELTA_AGB],
										s->value[DELTA_BGB],
										s->value[AGB],
										s->value[BGB],
										s->value[BGB]/s->value[AGB],
										s->value[DELTA_TREE_AGB],
										s->value[DELTA_TREE_BGB],
										s->value[C_HWP],
										s->value[VOLUME_HWP],
										s->value[YEARLY_STEM_AUT_RESP],
										s->value[YEARLY_LEAF_AUT_RESP],
										s->value[YEARLY_FROOT_AUT_RESP],
										s->value[YEARLY_CROOT_AUT_RESP],
										s->value[YEARLY_BRANCH_AUT_RESP]);
							}
							if ( c->heights[height].dbhs[dbh].ages[age].species_count > 1 ) {
								logger(g_annual_log,",*");
							}
						}
						if ( c->heights[height].dbhs[dbh].ages_count > 1 ) {
							logger(g_annual_log,",**");
						}
					}
					if ( c->heights[height].dbhs_count > 1 ) {
						logger(g_annual_log,",***");
					}
				}
				if ( c->tree_layers[layer].layer_n_height_class > 1 ) {
					logger(g_annual_log,",****");
				}
			}
			if ( c->tree_layers_count > 1 )
			{
				logger(g_annual_log,",*****");
			}
		}
	}
	else
	{
		//ALESSIOC TO ALLESSIOR PRINT EMPTY SPACES WHEN N_TREE = 0
	}
	/************************************************************************/
	/* printing variables at cell level only if there's more than one layer */
	logger(g_annual_log, ",%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f"
			",%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f",
			c->annual_gpp,
			c->annual_npp,
			c->annual_aut_resp,
			c->annual_het_resp,
			c->annual_soil_resp,
			c->annual_soil_respCO2,
			c->annual_r_eco,
			c->annual_nee,
			c->annual_nep,
			c->annual_et,
			c->annual_lh_flux,
			c->annual_soil_evapo,
                        c->years[year].yearly_mean.asw,
			//c->asw,   //  comment: this is actually the value at the 31 december, it has more sense to have the mean value
			c->annual_iwue,
			c->volume,
			c->cum_volume,
			c->annual_out_flow,
			c->litrC,
			c->litr1C,
			c->litr2C,
			c->litr3C,
			c->litr4C,
			c->cwd_C,
			c->cwd_2C,
			c->cwd_3C,
			c->cwd_4C,
			c->soilC,
			c->soil1C,
			c->soil2C,
			c->soil3C,
			c->soil4C,
			c->litrN,
			c->litr1N,
			c->litr2N,
			c->litr3N,
			c->litr4N,
			c->cwd_N,
			c->cwd_2N,
			c->cwd_3N,
			c->cwd_4N,
			c->soilN,
			c->soil1N,
			c->soil2N,
			c->soil3N,
			c->soil4N,
			c->soil_depth);
	/************************************************************************/
	/* print meteo variables at cell level */
	logger(g_annual_log, ",%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f",   // removed yearly mean asw, already written 
			c->years[year].yearly_mean.solar_rad     ,
			c->years[year].yearly_mean.tavg          ,
			c->years[year].yearly_mean.tmax          ,
			c->years[year].yearly_mean.tmin          ,
			c->years[year].yearly_mean.tday          ,
			c->years[year].yearly_mean.tnight        ,
			c->years[year].yearly_mean.vpd           ,
			c->years[year].yearly_mean.prcp          ,
			c->years[year].yearly_mean.tsoil         ,
			c->years[year].yearly_mean.rh_f          ,
			//c->years[year].yearly_mean.asw           ,
			c->years[year].co2Conc);
	/************************************************************************/
	/* end print */
	logger(g_annual_log,"\n");
}
#else

// modified august 2021
void EOY_print_output_cell_level(cell_t *const c, const int year, const int years_of_simulation )
{
	int layer;
	int height;
	int dbh;
	int age;
	int species;

	static int print_header = 0;

	species_t *s;

	/* return if annual logging is off*/
	if ( ! g_annual_log ) return;

	/* test */
	//note it can be used only if no other classes are added!!
	/*
	if ( !year)
	{
		c->initial_tree_layers_count = c->tree_layers_count;
		c->initial_heights_count = c->heights_count;

		for ( height = c->heights_count - 1; height >= 0 ; --height )
		{
			c->heights[height].initial_dbhs_count = c->heights[height].dbhs_count;
			for ( dbh = c->heights[height].dbhs_count - 1; dbh >= 0; --dbh )
			{
				c->heights[height].dbhs[dbh].initial_ages_count = c->heights[height].dbhs[dbh].ages_count;

				for ( age = 0; age < c->heights[height].dbhs[dbh].ages_count ; ++age )
				{
					c->heights[height].dbhs[dbh].ages[age].initial_species_count = c->heights[height].dbhs[dbh].ages[age].species_count;
				}
			}
		}
	}
	 */

	/* heading */
	if ( ! print_header )
	{
		print_header = 1;

		logger(g_annual_log, "X,Y,YEAR");

		for ( layer = c->tree_layers_count - 1; layer >= 0; --layer )

		{
			for ( height = c->heights_count - 1; height >= 0 ; --height )
			{
				if( layer == c->heights[height].height_z )
				{
					/* heading for layers */
					logger(g_annual_log,",LAYER");

					/* heading for heights value */
					logger(g_annual_log,",HEIGHT");

					for ( dbh = c->heights[height].dbhs_count - 1; dbh >= 0; --dbh )
					{
						/* heading for dbhs value */
						logger(g_annual_log, ",DBH");

						for ( age = 0; age < c->heights[height].dbhs[dbh].ages_count ; ++age )
						{
							/* heading for ages */
							logger(g_annual_log,",AGE");

							 for ( species = 0; species < c->heights[height].dbhs[dbh].ages[age].species_count; ++species )
                           
							{
								/* heading for species name */
								logger(g_annual_log,",SPECIES");
                                                             
								/* heading for management name */
								logger(g_annual_log,",MANAGEMENT");

								logger(g_annual_log,
										",GPP"
										",GPP_SUN:GPP"
										",GPP_SHADE:GPP"
										",Av_SUN:A_SUN"
										",Aj_SUN:A_SUN"
										",Av_SHADE:A_SHADE"
										",Aj_SHADE:A_SHADE"
										",Av_TOT:A_TOT"
										",Aj_TOT:A_TOT"
										",GR"
										",MR"
										",RA"
										",NPP"
										",BP"
										",reser_as_diff"
										",ResAlloc"
										",ResDeple"
										",ResUsage"
										",BP/NPP"
										",ResAlloc/NPP"
										",ResAlloc/BP"
										",ResDeple/NPP"
										",ResDeple/BP"
										",ResUsage/NPP"
										",ResUsage/BP"
										",CUE"
										",BPE"
										",diffCUE-BPE"
										",Y(perc)"
										",MAX_NSC_CONC"
										",MIN_NSC_CONC"
										",PeakLAI"
										",MaxLAI"
										",SLA"
										",SAPWOOD-AREA"
										",HEARTWOOD-AREA"
										",CC-Proj"
										",DBHDC"
										",CROWN_DIAMETER"
										",CROWN_HEIGHT"
										",CROWN_AREA_PROJ"
										",APAR"
										",LIVETREE"
										",DEADTREE"
										",THINNEDTREE"
										",VEG_D"
										",FIRST_VEG_DAY"
										",CTRANSP"
										",CINT"
										",CLE"
										",WUE"
										",MAX_ANN_RESERVE_C"
                                                                                ",MIN_ANN_RESERVE_C"
                                                                                ",TREE_MAX_ANN_RESERVE_C"
                                                                                ",TREE_MIN_ANN_RESERVE_C"
										",MIN_RESERVE_C"
										",RESERVE_C"
										",STEM_C"
										",STEMSAP_C"
										",STEMHEART_C"
										",STEMSAP_PERC"
										",STEMLIVE_C"
										",STEMDEAD_C"
										",STEMLIVE_PERC"
										",MAX_LEAF_C"
										",MAX_FROOT_C"
										",CROOT_C"
										",CROOTLIVE_C"
										",CROOTDEAD_C"
										",CROOTLIVE_PERC"
										",BRANCH_C"
										",BRANCHLIVE_C"
										",BRANCHDEAD_C"
										",BRANCHLIVE_PERC"
										",FRUIT_C"
										",MAX_FRUIT_C"
										",TOT_SAPWOOD_C"
										",TOT_HEARTWOOD_C"
										",RESERVE_N"
										",STEM_N"
										",STEMLIVE_N"
										",STEMDEAD_N"
										",CROOT_N"
										",CROOTLIVE_N"
										",CROOTDEAD_N"
										",BRANCH_N"
										",BRANCHLIVE_N"
										",BRANCHDEAD_N"
										",FRUIT_N"
										",STANDING_WOOD"
										",DELTA_WOOD"
										",CUM_DELTA_WOOD"
										",BASAL_AREA"
										",TREE_CAI"
										",TREE_MAI"
										",CAI"
										",MAI"
										",VOLUME"
										",TREE_VOLUME"
										",DELTA_TREE_VOL(perc)"
										",DELTA_AGB"
										",DELTA_BGB"
										",AGB"
										",BGB"
										",BGB:AGB"
										",DELTA_TREE_AGB"
										",DELTA_TREE_BGB"
										",C_HWP"
										",VOLUME_HWP"
										",STEM_RA"
										",LEAF_RA"
										",FROOT_RA"
										",CROOT_RA"
										",BRANCH_RA");
							}
							if ( c->heights[height].dbhs[dbh].ages[age].species_count > 1 ) {
								logger(g_annual_log,",*");
							}
						}
						if ( c->heights[height].dbhs[dbh].ages_count > 1 ) {
							logger(g_annual_log,",**");
						}
					}
					if ( c->heights[height].dbhs_count > 1 ) {
						logger(g_annual_log,",***");
					}
				}
				if ( c->tree_layers[layer].layer_n_height_class > 1 ) {
					logger(g_annual_log,",****");
				}
			}
			if ( c->tree_layers_count > 1 ) {
				logger(g_annual_log,",*****");
			}
		}
		/************************************************************************/

                 if (!g_settings->regeneration)
                 { 
		/* heading cell variables */
		logger(g_annual_log,",gpp,npp,ar,hr,rsoil,rsoilCO2,reco,nee,nep,et,le,soil-evapo,avg_asw,iWue,vol,cum_vol,run_off,"   // put avg_asw
				"litrC,litr1C,litr2C,litr3C,litr4C,cwdC,cwd2C,cwd3C,cwd4C,soilC,soil1C,soil2C,soil3C,soil4C,"
				"litrN,litr1N,litr2N,litr3N,litr4N,cwdN,cwd2N,cwd3N,cwd4N,soilN,soil1N,soil2N,soil3N,soil4N,soil_depth"); 
		/************************************************************************/
		/* heading meteo variables */
		logger(g_annual_log,",solar_rad,tavg,tmax,tmin,tday,tnight,vpd,prcp,tsoil,rh,[CO2]"); // removed avg_asw
		/************************************************************************/
                }
		/* end heading */
		logger(g_annual_log,"\n");
		/************************************************************************/

	}

	/*****************************************************************************************************/


	/* values */
	logger(g_annual_log, "%d,%d,%d", c->x, c->y, c->years[year].year);

	/* print class level values */
	if ( c->heights_count )
	{
		for ( layer = c->tree_layers_count - 1; layer >= 0; --layer )
		{
			qsort(c->heights, c->heights_count, sizeof(height_t), sort_by_heights_desc);

			for ( height = 0; height < c->heights_count; ++height )
			{
				if( layer == c->heights[height].height_z )
				{
					/* print layer */
					logger(g_annual_log,",%d", layer);

					/* print height */
					logger(g_annual_log,",%g", c->heights[height].value);

					for ( dbh = c->heights[height].dbhs_count - 1; dbh >= 0; --dbh )
					{
						/* print dbh */
						logger(g_annual_log,",%g", c->heights[height].dbhs[dbh].value);

						// ALESSIOR TO ALESSIOC
						// SHOULD THIS BE CHANGED TO
						// start from c->heights[height].dbhs[dbh].ages_count-1 ?
						for ( age = 0; age < c->heights[height].dbhs[dbh].ages_count ; ++age )
						{
							/* print age */
							logger(g_annual_log,",%d", c->heights[height].dbhs[dbh].ages[age].value);

							for ( species = 0; species < c->heights[height].dbhs[dbh].ages[age].species_count; ++species )
							{
								s  = &c->heights[height].dbhs[dbh].ages[age].species[species];

								/* print species name */
								logger(g_annual_log,",%s", c->heights[height].dbhs[dbh].ages[age].species[species].name);

								/* print management */
								logger(g_annual_log,",%c", sz_management[c->heights[height].dbhs[dbh].ages[age].species[species].management]);

								/* print variables at layer-class level */
								logger(g_annual_log,",%6.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%d,%d,%d,%d,%d,%3.4f"
										",%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f"
										",%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f"
										",%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f",
										s->value[YEARLY_GPP],
										s->value[YEARLY_GPP_SUN]   / s->value[YEARLY_GPP],
										s->value[YEARLY_GPP_SHADE] / s->value[YEARLY_GPP],
										s->value[YEARLY_Av_SUN]    / s->value[YEARLY_A_SUN],
										s->value[YEARLY_Aj_SUN]    / s->value[YEARLY_A_SUN],
										s->value[YEARLY_Av_SHADE]  / s->value[YEARLY_A_SHADE],
										s->value[YEARLY_Aj_SHADE]  / s->value[YEARLY_A_SHADE],
										s->value[YEARLY_Av_TOT]    / s->value[YEARLY_A_TOT],
										s->value[YEARLY_Aj_TOT]    / s->value[YEARLY_A_TOT],
										s->value[YEARLY_TOTAL_GROWTH_RESP],
										s->value[YEARLY_TOTAL_MAINT_RESP],
										s->value[YEARLY_TOTAL_AUT_RESP],
										s->value[YEARLY_NPP],
										s->value[YEARLY_BP],
										(s->value[YEARLY_NPP] - s->value[YEARLY_BP]),
										s->value[YEARLY_RESERVE_ALLOC],
										s->value[YEARLY_RESERVE_DEPLE],
										s->value[YEARLY_RESERVE_USAGE],
										(s->value[YEARLY_BP] / s->value[YEARLY_NPP]),
										(s->value[YEARLY_RESERVE_ALLOC] / s->value[YEARLY_NPP]),
										(s->value[YEARLY_RESERVE_ALLOC] / s->value[YEARLY_BP]),
										(s->value[YEARLY_RESERVE_DEPLE] / s->value[YEARLY_NPP]),
										(s->value[YEARLY_RESERVE_DEPLE] / s->value[YEARLY_BP]),
										(s->value[YEARLY_RESERVE_USAGE] / s->value[YEARLY_NPP]),
										(s->value[YEARLY_RESERVE_USAGE] / s->value[YEARLY_BP]),
										s->value[YEARLY_CUE],
										s->value[YEARLY_BPE],
										(s->value[YEARLY_CUE] - s->value[YEARLY_BPE]),
										s->value[YEARLY_TOTAL_AUT_RESP] / s->value[YEARLY_GPP] * 100.,
										s->value[MAX_RESERVE_C_CONC],
										s->value[MIN_RESERVE_C_CONC],
										s->value[PEAK_LAI_PROJ],
										s->value[MAX_LAI_PROJ],
										s->value[SLA_PROJ],
										s->value[SAPWOOD_AREA],
										s->value[HEARTWOOD_AREA],
										s->value[CANOPY_COVER_PROJ],
										s->value[DBHDC_EFF],
										s->value[CROWN_DIAMETER],
										s->value[CROWN_HEIGHT],
										s->value[CROWN_AREA_PROJ],
										s->value[YEARLY_APAR],
										s->counter[N_TREE],
										s->counter[DEAD_TREE],
										s->counter[THINNED_TREE],
										s->counter[YEARLY_VEG_DAYS],
										s->counter[FIRST_VEG_DAYS],
										s->value[YEARLY_CANOPY_TRANSP],
										s->value[YEARLY_CANOPY_INT],
										s->value[YEARLY_CANOPY_LATENT_HEAT],
										s->value[YEARLY_WUE],
										s->value[MAX_ANN_RESERVE_C],
                                                                               s->value[MIN_ANN_RESERVE_C],      //5p6
                                                                               s->value[TREE_MAX_ANN_RESERVE_C], //5p6
                                                                               s->value[TREE_MIN_ANN_RESERVE_C], //5p6 
										s->value[MIN_RESERVE_C],
										s->value[RESERVE_C],
										s->value[STEM_C],
										s->value[STEM_SAPWOOD_C],
										s->value[STEM_HEARTWOOD_C],
										(s->value[STEM_SAPWOOD_C] * 100. ) / s->value[STEM_C],
										s->value[STEM_LIVEWOOD_C],
										s->value[STEM_DEADWOOD_C],
										(s->value[STEM_LIVEWOOD_C] * 100. ) / s->value[STEM_C],
										s->value[MAX_LEAF_C],
										s->value[MAX_FROOT_C],
										s->value[CROOT_C],
										s->value[CROOT_LIVEWOOD_C],
										s->value[CROOT_DEADWOOD_C],
										(s->value[CROOT_LIVEWOOD_C] * 100. ) / s->value[CROOT_C],
										s->value[BRANCH_C],
										s->value[BRANCH_LIVEWOOD_C],
										s->value[BRANCH_DEADWOOD_C],
										(s->value[BRANCH_LIVEWOOD_C] * 100. ) / s->value[BRANCH_C],
										s->value[FRUIT_C],
										s->value[MAX_FRUIT_C],
										s->value[TOT_SAPWOOD_C],
										s->value[TOT_HEARTWOOD_C],
										s->value[RESERVE_N],
										s->value[STEM_N],
										s->value[STEM_LIVEWOOD_N],
										s->value[STEM_DEADWOOD_N],
										s->value[CROOT_N],
										s->value[CROOT_LIVEWOOD_N],
										s->value[CROOT_DEADWOOD_N],
										s->value[BRANCH_N],
										s->value[BRANCH_LIVEWOOD_N],
										s->value[BRANCH_DEADWOOD_N],
										s->value[FRUIT_N],
										s->value[STANDING_WOOD],
										s->value[YEARLY_C_TO_WOOD],
										s->value[CUM_YEARLY_C_TO_WOOD],
										s->value[STAND_BASAL_AREA_m2],
										s->value[TREE_CAI],
										s->value[TREE_MAI],
										s->value[CAI],
										s->value[MAI],
										s->value[VOLUME],
										s->value[TREE_VOLUME],
										(s->value[TREE_CAI]*100.)/s->value[TREE_VOLUME],
										s->value[DELTA_AGB],
										s->value[DELTA_BGB],
										s->value[AGB],
										s->value[BGB],
										s->value[BGB]/s->value[AGB],
										s->value[DELTA_TREE_AGB],
										s->value[DELTA_TREE_BGB],
										s->value[C_HWP],
										s->value[VOLUME_HWP],
										s->value[YEARLY_STEM_AUT_RESP],
										s->value[YEARLY_LEAF_AUT_RESP],
										s->value[YEARLY_FROOT_AUT_RESP],
										s->value[YEARLY_CROOT_AUT_RESP],
										s->value[YEARLY_BRANCH_AUT_RESP]);
							}
							if ( c->heights[height].dbhs[dbh].ages[age].species_count > 1 ) {
								logger(g_annual_log,",*");
							}
						}
						if ( c->heights[height].dbhs[dbh].ages_count > 1 ) {
							logger(g_annual_log,",**");
						}
					}
					if ( c->heights[height].dbhs_count > 1 ) {
						logger(g_annual_log,",***");
					}
				}
				if ( c->tree_layers[layer].layer_n_height_class > 1 ) {
					logger(g_annual_log,",****");
				}
			}
			if ( c->tree_layers_count > 1 )
			{
				logger(g_annual_log,",*****");
			}
		}
	}
	else
	{
		//ALESSIOC TO ALLESSIOR PRINT EMPTY SPACES WHEN N_TREE = 0
	}
	/************************************************************************/
	/* printing variables at cell level only if there's more than one layer */  //Note: this is actually not true
 if (!g_settings->regeneration)
                 {
        logger(g_annual_log, ",%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f"
			",%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f",
			c->annual_gpp,
			c->annual_npp,
			c->annual_aut_resp,
			c->annual_het_resp,
			c->annual_soil_resp,
			c->annual_soil_respCO2,
			c->annual_r_eco,
			c->annual_nee,
			c->annual_nep,
			c->annual_et,
			c->annual_lh_flux,
			c->annual_soil_evapo,
                        c->years[year].yearly_mean.asw,
			//c->asw,   // comment: this is actually the value at the 31 december, it has more sense to have the mean value
			c->annual_iwue,
			c->volume,
			c->cum_volume,
			c->annual_out_flow,
			c->litrC,
			c->litr1C,
			c->litr2C,
			c->litr3C,
			c->litr4C,
			c->cwd_C,
			c->cwd_2C,
			c->cwd_3C,
			c->cwd_4C,
			c->soilC,
			c->soil1C,
			c->soil2C,
			c->soil3C,
			c->soil4C,
			c->litrN,
			c->litr1N,
			c->litr2N,
			c->litr3N,
			c->litr4N,
			c->cwd_N,
			c->cwd_2N,
			c->cwd_3N,
			c->cwd_4N,
			c->soilN,
			c->soil1N,
			c->soil2N,
			c->soil3N,
			c->soil4N,
			c->soil_depth);
	/************************************************************************/
	/* print meteo variables at cell level */
	logger(g_annual_log, ",%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f",   // removed yearly mean asw, already written 
			c->years[year].yearly_mean.solar_rad     ,
			c->years[year].yearly_mean.tavg          ,
			c->years[year].yearly_mean.tmax          ,
			c->years[year].yearly_mean.tmin          ,
			c->years[year].yearly_mean.tday          ,
			c->years[year].yearly_mean.tnight        ,
			c->years[year].yearly_mean.vpd           ,
			c->years[year].yearly_mean.prcp          ,
			c->years[year].yearly_mean.tsoil         ,
			c->years[year].yearly_mean.rh_f          ,
			//c->years[year].yearly_mean.asw           ,
			c->years[year].co2Conc);
	/************************************************************************/

        }
        /* end print */
	logger(g_annual_log,"\n");
}
#endif

void EOY_print_output_cell_level_ddalmo(cell_t *const c, const int year, const int years_of_simulation )
{
	int layer;
	int height;
	int dbh;
	int age;
	int species;

	static int print_header = 0;

	species_t *s;

	/* return if annual logging is off*/
	if ( ! g_annual_log ) return;

	/* test */
	//note it can be used only if no other classes are added!!
	/*
	if ( !year)
	{
		c->initial_tree_layers_count = c->tree_layers_count;
		c->initial_heights_count = c->heights_count;

		for ( height = c->heights_count - 1; height >= 0 ; --height )
		{
			c->heights[height].initial_dbhs_count = c->heights[height].dbhs_count;
			for ( dbh = c->heights[height].dbhs_count - 1; dbh >= 0; --dbh )
			{
				c->heights[height].dbhs[dbh].initial_ages_count = c->heights[height].dbhs[dbh].ages_count;

				for ( age = 0; age < c->heights[height].dbhs[dbh].ages_count ; ++age )
				{
					c->heights[height].dbhs[dbh].ages[age].initial_species_count = c->heights[height].dbhs[dbh].ages[age].species_count;
				}
			}
		}
	}
	 */

	/* heading */
	if ( ! print_header )
	{
	  print_header = 1;

           for ( layer = c->tree_layers_count - 1; layer >= 0; --layer )

		{

		logger(g_annual_log, "X,Y,YEAR");

		//for ( layer = c->tree_layers_count - 1; layer >= 0; --layer )

		//{
			for ( height = c->heights_count - 1; height >= 0 ; --height )
			{
				if( layer == c->heights[height].height_z )
				{
					/* heading for layers */
					logger(g_annual_log,",LAYER");

					/* heading for heights value */
					logger(g_annual_log,",HEIGHT");

					for ( dbh = c->heights[height].dbhs_count - 1; dbh >= 0; --dbh )
					{
						/* heading for dbhs value */
						logger(g_annual_log, ",DBH");

						for ( age = 0; age < c->heights[height].dbhs[dbh].ages_count ; ++age )
						{
							/* heading for ages */
							logger(g_annual_log,",AGE");

							 for ( species = 0; species < c->heights[height].dbhs[dbh].ages[age].species_count; ++species )
                           
							{
								/* heading for species name */
								logger(g_annual_log,",SPECIES");
                                                             
								/* heading for management name */
								logger(g_annual_log,",MANAGEMENT");

								logger(g_annual_log,
										",GPP"
										",GPP_SUN:GPP"
										",GPP_SHADE:GPP"
										",Av_SUN:A_SUN"
										",Aj_SUN:A_SUN"
										",Av_SHADE:A_SHADE"
										",Aj_SHADE:A_SHADE"
										",Av_TOT:A_TOT"
										",Aj_TOT:A_TOT"
										",GR"
										",MR"
										",RA"
										",NPP"
										",BP"
										",reser_as_diff"
										",ResAlloc"
										",ResDeple"
										",ResUsage"
										",BP/NPP"
										",ResAlloc/NPP"
										",ResAlloc/BP"
										",ResDeple/NPP"
										",ResDeple/BP"
										",ResUsage/NPP"
										",ResUsage/BP"
										",CUE"
										",BPE"
										",diffCUE-BPE"
										",Y(perc)"
										",MAX_NSC_CONC"
										",MIN_NSC_CONC"
										",PeakLAI"
										",MaxLAI"
										",SLA"
										",SAPWOOD-AREA"
										",HEARTWOOD-AREA"
										",CC-Proj"
										",DBHDC"
										",CROWN_DIAMETER"
										",CROWN_HEIGHT"
										",CROWN_AREA_PROJ"
										",APAR"
										",LIVETREE"
										",DEADTREE"
										",THINNEDTREE"
										",VEG_D"
										",FIRST_VEG_DAY"
										",CTRANSP"
										",CINT"
										",CLE"
										",WUE"
										",MAX_ANN_RESERVE_C"
                                                                                ",MIN_ANN_RESERVE_C"
                                                                                ",TREE_MAX_ANN_RESERVE_C"
                                                                                ",TREE_MIN_ANN_RESERVE_C"
										",MIN_RESERVE_C"
										",RESERVE_C"
										",STEM_C"
										",STEMSAP_C"
										",STEMHEART_C"
										",STEMSAP_PERC"
										",STEMLIVE_C"
										",STEMDEAD_C"
										",STEMLIVE_PERC"
										",MAX_LEAF_C"
										",MAX_FROOT_C"
										",CROOT_C"
										",CROOTLIVE_C"
										",CROOTDEAD_C"
										",CROOTLIVE_PERC"
										",BRANCH_C"
										",BRANCHLIVE_C"
										",BRANCHDEAD_C"
										",BRANCHLIVE_PERC"
										",FRUIT_C"
										",MAX_FRUIT_C"
										",TOT_SAPWOOD_C"
										",TOT_HEARTWOOD_C"
										",RESERVE_N"
										",STEM_N"
										",STEMLIVE_N"
										",STEMDEAD_N"
										",CROOT_N"
										",CROOTLIVE_N"
										",CROOTDEAD_N"
										",BRANCH_N"
										",BRANCHLIVE_N"
										",BRANCHDEAD_N"
										",FRUIT_N"
										",STANDING_WOOD"
										",DELTA_WOOD"
										",CUM_DELTA_WOOD"
										",BASAL_AREA"
										",TREE_CAI"
										",TREE_MAI"
										",CAI"
										",MAI"
										",VOLUME"
										",TREE_VOLUME"
										",DELTA_TREE_VOL(perc)"
										",DELTA_AGB"
										",DELTA_BGB"
										",AGB"
										",BGB"
										",BGB:AGB"
										",DELTA_TREE_AGB"
										",DELTA_TREE_BGB"
										",C_HWP"
										",VOLUME_HWP"
										",STEM_RA"
										",LEAF_RA"
										",FROOT_RA"
										",CROOT_RA"
										",BRANCH_RA");
							}
							if ( c->heights[height].dbhs[dbh].ages[age].species_count > 1 ) {
								logger(g_annual_log,",*");
							}
						}
						if ( c->heights[height].dbhs[dbh].ages_count > 1 ) {
							logger(g_annual_log,",**");
						}
					}
					if ( c->heights[height].dbhs_count > 1 ) {
						logger(g_annual_log,",***");
					}
				}
				if ( c->tree_layers[layer].layer_n_height_class > 1 ) {
					logger(g_annual_log,",****");
				}
			}
			//if ( c->tree_layers_count > 1 ) {
			//	logger(g_annual_log,",*****");
			//}
		//} //loop on layer
		/************************************************************************/

                 if (!g_settings->regeneration)
                 { 
		  /* heading cell variables */
		  logger(g_annual_log,",gpp,npp,ar,hr,rsoil,rsoilCO2,reco,nee,nep,et,le,soil-evapo,avg_asw,iWue,vol,cum_vol,run_off,"   // put avg_asw
				"litrC,litr1C,litr2C,litr3C,litr4C,cwdC,cwd2C,cwd3C,cwd4C,soilC,soil1C,soil2C,soil3C,soil4C,"
				"litrN,litr1N,litr2N,litr3N,litr4N,cwdN,cwd2N,cwd3N,cwd4N,soilN,soil1N,soil2N,soil3N,soil4N,soil_depth"); 
		  /************************************************************************/
		  /* heading meteo variables */
		  logger(g_annual_log,",solar_rad,tavg,tmax,tmin,tday,tnight,vpd,prcp,tsoil,rh,[CO2]"); // removed avg_asw
		  /************************************************************************/
                 }
		/* end heading */
		logger(g_annual_log,"\n");

                } // of loop on layer

		// /* end heading */
		//logger(g_annual_log,"\n");
		/************************************************************************/

	}

	/*****************************************************************************************************/


	/* values */
	//logger(g_annual_log, "%d,%d,%d", c->x, c->y, c->years[year].year);

	/* print class level values */
	if ( c->heights_count )    
	{
		for ( layer = c->tree_layers_count - 1; layer >= 0; --layer )
		{
			qsort(c->heights, c->heights_count, sizeof(height_t), sort_by_heights_desc);

			for ( height = 0; height < c->heights_count; ++height )
			{
				if( layer == c->heights[height].height_z )
				{
 
                                        /* values */
                                	logger(g_annual_log, "%d,%d,%d", c->x, c->y, c->years[year].year);

					/* print layer */
					logger(g_annual_log,",%d", layer);

					/* print height */
					logger(g_annual_log,",%g", c->heights[height].value);

					for ( dbh = c->heights[height].dbhs_count - 1; dbh >= 0; --dbh )
					{
						/* print dbh */
						logger(g_annual_log,",%g", c->heights[height].dbhs[dbh].value);

						// ALESSIOR TO ALESSIOC
						// SHOULD THIS BE CHANGED TO
						// start from c->heights[height].dbhs[dbh].ages_count-1 ?
						for ( age = 0; age < c->heights[height].dbhs[dbh].ages_count ; ++age )
						{
							/* print age */
							logger(g_annual_log,",%d", c->heights[height].dbhs[dbh].ages[age].value);

							for ( species = 0; species < c->heights[height].dbhs[dbh].ages[age].species_count; ++species )
							{
								s  = &c->heights[height].dbhs[dbh].ages[age].species[species];

								/* print species name */
								logger(g_annual_log,",%s", c->heights[height].dbhs[dbh].ages[age].species[species].name);

								/* print management */
								logger(g_annual_log,",%c", sz_management[c->heights[height].dbhs[dbh].ages[age].species[species].management]);

								/* print variables at layer-class level */
								logger(g_annual_log,",%6.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%d,%d,%d,%d,%d,%3.4f"
										",%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f"
										",%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f"
										",%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f",
										s->value[YEARLY_GPP],
										s->value[YEARLY_GPP_SUN]   / s->value[YEARLY_GPP],
										s->value[YEARLY_GPP_SHADE] / s->value[YEARLY_GPP],
										s->value[YEARLY_Av_SUN]    / s->value[YEARLY_A_SUN],
										s->value[YEARLY_Aj_SUN]    / s->value[YEARLY_A_SUN],
										s->value[YEARLY_Av_SHADE]  / s->value[YEARLY_A_SHADE],
										s->value[YEARLY_Aj_SHADE]  / s->value[YEARLY_A_SHADE],
										s->value[YEARLY_Av_TOT]    / s->value[YEARLY_A_TOT],
										s->value[YEARLY_Aj_TOT]    / s->value[YEARLY_A_TOT],
										s->value[YEARLY_TOTAL_GROWTH_RESP],
										s->value[YEARLY_TOTAL_MAINT_RESP],
										s->value[YEARLY_TOTAL_AUT_RESP],
										s->value[YEARLY_NPP],
										s->value[YEARLY_BP],
										(s->value[YEARLY_NPP] - s->value[YEARLY_BP]),
										s->value[YEARLY_RESERVE_ALLOC],
										s->value[YEARLY_RESERVE_DEPLE],
										s->value[YEARLY_RESERVE_USAGE],
										(s->value[YEARLY_BP] / s->value[YEARLY_NPP]),
										(s->value[YEARLY_RESERVE_ALLOC] / s->value[YEARLY_NPP]),
										(s->value[YEARLY_RESERVE_ALLOC] / s->value[YEARLY_BP]),
										(s->value[YEARLY_RESERVE_DEPLE] / s->value[YEARLY_NPP]),
										(s->value[YEARLY_RESERVE_DEPLE] / s->value[YEARLY_BP]),
										(s->value[YEARLY_RESERVE_USAGE] / s->value[YEARLY_NPP]),
										(s->value[YEARLY_RESERVE_USAGE] / s->value[YEARLY_BP]),
										s->value[YEARLY_CUE],
										s->value[YEARLY_BPE],
										(s->value[YEARLY_CUE] - s->value[YEARLY_BPE]),
										s->value[YEARLY_TOTAL_AUT_RESP] / s->value[YEARLY_GPP] * 100.,
										s->value[MAX_RESERVE_C_CONC],
										s->value[MIN_RESERVE_C_CONC],
										s->value[PEAK_LAI_PROJ],
										s->value[MAX_LAI_PROJ],
										s->value[SLA_PROJ],
										s->value[SAPWOOD_AREA],
										s->value[HEARTWOOD_AREA],
										s->value[CANOPY_COVER_PROJ],
										s->value[DBHDC_EFF],
										s->value[CROWN_DIAMETER],
										s->value[CROWN_HEIGHT],
										s->value[CROWN_AREA_PROJ],
										s->value[YEARLY_APAR],
										s->counter[N_TREE],
										s->counter[DEAD_TREE],
										s->counter[THINNED_TREE],
										s->counter[YEARLY_VEG_DAYS],
										s->counter[FIRST_VEG_DAYS],
										s->value[YEARLY_CANOPY_TRANSP],
										s->value[YEARLY_CANOPY_INT],
										s->value[YEARLY_CANOPY_LATENT_HEAT],
										s->value[YEARLY_WUE],
										s->value[MAX_ANN_RESERVE_C],
                                                                               s->value[MIN_ANN_RESERVE_C],       //5p6
                                                                               s->value[TREE_MAX_ANN_RESERVE_C],  //5p6
                                                                               s->value[TREE_MIN_ANN_RESERVE_C],  //5p6 
										s->value[MIN_RESERVE_C],
										s->value[RESERVE_C],
										s->value[STEM_C],
										s->value[STEM_SAPWOOD_C],
										s->value[STEM_HEARTWOOD_C],
										(s->value[STEM_SAPWOOD_C] * 100. ) / s->value[STEM_C],
										s->value[STEM_LIVEWOOD_C],
										s->value[STEM_DEADWOOD_C],
										(s->value[STEM_LIVEWOOD_C] * 100. ) / s->value[STEM_C],
										s->value[MAX_LEAF_C],
										s->value[MAX_FROOT_C],
										s->value[CROOT_C],
										s->value[CROOT_LIVEWOOD_C],
										s->value[CROOT_DEADWOOD_C],
										(s->value[CROOT_LIVEWOOD_C] * 100. ) / s->value[CROOT_C],
										s->value[BRANCH_C],
										s->value[BRANCH_LIVEWOOD_C],
										s->value[BRANCH_DEADWOOD_C],
										(s->value[BRANCH_LIVEWOOD_C] * 100. ) / s->value[BRANCH_C],
										s->value[FRUIT_C],
										s->value[MAX_FRUIT_C],
										s->value[TOT_SAPWOOD_C],
										s->value[TOT_HEARTWOOD_C],
										s->value[RESERVE_N],
										s->value[STEM_N],
										s->value[STEM_LIVEWOOD_N],
										s->value[STEM_DEADWOOD_N],
										s->value[CROOT_N],
										s->value[CROOT_LIVEWOOD_N],
										s->value[CROOT_DEADWOOD_N],
										s->value[BRANCH_N],
										s->value[BRANCH_LIVEWOOD_N],
										s->value[BRANCH_DEADWOOD_N],
										s->value[FRUIT_N],
										s->value[STANDING_WOOD],
										s->value[YEARLY_C_TO_WOOD],
										s->value[CUM_YEARLY_C_TO_WOOD],
										s->value[STAND_BASAL_AREA_m2],
										s->value[TREE_CAI],
										s->value[TREE_MAI],
										s->value[CAI],
										s->value[MAI],
										s->value[VOLUME],
										s->value[TREE_VOLUME],
										(s->value[TREE_CAI]*100.)/s->value[TREE_VOLUME],
										s->value[DELTA_AGB],
										s->value[DELTA_BGB],
										s->value[AGB],
										s->value[BGB],
										s->value[BGB]/s->value[AGB],
										s->value[DELTA_TREE_AGB],
										s->value[DELTA_TREE_BGB],
										s->value[C_HWP],
										s->value[VOLUME_HWP],
										s->value[YEARLY_STEM_AUT_RESP],
										s->value[YEARLY_LEAF_AUT_RESP],
										s->value[YEARLY_FROOT_AUT_RESP],
										s->value[YEARLY_CROOT_AUT_RESP],
										s->value[YEARLY_BRANCH_AUT_RESP]);
							}
							if ( c->heights[height].dbhs[dbh].ages[age].species_count > 1 ) {
								logger(g_annual_log,",*");
							}
						}
						if ( c->heights[height].dbhs[dbh].ages_count > 1 ) {
							logger(g_annual_log,",**");
						}
					}
					if ( c->heights[height].dbhs_count > 1 ) {
						logger(g_annual_log,",***");
					}
				}
				if ( c->tree_layers[layer].layer_n_height_class > 1 ) {
					logger(g_annual_log,",****");
				}
			}
			//if ( c->tree_layers_count > 1 )
			//{
				//logger(g_annual_log,",*****");
			//}
                
                        /************************************************************************/
	          /* printing variables at cell level only if there's more than one layer */  //Note: this is actually not true
                  if (!g_settings->regeneration)
                   {
                        logger(g_annual_log, ",%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f"
			",%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f",
			c->annual_gpp,
			c->annual_npp,
			c->annual_aut_resp,
			c->annual_het_resp,
			c->annual_soil_resp,
			c->annual_soil_respCO2,
			c->annual_r_eco,
			c->annual_nee,
			c->annual_nep,
			c->annual_et,
			c->annual_lh_flux,
			c->annual_soil_evapo,
                        c->years[year].yearly_mean.asw,
			//c->asw,   // comment: this is actually the value at the 31 december, it has more sense to have the mean value
			c->annual_iwue,
			c->volume,
			c->cum_volume,
			c->annual_out_flow,
			c->litrC,
			c->litr1C,
			c->litr2C,
			c->litr3C,
			c->litr4C,
			c->cwd_C,
			c->cwd_2C,
			c->cwd_3C,
			c->cwd_4C,
			c->soilC,
			c->soil1C,
			c->soil2C,
			c->soil3C,
			c->soil4C,
			c->litrN,
			c->litr1N,
			c->litr2N,
			c->litr3N,
			c->litr4N,
			c->cwd_N,
			c->cwd_2N,
			c->cwd_3N,
			c->cwd_4N,
			c->soilN,
			c->soil1N,
			c->soil2N,
			c->soil3N,
			c->soil4N,
			c->soil_depth);
	            /************************************************************************/
	            /* print meteo variables at cell level */
	            logger(g_annual_log, ",%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f",   // removed yearly mean asw, already written 
			c->years[year].yearly_mean.solar_rad     ,
			c->years[year].yearly_mean.tavg          ,
			c->years[year].yearly_mean.tmax          ,
			c->years[year].yearly_mean.tmin          ,
			c->years[year].yearly_mean.tday          ,
			c->years[year].yearly_mean.tnight        ,
			c->years[year].yearly_mean.vpd           ,
			c->years[year].yearly_mean.prcp          ,
			c->years[year].yearly_mean.tsoil         ,
			c->years[year].yearly_mean.rh_f          ,
			//c->years[year].yearly_mean.asw           ,
			c->years[year].co2Conc);
	            /************************************************************************/

                   }
                   
                   /* end print */
	           logger(g_annual_log,"\n");

		} //end layer 

	}
	else   // no height_count
	{
		//ALESSIOC TO ALLESSIOR PRINT EMPTY SPACES WHEN N_TREE = 0
	}
}
void EOY_print_output_cell_level_mc(cell_t *const c, const int year, const int years_of_simulation )
{

// this version allow to print for each output date, for each row the output according to each class
	int layer;
	int height;
	int dbh;
	int age;
	int species;

	static int print_header = 0;
	
	int print_cell = 0;   // print cell level value

	species_t *s;

	/* return if annual logging is off*/
	if ( ! g_annual_log ) return;

	/* test */
	//note it can be used only if no other classes are added!!
	/*
	if ( !year)
	{
		c->initial_tree_layers_count = c->tree_layers_count;
		c->initial_heights_count = c->heights_count;

		for ( height = c->heights_count - 1; height >= 0 ; --height )
		{
			c->heights[height].initial_dbhs_count = c->heights[height].dbhs_count;
			for ( dbh = c->heights[height].dbhs_count - 1; dbh >= 0; --dbh )
			{
				c->heights[height].dbhs[dbh].initial_ages_count = c->heights[height].dbhs[dbh].ages_count;

				for ( age = 0; age < c->heights[height].dbhs[dbh].ages_count ; ++age )
				{
					c->heights[height].dbhs[dbh].ages[age].initial_species_count = c->heights[height].dbhs[dbh].ages[age].species_count;
				}
			}
		}
	}
	 */

	/* heading */
	if ( ! print_header )
	{
	  print_header = 1;
		
		// 5p6 The output is written for each class in separated rows.  
                // at the end the information at cell level is provided in each row. This of course is redundant
                // TODO: save in a separate file 
		
                logger(g_annual_log, "X,Y,YEAR,LAYER,HEIGHT,DBH,AGE,SPECIE,MANAGEMENT");
                logger(g_annual_log,
										",GPP"
										",GPP_SUN:GPP"
										",GPP_SHADE:GPP"
										",Av_SUN:A_SUN"
										",Aj_SUN:A_SUN"
										",Av_SHADE:A_SHADE"
										",Aj_SHADE:A_SHADE"
										",Av_TOT:A_TOT"
										",Aj_TOT:A_TOT"
										",GR"
										",MR"
										",RA"
										",NPP"
										",BP"
										",reser_as_diff"
										",ResAlloc"
										",ResDeple"
										",ResUsage"
										",BP/NPP"
										",ResAlloc/NPP"
										",ResAlloc/BP"
										",ResDeple/NPP"
										",ResDeple/BP"
										",ResUsage/NPP"
										",ResUsage/BP"
										",CUE"
										",BPE"
										",diffCUE-BPE"
										",Y(perc)"
										",MAX_NSC_CONC"
										",MIN_NSC_CONC"
										",PeakLAI"
										",MaxLAI"
										",SLA"
										",SAPWOOD-AREA"
										",HEARTWOOD-AREA"
										",CC-Proj"
										",DBHDC"
										",CROWN_DIAMETER"
										",CROWN_HEIGHT"
										",CROWN_AREA_PROJ"
										",APAR"
										",LIVETREE"
										",DEADTREE"
										",THINNEDTREE"
										",VEG_D"
										",FIRST_VEG_DAY"
										",CTRANSP"
										",CINT"
										",CLE"
										",WUE"
										",MAX_ANN_RESERVE_C"
                                                                                ",MIN_ANN_RESERVE_C"
                                                                                ",TREE_MAX_ANN_RESERVE_C"
                                                                                ",TREE_MIN_ANN_RESERVE_C"
										",MIN_RESERVE_C"
										",RESERVE_C"
										",STEM_C"
										",STEMSAP_C"
										",STEMHEART_C"
										",STEMSAP_PERC"
										",STEMLIVE_C"
										",STEMDEAD_C"
										",STEMLIVE_PERC"
										",MAX_LEAF_C"
										",MAX_FROOT_C"
										",CROOT_C"
										",CROOTLIVE_C"
										",CROOTDEAD_C"
										",CROOTLIVE_PERC"
										",BRANCH_C"
										",BRANCHLIVE_C"
										",BRANCHDEAD_C"
										",BRANCHLIVE_PERC"
										",FRUIT_C"
										",MAX_FRUIT_C"
										",TOT_SAPWOOD_C"
										",TOT_HEARTWOOD_C"
										",RESERVE_N"
										",STEM_N"
										",STEMLIVE_N"
										",STEMDEAD_N"
										",CROOT_N"
										",CROOTLIVE_N"
										",CROOTDEAD_N"
										",BRANCH_N"
										",BRANCHLIVE_N"
										",BRANCHDEAD_N"
										",FRUIT_N"
										",STANDING_WOOD"
										",DELTA_WOOD"
										",CUM_DELTA_WOOD"
										",BASAL_AREA"
										",TREE_CAI"
										",TREE_MAI"
										",CAI"
										",MAI"
										",VOLUME"
										",TREE_VOLUME"
										",DELTA_TREE_VOL(perc)"
										",DELTA_AGB"
										",DELTA_BGB"
										",AGB"
										",BGB"
										",BGB:AGB"
										",DELTA_TREE_AGB"
										",DELTA_TREE_BGB"
										",C_HWP"
										",VOLUME_HWP"
										",STEM_RA"
										",LEAF_RA"
										",FROOT_RA"
										",CROOT_RA"
										",BRANCH_RA");

		/************************************************************************/

               if (!print_cell)
                 { 
		  /* heading cell variables */
		  logger(g_annual_log,",gpp,npp,ar,hr,rsoil,rsoilCO2,reco,nee,nep,et,le,soil-evapo,avg_asw,iWue,vol,cum_vol,run_off,"  
				"litrC,litr1C,litr2C,litr3C,litr4C,cwdC,cwd2C,cwd3C,cwd4C,soilC,soil1C,soil2C,soil3C,soil4C,"
				"litrN,litr1N,litr2N,litr3N,litr4N,cwdN,cwd2N,cwd3N,cwd4N,soilN,soil1N,soil2N,soil3N,soil4N,soil_depth"); 
		  /************************************************************************/
		  /* heading meteo variables */
		  logger(g_annual_log,",solar_rad,tavg,tmax,tmin,tday,tnight,vpd,prcp,tsoil,rh,[CO2]"); 
		  /************************************************************************/
                 }
		/* end heading */
		logger(g_annual_log,"\n");
		/************************************************************************/

	}

	/*****************************************************************************************************/


	/* values */
	//logger(g_annual_log, "%d,%d,%d", c->x, c->y, c->years[year].year);

	/* print class level values */
	
	//if ( c->heights_count )    //
	//{
		for ( layer = c->tree_layers_count - 1; layer >= 0; --layer )
		{
			qsort(c->heights, c->heights_count, sizeof(height_t), sort_by_heights_desc);

			for ( height = 0; height < c->heights_count; ++height )
			{
				if( layer == c->heights[height].height_z )     // check if height is in the layer considered
				{
 
                                        /* values */
                                	//logger(g_annual_log, "%d,%d,%d", c->x, c->y, c->years[year].year);

					/* print layer */
					//logger(g_annual_log,",%d", layer);

					/* print height */
					//logger(g_annual_log,",%g", c->heights[height].value);

					for ( dbh = c->heights[height].dbhs_count - 1; dbh >= 0; --dbh )
					{
						/* print dbh */
						//logger(g_annual_log,",%g", c->heights[height].dbhs[dbh].value);

						// ALESSIOR TO ALESSIOC
						// SHOULD THIS BE CHANGED TO
						// start from c->heights[height].dbhs[dbh].ages_count-1 ?
						for ( age = 0; age < c->heights[height].dbhs[dbh].ages_count ; ++age )
						{
							/* print age */
							//logger(g_annual_log,",%d", c->heights[height].dbhs[dbh].ages[age].value);

                                                       // LOOP PIU' INTERNO
							for ( species = 0; species < c->heights[height].dbhs[dbh].ages[age].species_count; ++species )
							{
								s  = &c->heights[height].dbhs[dbh].ages[age].species[species];

  
                                                                /* values */
                                	                 logger(g_annual_log, "%d,%d,%d", c->x, c->y, c->years[year].year);

					                 /* print layer */
					                 logger(g_annual_log,",%d", layer);

					                 /* print height */
					                 logger(g_annual_log,",%g", c->heights[height].value);
					                 
                                                        /* print dbh */
						         logger(g_annual_log,",%g", c->heights[height].dbhs[dbh].value);
						         
						          /* print age */
							logger(g_annual_log,",%d", c->heights[height].dbhs[dbh].ages[age].value);
             

								/* print species name */
								logger(g_annual_log,",%s", c->heights[height].dbhs[dbh].ages[age].species[species].name);

								/* print management */
								logger(g_annual_log,",%c", sz_management[c->heights[height].dbhs[dbh].ages[age].species[species].management]);

								/* print variables at layer-class level */
								logger(g_annual_log,",%6.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%d,%d,%d,%d,%d,%3.4f"
										",%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f"
										",%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f"
										",%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f",
										s->value[YEARLY_GPP],
										s->value[YEARLY_GPP_SUN]   / s->value[YEARLY_GPP],
										s->value[YEARLY_GPP_SHADE] / s->value[YEARLY_GPP],
										s->value[YEARLY_Av_SUN]    / s->value[YEARLY_A_SUN],
										s->value[YEARLY_Aj_SUN]    / s->value[YEARLY_A_SUN],
										s->value[YEARLY_Av_SHADE]  / s->value[YEARLY_A_SHADE],
										s->value[YEARLY_Aj_SHADE]  / s->value[YEARLY_A_SHADE],
										s->value[YEARLY_Av_TOT]    / s->value[YEARLY_A_TOT],
										s->value[YEARLY_Aj_TOT]    / s->value[YEARLY_A_TOT],
										s->value[YEARLY_TOTAL_GROWTH_RESP],
										s->value[YEARLY_TOTAL_MAINT_RESP],
										s->value[YEARLY_TOTAL_AUT_RESP],
										s->value[YEARLY_NPP],
										s->value[YEARLY_BP],
										(s->value[YEARLY_NPP] - s->value[YEARLY_BP]),
										s->value[YEARLY_RESERVE_ALLOC],
										s->value[YEARLY_RESERVE_DEPLE],
										s->value[YEARLY_RESERVE_USAGE],
										(s->value[YEARLY_BP] / s->value[YEARLY_NPP]),
										(s->value[YEARLY_RESERVE_ALLOC] / s->value[YEARLY_NPP]),
										(s->value[YEARLY_RESERVE_ALLOC] / s->value[YEARLY_BP]),
										(s->value[YEARLY_RESERVE_DEPLE] / s->value[YEARLY_NPP]),
										(s->value[YEARLY_RESERVE_DEPLE] / s->value[YEARLY_BP]),
										(s->value[YEARLY_RESERVE_USAGE] / s->value[YEARLY_NPP]),
										(s->value[YEARLY_RESERVE_USAGE] / s->value[YEARLY_BP]),
										s->value[YEARLY_CUE],
										s->value[YEARLY_BPE],
										(s->value[YEARLY_CUE] - s->value[YEARLY_BPE]),
										s->value[YEARLY_TOTAL_AUT_RESP] / s->value[YEARLY_GPP] * 100.,
										s->value[MAX_RESERVE_C_CONC],
										s->value[MIN_RESERVE_C_CONC],
										s->value[PEAK_LAI_PROJ],
										s->value[MAX_LAI_PROJ],
										s->value[SLA_PROJ],
										s->value[SAPWOOD_AREA],
										s->value[HEARTWOOD_AREA],
										s->value[CANOPY_COVER_PROJ],
										s->value[DBHDC_EFF],
										s->value[CROWN_DIAMETER],
										s->value[CROWN_HEIGHT],
										s->value[CROWN_AREA_PROJ],
										s->value[YEARLY_APAR],
										s->counter[N_TREE],
										s->counter[DEAD_TREE],
										s->counter[THINNED_TREE],
										s->counter[YEARLY_VEG_DAYS],
										s->counter[FIRST_VEG_DAYS],
										s->value[YEARLY_CANOPY_TRANSP],
										s->value[YEARLY_CANOPY_INT],
										s->value[YEARLY_CANOPY_LATENT_HEAT],
										s->value[YEARLY_WUE],
										s->value[MAX_ANN_RESERVE_C],
                                                                               s->value[MIN_ANN_RESERVE_C],       //5p6
                                                                               s->value[TREE_MAX_ANN_RESERVE_C],  //5p6
                                                                               s->value[TREE_MIN_ANN_RESERVE_C],  //5p6
                                                                               s->value[MIN_RESERVE_C],
										s->value[RESERVE_C],
										s->value[STEM_C],
										s->value[STEM_SAPWOOD_C],
										s->value[STEM_HEARTWOOD_C],
										(s->value[STEM_SAPWOOD_C] * 100. ) / s->value[STEM_C],
										s->value[STEM_LIVEWOOD_C],
										s->value[STEM_DEADWOOD_C],
										(s->value[STEM_LIVEWOOD_C] * 100. ) / s->value[STEM_C],
										s->value[MAX_LEAF_C],
										s->value[MAX_FROOT_C],
										s->value[CROOT_C],
										s->value[CROOT_LIVEWOOD_C],
										s->value[CROOT_DEADWOOD_C],
										(s->value[CROOT_LIVEWOOD_C] * 100. ) / s->value[CROOT_C],
										s->value[BRANCH_C],
										s->value[BRANCH_LIVEWOOD_C],
										s->value[BRANCH_DEADWOOD_C],
										(s->value[BRANCH_LIVEWOOD_C] * 100. ) / s->value[BRANCH_C],
										s->value[FRUIT_C],
										s->value[MAX_FRUIT_C],
										s->value[TOT_SAPWOOD_C],
										s->value[TOT_HEARTWOOD_C],
										s->value[RESERVE_N],
										s->value[STEM_N],
										s->value[STEM_LIVEWOOD_N],
										s->value[STEM_DEADWOOD_N],
										s->value[CROOT_N],
										s->value[CROOT_LIVEWOOD_N],
										s->value[CROOT_DEADWOOD_N],
										s->value[BRANCH_N],
										s->value[BRANCH_LIVEWOOD_N],
										s->value[BRANCH_DEADWOOD_N],
										s->value[FRUIT_N],
										s->value[STANDING_WOOD],
										s->value[YEARLY_C_TO_WOOD],
										s->value[CUM_YEARLY_C_TO_WOOD],
										s->value[STAND_BASAL_AREA_m2],
										s->value[TREE_CAI],
										s->value[TREE_MAI],
										s->value[CAI],
										s->value[MAI],
										s->value[VOLUME],
										s->value[TREE_VOLUME],
										(s->value[TREE_CAI]*100.)/s->value[TREE_VOLUME],
										s->value[DELTA_AGB],
										s->value[DELTA_BGB],
										s->value[AGB],
										s->value[BGB],
										s->value[BGB]/s->value[AGB],
										s->value[DELTA_TREE_AGB],
										s->value[DELTA_TREE_BGB],
										s->value[C_HWP],
										s->value[VOLUME_HWP],
										s->value[YEARLY_STEM_AUT_RESP],
										s->value[YEARLY_LEAF_AUT_RESP],
										s->value[YEARLY_FROOT_AUT_RESP],
										s->value[YEARLY_CROOT_AUT_RESP],
										s->value[YEARLY_BRANCH_AUT_RESP]);
										
							 /************************************************************************/
	      /* printing variables at cell level only if there's more than one layer */  //Note: in the 5p6 we write the data anyhow
           
              if (!print_cell)
                   {
                        logger(g_annual_log, ",%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f"
			",%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.2f",
			c->annual_gpp,
			c->annual_npp,
			c->annual_aut_resp,
			c->annual_het_resp,
			c->annual_soil_resp,
			c->annual_soil_respCO2,
			c->annual_r_eco,
			c->annual_nee,
			c->annual_nep,
			c->annual_et,
			c->annual_lh_flux,
			c->annual_soil_evapo,
                        c->years[year].yearly_mean.asw,
			//c->asw,   //  comment: this is actually the value at the 31 december, it has more sense to have the mean value
			c->annual_iwue,
			c->volume,
			c->cum_volume,
			c->annual_out_flow,
			c->litrC,
			c->litr1C,
			c->litr2C,
			c->litr3C,
			c->litr4C,
			c->cwd_C,
			c->cwd_2C,
			c->cwd_3C,
			c->cwd_4C,
			c->soilC,
			c->soil1C,
			c->soil2C,
			c->soil3C,
			c->soil4C,
			c->litrN,
			c->litr1N,
			c->litr2N,
			c->litr3N,
			c->litr4N,
			c->cwd_N,
			c->cwd_2N,
			c->cwd_3N,
			c->cwd_4N,
			c->soilN,
			c->soil1N,
			c->soil2N,
			c->soil3N,
			c->soil4N,
			c->soil_depth);
	            /************************************************************************/
	            /* print meteo variables at cell level */
	           logger(g_annual_log, ",%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f",  
			c->years[year].yearly_mean.solar_rad     ,
			c->years[year].yearly_mean.tavg          ,
			c->years[year].yearly_mean.tmax          ,
			c->years[year].yearly_mean.tmin          ,
			c->years[year].yearly_mean.tday          ,
			c->years[year].yearly_mean.tnight        ,
			c->years[year].yearly_mean.vpd           ,
			c->years[year].yearly_mean.prcp          ,
			c->years[year].yearly_mean.tsoil         ,
			c->years[year].yearly_mean.rh_f          ,
		//	c->years[year].yearly_mean.asw           ,
			c->years[year].co2Conc);
	            /************************************************************************/

                   }
                   
                   /* end print */
	           logger(g_annual_log,"\n");
										
							
							}  // specie
							
						}  // age
						
					 }   // dbh
				
				}  // check on the layer
				
			}   // height
			

		} // layer

	//}    
	//else   // no height_count
	//{
		//TODO ALESSIOC TO ALLESSIOR PRINT EMPTY SPACES WHEN N_TREE = 0

	//}
}

void EOY_cell_msg(void)
{
	if  ( g_annual_log )
	{
		g_annual_log->std_output = 1;
		logger(g_annual_log, sz_launched, netcdf_get_version(), datetime_current());
		print_model_paths(g_annual_log);
		//const char* p;
		//p = file_get_name_only(g_annual_log->filename);
		logger(g_annual_log, "#output file = %s\n", g_annual_log->filename);
		print_model_settings(g_annual_log);
	}
}

/*************************************************************************************************************************************/
/*************************************************************************************************************************************/

void EOD_print_output_soil_cell_level(cell_t *const c, const int day, const int month, const int year, const int years_of_simulation )
{
	static int print_header = 0;

	/* return if monthly logging is off*/
	if ( ! g_daily_soil_log ) return;

	/* heading */
	if ( ! print_header )
	{
		print_header = 1;

		logger(g_daily_soil_log, "X,Y,YEAR,MONTH,DAY");

		/************************************************************************/
		/* heading variables at cell level only if there's more than one layer */

		/* carbon pools */
		logger(g_daily_soil_log,
				",***,"
				"leaf_litr1C,"
				"leaf_litr2C,"
				"leaf_litr3C,"
				"leaf_litr4C,"
				"froot_litr1C,"
				"froot_litr2C,"
				"froot_litr3C,"
				"froot_litr4C,"
				"litrC,"
				"litr1C,"
				"litr2C,"
				"litr3C,"
				"litr4C,"
				"soilC,"
				"soil1C,"
				"soil2C,"
				"soil3C,"
				"soil4C");

		/* nitrogen pools */
		logger(g_daily_soil_log,
				",***,"
				"leaf_litr1N,"
				"leaf_litr2N,"
				"leaf_litr3N,"
				"leaf_litr4N,"
				"froot_litr1N,"
				"froot_litr2N,"
				"froot_litr3N,"
				"froot_litr4N,"
				"litrN,"
				"litr1N,"
				"litr2N,"
				"litr3N,"
				"litr4N,"
				"soilN,"
				"soil1N,"
				"soil2N,"
				"soil3N,"
				"soil4N\n");
	}

	/* values */
	logger(g_daily_soil_log, "%d,%d,%d,%d,%d", c->x, c->y, c->years[year].year, month + 1, day + 1);

	/************************************************************************/

	/* carbon pools */
	logger(g_daily_soil_log, ",***,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f",
			c->leaf_litr1C,
			c->leaf_litr2C,
			c->leaf_litr3C,
			c->leaf_litr4C,
			c->froot_litr1C,
			c->froot_litr2C,
			c->froot_litr3C,
			c->froot_litr4C,
			c->litrC,
			c->litr1C,
			c->litr2C,
			c->litr3C,
			c->litr4C,
			c->soilC,
			c->soil1C,
			c->soil2C,
			c->soil3C,
			c->soil4C);

	/* nitrogen pools */
	logger(g_daily_soil_log, ",***,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f\n",
			c->leaf_litr1N,
			c->leaf_litr2N,
			c->leaf_litr3N,
			c->leaf_litr4N,
			c->froot_litr1N,
			c->froot_litr2N,
			c->froot_litr3N,
			c->froot_litr4N,
			c->litrN,
			c->litr1N,
			c->litr2N,
			c->litr3N,
			c->litr4N,
			c->soilN,
			c->soil1N,
			c->soil2N,
			c->soil3N,
			c->soil4N);
}

void EOD_soil_msg(void)
{
	if ( g_daily_soil_log )
	{
		g_daily_soil_log->std_output = 1;
		logger(g_daily_soil_log, sz_launched, netcdf_get_version(), datetime_current());
		print_model_paths(g_daily_soil_log);
		//const char* p;
		//p = file_get_name_only(g_daily_soil_log->filename);
		logger(g_daily_soil_log, "#output file = %s\n", g_daily_soil_log->filename);
		print_model_settings(g_daily_soil_log);
	}
}

void EOM_print_output_soil_cell_level(cell_t *const c, const int month, const int year, const int years_of_simulation )
{
	static int print_header = 0;

	/* return if monthly logging is off*/
	if ( ! g_monthly_soil_log ) return;

	/* heading */
	if ( ! print_header )
	{
		print_header = 1;

		logger(g_monthly_soil_log, "X,Y,YEAR,MONTH");

		/************************************************************************/
		/* heading variables at cell level only if there's more than one layer */

		/* carbon pools */
		logger(g_monthly_soil_log,
				",***,"
				"leaf_litr1C,"
				"leaf_litr2C,"
				"leaf_litr3C,"
				"leaf_litr4C,"
				"froot_litr1C,"
				"froot_litr2C,"
				"froot_litr3C,"
				"froot_litr4C,"
				"litrC,"
				"litr1C,"
				"litr2C,"
				"litr3C,"
				"litr4C,"
				"soilC,"
				"soil1C,"
				"soil2C,"
				"soil3C,"
				"soil4C");

		/* nitrogen pools */
		logger(g_monthly_soil_log,
				",***,"
				"leaf_litr1N,"
				"leaf_litr2N,"
				"leaf_litr3N,"
				"leaf_litr4N,"
				"froot_litr1N,"
				"froot_litr2N,"
				"froot_litr3N,"
				"froot_litr4N,"
				"litrN,"
				"litr1N,"
				"litr2N,"
				"litr3N,"
				"litr4N,"
				"soilN,"
				"soil1N,"
				"soil2N,"
				"soil3N,"
				"soil4N\n");
	}

	/* values */
	logger (g_monthly_soil_log, "%d,%d,%d,%d", c->x, c->y, c->years[year].year, month +1 );

	/************************************************************************/

	/* carbon pools */
	logger(g_monthly_soil_log, ",***,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f",
			c->leaf_litr1C,
			c->leaf_litr2C,
			c->leaf_litr3C,
			c->leaf_litr4C,
			c->froot_litr1C,
			c->froot_litr2C,
			c->froot_litr3C,
			c->froot_litr4C,
			c->litrC,
			c->litr1C,
			c->litr2C,
			c->litr3C,
			c->litr4C,
			c->soilC,
			c->soil1C,
			c->soil2C,
			c->soil3C,
			c->soil4C);

	/* nitrogen pools */
	logger(g_monthly_soil_log, ",***,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f\n",
			c->leaf_litr1N,
			c->leaf_litr2N,
			c->leaf_litr3N,
			c->leaf_litr4N,
			c->froot_litr1N,
			c->froot_litr2N,
			c->froot_litr3N,
			c->froot_litr4N,
			c->litrN,
			c->litr1N,
			c->litr2N,
			c->litr3N,
			c->litr4N,
			c->soilN,
			c->soil1N,
			c->soil2N,
			c->soil3N,
			c->soil4N);
}

void EOM_soil_msg(void)
{
	if ( g_monthly_soil_log )
	{
		g_monthly_soil_log->std_output = 1;
		logger(g_monthly_soil_log, sz_launched, netcdf_get_version(), datetime_current());
		print_model_paths(g_monthly_soil_log);
		//const char* p;
		//p = file_get_name_only(g_monthly_soil_log->filename);
		logger(g_monthly_soil_log, "#output file = %s\n", g_monthly_soil_log->filename);
		print_model_settings(g_monthly_soil_log);
	}
}

void EOY_print_output_soil_cell_level(cell_t *const c, const int year, const int years_of_simulation )
{
	static int print_header = 0;

	/* return if annual logging is off*/
	if ( ! g_annual_soil_log ) return;

	/* heading */
	if ( ! print_header )
	{
		print_header = 1;

		logger(g_annual_soil_log, "X,Y,YEAR");

		/************************************************************************/
		/* heading variables at cell level only if there's more than one layer */

		/* carbon pools */
		logger(g_annual_soil_log,
				",***,"
				"leaf_litr1C,"
				"leaf_litr2C,"
				"leaf_litr3C,"
				"leaf_litr4C,"
				"froot_litr1C,"
				"froot_litr2C,"
				"froot_litr3C,"
				"froot_litr4C,"
				"litrC,"
				"litr1C,"
				"litr2C,"
				"litr3C,"
				"litr4C,"
				"soilC,"
				"soil1C,"
				"soil2C,"
				"soil3C,"
				"soil4C");

		/* nitrogen pools */
		logger(g_annual_soil_log,
				",***,"
				"leaf_litr1N,"
				"leaf_litr2N,"
				"leaf_litr3N,"
				"leaf_litr4N,"
				"froot_litr1N,"
				"froot_litr2N,"
				"froot_litr3N,"
				"froot_litr4N,"
				"litrN,"
				"litr1N,"
				"litr2N,"
				"litr3N,"
				"litr4N,"
				"soilN,"
				"soil1N,"
				"soil2N,"
				"soil3N,"
				"soil4N\n");
	}
	/*****************************************************************************************************/

	/* values */
	logger(g_annual_soil_log, "%d,%d,%d", c->x, c->y, c->years[year].year);

	/************************************************************************/

	/* carbon pools */
	logger(g_annual_soil_log, ",***,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f",
			c->leaf_litr1C,
			c->leaf_litr2C,
			c->leaf_litr3C,
			c->leaf_litr4C,
			c->froot_litr1C,
			c->froot_litr2C,
			c->froot_litr3C,
			c->froot_litr4C,
			c->litrC,
			c->litr1C,
			c->litr2C,
			c->litr3C,
			c->litr4C,
			c->soilC,
			c->soil1C,
			c->soil2C,
			c->soil3C,
			c->soil4C);
	/* nitrogen pools */
	logger(g_annual_soil_log, ",***,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f,%3.4f\n",
			c->leaf_litr1N,
			c->leaf_litr2N,
			c->leaf_litr3N,
			c->leaf_litr4N,
			c->froot_litr1N,
			c->froot_litr2N,
			c->froot_litr3N,
			c->froot_litr4N,
			c->litrN,
			c->litr1N,
			c->litr2N,
			c->litr3N,
			c->litr4N,
			c->soilN,
			c->soil1N,
			c->soil2N,
			c->soil3N,
			c->soil4N);
}

void EOY_soil_msg(void)
{
	if ( g_annual_soil_log )
	{
		g_annual_soil_log->std_output = 1;
		logger(g_annual_soil_log, sz_launched, netcdf_get_version(), datetime_current());
		print_model_paths(g_annual_soil_log);
		//const char* p;
		//p = file_get_name_only(g_annual_soil_log->filename);
		logger(g_annual_log, "#output file = %s\n", g_annual_soil_log->filename);
		print_model_settings(g_annual_soil_log);
	}
}

#endif
