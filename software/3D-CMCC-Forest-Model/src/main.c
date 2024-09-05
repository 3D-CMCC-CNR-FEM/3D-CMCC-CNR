/* main.c */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#ifdef _WIN32
#ifdef _DEBUG
#define _CRTDBG_MAP_ALLOC
//#define  NC_USE  // to define if NC related functions have to be used (and if not, allowing to build an exe for windows). To use ncdf library,just un-comment!
#include <crtdbg.h>
#endif
#endif
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <errno.h>
#include <stdarg.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "compiler.h"
#include "soil_settings.h"
#include "print_output.h"
#include "common.h"
#include "constants.h"
#include "topo.h"
#include "settings.h"
#include "logger.h"
#include "output.h"
#include "matrix.h"
#include "meteo.h"
#include "nc.h"
#include "met_data.h"
#include "cropmodel_daily.h"
#include "g-function.h"
#include "leaf_fall.h"
#include "tree_model.h"
#include "soil_dndc.h"
#include "fluxes.h"
#include "print.h"
#include "utility.h"
#include "cell_model.h"
#include "soil_model.h"
#include "litter_model.h"
#include "compare.h"
#include "management.h"
#include "mortality.h"
#include "soil_water_balance.h"

//#define BENCHMARK_ONLY

/* Last cumulative days in months in non Leap years */
int MonthLength [] = { 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365};

/* Last cumulative days in months in Leap years */
int MonthLength_Leap [] = { 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366};

/* Days in Months */
int DaysInMonth [] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };

const char *szMonth[MONTHS_COUNT] = { "JANUARY", "FEBRUARY", "MARCH", "APRIL", "MAY", "JUNE", "JULY",
		"AUGUST", "SEPTEMBER", "OCTOBER", "NOVEMBER", "DECEMBER" };

/* do not change order in logger.h */
const char* log_types[LOG_TYPES_COUNT] = { "debug" , "daily" , "monthly" , "annual" , "soil_daily", "soil_monthly", "soil_annual" };

/* global vars */
soil_settings_t* g_soil_settings = NULL;
topo_t* g_topo = NULL;
settings_t* g_settings = NULL;
dataset_t* g_dataset = NULL;
management_t* g_management = NULL;

/* DO NOT REMOVE INITIALIZATION TO NULL, IT IS REQUIRED !! */
char	*g_sz_parameterization_path = NULL
		, *g_sz_input_path = NULL
		, *g_sz_output_path = NULL
		, *g_sz_dataset_file = NULL
		, *g_sz_input_met_file = NULL
		, *g_sz_soil_file = NULL
		, *g_sz_topo_file = NULL
		, *g_sz_settings_file =	NULL
                , *g_sz_manag_file =	NULL
		, *g_sz_ndep_file = NULL
		, *g_sz_co2_conc_file = NULL
		, *g_sz_output_vars_file = NULL
		, *g_sz_benchmark_path = NULL
		;

int g_year_start_index;
char g_sz_input_data_path[256];
char g_sz_output_fullpath[256];

static int years_of_simulation;	// default is none

/* extern variables */
extern logger_t* g_debug_log;
extern logger_t* g_daily_log;
extern logger_t* g_monthly_log;
extern logger_t* g_annual_log;
extern logger_t* g_daily_soil_log;
extern logger_t* g_monthly_soil_log;
extern logger_t* g_annual_soil_log;

/* strings */
const char sz_launched[] = "\n#--------------------------------------------------------------------------------\n"
		"#"PROGRAM_FULL_NAME"\n"
		"#compiled using "COMPILER" on "__DATE__" at "__TIME__"\n"
		"#using NetCDF %s\n"
		"#launched: %s\n"
		"#--------------------------------------------------------------------------------\n";

static const char banner[] = "\n#"PROGRAM_FULL_NAME"\n"
		"#compiled using "COMPILER" on "__DATE__" at "__TIME__"\n"
		"#using NetCDF %s\n"
		/*
#if defined _WIN32
		"previous commit hash: "COMMIT_HASH"\n"
#endif
		 */
		"(use -h parameter for more information)\n\n"
		"The "PROGRAM_FULL_NAME" has been developed by:\n"
		"Alessio Collalti [alessio.collalti@cnr.it, alessio.collalti@isafom.cnr.it],\n"
		"Alessio Ribeca [alessio.ribeca@cmcc.it]\n"
		"Carlo Trotta [trottacarlo@unitus.it]\n"
		"Corrado Biondo [@CMCC]\n"
		"Gina Marano [gina.marano@unina.it]\n"
		"Giorgio Matteucci [giorgio.matteucci@isafom.cnr.it]\n"
		"Daniela Dalmonech [daniela.dalmonech@isafom.cnr.it]\n"
		"National Research Council of Italy (CNR),\n"
		"Institute for Agricultural and Forestry Systems in the Mediterranean(ISAFOM),\n"
		"Via della Madonna Alta, 128, 06128 - Perugia (PG), Italy\n"
		"\n"
		"and \n"
		"Tuscia University (UNITUS),\n"
		"Department for innovation in biological, agro-food and forest systems (DIBAF),\n"
		"Forest Ecology Lab\n"
		"Programmers: Alessio Collalti - Alessio Ribeca - Carlo Trotta  \n"
		"\n"
		"\"DISCLAIMER\"\n"
		"CNR-UNITUS\n"
		"accepts no responsibility for the use of the 3D-CMCC-CNR FEM in\n"
		"the form supplied or as subsequently modified by third parties.\n"
		"CNR and UNITUS disclaims liability for all losses,\n"
		"damages and costs incurred by any person as a result of relying on this software.\n"
		"Use of this software assumes agreement to this condition of use.\n"
		"Removal of this statement violates the spirit in which 3D-CMCC-FEM,\n"
		"was released by CNR-UNITUS.\n"
		"for more information see:\n"
		"-Collalti et al., 2014 Ecological Modelling,\n"
		"-Collalti et al., 2016 Geoscientific Model Development\n"
		"-Marconi  et al., 2017 Forests\n"
		"-Collalti et al., 2017 Forest@\n"
		"-Collalti et al., 2018 Journal of Advances in Earth System Modeling\n"
		"-Collalti et al., 2019a Ecological Applications\n"
		"-Collalti et al., 2019b Global Change Biology\n"
		"--------------------------------------------------------------------------------\n";
static const char msg_input_path[]					=	"input path = %s\n";
static const char msg_parameterization_path[]		=	"parameterization path = %s\n";
static const char msg_soil_file[]					=	"soil file = %s\n";
static const char msg_topo_file[]					=	"topo file = %s\n";
static const char msg_met_path[]					=	"met path = %s\n";
static const char msg_settings_path[]				=	"settings path = %s\n";
static const char msg_co2_conc_file[]				=	"co2 conc file path = %s\n";
static const char msg_ndep_file[]					=	"ndep file path = %s\n";
static const char msg_debug_output_file[]			=	"debug output file path = %s\n";
static const char msg_daily_output_file[]			=	"daily output file path = %s\n";
static const char msg_monthly_output_file[]			=	"monthly output file path = %s\n";
static const char msg_annual_output_file[]			=	"annual output file path = %s\n";
static const char msg_soil_daily_output_file[]		=	"soil daily output file path = %s\n";
static const char msg_soil_monthly_output_file[]	=	"soil monthly output file path = %s\n";
static const char msg_soil_annual_output_file[]		=	"soil annual output file path = %s\n";
static const char msg_ok[]							=	"ok ";
static const char msg_usage[]						=	"\nusage:\n"
		"  3D-CMCC-FEM -i INPUT_DIR -d DATASET_FILENAME -m MET_FILE_LIST -s SITE_FILENAME -c SETTINGS_FILENAME [-o OUTPUT_FILENAME] [-h]\n"
		"  required options:\n"
#ifdef _WIN32
		"    -i input path (i.e.: -i c:\\input\\directory\\)\n"
		"    -o output path (i.e.: -o c:\\output\\directory\\)\n"
		"    -p parameterization directory (i.e.: -i c:\\parameterization\\directory\\)\n"
#else
		"    -i input directory (i.e.: -i path/to/input/folder/)\n"
		"    -o output directory (i.e.: -i path/to/output/folder/)\n"
		"    -p parameterization directory (i.e.: -i path/to/parameterization/directory/)\n"
#endif
		"    -d dataset filename stored into input directory (i.e.: -d input.txt)\n"
		"    -m met filename list stored into input directory (i.e.: -m 1999.txt,2003.txt,2009.txt)\n"
		"    -s soil filename stored into input directory (i.e.: -s soil.txt or soil.nc)\n"
		"    -t topo filename stored into input directory (i.e.: -t topo.txt or topo.nc)\n"
		"    -c settings filename stored into input directory (i.e.: -c settings.txt)\n"
		"    -k co2 concentration file (i.e.: -k co2_conc.txt)\n"
		"    -q management scheme file stored into input directory (i.e.: -q management.txt)\n"
		"    -n ndep file (i.e.: -n ndep.txt)\n"
		"    -u benchmark path\n"   // not used in the current version
		"    -h print this help\n"
		;

// COMMENT: 5p6 Version used for regional simulations in Basilicata

/* error messages */
extern const char sz_err_out_of_memory[];
const char err_unable_open_file[] = "unable to open file.";
const char err_empty_file[] = "empty file ?";
const char err_window_size_too_big[] = "window size too big.";
//static const char err_unable_get_current_path[] = "unable to retrieve current path.\n";
static const char err_unable_to_register_atexit[] = "unable to register clean-up routine.\n";

/*
static int cells_sort(const void *a, const void *b)
{
	if ( ((cell_t*)a)->elev < ((cell_t*)b)->elev )
	{
		return -1;
	} else if ( ((cell_t*)a)->elev > ((cell_t*)b)->elev )
	{
		return 1;
	} else
	{
		return 0;
	}
}
 */

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

static void clean_up(void)
{
	if ( g_sz_benchmark_path ) free(g_sz_benchmark_path);
	if ( g_sz_output_vars_file ) free(g_sz_output_vars_file);
	if ( g_sz_ndep_file ) free(g_sz_ndep_file);
	if ( g_sz_co2_conc_file ) free(g_sz_co2_conc_file);
	if ( g_sz_settings_file ) free(g_sz_settings_file);
        if ( g_sz_manag_file ) free(g_sz_manag_file);
	if ( g_sz_parameterization_path ) free(g_sz_parameterization_path);
	if ( g_sz_output_path ) free(g_sz_output_path);
	if ( g_sz_input_path ) free(g_sz_input_path);
	//if ( g_topo ) free(g_topo);
	if ( g_sz_topo_file ) free(g_sz_topo_file);
	//if ( g_soil_settings ) free(g_soil_settings);
	if ( g_sz_soil_file ) free(g_sz_soil_file);
	if ( g_management ) management_free(g_management);
	if ( g_dataset) dataset_free(g_dataset);
	if ( g_settings ) settings_free(g_settings);
	if ( g_sz_input_met_file ) free(g_sz_input_met_file);
	if ( g_sz_dataset_file ) free(g_sz_dataset_file);

#ifdef _WIN32
#ifdef _DEBUG
	_CrtSetReportMode(_CRT_WARN, _CRTDBG_MODE_FILE);
	_CrtSetReportFile(_CRT_WARN, _CRTDBG_FILE_STDOUT);
	_CrtSetReportMode(_CRT_ERROR, _CRTDBG_MODE_FILE);
	_CrtSetReportFile(_CRT_ERROR, _CRTDBG_FILE_STDOUT);
	_CrtSetReportMode(_CRT_ASSERT, _CRTDBG_MODE_FILE);
	_CrtSetReportFile(_CRT_ASSERT, _CRTDBG_FILE_STDOUT);
	_CrtDumpMemoryLeaks();

	system("PAUSE");
#endif
#endif
}

static void show_usage(void)
{
	puts(msg_usage);
}

static int output_path_create(void)
{
#define BUFFER_SIZE	1024
	char date[32];
	char buf[BUFFER_SIZE+1];
	int i;

	/* get current date */
	{
		struct tm* ptm;
		time_t t;

		/* get local datetime */
		time(&t);
		ptm = localtime(&t);
		sprintf(date, "%04d_%s_%02d"
				, ptm->tm_year + 1900
				, szMonth[ptm->tm_mon]
						  , ptm->tm_mday
		);
	}

	if ( ! g_sz_output_path ) {
		g_sz_output_path = string_copy("output");
		if ( ! g_sz_output_path ) {
			puts(sz_err_out_of_memory);
			return 0;
		}
	}

	i = has_path_delimiter(g_sz_output_path);
	i = sprintf(buf, "%s%soutput_%s_%s%s"
			, g_sz_output_path
			, i ? "" : FOLDER_DELIMITER
					, PROGRAM_VERSION
					, date
					, FOLDER_DELIMITER
	);
	if ( i < 0 ) {
		puts("unable to format output path!");
		return 0;
	}
	free(g_sz_output_path);
	g_sz_output_path = string_copy(buf);
	if ( ! g_sz_output_path ) {
		puts(sz_err_out_of_memory);
		return 0;
	}

	return 1;
#undef BUFFER_SIZE
}

static int parameterization_output_create(void) {
	int i;

	i = sprintf(g_sz_input_data_path
			, "%sinput_data%s"
			, g_sz_output_path
			, FOLDER_DELIMITER
	);
	if ( i < 0 )
	{
		return 0;
	}

	return path_create(g_sz_input_data_path);
}

static int log_start(const char* const sitename)
{
	char *p;
	char buffer[512];	/* should be enough */
	int len;

	/* build log filename */

	/* program version */
	len = sprintf(buffer, "_%s", PROGRAM_VERSION);

	/* sitename */
	if ( sitename && sitename[0] ) {
		len += sprintf(buffer+len, "_%s", sitename);
	}

	/* met file */
	len += sprintf(buffer+len, "_%s", get_filename(g_sz_input_met_file));

	/* co2_fixed */
	len += sprintf(buffer+len, "_(%d", g_settings->year_start);

	/* end simulation */
	len += sprintf(buffer+len, "-%d)", g_settings->year_end);

	/* co2_transient */
	switch ( g_settings->CO2_trans ) {
	case CO2_TRANS_ON:
		p = "ON";
		break;

	case CO2_TRANS_OFF:
		p = "OFF";
		break;

	case CO2_TRANS_VAR:
		p = "VAR";
		break;

	default:
		puts("bad CO2_trans value in settings file!\n");
		return 0;
	}
	len += sprintf(buffer+len, "_CO2_%s", p);

	/* CO2 file */
	len += sprintf(buffer+len, "_%s", get_filename(g_sz_co2_conc_file));

	/* management */
	switch ( g_settings->management ) {
	case MANAGEMENT_ON:
		p = "ON";
		break;

	case MANAGEMENT_OFF:
		p = "OFF";
		break;

	case MANAGEMENT_VAR:
		p = "VAR";
		break;

        case MANAGEMENT_VAR1:   //in case an additional management type will be used.
		p = "VAR1";
		break;

	default:
		puts("bad management value in settings file!\n");
		return 0;
	}
	len += sprintf(buffer+len, "_Man_%s", p);

	/* time */
	len += sprintf(buffer+len, "_%c", g_settings->time);

	/* sizeCell */
	len += sprintf(buffer+len, "_%d", (int)g_settings->sizeCell);

	/* date */
	//len += sprintf(buffer+len, "_%s", date);

	/* files type */
	{
		char *p;
		char *ext;
		p = strrchr(g_sz_input_met_file, '.');
		if ( ! p ) {
			puts("bad met file!\n");
			return -1;
		}
		++p;
		if ( ! string_compare_i(p, "nc") || ! string_compare_i(p, "nc4") ) {
			ext = "_nc";
		} else if ( ! string_compare_i(p, "lst") ) {
			ext = "_lst";
		} else if ( ! string_compare_i(p, "txt") ) {
			ext = "_txt";
		} else {
			puts("bad met file!\n");
			return 0;
		}
		len += sprintf(buffer+len, "%s", ext);
	}

	/* extension */
	sprintf(buffer+len, ".txt");

	/* create log files and input_data folder */
	{
		int i;
		int log_flag[LOG_TYPES_COUNT];
		logger_t** logs[LOG_TYPES_COUNT] = {
				&g_debug_log
				, &g_daily_log
				, &g_monthly_log
				, &g_annual_log
				, &g_daily_soil_log
				, &g_monthly_soil_log
				, &g_annual_soil_log
		};

		log_flag[0] = g_settings->debug_output;
		log_flag[1] = g_settings->daily_output;
		log_flag[2] = g_settings->monthly_output;
		log_flag[3] = g_settings->yearly_output;
		log_flag[4] = g_settings->soil_output;
		log_flag[5] = g_settings->soil_output;
		log_flag[6] = g_settings->soil_output;

		for ( i = 0 ; i < LOG_TYPES_COUNT; ++i ) {
			if ( log_flag[i] ) {
				*logs[i] = logger_new("%s%s%s%s%s"
						, g_sz_output_path
						, log_types[i]
									, FOLDER_DELIMITER
									, log_types[i]
												, buffer
				);
				if ( ! *logs[i] ) {
					printf("Unable to create %s log!\n\n", log_types[i]);
					return 0;
				}
			}
		}
	}

	/* disable screen output when "off" */
	if ( ! g_settings->screen_output ) { 
		logger_disable_std(g_debug_log);
	}

	logger(g_debug_log, sz_launched, netcdf_get_version(), datetime_current());

	/* show paths */
	logger(g_debug_log,"\nFILE NAMES AND PATHS\n");
	if ( g_sz_input_path )
	{
		logger(g_debug_log, msg_input_path, g_sz_input_path);
	}
	logger(g_debug_log, msg_parameterization_path, g_sz_parameterization_path);
	logger(g_debug_log, msg_soil_file, g_sz_soil_file);
	logger(g_debug_log, msg_topo_file, g_sz_topo_file);
	logger(g_debug_log, msg_met_path, g_sz_input_met_file);

	if ( g_sz_co2_conc_file ) 
		logger(g_debug_log, msg_co2_conc_file, g_sz_co2_conc_file);

	if ( g_sz_ndep_file ) 
		logger(g_debug_log, msg_ndep_file, g_sz_ndep_file);

	logger(g_debug_log, msg_settings_path, g_sz_settings_file);
	if ( g_debug_log )						logger(g_debug_log, msg_debug_output_file, g_debug_log->filename);
	if ( g_debug_log && g_daily_log )		logger(g_debug_log, msg_daily_output_file, g_daily_log->filename);
	if ( g_debug_log && g_monthly_log )		logger(g_debug_log, msg_monthly_output_file, g_monthly_log->filename);
	if ( g_debug_log && g_annual_log )		logger(g_debug_log, msg_annual_output_file, g_annual_log->filename);
	if ( g_debug_log && g_daily_soil_log )		logger(g_debug_log, msg_soil_daily_output_file, g_daily_soil_log->filename);
	if ( g_debug_log && g_monthly_soil_log )	logger(g_debug_log, msg_soil_monthly_output_file, g_monthly_soil_log->filename);
	if ( g_debug_log && g_annual_soil_log )		logger(g_debug_log, msg_soil_annual_output_file, g_annual_soil_log->filename);

	if ( g_daily_log ) g_daily_log->std_output = 0;
	if ( g_monthly_log ) g_monthly_log->std_output = 0;
	if ( g_annual_log ) g_annual_log->std_output = 0;
	if ( g_daily_soil_log ) g_daily_soil_log->std_output = 0;
	if ( g_monthly_soil_log ) g_monthly_soil_log->std_output = 0;
	if ( g_annual_soil_log ) g_annual_soil_log->std_output = 0;

	return 1;
}

/*
	copy path and add a / at end if missing
 */

char* path_copy(const char *const s) {
	if ( s ) {
		int i;
		i = strlen(s);
		if ( ('/' == s[i-1]) || ('\\' == s[i-1]) ) {
			return string_copy(s);
		} else {
			char *p = malloc(i+1+1);
			if ( p ) {
				if ( strcpy(p, s) ) {
					strcat(p, "/");
					return p;
				} else {
					free(p);
					p = NULL;
				}
			}
		}
	}
	return NULL;
}

/*
	parse and check passed args
 */
static int parse_args(int argc, char *argv[])
{
	char* p;
	int i;

	g_sz_input_path = NULL;
	g_sz_output_path = NULL;
	g_sz_parameterization_path = NULL;
	g_sz_dataset_file = NULL;
	g_sz_input_met_file = NULL;
	g_sz_soil_file = NULL;
	g_sz_topo_file = NULL;
	g_sz_settings_file = NULL;
        g_sz_manag_file    = NULL; 
        g_sz_input_met_file = NULL;
	g_sz_output_vars_file = NULL;
	g_sz_benchmark_path = NULL;

           printf("check files");

	for ( i = 1; i < argc; ++i ) {
		if ( argv[i][0] != '-' ) {
			continue;
		}

		switch ( argv[i][1] ) {
		case 'i': /* folder where input files are stored */
			if ( ! argv[i+1] ) {
				puts("input path not specified!");
				goto err;
			}
			g_sz_input_path = path_copy(argv[i+1]);
			if ( ! g_sz_input_path ) {
				puts(sz_err_out_of_memory);
				goto err;
			}
			break;

		case 'p': /* folder where parameterization files are stored */
			if ( ! argv[i+1] ) {
				puts("parameterization path not specified!");
				goto err;
			}
			g_sz_parameterization_path = path_copy(argv[i+1]);
			if ( ! g_sz_parameterization_path ) {
				puts(sz_err_out_of_memory);
				goto err;
			}
			break;

		case 'o': /* output folder*/
			if ( ! argv[i+1] ) {
				puts("output folder not specified");
				goto err;
			}
			g_sz_output_path = string_copy(argv[i+1]);
			if ( ! g_sz_output_path ) {
				puts(sz_err_out_of_memory);
				goto err;
			}
			break;

		case 'd': /* dataset file path */
			if ( ! argv[i+1] ) {
				puts("dataset file not specified!");
				goto err;
			}
			g_sz_dataset_file = string_copy(argv[i+1]);
			if( ! g_sz_dataset_file ) {
				puts(sz_err_out_of_memory);
				goto err;
			}
			break;

		case 'm': /* met filename */
			if ( ! argv[i+1] ) {
				puts("met file not specified!");
				goto err;
			}
			g_sz_input_met_file = string_copy(argv[i+1]);
			if ( ! g_sz_input_met_file ) {
				puts(sz_err_out_of_memory);
				goto err;
			}
			break;

		case 's': /*soil file */
			if ( ! argv[i+1] ) {
				puts("soil file not specified!");
				goto err;
			}
			g_sz_soil_file = string_copy(argv[i+1]);
			if ( ! g_sz_soil_file ) {
				puts(sz_err_out_of_memory);
				goto err;
			}
			break;

		case 't': /* topo file */
			if ( ! argv[i+1] ) {
				puts("topo file not specified!");
				goto err;
			}
			g_sz_topo_file = string_copy(argv[i+1]);
			if( ! g_sz_topo_file ) {
				puts(sz_err_out_of_memory);
				goto err;
			}
			break;

		case 'c': /* settings file */
			if ( ! argv[i+1] ) {
				puts("settings file not specified!");
				goto err;
			}
			g_sz_settings_file = string_copy(argv[i+1]);
			if( ! g_sz_settings_file ) {
				puts(sz_err_out_of_memory);
				goto err;
			}
			break;

		case 'k': /* co2 file */
			if ( ! argv[i+1] ) {
				puts("co2 conc file not specified!");
				goto err;
			}
			g_sz_co2_conc_file = string_copy(argv[i+1]);
			if( ! g_sz_co2_conc_file ) {
				puts(sz_err_out_of_memory);
				goto err;
			}
			break;

		case 'n': /* ndep file */
			if ( ! argv[i+1] ) {
				puts("ndep file not specified!");
				goto err;
			}
			g_sz_ndep_file = string_copy(argv[i+1]);
			if( ! g_sz_ndep_file ) {
				puts(sz_err_out_of_memory);
				goto err;
			}
			break;

		case 'r': /* output filename */
			if ( ! argv[i+1] ) {
				puts("output path not specified!");
				goto err;
			}
			g_sz_output_vars_file = string_copy(argv[i+1]);
			if( ! g_sz_output_vars_file ) {
				puts(sz_err_out_of_memory);
				goto err;
			}
			break;

              //  prescribed-management file (e.g. SETTINGS, where settings file is also saved)

                 case 'q': // 
 			if ( ! argv[i+1] ) {
				puts("output path not specified!");
				goto err;
			}
			g_sz_manag_file = string_copy(argv[i+1]);


			if( ! g_sz_manag_file ) {
				puts(sz_err_out_of_memory);
				goto err;
			}
			break;

		case 'u': /* compare path */
			if ( ! argv[i+1] ) {
				puts("compare path not specified!");
				goto err;
			}
			g_sz_benchmark_path = string_copy(argv[i+1]);
			if( ! g_sz_benchmark_path ) {
				puts(sz_err_out_of_memory);
				goto err;
			}
			break;

		case 'h': /* show help */
			goto err_show_usage;
			break;

		default:
			printf("invalid option (%s)!\n", argv[i]);
			goto err_show_usage;
		}
	}

	/* check that each mandatory parameter has been used */

	/* output path specified ? */
	if ( ! g_sz_output_path ) {
		g_sz_output_path = string_copy("output");
		if ( ! g_sz_output_path ) {
			puts(sz_err_out_of_memory);
			goto err;
		}
	}

	if ( ! g_sz_parameterization_path ) {
		puts("parameterization path not specified!");
		goto err_show_usage;
	}

	if ( g_sz_dataset_file && g_sz_input_path ) {
		p = concatenate_path(g_sz_input_path, g_sz_dataset_file); 
		if ( ! p )
		{
			puts(sz_err_out_of_memory);
			goto err;
		}
		free(g_sz_dataset_file);
		g_sz_dataset_file = p;
	}

	if ( ! g_sz_input_met_file ) {
		puts("met file list is missing!");
		goto err_show_usage;
	}  else if ( g_sz_input_path ) {
		p = concatenate_path(g_sz_input_path, g_sz_input_met_file); 
		if ( ! p )
		{
			puts(sz_err_out_of_memory);
			goto err;
		}
		free(g_sz_input_met_file);
		g_sz_input_met_file = p;
	}

	if ( ! g_sz_soil_file ) {
		puts("soil filename not specified!");
		goto err_show_usage;
	} else if ( g_sz_input_path ) {
		p = concatenate_path(g_sz_input_path, g_sz_soil_file); 
		if ( ! p )
		{
			puts(sz_err_out_of_memory);
			goto err;
		}
		free(g_sz_soil_file);
		g_sz_soil_file = p;
	}

	if ( ! g_sz_topo_file ) {
		puts("topo filename option is missing!");
		goto err_show_usage;
	} else if ( g_sz_input_path ) {
		p = concatenate_path(g_sz_input_path, g_sz_topo_file); 
		if ( ! p )
		{
			puts(sz_err_out_of_memory);
			goto err;
		}
		free(g_sz_topo_file);
		g_sz_topo_file = p;
	}

	if ( ! g_sz_settings_file ) {
		puts("settings filename option is missing!");
		goto err_show_usage;
	} else if ( g_sz_input_path ) {
		p = concatenate_path(g_sz_input_path, g_sz_settings_file); 
		if ( ! p )
		{
			puts(sz_err_out_of_memory);
			goto err;
		}
		free(g_sz_settings_file);
		g_sz_settings_file = p;
	}

	if ( g_sz_ndep_file ) {
		p = concatenate_path(g_sz_input_path, g_sz_ndep_file);
		if ( ! p )
		{
			puts(sz_err_out_of_memory);
			goto err;
		}
		free(g_sz_ndep_file);
		g_sz_ndep_file = p;
	}

	if ( g_sz_co2_conc_file ) {
		p = concatenate_path(g_sz_input_path, g_sz_co2_conc_file);
		if ( ! p )
		{
			puts(sz_err_out_of_memory);
			goto err;
		}
		free(g_sz_co2_conc_file);
		g_sz_co2_conc_file = p;
	}

	if ( g_sz_output_vars_file ) {
		p = concatenate_path(g_sz_input_path, g_sz_output_vars_file); 
		if ( ! p )
		{
			puts(sz_err_out_of_memory);
			goto err;
		}
		free(g_sz_output_vars_file);
		g_sz_output_vars_file = p;
	}

	return 1;

	err_show_usage:
	show_usage();

	err:
	return 0;
}

static int check_soil_values(void)
{
	/** check for soil mandatory values **/
	/* soil latitude and longitude */
	if (IS_INVALID_VALUE(g_soil_settings->values[SOIL_LAT])
			|| IS_INVALID_VALUE(g_soil_settings->values[SOIL_LON]))
	{
		logger_error(g_debug_log, "NO SOIL DATA AVAILABLE (Latitude and Longitude)");
		goto err;
	}
	/* soil texture */
	if (IS_INVALID_VALUE(g_soil_settings->values[SOIL_SAND_PERC])
			|| IS_INVALID_VALUE(g_soil_settings->values[SOIL_CLAY_PERC])
			|| IS_INVALID_VALUE(g_soil_settings->values[SOIL_SILT_PERC])
			|| IS_INVALID_VALUE(g_soil_settings->values[SOIL_DEPTH]))
	{
		logger_error(g_debug_log, "NO SOIL DATA AVAILABLE (Percentage of texture)");
		goto err;
	}
	/* soil depth */
	if (IS_INVALID_VALUE(g_soil_settings->values[SOIL_DEPTH]))
	{
		logger_error(g_debug_log, "NO SOIL DATA AVAILABLE (Percentage of texture)");
		goto err;
	}
	/* soil fertility */
	if (IS_INVALID_VALUE(g_soil_settings->values[SOIL_FR])
			|| IS_INVALID_VALUE(g_soil_settings->values[SOIL_FN0])
			|| IS_INVALID_VALUE(g_soil_settings->values[SOIL_FNN])
			|| IS_INVALID_VALUE(g_soil_settings->values[SOIL_M0]))
	{
		logger_error(g_debug_log, "NO SOIL DATA AVAILABLE (Soil fertility value)");
		goto err;
	}
	/* litter carbon  */
	if (IS_INVALID_VALUE(g_soil_settings->values[LITTERC]))
	{
		/* initialize to zero value */
		g_soil_settings->values[LITTERC] = 0.;
		logger_error(g_debug_log, "NO SOIL DATA AVAILABLE (Carbon litter value)\n");
	}
	/* litter nitrogen */
	if (IS_INVALID_VALUE(g_soil_settings->values[LITTERN]))
	{
		/* initialize to zero value */
		g_soil_settings->values[LITTERN] = 0.;
		logger_error(g_debug_log, "NO SOIL DATA AVAILABLE (Nitrogen litter value)\n");
	}
	/* litter coarse woody debris carbon */
	if (IS_INVALID_VALUE(g_soil_settings->values[DEADWOODC]))
	{
		/* initialize to zero value */
		g_soil_settings->values[DEADWOODC] = 0.;
		logger_error(g_debug_log, "NO SOIL DATA AVAILABLE (Coarse woody debris carbon value)\n");
	}
	/* soil carbon */
	if (IS_INVALID_VALUE(g_soil_settings->values[SOILC]))
	{
		/* initialize to zero value */
		g_soil_settings->values[SOILC] = 0.;
		logger_error(g_debug_log, "NO SOIL DATA AVAILABLE (Carbon soil value)\n");
	}
	/* soil nitrogen */
	if (IS_INVALID_VALUE(g_soil_settings->values[SOILN]))
	{
		/* initialize to zero value */
		g_soil_settings->values[SOILN] = 0.;
		logger_error(g_debug_log, "NO SOIL DATA AVAILABLE (Nitrogen soil value)\n");
	}

	return 1;

	err:
	return 0;
}
static int check_topo_values(void)
{
	/** check for topo mandatory values **/
	/* topo elevation */
	if (IS_INVALID_VALUE(g_topo->values[TOPO_ELEV]))
	{
		logger_error(g_debug_log, "NO TOPO DATA AVAILABLE (elevation value)");
		goto err;
	}
	return 1;

	err:
	return 0;
}

static int restart(const matrix_t*const m, int restart_year) {
	logger_flush(g_daily_log);
	logger_flush(g_monthly_log);
	logger_flush(g_annual_log);
	logger_flush(g_daily_soil_log);
	logger_flush(g_monthly_soil_log);
	logger_flush(g_annual_soil_log);

	printf("\nprogram paused due to restart on %d\nhit enter when ready\n", restart_year);
	getchar();

	printf("re-import settings file %s...", g_sz_settings_file);
	settings_free(g_settings);
	g_settings = settings_import(g_sz_settings_file);
	if ( ! g_settings )
	{
		return 0;
	}
	puts(msg_ok);

	printf("re-import species...");
	{
		int cell, height, dbh, age, species;
		for ( cell = 0; cell < m->cells_count; ++cell )
		{
			for (height = 0; height < m->cells[cell].heights_count; ++height)
			{
				for (dbh = 0; dbh < m->cells[cell].heights[height].dbhs_count; ++dbh)
				{
					for (age = 0; age < m->cells[cell].heights[height].dbhs[dbh].ages_count; ++age )
					{
						for ( species = 0; species < m->cells[cell].heights[height].dbhs[dbh].ages[age].species_count; ++species )
						{
							if ( ! fill_species_from_file(&m->cells[cell].heights[height].dbhs[dbh].ages[age].species[species]) )
							{
								return 0;
							}
						}
					}
				}
			}
		}
	}
	puts(msg_ok);

	return 1;
}

static int log_rename(void)
{
	char end_year[4+1] = { 0 };
	int i;
	logger_t* logs[LOG_TYPES_COUNT] =
	{
			g_debug_log
			, g_daily_log
			, g_monthly_log
			, g_annual_log
			, g_daily_soil_log
			, g_monthly_soil_log
			, g_annual_soil_log
	};

	sprintf(end_year, "%d", g_settings->year_end);

	for ( i = 0; i < LOG_TYPES_COUNT; ++i )
	{
		if ( logs[i] && logs[i]->f )
		{
			char* p;
			char* new_filename;
			int ret;

			new_filename = string_copy(logs[i]->filename);
			if ( ! new_filename ) return 0;

			p = strchr(new_filename, '(');
			if ( ! p ) return 0;
			p = strchr(p, '-');
			if ( ! p ) return 0;
			++p;

			p[0] = end_year[0];
			p[1] = end_year[1];
			p[2] = end_year[2];
			p[3] = end_year[3];

			fclose(logs[i]->f);

			if ( file_exists(new_filename) )
			{
				remove(new_filename);
			}

			ret = rename(logs[i]->filename, new_filename);
			if ( -1 == ret )
			{
				free(new_filename);
				return 0;
			}
			free(logs[i]->filename);
			logs[i]->filename = new_filename;		
			logs[i]->f = fopen(logs[i]->filename, "a");
			if ( ! logs[i]->f ) return 0;
		}
	}

	return 1;
}

void sort_all(matrix_t* m)
{
	int cell;
	int i;

	assert(m);

	for ( cell = 0; cell < m->cells_count; ++cell )
	{
		int height;
		for ( height = 0; height < m->cells[cell].heights_count; ++height )
		{
			int dbh;
			for ( dbh = 0; dbh < m->cells[cell].heights[height].dbhs_count; ++dbh )
			{
				int age;
				for ( age = 0; age < m->cells[cell].heights[height].dbhs[dbh].ages_count; ++age )
				{
					/*
					int species;
					for ( species = 0; species < m->cells[cell].heights[height].dbhs[dbh].ages[age].species_count; ++species )
					{
						// re-index species
						m->cells[cell].heights[height].dbhs[dbh].ages[age].species[species].index = species;
					}
					 */
					// sort ages
					qsort(m->cells[cell].heights[height].dbhs[dbh].ages, m->cells[cell].heights[height].dbhs[dbh].ages_count, sizeof(age_t),sort_by_ages_desc);
					// re-index
					for ( i = 0; i < m->cells[cell].heights[height].dbhs[dbh].ages_count; ++i )
					{
						m->cells[cell].heights[height].dbhs[dbh].ages[i].index = i;
					}
				}
				// sort dbhs
				qsort(m->cells[cell].heights[height].dbhs, m->cells[cell].heights[height].dbhs_count, sizeof(dbh_t), sort_by_dbhs_desc);
				// re-index
				for ( i = 0; i < m->cells[cell].heights[height].dbhs_count; ++i )
				{
					m->cells[cell].heights[height].dbhs[i].index = i;
				}
			}
			// sort height
			qsort(m->cells[cell].heights, m->cells[cell].heights_count, sizeof(height_t), sort_by_heights_desc);
			// re-index
			for ( i = 0; i < m->cells[cell].heights_count; ++i )
			{
				m->cells[cell].heights[i].index = i;
			}
		}
	}
}




#if 1
//note: 02/february/2017
// the model runs for all days/months/years and then changes cell
// Currently testes extensively only with one cell
// TODO test with more cells in combination with netcdf files

int main(int argc, char *argv[]) {
	int year;
	int month;
	int day;
	int cell;
	int prog_ret;
	int soil_settings_count;
	int topos_count;
	double timer;
	double start_timer;
	double end_timer;
	matrix_t* matrix;
	soil_settings_t* s;
	topo_t* t;

	//_CrtSetBreakAlloc(); //do not remove, for debugging under msvc++
	assert(MONTHS_COUNT == METEO_MONTHS_COUNT); // just to be sure

	/* initialize */
	matrix = NULL;
	t = NULL;
	s = NULL;
	prog_ret = 1;

	/* start timer */
	timer_init();
	start_timer = timer_get();

	/* register atexit */
	if ( -1 == atexit(clean_up) ) {
		puts(err_unable_to_register_atexit);
		return 1;
	}

	/* show banner */
	printf(banner, netcdf_get_version());

	/* parse arguments */
	printf("parsing arguments...");
	if ( ! parse_args(argc, argv) ) {
		return 1;
	}
	puts(msg_ok);
      
	if ( ! output_path_create() ) {
		return 0;
	}

	if ( ! parameterization_output_create() ) {
		puts("Unable to create parameterization output path\n");
		return 0;
	}

	if ( g_sz_benchmark_path && ! g_sz_output_path ) {
		puts("unable to benchmark. output path not specified!");
		return 0;
	}

#ifdef BENCHMARK_ONLY
	goto benchmark;
#endif

	printf("import settings file %s...", g_sz_settings_file);
	g_settings = settings_import(g_sz_settings_file);
	if ( ! g_settings ) {
		goto err;
	}
        puts(msg_ok);

	if ( ! g_settings->Ndep_fixed && ! g_sz_ndep_file ) {
		puts("ndep file not specified for ndep_fixed off");
		goto err;
	}
	

       // below old part of the code v5.5

	 /* if ( MANAGEMENT_VAR == g_settings->management )   
          {    
          char *p;
		char buffer[256];

		sprintf(buffer, "ISIMIP/%s_management_ISIMIP.txt", g_settings->sitename);

		printf("import management file %s...", buffer);
		if ( g_sz_input_path )
		{
			p = concatenate_path(g_sz_input_path, buffer); 
			if ( ! p )
			{
				puts(sz_err_out_of_memory);
				goto err;
			}
		}
		else
		{
			p = buffer;
		}

		g_management = management_load(p);
		if ( g_sz_input_path ) free(p);
		//if ( ! g_management ) goto err;
		puts(msg_ok);
	}  */
        
           // 5p6 load management file (when MANAGEMENT == VAR and the thinning/clear cut is prescribed)

           if ( MANAGEMENT_VAR == g_settings->management )  
	   {
		char *p;
		char buffer[256];
 

               // if management file info is provided, check if the file exists
               // g_sz_manag_file  : gives already the path to the file
   
	        if ( g_sz_manag_file ) {      // not necessarily we do need the management file (e.s.if MAN=VAR and the simulation cover the historical only 
                                              // with all density data observed
		  p = concatenate_path(g_sz_input_path, g_sz_manag_file); 
                  printf("g_sz_manag_file concatenate P  %s\n", p);
		  if ( ! p ) 
	     	  {
			puts(sz_err_out_of_memory);
			goto err;
		  }
		  //free(g_sz_manag_file);
		  ////  g_sz_manag_file = p;  //
		  
		   if ( ! file_copy(p, g_sz_input_data_path) )
	           {
		    printf("warning: unable to copy management file %s to %s\n", g_sz_manag_file, g_sz_input_data_path);
	           }
	

                  g_management = management_load_presc(p);
 
                  if ( g_sz_input_path ) free(p);
                  if ( ! g_management ) goto err;

		  puts(msg_ok);

	        }            
	   }
 
	printf("soil import...");
	s = soil_settings_import(g_sz_soil_file, &soil_settings_count);
	if ( ! s ) {
		goto err;
	}
	puts(msg_ok);
       

	printf("topo import...");
	t = topo_import(g_sz_topo_file, &topos_count);
	if ( ! t ) {
		goto err;
	}
   
	puts(msg_ok);
       
	printf("build matrix");
	if ( g_sz_dataset_file )
                printf(" using %s...", g_sz_dataset_file);		
	else
		printf("...");
       
	matrix = matrix_create(s, soil_settings_count, g_sz_dataset_file, &g_dataset);
        
	if ( ! matrix ) goto err;
	puts(msg_ok);

	// save input data files
	if ( ! file_copy(g_sz_dataset_file, g_sz_input_data_path) )
	{
		printf("warning: unable to copy dataset file %s to %s\n", g_sz_dataset_file, g_sz_input_data_path);
	}

	if ( ! file_copy(g_sz_settings_file, g_sz_input_data_path) )
	{
		printf("warning: unable to copy settings file %s to %s\n", g_sz_settings_file, g_sz_input_data_path);
	}

	if ( ! file_copy(g_sz_soil_file, g_sz_input_data_path) )
	{
		printf("warning: unable to copy soil file %s to %s\n", g_sz_soil_file, g_sz_input_data_path);
	}

	if ( ! file_copy(g_sz_topo_file, g_sz_input_data_path) )
	{
		printf("warning: unable to copy topo file %s to %s\n", g_sz_topo_file, g_sz_input_data_path);
	}

	g_year_start_index = -1;

#if 0
	/* add elev to each cells */
	for ( cell = 0; cell < matrix->cells_count; ++cell )
	{
		int i;
		int index;

		for ( i = 0; i < topos_count; ++i )
		{
			if ( (matrix->cells[cell].x == (int)t[i].values[TOPO_X])
					&& (matrix->cells[cell].y == (int)t[i].values[TOPO_Y]) )
			{
				index = i;
				break;
			}
		}

		if ( -1 == index )
		{
			logger_error(g_debug_log, "no topo settings found for cell at %d,%d\n"
					, matrix->cells[cell].x
					, matrix->cells[cell].y
			);
			continue;
		}

		matrix->cells[cell].elev = t[index].values[TOPO_ELEV];
	}

	/* sort cells */
	qsort(matrix->cells, matrix->cells_count, sizeof*matrix->cells, cells_sort);
#endif

	for ( cell = 0; cell < matrix->cells_count; ++cell )
	{
		/********************************** IMPORT SOIL AND TOPO **********************************/
		/* set g_soil_settings and g_topo */
		{
			int i;
			int index = -1;

			/* soil */
			logger(g_debug_log, "Processing soil data files for cell at %d,%d...\n", matrix->cells[cell].x, matrix->cells[cell].y);
			logger(g_debug_log, "input_soil_path = %s\n", g_sz_soil_file);

			for ( i = 0; i < soil_settings_count; ++i )
			{
				if ( (matrix->cells[cell].x == (int)s[i].values[SOIL_X])
						&& (matrix->cells[cell].y == (int)s[i].values[SOIL_Y]) )
				{
					index = i;
					break;
				}
			}

			if ( -1 == index )
			{
				logger_error(g_debug_log, "no soil settings found for cell at %d,%d\n"
						, matrix->cells[cell].x
						, matrix->cells[cell].y
				);
				continue;
			}
			g_soil_settings = &s[index];

			/* topo */
			logger(g_debug_log, "Processing topo data files for cell at %d,%d...\n", matrix->cells[cell].x, matrix->cells[cell].y);
			logger(g_debug_log, "input_topo_path = %s\n", g_sz_topo_file);
			for ( i = 0; i < topos_count; ++i )
			{
				if ( (matrix->cells[cell].x == (int)t[i].values[TOPO_X])
						&& (matrix->cells[cell].y == (int)t[i].values[TOPO_Y]) )
				{
					index = i;
					break;
				}
			}

			if ( -1 == index )
			{
				logger_error(g_debug_log, "no topo settings found for cell at %d,%d\n"
						, matrix->cells[cell].x
						, matrix->cells[cell].y
				);
				continue;
			}
			g_topo = &t[index];
		}

		/* only for first cell */
		if ( 0 == cell ) {
			if ( ! log_start(g_settings->sitename) ) {
				goto err;
			}
			logger(g_debug_log, "\n3D-CMCC-FEM STARTs....\n\n");
		}

		/* check soil values */
		if ( ! check_soil_values() ) {
			goto err;
		}

		/* check topo values */
		if ( ! check_topo_values() ) {
			goto err;
		}

		/* check hemisphere */
		if ( g_soil_settings->values[SOIL_LAT] > 0 )
		{
			matrix->cells[cell].north = 0;
		} else {
			matrix->cells[cell].north = 1;
		}


		// TODO make it better!
		/* import SOIL data from soil.txt and assign to matrix cells variables */
		/* (assign global g_soil_setting_variables to single cells) */
		matrix->cells[cell].lat           = g_soil_settings->values[SOIL_LAT];
		matrix->cells[cell].lon           = g_soil_settings->values[SOIL_LON];
		matrix->cells[cell].clay_perc     = g_soil_settings->values[SOIL_CLAY_PERC];
		matrix->cells[cell].silt_perc     = g_soil_settings->values[SOIL_SILT_PERC];
		matrix->cells[cell].sand_perc     = g_soil_settings->values[SOIL_SAND_PERC];
		matrix->cells[cell].soil_depth    = g_soil_settings->values[SOIL_DEPTH];
		matrix->cells[cell].fr            = g_soil_settings->values[SOIL_FR];
		matrix->cells[cell].fn0           = g_soil_settings->values[SOIL_FN0];
		matrix->cells[cell].fnn           = g_soil_settings->values[SOIL_FNN];
		matrix->cells[cell].m0            = g_soil_settings->values[SOIL_M0];
		matrix->cells[cell].init_litter_C = g_soil_settings->values[LITTERC];
		matrix->cells[cell].init_litter_N = g_soil_settings->values[LITTERN];
		matrix->cells[cell].init_soil_C   = g_soil_settings->values[SOILC];
		matrix->cells[cell].init_soil_N   = g_soil_settings->values[SOILN];
		matrix->cells[cell].init_dead_C   = g_soil_settings->values[DEADWOODC];

		/* import TOPO data from topo.txt and assign to matrix cells variables */
		/* (assign global g_topo variables to single cells) */
		matrix->cells[cell].elev          = g_topo->values[TOPO_ELEV];


		/********************************** IMPORT MET DATA **********************************/

		logger(g_debug_log, "Processing met data files for cell at %d,%d...\n", matrix->cells[cell].x, matrix->cells[cell].y);
		logger(g_debug_log, "input_met_path = %s\n", g_sz_input_met_file);
		logger_error(g_debug_log, "importing met data...");

		if ( ! import_meteo_data(g_sz_input_met_file, &years_of_simulation, &matrix->cells[cell]) ) goto err;
		logger_error(g_debug_log, "ok\n");

		/* set start year index */
		if ( -1 == g_year_start_index )
		{
			int i;
			int ii;

			for ( i = 0; i < years_of_simulation; ++i ) {
				if ( g_settings->year_start == matrix->cells[0].years[i].year ) {
					g_year_start_index = i;
					break;
				}
			}
			if ( -1 == g_year_start_index ) {
				logger_error(g_debug_log, "start year (%d) not found. range is %d-%d.\n"
						, g_settings->year_start
						, matrix->cells[0].years[0].year
						, matrix->cells[0].years[years_of_simulation-1].year
				);
				goto err;
			}

			/* set end year (adjusting years_of_simulation) */
			ii = -1;
			for ( i = 0; i < years_of_simulation; ++i ) {
				if ( g_settings->year_end == matrix->cells[0].years[i].year ) {
					if ( g_year_start_index > i ) {
						logger_error(g_debug_log, "start year (%d) cannot be > end year (%d)\n"
								, g_settings->year_start
								, g_settings->year_end
						);
						g_year_start_index = -1;
						goto err;
					}
					ii = i;
					break;
				}
			}
			if ( -1 == ii ) {
				logger_error(g_debug_log, "end year (%d) not found. end year will be changed to %d\n"
						, g_settings->year_end
						, matrix->cells[0].years[years_of_simulation-1].year
				);
				g_settings->year_end = matrix->cells[0].years[years_of_simulation-1].year;
				ii = years_of_simulation-1;
				logger_error(g_debug_log, "renaming logs...");
				if ( ! log_rename() )
				{
					logger_error(g_debug_log, "error!\n");
					goto err;
				}
				else
				{
					logger_error(g_debug_log, "ok\n");
				}
			}
			/* i is year_end_index */
			i = ii + 1;

			years_of_simulation = i - g_year_start_index;
		}
		else
		{
			years_of_simulation -= g_year_start_index;
		}

		/* move pointer for year */
		matrix->cells[cell].years += g_year_start_index;
	}

	logger(g_debug_log, "Total years_of_simulation = %d\n", years_of_simulation);
	logger(g_debug_log, "***************************************************\n\n");

	// sort
#ifdef USE_NEW_OUTPUT
	sort_all(matrix);
#endif

	/* for monthly and yearly means */
	for ( cell = 0; cell < matrix->cells_count; ++cell )
	{
		for ( year = 0; year < years_of_simulation; ++year )
		{
			matrix->cells[cell].years[year].yearly_mean.solar_rad     = 0.;
			matrix->cells[cell].years[year].yearly_mean.tavg          = 0.;
			matrix->cells[cell].years[year].yearly_mean.tmax          = 0.;
			matrix->cells[cell].years[year].yearly_mean.tmin          = 0.;
			matrix->cells[cell].years[year].yearly_mean.tday          = 0.;
			matrix->cells[cell].years[year].yearly_mean.tnight        = 0.;
			matrix->cells[cell].years[year].yearly_mean.vpd           = 0.;
			matrix->cells[cell].years[year].yearly_mean.prcp          = 0.;
			matrix->cells[cell].years[year].yearly_mean.tsoil         = 0.;
			matrix->cells[cell].years[year].yearly_mean.rh_f          = 0.;
			matrix->cells[cell].years[year].yearly_mean.incoming_par  = 0.;
			matrix->cells[cell].years[year].yearly_mean.par           = 0.;

			for ( month = 0; month < MONTHS_COUNT; ++month )
			{
				matrix->cells[cell].years[year].monthly_mean[month].solar_rad     = 0.;
				matrix->cells[cell].years[year].monthly_mean[month].tavg          = 0.;
				matrix->cells[cell].years[year].monthly_mean[month].tmax          = 0.;
				matrix->cells[cell].years[year].monthly_mean[month].tmin          = 0.;
				matrix->cells[cell].years[year].monthly_mean[month].tday          = 0.;
				matrix->cells[cell].years[year].monthly_mean[month].tnight        = 0.;
				matrix->cells[cell].years[year].monthly_mean[month].vpd           = 0.;
				matrix->cells[cell].years[year].monthly_mean[month].prcp          = 0.;
				matrix->cells[cell].years[year].monthly_mean[month].tsoil         = 0.;
				matrix->cells[cell].years[year].monthly_mean[month].rh_f          = 0.;
				matrix->cells[cell].years[year].monthly_mean[month].incoming_par  = 0.;
				matrix->cells[cell].years[year].monthly_mean[month].par           = 0.;
			}
		}
	}


	for ( year = 0; year < years_of_simulation; ++year )
	{
		int current_doy;
		int days_per_month;
		int leap_year;

		current_doy = 0;
		leap_year = IS_LEAP_YEAR(g_settings->year_start + year);

		if ( g_settings->year_restart == g_settings->year_start+year )
		{
			if ( ! restart(matrix, g_settings->year_restart) ) goto err;
		}

		for ( month = 0; month < MONTHS_COUNT; ++month )
		{
			days_per_month = DaysInMonth[month];

			if ( leap_year && (FEBRUARY == month) )
			{
				++days_per_month;
			}

			for ( day = 0; day < 31; ++day )
			{
				if ( day >= days_per_month ) 
				{
					break;
				}

				if( !day && !month && !year )
				{
					/* summary on model simulation */
					simulation_summary ( matrix );

					/* summary on site data */
					site_summary       ( matrix );

					/* summary on topographic data */
					topo_summary       ( matrix );

					/* summary on soil data */
					soil_summary       ( matrix, matrix->cells );

					/* if spin up off */
					if ( ! g_settings->spinup )
					{
						if( LANDUSE_F == g_soil_settings->landuse )
						{
							/* forest initialization */
							forest_initialization ( matrix, day, month, year );
						}
						else
						{
							/* include initialization for other land uses */
						}
					}

					/* litter initialization */
					litter_initialization ( matrix, day, month, year );

					/* soil initialization */
					soil_initialization   ( matrix, day, month, year );

				}

				for ( cell = 0; cell < matrix->cells_count; ++cell )
				{

					/* compute daily climate variables not coming from met data */
					Daily_avg_temperature       ( matrix->cells[cell].years[year].m, day, month );
					Daylight_avg_temperature    ( matrix->cells[cell].years[year].m, day, month );
					Nightime_avg_temperature    ( matrix->cells[cell].years[year].m, day, month );
					Soil_temperature            ( &matrix->cells[cell], day, month, year );
					Thermic_sum                 ( &matrix->cells[cell], matrix->cells[cell].years[year].m, day, month, year );
					Air_density                 ( matrix->cells[cell].years[year].m, day, month );
					Day_Length                  ( &matrix->cells[cell], day, month, year );
					Latent_heat                 ( matrix->cells[cell].years[year].m, day, month );
					Air_pressure                ( matrix->cells[cell].years[year].m, day, month );
					Psychrometric               ( matrix->cells[cell].years[year].m, day, month );
					Sat_vapour_pressure         ( &matrix->cells[cell], day, month, year );
					Dew_temperature             ( matrix->cells[cell].years[year].m, day, month );
					Radiation                   ( &matrix->cells[cell], day, month, year );
					Check_prcp                  ( &matrix->cells[cell], day, month, year );
					Daily_Ndeposition           ( &matrix->cells[cell], day, month, year );
					Averaged_temperature        ( &matrix->cells[cell], AVERAGED_TAVG, day, month, year );
					Averaged_temperature        ( &matrix->cells[cell], AVERAGED_TDAY, day, month, year );
					Averaged_temperature        ( &matrix->cells[cell], AVERAGED_TNIGHT, day, month, year );
					Averaged_temperature        ( &matrix->cells[cell], AVERAGED_TSOIL, day, month, year );
					Weighted_average_temperature( &matrix->cells[cell], WEIGHTED_MEAN_TAVG, day, month, year );
					Weighted_average_temperature( &matrix->cells[cell], WEIGHTED_MEAN_TDAY, day, month, year );
					Weighted_average_temperature( &matrix->cells[cell], WEIGHTED_MEAN_TNIGHT, day, month, year );
					Weighted_average_temperature( &matrix->cells[cell], WEIGHTED_MEAN_TSOIL, day, month, year );

					if ( LANDUSE_F == g_soil_settings->landuse )
					{
						/* compute annually the days for the growing season BEFORE any other process */
						Veg_Days (&matrix->cells[cell], day, month, year);
					}
					else
					{
						/* include other land use ???? */
					}

					/* for monthly mean */
					matrix->cells[cell].years[year].monthly_mean[month].solar_rad     += matrix->cells[cell].years[year].m[month].d[day].solar_rad;
					matrix->cells[cell].years[year].monthly_mean[month].tavg          += matrix->cells[cell].years[year].m[month].d[day].tavg;
					matrix->cells[cell].years[year].monthly_mean[month].tmax          += matrix->cells[cell].years[year].m[month].d[day].tmax;
					matrix->cells[cell].years[year].monthly_mean[month].tmin          += matrix->cells[cell].years[year].m[month].d[day].tmin;
					matrix->cells[cell].years[year].monthly_mean[month].tday          += matrix->cells[cell].years[year].m[month].d[day].tday;
					matrix->cells[cell].years[year].monthly_mean[month].tnight        += matrix->cells[cell].years[year].m[month].d[day].tnight;
					matrix->cells[cell].years[year].monthly_mean[month].vpd           += matrix->cells[cell].years[year].m[month].d[day].vpd;
					matrix->cells[cell].years[year].monthly_mean[month].prcp          += matrix->cells[cell].years[year].m[month].d[day].prcp;
					matrix->cells[cell].years[year].monthly_mean[month].tsoil         += matrix->cells[cell].years[year].m[month].d[day].tsoil;
					matrix->cells[cell].years[year].monthly_mean[month].rh_f          += matrix->cells[cell].years[year].m[month].d[day].rh_f;
					matrix->cells[cell].years[year].monthly_mean[month].incoming_par  += matrix->cells[cell].years[year].m[month].d[day].incoming_par;
					matrix->cells[cell].years[year].monthly_mean[month].par           += matrix->cells[cell].years[year].m[month].d[day].par;
				}
			}

			for ( cell = 0; cell < matrix->cells_count; ++cell )
			{
				/* compute monthly mean */
				matrix->cells[cell].years[year].monthly_mean[month].solar_rad /= days_per_month;
				matrix->cells[cell].years[year].monthly_mean[month].tavg /= days_per_month;
				matrix->cells[cell].years[year].monthly_mean[month].tmax /= days_per_month;
				matrix->cells[cell].years[year].monthly_mean[month].tmin          /= days_per_month;
				matrix->cells[cell].years[year].monthly_mean[month].tday          /= days_per_month;
				matrix->cells[cell].years[year].monthly_mean[month].tnight        /= days_per_month;
				matrix->cells[cell].years[year].monthly_mean[month].vpd           /= days_per_month;
				//matrix->cells[cell].years[year].monthly_mean[month].prcp /= days_per_month; // we need accumul, not mean
				matrix->cells[cell].years[year].monthly_mean[month].tsoil         /= days_per_month;
				matrix->cells[cell].years[year].monthly_mean[month].rh_f          /= days_per_month;
				matrix->cells[cell].years[year].monthly_mean[month].incoming_par  /= days_per_month;
				matrix->cells[cell].years[year].monthly_mean[month].par           /= days_per_month;

				/* for yearly mean */
				matrix->cells[cell].years[year].yearly_mean.solar_rad += matrix->cells[cell].years[year].monthly_mean[month].solar_rad;
				matrix->cells[cell].years[year].yearly_mean.tavg += matrix->cells[cell].years[year].monthly_mean[month].tavg;
				matrix->cells[cell].years[year].yearly_mean.tmax += matrix->cells[cell].years[year].monthly_mean[month].tmax;
				matrix->cells[cell].years[year].yearly_mean.tmin += matrix->cells[cell].years[year].monthly_mean[month].tmin;
				matrix->cells[cell].years[year].yearly_mean.tday += matrix->cells[cell].years[year].monthly_mean[month].tday;
				matrix->cells[cell].years[year].yearly_mean.tnight += matrix->cells[cell].years[year].monthly_mean[month].tnight;
				matrix->cells[cell].years[year].yearly_mean.vpd += matrix->cells[cell].years[year].monthly_mean[month].vpd;
				matrix->cells[cell].years[year].yearly_mean.prcp += matrix->cells[cell].years[year].monthly_mean[month].prcp;
				matrix->cells[cell].years[year].yearly_mean.tsoil += matrix->cells[cell].years[year].monthly_mean[month].tsoil;
				matrix->cells[cell].years[year].yearly_mean.rh_f += matrix->cells[cell].years[year].monthly_mean[month].rh_f;
				matrix->cells[cell].years[year].yearly_mean.incoming_par += matrix->cells[cell].years[year].monthly_mean[month].incoming_par;
				matrix->cells[cell].years[year].yearly_mean.par += matrix->cells[cell].years[year].monthly_mean[month].par;
			}
		}

		for ( cell = 0; cell < matrix->cells_count; ++cell )
		{
			/* compute yearly mean */
			matrix->cells[cell].years[year].yearly_mean.solar_rad /= MONTHS_COUNT;
			matrix->cells[cell].years[year].yearly_mean.tavg /= MONTHS_COUNT;
			matrix->cells[cell].years[year].yearly_mean.tmax /= MONTHS_COUNT;
			matrix->cells[cell].years[year].yearly_mean.tmin /= MONTHS_COUNT;
			matrix->cells[cell].years[year].yearly_mean.tday /= MONTHS_COUNT;
			matrix->cells[cell].years[year].yearly_mean.tnight /= MONTHS_COUNT;
			matrix->cells[cell].years[year].yearly_mean.vpd /= MONTHS_COUNT;
			//matrix->cells[cell].years[year].yearly_mean.prcp /= MONTHS_COUNT; // we need accumul, not mean
			matrix->cells[cell].years[year].yearly_mean.tsoil /= MONTHS_COUNT;
			matrix->cells[cell].years[year].yearly_mean.rh_f /= MONTHS_COUNT;
			matrix->cells[cell].years[year].yearly_mean.incoming_par /= MONTHS_COUNT;
			matrix->cells[cell].years[year].yearly_mean.par /= MONTHS_COUNT;
		}

		for ( month = 0; month < MONTHS_COUNT; ++month )
		{
			days_per_month = DaysInMonth[month];

			if ( leap_year && (FEBRUARY == month) )
			{
				++days_per_month;
			}

			current_doy += days_per_month;

			for ( day = 0; day < 31; ++day )
			{
				if ( day >= days_per_month ) 
				{
					break;
				}

				for ( cell = 0; cell < matrix->cells_count; ++cell )
				{

					/* counter "day of the year" */
					if( !day && !month )matrix->cells[cell].doy = 1;
					else ++matrix->cells[cell].doy;

					/* counter "days of simulation" */
					if( !day && !month && !year )matrix->cells[cell].dos = 1;
					else ++matrix->cells[cell].dos;

					/* print daily met data */
					if (matrix->cells_count > 0.) print_daily_met_data (&matrix->cells[cell], day, month, year);

					/* print cell data */
					if (matrix->cells_count > 0.) print_daily_cell_data ( &matrix->cells[cell] );

					/************************************************************************/
					/* note: if spinup 'on' tree_model_daily doesn't run */
					if ( ! g_settings->spinup )
					{
						if ( g_sz_dataset_file )
						{
							/* run tree model daily */
							if ( ( LANDUSE_F == g_soil_settings->landuse ) && ( matrix->cells[cell].heights_count != 0 ) )
							{
								if ( 'f' == g_settings->version )
								{
									if ( matrix->cells[cell].n_trees > 0 )
									{
										if ( !Tree_model_daily( matrix, cell, day, month, year ) )
										{
											logger(g_debug_log, "tree model daily failed!!!\n");
											goto err;
										}
										else
										{
											//printf("ok tree_model (x=%d,y=%d) (%02d-%02d-%d)\n"
											//		, matrix->cells[cell].x
											//		, matrix->cells[cell].y
											//		, day+1, month+1, year+g_settings->year_start
											//);
										}
									}
								}
							}
							else
							{
								/* run for possible other Land Use versions */
							}
						}
					}
					/************************************************************************/
					/* pruning */
					// v 5p6 not yet active
					if ( g_settings->prunings_count ) {
						int i;

						for ( i = 0; i < g_settings->prunings_count; ++i ) {
							if ( (year+g_settings->year_start == g_settings->prunings[i].year)
									&& (month == g_settings->prunings[i].month-1)
									&& (day == g_settings->prunings[i].day-1) ) {
								pruning_daily(matrix, cell, &g_settings->prunings[i]);
							}
						}
					}

					/************************************************************************/
					/* run for litter model */
					if ( ! Litter_model_daily(matrix, cell, day, month, year) )
					{
						logger_error(g_debug_log, "litter model daily failed!!!\n");
						goto err;
					}
					else
					{
						//printf("ok litr_model (x=%d,y=%d) (%02d-%02d-%d)\n"
						//		, matrix->cells[cell].x
						//		, matrix->cells[cell].y
						//		, day+1, month+1, year+g_settings->year_start
						//);
					}
					/************************************************************************/
					/* run for soil model */
					if ( ! Soil_model_daily(matrix, cell, day, month, year) )
					{
						logger_error(g_debug_log, "soil model daily failed!!!\n");
						goto err;
					}
					else
					{
						//printf("ok soil_model (x=%d,y=%d) (%02d-%02d-%d)\n"
						//		, matrix->cells[cell].x
						//		, matrix->cells[cell].y
						//		, day+1, month+1, year+g_settings->year_start
						//);
					}
					/************************************************************************/
					/* run for cell model */
					if ( !Cell_model_daily(matrix, cell, day, month, year) )
					{
						logger_error(g_debug_log, "cell model daily failed!!!\n");
						goto err;
					}
					else
					{
						printf("ok cell_model (x=%d,y=%d) (%02d-%02d-%d)\n"
								, matrix->cells[cell].x
								, matrix->cells[cell].y
								, day+1, month+1, year+g_settings->year_start
						);
					}
					/*************************************************************************/
					/******************************************************************************/
					/* print daily output */
					//EOD_print_output_cell_level (&matrix->cells[cell], day, month, year, years_of_simulation );
                                       //EOD_print_output_cell_level_ddalmo (&matrix->cells[cell], day, month, year, years_of_simulation );
                    EOD_print_output_cell_level_mc(&matrix->cells[cell], day, month, year, years_of_simulation );
					 EOD_print_output_soil_cell_level (&matrix->cells[cell], day, month, year, years_of_simulation );

					/* reset daily variables once printed */
					reset_daily_class_variables ( &matrix->cells[cell] );
					reset_daily_layer_variables ( &matrix->cells[cell] );
					reset_daily_cell_variables  ( &matrix->cells[cell] );

					logger(g_debug_log, "******************* END OF DAY (%d) *******************\n\n\n", day + 1 );

					/* end of month */
					if ( current_doy == matrix->cells[cell].doy )
					{
						/* print monthly output */
						//EOM_print_output_cell_level( &matrix->cells[cell], month, year, years_of_simulation );
                                               //EOM_print_output_cell_level_ddalmo( &matrix->cells[cell], month, year, years_of_simulation );
                                               EOM_print_output_cell_level_mc( &matrix->cells[cell], month, year, years_of_simulation );
						EOM_print_output_soil_cell_level( &matrix->cells[cell], month, year, years_of_simulation );

						reset_monthly_class_variables ( &matrix->cells[cell] );
						reset_monthly_layer_variables ( &matrix->cells[cell] );
						reset_monthly_cell_variables  ( &matrix->cells[cell] );

						logger(g_debug_log, "******************* END OF MONTH (%d) *******************\n\n", month + 1 );
					}

					/* reset annual variables once printed */
					if ( matrix->cells[cell].doy == (leap_year ? 366 : 365)  )
					{

						/* print annual output */

						//EOY_print_output_cell_level( &matrix->cells[cell], year, years_of_simulation );
                                               //EOY_print_output_cell_level_ddalmo( &matrix->cells[cell], year, years_of_simulation );
                                               EOY_print_output_cell_level_mc( &matrix->cells[cell], year, years_of_simulation );
						EOY_print_output_soil_cell_level( &matrix->cells[cell], year, years_of_simulation );

						reset_annual_class_variables ( &matrix->cells[cell] );
						reset_annual_layer_variables ( &matrix->cells[cell] );
						reset_annual_cell_variables  ( &matrix->cells[cell] );

						logger(g_debug_log, "******************* END OF YEAR (%d) *******************\n\n\n\n\n\n", matrix->cells[cell].years[year].year );
					}
				}
			}
		}
	}

	EOD_cell_msg();
	EOD_soil_msg();
	EOM_cell_msg();
	EOM_soil_msg();
	EOY_cell_msg();
	EOY_soil_msg();

	for ( cell = 0; cell < matrix->cells_count; ++cell )
	{
		if ( ! matrix->cells_count ) break;

		if ( g_year_start_index != -1 ) {
			free(matrix->cells[cell].years-g_year_start_index);
		} else {
			free(matrix->cells[cell].years);
		}
		matrix->cells[cell].years = NULL; /* required */
	}

	/* free memory */
	matrix_free(matrix); matrix = NULL;

#ifdef NC_USE       
	/* create nc files */
	if ( g_daily_log && g_settings->netcdf_output ) {
		printf("creating nc daily files...");
		logger_flush(g_daily_log);
		
		if ( ! convert_to_nc(g_daily_log->filename) ) {
			goto err;
		}
		puts("ok");
	}

	if ( g_annual_log && g_settings->netcdf_output ) {
		printf("creating nc annual files...");
		logger_flush(g_annual_log);
		
		if ( ! convert_to_nc(g_annual_log->filename) ) {
			goto err;
		}
		puts("ok");
	}
#endif
	/* ok ! */
	prog_ret = 0;

#ifdef BENCHMARK_ONLY
	benchmark:
#endif

	/* benchmark ? */
	if ( g_sz_benchmark_path ) {
		logger_error(g_debug_log, "benchmark with '%s'...", g_sz_benchmark_path);
		if ( compare(g_sz_output_path, g_sz_benchmark_path) ) {
			puts("ok!");
		}
	}

	err:

	/* cleanup memory...
		do not remove null pointer to prevent
		double free with clean_up func
	 */

	/* close logger */
	logger_close(g_annual_soil_log); g_annual_soil_log = NULL;
	logger_close(g_monthly_soil_log); g_monthly_soil_log = NULL;
	logger_close(g_daily_soil_log); g_daily_soil_log = NULL;
	logger_close(g_annual_log); g_annual_log = NULL;
	logger_close(g_monthly_log); g_monthly_log = NULL;
	logger_close(g_daily_log); g_daily_log = NULL;

	/* free memory */
	if ( t ) free(t);
	if ( s ) free(s);
	if ( matrix ) matrix_free(matrix);

	free(g_sz_input_met_file); g_sz_input_met_file = NULL;
	free(g_sz_settings_file); g_sz_settings_file = NULL;

	/* this should be freed before */
	free(g_sz_topo_file); g_sz_topo_file = NULL;
	free(g_sz_soil_file); g_sz_soil_file = NULL;

	/* print elapsed time */
	end_timer = timer_get();
	timer = end_timer - start_timer;
	logger_error(g_debug_log, "%.2f secs / %.2f mins / %.2f hours elapsed\n", timer, timer / 60.0, timer / 3600.0);
	logger_close(g_debug_log); g_debug_log = NULL;

	return prog_ret;
}
#else
//note: 02/february/2017 
// the model runs for all days/months/years and then changes cell

int main(int argc, char *argv[]) {
	char sz_date[32]; // should be enough
	int ret;
	int year;
	int month;
	int day;
	int cell;
	int prog_ret;
	int soil_settings_count;
	int topos_count;
	double timer;
	double start_timer;
	double end_timer;
	matrix_t* matrix;
	output_t* output_vars;
	soil_settings_t* s;
	topo_t* t;

	//_CrtSetBreakAlloc();

	/* initialize */
	matrix = NULL;
	output_vars = NULL;
	t = NULL;
	s = NULL;
	prog_ret = 1;

	/* start timer */
	timer_init();
	start_timer = timer_get();

	/* get current date */
	{
		struct tm* ptm;
		time_t t;

		/* get current date */
		time(&t);
		ptm = gmtime(&t);
		sprintf(sz_date, "%04d_%s_%02d"
				, ptm->tm_year + 1900
				, szMonth[ptm->tm_mon]
						  , ptm->tm_mday
		);
	}

	/* register atexit */
	if ( -1 == atexit(clean_up) ) {
		puts(err_unable_to_register_atexit);
		return 1;
	}

	/* show banner */
	printf(banner, netcdf_get_version());

	/* parse arguments */
	printf("parsing arguments...");
	if ( ! parse_args(argc, argv) ) {
		return 1;
	}
	puts(msg_ok);

	if ( ! parameterization_output_create(sz_date) ) {
		puts("Unable to create parameterization output path\n");
		return 0;
	}

	/* import output vars file ? */
	if ( g_sz_output_vars_file ) {
		printf("import output file...");
		output_vars = output_import(g_sz_output_vars_file);
		free(g_sz_output_vars_file);
		g_sz_output_vars_file = NULL;
		if ( ! output_vars ) {
			return 1;
		}
		puts(msg_ok);
	}

	if ( g_sz_benchmark_path && ! g_sz_output_path ) {
		puts("unable to benchmark. output path not specified!");
		return 0;
	}

#ifdef BENCHMARK_ONLY
	goto benchmark;
#endif

	printf("import settings file %s...", g_sz_settings_file);
        
	g_settings = settings_import(g_sz_settings_file);
	if ( ! g_settings )
	{
		goto err;
	}
          
      
	if ( g_settings->time != 'd' ) {
		puts("uncorrect time step choiced!");
		goto err;
	} else {
		puts(msg_ok);
	}

	if ( ! g_settings->Ndep_fixed && ! g_sz_ndep_file ) {
		puts("ndep file not specified for ndep_fixed off");
		goto err;
	}
        
       

	printf("soil import...");
	s = soil_settings_import(g_sz_soil_file, &soil_settings_count);
	if ( ! s ) {
		goto err;
	}
	puts(msg_ok);
     
        
	printf("topo import...");
      
	t = topo_import(g_sz_topo_file, &topos_count);
	if ( ! t ) {
		goto err;
	}
        
	puts(msg_ok);
         
       

	printf("build matrix");
	if ( g_sz_dataset_file )
		printf(" using %s...", g_sz_dataset_file);
	else
		printf("...");
	matrix = matrix_create(s, soil_settings_count, g_sz_dataset_file);
	if ( ! matrix ) goto err;
	puts(msg_ok);

	g_year_start_index = -1;

	logger(g_debug_log, "\n3D-CMCC FEM START....\n\n");

#if 0
	/* add elev to each cells */
	for ( cell = 0; cell < matrix->cells_count; ++cell )
	{
		int i;
		int index;

		for ( i = 0; i < topos_count; ++i )
		{
			if ( (matrix->cells[cell].x == (int)t[i].values[TOPO_X])
					&& (matrix->cells[cell].y == (int)t[i].values[TOPO_Y]) )
			{
				index = i;
				break;
			}
		}

		if ( -1 == index )
		{
			logger_error(g_debug_log, "no topo settings found for cell at %d,%d\n"
					, matrix->cells[cell].x
					, matrix->cells[cell].y
			);
			continue;
		}

		matrix->cells[cell].elev = t[index].values[TOPO_ELEV];
	}

	/* sort cells */
	qsort(matrix->cells, matrix->cells_count, sizeof*matrix->cells, cells_sort);
#endif


	/********************************** IMPORT SOIL AND TOPO **********************************/
	for ( cell = 0; cell < matrix->cells_count; ++cell )
	{
		logger(g_debug_log, "Processing met data files for cell at %d,%d...\n", matrix->cells[cell].x, matrix->cells[cell].y);
		logger(g_debug_log, "input_met_path = %s\n", g_sz_input_met_file);

		/* set g_soil_settings and g_topo */
		{
			int i;
			int index = -1;

			/* soil settings */
			for ( i = 0; i < soil_settings_count; ++i )
			{
				if ( (matrix->cells[cell].x == (int)s[i].values[SOIL_X])
						&& (matrix->cells[cell].y == (int)s[i].values[SOIL_Y]) )
				{
					index = i;
					break;
				}
			}

			if ( -1 == index )
			{
				logger_error(g_debug_log, "no soil settings found for cell at %d,%d\n"
						, matrix->cells[cell].x
						, matrix->cells[cell].y
				);
				continue;
			}
			g_soil_settings = &s[index];

			/* topo */
			for ( i = 0; i < topos_count; ++i )
			{
				if ( (matrix->cells[cell].x == (int)t[i].values[TOPO_X])
						&& (matrix->cells[cell].y == (int)t[i].values[TOPO_Y]) )
				{
					index = i;
					break;
				}
			}

			if ( -1 == index )
			{
				logger_error(g_debug_log, "no topo settings found for cell at %d,%d\n"
						, matrix->cells[cell].x
						, matrix->cells[cell].y
				);
				continue;
			}
			g_topo = &t[index];
		}

		/* only for first cell */
		if ( 0 == cell ) {
			if ( ! log_start(sz_date, g_settings->sitename) ) {
				goto err;
			}
		}

		if (	IS_INVALID_VALUE(g_soil_settings->values[SOIL_SAND_PERC])
				|| IS_INVALID_VALUE(g_soil_settings->values[SOIL_CLAY_PERC])
				|| IS_INVALID_VALUE(g_soil_settings->values[SOIL_SILT_PERC])
				|| IS_INVALID_VALUE(g_soil_settings->values[SOIL_DEPTH]) ) {
			logger_error(g_debug_log, "NO SOIL DATA AVAILABLE");
			goto err;
		}

		/* check hemisphere */
		if ( g_soil_settings->values[SOIL_LAT] > 0 ) {
			matrix->cells[cell].north = 0;
		} else {
			matrix->cells[cell].north = 1;
		}

		logger_error(g_debug_log, "importing met data...");
		matrix->cells[cell].years = import_meteo_data(g_sz_input_met_file, &years_of_simulation, matrix->cells[cell].x, matrix->cells[cell].y);
		if ( ! matrix->cells[cell].years ) goto err;
		logger_error(g_debug_log, "ok\n");

		/* set start year index */
		if ( -1 == g_year_start_index ) {
			int i;
			int ii;

			for ( i = 0; i < years_of_simulation; ++i ) {
				if ( g_settings->year_start == matrix->cells[0].years[i].year ) {
					g_year_start_index = i;
					break;
				}
			}
			if ( -1 == g_year_start_index ) {
				logger_error(g_debug_log, "start year (%d) not found. range is %d-%d\n"
						, g_settings->year_start
						, matrix->cells[0].years[0].year
						, matrix->cells[0].years[years_of_simulation-1].year
				);
				goto err;
			}

			/* set end year (adjusting years_of_simulation) */
			ii = -1;
			for ( i = 0; i < years_of_simulation; ++i ) {
				if ( g_settings->year_end == matrix->cells[0].years[i].year ) {
					if ( g_year_start_index > i ) {
						logger_error(g_debug_log, "start year (%d) cannot be > end year (%d)\n"
								, g_settings->year_start
								, g_settings->year_end
						);
						goto err;
					}
					ii = i;
					break;
				}
			}
			if ( -1 == ii ) {
				logger_error(g_debug_log, "end year (%d) not found. range is %d-%d\n"
						, g_settings->year_end
						, matrix->cells[0].years[0].year
						, matrix->cells[0].years[years_of_simulation-1].year
				);
				goto err;
			}
			/* i is year_end_index */
			i = ii + 1;

			years_of_simulation = i - g_year_start_index;
		}

		/* move pointer for year */
		matrix->cells[cell].years += g_year_start_index;

		logger(g_debug_log, "Total years_of_simulation = %d\n", years_of_simulation);
		logger(g_debug_log, "***************************************************\n\n");

		for ( year = 0; year < years_of_simulation; ++year )
		{
			/* for handling leap years */
			int days_per_month;

			for ( month = 0; month < MONTHS_COUNT; ++month )
			{
				days_per_month = DaysInMonth[month];
				if ( (FEBRUARY == month) && IS_LEAP_YEAR(matrix->cells[cell].years[year].year) )
				{
					++days_per_month;
				}

				for ( day = 0; day < days_per_month; ++day )
				{
					if( !day && !month && ! year )
					{
						/* general summary on model simulation */
						simulation_summary(matrix);

						/* summary on site */
						site_summary(matrix);

						/* summary on topography */
						topo_summary(matrix);

						/* summary on soil */
						soil_summary(matrix, matrix->cells);

						if( LANDUSE_F == g_soil_settings->landuse )
						{
							/* note: this happens just the first day of simulation */
							/* forest initialization */
							forest_initialization(matrix, day, month, year);
						}
						else
						{
							/* include summary for other land uses */
						}
					}

					/* compute daily climate variables not coming from met data */
					Avg_temperature(matrix->cells[cell].years[year].m, day, month);
					Daylight_avg_temperature(matrix->cells[cell].years[year].m, day, month);
					Nightime_avg_temperature(matrix->cells[cell].years[year].m, day, month);
					Soil_temperature(matrix->cells[cell].years[year].m, day, month);
					Thermic_sum(matrix->cells[cell].years[year].m, day, month, year);
					Air_density(matrix->cells[cell].years[year].m, day, month);
					Day_Length(&matrix->cells[cell], day, month, year);
					Latent_heat(matrix->cells[cell].years[year].m, day, month);
					Air_pressure(matrix->cells[cell].years[year].m, day, month);
					Psychrometric(matrix->cells[cell].years[year].m, day, month);
					Sat_vapour_pressure(&matrix->cells[cell], day, month, year);
					Dew_temperature(matrix->cells[cell].years[year].m, day, month);
					Radiation(&matrix->cells[cell], day, month, year);
					Check_prcp(&matrix->cells[cell], day, month, year);

					if ( LANDUSE_F == g_soil_settings->landuse )
					{
						/* compute annually the days for the growing season before any other process */
						Veg_Days (&matrix->cells[cell], day, month, year);
					}
					else
					{
						/* include other land use ???? */
					}
				}
			}
			for ( month = 0; month < MONTHS_COUNT; ++month )
			{
				days_per_month = DaysInMonth[month];
				if ( (FEBRUARY == month) && IS_LEAP_YEAR(matrix->cells[cell].years[year].year) )
				{
					++days_per_month;
				}
				for ( day = 0; day < days_per_month; ++day )
				{
					/* print daily met data */
					print_daily_met_data (&matrix->cells[cell], day, month, year);

					/* print cell data */
					print_daily_cell_data ( &matrix->cells[cell] );

					/************************************************************************/
					if ( g_sz_dataset_file )
					{
						if ( (LANDUSE_F == g_soil_settings->landuse) && (matrix->cells[cell].heights_count != 0) )
						{
							if ( 'f' == g_settings->version )
							{
								if ( !Tree_model_daily( matrix, cell, day, month, year ) )
								{
									logger(g_debug_log, "tree model daily failed!!!");
								}
								else
								{
									//printf("ok tree_model (%02d-%02d-%d)\n", day+1, month+1, year+g_settings->year_start);
								}
							}
							else
							{
								/* run for BGC version */
							}
						}
						else
						{
							//if ( ! crop_model_D(matrix, cell, year, month, day, years_of_simulation) )
							//{
							//	logger(g_debug_log, "crop model failed.");
							//}
							//else
							//{
							//	puts(msg_ok);
							//	//look if put it here or move before tree_model at the beginning of each month simulation
							//	soil_model (matrix, yos, years, month, years_of_simulation);
							//}
						}
					}

					/************************************************************************/
					/* run for soil model */
					if ( !Soil_model_daily(matrix, cell, day, month, year) )
					{
						logger_error(g_debug_log, "soil model daily failed!!!");
					}
					else
					{
						//printf("ok soil_model (%02d-%02d-%d)\n", day+1, month+1, year+g_settings->year_start);
					}
					/************************************************************************/
					/* run for cell model */
					if ( !Cell_model_daily(matrix, cell, day, month, year) )
					{
						logger_error(g_debug_log, "cell model daily failed!!!");
					}
					else
					{
						//printf("ok cell_model (%02d-%02d-%d)\n", day+1, month+1, year+g_settings->year_start);
					}
					/*************************************************************************/

					/******************************************************************************/
					/* print daily output */
					EOD_print_output_cell_level (&matrix->cells[cell], day, month, year, years_of_simulation );

					/* reset daily variables once printed */
					reset_daily_class_variables ( &matrix->cells[cell] );
					reset_daily_layer_variables ( &matrix->cells[cell] );
					reset_daily_cell_variables  ( &matrix->cells[cell] );

					//					if ( ! string_compare_i(g_settings->dndc, "on") )
					//					{
					//						Get_EOD_soil_balance_cell_level (&matrix->cells[cell], year, month, day);
					//					}
					logger(g_debug_log, "****************END OF DAY (%d)*******************\n\n\n", day + 1 );
				}

				/******************************************************************************/
				/* print monthly output */
				EOM_print_output_cell_level( &matrix->cells[cell], month, year, years_of_simulation );

				/* reset monthly variables once printed */
				if ( ( IS_LEAP_YEAR( matrix->cells[cell].years[year].year ) ? (MonthLength_Leap[month]) : (MonthLength[month] )) == matrix->cells[cell].doy )
				{
					reset_monthly_class_variables ( &matrix->cells[cell] );
					reset_monthly_layer_variables ( &matrix->cells[cell] );
					reset_monthly_cell_variables  ( &matrix->cells[cell] );
				}
				logger(g_debug_log, "****************END OF MONTH (%d)*******************\n\n", month + 1 );
			}

			/******************************************************************************/
			/* print annual output */
			EOY_print_output_cell_level( &matrix->cells[cell], year, years_of_simulation );
                        

			/* reset annual variables once printed */
			if ( ( IS_LEAP_YEAR( matrix->cells[cell].years[year].year ) ? (MonthLength_Leap[DECEMBER]) : (MonthLength[DECEMBER] )) == matrix->cells[cell].doy )
			{
				reset_annual_class_variables ( &matrix->cells[cell] );
				reset_annual_layer_variables ( &matrix->cells[cell] );
				reset_annual_cell_variables  ( &matrix->cells[cell] );
			}
			logger(g_debug_log, "****************END OF YEAR (%d)*******************\n", matrix->cells[cell].years[year].year );
		}

		if ( ! matrix->cells_count ) break;

		if ( g_year_start_index != -1 ) {
			free(matrix->cells[cell].years-g_year_start_index);
		} else {
			free(matrix->cells[cell].years);
		}
		matrix->cells[cell].years = NULL; /* required */
	}

	/* free memory */
	matrix_free(matrix); matrix = NULL;

	/* ok ! */
	prog_ret = 0;

#ifdef BENCHMARK_ONLY
	benchmark:
#endif

	// TODO: FIX THIS
	/* close logger */
	logger_close(g_soil_log); g_soil_log = NULL;
	logger_close(g_annual_log); g_annual_log = NULL;
	logger_close(g_monthly_log); g_monthly_log = NULL;
	logger_close(g_daily_log); g_daily_log = NULL;

	/* benchmark ? */
	if ( g_sz_benchmark_path ) {
		int i;
		char temp[256];
		i = has_path_delimiter(g_sz_output_path);
		sprintf(temp, "%s%soutput_%s_%s"
				, g_sz_output_path
				, i ? "" : FOLDER_DELIMITER
						, PROGRAM_VERSION
						, sz_date
		);
		logger_error(g_debug_log, "benchmark with '%s'...", g_sz_benchmark_path);
		if ( compare(temp, g_sz_benchmark_path) ) {
			puts("ok!");
		}
	}

	err:

	/* cleanup memory...
		do not remove null pointer to prevent
		double free with clean_up func
	 */

	/* close logger */
	logger_close(g_soil_log); g_soil_log = NULL;
	logger_close(g_annual_log); g_annual_log = NULL;
	logger_close(g_monthly_log); g_monthly_log = NULL;
	logger_close(g_daily_log); g_daily_log = NULL;

	/* free memory */
	if ( t ) free(t);
	if ( s ) free(s);
	if ( output_vars ) output_free(output_vars);
	if ( matrix ) matrix_free(matrix);

	free(g_sz_input_met_file); g_sz_input_met_file = NULL;
	free(g_sz_settings_file); g_sz_settings_file = NULL;

	/* this should be freed before */
	free(g_sz_topo_file); g_sz_topo_file = NULL;
	free(g_sz_soil_file); g_sz_soil_file = NULL;

	/* print elapsed time */
	end_timer = timer_get();
	timer = end_timer - start_timer;
	logger_error(g_debug_log, "%.2f secs / %.2f mins / %.2f hours elapsed\n", timer, timer / 60.0, timer / 3600.0);
	logger_close(g_debug_log); g_debug_log = NULL;

	return prog_ret;
}
#endif
