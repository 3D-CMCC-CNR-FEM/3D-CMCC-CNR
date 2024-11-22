/* settings.c */
#include "settings.h"
#include "common.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include "matrix.h"

enum {
	SETTINGS_SITENAME = 0
	, SETTINGS_VERSION
	, SETTINGS_SPATIAL
	, SETTINGS_TIME
	, SETTINGS_SPINUP
	, SETTINGS_SPINUP_YEARS
	, SETTINGS_SCREEN_OUTPUT
	, SETTINGS_DEBUG_OUTPUT
	, SETTINGS_DAILY_OUTPUT
	, SETTINGS_MONTHLY_OUTPUT
	, SETTINGS_YEARLY_OUTPUT
	, SETTINGS_YEARLY_MAN_OUTPUT
	, SETTINGS_SOIL_OUTPUT
	, SETTINGS_NETCDF_OUTPUT
	, SETTINGS_YEAR_START
	, SETTINGS_YEAR_END
	, SETTINGS_YEAR_RESTART
	, SETTINGS_PSN_MOD
	, SETTINGS_CO2_TRANS
	, SETTINGS_YEAR_START_CO2_FIXED
	, SETTINGS_NDEP_FIXED
	, SETTINGS_TBASE_RESP        // note: currently this setting is not available in the ifle  OPTIONAL 
	, SETTINGS_PHOTO_ACCL
	, SETTINGS_RESP_ACCL
	, SETTINGS_COLD_ACCL         // 5p67
	, SETTINGS_N_REG             // 5p67
	, SETTINGS_REGENERATION
	, SETTINGS_MANAGEMENT
	, SETTINGS_MANAGEMENT_TYPE   // 5p7 this value/setting is actually not reported in the setting file OPTIONAL 
	, SETTINGS_YEAR_START_MANAGEMENT
	, SETTINGS_PROGN_AUT_RESP
	, SETTINGS_SIZECELL
	, SETTINGS_Y
	, SETTINGS_CO2CONC
	, SETTINGS_CO2_INCR
	, SETTINGS_INIT_FRAC_MAXASW
	, SETTINGS_TREE_LAYER_LIMIT
	, SETTINGS_SOIL_LAYER
	, SETTINGS_MAX_LAYER_COVER
	, SETTINGS_THINNING_REGIME
	, SETTINGS_BRANCH_REM       // 5p67
	, SETTINGS_REGENERATION_SPECIES
	, SETTINGS_REGENERATION_MANAGEMENT
	, SETTINGS_REGENERATION_N_TREE
	, SETTINGS_REGENERATION_AGE
	, SETTINGS_REGENERATION_AVDBH
	, SETTINGS_REGENERATION_LAI
	, SETTINGS_REGENERATION_HEIGHT
	, SETTINGS_PRUNING
	, SETTINGS_IRRIGATION

	, SETTINGS_COUNT
};

enum
{
	SETTINGS_REPLANTED_SPECIES
	, SETTINGS_REPLANTED_MANAGEMENT
	, SETTINGS_REPLANTED_TREE
	, SETTINGS_REPLANTED_AGE
	, SETTINGS_REPLANTED_AVDBH
	, SETTINGS_REPLANTED_LAI
	, SETTINGS_REPLANTED_HEIGHT

	, REPLANTED_SETTINGS_COUNT
};

typedef struct
{
	replanted_vars_t values;
	int assigned_values[REPLANTED_SETTINGS_COUNT];
	int index;

} replanted_temp_t;

extern char* g_sz_input_path;
extern const char sz_err_out_of_memory[];

const char* sz_settings[SETTINGS_COUNT] = {
	"SITENAME"
	, "VERSION"
	, "SPATIAL"
	, "TIME"
	, "SPINUP"
	, "SPINUP_YEARS"
	, "SCREEN_OUTPUT"
	, "DEBUG_OUTPUT"
	, "DAILY_OUTPUT"
	, "MONTHLY_OUTPUT"
	, "ANNUAL_OUTPUT"
	, "ANNUAL_MAN_OUTPUT"
	, "SOIL_OUTPUT"
	, "NETCDF_OUTPUT"
	, "YEAR_START"
	, "YEAR_END"
	, "YEAR_RESTART"
	, "PSN_MOD"
	, "CO2_TRANS"
	, "YEAR_START_CO2_FIXED"
	, "NDEP_FIXED"
	, "TBASE_RESP"
	, "PHOTO_ACCL"
	, "RESP_ACCL"
	, "COLD_ACCL"     // 5p7
	, "N_REG"         // 5p7
	, "REGENERATION"
	, "MANAGEMENT"
	, "MANAGEMENT_TYPE"
	, "YEAR_START_MANAGEMENT"
	, "PROGN_AUT_RESP"
	, "SIZECELL"
	, "Y"
	, "CO2CONC"
	, "CO2_INCR"
	, "INIT_FRAC_MAXASW"
	, "TREE_LAYER_LIMIT"
	, "SOIL_LAYER"
	, "MAX_LAYER_COVER"
	, "THINNING_REGIME"
	,  "BRANCH_REM"   // 5p7 percentage of branches left in the stand after thinning 
	, "REGENERATION_SPECIES"
	, "REGENERATION_MANAGEMENT"
	, "REGENERATION_N_TREE"
	, "REGENERATION_AGE"
	, "REGENERATION_AVDBH"
	, "REGENERATION_LAI"
	, "REGENERATION_HEIGHT"
	, "PRUNING"
	, "IRRIGATION"
};

const int optional[] = {
	SETTINGS_YEAR_RESTART
	, SETTINGS_MANAGEMENT_TYPE
	, SETTINGS_REGENERATION_SPECIES
	, SETTINGS_TBASE_RESP
};

static const char* sz_replanted_settings[REPLANTED_SETTINGS_COUNT] =
{
	"REPLANTED_SPECIES"
	, "REPLANTED_MANAGEMENT"
	, "REPLANTED_TREE"
	, "REPLANTED_AGE"
	, "REPLANTED_AVDBH"
	, "REPLANTED_LAI"
	, "REPLANTED_HEIGHT"
};

static int realloc_replanted_temp(replanted_temp_t** t, int* count) {
	replanted_temp_t *no_leak;

	no_leak = realloc(*t, ++*count*sizeof*no_leak);
	if ( ! no_leak )
	{
		--*count;
		return 0;
	}
	*t = no_leak;
	return 1;
}

// index is zero based
static replanted_temp_t* get_current_by_index(replanted_temp_t* p, int count, int index)
{
	int i;

	for ( i = 0; i < count; i++ )
	{
		if ( index == p[i].index )
		{
			return &p[i];
		}
	}
	return NULL;
}

static int settings_replanted_import(const char* const filename, settings_t* s)
{
#define BUFFER_SIZE	256

	char buffer[BUFFER_SIZE];
	int i;
	//int imported[REPLANTED_SETTINGS_COUNT] = { 0 };
	replanted_temp_t* temp;
	int temp_count;
	
	replanted_vars_t replanted = { 0 };
	FILE* f;

	const char delimiter[] = " /\"\t\r\n";

	assert(filename && s);

	temp = NULL;
	temp_count = 0;

	f = fopen(filename, "r");
	if ( ! f ) {
		puts("file not found!");
		return 0;
	}

	while ( fgets(buffer, BUFFER_SIZE, f) )
	{
		char *p;
		char *p2;
		char *token;
		int index;
		replanted_temp_t* current;

		/* remove initial spaces (if any) */
		p2 = buffer;
		while ( isspace(*p2) ) ++p2;

		/* skip empty lines and comments */
		if (	('\r' == p2[0])
				|| ('\n' == p2[0])
				|| ('/' == p2[0]) ) {
			continue;
		}

		/* get setting name */
		current = NULL;
		token = string_tokenizer(p2, delimiter, &p);
		if ( ! token ) continue;

		/* check for name */
		for ( i = 0; i < REPLANTED_SETTINGS_COUNT; i++ )
		{
			if ( ! string_compare_i(token, sz_replanted_settings[i]) )
			{ 
				// token found, without index

				current = get_current_by_index(temp, temp_count, 0);
				if ( ! current )
				{
					// add new
					replanted_temp_t* no_leak;
					no_leak = realloc(temp, (temp_count+1)*sizeof*no_leak);
					if ( ! no_leak )
					{
						printf("unable to retrieve index from %s\n", token);
						if ( temp ) free(temp);
						fclose(f);
						return 0;
					}
					temp = no_leak;
					current = &temp[temp_count++];
					current->index = 0;
					{
						int y;
						for ( y = 0; y < REPLANTED_SETTINGS_COUNT; y++ )
						{
							current->assigned_values[y] = 0;
						}
					}
				}
				break;
			}
			else
			{
				// if not found, look for the name of variable with index
				// e.g: REPLANTED_SPECIES_0  and so on
				p2 = strstr(token, sz_replanted_settings[i]);
				if ( 0 == p2 - token )
				{
					// get index
					if ( 1 != sscanf(token, "%*[^0-9]%d", &index) )
					{
						printf("unable to retrieve index from %s\n", token);
						if ( temp ) free(temp);
						fclose(f);
						return 0;
					}
					if ( index < 0 )
					{
						printf("bad index specified for %s\n", token);
						if ( temp ) free(temp);
						fclose(f);
						return 0;
					}
					if ( index > temp_count )
					{
						printf("bad index specified for %s: should be %d\n", token, temp_count);
						if ( temp ) free(temp);
						fclose(f);
						return 0;
					}
					current = get_current_by_index(temp, temp_count, index);
					if ( ! current )
					{
						// add new
						replanted_temp_t* no_leak;
						no_leak = realloc(temp, (temp_count+1)*sizeof*no_leak);
						if ( ! no_leak )
						{
							printf("unable to retrieve index from %s\n", token);
							if ( temp ) free(temp);
							fclose(f);
							return 0;
						}
						temp = no_leak;
						current = &temp[temp_count++];
						current->index = index;
						{
							int y;
							for ( y = 0; y < REPLANTED_SETTINGS_COUNT; y++ )
							{
								current->assigned_values[y] = 0;
							}
						}
					}
					break;
				}
			}
		}
		
		if ( REPLANTED_SETTINGS_COUNT == i )
		{
			continue;
		}

		if ( ! current )
		{
			fclose(f);
			return 0;
		}

		if ( current->assigned_values[i] )
		{
			if ( current->index )
				printf("%s_%d already specified\n", sz_replanted_settings[i], current->index);
			else
				printf("%s already specified\n", sz_replanted_settings[i]);
			free(temp);
			fclose(f);
			return 0;
		}

		/* get value */
		token = string_tokenizer(NULL, delimiter, &p);
		if ( ! token ) {
			if ( current->index )
				printf("no value specified for %s_%d\n", sz_replanted_settings[i], current->index);
			else
				printf("no value specified for %s\n", sz_replanted_settings[i]);
			free(temp);
			fclose(f);
			return 0;
		}

		if ( SETTINGS_REPLANTED_SPECIES == i )
		{
			strncpy(current->values.species, (const char*)token, SETTINGS_REPLANTED_SPECIES_MAX_SIZE-1);
		}
		else if ( i == SETTINGS_REPLANTED_MANAGEMENT )
		{
			if ( ('T' == token[0]) || ('t' == token[0]) ) {
				current->values.management = T;
			} else if ( ('C' == token[0]) || ('c' == token[0]) ) {
				current->values.management = C;
			} else if ( ('N' == token[0]) || ('n' == token[0]) ) {
				current->values.management = N;
			} else {
				printf("bad management habitus specified: %s\n", token);
				free(temp);
				fclose(f);
				return 0;
			}
		} else {
			int err;
			double value;

			value = convert_string_to_float(token, &err);
			if ( err )
			{
				printf("unable to convert value for %s: \"%s\"\n", sz_replanted_settings[i], token);
				free(temp);
				fclose(f);
				return 0;
			}
			switch ( i )
			{
				case SETTINGS_REPLANTED_TREE:
					current->values.n_tree = value;
				break;

				case SETTINGS_REPLANTED_AGE:
					current->values.age = value;
				break;

				case SETTINGS_REPLANTED_AVDBH:
					current->values.avdbh = value;
				break;

				case SETTINGS_REPLANTED_LAI:
					current->values.lai = value;
				break;

				case SETTINGS_REPLANTED_HEIGHT:
					current->values.height = value;
				break;
			}
		}
		current->assigned_values[i] = 1;
	}
	fclose(f);

	// parse imported stuff
	if ( ! temp_count )
	{
		puts("no replanted stuff found!");
		return 0;
	}

	// imported all stuff ?
	for  ( i = 0; i < temp_count; i++ )
	{
		int y;
		int flag = 0;
		for  ( y = 0; y < REPLANTED_SETTINGS_COUNT; y++ )
		{
			if ( ! temp[i].assigned_values[y] )
			{
				if ( temp[i].index )
				{
					printf("%s_%d not specified\n", sz_replanted_settings[y], temp[i].index);
				}
				else
				{
					printf("%s not specified\n", sz_replanted_settings[y]);
				}
				free(temp);
				return 0;
			}
		}
	}

	// assign values
	s->replanted_count = temp_count;
	s->replanted = malloc(s->replanted_count*sizeof*s->replanted);
	if ( ! s->replanted )
	{
		puts("out of memory!");
		free(temp);
		return 0;
	}

    


	for ( i = 0; i < s->replanted_count; i++ )
	{
		strcpy(s->replanted[i].species, temp[i].values.species);
		s->replanted[i].management = temp[i].values.management;
		s->replanted[i].n_tree = temp[i].values.n_tree;
		s->replanted[i].age = temp[i].values.age;
		s->replanted[i].avdbh = temp[i].values.avdbh;
		s->replanted[i].lai = temp[i].values.lai;
		s->replanted[i].height = temp[i].values.height;
	}

	free(temp);

	return 1;

#undef BUFFER_SIZE
}

static int setting_pruning_import(settings_t* s) {
#define PRUNING_FILENAME "pruning.csv"
	char* path;
	int ret;
	int rows_count;
	int prunings_count;
	FILE *f;
	pruning_t* prunings;

	assert(s);
	// must be empty!
	assert(!s->prunings);

	ret = 0; // defaults to err;
	path = NULL;
	f = NULL;

	// compute path
	if ( g_sz_input_path ) {
		path = concatenate_path(g_sz_input_path, PRUNING_FILENAME);
	} else {
		path = string_copy(PRUNING_FILENAME);
	}
	if ( ! path ) {
		puts(sz_err_out_of_memory);
		goto quit;
	}

	// TODO fix path
	// returns -1 if file not found
	// returns -2 if out of memory
	// returns -3 if read error
	rows_count = file_get_rows_count(path);
	if ( rows_count <= 0 ) {
		printf("%s", path);
		if (  0 == rows_count ) puts(" is empty");
		if ( -1 == rows_count ) puts(" not found");
		if ( -2 == rows_count ) printf(": %s\n", sz_err_out_of_memory);
		if ( -3 == rows_count ) puts(" read error");
		goto quit;
	}

	// remove header
	if ( ! --rows_count ) {
		printf("%s is empty\n", path);
		goto quit;
	}

	s->prunings = malloc(rows_count*sizeof*prunings);
	if ( ! s->prunings ) {
		puts(sz_err_out_of_memory);
		goto quit;
	}

	// import file
	{
	#define BUFFER_SIZE	256
	#define COLUMNS_COUNT 5
		char buffer[BUFFER_SIZE];
		char* token;
		char* p;
		int row;
		int imported_rows;

		const char header[] = "YEAR,MONTH,DAY,INTESITY,HEIGHT";

		f = fopen(path, "r");
		if ( ! f ) {
			printf("unable to open %s\n", path);
			goto quit;
		}

		// get header
		if ( ! fgets(buffer, BUFFER_SIZE, f) ) {
			printf("%s is empty\n", path);
			goto quit;
		}

		// remove endline
		for ( row = 0; buffer[row]; ++row ) {
			if ( ('\r' == buffer[row]) || ('\n' == buffer[row]) ) {
				buffer[row] = '\0';
				break;
			}
		}

		// check if header is compliant
		if ( string_compare_i(buffer, header) ) {
			printf("%s has an invalid header, must be: %s\n", path, header);
			goto quit;
		}

		row = -1;
		imported_rows = 0;
		prunings_count = 0;
		while ( fgets(buffer, BUFFER_SIZE, f) ) {
			int i;
			int imported_columns;

			const char delimiter[] = ",";

			// on crash there's an error during rows_count computation!
			// ( file_get_rows_count func inside common.c )
			++row; // DO NOT INCREMENT INSIDE CHECK CONDITION
			CHECK_CONDITION(row, >=, rows_count);

			// remove endlines
			for ( i = 0; buffer[i]; ++i ) {
				if ( ('\r' == buffer[i]) || ('\n' == buffer[i]) ) {
					buffer[i] = '\0';
					break;
				}
			}
			// skip empty row
			if ( ! buffer[0] )
				continue;

			imported_columns = 0;
			for ( i = 0, token = string_tokenizer(buffer, delimiter, &p); token; token = string_tokenizer(NULL, delimiter, &p), ++i ) {
				int err;
				double v;

				if ( i >= COLUMNS_COUNT ) {
					printf("%s has too many columns at row %d\n", path, row+1);
					goto quit;
				}

				v = convert_string_to_float(token, &err);
				if ( err ) {
					//printf("%s: unable to convert value '%s' at row %d, column %d\n", token, row+1, i+1);
					//format correction
					 printf(" unable to convert value %s at row %d, column %d\n", token, row+1, i+1);
				
					goto quit;
				}

				switch ( i )
				{
					// year
					case 0:
						s->prunings[row].year = (int)v;
					break;

					// month
					case 1:
						s->prunings[row].month = (int)v;
					break;

					// day
					case 2:
						s->prunings[row].day = (int)v;
					break;

					// intensity
					case 3:
						s->prunings[row].intensity = v;
					break;

					// height
					case 4:
						s->prunings[row].height = v;
					break;
				}
				++imported_columns;
			}

			if ( imported_columns != COLUMNS_COUNT ) {
				printf("%s has %d columns instead of %d columns at row %d\n", path, imported_columns, COLUMNS_COUNT, row+1);
				goto quit;
			}
			++imported_rows;
		}

		s->prunings_count = imported_rows;

	#undef COLUMNS_COUNT
	#undef BUFFER_SIZE
	}

	if ( ! s->prunings_count ) {
		printf("%s is empty\n", path);
		goto quit;
	}

	// check valid date inside struct
	{
		int i;

		const int days_per_month[] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };

		for ( i = 0; i < s->prunings_count; ++i )
		{
			int year;
			int month;
			int day;
			int days;
						
			year = s->prunings[i].year;
			if ( year < 0 ) {
				printf("%s has invalid year (%d) at row %d\n", path, year, i+1+1); // we start from zero and add header too
				goto quit;
			}
			
			month = s->prunings[i].month;
			if ( (month < 1) || (month > 12) ) {
				printf("%s has invalid month (%d) at row %d\n", path, month, i+1+1); // we start from zero and add header too
				goto quit;
			}

			days = days_per_month[month-1];
			// fix for leap year
			if ( (2 == month) && IS_LEAP_YEAR(year) ) {
				++days;
			}
			day = s->prunings[i].day;
			if ( (day < 1) || (day > days) ) {
				printf("%s has invalid day (%d) at row %d\n", path, day, i+1+1); // we start from zero and add header too
				goto quit;
			}
		}
	}

	ret = 1; // OK

quit:
	if ( f )
		fclose(f);

	if ( path )
		free(path);

	if ( ! ret ) {
		s->prunings_count = 0;
		if ( s->prunings ) {
			free(s->prunings);
			s->prunings = NULL;
		}
	}

	return ret;
#undef PRUNING_FILENAME
}

static int setting_irrigation_import(settings_t* s) {
#define IRRIGATION_FILENAME "irrigation.csv"
	char* path;
	int ret;
	int rows_count;
	int irrigations_count;
	FILE *f;
	irrigation_t* irrigations;

	assert(s);
	// must be empty!
	assert(!s->irrigations);

	ret = 0; // defaults to err;
	path = NULL;
	f = NULL;

	// compute path
	if ( g_sz_input_path ) {
		path = concatenate_path(g_sz_input_path, IRRIGATION_FILENAME);
	} else {
		path = string_copy(IRRIGATION_FILENAME);
	}
	if ( ! path ) {
		puts(sz_err_out_of_memory);
		goto quit;
	}

	// TODO fix path
	// returns -1 if file not found
	// returns -2 if out of memory
	// returns -3 if read error
	rows_count = file_get_rows_count(path);
	if ( rows_count <= 0 ) {
		printf("%s", path);
		if (  0 == rows_count ) puts(" is empty");
		if ( -1 == rows_count ) puts(" not found");
		if ( -2 == rows_count ) printf(": %s\n", sz_err_out_of_memory);
		if ( -3 == rows_count ) puts(" read error");
		goto quit;
	}

	// remove header
	if ( ! --rows_count ) {
		printf("%s is empty\n", path);
		goto quit;
	}

	s->irrigations = malloc(rows_count*sizeof*irrigations);
	if ( ! s->irrigations ) {
		puts(sz_err_out_of_memory);
		goto quit;
	}

	// import file
	{
	#define BUFFER_SIZE	256
	#define COLUMNS_COUNT 4
		char buffer[BUFFER_SIZE];
		char* token;
		char* p;
		int row;
		int imported_rows;

		const char header[] = "YEAR,MONTH,DAY,AMOUNT";

		f = fopen(path, "r");
		if ( ! f ) {
			printf("unable to open %s\n", path);
			goto quit;
		}

		// get header
		if ( ! fgets(buffer, BUFFER_SIZE, f) ) {
			printf("%s is empty\n", path);
			goto quit;
		}

		// remove endline
		for ( row = 0; buffer[row]; ++row ) {
			if ( ('\r' == buffer[row]) || ('\n' == buffer[row]) ) {
				buffer[row] = '\0';
				break;
			}
		}

		// check if header is compliant
		if ( string_compare_i(buffer, header) ) {
			printf("%s has an invalid header, must be: %s\n", path, header);
			goto quit;
		}

		row = -1;
		imported_rows = 0;
		irrigations_count = 0;
		while ( fgets(buffer, BUFFER_SIZE, f) ) {
			int i;
			int imported_columns;

			const char delimiter[] = ",";

			// on crash there's an error during rows_count computation!
			// ( file_get_rows_count func inside common.c )
			++row; // DO NOT INCREMENT INSIDE CHECK CONDITION
			CHECK_CONDITION(row, >=, rows_count);

			// remove endlines
			for ( i = 0; buffer[i]; ++i ) {
				if ( ('\r' == buffer[i]) || ('\n' == buffer[i]) ) {
					buffer[i] = '\0';
					break;
				}
			}
			// skip empty row
			if ( ! buffer[0] )
				continue;

			imported_columns = 0;
			for ( i = 0, token = string_tokenizer(buffer, delimiter, &p); token; token = string_tokenizer(NULL, delimiter, &p), ++i ) {
				int err;
				double v;

				if ( i >= COLUMNS_COUNT ) {
					printf("%s has too many columns at row %d\n", path, row+1);
					goto quit;
				}

				v = convert_string_to_float(token, &err);
				if ( err ) {
					//printf("%s: unable to convert value '%s' at row %d, column %d\n", token, row+1, i+1);
					// format correction
					printf(" unable to convert value %s at row %d, column %d\n", token, row+1, i+1);
					goto quit;
				}

				switch ( i )
				{
					// year
					case 0:
						s->irrigations[row].year = (int)v;
					break;

					// month
					case 1:
						s->irrigations[row].month = (int)v;
					break;

					// day
					case 2:
						s->irrigations[row].day = (int)v;
					break;

					// amount
					case 3:
						s->irrigations[row].amount = v;
					break;
				}
				++imported_columns;
			}

			if ( imported_columns != COLUMNS_COUNT ) {
				printf("%s has %d columns instead of %d columns at row %d\n", path, imported_columns, COLUMNS_COUNT, row+1);
				goto quit;
			}

			++imported_rows;
		}

		s->irrigations_count = imported_rows;

	#undef COLUMNS_COUNT
	#undef BUFFER_SIZE
	}

	if ( ! s->irrigations_count ) {
		printf("%s is empty\n", path);
		goto quit;
	}

	// check valid date inside struct
	{
		int i;

		const int days_per_month[] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };

		for ( i = 0; i < s->irrigations_count; ++i )
		{
			int year;
			int month;
			int day;
			int days;
						
			year = s->irrigations[i].year;
			if ( year < 0 ) {
				printf("%s has invalid year (%d) at row %d\n", path, year, i+1+1); // we start from zero and add header too
				goto quit;
			}
			
			month = s->irrigations[i].month;
			if ( (month < 1) || (month > 12) ) {
				printf("%s has invalid month (%d) at row %d\n", path, month, i+1+1); // we start from zero and add header too
				goto quit;
			}

			days = days_per_month[month-1];
			// fix for leap year
			if ( (2 == month) && IS_LEAP_YEAR(year) ) {
				++days;
			}
			day = s->irrigations[i].day;
			if ( (day < 1) || (day > days) ) {
				printf("%s has invalid day (%d) at row %d\n", path, day, i+1+1); // we start from zero and add header too
				goto quit;
			}
		}
	}

	ret = 1; // OK

quit:
	if ( f )
		fclose(f);

	if ( path )
		free(path);

	if ( ! ret ) {
		s->irrigations_count = 0;
		if ( s->irrigations ) {
			free(s->irrigations);
			s->irrigations = NULL;
		}
	}

	return ret;
#undef IRRIGATION_FILENAME
}

settings_t* settings_import(const char *const filename) {
#define BUFFER_SIZE	256
	char buffer[BUFFER_SIZE];
	int i;
	int index;
	int err;
	settings_t* s;
	int pruning;
	int irrigation;
	int imported[SETTINGS_COUNT] = { 0 };
	FILE *f;

	const char delimiter[] = " /\"\t\r\n";

	f = fopen(filename, "r");
	if ( ! f ) {
		puts("file not found!");
		return 0;
	}

	s = malloc(sizeof*s);
	if ( ! s  ) {
		puts(sz_err_out_of_memory);
		fclose(f);
		return 0;
	}
	/* all settings values defaults to 0 ( off ) */
	memset(s, 0, sizeof*s);
	memset(imported, 0, (sizeof*imported)*SETTINGS_COUNT);

	pruning = 0;
	irrigation = 0;
	while ( fgets(buffer, BUFFER_SIZE, f) ) {
		char *p;
		char *p2;
		char *token;
		double value;

		/* remove initial spaces (if any) */
		p2 = buffer;
		while ( isspace(*p2) ) ++p2;

		/* skip empty lines and comments */
		if (	('\r' == p2[0])
				|| ('\n' == p2[0])
				|| ('/' == p2[0]) ) {
			continue;
		}

		/* get setting name */
		token = string_tokenizer(p2, delimiter, &p);
		if ( ! token ) continue;

		/* check for name */
		index = -1;
		for ( i = 0; i < SETTINGS_COUNT; i++ ) {
			if ( ! string_compare_i(token, sz_settings[i]) ) {
				index = i;
				break;
			}
		}
		// not found ?
		if ( -1 == index )
		{
			//printf("unknown parameter specified in settings: %s. skipped.\n", token);
			continue;
		}

		if ( imported[index] ) {
			printf("%s already specified. skipped.\n", token);
			continue;
		}

		/* get value */
		token = string_tokenizer(NULL, delimiter, &p);
		if ( ! token ) {
			int flag;

			printf("no value specified for %s\n", sz_settings[index]);
			/* check if is an optional parameter */
			flag = 0;
			for ( i = 0; i < (int)SIZE_OF_ARRAY(optional); ++i ) {
				if ( index == optional[i] ) {
					flag = 1;
					break;
				}
			}
			if ( ! flag ) {
				settings_free(s);
				fclose(f);
				return NULL;
			} else {
				puts(". skipped.");
				continue;
			}
		}

		switch ( index ) {
			case SETTINGS_SITENAME:
				strncpy(s->sitename, token, SETTINGS_SITENAME_MAX_SIZE-1);
			break;

			case SETTINGS_VERSION:
				s->version = *token;
			break;

			case SETTINGS_SPATIAL:
				s->spatial = *token;
			break;

			case SETTINGS_TIME:
				s->time = *token;
				if ( s->time != 'd' ) {
					puts("uncorrect time step choiced!");
					settings_free(s);
					fclose(f);
					return 0;
				}
			break;

			case SETTINGS_SPINUP:
				if ( ! string_compare_i(token, "on") ) {
					s->spinup = 1;
				}
			break;

			case SETTINGS_SCREEN_OUTPUT:
				if ( ! string_compare_i(token, "on") ) {
					s->screen_output = 1;
				}
			break;

			case SETTINGS_DEBUG_OUTPUT:
				if ( ! string_compare_i(token, "on") ) {
					s->debug_output = 1;
				}
			break;

			case SETTINGS_DAILY_OUTPUT:
				if ( ! string_compare_i(token, "on") ) {
					s->daily_output = 1;
				}
			break;

			case SETTINGS_MONTHLY_OUTPUT:
				if ( ! string_compare_i(token, "on") ) {
					s->monthly_output = 1;
				}
			break;

			case SETTINGS_YEARLY_OUTPUT:
				if ( ! string_compare_i(token, "on") ) {
					s->yearly_output = 1;
				}
			break;

			case SETTINGS_YEARLY_MAN_OUTPUT:
				if ( ! string_compare_i(token, "on") ) {
					s->yearly_man_output = 1;
				}
			break;

			case SETTINGS_SOIL_OUTPUT:
				if ( ! string_compare_i(token, "on") ) {
					s->soil_output = 1;
				}
			break;

			case SETTINGS_NETCDF_OUTPUT:
				if ( ! string_compare_i(token, "on") ) {
					s->netcdf_output = 1;
				}
			break;

			case SETTINGS_CO2_TRANS:
				if ( ! string_compare_i(token, "on") ) {
					s->CO2_trans = CO2_TRANS_ON;
				} else if ( ! string_compare_i(token, "var") ) {
					s->CO2_trans = CO2_TRANS_VAR;
				}
			break;

			case SETTINGS_NDEP_FIXED:
				if ( ! string_compare_i(token, "on") ) {
					s->Ndep_fixed = 1;
				}
			break;

			case SETTINGS_PHOTO_ACCL:
				if ( ! string_compare_i(token, "on") ) {
					s->Photo_accl = 1;
				}
			break;

			case SETTINGS_RESP_ACCL:
				if ( ! string_compare_i(token, "on") ) {
					s->Resp_accl = 1;
				}
			break;
			
			//5p7 
			
			case SETTINGS_COLD_ACCL:
				if ( ! string_compare_i(token, "on") ) {
				s->Cold_accl = 1;
				}
		    break;
        
		    //5p7 
			case SETTINGS_N_REG:
				if ( ! string_compare_i(token, "on") ) {
				s->N_reg = 1;
				}
		    break;
			

			case SETTINGS_REGENERATION:
				if ( ! string_compare_i(token, "on") ) {
					s->regeneration = 1;
				}
			break;
                                 
			case SETTINGS_MANAGEMENT:
				if ( ! string_compare_i(token, "on") ) {
					s->management = MANAGEMENT_ON;
				} else if ( ! string_compare_i(token, "var") ) {
					s->management = MANAGEMENT_VAR;
				} else if ( ! string_compare_i(token, "var1") ) {
					s->management = MANAGEMENT_VAR1;
				}
			break;
                             
			case SETTINGS_PROGN_AUT_RESP:
				if ( ! string_compare_i(token, "on") ) {
					s->Prog_Aut_Resp = 1;
				}
			break;

			case SETTINGS_THINNING_REGIME:
				if ( ! string_compare_i(token, "below") ) {
					s->thinning_regime = THINNING_REGIME_BELOW;
				}
			break;

			case SETTINGS_REGENERATION_SPECIES:
				strncpy(s->regeneration_species, (const char*)token, SETTINGS_REGENERATION_SPECIES_MAX_SIZE-1);
			break;

			case SETTINGS_REGENERATION_MANAGEMENT:
				if ( ('T' == token[0]) || ('t' == token[0]) ) {
					s->regeneration_management = T;
				} else if ( ('C' == token[0]) || ('c' == token[0]) ) {
					s->regeneration_management = C;
				} else if ( ('N' == token[0]) || ('n' == token[0]) ) {
					s->regeneration_management = N;
				} else {
					printf("bad regeneration habitus specified: %s\n", token);
					settings_free(s);
					fclose(f);
					return 0;
				}
			break;
                                  
			case SETTINGS_PRUNING:
				if ( ! string_compare_i(token, "on") ) {
					pruning = 1;
				}
			break;

			case SETTINGS_IRRIGATION:
				if ( ! string_compare_i(token, "on") ) {
					irrigation = 1;
				}
			break;

			default:
				value = convert_string_to_float(token, &err);
				if ( err ) {
					if ( ! ((SETTINGS_YEAR_RESTART == index) && ! string_compare_i(token, "off")) ) {
						printf("unable to convert value for %s: \"%s\"\n", sz_settings[index], token);
						free(s);
						fclose(f);
						return 0;
					}
				}
				switch ( index ) {
					case SETTINGS_SPINUP_YEARS:
						s->spinup_years = (int)value;
					break;

					case SETTINGS_YEAR_START:
						s->year_start = (int)value;
					break;

					case SETTINGS_YEAR_END:
						s->year_end = (int)value;
					break;

					case SETTINGS_YEAR_RESTART:
						s->year_restart = (int)value;
					break;

					case SETTINGS_PSN_MOD:
						s->PSN_mod = (int)value;
					break;

					case SETTINGS_YEAR_START_CO2_FIXED:
						s->year_start_co2_fixed = (int)value;
					break;

					case SETTINGS_TBASE_RESP:
						s->Tbase_resp = value;
					break;

					case SETTINGS_YEAR_START_MANAGEMENT:
						s->year_start_management = (int)value;
					break;

					case SETTINGS_MANAGEMENT_TYPE:
						s->management_type = (int)value;
					break;

					case SETTINGS_SIZECELL:
						s->sizeCell = value;
						// ALESSIOC fill me
						s->sizeCell *= s->sizeCell;
					break;
					
					case SETTINGS_Y:
						s->Fixed_Aut_Resp_rate = value;
					break;

					case SETTINGS_CO2CONC:
						s->co2Conc = value;
					break;

					case SETTINGS_CO2_INCR:
						s->co2_incr = value;
					break;

					case SETTINGS_INIT_FRAC_MAXASW:
						s->init_frac_maxasw = value;
					break;

					case SETTINGS_TREE_LAYER_LIMIT:
						s->tree_layer_limit = value;
					break;

					case SETTINGS_SOIL_LAYER:
						s->number_of_soil_layer = value;
					break;

					case SETTINGS_MAX_LAYER_COVER:
						s->max_layer_cover = value;
					break;
                    //5p7 
					case SETTINGS_BRANCH_REM:
						s->thinning_branch_rem = (int)value;
					break;

					case SETTINGS_REGENERATION_N_TREE:
						s->regeneration_n_tree = value;
					break;

					case SETTINGS_REGENERATION_AGE:
						s->regeneration_age = value;
					break;

					case SETTINGS_REGENERATION_AVDBH:
						s->regeneration_avdbh = value;
					break;

					case SETTINGS_REGENERATION_LAI:
						s->regeneration_lai = value;
					break;

					case SETTINGS_REGENERATION_HEIGHT:
						s->regeneration_height = value;
					break;
				}
		}
		imported[index] = 1;
	}
	fclose(f);

	/* check if all setting values are imported */
	for (  i = 0; i < SETTINGS_COUNT; i++ ) {
		if ( ! imported[i] ) {
			int y;
			int flag;

			/* check for optional parameter */
			flag = 0;
			for ( y = 0; y < (int)SIZE_OF_ARRAY(optional); ++y ) {
				if ( i == optional[y] ) {
					flag = 1;
					break;
				}
			}

			printf("%s was not specified.", sz_settings[i]);

			if ( ! flag ) {
				puts("");
				settings_free(s);
				return 0;
			} else {
				puts(" optional.");
			}
		}
	}

       

	/* check for restart year */
	if ( 0 == s->year_restart ) {
		s->year_restart = -1;
	} else if ( (s->year_restart <= s->year_start) || (s->year_restart >= s->year_end) ) {
		printf("%s must be between %d and %d not %d\n", sz_settings[SETTINGS_YEAR_RESTART]
					, s->year_start+1
					, s->year_end-1
					, s->year_restart
		);
		settings_free(s);
		return NULL;
	}

	if ( ! settings_replanted_import(filename, s) )
	{
		settings_free(s);
		s = NULL;
	}

	if ( pruning && ! setting_pruning_import(s) )
	{
		settings_free(s);
		s = NULL;
	}
	else if ( irrigation && ! setting_irrigation_import(s) )
	{
		settings_free(s);
		s = NULL;
	}

	return s;
#undef BUFFER_SIZE
}

void settings_free(settings_t* s) {
	if ( s )
	{
		if ( s->irrigations_count )
		{
			free(s->irrigations);
		}

		if ( s->prunings_count )
		{
			free(s->prunings);
		}

		if ( s->replanted_count )
		{
			free(s->replanted);
		}

		free(s);
	}
}
