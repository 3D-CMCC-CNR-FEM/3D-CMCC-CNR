/* logger.c */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <assert.h>
#include "common.h"
#include "logger.h"
#include "settings.h"

/* DO NOT REMOVE INITIALIZATION TO NULL, IT IS REQUIRED !! */
logger_t* g_debug_log = NULL;
logger_t* g_daily_log = NULL;
logger_t* g_monthly_log = NULL;
logger_t* g_annual_log = NULL;
logger_t* g_annual_man_log = NULL;
logger_t* g_daily_soil_log = NULL;
logger_t* g_monthly_soil_log = NULL;
logger_t* g_annual_soil_log = NULL;

extern settings_t* g_settings;

logger_t* logger_new(const char* const path, ...) {
#define BUFFER_SIZE	1024
	char buffer[BUFFER_SIZE];

	va_list va;
	logger_t* p;

	va_start(va, path);
	vsnprintf(buffer, BUFFER_SIZE, path, va);
	va_end(va);

	if ( ! path_create(buffer) ) {
		return NULL;
	}

	p = malloc(sizeof*p);

	if ( ! p ) return NULL;
	p->filename = string_copy(buffer);
	if ( ! p->filename ) {
		free(p);
		return NULL;
	}
	p->file_output = 1;
	p->std_output = 1;

	p->f = fopen(buffer, "w");
	if ( ! p->f ) {
		free(p);
		p = NULL;
	}
	return p;
#undef BUFFER_SIZE
}

void logger(logger_t *p, const char *text, ...) {
#define LOGGER_BUFFER_SIZE	4096
	char buffer[LOGGER_BUFFER_SIZE];
	va_list va;

	va_start(va, text);
	vsnprintf(buffer, LOGGER_BUFFER_SIZE, text, va);
	va_end(va);

	if ( (p && p->std_output) || (g_settings && g_settings->screen_output) ) {
		fputs(buffer, stdout);
	}

	if ( p && p->file_output && p->f ) {
		fputs(buffer, p->f);
	}
#undef LOGGER_BUFFER_SIZE
}

void logger_error(logger_t *p, const char *text, ...) {
#define LOGGER_BUFFER_SIZE	4096
	char buffer[LOGGER_BUFFER_SIZE];
	va_list va;

	va_start(va, text);
	vsnprintf(buffer, LOGGER_BUFFER_SIZE, text, va);
	va_end(va);

	fputs(buffer, stdout);

	if ( p && p->file_output && p->f ) {
		fputs(buffer, p->f);
	}
#undef LOGGER_BUFFER_SIZE
}

void error_log(const char *text, ...) {
#define LOGGER_BUFFER_SIZE	4096
	char buffer[LOGGER_BUFFER_SIZE];
	int std_out[LOG_TYPES_COUNT] = { 0 };
	va_list va;

	if ( g_debug_log && g_debug_log->std_output ) { g_debug_log->std_output = 0; std_out[DEBUG_LOG] = 1; }
	if ( g_daily_log && g_daily_log->std_output ) { g_daily_log->std_output = 0; std_out[DAILY_LOG] = 1; }
	if ( g_monthly_log && g_monthly_log->std_output ) { g_monthly_log->std_output = 0; std_out[MONTHLY_LOG] = 1; }
	if ( g_annual_log && g_annual_log->std_output ) { g_annual_log->std_output = 0; std_out[YEARLY_LOG] = 1; }
	
	if ( g_annual_man_log && g_annual_man_log->std_output ) { g_annual_man_log->std_output = 0; std_out[YEARLY_MAN_LOG] = 1; }
	
	if ( g_daily_soil_log && g_daily_soil_log->std_output ) { g_daily_soil_log->std_output = 0; std_out[SOIL_DAILY_LOG] = 1; }
	if ( g_monthly_soil_log && g_monthly_soil_log->std_output ) { g_monthly_soil_log->std_output = 0; std_out[SOIL_MONTHLY_LOG] = 1; }
	if ( g_annual_soil_log && g_annual_soil_log->std_output ) { g_annual_soil_log->std_output = 0; std_out[SOIL_YEARLY_LOG] = 1; }

	va_start(va, text);
	vsnprintf(buffer, LOGGER_BUFFER_SIZE, text, va);
	va_end(va);

	fputs(buffer, stdout);

	if ( g_debug_log && g_debug_log->f ) fputs(buffer, g_debug_log->f);
	if ( g_daily_log && g_daily_log->f ) fputs(buffer, g_daily_log->f);
	if ( g_monthly_log && g_monthly_log->f ) fputs(buffer, g_monthly_log->f);
	if ( g_annual_log && g_annual_log->f ) fputs(buffer, g_annual_log->f);

	if ( g_annual_man_log && g_annual_man_log->f ) fputs(buffer, g_annual_man_log->f);

	if ( g_daily_soil_log && g_daily_soil_log->f ) fputs(buffer, g_daily_soil_log->f);
	if ( g_monthly_soil_log && g_monthly_soil_log->f ) fputs(buffer, g_monthly_soil_log->f);
	if ( g_annual_soil_log && g_annual_soil_log->f ) fputs(buffer, g_annual_soil_log->f);

	if ( std_out[DEBUG_LOG] ) g_debug_log->std_output = 1;
	if ( std_out[DAILY_LOG] ) g_daily_log->std_output = 1;
	if ( std_out[MONTHLY_LOG] ) g_monthly_log->std_output = 1;
	if ( std_out[YEARLY_LOG] ) g_annual_log->std_output = 1;

	if ( std_out[YEARLY_MAN_LOG] ) g_annual_man_log->std_output = 1;

	if ( std_out[SOIL_DAILY_LOG] ) g_daily_soil_log->std_output = 1;
	if ( std_out[SOIL_MONTHLY_LOG] ) g_monthly_soil_log->std_output = 1;
	if ( std_out[SOIL_YEARLY_LOG] ) g_annual_soil_log->std_output = 1;

#undef LOGGER_BUFFER_SIZE
}

extern logger_t* g_debug_log;
extern logger_t* g_daily_log;
extern logger_t* g_monthly_log;
extern logger_t* g_annual_log;

extern logger_t* g_annual_man_log;

extern logger_t* g_daily_soil_log;
extern logger_t* g_monthly_soil_log;
extern logger_t* g_annual_soil_log;

void logger_close(logger_t* p) {
	if ( p ) {
		if ( p->filename ) {
			free(p->filename);
		}
		if ( p->f ) {
			fclose(p->f);
		}
		free(p);
	}
}

void logger_flush(logger_t *p) {
	if ( p ) {
		if ( p->f ) {
			fflush(p->f);
		}
	}
}
