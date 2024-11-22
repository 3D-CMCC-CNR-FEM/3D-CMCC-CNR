/* common.h */
#ifndef COMMON_H
#define COMMON_H

#define PROGRAM_NAME		"3D (Three Dimensional) CMCC (Coupled Model Carbon Cycle) FEM (Forest Ecosystem Module)"
#define PROGRAM_VERSION		"5.606"
#define PROGRAM_FULL_NAME	PROGRAM_NAME" v."PROGRAM_VERSION""

#ifdef _WIN32
#define FOLDER_DELIMITER	"\\"
#define FOLDER_DELIMITER_C	'\\'
#else
#define FOLDER_DELIMITER	"/"
#define FOLDER_DELIMITER_C	'/'
#endif

#define INVALID_VALUE			-9999
#define IS_INVALID_VALUE(x)		((INVALID_VALUE==(x)))
#define IS_LEAP_YEAR(x)			(((x) % 4 == 0 && (x) % 100 != 0) || (x) % 400 == 0)
#define SIZE_OF_ARRAY(a)		sizeof((a))/sizeof((a)[0])
#define ARE_FLOATS_EQUAL(a,b)	(((a)==(b)))

#ifndef MAX
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a, b)  (((a) < (b)) ? (a) : (b))
#endif

#ifndef ROUND
/* we add 0.5 so if x is > 0.5 we truncate to next integer */
#define ROUND(x)	((x)>=0?(long)((x)+0.5):(long)((x)-0.5))
#else
#define ROUND		round
#endif

#define XSTR(a)		STR(a)
#define STR(a)		#a

#include "logger.h"

// ALESSIOR:
// cannot include header 'cause "print_output.h" include "matrix.h"
// that define CAI that is already defined by windows SDK
//#include "print_output.h"
void print_model_paths(logger_t *const _log);
void print_model_settings(logger_t*const log);

extern logger_t* g_daily_log;
extern logger_t* g_monthly_log;
extern logger_t* g_annual_log;
extern logger_t* g_annual_man_log;
//extern logger_t* g_soil_log;

#define CHECK_CONDITION(x,c,y) {					\
		if ( is_nan((x)) ) {						\
			puts(XSTR(x)" is NAN in " XSTR(__FILE__) " in " XSTR(__LINE__));					\
			exit(1);								\
		} else if ( is_inf((x)) ) {					\
			puts(XSTR(x)" is INF in " XSTR(__FILE__) " in " XSTR(__LINE__));					\
			exit(1);								\
		} else if ( (x)c(y) ) { 					\
			char buf[256]; 							\
			sprintf(buf, "\n#error:\n#condition \"%s %s %s\" is true\n#value of %s is %g\n#value of %s is %g\n#in %s on line %d\n"	\
					, XSTR(x)				\
					, XSTR(c)				\
					, XSTR(y)				\
					, XSTR(x)				\
					, (double)(x)			\
					, XSTR(y)				\
					, (double)(y)			\
					, __FILE__				\
					, __LINE__				\
			);										\
			logger_error(NULL, buf);				\
			if ( g_daily_log ) {					\
				g_daily_log->std_output = 0;		\
				g_daily_log->file_output = 1;		\
				logger(g_daily_log, buf);			\
				print_model_paths(g_daily_log);		\
				print_model_settings(g_daily_log);	\
			}										\
			if ( g_monthly_log ) {					\
				g_monthly_log->std_output = 0;		\
				g_monthly_log->file_output = 1;		\
				logger(g_monthly_log, buf);			\
				print_model_paths(g_monthly_log);	\
				print_model_settings(g_monthly_log);\
			}										\
			if ( g_annual_log ) {					\
				g_annual_log->std_output = 0;		\
				g_annual_log->file_output = 1;		\
				logger(g_annual_log, buf);			\
				print_model_paths(g_annual_log);	\
				print_model_settings(g_annual_log);	\
			}										\
			if ( g_annual_man_log ) {					\
				g_annual_man_log->std_output = 0;		\
				g_annual_man_log->file_output = 1;		\
				logger(g_annual_man_log, buf);			\
				print_model_paths(g_annual_man_log);	\
				print_model_settings(g_annual_man_log);	\
			}	                                    \
			exit(1);								\
		}											\
}




#define convert_string_to_int(s, err) ((int)convert_string_to_float((s),(err)))

double convert_string_to_float(const char *const string, int *const error);
int string_compare_i(const char *str1, const char *str2);
char* string_copy(const char *const string);
char* string_tokenizer(char *string, const char *delimiters, char **p);
char* strstr_i(char* str1, const char* str2);
char* get_current_path(void);
double QuadM (const double a, const double b, const double c, int* const err);
double QuadP (const double a, const double b, const double c, int* const err);

int file_get_rows_count(const char* const filename);
int file_load_in_memory(const char* const filename, char** result);

int file_copy(const char* const filename, const char* const path);
const char* file_get_name_only(const char* const filename);
int path_create(const char *const path);
char* concatenate_path(char* s1, char* s2);
int file_exists(const char *const file);
int has_path_delimiter(const char* const s);
const char* datetime_current(void);
void timer_init(void);
double timer_get(void);
int istab(const int c);
int is_nan(double x);
int is_inf(double x);

#endif /* COMMON_H */
