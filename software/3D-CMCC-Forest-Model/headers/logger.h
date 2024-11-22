/* logger.h */
#ifndef LOGGER_H_
#define LOGGER_H_

#include <stdio.h> // for FILE

/* do not change order */
enum {
	DEBUG_LOG
	, DAILY_LOG
	, MONTHLY_LOG
	, YEARLY_LOG
	, YEARLY_MAN_LOG
	, SOIL_DAILY_LOG
	, SOIL_MONTHLY_LOG
	, SOIL_YEARLY_LOG

	, LOG_TYPES_COUNT
};

typedef struct {
	FILE *f;
	char *filename;
	int file_output;
	int std_output;
} logger_t;

logger_t* logger_new(const char* const text, ...);
#define logger_enable_file(p)	{if((p))(p)->file_output=1;}
#define logger_disable_file(p)	{if((p))(p)->file_output=0;}
#define logger_enable_std(p)	{if((p))(p)->std_output=1;}
#define logger_disable_std(p)	{if((p))(p)->std_output=0;}
void logger(logger_t *p, const char *text, ...);
void logger_error(logger_t *p, const char *text, ...);
void error_log(const char *text, ...);
void logger_flush(logger_t *p);
void logger_close(logger_t* p);


#endif /* LOGGER_H_ */
