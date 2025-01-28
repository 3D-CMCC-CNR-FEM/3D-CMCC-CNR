/* print_output.h */
#ifndef PRINT_OUTPUT_H_
#define PRINT_OUTPUT_H_

#include "matrix.h"
#include "logger.h"


void EOD_print_output_cell_level(cell_t *const, const int, const int, const int, const int);
void EOD_cell_msg(void);

void EOM_print_output_cell_level(cell_t *const, const int, const int, const int);
void EOM_cell_msg(void);

void EOY_print_output_cell_level(cell_t *const, const int, const int);
void EOY_cell_msg(void);

void EOY_print_output_class_level_mortality(cell_t *const c, const int height, const int dbh , const int age , const int species, const int year );
//void EOY_cell_msg(void);


void EOD_print_output_soil_cell_level(cell_t *const, const int, const int, const int, const int);
void EOD_soil_msg(void);

void EOM_print_output_soil_cell_level(cell_t *const, const int, const int, const int);
void EOM_soil_msg(void);

void EOY_print_output_soil_cell_level(cell_t *const, const int, const int);
void EOY_soil_msg(void);

void print_model_paths(logger_t *const _log);

void print_model_settings(logger_t*const log);

#endif /* PRINT_OUTPUT_H_ */
