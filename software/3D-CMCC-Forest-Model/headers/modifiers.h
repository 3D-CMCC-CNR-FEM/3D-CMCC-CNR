/* modifiers.h */
#ifndef MODIFIERS_H_
#define MODIFIERS_H_

#include "matrix.h"

void modifiers(cell_t *const c, const int layer, const int height, const int dbh, const int age, const int species, const meteo_daily_t *const meteo_daily,	const meteo_annual_t *const meteo_annual ,  const int month , const int day);
//5p606 added month and doy

#endif /* MODIFIERS_H_ */
