/* met_data.h */
#ifndef MET_DATA_H_
#define MET_DATA_H_

#include "matrix.h"

typedef enum {
	WEIGHTED_MEAN_TAVG
	, WEIGHTED_MEAN_TDAY
	, WEIGHTED_MEAN_TNIGHT
	, WEIGHTED_MEAN_TSOIL

	, WEIGHTED_MEAN_COUNT
} e_weighted_average_var;
typedef enum {
	AVERAGED_TAVG
	, AVERAGED_TDAY
	, AVERAGED_TNIGHT
	, AVERAGED_TSOIL

	, AVERAGED_COUNT
} e_averaged_var;

void Radiation (cell_t *const c, const int day, const int month, const int year);
void Daily_avg_temperature(meteo_t *met, const int day, const int month);
void Daylight_avg_temperature(meteo_t *const met, const int day, const int month);
void Nightime_avg_temperature(meteo_t *const met, const int day, const int month);
void Thermic_sum (cell_t *c, meteo_t *met, const int day, const int month, const int year);
void Air_pressure(meteo_t *met, const int day, const int month);
void Air_density (meteo_t *met, const int day, const int month);
void Day_Length(cell_t *c, const int day, const int month, const int year);
void Latent_heat(meteo_t *met, const int day, const int month);
void Psychrometric(meteo_t *met, const int day, const int month);
void Sat_vapour_pressure(cell_t *const c, const int day, const int month, const int year);
void Dew_temperature(meteo_t *const met, const int day, const int month);
void Check_prcp(cell_t *c, const int day, const int month, const int year);
void Averaged_temperature(const cell_t *const c, const e_averaged_var var, int day, int month, int year);
void Weighted_average_temperature(const cell_t *const c, const e_weighted_average_var var, int day, int month, int year);
void Soil_temperature(const cell_t *const c, int day, int month, int year);
void Daily_Ndeposition(const cell_t *const c, int day, int month, int year_iyearndex);
// 5.7
//void Seedling_soil_par (const cell_t *const c, const meteo_daily_t *const meteo_daily, int day, int month, int year);
void Seedling_soil_par (const cell_t *const c, const meteo_daily_t *const meteo_daily, int month, int year);


#endif /* MET_DATA_H_ */
