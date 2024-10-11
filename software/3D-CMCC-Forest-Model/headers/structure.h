/* structure.h */
#ifndef STRUCTURE_H_
#define STRUCTURE_H_

int alloc_struct(void** t, int* count, int* avail, unsigned int size);

int daily_forest_structure ( cell_t *const c, const meteo_daily_t *const meteo_daily , const int  );

int monthly_forest_structure ( cell_t* const c );

int annual_forest_structure( cell_t *const, const int );

//void potential_max_min_canopy_cover (cell_t *const c);

void potential_max_min_density ( cell_t *const c );

#endif /* STRUCTURE_H_ */
