/*peak_lai.c*/

/* includes */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "common.h"
#include "peak_lai.h"
#include "constants.h"
#include "settings.h"
#include "logger.h"
#include "lai.h"

extern settings_t* g_settings;

void peak_lai( age_t *const a, species_t *const s, const int day, const int month, const int years)
{
	/* compute age-related sla */
	specific_leaf_area      ( a, s );

	/* compute annual Peak Projected Lai (m2/m2) (tree level) */
	s->value[PEAK_LAI_PROJ] = ( ( s->value[SAPWOOD_AREA] / 10000. ) * s->value[SAP_LEAF]) / s->value[CROWN_AREA_PROJ];
   //   printf(" IN PEAK LAI 5 s->value[PEAK_LAI_PROJ] %f!!!\n",s->value[PEAK_LAI_PROJ]);
	/* check if Peak Lai exceeds Maximum prescribed Peak Lai (this shouldn't happens) */
	if (s->value[PEAK_LAI_PROJ] > MAX_PEAK_LAI_PROJ) s->value[PEAK_LAI_PROJ] = MAX_PEAK_LAI_PROJ;

	/* compute max leaf carbon (tC/cell) at Peak Projected Lai */
	s->value[MAX_LEAF_C]   = ( ( s->value[PEAK_LAI_PROJ] / s->value[SLA_AVG] ) / 1e3 ) * ( s->value[CANOPY_COVER_PROJ] * g_settings->sizeCell );

	/* then compute fine root carbon (tC/cell) at Peak Projected Lai */
	s->value[MAX_FROOT_C]  = s->value[MAX_LEAF_C] * s->value[FINE_ROOT_LEAF];

	/***************************************************************************************************************/

	/* compute annual Peak Exposed Lai (m2/m2) */
	s->value[PEAK_LAI_EXP] = ( ( s->value[SAPWOOD_AREA] / 10000. ) * s->value[SAP_LEAF]) / s->value[CROWN_AREA_PROJ];

	/***************************************************************************************************************/
	/* note: special case for evergreen */

	if ((s->value[PHENOLOGY] == 1.1 || s->value[PHENOLOGY] == 1.2) && ( !day && !month && !years ))
	{
		if ((s->value[LAI_PROJ] > s->value[PEAK_LAI_PROJ]) || (s->value[LEAF_C] > s->value[MAX_LEAF_C]))
		{
			/* force current LAI to Peak Lai */
			s->value[LAI_PROJ] = s->value[PEAK_LAI_PROJ];

			/* then recompute leaf biomass */
			s->value[LEAF_C]   = s->value[MAX_LEAF_C];

			/* then recompute fine root biomass */
			s->value[FROOT_C]  = s->value[MAX_FROOT_C];
		}
	}
}
