/*lai.c*/

/* includes */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lai.h"
#include "common.h"
#include "constants.h"
#include "settings.h"
#include "logger.h"

extern settings_t* g_settings;
extern logger_t* g_debug_log;

void specific_leaf_area ( const age_t *const a, species_t *const s )
{

	double n = 2;

	logger(g_debug_log, "\n**SPECIFIC LEAF AREA**\n");

	/* ager related SLA based on 3-PG model approach */
	s->value[SLA_AVG] = s->value[SLA_AVG1] + ( s->value[SLA_AVG0] - s->value[SLA_AVG1] )
									* exp( -LN2 * pow(((double)a->value / s->value[TSLA]), n));

	/* convert SLA from m2/kgDM -> m2/kgC */
	s->value[SLA_AVG] *= GC_GDM;

     //printf("s->value[SLA_AVG]   = %g\n", s->value[SLA_AVG1]  );
}

void daily_lai (cell_t *const c, age_t *const a, species_t *const s)
{
	/* NOTE: it mainly follows rationale of TREEDYN 3, Bossel, 1996, Ecological Modelling (eq. 30) */
	/* no consideration of leaf angle is taken into account at daily scale as described by Thornton (1998)
	 * it is approximates over the full day */

	/* compute age-related sla */
	specific_leaf_area      ( a, s );

	/**************************************************************************************************/

	/* total (m2) */
	s->value[LEAF_AREA]       = ( ( s->value[LEAF_C] * 1e3 / g_settings->sizeCell) * s->value[SLA_AVG] ) ;

	/* sun */
	s->value[LEAF_AREA_SUN]   = 1. - exp ( - s->value[LEAF_AREA] );

	/* shade */
	s->value[LEAF_AREA_SHADE] = s->value[LEAF_AREA] - s->value[LEAF_AREA_SUN];

	/**************************************************************************************************/

	/* compute LAI for Projected Area "Big-leaf sun/shade-model" as in DePury and Farquhar (1997)*/
	/* see Campbell and Norman, Environmental Biophysics 2nd Edition pg 259 */

	/* total */
	s->value[LAI_PROJ]       = s->value[LEAF_AREA] / s->value[CANOPY_COVER_PROJ];

	/* sun */
	s->value[LAI_SUN_PROJ]   = 1. - exp ( - s->value[LAI_PROJ] );

	/* shade */
	s->value[LAI_SHADE_PROJ] = s->value[LAI_PROJ] - s->value[LAI_SUN_PROJ];

	/* assign max annual LAI */
	s->value[MAX_LAI_PROJ]   = MAX(s->value[MAX_LAI_PROJ], s->value[LAI_PROJ]);

	/**************************************************************************************************/
	/* compute LAI for Projected Area "Big-leaf sun/shade-model" as in DePury and Farquhar (1997)*/
	/* see Campbell and Norman, Environmental Biophysics 2nd Edition pg 259 */

	/* total */
	s->value[TREE_LAI_PROJ]       = ( ( s->value[LEAF_C] * 1e3 / (double)s->counter[N_TREE] ) * s->value[SLA_AVG] ) / s->value[CROWN_AREA];

	/* sun */
	s->value[TREE_LAI_SUN_PROJ]   = 1. - exp ( - s->value[TREE_LAI_PROJ] );

	/* shade */
	s->value[TREE_LAI_SHADE_PROJ] = s->value[TREE_LAI_PROJ] - s->value[TREE_LAI_SUN_PROJ];

	/**************************************************************************************************/

	/* compute all-sided Leaf Area */
	s->value[ALL_LAI_PROJ] = s->value[LAI_PROJ] * s->value[CANOPY_COVER_PROJ];

	/* control it should anyway never happens */
	if ( s->value[LAI_PROJ] > s->value[PEAK_LAI_PROJ] )
	{
		s->value[LAI_PROJ]       =  s->value[PEAK_LAI_PROJ];

		s->value[LAI_SUN_PROJ]   = 1. - exp ( - s->value[LAI_PROJ] );

		s->value[LAI_SHADE_PROJ] = s->value[LAI_PROJ] - s->value[LAI_SUN_PROJ];
	}

	/**************************************************************************************************/

	/* compute specific leaf area projected for sun and shaded (m2/kgC) */
	if ( s->value[LEAF_C] > 0. )
	{
		/* for total leaves */
		s->value[SLA_PROJ]           =  ( s->value[LAI_PROJ] * ( s->value[CANOPY_COVER_PROJ] * g_settings->sizeCell ) ) / ( s->value[LEAF_C] * 1e3 ) ;

		/* sun leaves */
		s->value[SLA_SUN_PROJ]       = ( ( s->value[LAI_SUN_PROJ] + ( s->value[LAI_SHADE_PROJ] / s->value[SLA_RATIO] ) ) *
				( s->value[CANOPY_COVER_PROJ] * g_settings->sizeCell ) ) / ( s->value[LEAF_C] * 1e3 );

		/* shaded leaves */
		s->value[SLA_SHADE_PROJ]     = s->value[SLA_SUN_PROJ] * s->value[SLA_RATIO];
	}

	/**************************************************************************************************/

	/* computing Leaf sun and shaded Carbon and Nitrogen pools (tC-tN/ha) */
	if ( ! s->value[LEAF_C] )
	{
		s->value[LEAF_SUN_C]   = 0.;
		s->value[LEAF_SHADE_C] = 0.;

		s->value[LEAF_SUN_N]   = 0.;
		s->value[LEAF_SHADE_N] = 0.;
	}
	else
	{
		s->value[LEAF_SUN_C]   = ( ( s->value[LAI_SUN_PROJ]   * ( s->value[CANOPY_COVER_PROJ] * g_settings->sizeCell ) ) / s->value[SLA_SUN_PROJ] ) / 1e3;
		s->value[LEAF_SHADE_C] = s->value[LEAF_C] - s->value[LEAF_SUN_C];

		s->value[LEAF_SUN_N]   = s->value[LEAF_SUN_C] / s->value[CN_LEAVES];
		s->value[LEAF_SHADE_N] = s->value[LEAF_N] - s->value[LEAF_SUN_N];
	}

	/**************************************************************************************************/

	/* checks */


	CHECK_CONDITION(s->value[LAI_PROJ],       <, ZERO);
	CHECK_CONDITION(s->value[LAI_SUN_PROJ],   <, ZERO);
	CHECK_CONDITION(s->value[LAI_SHADE_PROJ], <, ZERO);
	CHECK_CONDITION(s->value[ALL_LAI_PROJ],   <, ZERO);
	CHECK_CONDITION(s->value[SLA_SUN_PROJ],   <, ZERO);
	CHECK_CONDITION(s->value[SLA_SHADE_PROJ], <, ZERO);
	CHECK_CONDITION(s->value[LAI_PROJ], >, (s->value[PEAK_LAI_PROJ] + eps));
	CHECK_CONDITION(fabs((s->value[LAI_SUN_PROJ] + s->value[LAI_SHADE_PROJ]) - s->value[LAI_PROJ]), >, eps );
	CHECK_CONDITION(fabs( ( s->value[LEAF_SUN_C] + s->value[LEAF_SHADE_C] ) - s->value[LEAF_C]) , > , eps);
	CHECK_CONDITION(fabs( ( s->value[LEAF_SUN_N] + s->value[LEAF_SHADE_N] ) - s->value[LEAF_N]) , > , eps);

	/* test: computing cell level max_lai */
	//FIXME in case of multi-class!!!!!
        // note: it should be maybe the max across class-level max value?

	c->max_lai_proj = s->value[MAX_LAI_PROJ];

	s->value[MONTHLY_LAI_PROJ] += s->value[LAI_PROJ];

}


