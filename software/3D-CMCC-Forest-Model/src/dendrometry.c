/* dendrometry.c */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dendometry.h"
#include "common.h"
#include "constants.h"
#include "logger.h"
#include "settings.h"

extern logger_t* g_debug_log;

void dendrometry_old(cell_t *const c, const int layer, const int height, const int dbh, const int age, const int species, const int year)
{
	double oldavDBH;
	double oldTreeHeight;
	double oldBasalArea;

	double max_diff = 2.;          /* maximum allowed differences in meter between initialized tree height and computed at the first year */

	height_t *h;
	dbh_t *d;
	species_t *s;

	/* this function compute at the temporal scale at which is called:
	 * -mass density
	 * -average dbh
	 * -tree height
	 * -basal area
	 * -recompute fractions of live and dead wood
	 * */

	h = &c->heights[height];
	d = &h->dbhs[dbh];
	s = &c->heights[height].dbhs[dbh].ages[age].species[species];

	/* assign previous month values */
	oldavDBH      = d->value;
	oldTreeHeight = h->value;
	oldBasalArea  = s->value[BASAL_AREA];

	logger(g_debug_log, "\n**DENDROMETRY**\n");

	/* compute tree AVDBH */
	if (s->value[STEMCONST_P] == NO_DATA && s->value[STEMPOWER_P] == NO_DATA)
	{
		/* use generic stemconst stempower values */
		if (oldavDBH < 9)
		{
			d->value = pow(((s->value[TREE_STEM_C] * 1e3) * GC_GDM / s->value[STEMCONST]), (1. / STEMPOWER_A));
		}
		else if (oldavDBH > 9 && oldavDBH < 15)
		{
			d->value = pow(((s->value[TREE_STEM_C] * 1e3) * GC_GDM / s->value[STEMCONST]), (1. / STEMPOWER_B));
		}
		else
		{
			d->value = pow(((s->value[TREE_STEM_C] * 1e3) * GC_GDM / s->value[STEMCONST]), (1. / STEMPOWER_C));
		}
	}
	else
	{
           // 5p6 recompute tree_stem_C based on the new stem_c and the new N_tree!
           // otherwise it  use the value of the day before. The differences would be anyhow minimal.
      //     printf("DENDROMETRY   s->value[STEM_C]              = %g  \n", s->value[STEM_C] );
		  

                s->value[TREE_STEM_C]                = (s->value[STEM_C]             / (double)s->counter[N_TREE]);
  //printf("DENDROMETRY   s->value[TREE_STEM_C]             = %g \n",  s->value[TREE_STEM_C]);
	//	   printf("DENDROMETRY    s->value[STEMCONST_P]             = %g  \n",  s->value[STEMCONST_P]);
	//	   printf("DENDROMETRY s->value[STEMPOWER_P]            = %g  \n",s->value[STEMPOWER_P]);
		/* use site specific stemconst stempower values */
		d->value = pow(((s->value[TREE_STEM_C] * 1e3) * GC_GDM) / s->value[STEMCONST_P], (1. / s->value[STEMPOWER_P]));
      //   printf("DENDROMETRY  d->value            = %g \n",  d->value );    
	}
	logger(g_debug_log, "-STEM_C          = %f tC/cell\n", s->value[STEM_C]);
	logger(g_debug_log, "-TREE_STEM_C     = %f tC/tree\n", s->value[TREE_STEM_C]);
	logger(g_debug_log, "-Old AVDBH       = %f cm\n", oldavDBH);
	logger(g_debug_log, "-New Average DBH = %f cm\n", d->value);

	/* check */
	CHECK_CONDITION( d->value, <, oldavDBH - eps );

	/*************************************************************************************************************************/
	/* compute Tree Height */
	/* using Chapman_Richards Function */
	/*
	 * CRA, CRB, CRC are species-specific regression coefficients
	 *-CRA coefficient is the asymptotic maximum height this coefficient represents the theoretic maximum height obtainable
	 * for a given stand or plot, this parameter could be used to represent maximum tree height.
	 * By modifying equation so that CRA changes as a function of stand age, site index, dominant height,
	 * or stand density, equation can be generalized to be more applicable to a wider range of	sites and stand ages,
	 * or more regionally applied
	 *-CRB represents exponential decay parameter
	 *-CRC represents shape parameters
	 * for references see also: M. Sharma, J. Parton, Forest Ecology and Management 249 (2007) 187–198
	 */

	logger(g_debug_log, "\n**Average Tree Height**\n");

	/* Chapman-Richards function */
	/* for references see also: R. Pilli et al. Forest Ecology and Management 237 (2006) 583–593 */
	/* note: this shouldn't be applied to saplings that are lower than 1.3 meter */
	h->value = DBH_ref + s->value[CRA] * pow (1. - exp ( - s->value[CRB] * d->value) , s->value[CRC]);

       // 5p6 Check

       if ( h->value > (s->value[CRA] + DBH_ref) )    // according to the RC equation, the max value of HA is DBH_ref + CRA (asymptotic value)
	{
		h->value = (s->value[CRA]+ DBH_ref);
	}
	logger(g_debug_log, "-Tree Height using Chapman-Richard function = %f m\n", h->value);

        /* check */
	CHECK_CONDITION (h->value, >, (s->value[CRA] + DBH_ref + eps) )

         

	logger(g_debug_log, "-Old Tree Height = %f m\n", oldTreeHeight);
	logger(g_debug_log, "-New Tree Height = %f m\n", h->value);

	/** Weibull function **/
	/* for references see also: R. Pilli et al. Forest Ecology and Management 237 (2006) 583–593 */
	/*
	 h->value = DBH_ref + s->value[WA] *(1.0 - exp ( - s->value[WB] * d->value * s->value[WC]));
	logger(g_debug_log, "-Old Tree Height = %f m\n", oldTreeHeight);
	logger(g_debug_log, "-New Tree Height = %f m\n", h->value);
	 */

	/* check */
	if ( ! s->counter[YOS] ){ CHECK_CONDITION( oldTreeHeight - h->value , > , max_diff );}
	else {CHECK_CONDITION( h->value, < , oldTreeHeight - eps );}

	/*************************************************************************************************************************/

	/* compute Basal Area and sapwood-heartwood area */

	logger(g_debug_log, "\n**Basal Area and sapwood-heartwood area**\n");

	/* basal area */
	s->value[BASAL_AREA]          = ( ( pow( ( d->value / 2.), 2. ) ) * Pi);
	s->value[BASAL_AREA_m2]       = s->value[BASAL_AREA]    * 0.0001;
	
	s->value[STAND_BASAL_AREA]    = s->value[BASAL_AREA]    * s->counter[N_TREE];
	s->value[STAND_BASAL_AREA_m2] = s->value[BASAL_AREA_m2] * s->counter[N_TREE];

	logger(g_debug_log, " BASAL AREA = %f cm^2\n", s->value[BASAL_AREA]);
	logger(g_debug_log, " BASAL BASAL_AREA_m2 = %f m^2\n", s->value[BASAL_AREA_m2]);
	logger(g_debug_log, " Stand level class basal area = %f cm2/class cell\n", s->value[STAND_BASAL_AREA]);
	logger(g_debug_log, " Stand level class basal area = %f m2/class cell\n", s->value[STAND_BASAL_AREA_m2]);

	logger(g_debug_log, "-Old Basal Area = %f cm^2\n", oldBasalArea);
	logger(g_debug_log, "-New Basal Area = %f cm^2\n", s->value[BASAL_AREA]);

	/* check */
	CHECK_CONDITION( s->value[BASAL_AREA], <, oldBasalArea - eps );

	/*************************************************************************************************************************/
#if 1
	// 5p5
	/* compute sapwood and heartwood area as fraction proportional to stem/sapwood/heartwood carbon */
	s->value[SAPWOOD_AREA]   = ( s->value[BASAL_AREA] * s->value[STEM_SAPWOOD_C] ) / s->value[STEM_C];
	s->value[HEARTWOOD_AREA] = ( s->value[BASAL_AREA] - s->value[SAPWOOD_AREA]) ;


#else
	//
	/* sapwood and heartwood are based on  Vertessy  et  al.  (1995)  and  Meinzer  et  al.(2001, 2005)*/
	s->value[SAPWOOD_AREA]     = s->value[SAP_A] * pow (d->value, s->value[SAP_B]);
	s->value[HEARTWOOD_AREA]   = s->value[BASAL_AREA] - s->value[SAPWOOD_AREA];
	s->value[SAPWOOD_PERC]     = s->value[SAPWOOD_AREA] / s->value[BASAL_AREA];

//	/* sapwood and heartwood pools */
//	s->value[STEM_SAPWOOD_C]   = s->value[STEM_C]   * s->value[SAPWOOD_PERC];
//	s->value[CROOT_SAPWOOD_C]  = s->value[CROOT_C]  * s->value[SAPWOOD_PERC];
//	s->value[BRANCH_SAPWOOD_C] = s->value[BRANCH_C] * s->value[SAPWOOD_PERC];
//	s->value[TOT_SAPWOOD_C]    = s->value[STEM_SAPWOOD_C] + s->value[CROOT_SAPWOOD_C] + s->value[BRANCH_SAPWOOD_C];
#endif

	/**/	
	// NOTE: this should not anyhow happen as the sapwood area eventually 
	// is going to decrease with age or stress events

	if ( s->value[SAPWOOD_AREA] > s->value[BASAL_AREA] + eps3 )
	{
		s->value[SAPWOOD_AREA] = s->value[BASAL_AREA];
		printf("Warning: s->value[SAPWOOD_AREA] > s->value[BASAL_AREA]\n");
	}
	/**/
	//NOTE: this should not anyhow happen 
	if ( s->value[HEARTWOOD_AREA] > s->value[BASAL_AREA] + eps3 )
	{
		s->value[HEARTWOOD_AREA] = s->value[BASAL_AREA];
		printf("Warning: s->value[HEARTWOOD_AREA] > s->value[BASAL_AREA]\n");
	}

	/* check */
	CHECK_CONDITION( s->value[SAPWOOD_AREA]   , > , s->value[BASAL_AREA] + eps3);
	CHECK_CONDITION( s->value[HEARTWOOD_AREA] , > , s->value[BASAL_AREA] + eps3);
	CHECK_CONDITION( (fabs(s->value[SAPWOOD_AREA] + s->value[HEARTWOOD_AREA]) - s->value[BASAL_AREA]), > , eps);


	logger(g_debug_log, "SAPWOOD_AREA     = %f cm^2\n",    s->value[SAPWOOD_AREA]);
	logger(g_debug_log, "HEART_WOOD_AREA  = %f cm^2\n",    s->value[HEARTWOOD_AREA]);
	logger(g_debug_log, "SAPWOOD_PERC     = %f%%\n",       s->value[SAPWOOD_PERC] * 100.);
	logger(g_debug_log, "CROOT_SAPWOOD_C  = %f tC cell\n", s->value[CROOT_SAPWOOD_C]);
	logger(g_debug_log, "BRANCH_SAPWOOD_C = %f tC cell\n", s->value[BRANCH_SAPWOOD_C]);
	logger(g_debug_log, "TOT_SAPWOOD_C    = %f tc cell\n", s->value[TOT_SAPWOOD_C]);
}


void annual_tree_increment(cell_t *const c, const int height, const int dbh, const int age, const int species, const int year)
{
	double prev_vol;
	double tree_prev_vol;

    double prev_vol2;
	double tree_prev_vol2;
	double max_diff = 0.05;           /* maximum allowed differences in m3/tree between initialized volume and computed at the first year */

	height_t *h;
	dbh_t *d;
	age_t *a;
	species_t *s;

	h = &c->heights[height];
	d = &c->heights[height].dbhs[dbh];
	a = &c->heights[height].dbhs[dbh].ages[age];
	s = &c->heights[height].dbhs[dbh].ages[age].species[species];

	/* in m^3/cell/yr */
	/* assumption: CAI = Volume t1 - Volume t0 */
	/* assumption: MAI = Volume t1 / Age t1 */

	/* CURRENT ANNUAL INCREMENT-CAI */

	logger(g_debug_log, "***VOLUME & CAI & MAI***\n");

	/* STAND VOLUME-(STEM VOLUME) */
	/* assign previous volume to temporary variables */
	prev_vol              = s->value[VOLUME];
	logger(g_debug_log, "Previous stand volume        = %f m^3/cell\n", prev_vol );

	tree_prev_vol         = s->value[TREE_VOLUME];
	logger(g_debug_log, "Previous single tree volume  = %f m^3DM/tree\n", tree_prev_vol );

	/* compute single tree volume */
	s->value[TREE_VOLUME] = (Pi * s->value[FORM_FACTOR] * pow((d->value / 100.) , 2.) * h->value) / 4.;
	logger(g_debug_log, "-Single tree volume = %g m3/tree\n", s->value[TREE_VOLUME]);

	/* compute class volume */
	s->value[VOLUME]      = s->value[TREE_VOLUME] * s->counter[N_TREE];
	logger(g_debug_log, "-Class volume = %g m3/sizeCell\n", s->value[VOLUME]);

	/* CAI-Current Annual Increment */
	s->value[CAI]         = s->value[VOLUME]      - prev_vol;
	logger(g_debug_log, "CAI-Current Annual Increment = %f m^3DM/cell/yr\n", s->value[CAI]);

	s->value[TREE_CAI]    = s->value[TREE_VOLUME] - tree_prev_vol;
	logger(g_debug_log, "CAI-Current Annual Increment = %f m^3DM/tree/yr\n", s->value[TREE_CAI]);

	/* MAI-Mean Annual Increment */
	s->value[MAI]         = s->value[VOLUME]      / (double)a->value;
	logger(g_debug_log, "MAI-Mean Annual Increment    = %f m^3DM/cell/yr \n", s->value[MAI]);

	s->value[TREE_MAI]    = s->value[TREE_VOLUME] / (double)a->value;
	logger(g_debug_log, "MAI-Mean Annual Increment    = %f m^3DM/tree/yr \n", s->value[TREE_MAI]);

	c->volume     += s->value[VOLUME];
	c->cum_volume += s->value[VOLUME];

	/* check every year after the first */
	if ( ! s->counter[YOS] ) { CHECK_CONDITION( ( tree_prev_vol - s->value[TREE_VOLUME] ), > , max_diff ); }
	else { CHECK_CONDITION( s->value[TREE_VOLUME], < , tree_prev_vol - eps ); }

    // compute data from volume2 (Vangi Elia's conversion factors)

		/* in m^3/cell/yr */
	/* assumption: CAI = Volume t1 - Volume t0 */
	/* assumption: MAI = Volume t1 / Age t1 */

	/* CURRENT ANNUAL INCREMENT-CAI */

	logger(g_debug_log, "***VOLUME2 & CAI2 & MAI2***\n");

	/* STAND VOLUME2-(STEM VOLUME2) */
	/* assign previous volume to temporary variables */
	prev_vol2              = s->value[VOLUME2];
	logger(g_debug_log, "Previous stand volume2        = %f m^3/cell\n", prev_vol2 );

	tree_prev_vol2         = s->value[TREE_VOLUME2];
	logger(g_debug_log, "Previous single tree volume2  = %f m^3DM/tree\n", tree_prev_vol2 );

	/* compute single tree volume */
	s->value[TREE_VOLUME2] = (s->value[TREE_STEM_C]*GC_GDM)/s->value[CONV_VOL_FACTOR] ; //(Pi * s->value[FORM_FACTOR] * pow((d->value / 100.) , 2.) * h->value) / 4.;
	logger(g_debug_log, "-Single tree-stem volume2 = %g m3/tree\n", s->value[TREE_VOLUME2]);

	/* compute class volume */
	s->value[VOLUME2]      = s->value[TREE_VOLUME2] * s->counter[N_TREE];
	logger(g_debug_log, "-Class volume2 = %g m3/sizeCell\n", s->value[VOLUME2]);

	/* CAI2-Current Annual Increment */
	s->value[CAI2]         = s->value[VOLUME2]      - prev_vol2;
	logger(g_debug_log, "CAI2-Current Annual Increment = %f m^3DM/cell/yr\n", s->value[CAI2]);

	s->value[TREE_CAI2]    = s->value[TREE_VOLUME2] - tree_prev_vol2;
	logger(g_debug_log, "CAI2-Current Annual Increment = %f m^3DM/tree/yr\n", s->value[TREE_CAI2]);

	/* MAI2-Mean Annual Increment */
	s->value[MAI2]         = s->value[VOLUME2]      / (double)a->value;
	logger(g_debug_log, "MAI2-Mean Annual Increment    = %f m^3DM/cell/yr \n", s->value[MAI2]);

	s->value[TREE_MAI2]    = s->value[TREE_VOLUME2] / (double)a->value;
	logger(g_debug_log, "MAI2-Mean Annual Increment    = %f m^3DM/tree/yr \n", s->value[TREE_MAI2]);

	c->volume2     += s->value[VOLUME2];
	c->cum_volume2 += s->value[VOLUME2];

	/* check every year after the first */
	if ( ! s->counter[YOS] ) { CHECK_CONDITION( ( tree_prev_vol2 - s->value[TREE_VOLUME2] ), > , max_diff ); }
	else { CHECK_CONDITION( s->value[TREE_VOLUME2], < , tree_prev_vol2 - eps ); }
}


void annual_gross_tree_increment(cell_t *const c, const int height, const int dbh, const int age, const int species, const int year)
{
	double prev_vol;
    double prev_vol2;

	double current_vol;       
	double tree_current_vol;

	double current_vol2;
	double tree_current_vol2;

	height_t *h;
	dbh_t *d;
	age_t *a;
	species_t *s;

	h = &c->heights[height];
	d = &c->heights[height].dbhs[dbh];
	a = &c->heights[height].dbhs[dbh].ages[age];
	s = &c->heights[height].dbhs[dbh].ages[age].species[species];

	/* in m^3/cell/yr */

	/* CURRENT GROSS ANNUAL INCREMENT */

	logger(g_debug_log, "*** GROSS ANNUAL INCREMENT***\n");

	/* STAND VOLUME-(STEM VOLUME) */
	/* assign previous volume to temporary variables */
	prev_vol              = s->value[VOLUME];
	logger(g_debug_log, "Previous stand volume        = %f m^3/cell\n", prev_vol );

	/* compute single tree volume */
	tree_current_vol = (Pi * s->value[FORM_FACTOR] * pow((d->value / 100.) , 2.) * h->value) / 4.;
	logger(g_debug_log, "-Single tree volume = %g m3/tree\n", tree_current_vol);

	/* compute class volume */
	current_vol      = tree_current_vol * s->counter[N_TREE];
	logger(g_debug_log, "-Class volume = %g m3/sizeCell\n", current_vol);

	/* GROSS Annual Increment */
	s->value[GROSS_INCR]         = current_vol      - prev_vol;
	logger(g_debug_log, "CAI-Current Annual Increment = %f m^3DM/cell/yr\n", s->value[GROSS_INCR]);

    // compute data from volume2 (Vangi Elia's conversion factors)

		/* in m^3/cell/yr */
	
	/* CURRENT GROSS ANNUAL INCREMENT */

	logger(g_debug_log, "*** GROSS ANNUAL INCREMENT 2 (Vangi's factors) ***\n");

	/* STAND VOLUME2-(STEM VOLUME2) */
	/* assign previous volume to temporary variables */
	prev_vol2              = s->value[VOLUME2];
	logger(g_debug_log, "Previous stand volume2        = %f m^3/cell\n", prev_vol2 );

	/* compute single tree volume */
	tree_current_vol2 = (s->value[TREE_STEM_C]*GC_GDM)/s->value[CONV_VOL_FACTOR] ; //(Pi * s->value[FORM_FACTOR] * pow((d->value / 100.) , 2.) * h->value) / 4.;
	logger(g_debug_log, "-Single tree-stem volume2 = %g m3/tree\n", tree_current_vol2 );

	/* compute class volume */
	current_vol2      = tree_current_vol2 * s->counter[N_TREE];
	logger(g_debug_log, "-Class volume2 = %g m3/sizeCell\n", current_vol2);

	/* Current GROSS Annual Increment */
	s->value[GROSS_INCR2]         = current_vol2      - prev_vol2;
	logger(g_debug_log, "CAI2-Current Annual Increment = %f m^3DM/cell/yr\n", s->value[GROSS_INCR2]);


}


void annual_minimum_reserve (species_t *const s)
{

        double actuall_sapwood;    // effective sapwood in the stand at the beginning of the year (in case that year thinning happened)

	/* recompute annual minimum reserve pool for  year allocation */
	/* these values are taken from: Schwalm and Ek, 2004 Ecological Modelling */
	/* following Krinner et al., 2005 */

	/* IMPORTANT! reserve computation if not in init data are computed from DM */
	logger(g_debug_log, "\n*annual minimum reserve*\n");


	/* compute total sapwood biomass */
	/* note: since SAP_WRES is parameterized for DryMatter it must be converted into DryMatter */
	/* Generally values range from:
	 * FOR BRANCH: 5.8 (+-0.6 s.d.) % of total DM (evergreen) and 12.8 (+-0.6 s.d.) (deciduous)
	 * FOR STEM:  1.8 (+-0.1 s.d.) % of total DM (evergreen) and 4.7 (+-0.1 s.d.) (deciduous)
	 * see: Hoch G. et al., (2003), Plant Cell and Environment (26: 1067-1081) */

       // 5p6: after the thinning is performed, the number of living trees is updated but not the actuall standing biomass (
       // the update of the biomass is done in CN_allocation). This means that in order to correctly compute the MIN_RESERVE_C
       // we need to consider the total sapwood_C after thinning as:

        //5p6   
        actuall_sapwood = s->value[TREE_SAPWOOD_C] * s->counter[N_TREE];  //N tree is updated if thinning happend 

        s->value[TOT_SAPWOOD_DM] = actuall_sapwood * GC_GDM;
	logger(g_debug_log, "--WTOT_sap_tDM = %f tDM/class \n", s->value[TOT_SAPWOOD_DM]);

	/* compute minimum annual reserve */
	s->value[MIN_RESERVE_C]  = s->value[TOT_SAPWOOD_DM] * s->value[SAP_WRES];
	logger(g_debug_log, "--MINIMUM Reserve Biomass = %f t res/class \n", s->value[MIN_RESERVE_C]);

	s->value[TREE_MIN_RESERVE_C] = s->value[MIN_RESERVE_C] / (double)s->counter[N_TREE];
	logger(g_debug_log, "--Average MINIMUM Reserve Biomass = %f tC/class tree \n", s->value[TREE_MIN_RESERVE_C]);


        // 5p6
        // set min and max reserve annual value as the previous day value (or the value set in the initialization)
        // in case it is the 
       
        s->value[MAX_ANN_RESERVE_C]          = s->value[RESERVE_C];
        s->value[MIN_ANN_RESERVE_C]          = s->value[RESERVE_C]; 

        s->value[TREE_MAX_ANN_RESERVE_C]          = s->value[RESERVE_C] / (double)s->counter[N_TREE];
        s->value[TREE_MIN_ANN_RESERVE_C]          = s->value[RESERVE_C] / (double)s->counter[N_TREE]; 

        // set min and max reserve concentration annual value at the beginning of the year 
        s->value[MAX_RESERVE_C_CONC] = ( s->value[MAX_ANN_RESERVE_C] / (s->value[TOT_SAPWOOD_C] * GC_GDM ) ) * 100.;
        s->value[MIN_RESERVE_C_CONC] = ( s->value[MIN_ANN_RESERVE_C] / (s->value[TOT_SAPWOOD_C] * GC_GDM ) ) * 100.;

}

void dendrometry ( cell_t *const c, const int layer, const int height, const int dbh, const int age, const int species, const meteo_daily_t *const meteo_daily, const int month, const int year )
{
	//	double oldavDBH;
	//	double oldTreeHeight;
	//	double oldBasalArea;
	//	double mass_density_kg;
	//
	//	double pot_par;                    /* potential absorbable incoming par */
	//	double pot_apar;                   /* potential absorbed par */
	//	double pot_apar_sun;               /* potential absorbed par sun */
	//	double pot_apar_shade;             /* potential absorbed par shade */
	//	double Light_refl_par_frac;
	//	double current_ccf;                /* crown competition factor */
	//	double current_hdf;                /* height-diameter competition factor */
	//	double current_lcf;                /* current light competition factor */
	//	double delta_stem;                 /* stem increment (kgDM/tree) */
	//	double delta_dbh;                  /* dbh increment (cm) */
	//	double delta_height;               /* height increment (m) */
	//	double dbh_m;                      /* dbh (m) */
	//	double hd_factor;                  /* HD factor based on minimum between light and crown competition */
	//
	//	double slope;
	//	double phi;
	//
	//
	//	static double annual_lcf;
	//	static int counter_annual_lcf;
	//	static double annual_ccf;
	//	static int counter_annual_ccf;
	//	static double old_stem;
	//
	//	height_t *h;
	//	dbh_t *d;
	//	age_t *a;
	//	species_t *s;
	//
	//	/* this function compute at the temporal scale at which is called
	//	 * globally:
	//	 * -mass density
	//	 * -delta dbh
	//	 * -delta tree height
	//	 * -basal area
	//	 * -sapwood and heartwood
	//	 * -recompute fractions of live and dead wood
	//	 * locally:
	//	 * -mass density in kg
	//	 * -volume
	//	 * */
	//
	//	h = &c->heights[height];
	//	d = &h->dbhs[dbh];
	//	a = &c->heights[height].dbhs[dbh].ages[age];
	//	s = &c->heights[height].dbhs[dbh].ages[age].species[species];
	//
	//	/* assign previous month values */
	//	oldavDBH = d->value;
	//	oldTreeHeight = h->value;
	//	oldBasalArea = s->value[BASAL_AREA];
	//
	//	logger(g_debug_log, "\n**DENDROMETRY**\n");
	//
	//	/*************************************************************************************************************************/
	//
	//	logger(g_debug_log, "\n**Mass density**\n");
	//
	//	/* compute mass density tDM/m3 */
	//	s->value[MASS_DENSITY] = s->value[RHO1] + (s->value[RHO0] - s->value[RHO1]) * exp(-LN2 * ((double)a->value / s->value[TRHO]));
	//	logger(g_debug_log, "-Mass Density = %f tDM/m3\n", s->value[MASS_DENSITY]);
	//
	//	/* convert to tDM-->kgDM m3 */
	//	mass_density_kg = s->value[MASS_DENSITY] * 1000.;
	//
	//	/*************************************************************************************************************************/
	//	/*
	//	 - Bossel 1996 Ecological Modelling 90, 187-227
	//	 - Peng et al., 2002, Ecological Modelling 153, 109-130
	//	 - Seidl et al., 2012, Ecological Modelling 231, 87-100
	//	 */
	//
	//	logger(g_debug_log, "*Physically based height diameter approach*\n");
	//
	//	/* compute phi (kg) */
	//	phi = (mass_density_kg * s->value[FORM_FACTOR] * Pi ) / 4.;
	//
	//	/* convert dbh cm --> m as in Bossel et al., 1996; Peng et al., 2002; Seidl et al., 2012 */
	//	dbh_m = a->value / 100.;
	//
	//	/*******************************************************************************************/
	//
	//	//	logger(g_debug_log, "Ind old stem mass = %f kgDM/tree\n", (s->value[STEM_C] * GC_GDM * 1000.)/s->counter[N_TREE]);
	//	//	logger(g_debug_log, "Ind old stem volume = %f m3/tree\n", ((s->value[STEM_C] * GC_GDM)/s->counter[N_TREE])/s->value[MASS_DENSITY]);
	//	//
	//	//	/* old stem mass from volume (cylinder) */
	//	//	old_stem_mass_from_volume = pow(dbh_m,2.) * h->value * phi;
	//	//	logger(g_debug_log, "Ind old stem mass from volume = %f kgDM/tree\n", old_stem_mass_from_volume);
	//	//
	//	//	/* old stem volume (cylinder) */
	//	//	old_stem_volume = s->value[FORM_FACTOR] * pow(d->value,2.) * h->value * 0.0001;
	//	//	logger(g_debug_log, "Ind stem volume = %f m3/tree\n", old_stem_volume);
	//
	//	/*******************************************************************************************/
	//
	//	/* compute dbh related HDMAX and HDMIN following Seidl et al., 2012 */
	//	/* compute HD_MAX */
	//	s->value[HD_MAX] = s->value[HDMAX_A] * pow(d->value, s->value[HDMAX_B]);
	//	logger(g_debug_log, "HD_MAX = %f \n", s->value[HD_MAX]);
	//	/*compute HD_MIN */
	//	s->value[HD_MIN] = s->value[HDMIN_A] * pow(d->value, s->value[HDMIN_B]);
	//	logger(g_debug_log, "HD_MIN = %f \n", s->value[HD_MIN]);
	//
	//	/*******************************************************************************************/
	//
	//	/* current height diameter factor (current_ccf) following Seidl et al., 2012 (eq.1b) */
	//	current_hdf = s->value[STEMCONST_P]*pow(dbh_m, s->value[STEMPOWER_P]);
	//	logger(g_debug_log, "height/diameter factor (Seidl method) = %f\n", current_hdf);
	//
	//	/* current height diameter factor (current_ccf) */
	//	//fixme????
	//	//current_hdf = h->value / dbh_m;
	//	//logger(g_debug_log, "height/diameter factor = %f\n", current_hdf);
	//
	//	/*******************************************************************************************/
	//
	//	/* compute current light competition factor */
	//
	//	if( s->value[LAI_PROJ] >= 1.0 )
	//	{
	//		Light_refl_par_frac = s->value[ALBEDO]/3.0;
	//	}
	//	else if ( ! s->value[LAI_PROJ] )
	//	{
	//		Light_refl_par_frac = 0.0;
	//	}
	//	else
	//	{
	//		Light_refl_par_frac = (s->value[ALBEDO]/3.0) * s->value[LAI_PROJ];
	//	}
	//
	//	if ( s->value[LAI_PROJ] > 0.)
	//	{
	//		/* compute potential absorbable incoming par less reflected */
	//		pot_par = meteo_daily->incoming_par - (meteo_daily->incoming_par * Light_refl_par_frac * s->value[DAILY_CANOPY_COVER_EXP]);
	//
	//		/* compute potential absorbed incoming par */
	//		pot_apar_sun   = pot_par * (1. - (exp(- s->value[K] * s->value[LAI_SUN_PROJ]))) * s->value[DAILY_CANOPY_COVER_EXP];
	//		pot_apar_shade = (pot_par - pot_apar_sun) * (1. - (exp(- s->value[K] * s->value[LAI_SHADE_PROJ]))) * s->value[DAILY_CANOPY_COVER_EXP];
	//		pot_apar       = pot_apar_sun + pot_apar_shade;
	//
	//		/* current light competition factor */
	//		current_lcf = s->value[APAR] / pot_apar;
	//		logger(g_debug_log, "light_competition factor = %f\n", current_lcf);
	//
	//		/*check */
	//		CHECK_CONDITION(current_lcf, >, 1);
	//		CHECK_CONDITION(current_lcf, <, ZERO );
	//
	//		/* cumulate annual light competition factor */
	//		annual_lcf += current_lcf;
	//		++counter_annual_lcf;
	//
	//		/*******************************************************************************************/
	//
	//		/* current crown competition factor (current_ccf) */
	//		slope       = s->value[DBHDCMAX]- s->value[DBHDCMIN];
	//		current_ccf = (s->value[DBHDC_EFF] - s->value[DBHDCMIN])/slope;
	//		logger(g_debug_log, "crown_competition factor = %f\n", current_ccf);
	//
	//		/*check */
	//		CHECK_CONDITION(current_ccf, >, 1);
	//		CHECK_CONDITION(current_ccf, <, ZERO );
	//
	//		/* cumulate annual crown competition factor */
	//		annual_ccf += current_ccf;
	//		++counter_annual_ccf;
	//
	//	}
	//
	//	/* annual computation */
	//	if ( c->doy == ( IS_LEAP_YEAR( c->years[year].year ) ? 366 : 365) )
	//	{
	//		logger(g_debug_log, "\n**ANNUAL DENDROMETRY**\n");
	//
	//		/* compute average annual light competition factor */
	//		annual_lcf /= (double)counter_annual_lcf;
	//
	//		/* check */
	//		CHECK_CONDITION(annual_lcf, >, 1);
	//		CHECK_CONDITION(annual_lcf, <, ZERO);
	//
	//		/* compute average annual crown competition factor */
	//		annual_ccf /= (double)counter_annual_ccf;
	//
	//		/* check */
	//		CHECK_CONDITION(annual_ccf, >, 1);
	//		CHECK_CONDITION(annual_ccf, <, ZERO);
	//
	//		/*******************************************************************************************/
	//
	//		/* partially inspired to Seidl et al., 2012 method */
	//		/* HD factor based on minimum value between light and crown competition */
	//		hd_factor = MIN (annual_ccf, annual_lcf);
	//		logger(g_debug_log, "minimum factor = %f\n", hd_factor);
	//
	//		/* reset values */
	//		annual_lcf         = 0.;
	//		counter_annual_lcf = 0;
	//		annual_ccf         = 0.;
	//		counter_annual_ccf = 0;
	//
	//		/*******************************************************************************************/
	//
	//		/* compute HD effective using a partially inspired to Seidl et al., 2012 */
	//		s->value[HD_EFF] = (s->value[HD_MAX] * (1. - hd_factor)) + (s->value[HD_MIN] * hd_factor);
	//		logger(g_debug_log, "Effective H/D ratio = %f\n", s->value[HD_EFF]);
	//
	//		/* check */
	//		CHECK_CONDITION( s->value[HD_EFF], <, s->value[HD_MIN]);
	//		CHECK_CONDITION( s->value[HD_EFF], >, s->value[HD_MAX]);
	//
	//		/*******************************************************************************************/
	//
	//		//	/* convert individual delta C stem and tDM/sizeCell --> kgDM/ tree */
	//		//	delta_stem = (s->value[C_TO_STEM] * GC_GDM * 1000. ) / s->counter[N_TREE];
	//		//	logger(g_debug_log, "delta_carbon stem = %f KgDM/month/tree \n", delta_stem);
	//		//
	//		//	/* compute individual dbh and height increment */
	//		//
	//		//	/* compute diameter increment (in m) Peng et al., 2002; Seidl et al., 2012 */
	//		//	delta_dbh_m = delta_stem / (phi * pow(dbh_m,2.) * ( ( 2. * current_hdf ) + s->value[HD_EFF]) );
	//		//	logger(g_debug_log, "delta_dbh = %f m \n", delta_dbh_m);
	//		//
	//		//	/* convert dbh from m --> cm */
	//		//	delta_dbh = delta_dbh_m * 100.;
	//		//	logger(g_debug_log, "delta_dbh = %f cm \n", delta_dbh);
	//		//
	//		//	/* compute tree height increment (in m) Peng et al., 2002; Seidl et al., 2012 */
	//		//	delta_height = delta_dbh_m * s->value[HD_EFF];
	//		//	logger(g_debug_log, "delta_height = %f m \n", delta_height);
	//
	//		/* following Peng et al., 2002*/
	//		s->value[HD_EFF] /= 100.;
	//		delta_stem = ( s->value[STEM_C] - old_stem ) * GC_GDM;
	//
	//		logger(g_debug_log, "-Annual stem increment = %f tC/year\n", delta_stem);
	//
	//		//fixme old_stem MUST be in matrix since every class have to be their own!!!
	//		old_stem = s->value[STEM_C];
	//
	//		delta_dbh = ( 4. * delta_stem ) / ( Pi * s->value[FORM_FACTOR] * s->value[MASS_DENSITY] * pow ( d->value , 2. ) *
	//				( ( 2 * ( h->value / d->value ) + s->value[HD_EFF])));
	//
	//		delta_height = delta_dbh * s->value[HD_EFF];
	//
	//		/*******************************************************************************************/
	//
	//		logger(g_debug_log, "\n**Average DBH**\n");
	//
	//		/* compute tree dbh */
	//		d->value += delta_dbh;
	//		logger(g_debug_log, "-Old AVDBH = %f cm\n", oldavDBH);
	//		logger(g_debug_log, "-New Average DBH = %f cm\n", d->value);
	//
	//		/* convert dbh cm --> m as in Bossel 1996; Peng et al., 2002; Seidl et al., 2012 */
	//		dbh_m = a->value / 100.;
	//
	//		/* check */
	//		CHECK_CONDITION( d->value, <, oldavDBH - eps );
	//
	//		/*******************************************************************************************/
	//
	//		logger(g_debug_log, "\n**Average Tree Height**\n");
	//
	//		/* compute tree height */
	//		h->value += delta_height;
	//		logger(g_debug_log, "-Old Tree Height = %f m\n", oldTreeHeight);
	//		logger(g_debug_log, "-New Tree Height = %f m\n", h->value);
	//
	//		/* check */
	//		CHECK_CONDITION( h->value, <, oldTreeHeight - eps );
	//		//fixme once change CRA with HMAX
	//		//CHECK_CONDITION( h->value, > s->value[CRA] - eps );
	//
	//		/*************************************************************************************************************************/
	//
	//		/* compute Basal Area and sapwood-heartwood area */
	//
	//		logger(g_debug_log, "\n**Basal Area and sapwood-heartwood area**\n");
	//
	//		s->value[BASAL_AREA]    = ( ( pow ( ( d->value / 2.), 2. ) ) * Pi );
	//		s->value[BASAL_AREA_m2] = s->value[BASAL_AREA] * 0.0001;
	//
	//		logger(g_debug_log, " BASAL AREA = %f cm^2\n", s->value[BASAL_AREA]);
	//		logger(g_debug_log, " BASAL BASAL_AREA_m2 = %f m^2\n", s->value[BASAL_AREA_m2]);
	//#if 1
	//		//test_new
	//		s->value[SAPWOOD_AREA]     = (s->value[SAP_A] * pow (d->value, s->value[SAP_B]));
	//		s->value[HEARTWOOD_AREA]   = (s->value[BASAL_AREA] - s->value[SAPWOOD_AREA]);
	//		s->value[SAPWOOD_PERC]     = (s->value[SAPWOOD_AREA]) / s->value[BASAL_AREA];
	//		s->value[STEM_SAPWOOD_C]   = (s->value[STEM_C]   * s->value[SAPWOOD_PERC]);
	//		s->value[CROOT_SAPWOOD_C]  = (s->value[CROOT_C]  * s->value[SAPWOOD_PERC]);
	//		s->value[BRANCH_SAPWOOD_C] = (s->value[BRANCH_C] * s->value[SAPWOOD_PERC]);
	//
	//		logger(g_debug_log, " SAPWOOD_AREA                    = %f cm2/tree\n", s->value[SAPWOOD_AREA]);
	//		logger(g_debug_log, " HEART_WOOD_AREA                 = %f cm2/tree\n", s->value[HEARTWOOD_AREA]);
	//		logger(g_debug_log, " sapwood perc                    = %f%%\n", s->value[SAPWOOD_PERC]*100);
	//		logger(g_debug_log, " Sapwood stem biomass            = %f tC class cell \n", s->value[STEM_SAPWOOD_C]);
	//		logger(g_debug_log, " Sapwood coarse root biomass     = %f tC class cell \n", s->value[CROOT_SAPWOOD_C]);
	//		logger(g_debug_log, " Sapwood branch and bark biomass = %f tC class cell \n", s->value[BRANCH_SAPWOOD_C]);
	//#endif
	//		s->value[TOT_SAPWOOD_C]       = s->value[STEM_SAPWOOD_C] + s->value[CROOT_SAPWOOD_C] + s->value[BRANCH_SAPWOOD_C];
	//		s->value[STAND_BASAL_AREA]    = s->value[BASAL_AREA] * s->counter[N_TREE];
	//		s->value[STAND_BASAL_AREA_m2] = s->value[BASAL_AREA_m2] * s->counter[N_TREE];
	//
	//		logger(g_debug_log, " Total Sapwood biomass = %f tc class cell \n", s->value[TOT_SAPWOOD_C]);
	//		logger(g_debug_log, " Stand level class basal area = %f cm^2/class cell\n", s->value[STAND_BASAL_AREA]);
	//		logger(g_debug_log, " Stand level class basal area (meters) = %f m^2/class cell\n", s->value[STAND_BASAL_AREA_m2]);
	//		logger(g_debug_log, "-Old Basal Area = %f cm^2\n", oldBasalArea);
	//		logger(g_debug_log, "-New Basal Area = %f cm^2\n", s->value[BASAL_AREA]);
	//
	//		/* check */
	//		CHECK_CONDITION( s->value[BASAL_AREA], <, oldBasalArea - eps );
	//	}
}
