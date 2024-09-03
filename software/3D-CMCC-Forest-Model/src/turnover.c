/*turnover.c*/

/* includes */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "turnover.h"
#include "common.h"
#include "constants.h"
#include "logger.h"

void sapwood_turnover( cell_t *const c, age_t *const a, species_t *const s, const int day, const int month, const int year )
{
	double effective_sapwood_turnover;       /* effective sapwood turnover */

	/*****************************************************************************************************************/
	/*****************************************************************************************************************/

#if 1

	double juv_sapwood_turnover_frac;       /* juvenile fraction of turnover */
	double mat_sapwood_turnover_frac;       /* mature fraction of turnover */
	int t_age;                              /* age at which TURNOVER = (juvenile_sapwood_turnover_frac/mature_sapwood_turnover_frac)/2 */
	double n = 2.;


	/** age related function for sapwood turnover **/

	/* "Similar analyses of Norway spruce, Scots pine, Sessile oak, and European beech
		show a species-specific decrease in the proportion of sapwood from f sw = 1.0 in
		the juvenile phase to 0.25–0.50 in the mature phase. Forest Dynamics, Growth and Yield book Hans Pretzsch (pg. 99) "

	/* assumption: this function assumes that turnover rate increases with increasing tree age */

        juv_sapwood_turnover_frac = 0.016;   // 60 year of MRT (ref?LPJ Poulter?) to mimic that very young saplings have high amount of sapwood
                                             // FIXME it should be moved to species.txt
	mat_sapwood_turnover_frac = s->value[SAPWOOD_TURNOVER]; 
	t_age                     = 30;                         //fixme it should be moved to species.txt

	/* note: it recomputes variables in species.txt */
	effective_sapwood_turnover = mat_sapwood_turnover_frac + ( juv_sapwood_turnover_frac - mat_sapwood_turnover_frac ) *
			exp( -LN2 * pow( ( (double)a->value / t_age ), n ) );
       
#else

	/* note: assumes fixed value in species.txt */
	effective_sapwood_turnover = s->value[SAPWOOD_TURNOVER];

#endif
      
	s->value[DAILY_SAPWOOD_TURNOVER] = effective_sapwood_turnover / ( IS_LEAP_YEAR ( c->years[year].year ) ? 366 : 365 );

	/*****************************************************************************************************************/
	/************************************************ sapwood turnover ***********************************************/

	if ( ! day &&  ! month )  // first january
	{

     
		//		if ( s->counter[YOS] == 1 )
		//		{
		//			s->value[YEARLY_C_TO_STEM_SAPWOOD]   = s->value[STEM_SAPWOOD_C];
		//			s->value[YEARLY_C_TO_CROOT_SAPWOOD]  = s->value[CROOT_SAPWOOD_C];
		//			s->value[YEARLY_C_TO_BRANCH_SAPWOOD] = s->value[BRANCH_SAPWOOD_C];
		//		}

#if 0
             
            

		/* as in: LPJ, Sitch et al., 2003 and  Orchidee, Krinner et al., 2005; all sapwood tunovers each year */

		/* stem sapwood */
		s->value[C_STEM_SAPWOOD_TO_HEARTWOOD]   = s->value[YEARLY_C_TO_STEM_SAPWOOD]   * s->value[DAILY_SAPWOOD_TURNOVER];

		/* coarse root sapwood */
		s->value[C_CROOT_SAPWOOD_TO_HEARTWOOD]  = s->value[YEARLY_C_TO_CROOT_SAPWOOD]  * s->value[DAILY_SAPWOOD_TURNOVER];

		/* branch sapwood */
		s->value[C_BRANCH_SAPWOOD_TO_HEARTWOOD] = s->value[YEARLY_C_TO_BRANCH_SAPWOOD] * s->value[DAILY_SAPWOOD_TURNOVER];

                /* reset annual sapwood values once used */
		s->value[YEARLY_C_TO_STEM_SAPWOOD]   = 0.;
		s->value[YEARLY_C_TO_CROOT_SAPWOOD]  = 0.;
		s->value[YEARLY_C_TO_BRANCH_SAPWOOD] = 0.;

#else

               // 5p6 according to GCB 2019 paper 
                /* stem sapwood */
		s->value[C_STEM_SAPWOOD_TO_HEARTWOOD]   = s->value[STEM_SAPWOOD_C]   * s->value[DAILY_SAPWOOD_TURNOVER];

		/* coarse root sapwood */
		s->value[C_CROOT_SAPWOOD_TO_HEARTWOOD]  = s->value[CROOT_SAPWOOD_C]  * s->value[DAILY_SAPWOOD_TURNOVER];

		/* branch sapwood */
		s->value[C_BRANCH_SAPWOOD_TO_HEARTWOOD] = s->value[BRANCH_SAPWOOD_C] * s->value[DAILY_SAPWOOD_TURNOVER];


		/* reset annual sapwood values once used */
		s->value[YEARLY_C_TO_STEM_SAPWOOD]   = 0.;
		s->value[YEARLY_C_TO_CROOT_SAPWOOD]  = 0.;
		s->value[YEARLY_C_TO_BRANCH_SAPWOOD] = 0.;

#endif

	}

       

	/* daily stem turnover sapwood to heartwood*/

	s->value[STEM_SAPWOOD_C]     -= s->value[C_STEM_SAPWOOD_TO_HEARTWOOD];
	s->value[STEM_HEARTWOOD_C]   += s->value[C_STEM_SAPWOOD_TO_HEARTWOOD];

	/* daily coarse root turnover sapwood to heartwood */

	s->value[CROOT_SAPWOOD_C]    -= s->value[C_CROOT_SAPWOOD_TO_HEARTWOOD];
	s->value[CROOT_HEARTWOOD_C]  += s->value[C_CROOT_SAPWOOD_TO_HEARTWOOD];

	/* daily branch turnover sapwood to heartwood */

	s->value[BRANCH_SAPWOOD_C]   -= s->value[C_BRANCH_SAPWOOD_TO_HEARTWOOD];
	s->value[BRANCH_HEARTWOOD_C] += s->value[C_BRANCH_SAPWOOD_TO_HEARTWOOD];

	/* Overall sapwood and heartwood stocks */

	s->value[TOT_SAPWOOD_C]       = s->value[STEM_SAPWOOD_C]   + s->value[CROOT_SAPWOOD_C]   + s->value[BRANCH_SAPWOOD_C];
	s->value[TOT_HEARTWOOD_C]     = s->value[STEM_HEARTWOOD_C] + s->value[CROOT_HEARTWOOD_C] + s->value[BRANCH_HEARTWOOD_C];

      // update tree level pool

	s->value[TREE_STEM_SAPWOOD_C]        = (s->value[STEM_SAPWOOD_C]     / (double)s->counter[N_TREE]);
	s->value[TREE_STEM_HEARTWOOD_C]      = (s->value[STEM_HEARTWOOD_C]   / (double)s->counter[N_TREE]);
	s->value[TREE_CROOT_SAPWOOD_C]       = (s->value[CROOT_SAPWOOD_C]    / (double)s->counter[N_TREE]);
	s->value[TREE_CROOT_HEARTWOOD_C]     = (s->value[CROOT_HEARTWOOD_C]  / (double)s->counter[N_TREE]);
	s->value[TREE_BRANCH_SAPWOOD_C]      = (s->value[BRANCH_SAPWOOD_C]   / (double)s->counter[N_TREE]);
	s->value[TREE_BRANCH_HEARTWOOD_C]    = (s->value[BRANCH_HEARTWOOD_C] / (double)s->counter[N_TREE]);
	s->value[TREE_SAPWOOD_C]             = (s->value[TOT_SAPWOOD_C]      / (double)s->counter[N_TREE]);
	s->value[TREE_HEARTWOOD_C]           = (s->value[TOT_HEARTWOOD_C]    / (double)s->counter[N_TREE]);


	/******************************************************************************************************************/
	/******************************************************************************************************************/

}

void livewood_turnover( cell_t *const c, age_t *const a, species_t *const s, const int day, const int month, const int year )
{
	double effective_livewood_turnover;      /* effective livewood turnover */

	/*****************************************************************************************************************/
	/*****************************************************************************************************************/

#if 1

	double juv_livewood_turnover_frac;       /* juvenile fraction of turnover */
	double mat_livewood_turnover_frac;       /* mature fraction of turnover */
	int t_age;                               /* age at which TURNOVER = (juvenile_livewood_turnover_frac/mature_livewood_turnover_frac)/2 */
	double n = 2.;

	/** age related function for livewood turnover **/

	/* assumption: this function assumes that turnover rate increases with increasing tree age */

        
        juv_livewood_turnover_frac = 0.016;  	//fixme it should be moved to species.txt

	mat_livewood_turnover_frac = s->value[LIVEWOOD_TURNOVER]; //fixme it should be moved to species.txt
	t_age                      = 30;                          //fixme it should be moved to species.txt

	/* note: it recomputes variables in species.txt */
	effective_livewood_turnover = mat_livewood_turnover_frac + ( juv_livewood_turnover_frac - mat_livewood_turnover_frac ) *
			exp( -LN2 * pow( ( (double)a->value / t_age ), n ) );

#else

	/* note: assumes fixed value in species.txt */
	effective_livewood_turnover = s->value[LIVEWOOD_TURNOVER];

#endif
     
	/* test it accounts for fractions for previous year carbon biomass for daily turnover */
	/* following Krinner et al., 2005 turnover occurs every day of the year */

	/* daily fraction of biomass to remove */
	s->value[DAILY_LIVEWOOD_TURNOVER] = effective_livewood_turnover / ( IS_LEAP_YEAR ( c->years[year].year ) ? 366 : 365 );

	/*****************************************************************************************************************/
	/************************************************ livewood turnover ***********************************************/

	if ( ! day &&  ! month ) // first january
	{
		//		if ( s->counter[YOS] == 1 )
		//		{
		//			s->value[YEARLY_C_TO_STEM_SAPWOOD]   = s->value[STEM_LIVEWOOD_C];
		//			s->value[YEARLY_C_TO_CROOT_SAPWOOD]  = s->value[CROOT_LIVEWOOD_C];
		//			s->value[YEARLY_C_TO_BRANCH_SAPWOOD] = s->value[BRANCH_LIVEWOOD_C];
		//		}

#if 0
           

		/* stem livewood */
		s->value[C_STEM_LIVEWOOD_TO_DEADWOOD]   = s->value[YEARLY_C_TO_STEM_LIVEWOOD]   * s->value[DAILY_LIVEWOOD_TURNOVER];

		/* coarse root livewood */
		s->value[C_CROOT_LIVEWOOD_TO_DEADWOOD]  = s->value[YEARLY_C_TO_CROOT_LIVEWOOD]  * s->value[DAILY_LIVEWOOD_TURNOVER];

		/* branch livewood */
		s->value[C_BRANCH_LIVEWOOD_TO_DEADWOOD] = s->value[YEARLY_C_TO_BRANCH_LIVEWOOD] * s->value[DAILY_LIVEWOOD_TURNOVER];

             //  printf("Computing livewood to deadwood \n");
             // printf("C_STEM_LIVEWOOD_TO_DEADWOOD %f\n",s->value[C_STEM_LIVEWOOD_TO_DEADWOOD] );


                /* reset annual  livewood values once used */
		s->value[YEARLY_C_TO_STEM_LIVEWOOD]   = 0.;
		s->value[YEARLY_C_TO_CROOT_LIVEWOOD]  = 0.;
		s->value[YEARLY_C_TO_BRANCH_LIVEWOOD] = 0.;
     
#else
               // 5p6 according to GCB 2019 paper 

                s->value[C_STEM_LIVEWOOD_TO_DEADWOOD]   = s->value[STEM_LIVEWOOD_C]   * s->value[DAILY_LIVEWOOD_TURNOVER];

		/* coarse root livewood */
		s->value[C_CROOT_LIVEWOOD_TO_DEADWOOD]  = s->value[CROOT_LIVEWOOD_C] * s->value[DAILY_LIVEWOOD_TURNOVER];

		/* branch livewood */
		s->value[C_BRANCH_LIVEWOOD_TO_DEADWOOD] = s->value[BRANCH_LIVEWOOD_C] * s->value[DAILY_LIVEWOOD_TURNOVER];

		/* reset annual  livewood values once used */
		s->value[YEARLY_C_TO_STEM_LIVEWOOD]   = 0.;
		s->value[YEARLY_C_TO_CROOT_LIVEWOOD]  = 0.;
		s->value[YEARLY_C_TO_BRANCH_LIVEWOOD] = 0.;

#endif

	}

	/* daily stem turnover livewood to heartwood */

	s->value[STEM_LIVEWOOD_C]   -= s->value[C_STEM_LIVEWOOD_TO_DEADWOOD];
	s->value[STEM_DEADWOOD_C]   += s->value[C_STEM_LIVEWOOD_TO_DEADWOOD];

	/* daily coarse root turnover livewood to heartwood */

	s->value[CROOT_LIVEWOOD_C]  -= s->value[C_CROOT_LIVEWOOD_TO_DEADWOOD];
	s->value[CROOT_DEADWOOD_C]  += s->value[C_CROOT_LIVEWOOD_TO_DEADWOOD];

	/* daily branch turnover livewood to heartwood */

	s->value[BRANCH_LIVEWOOD_C] -= s->value[C_BRANCH_LIVEWOOD_TO_DEADWOOD];
	s->value[BRANCH_DEADWOOD_C] += s->value[C_BRANCH_LIVEWOOD_TO_DEADWOOD];


	/* Overall livewood and heartwood stocks */

	s->value[TOT_LIVEWOOD_C] = s->value[STEM_LIVEWOOD_C] + s->value[CROOT_LIVEWOOD_C] + s->value[BRANCH_LIVEWOOD_C];
	s->value[TOT_DEADWOOD_C] = s->value[STEM_DEADWOOD_C] + s->value[CROOT_DEADWOOD_C] + s->value[BRANCH_DEADWOOD_C];

    // update tree level pool

	s->value[TREE_STEM_LIVEWOOD_C]       = (s->value[STEM_LIVEWOOD_C]    / (double)s->counter[N_TREE]);
	s->value[TREE_STEM_DEADWOOD_C]       = (s->value[STEM_DEADWOOD_C]    / (double)s->counter[N_TREE]);
	s->value[TREE_CROOT_LIVEWOOD_C]      = (s->value[CROOT_LIVEWOOD_C]   / (double)s->counter[N_TREE]);
	s->value[TREE_CROOT_DEADWOOD_C]      = (s->value[CROOT_DEADWOOD_C]   / (double)s->counter[N_TREE]);
	s->value[TREE_BRANCH_LIVEWOOD_C]     = (s->value[BRANCH_LIVEWOOD_C]  / (double)s->counter[N_TREE]);
	s->value[TREE_BRANCH_DEADWOOD_C]     = (s->value[BRANCH_DEADWOOD_C]  / (double)s->counter[N_TREE]);
	s->value[TREE_TOT_LIVEWOOD_C]        = (s->value[TOT_LIVEWOOD_C]     / (double)s->counter[N_TREE]);
	s->value[TREE_TOT_DEADWOOD_C]        = (s->value[TOT_DEADWOOD_C]     / (double)s->counter[N_TREE]);



	/******************************************************************************************************************/
	/******************************************************************************************************************/
}


//this is the old version used in the 5.4 version (and lower)
void turnover( cell_t *const c, age_t *const a, species_t *const s, const int day, const int month, const int year )
{
	double effective_livewood_turnover;      /* effective livewood turnover */

	/*****************************************************************************************************************/
	/*****************************************************************************************************************/

#if 1

	double juv_livewood_turnover_frac;       /* juvenile fraction of turnover */
	double mat_livewood_turnover_frac;       /* mature fraction of turnover */
	int t_age;                               /* age at which TURNOVER = (juvenile_livewood_turnover_frac/mature_livewood_turnover_frac)/2 */
	double n = 2.;

	/** age related function for turnover **/

	/* assumption: this function assumes that turnover rate increases with increasing tree age */
        
        juv_livewood_turnover_frac = 0.1;                         //fixme it should be moved to species.txt
	mat_livewood_turnover_frac = s->value[LIVEWOOD_TURNOVER]; //fixme it should be moved to species.txt
	t_age                      = 30;                          //fixme it should be moved to species.txt

	/* note: it recomputes variables in species.txt */
	effective_livewood_turnover = mat_livewood_turnover_frac + ( juv_livewood_turnover_frac - mat_livewood_turnover_frac ) *
			exp( -LN2 * pow( ( (double)a->value / t_age ), n ) );

#else

	/* note: assumes fixed value in species.txt */
	effective_livewood_turnover = s->value[LIVEWOOD_TURNOVER];

#endif

	/**************************************************************************************************************************/

	/* test it accounts for fractions for previous year carbon biomass for daily turnover */
	/* following Krinner et al., 2005 turnover occurs every day of the year */
	/* daily fraction of biomass to remove */

	s->value[DAILY_LIVEWOOD_TURNOVER] = effective_livewood_turnover / ( IS_LEAP_YEAR ( c->years[year].year ) ? 366 : 365 );


	if ( ! day &&  ! month && s->counter[YOS] > 1 )
	{
		/* stem */
		s->value[C_STEM_LIVEWOOD_TO_DEADWOOD]   = s->value[YEARLY_C_TO_STEM]   * s->value[DAILY_LIVEWOOD_TURNOVER];

		/* coarse root */
		s->value[C_CROOT_LIVEWOOD_TO_DEADWOOD]  = s->value[YEARLY_C_TO_CROOT]  * s->value[DAILY_LIVEWOOD_TURNOVER];

		/* branch */
		s->value[C_BRANCH_LIVEWOOD_TO_DEADWOOD] = s->value[YEARLY_C_TO_BRANCH] * s->value[DAILY_LIVEWOOD_TURNOVER];

		/* reset annual values once used */
		s->value[YEARLY_C_TO_STEM]   = 0.;
		s->value[YEARLY_C_TO_CROOT]  = 0.;
		s->value[YEARLY_C_TO_BRANCH] = 0.;
	}

	/* check */
	CHECK_CONDITION (s->value[C_STEM_LIVEWOOD_TO_DEADWOOD],   <, 0.);
	CHECK_CONDITION (s->value[C_CROOT_LIVEWOOD_TO_DEADWOOD],  <, 0.);
	CHECK_CONDITION (s->value[C_BRANCH_LIVEWOOD_TO_DEADWOOD], <, 0.);
	CHECK_CONDITION (s->value[C_STEM_LIVEWOOD_TO_DEADWOOD],   >, s->value[STEM_SAPWOOD_C]);
	CHECK_CONDITION (s->value[C_CROOT_LIVEWOOD_TO_DEADWOOD],  >, s->value[CROOT_SAPWOOD_C]);
	CHECK_CONDITION (s->value[C_BRANCH_LIVEWOOD_TO_DEADWOOD], >, s->value[BRANCH_SAPWOOD_C]);

	/*******************************************************************************************************************/

	/* daily stem turnover live to dead wood*/

	s->value[STEM_SAPWOOD_C]   -= s->value[C_STEM_LIVEWOOD_TO_DEADWOOD];
	s->value[STEM_HEARTWOOD_C] += s->value[C_STEM_LIVEWOOD_TO_DEADWOOD];

	/*******************************************************************************************************************/

	/* daily coarse root turnover live to dead wood */

	s->value[CROOT_SAPWOOD_C]   -= s->value[C_CROOT_LIVEWOOD_TO_DEADWOOD];
	s->value[CROOT_HEARTWOOD_C] += s->value[C_CROOT_LIVEWOOD_TO_DEADWOOD];

	/*******************************************************************************************************************/

	/* daily branch turnover live to dead wood */

	s->value[BRANCH_SAPWOOD_C]   -= s->value[C_BRANCH_LIVEWOOD_TO_DEADWOOD];
	s->value[BRANCH_HEARTWOOD_C] += s->value[C_BRANCH_LIVEWOOD_TO_DEADWOOD];

	/*******************************************************************************************************************/

	/* Overall sapwood and heartwood stocks */

	s->value[TOT_SAPWOOD_C]       = s->value[STEM_SAPWOOD_C]   + s->value[CROOT_SAPWOOD_C]   + s->value[BRANCH_SAPWOOD_C];
	s->value[TOT_HEARTWOOD_C]     = s->value[STEM_HEARTWOOD_C] + s->value[CROOT_HEARTWOOD_C] + s->value[BRANCH_HEARTWOOD_C];

	/*******************************************************************************************************************/
  

}

