/*
 * tree_model_check_balance.c
 *
 *  Created on: 03 apr 2017
 *      Author: alessio
 */
/* includes */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "common.h"
#include "check_balance.h"
#include "constants.h"
#include "settings.h"
#include "logger.h"

extern logger_t* g_debug_log;

/* note: when model checks fluxes there's no need to use global variables */
/* note: when model checks for stocks and it uses at least one stock for balance or it
 *  simulates more than one class it needs to use global variables */


int check_tree_class_radiation_flux_balance(cell_t *const c, const int layer, const int height, const int dbh, const int age, const int species)
{
	double in;
	double out;
	double store;
	double balance;

	species_t *s;
	s = &c->heights[height].dbhs[dbh].ages[age].species[species];

	/* DAILY CHECK ON CLASS LEVEL CANOPY POOL-ATMOSPHERE RADIATIVE BALANCE */

	/* PAR balance */
	/* sum of sources */
	in = s->value[PAR];

	/* sum of sinks */
	out = s->value[PAR_REFL] + s->value[TRANSM_PAR];

	/* sum of current storage */
	store = s->value[APAR_SUN] + s->value[APAR_SHADE];

	/* check canopy water pool balance */
	balance = in - out - store;

	logger(g_debug_log, "\nCLASS LEVEL PAR BALANCE\n");

	/* check for PAR balance closure*/
	if ( fabs( balance ) > eps )
	{
		error_log("DOY = %d\n", c->doy);
		error_log("PAR in = %g\n", in);
		error_log("PAR out = %g\n", out);
		error_log("PAR store = %g\n", store);
		error_log("PAR balance = %g\n", balance);
		error_log("...FATAL ERROR in 'Tree_model_daily' PAR balance (exit)\n");
		CHECK_CONDITION(fabs( balance ), > , eps);

		return 0;
	}
	else
	{
		logger(g_debug_log, "...ok 'Tree_model_daily' PAR balance\n");
	}
	/****************************************************************************************************************/
	/* Net Short-Wave radiation balance */
	/* sum of sources */
	in      = s->value[SW_RAD];

	/* sum of sinks */
	/* it must take into account the overall transmitted NET_RAD (reflected is yet computed for net radiation) */
	out     = s->value[SW_RAD_REFL] + s->value[SW_RAD_TRANSM];

	/* sum of current storage */
	store   = s->value[SW_RAD_ABS_SUN] + s->value[SW_RAD_ABS_SHADE];

	/* check canopy water pool balance */
	balance = in - out - store;

	logger(g_debug_log, "\nCLASS LEVEL Short Wave BALANCE\n");

	/* check for NET_RAD balance closure*/
	if ( fabs( balance ) > eps )
	{
		error_log("DOY = %d\n", c->doy);
		error_log("NET_SW_RAD in = %g\n", in);
		error_log("NET_RAD out = %g\n", out);
		error_log("NET_RAD store = %g\n", store);
		error_log("NET_RAD balance = %g\n", balance);
		error_log("...FATAL ERROR AT CLASS LEVEL NET_RAD balance (exit)\n");
		CHECK_CONDITION(fabs( balance ), > , eps);

		return 0;
	}
	else
	{
		logger(g_debug_log, "...ok in 'Tree_model_daily' Short Wave radiation balance \n");
	}
	/****************************************************************************************************************/
	/* PPFD balance */
	/* sum of sources */
	in      = s->value[PPFD];

	/* sum of sinks */
	out     = s->value[PPFD_REFL] +  s->value[PPFD_TRANSM];

	/* sum of current storage */
	store   = s->value[PPFD_ABS_SUN] + s->value[PPFD_ABS_SHADE];

	/* check canopy water pool balance */
	balance = in - out - store;

	logger(g_debug_log, "\nCLASS LEVEL PPFD BALANCE\n");

	/* check for PPFD balance closure*/
	if ( fabs( balance ) > eps )
	{
		error_log("DOY = %d\n", c->doy);
		error_log("PPFD in = %g\n", in);
		error_log("PPFD out = %g\n", out);
		error_log("PPFD store = %g\n", store);
		error_log("PPFD balance = %g\n", balance);
		error_log("...FATAL ERROR in 'Tree_model_daily' PPFD balance (exit)\n");
		CHECK_CONDITION(fabs( balance ), > , eps);

		return 0;
	}
	else
	{
		logger(g_debug_log, "...ok in 'Tree_model_daily' PPFD balance\n");
	}
	/* ok */
	return 1;
}

int check_tree_class_carbon_flux_balance(cell_t *const c, const int layer, const int height, const int dbh, const int age, const int species)
{
	double in;
	double out;
	double store;
	double balance;

	species_t *s;
	s = &c->heights[height].dbhs[dbh].ages[age].species[species];

	/* DAILY CHECK ON CLASS LEVEL CARBON FLUX BALANCE */
	/* check complete tree level carbon flux balance */

	/* sum of sources */
	in      = s->value[GPP];

	/* sum of sinks */
	out     = s->value[TOTAL_GROWTH_RESP] + s->value[TOTAL_MAINT_RESP];

	/* sum of current storage */
	store   = s->value[NPP];

	/* check carbon flux balance */
	balance = in - out - store;

	logger(g_debug_log, "\nCLASS LEVEL CARBON FLUX BALANCE (first)\n");

	/* check for carbon flux balance closure */
	if ( fabs( balance ) > eps )
	{
		error_log("DOY = %d\n",                      c->doy);
		error_log("\nin = %g gC/m2/day\n",           in);
		error_log("GPP = %g gC/m2/day\n",            s->value[GPP]);
		error_log("\nout = %g gC/m2/day\n",          out);
		error_log("TOTAL_AUT_RESP = %g gC/m2/day\n", s->value[TOTAL_AUT_RESP]);
		error_log("\nstore = %g gC/m2/day\n",        store);
		error_log("NPP = %g gC/m2/day\n",            s->value[NPP]);
		error_log("\nbalance = %g gC/m2/day\n",      balance);
		error_log("...FATAL ERROR in 'Tree_model_daily' carbon flux balance (first) (exit)\n");
		CHECK_CONDITION(fabs( balance ), > , eps);

		return 0;
	}
	else
	{
		logger(g_debug_log, "...ok in 'Tree_model_daily' carbon flux balance (first)\n");
	}

	/* DAILY CHECK ON CLASS LEVEL CARBON BALANCE */
	/* check complete tree level carbon balance */

	/* sum of sources */
	in      = s->value[GPP_tC];

	/* sum of sinks */
	out     = s->value[TOTAL_MAINT_RESP_tC] +
			s->value[TOTAL_GROWTH_RESP_tC]  +
			s->value[C_LEAF_TO_LITR]        +
			s->value[C_FROOT_TO_LITR]       +
			s->value[C_BRANCH_TO_CWD]       +
			s->value[DEAD_LEAF_C]           +
			s->value[DEAD_FROOT_C]          +
			s->value[DEAD_STEM_C]           +
			s->value[DEAD_CROOT_C]          +
			s->value[DEAD_BRANCH_C]         +
			s->value[DEAD_RESERVE_C]        +
			s->value[DEAD_FRUIT_C]          ;

	/* sum of current storage */
	store   = s->value[C_TO_LEAF]  +
			s->value[C_TO_STEM]    +
			s->value[C_TO_FROOT]   +
			s->value[C_TO_CROOT]   +
			s->value[C_TO_BRANCH]  +
			s->value[C_TO_RESERVE] +
			s->value[C_TO_FRUIT]   ;

	/* check carbon flux balance */
	balance = in - out - store;

	logger(g_debug_log, "\nCLASS LEVEL CARBON FLUX BALANCE (second)\n");

	/* check for carbon flux balance closure */
	if ( fabs( balance ) > eps )
	{
		error_log("DOY = %d\n",                           c->doy);
		error_log("\nin = %g tC/cell/day\n",              in);
		error_log("GPP_gC = %g tC/cell/day\n",            s->value[GPP_tC]);
		error_log("\nout = %g tC/cell/day\n",             out);
		error_log("TOTAL_MAINT_RESP = %g tC/cell/day\n",  s->value[TOTAL_MAINT_RESP_tC]);
		error_log("TOTAL_GROWTH_RESP = %g tC/cell/day\n", s->value[TOTAL_GROWTH_RESP_tC]);
		error_log("C_LEAF_TO_LITR = %g tC/cell/day\n",    s->value[C_LEAF_TO_LITR]);
		error_log("C_FROOT_TO_LITR = %g tC/cell/day\n",   s->value[C_FROOT_TO_LITR]);
		error_log("C_BRANCH_TO_CWD = %g tC/cell/day\n",   s->value[C_BRANCH_TO_CWD]);
		error_log("DEAD_LEAF_C = %g tC/cell/day\n",       s->value[DEAD_LEAF_C]);
		error_log("DEAD_FROOT_C = %g tC/cell/day\n",      s->value[DEAD_FROOT_C]);
		error_log("DEAD_STEM_C = %g tC/cell/day\n",       s->value[DEAD_STEM_C]);
		error_log("DEAD_CROOT_C = %g tC/cell/day\n",      s->value[DEAD_CROOT_C]);
		error_log("DEAD_BRANCH_C = %g tC/cell/day\n",     s->value[DEAD_BRANCH_C]);
		error_log("DEAD_RESERVE_C = %g tC/cell/day\n",    s->value[DEAD_RESERVE_C]);
		error_log("DEAD_FRUIT_C = %g tC/cell/day\n",      s->value[DEAD_FRUIT_C]);
		error_log("\nstore = %g tC/cell/day\n",           store);
		error_log("C_TO_LEAF = %g tC/cell/day\n",         s->value[C_TO_LEAF]);
		error_log("C_TO_FROOT = %g tC/cell/day\n",        s->value[C_TO_FROOT]);
		error_log("C_TO_STEM = %g tC/cell/day\n",         s->value[C_TO_STEM]);
		error_log("C_TO_CROOT = %g tC/cell/day\n",        s->value[C_TO_CROOT]);
		error_log("C_TO_BRANCH = %g tC/cell/day\n",       s->value[C_TO_BRANCH]);
		error_log("C_TO_RESERVE = %g tC/cell/day\n",      s->value[C_TO_RESERVE]);
		error_log("C_TO_FRUIT = %g tC/cell/day\n",        s->value[C_TO_FRUIT]);
		error_log("\nbalance = %g tC/cell/day\n",         balance);
		error_log("...FATAL ERROR in 'Tree_model_daily'  carbon flux balance (second) (exit)\n");
		CHECK_CONDITION(fabs( balance ), > , eps);

		return 0;
	}
	else
	{
		logger(g_debug_log, "...ok in 'Tree_model_daily' carbon flux balance (second)\n");
	}
	/* ok */
	return 1;
}

int check_tree_class_carbon_mass_balance ( cell_t *const c, const int layer, const int height, const int dbh, const int age, const int species )
{
	species_t *s;
	s = &c->heights[height].dbhs[dbh].ages[age].species[species];

	/* DAILY CHECK ON CLASS LEVEL CARBON MASS BALANCE */
	/* check complete tree level carbon mass balance */

	/* sum of sources */
	s->value[TREEC_IN]    = s->value[GPP_tC];

	/* sum of sinks */
	s->value[TREEC_OUT]   = s->value[TOTAL_MAINT_RESP_tC] +
			s->value[TOTAL_GROWTH_RESP_tC]                +
			s->value[C_LEAF_TO_LITR]                      +
			s->value[C_FROOT_TO_LITR]                     +
			s->value[C_BRANCH_TO_CWD]                     +
			s->value[DEAD_LEAF_C]                         +
			s->value[DEAD_FROOT_C]                        +
			s->value[DEAD_STEM_C]                         +
			s->value[DEAD_CROOT_C]                        +
			s->value[DEAD_BRANCH_C]                       +
			s->value[DEAD_RESERVE_C]                      +
			s->value[DEAD_FRUIT_C]                        ;

	/* sum of current storage */
	s->value[TREEC_STORE] = s->value[LEAF_C] +
			s->value[FROOT_C]                +
			s->value[CROOT_C]                +
			s->value[STEM_C]                 +
			s->value[BRANCH_C]               +
			s->value[RESERVE_C]              +
			s->value[FRUIT_C]                ;

	/* check carbon pool balance */
	s->value[TREEC_BALANCE] = s->value[TREEC_IN] - s->value[TREEC_OUT] - (s->value[TREEC_STORE] - s->value[TREEC_OLDSTORE]);

	logger(g_debug_log, "\nCLASS LEVEL CARBON MASS BALANCE\n");

	/* check for carbon mass balance closure */
	if ( ( fabs( s->value[TREEC_BALANCE] ) > eps ) && ( s->counter[DOS] > 1 ) )
	{
		error_log("TREE DOS = %d\n",                         s->counter[DOS]);
		error_log("\nin = %g tC/cell/day\n",                 s->value[TREEC_IN]);
		error_log("GPP = %g tC/cell\n",                      s->value[GPP_tC]);
		error_log("\nout = %g tC/cell/day\n",                s->value[TREEC_OUT]);
		error_log("TOTAL_MAINT_RESP_tC = %g tC/cell/day\n",  s->value[TOTAL_MAINT_RESP_tC]);
		error_log("TOTAL_GROWTH_RESP_tC = %g tC/cell/day\n", s->value[TOTAL_GROWTH_RESP_tC]);
		error_log("C_LEAF_TO_LITR = %g tC/cell/day\n",       s->value[C_LEAF_TO_LITR]);
		error_log("C_FROOT_TO_LITR = %g tC/cell/day\n",      s->value[C_FROOT_TO_LITR]);
		error_log("C_BRANCH_TO_CWD = %g tC/cell/day\n",      s->value[C_BRANCH_TO_CWD]);
		error_log("DEAD_LEAF_C = %g tC/cell/day\n",          s->value[DEAD_LEAF_C]);
		error_log("DEAD_FROOT_C = %g tC/cell/day\n",         s->value[DEAD_FROOT_C]);
		error_log("DEAD_STEM_C = %g tC/cell/day\n",          s->value[DEAD_STEM_C]);
		error_log("DEAD_CROOT_C = %g tC/cell/day\n",         s->value[DEAD_CROOT_C]);
		error_log("DEAD_BRANCH_C = %g tC/cell/day\n",        s->value[DEAD_BRANCH_C]);
		error_log("DEAD_RESERVE_C = %g tC/cell/day\n",       s->value[DEAD_RESERVE_C]);
		error_log("DEAD_FRUIT_C = %g tC/cell/day\n",         s->value[DEAD_FRUIT_C]);
		error_log("\nold_store = %g tC/cell\n",              s->value[TREEC_OLDSTORE]);
		error_log("store = %g tC/cell\n",                    s->value[TREEC_STORE]);
		error_log("store - old_store = %g tC/cell\n",        s->value[TREEC_STORE] - s->value[TREEC_OLDSTORE]);
		error_log("LEAF_C = %g tC/cell/day\n",               s->value[LEAF_C]);
		error_log("FROOT_C = %g tC/cell/day\n",              s->value[FROOT_C]);
		error_log("CROOT_C = %g tC/cell/day\n",              s->value[CROOT_C]);
		error_log("STEM_C = %g tC/cell/day\n",               s->value[STEM_C]);
		error_log("BRANCH_C = %g tC/cell/day\n",             s->value[BRANCH_C]);
		error_log("RESERVE_C = %g tC/cell/day\n",            s->value[RESERVE_C]);
		error_log("FRUIT_C = %g tC/cell/day\n",              s->value[FRUIT_C]);
		error_log("\nbalance = %g tC/cell\n",                s->value[TREEC_BALANCE]);
		error_log("...FATAL ERROR in 'Tree_model_daily' carbon mass balance (exit)\n");
		CHECK_CONDITION(fabs( s->value[TREEC_BALANCE] ), > , eps);

		return 0;
	}
	else
	{
		s->value[TREEC_OLDSTORE] = s->value[TREEC_STORE];
		logger(g_debug_log, "...ok in 'Tree_model_daily' carbon mass balance\n");
	}
	/* ok */
	return 1;
}
int check_tree_class_nitrogen_flux_balance ( cell_t *const c, const int layer, const int height, const int dbh, const int age, const int species )
{
	double in;
	double out;
	double store;
	double balance;

	species_t *s;
	s = &c->heights[height].dbhs[dbh].ages[age].species[species];

	/* DAILY CHECK ON CLASS LEVEL NITROGEN BALANCE */
	/* check complete tree level nitrogen balance */

	/* sum of sources */
	in      = s->value[NPP_tN];

	/* sum of sinks */
	out     =  s->value[N_LEAF_TO_LITR] +
			s->value[N_FROOT_TO_LITR]   +
			s->value[N_STEM_TO_CWD]     +
			s->value[N_CROOT_TO_CWD]    +
			s->value[N_BRANCH_TO_CWD]   +
			s->value[N_RESERVE_TO_CWD]  +
			s->value[N_FRUIT_TO_CWD]    ;

	/* sum of current storage */
	store   = s->value[N_TO_LEAF] +
			s->value[N_TO_STEM]   +
			s->value[N_TO_FROOT]  +
			s->value[N_TO_CROOT]  +
			s->value[N_TO_BRANCH] +
			s->value[N_TO_RESERVE]+
			s->value[N_TO_FRUIT]  ;

	/* check nitrogen flux balance */
	balance = in - out - store;

	logger(g_debug_log, "\nCLASS LEVEL NITROGEN FLUX BALANCE\n");

	/* check for nitrogen flux balance closure */
	if ( fabs( balance ) > eps )
	{
		error_log("DOY = %d\n",                          c->doy);
		error_log("\nin = %g tN/cell/day\n",             in);
		error_log("NPP_tN = %g tN/cell/day\n",           s->value[NPP_tN]);
		error_log("\nout = %g tN/cell/day\n",            out);
		error_log("N_LEAF_TO_LITR = %g tN/cell/day\n",   s->value[N_LEAF_TO_LITR]);
		error_log("N_FROOT_TO_LITR = %g tN/cell/day\n",  s->value[N_FROOT_TO_LITR]);
		error_log("N_STEM_TO_CWD = %g tN/cell/day\n",    s->value[N_STEM_TO_CWD]);
		error_log("N_CROOT_TO_CWD = %g tN/cell/day\n",   s->value[N_CROOT_TO_CWD]);
		error_log("N_BRANCH_TO_CWD = %g tN/cell/day\n",  s->value[N_BRANCH_TO_CWD]);
		error_log("N_RESERVE_TO_CWD = %g tN/cell/day\n", s->value[N_RESERVE_TO_CWD]);
		error_log("N_FRUIT_TO_CWD = %g tN/cell/day\n",   s->value[N_FRUIT_TO_CWD]);
		error_log("\nstore = %g tN/cell/day\n",          store);
		error_log("N_TO_LEAF = %g tN/cell/day\n",        s->value[N_TO_LEAF]);
		error_log("N_TO_FROOT = %g tN/cell/day\n",       s->value[N_TO_FROOT]);
		error_log("N_TO_CROOT = %g tN/cell/day\n",       s->value[N_TO_CROOT]);
		error_log("N_TO_STEM = %g tN/cell/day\n",        s->value[N_TO_STEM]);
		error_log("N_TO_RESERVE = %g tN/cell/day\n",     s->value[N_TO_RESERVE]);
		error_log("N_TO_BRANCH = %g tN/cell/day\n",      s->value[N_TO_BRANCH]);
		error_log("N_TO_FRUIT = %g tN/cell/day\n",       s->value[N_TO_FRUIT]);
		error_log("\nbalance = %g tN/cell/day\n",        balance);
		error_log("...FATAL ERROR in 'Tree_model_daily' nitrogen flux balance (exit)\n");
		CHECK_CONDITION(fabs( balance ), > , eps);

		return 0;
	}
	else
	{
		logger(g_debug_log, "...ok in 'Tree_model_daily' nitrogen flux balance\n");
	}
	/* ok */
	return 1;
}
int check_tree_class_nitrogen_mass_balance ( cell_t *const c, const int layer, const int height, const int dbh, const int age, const int species )
{
	species_t *s;
	s = &c->heights[height].dbhs[dbh].ages[age].species[species];

	/* DAILY CHECK ON CLASS LEVEL NITROGEN MASS BALANCE */
	/* check complete tree level nitrogen mass balance */

	/* sum of sources */
	s->value[TREEN_IN]    = /*s->value[GPP_tN]*/0.;

	/* sum of sinks */
	s->value[TREEN_OUT]   = s->value[C_LEAF_TO_LITR]      +
			s->value[N_FROOT_TO_LITR]                     +
			s->value[N_BRANCH_TO_CWD]                     +
			s->value[DEAD_LEAF_N]                         +
			s->value[DEAD_FROOT_N]                        +
			s->value[DEAD_STEM_N]                         +
			s->value[DEAD_CROOT_N]                        +
			s->value[DEAD_BRANCH_N]                       +
			s->value[DEAD_RESERVE_N]                      +
			s->value[DEAD_FRUIT_N]                        ;

	/* sum of current storage */
	s->value[TREEC_STORE] = s->value[LEAF_N] +
			s->value[FROOT_N]                +
			s->value[CROOT_N]                +
			s->value[STEM_N]                 +
			s->value[BRANCH_N]               +
			s->value[RESERVE_N]              +
			s->value[FRUIT_N]                ;

	/* check nitrogen pool balance */
	s->value[TREEN_BALANCE] = s->value[TREEN_IN] - s->value[TREEN_OUT] - (s->value[TREEN_STORE] - s->value[TREEN_OLDSTORE]);

	logger(g_debug_log, "\nCLASS LEVEL CARBON MASS BALANCE\n");

	/* check for nitrogen mass balance closure */
	if ( ( fabs( s->value[TREEN_BALANCE] ) > eps ) && ( s->counter[DOS] > 1 ) )
	{
		error_log("TREE DOS = %d\n",                         s->counter[DOS]);
		error_log("\nin = %g tN/cell/day\n",                 s->value[TREEN_IN]);
		error_log("\nout = %g tN/cell/day\n",                s->value[TREEN_OUT]);
		error_log("N_LEAF_TO_LITR = %g tN/cell/day\n",       s->value[N_LEAF_TO_LITR]);
		error_log("N_FROOT_TO_LITR = %g tN/cell/day\n",      s->value[N_FROOT_TO_LITR]);
		error_log("N_BRANCH_TO_CWD = %g tN/cell/day\n",      s->value[N_BRANCH_TO_CWD]);
		error_log("DEAD_LEAF_N = %g tN/cell/day\n",          s->value[DEAD_LEAF_N]);
		error_log("DEAD_FROOT_N = %g tN/cell/day\n",         s->value[DEAD_FROOT_N]);
		error_log("DEAD_STEM_N = %g tN/cell/day\n",          s->value[DEAD_STEM_N]);
		error_log("DEAD_NROOT_N = %g tN/cell/day\n",         s->value[DEAD_CROOT_N]);
		error_log("DEAD_BRANCH_N = %g tN/cell/day\n",        s->value[DEAD_BRANCH_N]);
		error_log("DEAD_RESERVE_N = %g tN/cell/day\n",       s->value[DEAD_RESERVE_N]);
		error_log("DEAD_FRUIT_N = %g tN/cell/day\n",         s->value[DEAD_FRUIT_N]);
		error_log("\nold_store = %g tN/cell\n",              s->value[TREEN_OLDSTORE]);
		error_log("store = %g tN/cell\n",                    s->value[TREEN_STORE]);
		error_log("store - old_store = %g tN/cell\n",        s->value[TREEN_STORE] - s->value[TREEN_OLDSTORE]);
		error_log("LEAF_N = %g tN/cell/day\n",               s->value[LEAF_N]);
		error_log("FROOT_N = %g tN/cell/day\n",              s->value[FROOT_N]);
		error_log("CROOT_N = %g tN/cell/day\n",              s->value[CROOT_N]);
		error_log("STEM_N = %g tN/cell/day\n",               s->value[STEM_N]);
		error_log("BRANCH_N = %g tN/cell/day\n",             s->value[BRANCH_N]);
		error_log("RESERVE_N = %g tN/cell/day\n",            s->value[RESERVE_N]);
		error_log("FRUIT_N = %g tN/cell/day\n",              s->value[FRUIT_N]);
		error_log("\nbalance = %g tN/cell\n",                s->value[TREEN_BALANCE]);
		error_log("...FATAL ERROR in 'Tree_model_daily' nitrogen mass balance (exit)\n");
		CHECK_CONDITION(fabs( s->value[TREEN_BALANCE] ), > , eps);

		return 0;
	}
	else
	{
		s->value[TREEN_OLDSTORE] = s->value[TREEN_STORE];
		logger(g_debug_log, "...ok in 'Tree_model_daily' nitrogen mass balance\n");
	}
	/* ok */
	return 1;
}

int check_tree_class_water_flux_balance ( cell_t *const c, const int layer, const int height, const int dbh, const int age, const int species )
{
	species_t *s;
	s = &c->heights[height].dbhs[dbh].ages[age].species[species];

	/* DAILY CHECK ON CLASS LEVEL CANOPY POOL-ATMOSPHERE WATER BALANCE */

	/* sum of sources */
	s->value[TREEW_IN]      = s->value[CANOPY_INT];

	/* sum of sinks */
	s->value[TREEW_OUT]     = s->value[CANOPY_EVAPO];

	/* sum of current storage */
	s->value[TREEW_STORE]   = s->value[CANOPY_WATER];

	/* check canopy water pool balance */
	s->value[TREEW_BALANCE] = s->value[TREEW_IN] - s->value[TREEW_OUT] - ( s->value[TREEW_STORE] - s->value[TREEW_OLDSTORE] );

	logger(g_debug_log, "\nCLASS LEVEL WATER BALANCE\n");

	if ( ( fabs( s->value[TREEW_BALANCE] ) > eps ) && ( s->counter[DOS] > 1 ) )
	{
		error_log("TREE DOS = %d\n",             s->counter[DOS]);
		error_log("canopy water in = %g\n",      s->value[TREEW_IN]);
		error_log("canopy water out = %g\n",     s->value[TREEW_OUT]);
		error_log("canopy water store = %g\n",   s->value[TREEW_STORE]);
		error_log("canopy water balance = %g\n", s->value[TREEW_BALANCE]);
		error_log("...FATAL ERROR in 'Tree_model_daily' canopy water balance (exit)\n");
		CHECK_CONDITION(fabs( s->value[TREEW_BALANCE] ), > , eps);

		return 0;
	}
	else
	{
		s->value[TREEW_OLDSTORE] = s->value[TREEW_STORE];
		logger(g_debug_log, "...ok in 'Tree_model_daily' canopy water balance\n");
	}
	/* ok */
	return 1;
}

