/*utility.c*/

/* includes */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utility.h"
#include "constants.h"
#include "logger.h"

extern logger_t* g_log;

void Reset_daily_variables (cell_t *const c, const int count)
{
	int height;
	int age;
	int species;
	int i;

	logger(g_log, "...resetting cell level daily variables...\n");

	/* reset daily radiative variables */
	c->apar = 0.;
	c->par_transm = 0.;
	c->par_reflected = 0.;
	c->par_for_soil = 0.;
	c->par_reflected_soil = 0.;

	c->sw_rad_refl = 0.;
	c->sw_rad_for_soil_refl = 0.;
	c->net_sw_rad_abs = 0.;
	c->net_sw_rad_transm = 0.;
	c->net_sw_rad_for_soil = 0.;

	c->ppfd_abs = 0.;
	c->ppfd_transm = 0.;
	c->ppfd_reflected = 0.;
	c->ppfd_for_soil = 0.;
	c->ppfd_reflected_soil = 0.;

	c->canopy_temp_k = 0.;


	/*reset daily carbon variables*/
	c->daily_gpp = 0.;
	c->daily_npp_tDM = 0.;
	c->daily_npp_gC = 0.;
	c->daily_leaf_maint_resp = 0.;
	c->daily_stem_maint_resp = 0.;
	c->daily_fine_root_maint_resp = 0.;
	c->daily_branch_maint_resp = 0.;
	c->daily_maint_resp = 0.;
	c->daily_coarse_root_maint_resp = 0.;
	c->daily_leaf_growth_resp = 0.;
	c->daily_stem_growth_resp = 0.;
	c->daily_fine_root_growth_resp = 0.;
	c->daily_branch_growth_resp = 0.;
	c->daily_coarse_root_growth_resp = 0.;
	c->daily_growth_resp = 0.;
	c->daily_aut_resp = 0.;
	c->daily_aut_resp_tC = 0.;
	c->daily_leaf_aut_resp = 0.;
	c->daily_stem_aut_resp = 0.;
	c->daily_branch_aut_resp = 0.;
	c->daily_fine_root_aut_resp = 0.;
	c->daily_coarse_root_aut_resp = 0.;
	c->daily_leaf_carbon = 0.;
	c->daily_stem_carbon = 0.;
	c->daily_fine_root_carbon = 0.;
	c->daily_coarse_root_carbon = 0.;
	c->daily_root_carbon = 0.;
	c->daily_branch_carbon = 0.;
	c->daily_reserve_carbon = 0.;
	c->daily_litter_carbon = 0.;
	c->daily_fruit_carbon = 0.;
	c->daily_litterfall = 0.;

	/*reset daily water variables*/
	c->prcp_rain = 0.;
	c->prcp_snow = 0.;
	c->snow_melt = 0.;
	c->snow_subl = 0.;
	c->out_flow = 0.;

	c->daily_c_transp = 0.;
	c->daily_c_int = 0.;
	c->daily_c_evapo = 0.;
	c->daily_soil_evapo = 0.;
	c->daily_c_bl_cond = 0.;
	c->daily_c_sensible_heat_flux = 0.;
	c-> daily_soil_latent_heat_flux = 0.;
	c->daily_latent_heat_flux = 0.;
	c->daily_sensible_heat_flux = 0.;

	c->daily_litterfall = 0.;
	c->cell_cover = 0.;
	//ALESSIOC
	//c->dominant_veg_counter = 0;
	//c->dominated_veg_counter = 0;
	//c->subdominated_veg_counter = 0;

	//c->layer_daily_dead_tree[0] = 0;
	//c->layer_daily_dead_tree[1] = 0;
	//c->layer_daily_dead_tree[2] = 0;
	//c->daily_dead_tree = 0;

	logger(g_log, "...resetting class level daily variables...\n");

	for ( height = count - 1; height >= 0; height-- )
	{
		i = c->heights[height].z;

		for (age = c->heights[height].ages_count - 1; age >= 0; age --)
		{
			for (species = c->heights[height].ages[age].species_count - 1; species >= 0; species -- )
			{
				/* reset daily radiative variables */
				c->heights[height].ages[age].species[species].value[PAR] = 0.;
				c->heights[height].ages[age].species[species].value[APAR] = 0.;
				c->heights[height].ages[age].species[species].value[APAR_SUN] = 0.;
				c->heights[height].ages[age].species[species].value[APAR_SHADE] = 0.;
				c->heights[height].ages[age].species[species].value[TRANSM_PAR] = 0.;
				c->heights[height].ages[age].species[species].value[TRANSM_PAR_SUN] = 0.;
				c->heights[height].ages[age].species[species].value[TRANSM_PAR_SHADE] = 0.;
				c->heights[height].ages[age].species[species].value[NET_SW_RAD] = 0.;
				c->heights[height].ages[age].species[species].value[NET_SW_RAD_ABS] = 0.;
				c->heights[height].ages[age].species[species].value[NET_SW_RAD_ABS_SUN] = 0.;
				c->heights[height].ages[age].species[species].value[NET_SW_RAD_ABS_SHADE] = 0.;
				c->heights[height].ages[age].species[species].value[NET_SW_RAD_TRANSM] = 0.;
				c->heights[height].ages[age].species[species].value[NET_SW_RAD_TRANSM_SUN] = 0.;
				c->heights[height].ages[age].species[species].value[NET_SW_RAD_TRANSM_SHADE] = 0.;
				c->heights[height].ages[age].species[species].value[PPFD] = 0.;
				c->heights[height].ages[age].species[species].value[PPFD_ABS] = 0.;
				c->heights[height].ages[age].species[species].value[PPFD_ABS_SUN] = 0.;
				c->heights[height].ages[age].species[species].value[PPFD_ABS_SHADE] = 0.;
				c->heights[height].ages[age].species[species].value[PPFD_TRANSM] = 0.;
				c->heights[height].ages[age].species[species].value[PPFD_TRANSM_SUN] = 0.;
				c->heights[height].ages[age].species[species].value[PPFD_TRANSM_SHADE] = 0.;

				/* reset daily carbon fluxes */
				c->heights[height].ages[age].species[species].value[DAILY_GPP_gC] = 0.;
				c->heights[height].ages[age].species[species].value[DAILY_POINT_GPP_gC] = 0.;
				c->heights[height].ages[age].species[species].value[NPP_gC] = 0.;
				c->heights[height].ages[age].species[species].value[NPP_tDM] = 0.;
				c->heights[height].ages[age].species[species].value[C_FLUX] = 0.;

				/* reset daily water fluxes */
				c->heights[height].ages[age].species[species].value[CANOPY_INT] = 0.;
				c->heights[height].ages[age].species[species].value[CANOPY_EVAPO] = 0.;
				c->heights[height].ages[age].species[species].value[CANOPY_TRANSP] = 0.;
				c->heights[height].ages[age].species[species].value[CANOPY_EVAPO_TRANSP] = 0.;

				/* reset daily multipliers */
				c->heights[height].ages[age].species[species].value[F_CO2] = 0.;
				c->heights[height].ages[age].species[species].value[F_LIGHT] = 0.;
				c->heights[height].ages[age].species[species].value[F_T] = 0.;
				c->heights[height].ages[age].species[species].value[F_FROST] = 0.;
				c->heights[height].ages[age].species[species].value[F_VPD] = 0.;
				c->heights[height].ages[age].species[species].value[F_AGE] = 0.;
				c->heights[height].ages[age].species[species].value[F_NUTR] = 0.;
				c->heights[height].ages[age].species[species].value[F_SW] = 0.;
				c->heights[height].ages[age].species[species].value[F_PSI] = 0.;
				c->heights[height].ages[age].species[species].value[PHYS_MOD] = 0.;

				/* reset daily carbon fluxes among pools */
				c->heights[height].ages[age].species[species].value[C_TO_LEAF] = 0.;
				c->heights[height].ages[age].species[species].value[C_TO_ROOT] = 0.;
				c->heights[height].ages[age].species[species].value[C_TO_FINEROOT] = 0.;
				c->heights[height].ages[age].species[species].value[C_TO_COARSEROOT] = 0.;
				c->heights[height].ages[age].species[species].value[C_TO_TOT_STEM] = 0.;
				c->heights[height].ages[age].species[species].value[C_TO_STEM] = 0.;
				c->heights[height].ages[age].species[species].value[C_TO_BRANCH] = 0.;
				c->heights[height].ages[age].species[species].value[C_TO_RESERVE] = 0.;
				c->heights[height].ages[age].species[species].value[C_TO_FRUIT] = 0.;
				c->heights[height].ages[age].species[species].value[C_TO_LITTER] = 0.;
				c->heights[height].ages[age].species[species].value[C_LEAF_TO_RESERVE] = 0.;
				c->heights[height].ages[age].species[species].value[C_FINEROOT_TO_RESERVE] = 0.;

				/* reset daily maint and growth respiration */
				c->heights[height].ages[age].species[species].value[DAILY_LEAF_MAINT_RESP] = 0.;
				c->heights[height].ages[age].species[species].value[NIGHTLY_LEAF_MAINT_RESP] = 0.;
				c->heights[height].ages[age].species[species].value[TOT_DAY_LEAF_MAINT_RESP] = 0.;
				c->heights[height].ages[age].species[species].value[FINE_ROOT_MAINT_RESP] = 0.;
				c->heights[height].ages[age].species[species].value[STEM_MAINT_RESP] = 0.;
				c->heights[height].ages[age].species[species].value[BRANCH_MAINT_RESP] = 0.;
				c->heights[height].ages[age].species[species].value[COARSE_ROOT_MAINT_RESP] = 0.;
				c->heights[height].ages[age].species[species].value[TOTAL_MAINT_RESP] = 0.;
				c->heights[height].ages[age].species[species].value[LEAF_GROWTH_RESP] = 0.;
				c->heights[height].ages[age].species[species].value[FINE_ROOT_GROWTH_RESP] = 0.;
				c->heights[height].ages[age].species[species].value[COARSE_ROOT_GROWTH_RESP] = 0.;
				c->heights[height].ages[age].species[species].value[STEM_GROWTH_RESP] = 0.;
				c->heights[height].ages[age].species[species].value[BRANCH_GROWTH_RESP] = 0.;
				c->heights[height].ages[age].species[species].value[TOTAL_GROWTH_RESP] = 0.;

				/* reset layer level daily carbon biomass increment in tC/cell/day */
				//ALESSIOC
				//c->daily_delta_wts[i] = 0.;
				//c->daily_delta_ws[i] = 0.;
				//c->daily_delta_wf[i] = 0.;
				//c->daily_delta_wbb[i] = 0.;
				//c->daily_delta_wfr[i] = 0.;
				//c->daily_delta_wcr[i] = 0.;
				//c->daily_delta_wres[i] = 0.;
			}
		}
	}


}

void Reset_monthly_variables (cell_t *const c, const int count)
{
	int height;
	int age;
	int species;

	logger(g_log, "...resetting monthly variables...\n");

	c->basal_area = 0.;

	c->monthly_gpp = 0.;
	c->monthly_npp_gC = 0.;
	c->monthly_npp_tDM = 0.;
	c->monthly_aut_resp = 0.;
	c->monthly_aut_resp_tC = 0.;
	c->monthly_maint_resp = 0.;
	c->monthly_growth_resp = 0.;
	c->monthly_r_eco = 0.;
	c->monthly_het_resp = 0.;
	c->monthly_gpp = 0.;
	c->monthly_C_flux = 0.;
	c->monthly_nee = 0.;
	c->monthly_litterfall = 0.;
	c->monthly_tot_w_flux = 0.;
	c->monthly_c_int = 0.;
	c->monthly_c_transp = 0.;
	c->monthly_c_evapo = 0.;
	c->monthly_c_water_stored = 0.;
	c->monthly_c_evapotransp = 0.;
	c->monthly_soil_evapo = 0.;
	c->monthly_et = 0.;
	c->monthly_c_bl_cond = 0.;
	c->monthly_latent_heat_flux = 0.;
	c->monthly_sensible_heat_flux = 0.;

	//ALESSIOC
	//c->layer_monthly_gpp[2] = 0.;
	//tocontinue...

	for ( height = count - 1; height >= 0; height-- )
	{
		for (age = c->heights[height].ages_count - 1; age >= 0; age --)
		{
			for (species = c->heights[height].ages[age].species_count - 1; species >= 0; species -- )
			{
				c->heights[height].ages[age].species[species].value[MONTHLY_EVAPOTRANSPIRATION] = 0.;
			}
		}
	}
}


void Reset_annual_variables (cell_t *const c, const int count)
{
	int height;
	int age;
	int species;
	logger(g_log, "...resetting annual variables...\n");

	/*reset cell related variables*/
	// ALESSIOC
	//c->canopy_cover_dominant = 0.;
	//c->canopy_cover_dominated = 0.;
	//c->canopy_cover_subdominated = 0.;

	c->annual_gpp = 0.;
	c->annual_npp_gC = 0.;
	c->annual_npp_tDM = 0.;
	c->annual_aut_resp = 0.;
	c->annual_aut_resp_tC = 0.;
	c->annual_maint_resp = 0.;
	c->annual_growth_resp = 0.;
	c->annual_r_eco = 0.;
	c->annual_het_resp = 0.;
	c->annual_c_int = 0.;
	c->annual_c_transp = 0.;
	c->annual_c_evapo = 0.;
	c->annual_c_water_stored = 0.;
	c->annual_c_evapotransp = 0.;
	c->annual_soil_evapo = 0.;
	c->annual_et = 0.;
	c->annual_c_bl_cond = 0.;
	c->annual_latent_heat_flux = 0.;
	c->annual_sensible_heat_flux = 0.;


	c->stand_agb = 0.;
	c->stand_bgb = 0.;
	//c->dead_tree = 0;
	c->annual_soil_evapo = 0.;

	for ( height = count - 1; height >= 0; height-- )
	{
		for (age = c->heights[height].ages_count - 1; age >= 0; age --)
		{
			for (species = c->heights[height].ages[age].species_count - 1; species >= 0; species -- )
			{
				c->heights[height].ages[age].species[species].value[PEAK_LAI] = 0.;
				c->heights[height].ages[age].species[species].value[MAX_LEAF_C] = 0.;
				/*reset cumulative values*/

				c->heights[height].ages[age].species[species].counter[VEG_DAYS] = 0;
				c->heights[height].ages[age].species[species].value[YEARLY_PHYS_MOD] = 0;

				c->heights[height].ages[age].species[species].value[YEARLY_GPP_gC] = 0;
				c->heights[height].ages[age].species[species].value[YEARLY_POINT_GPP_gC] = 0;
				c->heights[height].ages[age].species[species].value[YEARLY_NPP_tDM] = 0;

				// ALESSIOR DEL_STEMS used instead of DEAD_STEMS
				c->heights[height].ages[age].species[species].counter[DEAD_STEMS] = 0;
				c->heights[height].ages[age].species[species].counter[N_TREE_SAP] = 0;

				//INITIALIZE AVERAGE YEARLY MODIFIERS
				c->heights[height].ages[age].species[species].value[AVERAGE_F_VPD]  = 0.;
				c->heights[height].ages[age].species[species].value[AVERAGE_F_T]  = 0.;
				c->heights[height].ages[age].species[species].value[AVERAGE_F_SW]  = 0.;
				c->heights[height].ages[age].species[species].value[F_AGE]  = 0.;

				c->heights[height].ages[age].species[species].value[DEL_Y_WS] = 0.;
				c->heights[height].ages[age].species[species].value[DEL_Y_WF] = 0.;
				c->heights[height].ages[age].species[species].value[DEL_Y_WFR] = 0.;
				c->heights[height].ages[age].species[species].value[DEL_Y_WCR] = 0.;
				c->heights[height].ages[age].species[species].value[DEL_Y_WRES] = 0.;
				c->heights[height].ages[age].species[species].value[DEL_Y_WR] = 0.;
				c->heights[height].ages[age].species[species].value[DEL_Y_BB] = 0.;

				//SERGIO
				c->fineRootBiomass = c->heights[height].ages[age].species[species].value[BIOMASS_FINE_ROOT_tDM];
				//fixme
				if (c->heights[height].ages[age].species[species].value[PHENOLOGY] == 0.1 || c->heights[height].ages[age].species[species].value[PHENOLOGY] == 0.2)
				{
					c->fineRootBiomass = 0.;
				}
				c->coarseRootBiomass =  c->heights[height].ages[age].species[species].value[BIOMASS_COARSE_ROOT_LIVE_WOOD_tDM];
				c->stemBranchBiomass =  c->heights[height].ages[age].species[species].value[BIOMASS_STEM_BRANCH_LIVE_WOOD_tDM];
				c->stemBiomass =  c->heights[height].ages[age].species[species].value[BIOMASS_STEM_LIVE_WOOD_tDM];
				c->heights[height].ages[age].species[species].value[OLD_BIOMASS_ROOTS_COARSE] =
						c->heights[height].ages[age].species[species].value[BIOMASS_COARSE_ROOT_tDM];
				c->heights[height].ages[age].species[species].value[OLD_BIOMASS_FINE_ROOT_tDM] =
						c->heights[height].ages[age].species[species].value[BIOMASS_FINE_ROOT_tDM];
				c->heights[height].ages[age].species[species].value[OLD_BIOMASS_STEM] =
						c->heights[height].ages[age].species[species].value[BIOMASS_STEM_tDM];
				c->heights[height].ages[age].species[species].value[OLD_BIOMASS_BRANCH] =
						c->heights[height].ages[age].species[species].value[BIOMASS_BRANCH_tDM];
				c->heights[height].ages[age].species[species].value[OLD_BIOMASS_LEAVES] =
						c->heights[height].ages[age].species[species].value[BIOMASS_FOLIAGE_tDM];
				c->heights[height].ages[age].species[species].value[OLD_BIOMASS_STEM_LIVE_WOOD] =
						c->heights[height].ages[age].species[species].value[BIOMASS_STEM_LIVE_WOOD_tDM];
				c->heights[height].ages[age].species[species].value[OLD_BIOMASS_COARSE_ROOT_LIVE_WOOD] =
						c->heights[height].ages[age].species[species].value[BIOMASS_COARSE_ROOT_LIVE_WOOD_tDM];
				c->heights[height].ages[age].species[species].value[OLD_BIOMASS_STEM_BRANCH_LIVE_WOOD] =
						c->heights[height].ages[age].species[species].value[BIOMASS_STEM_BRANCH_LIVE_WOOD_tDM];

			}
		}
	}
}

void First_day(cell_t *const c, const int count)
{
	int height;
	int age;
	int species;

	logger(g_log, "..first day..\n");

	c->days_since_rain = 0.;

	for ( height = count - 1; height >= 0; height-- )
	{
		for (age = c->heights[height].ages_count - 1; age >= 0; age --)
		{
			for (species = c->heights[height].ages[age].species_count - 1; species >= 0; species -- )
			{
				/* compute cell level number of trees */
				c->n_tree += c->heights[height].ages[age].species[species].counter[N_TREE];

				c->heights[height].ages[age].species[species].turnover->FINERTOVER = 365 /
						c->heights[height].ages[age].species[species].value[LEAF_FINEROOT_TURNOVER];
				c->heights[height].ages[age].species[species].turnover->COARSERTOVER = 365 /
						c->heights[height].ages[age].species[species].value[COARSEROOT_TURNOVER];
				c->heights[height].ages[age].species[species].turnover->STEMTOVER = 365 /
						c->heights[height].ages[age].species[species].value[LIVE_WOOD_TURNOVER];
				c->heights[height].ages[age].species[species].turnover->BRANCHTOVER = 365 /
						c->heights[height].ages[age].species[species].value[BRANCHTTOVER];

				/* compute value for volume for next years comparisons (CAI-MAI) */
				c->heights[height].ages[age].species[species].value[MASS_DENSITY] = c->heights[height].ages[age].species[species].value[RHOMAX] +
						(c->heights[height].ages[age].species[species].value[RHOMIN] - c->heights[height].ages[age].species[species].value[RHOMAX]) *
						exp(-ln2 * (c->heights[height].ages[age].value / c->heights[height].ages[age].species[species].value[TRHO]));
				c->heights[height].ages[age].species[species].value[PREVIOUS_VOLUME] = c->heights[height].ages[age].species[species].value[STEM_C] * GC_GDM *
						(1 - c->heights[height].ages[age].species[species].value[FRACBB]) / c->heights[height].ages[age].species[species].value[MASS_DENSITY];
			}
		}
	}
}



