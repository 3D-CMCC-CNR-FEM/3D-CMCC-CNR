/* matrix.h */
#ifndef MATRIX_H_
#define MATRIX_H_

#include "meteo.h"
#include "soil_settings.h"

#define MAXTURNTIME 5000

typedef enum {
	T		/* timber */ // FIXME change with H (high forest)
	, C		/* coppice */
	, N		/* natural (no management) */
} e_management;

enum {

	/*********************************************** PARAMETERS **********************************************/
	/*********************************************************************************************************/

	/* SPECIES-SPECIFIC CONSTANT ECO-PHYSIOLOGICAL PARAMETER VALUES */
	/* NOTE: DON'T CHANGE THEIR ORDER!! */

	LIGHT_TOL,                           /* Light Tolerance 4 = very shade intolerant (cc = 90%), 3 = shade intolerant (cc = 100%), 2 = shade tolerant (cc = 110%), 1 = very shade tolerant (cc = 120%)*/
	PHENOLOGY,                           /* PHENOLOGY 0.1 = deciduous broadleaf, 0.2 = deciduous needle leaf, 1.1 = broad leaf evergreen, 1.2 = needle leaf evergreen*/
	ALPHA,                               /* Canopy quantum efficiency (molC/molPAR) */
	EPSILONgCMJ,                         /* Light Use Efficiency  (gC/MJ)(used if ALPHA is not available) */
	GAMMA_LIGHT,                         /* Empirical parameter for Light modifiers (see Makela et al., 2008) (m2 mol-1) */
	K,                                   /* Extinction coefficient for absorption of PAR by canopy */
	ALBEDO,                              /* Canopy albedo */
	INT_COEFF,                           /* Precipitation interception coefficient */
	SLA_AVG0,                            /* AVERAGE Specific Leaf Area m2/Kg for sunlit/shaded leaves (juvenile) */
	SLA_AVG1,                            /* AVERAGE Specific Leaf Area m2/Kg for sunlit/shaded leaves (mature) */
	TSLA,                                /* Age at which SLA_AVG = (SLA_AVG1 + SLA_AVG0 )/2 */
	SLA_RATIO,                           /* (DIM) ratio of shaded to sunlit projected SLA */
	LAI_RATIO,                           /* (DIM) all-sided to projected leaf area ratio */
	FRACBB0,                             /* Branch and Bark fraction at age 0 (m^2/kg) */
	FRACBB1,                             /* Branch and Bark fraction for mature stands (m^2/kg) */
	TBB,                                 /* Age at which fracBB = (FRACBB0 + FRACBB1 )/ 2 */
	RHO0,                                /* Minimum Basic Density for young Trees (tDM/m3) */
	RHO1,                                /* Maximum Basic Density for mature Trees (tDM/m3) */
	TRHO,                                /* Age at which rho = (RHO0 + RHO1 )/2 */
	FORM_FACTOR,                         /* Stem form factor (adim) */
	COEFFCOND,                           /* Define stomatal responsee to VPD in m/sec */
	BLCOND,                              /* Canopy Boundary Layer conductance */
	MAXCOND,                             /* Maximum leaf Conductance in m/sec */
	CUTCOND,                             /* Cuticular conductance in m/sec */
	MAXAGE,                              /* Maximum tree age */
	RAGE,                                /* Relative Age to give fAGE = 0.5 */
	NAGE,                                /* Power of relative Age in function for Age */
	GROWTHTMIN,                          /* Minimum temperature for growth */
	GROWTHTMAX,                          /* Maximum temperature for growth */
	GROWTHTOPT,                          /* Optimum temperature for growth */
	GROWTHSTART,                         /* Thermic sum  value for starting growth in °C */
	TAU,                                 // day   
	X0,                                  // °C 
	Smax,                                // °C 
	MINDAYLENGTH,                        /* Minimum day length for phenology */
	SWPOPEN,                             /* (MPa) Soil water potential open */
	SWPCLOSE,                            /* (MPa) Soil water potential close */
	OMEGA_CTEM,                          /* ALLOCATION PARAMETER */
	S0CTEM,                              /* PARAMETER CONTROLLING ALLOCATION TO STEM */
	R0CTEM,                              /* PARAMETER CONTROLLING ALLOCATION TO ROOT */
	F0CTEM,                              /* PARAMETER CONTROLLING ALLOCATION TO FOLIAGE */
	FRUIT_PERC,                          /* percentage of npp to fruit */
	CONES_LIFE_SPAN,                     /* cone life span */
	FINE_ROOT_LEAF,                      /* allocation new fine root C:new leaf (ratio) */
	STEM_LEAF,                           /* allocation new stem C:new leaf (ratio) */
	COARSE_ROOT_STEM,                    /* allocation new coarse root C:new stem (ratio) */
	LIVE_TOTAL_WOOD,                     /* allocation new live wood C:new total wood C (ratio) */
	N_RUBISCO,                           /* Fraction of leaf N in Rubisco (ratio) */
	CN_LEAVES,                           /* CN of leaves (kgC/kgN) */
	CN_FALLING_LEAVES,                   /* CN of leaf litter (kgC/kgN) */
	CN_FINE_ROOTS,                       /* CN of fine roots (kgC/kgN) */
	CN_LIVEWOOD,                         /* CN of live woods (kgC/kgN) */
	CN_DEADWOOD,                         /* CN of dead woods (kgC/kgN) */
	LEAF_LITR_LAB_FRAC,                  /* (DIM) leaf litter labile fraction */
	LEAF_LITR_CEL_FRAC,                  /* (DIM) leaf litter cellulose fraction */
	LEAF_LITR_LIGN_FRAC,                 /* (DIM) leaf litter lignin fraction */
	FROOT_LITR_LAB_FRAC,                 /* (DIM) fine root litter labile fraction */
	FROOT_LITR_CEL_FRAC,                 /* (DIM) fine root litter cellulose fraction */
	FROOT_LITR_LIGN_FRAC,                /* (DIM) fine root litter lignin fraction */
	DEADWOOD_CEL_FRAC,                   /* (DIM) dead wood litter cellulose fraction */
	DEADWOOD_LIGN_FRAC,                  /* (DIM) dead wood litter lignin fraction */
	BUD_BURST,                           /* days of bud burst at the beginning of growing season (only for deciduous) */
	LEAF_FALL_FRAC_GROWING,              /* proportions of the growing season of leaf fall */
	LEAF_FROOT_TURNOVER,                 /* Average yearly leaves and fine root turnover rate */
	LIVEWOOD_TURNOVER,                   /* Average yearly livewood turnover rate */
	SAPWOOD_TURNOVER,                    /* Average yearly sapwood turnover rate */
	DBHDCMAX,                            /* maximum dbh crown diameter relationship when minimum density */
	DBHDCMIN,                            /* minimum dbh crown diameter relationship when maximum density */
	SAP_A,                               /* a coefficient for sapwood */
	SAP_B,                               /* b coefficient for sapwood */
	SAP_LEAF,                            /* (m2/m2) sapwood_max leaf area ratio in pipe model */
	SAP_WRES,                            /* Sapwood-Reserve biomass ratio used if no Wres data are available */
	STEMCONST_P,                         /* constant in the stem mass vs. diameter relationship */
	STEMPOWER_P,                         /* power in the stem mass vs. diameter relationship */
	CRA,                                 /* Chapman-Richards a parameter (maximum height, meter) */
	CRB,                                 /* Chapman-Richards b parameter */
	CRC,                                 /* Chapman-Richards c parameter */
	HDMAX_A,                             /* A parameter for Height (m) to Base diameter (m) ratio MAX */
	HDMAX_B,                             /* B parameter for Height (m) to Base diameter (m) ratio MAX */
	HDMIN_A,                             /* A parameter for Height (m) to Base diameter (m) ratio MIN */
	HDMIN_B,                             /* B parameter for Height (m) to Base diameter (m) ratio MIN */
	CROWN_FORM_FACTOR,                   /* Crown form factor (0 = cylinder, 1 = cone, 2 = sphere, 3 = ellipsoid) */
	CROWN_A,                             /* Crown a parameter */
	CROWN_B,                             /* Crown b parameter */
	MAXSEED,                             /* maximum seeds number (see TREEMIG) */
	MASTSEED,                            /* Musting year (see TREEMIG) */
	WEIGHTSEED,                          /* single fruit weight in g */
	SEXAGE,                              /* Age for sexual maturity */
	GERMCAPACITY,                        /* Geminability rate */
	ROTATION,                            /* rotation for final harvest (based on tree age) */
	THINNING,                            /* thinning regime (based on year simulation) */
	THINNING_REGIME,                     /* thinning regime (0 = above, 1 = below) */
	THINNING_INTENSITY,                  /* thinning intensity (% of Basal Area/N-tree to remove) */
    CONV_VOL_FACTOR,                     // conversion factor from stem carbon sto stem volume (GSV)

	MINSTOCKGROW,
	THINNING_80,
	THINNING_60_80,
	THINNING_40_60,
	THINNING_40_20,
	THINNING_0,

	/*********************************************STATE VARIABLES*********************************************/
	/*********************************************************************************************************/

	DBHDC_EFF,                           /* Crown Projected Diameter from DBH in function of density*/
	DENMAX,                              /* (trees/cell) maximum stand density  */
	DENMIN,                              /* (trees/cell) minimum stand density */
	MAX_LAYER_COVER,                     /* (fraction) maximum fraction of layer covered */
	CROWN_RADIUS,                        /* (m) Crown Projected Radius */
	CROWN_DIAMETER,                      /* (m) Crown Projected Diameter */
	CROWN_HEIGHT,                        /* (m) Crown height */
	CROWN_AREA,                          /* (m2) Crown Area */
	CROWN_AREA_PROJ,                     /* (m2) Crown Projected Area (at zenith angle) */
	//CROWN_AREA_EXP,                      /* (m2) Crown Exposed Area (considered as the area that can intercepts light) */
	CROWN_VOLUME,                        /* (m3) Crown Volume */
	CROWN_DENSITY,                       /* (m) Crown density */
	CANOPY_COVER_PROJ,                   /* (DIM) Canopy Projected Cover % of pixel covered (at zenith angle) */
	//CANOPY_COVER_EXP,                    /* (DIM) Canopy Exposed Cover % of pixel exposed (integrated all over all viewing angles) */
	DAILY_CANOPY_COVER_PROJ,             /* (DIM) Daily Canopy Projected Cover % of pixel covered (at zenith angle) */
	DAILY_CANOPY_COVER_EXP,              /* (DIM) Daily Canopy Exposed Cover % of pixel exposed (integrated all over the daylength) */
	HD_EFF,                              /* (DIM) Effective Height/Diameter ratio */
	HD_MAX,                              /* (DIM) Height (m) to Base diameter (m) ratio MAX */
	HD_MIN,                              /* (DIM) Height (m) to Base diameter (m) ratio MIN */
	SAPWOOD_AREA,                        /* (cm2) tree sapwood area */
	SAPWOOD_PERC,                        /* sapwood area to total stem area */
	HEARTWOOD_AREA,                      /* (cm2) tree heartwood area */
	HEARTWOOD_PERC,                      /* (DIM) heartwood area to total stem area */
	BASAL_AREA,                          /* (cm2/tree) Individual Basal Area */
	BASAL_AREA_m2,                       /* (m2/tree) Individual Basal Area */
	STAND_BASAL_AREA,                    /* (cm2/class) Class Basal Area of overall class */
	STAND_BASAL_AREA_m2,                 /* (m2/class) Class Basal Area of overall class */
	MASS_DENSITY,                        /* (tDM/m3) mass density */
	VOLUME,                              /* (m3/class) Stem Volume */
	TREE_VOLUME,                         /* (m3/tree) Single Tree-stem  Volume */
	CAI,                                 /* (m3/class/year) Current Annual Volume Increment */
	TREE_CAI,                            /* (m3/tree/year) Single Tree-stem  Current Annual Volume Increment */
	MAI,                                 /* (m3/class/year) Mean Annual Volume Increment */
	TREE_MAI,                            /* (m3/tree/year) Single Tree-stem  Mean Annual Volume Increment */

	VOLUME2,                              /* (m3/class) Stem Volume using conversion factor  */
	TREE_VOLUME2,                         /* (m3/tree) Single Tree-stem  Volume using conversion factor*/
	CAI2,                                 /* (m3/class/year) Current Annual Volume Increment using VOLUME2 */
	TREE_CAI2,                            /* (m3/tree/year) Single Tree-stem Current Annual Volume Increment using VOLUME2 */
	MAI2,                                 /* (m3/class/year) Mean Annual Volume Increment using VOLUME2 */
	TREE_MAI2,                            /* (m3/tree/year) Single Tree-stem  Mean Annual Volume Increment using VOLUME2 */
	
	STEMCONST,                           /* constant in the stem mass vs. diameter relationship */
	MAX_SPECIES_COVER,                   /* (DIM) maximum species-specific cell cover */
    
	/* PAR */
	PAR,                                /* (molPAR/m2/day) Photosynthetically Active Radiation  */
	PAR_SUN,                            /* (molPAR/m2/day) Photosynthetically Active Radiation for sun leaves */
	PAR_SHADE,                          /* (molPAR/m2/day) Photosynthetically Active Radiation for shaded leaves */
	PAR_REFL,                           /* (molPAR/m2/day) Reflected Photosynthetically Active Radiation OVERALL CANOPY */
	PAR_REFL_SUN,                       /* (molPAR/m2/day) Reflected Photosynthetically Active Radiation for sun leaves */
	PAR_REFL_SHADE,                     /* (molPAR/m2/day) Reflected Photosynthetically Active Radiation for shaded leaves */
	APAR,                               /* (molPAR/m2/day) Absorbed Photosynthetically Active Radiation */
	APAR_SUN,                           /* (molPAR/m2/day) Absorbed Photosynthetically Active Radiation for sun leaves */
	APAR_SHADE,                         /* (molPAR/m2/day) Absorbed Photosynthetically Active Radiation for shaded leaves */
	TRANSM_PAR,                         /* (molPAR/m2/day) Transmitted Photosynthetically Active Radiation */
	TRANSM_PAR_SUN,                     /* (molPAR/m2/day) Transmitted Photosynthetically Active Radiation from sun leaves */
	TRANSM_PAR_SHADE,                   /* (molPAR/m2/day) Transmitted Photosynthetically Active Radiation from shaded leaves */
	fAPAR,                              /* (unitless) Fraction of absorbed Photosynthetically Active Radiation */

	YEARLY_APAR,

	/* short wave */
	SW_RAD,                             /* (W/m2) Short Wave Radiation  */
	SW_RAD_SUN,                         /* (W/m2) Short Wave Radiation for sun leaves */
	SW_RAD_SHADE,                       /* (W/m2) Short Wave Radiation for shaded leaves */
	SW_RAD_REFL,                        /* (W/m2) Reflected Short Wave Radiation OVERALL CANOPY */
	SW_RAD_REFL_SUN,                    /* (W/m2) Reflected Short Wave Radiation for sun leaves */
	SW_RAD_REFL_SHADE,                  /* (W/m2) Reflected Short Wave Radiation for shaded leaves */
	SW_RAD_ABS,                         /* (W/m2) Absorbed Net Short Wave radiation */
	SW_RAD_ABS_SUN,                     /* (W/m2) Absorbed Net Short Wave Radiation for sun leaves */
	SW_RAD_ABS_SHADE,                   /* (W/m2) Absorbed Net Short Wave Radiation for shaded leaves */
	SW_RAD_TRANSM,                      /* (W/m2) Transmitted Net Short Wave Radiation */
	SW_RAD_TRANSM_SUN,                  /* (W/m2) Transmitted Net Short Wave Radiation for sun leaves */
	SW_RAD_TRANSM_SHADE,                /* (W/m2) Transmitted Net Short Wave Radiation for shaded leaves */

	/* long wave */
	NET_LW_RAD,                         /* (W/m2) Long Wave Radiation */
	LW_RAD_REFL,                        /* (W/m2) Reflected Long Wave Radiation OVERALL CANOPY */
	LW_RAD_EMIT,                        /* (W/m2) Emitted Long Wave radiation */
	LW_RAD_EMIT_SUN,                    /* (W/m2) Emitted Long Wave Radiation for sun leaves */
	LW_RAD_EMIT_SHADE,                  /* (W/m2) Emitted Long Wave Radiation for shaded leaves */
	LW_RAD_ABS,                         /* (W/m2) Absorbed Long Wave radiation */
	LW_RAD_ABS_SUN,                     /* (W/m2) Absorbed Long Wave Radiation for sun leaves */
	LW_RAD_ABS_SHADE,                   /* (W/m2) Absorbed Long Wave Radiation for shaded leaves */
	LW_RAD_TRANSM,                      /* (W/m2) Transmitted Long Wave Radiation */
	LW_RAD_TRANSM_SUN,                  /* (W/m2) Transmitted Long Wave Radiation for sun leaves */
	LW_RAD_TRANSM_SHADE,                /* (W/m2) Transmitted Long Wave Radiation for shaded leaves */

	/* net radiation */
	NET_RAD,                            /* (W/m2) Short Wave Radiation */
	NET_RAD_SUN,
	NET_RAD_SHADE,
	NET_RAD_REFL,                       /* (W/m2) Reflected Radiation OVERALL CANOPY */
	NET_RAD_ABS,                        /* (W/m2) Absorbed Net radiation */
	NET_RAD_ABS_SUN,                    /* (W/m2) Absorbed Net Radiation for sun leaves */
	NET_RAD_ABS_SHADE,                  /* (W/m2) Absorbed Net Radiation for shaded leaves */
	NET_RAD_TRANSM,                     /* (W/m2) Transmitted Net Radiation */
	NET_RAD_TRANSM_SUN,                 /* (W/m2) Transmitted Net Radiation for sun leaves */
	NET_RAD_TRANSM_SHADE,               /* (W/m2) Transmitted Net Radiation for shaded leaves */
	NET_RAD_REFL_SUN,
	NET_RAD_REFL_SHADE,

	/* modifiers */
	F_VPD,                              /* (DIM) VPD modifier */
	F_LIGHT,                            /* (DIM) LIGHT modifier */
	F_LIGHT_SUN,                        /* (DIM) LIGHT modifier for Sun leaves */
	F_LIGHT_SHADE,                      /* (DIM) LIGHT modifier for Shaded leaves */
	F_AGE,                              /* (DIM) AGE modifier */
	FR,                                 /* (DIM) Nitrogen fertility rate */
	F_NUTR,                             /* (DIM) SOIL NUTRIENT Modifier */
	F_T,                                /* (DIM) TEMPERATURE modifier */
	F_SW,                               /* (DIM) SOIL WATER modifier */
	F_DROUGHT,                          /* (DIM) SOIL DROUGHT modifier (see Duursma et al., 2008) */
	F_PSI,                              /* (DIM) SOIL WATER modifier using PSI, see Biome */
	F_CO2,                              /* (DIM) CO2 fertilization effect (as choiced in script) */
	F_CO2_VER,                          /* (DIM) CO2 fertilization effect (Veroustraete's version) */
	F_CO2_WANG,                         /* (DIM) CO2 fertilization effect (Wang et al's version) */
	F_CO2_FRANKS,                       /* (DIM) CO2 fertilization effect (Franks et al's version) */
	F_CO2_TR,                           /* (DIM) CO2 reduction effect (for stomatal conductance) */
	PHYS_MOD,                           /* (DIM) PHYSIOLOGICAL modifier */
	F_LIGHT_VERT,

        old_X_acl,                          //
        F_ACCL,                             // cold-acclimation factor for assimilation
                
        m_corr,                             // correction factor for assimilation,when stomatal conductance is close to 0
	m_corr_sun,
	m_corr_shade,

	F_LIGHT_MAKELA,                            /* (DIM) LIGHT modifier */
	F_LIGHT_SUN_MAKELA,                        /* (DIM) LIGHT modifier for Sun leaves */
	F_LIGHT_SHADE_MAKELA,                      /* (DIM) LIGHT modifier for Shaded leaves */

	/* water */
	/* pools */
	CANOPY_WATER,                       /* (mm/m2) canopy water pool */
	CANOPY_SNOW,                        /* (mm/m2) canopy snow pool  */

	/* fluxes */
	STOMATAL_CONDUCTANCE,               /* (m/sec) stomatal conductance */
	STOMATAL_SUN_CONDUCTANCE,           /* (m/sec) stomatal conductance for sun leaves */
	STOMATAL_SHADE_CONDUCTANCE,         /* (m/sec) stomatal conductance for shade leaves */
	LEAF_CONDUCTANCE,                   /* (m/sec) leaf conductance */
	LEAF_SUN_CONDUCTANCE,               /* (m/sec) leaf conductance for sun leaves */
	LEAF_SHADE_CONDUCTANCE,             /* (m/sec) leaf conductance for shade leaves */
	CANOPY_BLCOND,                      /* (m/sec) Canopy Boundary Layer conductance */
	FRAC_DAYTIME_TRANSP,                /* (DIM) fraction of daytime for transpiration (and photosynthesis) */
	CANOPY_INT_RAIN,                    /* (mm/day) canopy interception of rainfall */
	CANOPY_INT_SNOW,                    /* (mm/day) canopy interception of snow */
	CANOPY_EVAPO,                       /* (mm/day) canopy evaporation */
	CANOPY_SUBLI,                       /* (mm/day) canopy sublimation */
	CANOPY_MELT,                        /* (mm/day) canopy melting */
	CANOPY_TRANSP,                      /* (mm/day) canopy transpiration */
	CANOPY_TRANSP_SUN,                  /* (mm/day) canopy transpiration for sun leaves */
	CANOPY_TRANSP_SHADE,                /* (mm/day) canopy transpiration for shade leaves */
	CANOPY_EVAPO_TRANSP,                /* (mm/day) canopy evapotranspiration */
	CANOPY_FRAC_DAY_TRANSP,             /* (DIM) Fraction of daily canopy transpiration */
	CANOPY_LATENT_HEAT,                 /* (W/m2) canopy latent heat */
	CANOPY_SENSIBLE_HEAT,               /* (W/m2) canopy sensible heat */

	MONTHLY_CANOPY_TRANSP,              /* (mm/month) monthly canopy transpiration */
	MONTHLY_CANOPY_TRANSP_SUN,          /* (mm/month) monthly canopy transpiration for sun leaves */
	MONTHLY_CANOPY_TRANSP_SHADE,        /* (mm/month) monthly canopy transpiration for shade leaves */
	MONTHLY_CANOPY_INT,                 /* (mm/month) monthly canopy interception */
	MONTHLY_CANOPY_EVAPO_TRANSP,        /* (mm/month) monthly canopy evapotranspiration */
	YEARLY_CANOPY_TRANSP,               /* (mm/year) annual canopy transpiration */
	YEARLY_CANOPY_TRANSP_SUN,           /* (mm/year) annual canopy transpiration for sun leaves */
	YEARLY_CANOPY_TRANSP_SHADE,         /* (mm/year) annual canopy transpiration for shade leaves */
	YEARLY_CANOPY_INT,                  /* (mm/year) annual canopy interception */
	YEARLY_CANOPY_EVAPO_TRANSP,         /* (mm/year) annual canopy evapotranspiration */
	MONTHLY_CANOPY_LATENT_HEAT,         /* (W/m2) monthly canopy latent heat */
	YEARLY_CANOPY_LATENT_HEAT,          /* (W/m2) annual canopy latent heat */

	/* carbon use efficiency */
	CUE,                                /* (gCNPP/gCGPP) daily Carbon Use Efficiency */
	MONTHLY_CUE,                        /* (gCNPP/gCGPP) monthly Carbon Use Efficiency */
	YEARLY_CUE,                         /* (gCNPP/gCGPP) annual Carbon Use Efficiency */
	/* biomass production efficiency */
	BPE,                                /* (gCBP/gCGPP) biomass production Efficiency */
	MONTHLY_BPE,                        /* (gCBP/gCGPP) monthly biomass production Efficiency */
	YEARLY_BPE,                         /* (gCBP/gCGPP) annual biomass production Efficiency */

	/* water use efficiency */
	WUE,                                /* daily Water Use Efficiency */
	MONTHLY_WUE,                        /* monthly Water Use Efficiency */
	YEARLY_WUE,                         /* annual Water Use Efficiency */
	iWUE,                               /* daily intrinsic Water Use Efficiency */
	MONTHLY_iWUE,                       /* monthly intrinsic Water Use Efficiency */
	YEARLY_iWUE,                        /* annual intrinsic Water Use Efficiency */
	iWUE_SUN,                           /* daily intrinsic Water Use Efficiency for sun leaves */
	MONTHLY_iWUE_SUN,                   /* monthly intrinsic Water Use Efficiency for sun leaves */
	YEARLY_iWUE_SUN,                    /* annual intrinsic Water Use Efficiency for sun leaves */
	iWUE_SHADE,                         /* daily intrinsic Water Use Efficiency for shade leaves */
	MONTHLY_iWUE_SHADE,                 /* monthly intrinsic Water Use Efficiency for shade leaves */
	YEARLY_iWUE_SHADE,                  /* annual intrinsic Water Use Efficiency for shade leaves */
	gsWUE,                              /* daily intrinsic Water Use Efficiency (stomatal) */
	gsWUE_SUN,                          /* daily intrinsic Water Use Efficiency (stomatal) for sun leaves */
	gsWUE_SHADE,                        /* daily intrinsic Water Use Efficiency (stomatal) for shade leaves */

	/* LAI */
	SLA_AVG,                            /* (m2/KgC) Age-related Average Specific Leaf Area */
	SLA_PROJ,                           /* (m2/KgC) Current Age-related Average Specific Leaf Area */
	SLA_SUN_PROJ,                       /* (m2/KgC) Current Age-related Porjected Average Specific Leaf Area for sun leaves*/
	SLA_SHADE_PROJ,                     /* (m2/KgC) Current Age-related Porjected Average Specific Leaf Area for shaded leaves*/
	LEAF_AREA,                          /* (m2) Leaf area */
	LEAF_AREA_SUN,                      /* (m2) Leaf Area for sun leaves */
	LEAF_AREA_SHADE,                    /* (m2) Leaf Area for shaded leaves */
	LAI_PROJ,                           /* (m2/m2) LAI for Projected Area covered (at zenith angle) */
	LAI_SUN_PROJ,                       /* (m2/m2) LAI for sun leaves for Projected Area covered (at zenith angle) */
	LAI_SHADE_PROJ,                     /* (m2/m2) LAI for shaded leaves for Projected Area covered (at zenith angle) */
	TREE_LAI_PROJ,                           /* (m2/m2) LAI for Projected Area covered (at zenith angle) */
	TREE_LAI_SUN_PROJ,                       /* (m2/m2) LAI for sun leaves for Projected Area covered (at zenith angle) */
	TREE_LAI_SHADE_PROJ,                     /* (m2/m2) LAI for shaded leaves for Projected Area covered (at zenith angle) */
	MAX_LAI_PROJ,                       /* (m2/m2) Maximum annual LAI */
	PEAK_LAI_PROJ,                      /* (m2/m2) Peak Projected LAI */
	ALL_LAI_PROJ,                       /* (m2) All side Projected LAI */
	LAI_EXP,                            /* (m2/m2) LAI for Exposed Area covered (all over all viewing angles) */
	LAI_SUN_EXP,                        /* (m2/m2) LAI for sun leaves for Exposed Area covered (all over all viewing angles) */
	LAI_SHADE_EXP,                      /* (m2/m2) LAI for shaded leaves for Exposed Area covered (all over all viewing angles) */
	MAX_LAI_EXP,                        /* (m2/m2) Maximum annual LAI exposed */
	PEAK_LAI_EXP,                       /* (m2/m2) Peak Projected LAI exposed */
	MAX_LAI_LEAFFALL_PROJ,              /* (m2/m2) LAI Projected at the very first day of senescence: parameter of the sigmoid function */
	LEAF_TEMP_K,                        /* (K) leaf temperature */
	CANOPY_TEMP_K,                      /* (K) canopy temperature */
	CANOPY_TEMP_K_OLD,                  /* (K) previous time canopy temperature */
    MONTHLY_LAI_PROJ,                  /* (m2/m2) monthly mean LAI */
	
	/* assimilation Farqhuar */
	A_SUN,
	Av_SUN,
	Aj_SUN,
	YEARLY_A_SUN,
	YEARLY_Av_SUN,
	YEARLY_Aj_SUN,
	A_SHADE,
	Av_SHADE,
	Aj_SHADE,
	YEARLY_A_SHADE,
	YEARLY_Av_SHADE,
	YEARLY_Aj_SHADE,
	A_TOT,
	Av_TOT,
	Aj_TOT,
	YEARLY_A_TOT,
	YEARLY_Av_TOT,
	YEARLY_Aj_TOT,


	ASSIMILATION,                       /* (gC/m2/day) Daily Gross Assimilation (carboxylation) (which include dark respiration */
	ASSIMILATION_SUN,                   /* (gC/m2/day) Daily Gross Assimilation (carboxylation) for sun leaves (which include dark respiration */
	ASSIMILATION_SHADE,                 /* (gC/m2/day) Daily Gross Assimilation (carboxylation) for shaded leaves (which include dark respiration */
	MONTHLY_ASSIMILATION,               /* (gC/m2/month) Daily Gross Assimilation (carboxylation) (which include dark respiration */
	MONTHLY_ASSIMILATION_SUN,           /* (gC/m2/month) Daily Gross Assimilation (carboxylation) for sun leaves (which include dark respiration */
	MONTHLY_ASSIMILATION_SHADE,         /* (gC/m2/month) Daily Gross Assimilation (carboxylation) for shaded leaves (which include dark respiration */
	YEARLY_ASSIMILATION,                /* (gC/m2/year) Daily Gross Assimilation (carboxylation) (which include dark respiration */
	YEARLY_ASSIMILATION_SUN,            /* (gC/m2/year) Daily Gross Assimilation (carboxylation) for sun leaves (which include dark respiration */
	YEARLY_ASSIMILATION_SHADE,          /* (gC/m2/year) Daily Gross Assimilation (carboxylation) for shaded leaves (which include dark respiration */

	/* carbon fluxes */
	C_FLUX,
	ALPHA_EFF,                          /* (molC/molPAR) Effective quantum canopy efficiency */
	ALPHA_EFF_SUN,                      /* (molC/molPAR) Effective quantum canopy efficiency */
	ALPHA_EFF_SHADE,                    /* (molC/molPAR) Effective quantum canopy efficiency */
	GPP,                                /* (gC/m2/day) Daily Gross Primary Production */
	GPP_SUN,                            /* (gC/m2/day) Daily Sun Leaves Gross Primary Production */
	GPP_SHADE,                          /* (gC/m2/day) Daily Shade Leaves Gross Primary Production */
	MONTHLY_GPP,                        /* (gC/m2/month) Monthly Gross Primary Production */
	MONTHLY_GPP_SUN,                    /* (gC/m2/month) Monthly Gross Primary Production for sun leaves */
	MONTHLY_GPP_SHADE,                  /* (gC/m2/month) Monthly Gross Primary Production for shaded leaves */
	YEARLY_GPP,                         /* (gC/m2/year) Yearly Gross Primary Production */
	YEARLY_GPP_SUN,                     /* (gC/m2/year) Yearly Gross Primary Production for sun leaves */
	YEARLY_GPP_SHADE,                   /* (gC/m2/year) Yearly Gross Primary Production for shaded leaves */
	GPP_tC,                             /* (tC/cell/day) Daily Gross Primary Production  */
	NPP,                                /* (gC/m2/day) Daily Net Primary Production */
	MONTHLY_NPP,                        /* (gC/m2/monthly) Monthly Net Primary Production */
	YEARLY_NPP,                         /* (gC/m2/year) Yearly Net Primary Production */
	NPP_tC,                             /* (tC/cell/day) Daily Net Primary Production */
	MONTHLY_NPP_tC,                     /* (tC/cell/month) Monthly Net Primary Production */
	YEARLY_NPP_tC,                      /* (tC/cell/day) Yearly Net Primary Production */
	BP,                                 /* (gC/m2/day) Daily Biomass Production */
	MONTHLY_BP,                         /* (gC/m2/monthly) Monthly Biomass Production */
	YEARLY_BP,                          /* (gC/m2/year) Yearly Biomass Production */

	/* nitrogen fluxes */
	NPP_gN_DEMAND,                      /* (gN/m2/day) Daily Net Primary Production demand */
	MONTHLY_NPP_gN,                     /* (gN/m2/month) Monthly Net Primary Production demand */
	YEARLY_NPP_gN,                      /* (gN/m2/year) Yearly Net Primary Production demand */
	NPP_tN_DEMAND,                      /* (tNC/cell/day) Daily Net Primary Production demand */
	MONTHLY_NPP_tN,                     /* (tNC/cell/month) Monthly Net Primary Production demand */
	YEARLY_NPP_tN,                      /* (tNC/cell/year) Yearly Net Primary Production demand */

	/* transfer pools carbon to carbon pools */
	C_TO_LEAF,                          /* (tC/cell/day) Daily Net Primary Production to Leaf pool */
	C_TO_FROOT,                         /* (tC/cell/day) Daily Net Primary Production to Fine root pool */
	C_TO_CROOT,                         /* (tC/cell/day) Daily Net Primary Production to Coarse root pool */
	C_TO_STEM,                          /* (tC/cell/day) Daily Net Primary Production to Stem pool */
	C_TO_BRANCH,                        /* (tC/cell/day) Daily Net Primary Production to Branch pool */
	C_TO_RESERVE,                       /* (tC/cell/day) Daily Net Primary Production to Reserve pool */
	C_TO_FRUIT,                         /* (tC/cell/day) Daily Net Primary Production to Fruit pool */
	LEAF_C_TO_REMOVE,                   /* (tC/cell/day) Daily Leaf to remove */
	FROOT_C_TO_REMOVE,                  /* (tC/cell/day) Daily Fine Root to remove tC/sizeCell */
	FRUIT_C_TO_REMOVE,                  /* (tC/cell/day) Daily Fruit to remove */
	BRANCH_C_TO_REMOVE,                 /* (tC/cell/day) Annual Branch to remove */
	CROOT_C_TO_REMOVE,                  /* (tC/cell/day) Annual Coarse root to remove */
	BRANCH_SAP_C_TO_REMOVE,                 /* (tC/cell/day) Annual Branch to remove */
	CROOT_SAP_C_TO_REMOVE,                  /* (tC/cell/day) Annual Coarse root to remove */

	LEAF_FALL_C,                        //(tC/cell/day) Daily Leaf to remove -evergreen - constant for each day 
	FROOT_FALL_C,                       //(tC/cell/day) Daily froot to remove -evergreen - constant for each day
	FRUIT_FALL_C,                       //(tC/cell/day) Daily fruit to remove -evergreen - constant for each day

	C_LEAF_TO_RESERVE,                  /* (tC/cell/day) Daily retranslocated C from Leaf pool to Reserve */
	C_FROOT_TO_RESERVE,                 /* (tC/cell/day) Daily retranslocated C from Fine root pool to Reserve */
	C_LEAF_TO_LITR,                     /* (tC/cell/day) Daily tC from Leaf pool to Litter */
	C_FROOT_TO_LITR,                    /* (tC/cell/day) Daily tC from Fine root pool to Soil */
	C_RESERVE_TO_LITR,                   /* (tC/cell/day) Daily tC from Reserve pool to litter pool */
	C_BRANCH_TO_RESERVE,                /* (tC/cell/day) Daily tC from Branch pool to Litter */
	C_CROOT_TO_RESERVE,                 /* (tC/cell/day) Daily tC from Coarse root pool to Litter */
	C_STEM_TO_CWD,                      /* (tC/cell/day) Daily tC from Stem pool to Coarse Woody Debris */
	C_CROOT_TO_CWD,                     /* (tC/cell/day) Daily tC from Coarse root pool to Coarse Woody Debris */
	C_BRANCH_TO_CWD,                    /* (tC/cell/day) Daily tC from Branch pool to Coarse Woody Debris */
	//C_RESERVE_TO_CWD,                   /* (tC/cell/day) Daily tC from Reserve pool to Coarse Woody Debris */
	C_FRUIT_TO_CWD,                     /* (tC/cell/day) Daily tC from Fruit pool to Coarse Woody Debris */
	C_STEM_LIVEWOOD_TO_CWD,             /* (tC/cell/day) Daily tC from Stem livewood to Coarse Woody Debris */
	C_CROOT_LIVEWOOD_TO_CWD,            /* (tC/cell/day) Daily tC from Coarse root livewood to Coarse Woody Debris */
	C_BRANCH_LIVEWOOD_TO_CWD,           /* (tC/cell/day) Daily tC from Branch livewood to Coarse Woody Debris */
	C_STEM_DEADWOOD_TO_CWD,             /* (tC/cell/day) Daily tC from Stem deadwood to Coarse Woody Debris */
	C_CROOT_DEADWOOD_TO_CWD,            /* (tC/cell/day) Daily tC from Coarse root deadwood to Coarse Woody Debris */
	C_BRANCH_DEADWOOD_TO_CWD,           /* (tC/cell/day) Daily tC from Branch deadwood to Coarse Woody Debris */
	C_STEM_SAPWOOD_TO_CWD,              /* (tC/cell/day) Daily tC from Stem Sapwood to Coarse Woody Debris */
	C_CROOT_SAPWOOD_TO_CWD,             /* (tC/cell/day) Daily tC from Coarse root Sapwood to Coarse Woody Debris */
	C_BRANCH_SAPWOOD_TO_CWD,            /* (tC/cell/day) Daily tC from Branch Sapwood to Coarse Woody Debris */
	C_STEM_HEARTWOOD_TO_CWD,            /* (tC/cell/day) Daily tC from Stem Heartwood to Coarse Woody Debris */
	C_CROOT_HEARTWOOD_TO_CWD,           /* (tC/cell/day) Daily tC from Coarse root Heartwood to Coarse Woody Debris */
	C_BRANCH_HEARTWOOD_TO_CWD,          /* (tC/cell/day) Daily tC from Branch Heartwood to Coarse Woody Debris */
	C_TO_LITR,                          /* (tC/cell/day) Daily Leaf and Fine root Carbon to Litter pool */
	C_TO_CWD,                           /* (tC/cell/day) Daily Wood Carbon to Coarse Woody Debirs */
	C_STEM_LIVEWOOD_TO_DEADWOOD,        /* (tC/cell/day) Daily tC from Stem live wood pool to Stem dead wood */
	C_CROOT_LIVEWOOD_TO_DEADWOOD,       /* (tC/cell/day) Daily tC from Coarse live wood pool to Coarse dead wood */
	C_BRANCH_LIVEWOOD_TO_DEADWOOD,      /* (tC/cell/day) Daily tC from Branch live wood pool to Branch dead wood */
	C_STEM_SAPWOOD_TO_HEARTWOOD,        /* (tC/cell/day) Daily tC from Stem sapwood pool to Stem heartwood */
	C_CROOT_SAPWOOD_TO_HEARTWOOD,       /* (tC/cell/day) Daily tC from Coarse sapwood pool to Coarse heartwood */
	C_BRANCH_SAPWOOD_TO_HEARTWOOD,      /* (tC/cell/day) Daily tC from Branch sapwood pool to Branch heartwood */
	C_RESERVE_TO_LEAF_BUDBURST,         /* (tC/cell/day) Daily tC from Reserve pool to Leaf Budburst */
	C_RESERVE_TO_FROOT_BUDBURST,        /* (tC/cell/day) Daily tC from Reserve pool to Fine root Budburst */
	C_RESERVE_TO_BUDBURST,              /* (tC/cell/day) Daily tC from Reserve pool to Budburst */
	TOT_C_RESERVE_TO_LEAF_BUDBURST,     /* (tC/cell/day) Total tC from Reserve pool to Leaf root Budburst */
	TOT_C_RESERVE_TO_FROOT_BUDBURST,    /* (tC/cell/day) Total tC from Reserve pool to Fine root Budburst */

	YEARLY_RESERVE_ALLOC,
	YEARLY_RESERVE_DEPLE,
	YEARLY_RESERVE_USAGE,

	YEARLY_C_TO_STEM,                   /* (tC/cell/year) Annual Stem biomass for turnover */
	YEARLY_C_TO_CROOT,                  /* (tC/cell/year) Annual Coarse root biomass for turnover */
	YEARLY_C_TO_BRANCH,                 /* (tC/cell/year) Annual Branch biomass for turnover */

	YEARLY_C_TO_STEM_LIVEWOOD,           /* (tC/cell/year) Annual Stem biomass for turnover */
	YEARLY_C_TO_CROOT_LIVEWOOD,          /* (tC/cell/year) Annual Coarse root biomass for turnover */
	YEARLY_C_TO_BRANCH_LIVEWOOD,         /* (tC/cell/year) Annual Branch biomass for turnover */

	YEARLY_STEM_LIVEWOOD_TO_DEADWOOD,   /* (tC/cell/year) Annual Stem biomass for turnover */
	YEARLY_CROOT_LIVEWOOD_TO_DEADWOOD,  /* (tC/cell/year) Annual Coarse root biomass for turnover */
	YEARLY_BRANCH_LIVEWOOD_TO_DEADWOOD, /* (tC/cell/year) Annual Branch biomass for turnover */

	YEARLY_C_TO_STEM_SAPWOOD,           /* (tC/cell/year) Annual Stem biomass for turnover */
	YEARLY_C_TO_CROOT_SAPWOOD,          /* (tC/cell/year) Annual Coarse root biomass for turnover */
	YEARLY_C_TO_BRANCH_SAPWOOD,         /* (tC/cell/year) Annual Branch biomass for turnover */

	YEARLY_STEM_SAPWOOD_TO_HEARTWOOD,   /* (tC/cell/year) Annual Stem biomass for turnover */
	YEARLY_CROOT_SAPWOOD_TO_HEARTWOOD,  /* (tC/cell/year) Annual Coarse root biomass for turnover */
	YEARLY_BRANCH_SAPWOOD_TO_HEARTWOOD, /* (tC/cell/year) Annual Branch biomass for turnover */

	YEARLY_C_TO_WOOD,                   /* (tC/cell/year) Annual Carbon stocked into wood pool */
	CUM_YEARLY_C_TO_WOOD,               /* (tC/cell) Cumulated Annual Carbon stocked into wood pool */

	/* carbon pools in tons of tC/sizeCell */
	LEAF_C,                             /* (tC/cell) Current Leaf carbon pool */
	LEAF_SUN_C,                         /* (tC/cell) Current Leaf sun carbon pool */
	LEAF_SHADE_C,                       /* (tC/cell) Current Leaf shade carbon pool */
	MAX_LEAF_C,                         /* (tC/cell/year) Maximum Leaf carbon pool */
	CROOT_C,                            /* (tC/cell) Current Coarse root carbon pool */
	FROOT_C,                            /* (tC/cell) Current Fine root carbon pool */
	MAX_FROOT_C,                        /* (tC/cell/year) Maximum Fine root carbon pool */
	MAX_BUD_BURST_C,                    /* (tC/cell) Maximum Current Leaf carbon pool for BudBurst */
	STEM_C,                             /* (tC/cell) Current Stem carbon pool */
	BRANCH_C,                           /* (tC/cell) Current Branch carbon pool */
	MIN_RESERVE_C,                      /* (tC/cell) Current Minimum reserve carbon pool */  // corr:  minimum reserve pool which has to be guaranteed that year.
        MIN_ANN_RESERVE_C,                  /* (tC/cell) Current Minimum reserve carbon pool */  // corr:  minimum reserve value for that year
	MAX_ANN_RESERVE_C,                  /* (tC/cell) Current Minimum reserve carbon pool */  // corr:  maximum reserve value for that year
	MAX_RESERVE_C_CONC,                 /* (%) Current Maximum reserve concentration carbon pool */  // max for the year
	MIN_RESERVE_C_CONC,                 /* (%) Current Minimum reserve concentration carbon pool */  // min for the year
	RESERVE_C,                          /* (tC/cell) Current Reserve carbon pool */
	FRUIT_C,                            /* (tC/cell) Current Fruit carbon pool */
	MAX_FRUIT_C,                        /* (tC/cell/year) Annual Fruit carbon pool */
	LITR_C,                             /* (tC/cell) Current Litter carbon pool */
	CWD_C,                              /* (tC/cell) Current Coarse Woody Debris carbon pool */
	TOT_LIVEWOOD_C,                     /* (tC/cell) Current Live wood carbon pool */
	TOT_DEADWOOD_C,                     /* (tC/cell) Current Dead wood carbon pool */
	TOT_WOOD_C,                         /* (tC/cell) Current Total wood carbon pool */
	TOTAL_C,                            /* (tC/cell) Current Total carbon pool */
	STEM_LIVEWOOD_C,                    /* (tC/cell) Current Stem live wood carbon pool */
	STEM_DEADWOOD_C,                    /* (tC/cell) Current Stem dead wood carbon pool */
	CROOT_LIVEWOOD_C,                   /* (tC/cell) Current Coarse root live wood carbon pool */
	CROOT_DEADWOOD_C,                   /* (tC/cell) Current Coarse root dead wood carbon pool */
	BRANCH_LIVEWOOD_C,                  /* (tC/cell) Current Branch live wood carbon pool */
	BRANCH_DEADWOOD_C,                  /* (tC/cell) Current Branch dead wood carbon pool */
	TOT_SAPWOOD_C,                      /* (tC/cell) Current total sapwood carbon pool */
	TOT_HEARTWOOD_C,                    /* (tC/cell) Current total heartwood carbon pool */
	STEM_SAPWOOD_C,                     /* (tC/cell) Current Stem sapwood carbon pool */
	STEM_HEARTWOOD_C,                   /* (tC/cell) Current Stem heartwood carbon pool */
	CROOT_SAPWOOD_C,                    /* (tC/cell) Current Coarse root sapwood carbon pool */
	CROOT_HEARTWOOD_C,                  /* (tC/cell) Current Coarse root heartwood carbon pool */
	BRANCH_SAPWOOD_C,                   /* (tC/cell) Current Branch sapwood carbon pool */
	BRANCH_HEARTWOOD_C,                 /* (tC/cell) Current Branch heartwood carbon pool */
	EFF_LIVE_TOTAL_WOOD_FRAC,           /* (DIM) Age-related fraction of Live biomass per Total biomass */
	DAILY_LIVEWOOD_TURNOVER,            /* (tC/cell/day) Daily live wood turnover rate */
	DAILY_SAPWOOD_TURNOVER,             /* (tC/cell/day) Daily sapwood turnover rate */
	STANDING_WOOD,                      /* (tC/cell) standing wood carbon */

	/* management C */
	C_HWP,                              /* (tC/cell/yr) annual harvested woody products removed from stand */
	CUM_C_HWP,                          /* (tC/cell) cumulated harvested woody products removed from stand */
	VOLUME_HWP,                         /* (m3/cell/yr) annual volume harvested woody products removed from stand */
	CUM_VOLUME_HWP,                     /* (m3/cell) annual volume harvested woody products removed from stand */
    
    VOLUME2_HWP,                         /* (m3/cell/yr) annual volume harvested woody products removed from stand using VOLUME2*/
	CUM_VOLUME2_HWP,                     /* (m3/cell) annual volume harvested woody products removed from stand using VOLUME2*/
    

	C_STEM_TO_HWP,                      /* (tC/cell/yr) annual harvested woody products removed from stand as stem */
	C_BRANCH_TO_HWP,                    /* (tC/cell/yr) annual harvested woody products removed from stand as branch*/

	OLD_RESERVE_C,    /* TO REMOVE */   /* (tC/cell) Old reserve carbon pool */

	/* per tree in tC */
	TREE_LEAF_C,                        /* (tC/tree) Single Tree Leaf carbon pool */
	TREE_STEM_C,                        /* (tC/tree) Single Tree Stem carbon pool */
	TREE_FROOT_C,                       /* (tC/tree) Single Tree Fine root carbon pool */
	TREE_CROOT_C,                       /* (tC/tree) Single Tree Coarse carbon pool */
	TREE_RESERVE_C,                     /* (tC/tree) Single Tree Reserve carbon pool */
	TREE_MIN_RESERVE_C,                 /* (tC/tree) Single Tree Minimum Reserve carbon pool */
        TREE_MIN_ANN_RESERVE_C,             /* (tC/tree) Minimum reserve carbon pool */  // ddalmo corr:  minimum reserve value for that year
	TREE_MAX_ANN_RESERVE_C,             /* (tC/tree) Minimum reserve carbon pool */  // ddalmo corr:  maximum reserve value for that year
	TREE_FRUIT_C,                       /* (tC/tree) Single Tree Fruit carbon pool */
	TREE_BRANCH_C,                      /* (tC/tree) Single Tree Branch carbon pool */
	TREE_STEM_SAPWOOD_C,                /* (tC/tree) Single Tree Sapwood Stem carbon pool */
	TREE_STEM_HEARTWOOD_C,              /* (tC/tree) Single Tree Heartwood Stem carbon pool */
	TREE_STEM_LIVEWOOD_C,               /* (tC/tree) Single Tree Live Stem carbon pool */
	TREE_STEM_DEADWOOD_C,               /* (tC/tree) Single Tree Dead Stem carbon pool */
	TREE_CROOT_SAPWOOD_C,               /* (tC/tree) Single Tree Sapwood Coarse root carbon pool */
	TREE_CROOT_HEARTWOOD_C,             /* (tC/tree) Single Tree Heartwood Coarse root carbon pool */
	TREE_CROOT_LIVEWOOD_C,              /* (tC/tree) Single Tree Live Coarse root carbon pool */
	TREE_CROOT_DEADWOOD_C,              /* (tC/tree) Single Tree Dead Coarse root carbon pool */
	TREE_BRANCH_SAPWOOD_C,              /* (tC/tree) Single Tree Sapwood Branch carbon pool */
	TREE_BRANCH_HEARTWOOD_C,            /* (tC/tree) Single Tree Heartwood Branch carbon pool */
	TREE_BRANCH_LIVEWOOD_C,             /* (tC/tree) Single Tree Live Branch carbon pool */
	TREE_BRANCH_DEADWOOD_C,             /* (tC/tree) Single Tree Dead Stem carbon pool */
	TREE_TOT_LIVEWOOD_C,                /* (tC/tree) Single Tree Live Wood carbon pool */
	TREE_TOT_DEADWOOD_C,                /* (tC/tree) Single Tree Dead Wood carbon pool */
	TREE_TOT_WOOD_C,                    /* (tC/tree) Single Tree Total Wood carbon pool */
	TREE_SAPWOOD_C,                     /* (tC/tree) Single Tree Sapwood carbon pool */
	TREE_HEARTWOOD_C,                   /* (tC/tree) Single Tree Heartwood carbon pool */

	/* carbon biomass pools in tDM/sizeCell */
	LEAF_DM,                            /* (tDM/cell) Current Leaf dry matter pool */
	FROOT_DM,                           /* (tDM/cell) Current Fine root dry matter pool */
	STEM_DM,                            /* (tDM/cell) Current Stem dry matter pool */
	CROOT_DM,                           /* (tDM/cell) Current Coarse root dry matter pool */
	BRANCH_DM,                          /* (tDM/cell) Current Branch dry matter pool */
	RESERVE_DM,                         /* (tDM/cell) Current Reserve dry matter pool */
	FRUIT_DM,                           /* (tDM/cell) Current Fruit dry matter pool */
	TOT_SAPWOOD_DM,                     /* (tDM/cell) Current Sapwood dry matter pool */
	TOT_HEARTWOOD_DM,                   /* (tDM/cell) Current Heartwood dry matter pool */
	STEM_SAPWOOD_DM,                    /* (tDM/cell) Current Live Stem dry matter pool */
	STEM_HEARTWOOD_DM,                  /* (tDM/cell) Current Dead Stem dry matter pool */
	CROOT_SAPWOOD_DM,                   /* (tDM/cell) Current Live Coarse root dry matter pool */
	CROOT_HEARTWOOD_DM,                 /* (tDM/cell) Current Dead Coarse root dry matter pool */
	BRANCH_SAPWOOD_DM,                  /* (tDM/cell) Current Live Branch dry matter pool */
	BRANCH_HEARTWOOD_DM,                /* (tDM/cell) Current Dead Branch dry matter pool */
	TOT_LIVEWOOD_DM,                    /* (tDM/cell) Current Live Wood dry matter pool */
	TOT_DEADWOOD_DM,                    /* (tDM/cell) Current Dead Wood dry matter pool */
	STEM_LIVEWOOD_DM,                   /* (tDM/cell) Current Live Stem dry matter pool */
	STEM_DEADWOOD_DM,                   /* (tDM/cell) Current Dead Stem dry matter pool */
	CROOT_LIVEWOOD_DM,                  /* (tDM/cell) Current Live Coarse root dry matter pool */
	CROOT_DEADWOOD_DM,                  /* (tDM/cell) Current Dead Coarse root dry matter pool */
	BRANCH_LIVEWOOD_DM,                 /* (tDM/cell) Current Live Branch dry matter pool */
	BRANCH_DEADWOOD_DM,                 /* (tDM/cell) Current Dead Branch dry matter pool */

	/* Above and Below Ground Biomass */
	AGB,                                /* (tC/cell) Above Ground Biomass pool */
	OLD_AGB,                            /* (tC/cell) Previous year Above Ground Biomass pool */
	BGB,                                /* (tC/cell) Below Ground Biomass pool */
	OLD_BGB,                            /* (tC/cell) Previous year Below Ground Biomass pool */
	DELTA_AGB,                          /* (tC/cell/year) Above Ground Biomass increment */
	DELTA_BGB,                          /* (tC/cell/year) Below Ground Biomass increment */
	TREE_AGB,                           /* (tC/tree) Above Ground Biomass pool */
	TREE_OLD_AGB,                       /* (tC/cell) Previous year Above Ground Biomass pool */
	TREE_BGB,                           /* (tC/tree) Below Ground Biomass pool */
	TREE_OLD_BGB,                       /* (tC/cell) Previous year Below Ground Biomass pool */
	DELTA_TREE_AGB,                     /* (tC/tree/year) Above Ground Biomass increment */
	DELTA_TREE_BGB,                     /* (tC/tree/year) Below Ground Biomass increment */

	/* converted fraction for allocation */
	FROOT_LEAF_FRAC,                    /* (DIM) Fraction fine root to leaf */
	STEM_LEAF_FRAC,                     /* (DIM) Fraction stem to leaf */
	CROOT_STEM_FRAC,                    /* (DIM) Fraction coarse root to stem */
	LIVE_TOTAL_WOOD_FRAC,               /* (DIM) Fraction live to total stem wood */
	FINE_COARSE_ROOT,                   /* (DIM) Fraction fine to coarse root */
	FRACBB,                             /* (DIM) Fraction branch and bark to stem */

	/* Maintenance respiration */
	DAILY_LEAF_MAINT_RESP,              /* (gC/m2/day) Day-time leaf maintenance respiration */
	DAILY_LEAF_SUN_MAINT_RESP,          /* (gC/m2/day) Day-time leaf sun maintenance respiration */
	DAILY_LEAF_SHADE_MAINT_RESP,        /* (gC/m2/day) Day-time leaf shade maintenance respiration */
	NIGHTLY_LEAF_MAINT_RESP,            /* (gC/m2/day) Night-time leaf maintenance respiration */
	NIGHTLY_LEAF_SUN_MAINT_RESP,        /* (gC/m2/day) Night-time leaf sun maintenance respiration */
	NIGHTLY_LEAF_SHADE_MAINT_RESP,      /* (gC/m2/day) Night-time leaf shade maintenance respiration */
	TOT_LEAF_MAINT_RESP,                /* (gC/m2/day) Leaf maintenance respiration */
	FROOT_MAINT_RESP,                   /* (gC/m2/day) Fine root maintenance respiration */
	STEM_MAINT_RESP,                    /* (gC/m2/day) Stem maintenance respiration */
	CROOT_MAINT_RESP,                   /* (gC/m2/day) Coarse root maintenance respiration */
	BRANCH_MAINT_RESP,                  /* (gC/m2/day) Branch maintenance respiration */
	TOTAL_MAINT_RESP,                   /* (gC/m2/day) Total maintenance respiration */
	MONTHLY_TOTAL_MAINT_RESP,           /* (gC/m2/month) Monthly Total maintenance respiration */
	YEARLY_TOTAL_MAINT_RESP,            /* (gC/m2/year) Yearly Total maintenance respiration */
	TOTAL_MAINT_RESP_tC,                /* (tC/cell/day) Total maintenance respiration */

	/* Growth respiration */
	EFF_GRPERC,                         /* (DIM) Growth respiration ratio based on age */
	LEAF_GROWTH_RESP,                   /* (gC/m2/day) Leaf growth respiration */
	FROOT_GROWTH_RESP,                  /* (gC/m2/day) Fine root growth respiration */
	STEM_GROWTH_RESP,                   /* (gC/m2/day) Stem growth respiration */
	BRANCH_GROWTH_RESP,                 /* (gC/m2/day) Branch growth respiration */
	CROOT_GROWTH_RESP,                  /* (gC/m2/day) Coarse root growth respiration */
	TOTAL_GROWTH_RESP,                  /* (gC/m2/day) Total growth respiration */
	MONTHLY_TOTAL_GROWTH_RESP,          /* (gC/m2/monthly) Monthly Total growth respiration */
	YEARLY_TOTAL_GROWTH_RESP,           /* (gC/m2/year) Yearly Total growth respiration */
	TOTAL_GROWTH_RESP_tC,               /* (tC/cell/day) Total growth respiration */

	/* Autotrophic respiration */
	LEAF_AUT_RESP,                      /* (gC/m2/day) Leaf autotrophic respiration */
	FROOT_AUT_RESP,                     /* (gC/m2/day) Fine root autotrophic respiration */
	STEM_AUT_RESP,                      /* (gC/m2/day) Stem autotrophic respiration */
	CROOT_AUT_RESP,                     /* (gC/m2/day) Coarse root autotrophic respiration */
	BRANCH_AUT_RESP,                    /* (gC/m2/day) Branch autotrophic respiration */
	TOTAL_AUT_RESP,                     /* (gC/m2/day) Total autotrophic respiration */
	MONTHLY_LEAF_AUT_RESP,              /* (gC/m2/month) Monthly Leaf autotrophic respiration */
	MONTHLY_FROOT_AUT_RESP,             /* (gC/m2/month) Monthly Fine root autotrophic respiration */
	MONTHLY_STEM_AUT_RESP,              /* (gC/m2/month) Monthly Stem autotrophic respiration */
	MONTHLY_CROOT_AUT_RESP,             /* (gC/m2/month) Monthly Coarse root autotrophic respiration */
	MONTHLY_BRANCH_AUT_RESP,            /* (gC/m2/month) Monthly Branch autotrophic respiration */
	MONTHLY_TOTAL_AUT_RESP,             /* (gC/m2/month) Monthly Total autotrophic respiration */
	YEARLY_LEAF_AUT_RESP,               /* (gC/m2/year) Annual Leaf autotrophic respiration */
	YEARLY_FROOT_AUT_RESP,              /* (gC/m2/year) Annual Fine root autotrophic respiration */
	YEARLY_STEM_AUT_RESP,               /* (gC/m2/year) Annual Stem autotrophic respiration */
	YEARLY_CROOT_AUT_RESP,              /* (gC/m2/year) Annual Coarse root autotrophic respiration */
	YEARLY_BRANCH_AUT_RESP,             /* (gC/m2/year) Annual Branch autotrophic respiration */
	YEARLY_TOTAL_AUT_RESP,              /* (gC/m2/year) Annual Total autotrophic respiration */
	TOTAL_AUT_RESP_tC,                  /* (tC/cell/day) Total autotrophic respiration */

	/* nitrogen to nitrogen pools */
	N_TO_LEAF,                          /* (tN/cell/day) Daily Net Primary Production to Leaf pool */
	N_TO_ROOT,                          /* (tN/cell/day) Daily Net Primary Production to Root pool */
	N_TO_FROOT,                         /* (tN/cell/day) Daily Net Primary Production to Fine root pool */
	N_TO_CROOT,                         /* (tN/cell/day) Daily Net Primary Production to Coarse root pool */
	N_TO_TOT_STEM,                      /* (tN/cell/day) Daily Net Primary Production to Stem and Branch pool */
	N_TO_STEM,                          /* (tN/cell/day) Daily Net Primary Production to Stem pool */
	N_TO_BRANCH,                        /* (tN/cell/day) Daily Net Primary Production to Branch pool */
	N_TO_RESERVE,                       /* (tN/cell/day) Daily Net Primary Production to Reserve pool */
	N_TO_FRUIT,                         /* (tN/cell/day) Daily Net Primary Production to Fruit pool */
	N_TO_LITTER,                        /* (tN/cell/day) Daily Net Primary Production to Litter pool */
	N_TO_LITR,                          /* (tN/cell/day) Daily Litter to Litter nitrogen pool */
	LEAF_N_TO_REMOVE,                   /* (tN/cell/day) Daily Leaf to remove */
	FROOT_N_TO_REMOVE,                  /* (tN/cell/day) Daily Fine Root to remove */
	FRUIT_N_TO_REMOVE,                  /* (tN/cell/day) Daily Fruit to remove */
	BRANCH_N_TO_REMOVE,                 /* (tN/cell/day) Annual Branch to remove */
	CROOT_N_TO_REMOVE,                  /* (tN/cell/day) Annual Coarse root to remove */

	LEAF_FALL_N,                        //(tN/cell/day) Daily Leaf to remove -evergreen - constant for each day 
	FROOT_FALL_N,                       //(tN/cell/day) Daily froot to remove -evergreen - constant for each day
	FRUIT_FALL_N,                       //(tN/cell/day) Daily fruit to remove -evergreen - constant for each day

	N_TO_CWD,                           /* (tN/cell/day) Daily Woody Nitrogen to Coarse Woody Debris */
	N_LEAF_TO_LITR,                     /* (tN/cell/day) Daily Leaf Litter to Litter labile nitrogen pool */
	N_FROOT_TO_LITR,                    /* (tN/cell/day) Daily Leaf Litter to Litter labile nitrogen pool */
	//N_FRUIT_TO_LITR,                    /* (tN/cell/day) Daily Fruit to Litter labile nitrogen pool */
	N_LEAF_TO_LITR1N,                   /* (tN/cell/day) Daily Leaf Litter to Litter labile nitrogen pool */
	N_LEAF_TO_LITR2N,                   /* (tN/cell/day) Daily Leaf Litter to Litter unshielded cellulose nitrogen pool */
	N_LEAF_TO_LITR3N,                   /* (tN/cell/day) Daily Leaf Litter to Litter shielded cellulose nitrogen pool */
	N_LEAF_TO_LITR4N,                   /* (tN/cell/day) Daily Leaf Litter to Litter lignin nitrogen pool */
	N_FROOT_TO_LITR1N,                  /* (tN/cell/day) Daily Fine Root Litter to Litter labile nitrogen pool */
	N_FROOT_TO_LITR2N,                  /* (tN/cell/day) Daily Fine Root Litter to Litter unshielded cellulose nitrogen pool */
	N_FROOT_TO_LITR3N,                  /* (tN/cell/day) Daily Fine Root Litter to Litter shielded cellulose nitrogen pool */
	N_FROOT_TO_LITR4N,                  /* (tN/cell/day) Daily Fine Root Litter to Litter lignin nitrogen pool */
	N_STEM_TO_CWD,                      /* (tN/cell/day) tN from Stem pool to Coarse Woody Debris */
	N_CROOT_TO_CWD,                     /* (tN/cell/day) tN from Coarse root pool to Coarse Woody Debris */
	N_BRANCH_TO_CWD,                    /* (tN/cell/day) tN from Branch pool to Coarse Woody Debris */
	N_RESERVE_TO_CWD,                   /* (tN/cell/day) tN from Reserve pool to Coarse Woody Debris */
	N_FRUIT_TO_CWD,                     /* (tN/cell/day) tN from Fruit pool to Coarse Woody Debris */
	N_LEAF_TO_RESERVE,                  /* (tN/cell/day) Daily retranslocated C from Leaf pool to Reserve */
	N_FROOT_TO_RESERVE,                 /* (tN/cell/day) Daily retranslocated C from Fine root pool to Reserve */
	N_BRANCH_TO_RESERVE,                /* (tN/cell/day) Annual retranslocated C from Branch pool to Reserve */
	N_CROOT_TO_RESERVE,                 /* (tN/cell/day) Annual retranslocated C from Coarse root pool to Reserve */
	N_STEM_LIVEWOOD_TO_DEADWOOD,        /* (tN/cell/day) Daily tN from Stem live wood pool to Stem dead wood */
	N_CROOT_LIVE_WOOD_TO_DEADWOOD,      /* (tN/cell/day) Daily tN from Coarse live wood pool to Coarse dead wood */
	N_BRANCH_LIVE_WOOD_TO_DEAD_WOOD,    /* (tN/cell/day) Daily tN from Branch live wood pool to Branch dead wood */

	YEARLY_N_TO_WOOD,                   /* (tN/cell/year) Annual Nitrogen stocked into wood pool */
	CUM_YEARLY_N_TO_WOOD,               /* (tN/cell) Cumulated Annual Nitrogen stocked into wood pool */

	/* nitrogen pools (tN/area) */
	LEAF_N,                             /* (tN/cell) Current Leaf nitrogen pool */
	LEAF_SUN_N,                         /* (tN/cell) Current Leaf sun nitrogen pool */
	LEAF_SHADE_N,                       /* (tN/cell) Current Leaf shade nitrogen pool */
	LEAF_FALLING_N,                     /* (tN/cell) Current Leaf falling nitrogen pool */
	FROOT_N,                            /* (tN/cell) Current Fine root nitrogen pool */
	DEADWOODN,                          /* (tN/cell) Current Coarse woody debris falling nitrogen pool */
	CROOT_N,                            /* (tN/cell) Current Coarse root nitrogen pool */
	CROOT_LIVEWOOD_N,                   /* (tN/cell) Current Live Coarse root nitrogen pool */
	CROOT_DEADWOOD_N,                   /* (tN/cell) Current Dead Coarse root nitrogen pool */
	STEM_N,                             /* (tN/cell) Current Stem nitrogen pool */
	STEM_LIVEWOOD_N,                    /* (tN/cell) Current Live Stem nitrogen pool */
	STEM_DEADWOOD_N,                    /* (tN/cell) Current Dead Stem nitrogen pool */
	BRANCH_N,                           /* (tN/cell) Current Branch nitrogen pool */
	BRANCH_LIVEWOOD_N,                  /* (tN/cell) Current Live Branch nitrogen pool */
	BRANCH_DEADWOOD_N,                  /* (tN/cell) Current Dead Branch nitrogen pool */
	RESERVE_N,                          /* (tN/cell) Current Reserve nitrogen pool */
	FRUIT_N,                            /* (tN/cell) Current Fruit nitrogen pool */
	MAX_FRUIT_N,                        /* (tN/cell/year) Annual Fruit nitrogen pool */
	TREE_N_DEMAND,                      /* (tN/cell) Current Nitrogen demand for new plant tissues */
	LITR_N,                             /* (tN/cell) Current Litter nitrogen pool */
	CWD_N,                              /* (tN/cell) Current Coarse Woody Debris nitrogen pool */
	TOT_STEM_N,                         /* (tN/cell) Current Total Stem nitrogen pool */
	TOT_ROOT_N,                         /* (tN/cell) Current Total root nitrogen pool */
	DEADWOOD_LITTER_N,                  /* (tN/cell) Current deadwood litter nitrogen pool */
	TOTAL_N,                            /* (tN/cell) Current Total nitrogen pool */

	/* per tree in tN */
	TREE_LEAF_N,                        /* (tN/tree) Single Tree Leaf nitrogen pool */
	TREE_STEM_N,                        /* (tN/tree) Single Tree Stem nitrogen pool */
	TREE_TOT_STEM_N,                    /* (tN/tree) Single Tree Stem + Branch nitrogen pool */
	TREE_ROOT_N,                        /* (tN/tree) Single Tree Total root nitrogen pool */
	TREE_FROOT_N,                       /* (tN/tree) Single Tree Fine root nitrogen pool */
	TREE_CROOT_N,                       /* (tN/tree) Single Tree Coarse nitrogen pool */
	TREE_BRANCH_N,                      /* (tN/tree) Single Tree Branch nitrogen pool */
	TREE_RESERVE_N,                     /* (tN/tree) Single Tree Reserve nitrogen pool */
	TREE_FRUIT_N,                       /* (tN/tree) Single Tree Fruit nitrogen pool */

	/* litter and soil flux */
	LEAF_TO_LITRC,                      /* (tC/cell/day) leaf litter  C flux */
	LEAF_TO_LITR1C,                     /* (tC/cell/day) leaf litter labile C flux */
	RESERVE_TO_LITRC,                      /* (tC/cell/day) reserve litter  C flux */
	RESERVE_TO_LITR1C,                     /* (tC/cell/day) reserve litter labile C flux */
	LEAF_TO_LITR2C,                     /* (tC/cell/day) leaf litter unshielded cell/dayulose C flux */
	LEAF_TO_LITR3C,                     /* (tC/cell/day) leaf litter shielded cell/dayulose C flux */
	LEAF_TO_LITR4C,                     /* (tC/cell/day) leaf litter lignin C flux */
	LEAF_TO_LITRN,                      /* (tN/cell/day) leaf litter  N flux */
	LEAF_TO_LITR1N,                     /* (tN/cell/day) leaf litter labile N flux */
	LEAF_TO_LITR2N,                     /* (tN/cell/day) leaf litter unshielded Cellulose N flux */
	LEAF_TO_LITR3N,                     /* (tN/cell/day) leaf litter shielded Cellulose N flux */
	LEAF_TO_LITR4N,                     /* (tN/cell/day) leaf litter lignin N flux */
	FROOT_TO_LITRC,                     /* (tC/cell/day) froot litter labile C flux */
	FROOT_TO_LITR1C,                    /* (tC/cell/day) froot litter labile C flux */
	FROOT_TO_LITR2C,                    /* (tC/cell/day) froot litter unshielded celluse C flux */
	FROOT_TO_LITR3C,                    /* (tC/cell/day) froot litter shielded cell/dayulose C flux */
	FROOT_TO_LITR4C,                    /* (tC/cell/day) froot litter lignin C flux */
	FROOT_TO_LITRN,                     /* (tN/cell/day) froot litter labile N flux */
	FROOT_TO_LITR1N,                    /* (tN/cell/day) froot litter labile N flux */
	FROOT_TO_LITR2N,                    /* (tN/cell/day) froot litter unshielded Cellulose N flux */
	FROOT_TO_LITR3N,                    /* (tN/cell/day) froot litter shielded Cellulose N flux */
	FROOT_TO_LITR4N,                    /* (tN/cell/day) froot litter lignin N flux */
	CWD_TO_LITRC,                       /* (tC/cell/day) Coarse Woody Debris flux */
	CWD_TO_LITR2C,                      /* (tC/cell/day) Coarse Woody Debris unshielded flux */
	CWD_TO_LITR3C,                      /* (tC/cell/day) Coarse Woody Debris shielded flux */
	CWD_TO_LITR4C,                      /* (tC/cell/day) Coarse Woody Debris lignin flux */
	CWD_TO_LITRN ,                      /* (tN/cell/day) Coarse Woody Debris flux */
	CWD_TO_LITR2N,                      /* (tN/cell/day) Coarse Woody Debris unshielded flux */
	CWD_TO_LITR3N,                      /* (tN/cell/day) Coarse Woody Debris shielded flux */
	CWD_TO_LITR4N,                      /* (tN/cell/day) Coarse Woody Debris lignin flux */

	/* litter and soil pools */
	LEAF_LITRC,                         /* (tC/cell) leaf litter C pool */
	LEAF_LITR1C,                        /* (tC/cell) leaf litter labile C pool */
	LEAF_LITR2C,                        /* (tC/cell) leaf litter unshielded cellulose C pool */
	LEAF_LITR3C,                        /* (tC/cell) leaf litter shielded cellulose C pool */
	LEAF_LITR4C,                        /* (tC/cell) leaf litter lignin C pool */
	LEAF_LITRN,                         /* (tN/cell) leaf litter N pool */
	LEAF_LITR1N,                        /* (tN/cell) leaf litter labile N pool */
	LEAF_LITR2N,                        /* (tN/cell) leaf litter unshielded Cellulose N pool */
	LEAF_LITR3N,                        /* (tN/cell) leaf litter shielded Cellulose N pool */
	LEAF_LITR4N,                        /* (tN/cell) leaf litter lignin N pool */
	FROOT_LITRC,                        /* (tC/cell) froot litter C pool */
	FROOT_LITR1C,                       /* (tC/cell) froot litter labile C pool */
	FROOT_LITR2C,                       /* (tC/cell) froot litter unshielded cellulose C pool */
	FROOT_LITR3C,                       /* (tC/cell) froot litter shielded cellulose C pool */
	FROOT_LITR4C,                       /* (tC/cell) froot litter lignin C pool */
	FROOT_LITRN,                        /* (tN/cell) froot litter N pool */
	FROOT_LITR1N,                       /* (tN/cell) froot litter labile N pool */
	FROOT_LITR2N,                       /* (tN/cell) froot litter unshielded Cellulose N pool */
	FROOT_LITR3N,                       /* (tN/cell) froot litter shielded Cellulose N pool */
	FROOT_LITR4N,                       /* (tN/cell) froot litter lignin N pool */
	CWD_LITRC ,                         /* (tC/cell) Coarse Woody Debris pool */
	CWD_LITR2C,                         /* (tC/cell) Coarse Woody Debris unshielded pool */
	CWD_LITR3C,                         /* (tC/cell) Coarse Woody Debris shielded pool */
	CWD_LITR4C,                         /* (tC/cell) Coarse Woody Debris lignin pool */
	CWD_LITRN,                          /* (tN/cell) Coarse Woody Debris pool */
	CWD_LITR2N,                         /* (tN/cell) Coarse Woody Debris unshielded pool */
	CWD_LITR3N,                         /* (tN/cell) Coarse Woody Debris shielded pool */
	CWD_LITR4N,                         /* (tN/cell) Coarse Woody Debris lignin pool */

	/* fractions */
	LEAF_LITR_SCEL_FRAC,                /* (DIM) leaf litter shielded cellulose fraction */
	LEAF_LITR_USCEL_FRAC,               /* (DIM) leaf litter unshielded cellulose fraction */
	FROOT_LITR_SCEL_FRAC,               /* (DIM) fine root litter shielded cellulose fraction */
	FROOT_LITR_USCEL_FRAC,              /* (DIM) fine root litter unshielded fraction */
	DEADWOOD_SCEL_FRAC,                 /* (DIM) dead wood litter shielded cellulose fraction */
	DEADWOOD_USCEL_FRAC,                /* (DIM) dead wood litter unshielded fraction */

	/*** BALANCES VALUES ***/
	TREEC_IN,                           /* tree carbon balance in */
	TREEC_OUT,                          /* tree carbon balance out */
	TREEC_STORE,                        /* tree carbon balance store */
	TREEC_OLDSTORE,                     /* tree carbon balance old_store */
	TREEC_BALANCE,                      /* tree carbon balance balance */

	TREEW_IN,                           /* tree water balance in */
	TREEW_OUT,                          /* tree water balance out */
	TREEW_STORE,                        /* tree water balance store */
	TREEW_OLDSTORE,                     /* tree water balance old_store */
	TREEW_BALANCE,                      /* tree water balance balance */

	TREEN_IN,                           /* tree nitrogen balance in */
	TREEN_OUT,                          /* tree nitrogen balance out */
	TREEN_STORE,                        /* tree nitrogen balance store */
	TREEN_OLDSTORE,                     /* tree nitrogen balance old_store */
	TREEN_BALANCE,                      /* tree nitrogen balance balance */

	VALUES
};

enum {

	N_TREE,                             /* (tree/cell) Live Tree per class */
	DEAD_TREE,                          /* (tree/cell) Dead Tree per class */
	THINNED_TREE,                       /* (tree/cell) Tree removed for thinning or harvesting */
	tree_remove_crowded ,           // (tree/cell) number to be remove when self-thinning is activated
	N_STUMP,                            /* (stump/cell) Number of stumps for coppice per class */
	N_SEED,                             /* (seeds/cell) Seeds number per class */
	N_TREE_SAP,                         /* (saplings/cell) Numbers of Saplings per class */
	VEG_UNVEG,                          /* (DIM) vegetative-un vegetative period per class */
	VEG_DAYS,                           /* (days/year) day of vegetative period for class if != 0 is in veg period */
	FIRST_VEG_DAYS,                     /* (DIM) first annual day of veg period */
	YEARLY_VEG_DAYS,                    /* (days/year) annual number of vegetative days */
	BUD_BURST_COUNTER,                  /* (days/year) counter days of budburst per class */
	DAYS_FOR_BUDBURST,                  /* (days/year) days of budburst per class */
	BUD_BURST_DAY_COUNTER,              /* (days/year) days of budburst per class */
	BUD_BURST_WEIGHT,                   /* */
	BUDBURST_A,                         /* */
	BUDBURST_B,                         /* */
	DAYS_LEAFFALL,                      /* (days/year) days of leaf fall per class */
	DAY_VEG_FOR_LEAF_FALL,              /* (days/year) days for leaf fall */
	LEAF_FALL_COUNTER,                  /* (DIM) counter for leaf fall */
	SENESCENCE_DAY_ONE,                 /* (DIM) doy at first day of the senescence */
	DAYS_WITH_CANOPY_WET,               /* (days/year) with canopy wet */
	DAYS_WITH_CANOPY_SNOW,              /* (days/year) with canopy snow */
	DOS,                                /* (days) day of simulation */
	MOS,                                /* (months) month of simulation */
	YOS,                                /* (years) year of simulation */
	YEARS_THINNING,
	                                    // FIXME the followings should be moved to a different ENUM
	THINNING_HAPPENS,                   /* 0 = thinning off, 1 = thinning on */
	HARVESTING_HAPPENS,                 /* 0 = harvesting off, 1 = harvesting on */
	                                     
	COUNTERS
};

typedef struct
{
	int year_stand;
	int x;
	int y;
	int age;
	char *species;
	e_management management;
	int n;
	int stool;
	double avdbh;
	double height;
	double wf;
	double wrc;
	double wrf;
	double ws;
	double wbb;
	double wres;
	double lai;
} row_t;

typedef struct
{
	row_t* rows;
	int rows_count;

} dataset_t;

typedef struct
{
	int index;

	e_management management;

	char *name;                         /** species name **/
	int phenology_phase;                /* species-specific phenological phase */
	int counter[COUNTERS];
	double value[VALUES];
} species_t;

typedef struct
{
	int index;

	species_t *species;

	int value;                          /** age value (year) **/
	int species_count;                  /* number of different species */
	int species_avail;

	/* for logger */
	int initial_species_count;
} age_t;

typedef struct
{
	int index;

	age_t* ages;

	double value;                       /** dbh value (cm) **/
	int ages_count;                     /* number of different cohorts (ages) */
	int ages_avail;

	/* for logger */
	int initial_ages_count;

	int dbh_n_trees;                    /* number of trees per dbh class */
	double dbh_density;                 /* density of trees per dbh class */
} dbh_t;

typedef struct
{
	int layer_z;                        //fixme remove
	int layer_n_height_class;           /* number of height class per layer */
	int layer_n_trees;                  /* number of trees per layer */
	int layer_height_class_counter;
	int canopy_int_layer_height_class_counter;
	

	double layer_density;               /* tree density per layer (n_tree/sizecell) */
	double layer_cover_proj;            /* layer canopy cover projected per layer */
	double layer_cover_exp;             /* layer canopy cover exposed per layer */
	double daily_layer_cover_proj;      /* layer canopy cover projected per layer */
	//double daily_layer_cover_exp;       /* layer canopy cover exposed per layer */
	double layer_avg_tree_height;       /* average tree layer height */
	double layer_tree_height_modifier;  /* layer level tree height modifier (see Wallace et al., 1991; Cannell & Grace 1993) */

} tree_layer_t;

typedef struct
{
	int index;

	dbh_t *dbhs;

	double value;                       /** height value (m) **/
	int dbhs_count;                     /* number of different dbh classes */
	int dbhs_avail;

	/* for logger */
	int initial_dbhs_count;
	int height_z;
	double height_density;              /* tree density per height class (n_tree/sizecell) */
	double height_cover;                /* height class cover per height class (n_tree/sizecell) */
} height_t;

enum
{
	SOIL_CARBON,                    /* SOIL CARBON (gC/m2 soil) */
	PREVIOUS_ASW,
	WATER_BALANCE,
	OLD_WATER_BALANCE,
	MOIST_RATIO,
	VSWC,                           /* VOLUMETRIC SOIL WATER CONTENT (%VOL) */
	PSI_SAT,                        /* SOIL SATURATED MATRIC POTENTIAL */
	VWC,
	VWC_FC,                         /* SOIL FIELD CAPACITY VOLUMETRIC WATER CONTENT */
	VWC_SAT,                        /* SOIL SATURATED VOLUMETRIC WATER CONTENT */
	SOIL_B,                         /* SOIL MOISTURE PARAMETER */
	SOIL_W_SAT,                     /* (KGH2O/M2) SOILWATER AT SATURATION */
	SOIL_W_FC,                      /* (KGH2O/M2) SOILWATER AT FIELD CAPACITY */
	DAILYSOILEVAPORATIONWATT,
	DAILYSOILLATENTHEATFLUX,
	DAILYSOILSENSIBLEHEATFLUX,


	SOIL_COUNT
};

typedef struct
{
	double value[SOIL_COUNT];

} soil_layer_s;

typedef struct
{
	/* tree height struct */
	height_t* heights;
	int heights_count;
	int heights_avail;

	/* tree layer struct */
	tree_layer_t* tree_layers;
	int tree_layers_count;
	int tree_layers_count_temp;
	int t_layers_avail;

	/* soil layer struct */
	soil_layer_s* soil_layers;
	int soil_layers_count;
	int s_layers_avail;

	int years_count;
	meteo_annual_t *years;

	int doy;                                                              /* day of the year */
	int dos;                                                              /* day of simulation */

    // FLAGS 
	int harvesting;                                                       /* 0 = off, 1 = on */
    int removal;                                                          // 0=off, 1 =0n (set if removal in the stand is performed)
    int GREFFMORT_HAPPENS;                 /* 0 = no classes removed because of growth efficiency mortality, 1 = at least one class removed  */
   
	/* for logger */
	int initial_tree_layers_count;
	int initial_heights_count;

	/* overall cell counter */
	int cell_tree_layers_count;
	int cell_heights_count;
	int cell_dbhs_count;                                                  /* (dbh/cell) number of different dbh */
	int cell_ages_count;                                                  /* (ages/cell) number of different ages */
	int cell_species_count;                                               /* (species/cell) number of different species */
	int cell_height_class_counter;
	int canopy_int_cell_height_class_counter;
    int cell_age;                                                         // age of the stand (abovegraound/soprassuolo)
 
    int cell_coppice_count;                                             //  reduction counter for coppice density
   
	int year_stand;                                                       /* input stand.txt row stand year : used when man = var */
 
	int x;                                                                /* cell index within the matrix */
	int y;                                                                /* cell index within the matrix */

	/* topo file */
	//double elev;                                                          /* cell elevation (m) */

	int north;                                                            /* northern hemisphere north = 0, south hemisphere south = 1 */

	int doy_daylength;

	//fixme ALESSIOR move to meteo struct!
	double previous_thermic_sum;
	double previous_solar_rad;
	double previous_tavg;
	double previous_tmax;
	double previous_tmin;
	double previous_vpd;
	double previous_ts_f;
	double previous_prcp;
	double previous_swc;
	double previous_ndvi_lai;

	/* annual met values */
	double abscission_daylength;                                          /* (hrs/day) day length for starting abscission (not used) */

	/*forest structure variables*/	
	int Veg_Counter;                                                      /* number of class in vegetative period  in the cell */
	double daily_cell_cover_proj;                                         /* (ratio) Daily canopy cover projected */
	double daily_cell_cover_exp;                                          /* (ratio) Daily canopy cover exposed */
	double cell_cover_proj;                                               /* (ratio) Canopy cover projected */
	double cell_cover_exp;                                                /* (ratio) Canopy cover exposed */
	double bare_soil_cover;                                               /* (ratio) Bare soil cover */
	int n_trees;			                                              /* (trees/cell/day) number of trees */
	int daily_dead_tree;                                                  /* (trees/cell/day) daily number of dead tree */
	int monthly_dead_tree;                                                /* (trees/cell/month) monthly number of dead tree */
	int annual_dead_tree ;                                                /* (trees/cell/year) annual number of dead tree */
	int saplings_counter;                                                 /* (trees/cell) number of class as sapling */
	
	double basal_area;                                                    /* (m2/cell) cumulated basal area at cell level */
	double agb;                                                           /* (tC/cell) current above ground biomass at cell level */
	double bgb;                                                           /* (tC/cell) current below ground biomass at cell level */
	double volume;                                                        /* (m3/cell) current stem volume */
	double cum_volume;                                                    /* (m3/cell) cumulated stem volume */
	double volume2;                                                        /* (m3/cell) current stem volume using conversion factors for volume*/
	double cum_volume2;                                                    /* (m3/cell) cumulated stem volume using conversion factorsfor volume */
	double standing_wood;                                                 /* (tC/cell) standing wood biomass */
	double yearly_C_to_wood;                                              /* (tC/cell/year) Annual Carbon stocked into wood pool */
	double cum_yearly_C_to_wood;                                          /* (tC/cell) Cumulated Annual Carbon stocked into wood pool */
	double cum_npp;                                                       /* (tC/cell) Cumulated Annual npp */

	double lai;                                                           //fixme
	double max_lai_proj;                                                  //fixme

	/* radiation variables */
	/* short wave radiation */
	double sw_rad_down_MJ;                                                /* (MJ/m2/day) Downward short wave radiation flux */
	double short_wave_radiation_upward_MJ;                                /* (MJ/m2/day) Upward short wave radiation flux */
	double short_wave_radiation_upward_W;                                 /* (W/m2) Upward short wave radiation flux */
	double sw_rad_refl;                                                   /* (W/m2) Short Wave radiation reflected */
	double sw_rad_refl_soil;                                              /* (W/m2) Short Wave radiation reflected from soil level */
	double sw_rad_refl_snow;                                              /* (W/m2) (W/day)  Wave radiation reflected from snow level */
	double sw_rad_abs;                                                    /* (W/m2) Short Wave radiation absorbed tree level */
	double sw_rad_abs_soil;                                               /* (W/m2) Short Wave radiation absorbed soil level */
	double sw_rad_abs_snow;                                               /* (W/m2) Short Wave radiation absorbed snow level */

	double net_rad_abs;
	double net_rad_refl;
	double net_rad_refl_soil;
	double net_rad_refl_snow;
	double net_rad_abs_soil;
	double net_rad_abs_snow;
	double temp_net_rad_abs;                                              /* (W/m2) temporary absorbed Net radiation for layer */
	double temp_net_rad_refl;                                             /* (W/m2) temporary reflected Net radiation for layer */

	/* long wave radiation */
	double long_wave_radiation_upward_MJ;                                 /* (MJ/m2/day) Upward long wave radiation flux */
	double long_wave_radiation_upward_W;                                  /* (W/m2) Upward long wave radiation flux */
	double net_long_wave_radiation_MJ;                                    /* (MJ/m2/day) Net long wave radiation flux  */
	double net_long_wave_radiation_W;                                     /* (W/m2) Net long wave radiation flux */
	double long_wave_absorbed;                                            /* (W/m2) Total absorbed long wave wave radiation flux */
	double long_wave_reflected;                                           /* (W/m2) Total reflected long wave radiation flux */
	double long_wave_emitted;                                             /* (W/m2) Total emitted long wave radiation flux */
	double soil_long_wave_emitted;                                        /* (W/m2) Soil emitted long wave radiation flux */
	double net_lw_rad_for_soil;                                           /* (W/m2) Net Long Wave radiation to soil level */
	double net_rad_for_soil;                                              /* (W/m2) Net radiation to soil level */
	double temp_sw_rad_abs;                                               /* (W/m2) temporary absorbed short wave for layer */
	double temp_sw_rad_refl;                                              /* (W/m2) temporary reflected short wave for layer */

	/* PAR radiation */
	double apar;                                                          /* (molPAR/m2/day) cumulated absorbed PAR at tree level */
	double apar_soil;                                                     /* (molPAR/m2/day) cumulated absorbed PAR at soil level */
	double apar_snow;                                                     /* (molPAR/m2/day) cumulated absorbed PAR at snow level */
	double par_transm;                                                    /* (molPAR/m2/day) cumulated transmitted PAR at cell level */
	double par_refl;                                                      /* (molPAR/m2/day) cumulated reflected PAR at cell level */
	double par_refl_soil;                                                 /* (molPAR/m2/day) reflected PAR by the soil */
	double par_refl_snow;                                                 /* (molPAR/m2/day) reflected PAR by the snow */
	double temp_apar;                                                     /* (molPAR/m2/day) temporary absorbed PAR for layer */
	double temp_par_refl;                                                 /* (molPAR/m2/day) temporary reflected PAR for layer */


	/************************************************** forest trees **************************************************/

	double canopy_temp;                                                   /* (°C) daily canopy temperature */
	double canopy_temp_k;                                                 /* (K) daily canopy temperature */
	double canopy_temp_diff;                                              /* (K) daily canopy temperature */

	/* carbon fluxes */
	double daily_gpp, monthly_gpp, annual_gpp;                            /* (gC/m2/time) daily, monthly and annual GPP at cell level */
	double daily_ass, monthly_ass, annual_ass;                            /* (gC/m2/time) daily, monthly and annual assimilation at cell level */
	double daily_gpp_tC, monthly_gpp_tC, annual_gpp_tC;                   /* (tC/m2/time) daily, monthly and annual GPP at cell level */
	double daily_npp, monthly_npp, annual_npp;                            /* (gC/m2/time) daily, monthly and annual NPP at cell level */
	double daily_npp_tC, monthly_npp_tC, annual_npp_tC;                   /* (tC/m2/time) daily, monthly and annual NPP at cell level */
	double daily_npp_tDM, monthly_npp_tDM, annual_npp_tDM;                /* (tDM/m2/time) daily, monthly and annual NPP at cell level */
	double daily_aut_resp, monthly_aut_resp, annual_aut_resp;             /* (gC/m2/time) daily, monthly and annual aut resp at cell level */
	double daily_aut_resp_tC, monthly_aut_resp_tC, annual_aut_resp_tC;    /* (tC/m2/time) daily, monthly and annual aut resp at cell level */
	double daily_maint_resp, monthly_maint_resp, annual_maint_resp;       /* (gC/m2/time) daily, monthly and annual maint resp at cell level */
	double daily_growth_resp, monthly_growth_resp, annual_growth_resp;    /* (gC/m2/time) daily, monthly and annual growth resp at cell level */
	double daily_C_flux, monthly_C_flux, annual_C_flux;                   /* (gC/m2/time) daily, monthly and annual carbon fluxes at cell level */
	double daily_leaf_carbon;                                             /* (gC/m2/day) daily carbon assimilated to c pool at cell level */
	double daily_stem_carbon;                                             /* (gC/m2/day) daily carbon assimilated to c pool at cell level */
	double daily_froot_carbon;                                            /* (gC/m2/day) daily carbon assimilated to c pool at cell level */
	double daily_croot_carbon;                                            /* (gC/m2/day) daily carbon assimilated to c pool at cell level */
	double daily_branch_carbon;                                           /* (gC/m2/day) daily carbon assimilated to c pool at cell level */
	double daily_reserve_carbon;                                          /* (gC/m2/day) daily carbon assimilated to c pool at cell level */
	double daily_fruit_carbon;                                            /* (gC/m2/day) daily carbon assimilated to fruit c pool at cell level */
	double daily_leaf_maint_resp;                                         /* (gC/m2/day) daily leaf maint resp at cell level */
	double daily_stem_maint_resp;                                         /* (gC/m2/day) daily stem maint resp at cell level */
	double daily_froot_maint_resp;                                        /* (gC/m2/day) daily fine root maint resp at cell level */
	double daily_branch_maint_resp;                                       /* (gC/m2/day) daily branch and bark maint resp at cell level */
	double daily_croot_maint_resp;                                        /* (gC/m2/day) daily coarse root maint resp at cell level */
	double daily_leaf_growth_resp;                                        /* (gC/m2/day) daily leaf growth resp at cell level */
	double daily_stem_growth_resp;                                        /* (gC/m2/day) daily stem growth resp at cell level */
	double daily_froot_growth_resp;                                       /* (gC/m2/day) daily fine root growth resp at cell level */
	double daily_branch_growth_resp;                                      /* (gC/m2/day) daily branch and bark growth resp at cell level */
	double daily_croot_growth_resp;                                       /* (gC/m2/day) daily coarse root growth resp at cell level */
	double daily_fruit_growth_resp;                                       /* (gC/m2/day) daily fruit growth resp at cell level */
	double daily_leaf_aut_resp;                                           /* (gC/m2/day) daily leaf aut resp at cell level */
	double daily_stem_aut_resp;                                           /* (gC/m2/day) daily stem aut resp at cell level */
	double daily_branch_aut_resp;                                         /* (gC/m2/day) daily branch and bark aut resp at cell level */
	double daily_froot_aut_resp;                                          /* (gC/m2/day) daily fine root aut resp at cell level */
	double daily_croot_aut_resp;                                          /* (gC/m2/day) daily coarse root aut resp at cell level */
    double annual_hwp;                         // (tC/cell/yr) cell level annual harvested woody products removed from stand - at DOY == 1  */
    double annual_Ctree_repl;                         // (tC/cell/yr) cell level annual replanted trees total C - at DOY ==1 */

	/* tree carbon pools */
	double leaf_carbon;                                                   /* (gC/m2) leaf carbon at cell level */
	double froot_carbon;                                                  /* (gC/m2) fine root carbon at cell level */
	double stem_carbon;                                                   /* (gC/m2) stem carbon at cell level */
	double stem_live_wood_carbon;                                         /* (gC/m2) stem live wood carbon at cell level */
	double stem_dead_wood_carbon;                                         /* (gC/m2) stem dead wood carbon at cell level */
	double tot_stem_carbon;                                               /* (gC/m2) total stem carbon at cell level */
	double croot_carbon;                                                  /* (gC/m2) coarse root carbon at cell level */
	double croot_live_wood_carbon;                                        /* (gC/m2) coarse root live wood carbon at cell level */
	double croot_dead_wood_carbon;                                        /* (gC/m2) coarse root dead wood carbon at cell level */
	double tot_root_carbon;                                               /* (gC/m2) total root carbon at cell level */
	double branch_carbon;                                                 /* (gC/m2) branch carbon at cell level */
	double branch_live_wood_carbon;                                       /* (gC/m2) branch live wood carbon at cell level */
	double branch_dead_wood_carbon;                                       /* (gC/m2) branch dead wood carbon at cell level */
	double reserve_carbon;                                                /* (gC/m2) reserve at cell level */
	double fruit_carbon;                                                  /* (gC/m2) fruit at cell level */
  

	/* old tree carbon pools (used for mass balance) */
	double old_leaf_carbon;
	double old_froot_carbon;
	double old_croot_carbon;
	double old_stem_carbon;
	double old_branch_carbon;
	double old_reserve_carbon;
	double old_fruit_carbon;
	double old_litrC;
	double old_soilC;
	double old_cwd_C;

	/* nitrogen fluxes */
	double daily_leaf_nitrogen;                                           /* (gN/m2/day) daily nitrogen assimilated to c pool at cell level */
	double daily_stem_nitrogen;                                           /* (gN/m2/day) daily nitrogen assimilated to c pool at cell level */
	double daily_froot_nitrogen;                                          /* (gN/m2/day) daily nitrogen assimilated to c pool at cell level */
	double daily_croot_nitrogen;                                          /* (gN/m2/day) daily nitrogen assimilated to c pool at cell level */
	double daily_root_nitrogen;                                           /* (gN/m2/day) daily nitrogen assimilated to c pool at cell level */
	double daily_branch_nitrogen;                                         /* (gN/m2/day) daily nitrogen assimilated to c pool at cell level */
	double daily_reserve_nitrogen;                                        /* (gN/m2/day) daily nitrogen assimilated to c pool at cell level */
	double daily_fruit_nitrogen;                                          /* (gN/m2/day) daily nitrogen assimilated to fruit c pool at cell level */

	/* tree nitrogen pools */
	double leaf_nitrogen;                                                 /* (gN/m2) leaf carbon at cell level */
	double froot_nitrogen;                                                /* (gN/m2) fine root carbon at cell level */
	double stem_nitrogen;                                                 /* (gN/m2) stem carbon at cell level */
	double stem_live_wood_nitrogen;                                       /* (gN/m2) stem live wood carbon at cell level */
	double stem_dead_wood_nitrogen;                                       /* (gN/m2) stem dead wood carbon at cell level */
	double tot_stem_nitrogen;                                             /* (gN/m2) total stem carbon at cell level */
	double croot_nitrogen;                                                /* (gN/m2) coarse root carbon at cell level */
	double croot_live_wood_nitrogen;                                      /* (gN/m2) coarse root live wood carbon at cell level */
	double croot_dead_wood_nitrogen;                                      /* (gN/m2) coarse root dead wood carbon at cell level */
	double tot_root_nitrogen;                                             /* (gN/m2) total root carbon at cell level */
	double branch_nitrogen;                                               /* (gN/m2) branch carbon at cell level */
	double branch_live_wood_nitrogen;                                     /* (gN/m2) branch live wood carbon at cell level */
	double branch_dead_wood_nitrogen;                                     /* (gN/m2) branch dead wood carbon at cell level */
	double reserve_nitrogen;                                              /* (gN/m2) reserve at cell level */
	double fruit_nitrogen;                                                /* (gN/m2) fruit at cell level */

	/* old tree nitrogen pools (used for mass balance) */
	double old_leaf_nitrogen;
	double old_froot_nitrogen;
	double old_croot_nitrogen;
	double old_stem_nitrogen;
	double old_branch_nitrogen;
	double old_reserve_nitrogen;
	double old_fruit_nitrogen;
	double old_litrN;
	double old_soilN;

	/* carbon use efficiency */
	double daily_cue;                                                     /* (DIM) daily carbon use efficiency */
	double monthly_cue;                                                   /* (DIM) monthly carbon use efficiency */
	double annual_cue;                                                    /* (DIM) annual carbon use efficiency */

	/* water use efficiency */
	double daily_wue;                                                     /* (DIM) daily water use efficiency */
	double daily_iwue;                                                    /* (DIM) daily intrinsic water use efficiency */
	double monthly_wue;                                                   /* (DIM) monthly water use efficiency */
	double monthly_iwue;                                                  /* (DIM) monthly intrinsic water use efficiency */
	double annual_wue;                                                    /* (DIM) annual water use efficiency */
	double annual_iwue;                                                   /* (DIM) annual intrinsic water use efficiency */

	/* energy fluxes */
	double temp_int_rain;                                                 /* (mm/m2) temporary intercepted rain */
	double temp_int_snow;                                                 /* (mm/m2) temporary intercepted snow */
	double daily_canopy_transp_watt;                                      /* (W/m2) daily canopy transpiration at cell level */
	double daily_canopy_evapo_watt;                                       /* (W/m2) daily canopy evaporation at cell level */
	double daily_canopy_et_watt;                                          /* (W/m2) daily canopy evapotranspiration at cell level */
	double daily_canopy_lh_flux;                                          /* (W/m2) daily canopy latent heat flux at cell level  */
	double daily_canopy_sh_flux;                                          /* (W/m2) current daily canopy sensible heat flux */
	double daily_lh_flux, monthly_lh_flux, annual_lh_flux;                /* (W/m2) daily, monthly and annual latent heat flux at cell level */
	double daily_sh_flux, monthly_sh_flux, annual_sh_flux;                /* (W/m2) daily, monthly and annual sensible heat flux at cell level */


	/************************************************** LITTER AND SOIL **************************************************/

	/* SOIL AND LITETR INITIALIZATION */
	/* from "soil.txt" file */
	double lat;                                                           /* latitude */
	double lon;                                                           /* longitude */
	double clay_perc;                                                     /* (%) clay percentage */
	double silt_perc;   												  /* (%) silt percentage */
	double sand_perc;                                                     /* (%) sand percentage */
	double soil_depth;                                                    /* (cm) soil depth */
	double fr;                                                            /* fertility rate */
	double fn0;                                                           /* */
	double fnn;                                                           /* */
	double m0;                                                            /* */
	double init_litter_C;                                                 /* initial litter carbon */
	double init_litter_N;                                                 /* initial litter nitrogen */
	double init_soil_C;                                                   /* initial soil carbon */
	double init_soil_N;                                                   /* initial soil nitrogen */
	double init_dead_C;                                                   /* initial dead carbon */

	/* computed internally during initialization */
	double bulk_density;                                                  /* (g/cm3) soil bulk density */
	double wilting_point;                                                 /* (mm/day) volumetric water content at wilting point */
	double field_capacity;                                                /* volumetric water content at field capacity */
	double sat_hydr_conduct;                                              /* (mm/day) saturated hydraulic conductivity (mm/day) */
	double psi;                                                           /* (MPa) water potential of soil and leaves */
	double soil_moist_ratio;                                              /* (DIM) soil moisture ratio */
	double vwc;                                                           /* (mm/day) soil volumetric water content */
	double vwc_fc;                                                        /* (mm/day) soil volumetric water content at field capacity */
	double vwc_sat;                                                       /* (mm/day) soil volumetric water content at saturation */
	double psi_sat;                                                       /* (MPa) soil matric potential */
	double soil_b;                                                        /* (DIM) slope of log(psi) vs log(rwc) */
	double soilw_fc;                                                      /* (kgH2O/day) soil water at field capacity */
	double soilw_sat;                                                     /* (kgH2O/day) soil water at saturation */

	/* soil scalars */
	double tsoil_scalar;                                                  /* (DIM) soil temperature scalar for decomposition */
	double wsoil_scalar;                                                  /* (DIM) soil water scalar  for decomposition */
	double rate_scalar;                                                   /* (DIM) soil (temperature * water) scalar for decomposition */
	double daily_gross_nmin;                                              /* (tN/cell/m2) daily gross N mineralization */
	double daily_gross_nimmob;                                            /* (tN/cell/m2) daily gross N immobilization */
	double daily_net_nmin;                                                /* (tN/cell/m2) daily net N mineralization */
	double fpi;                                                           /* (DIM) fraction of potential immobilization */

	/* litter and soil decomposition carbon fluxes */
	/* cwd + litter */
	double daily_cwd_to_litrC;                                            /* (gC/m2/day) daily cwd carbon to litter carbon pool cell level */
	double daily_cwd_to_litr2C;                                           /* (gC/m2/day) daily decomposition fluxes cwd to litter */
	double daily_cwd_to_litr3C;                                           /* (gC/m2/day) daily decomposition fluxes cwd to litter */
	double daily_cwd_to_litr4C;                                           /* (gC/m2/day) daily decomposition fluxes cwd to litter */
	double daily_leaf_to_litrC;                                           /* (gC/m2/day) daily leaf to litter carbon pool at cell level */
	double daily_leaf_to_litr1C;                                          /* (gC/m2/day) daily leaf to litter labile carbon pool at cell level */
	double daily_leaf_to_litr2C;                                          /* (gC/m2/day) daily leaf to litter unshelded cellulose carbon pool at cell level */
	double daily_leaf_to_litr3C;                                          /* (gC/m2/day) daily leaf to litter  shelded cellulose carbon pool at cell level */
	double daily_leaf_to_litr4C;                                          /* (gC/m2/day) daily leaf to litter lignin carbon pool at cell level */
	double daily_froot_to_litrC;                                          /* (gC/m2/day) daily fine root to litter carbon pool at cell level */
	double daily_froot_to_litr1C;                                         /* (gC/m2/day) daily fine root to litter labile carbon pool at cell level */
	double daily_froot_to_litr2C;                                         /* (gC/m2/day) daily fine root to litter unshelded cellulose carbon pool at cell level */
	double daily_froot_to_litr3C;                                         /* (gC/m2/day) daily fine root to litter shelded cellulose carbon pool at cell level */
	double daily_froot_to_litr4C;                                         /* (gC/m2/day) daily fine root to litter lignin carbon pool at cell level */
	double daily_reserve_to_litrC;                                          /* (gC/m2/day) daily reserve to litter  carbon pool at cell level */
	double daily_reserve_to_litr1C;                                          /* (gC/m2/day) daily reserve to litter labile carbon pool at cell level */
	double daily_to_litrC;                                                /* (gC/m2/day) daily balance of carbon to litter carbon pool at cell level */
	double daily_to_litr1C;                                               /* (gC/m2/day) daily balance of carbon to litter labile carbon pool at cell level */
	double daily_to_litr2C;                                               /* (gC/m2/day) daily balance of carbon to litter unshelded cellulose carbon pool at cell level */
	double daily_to_litr3C;                                               /* (gC/m2/day) daily balance of carbon to litter shelded cellulose carbon pool at cell level */
	double daily_to_litr4C;                                               /* (gC/m2/day) daily balance of carbon to litter lignin carbon pool at cell level */
	double daily_litr_to_soilC;                                           /* (gC/m2/day) daily litter carbon pool to soil pool */
	double daily_litr1C_to_soil1C;                                        /* (gC/m2/day) daily decomposition fluxes litter to soil */
	double daily_litr2C_to_soil2C;                                        /* (gC/m2/day) daily decomposition fluxes litter to soil */
	double daily_litr3C_to_litr2C;                                        /* (gC/m2/day) daily decomposition fluxes litter to soil */
	double daily_litr4C_to_soil3C;                                        /* (gC/m2/day) daily decomposition fluxes litter to soil */
	/* soil */
	double daily_to_soilC;                                                /* (gC/m2/day) daily litter carbon to soil */
	double daily_to_soil1C;                                               /* (gC/m2/day) daily microbial recycling pool C (fast) */
	double daily_to_soil2C;                                               /* (gC/m2/day) daily microbial recycling pool C (medium) */
	double daily_to_soil3C;                                               /* (gC/m2/day) daily microbial recycling pool C (slow) */
	double daily_to_soil4C;                                               /* (gC/m2/day) daily recalcitrant SOM C (humus, slowest) */
	double daily_soil1C_to_soil2C;                                        /* (gC/m2/day) daily decomposition fluxes soil to soil */
	double daily_soil2C_to_soil3C;                                        /* (gC/m2/day) daily decomposition fluxes soil to soil */
	double daily_soil3C_to_soil4C;                                        /* (gC/m2/day) daily decomposition fluxes soil to soil */

	/* litter and soil nitrogen fluxes */
	double daily_N_miner;                                                 /*(gN/m2/day) daily nitrogen mineralized */
	double daily_N_pot_imm;                                               /*(gN/m2/day) daily nitrogen immobilized */
	/* litter */
	double daily_to_litrN;                                                /* (gN/m2/day) daily nitrogen to litter carbon pool at cell level */
	double daily_to_litr1N;                                               /* (gN/m2/day) daily nitrogen to litter labile carbon pool at cell level */
	double daily_to_litr2N;                                               /* (gN/m2/day) daily nitrogen to litter unshelded cellulose carbon pool at cell level */
	double daily_to_litr3N;                                               /* (gN/m2/day) daily nitrogen to litter shelded cellulose carbon pool at cell level */
	double daily_to_litr4N;                                               /* (gN/m2/day) daily nitrogen to litter lignin carbon pool at cell level */
	double daily_cwd_to_litrN;                                            /* (gN/m2/day) daily cwd to litter nitrogen pool at cell level */
	double daily_cwd_to_litr2N;                                           /* (gN/m2/day) daily decomposition fluxes cwd to litter */
	double daily_cwd_to_litr3N;                                           /* (gN/m2/day) daily decomposition fluxes cwd to litter */
	double daily_cwd_to_litr4N;                                           /* (gN/m2/day) daily decomposition fluxes cwd to litter */
	double daily_leaf_to_litrN;                                           /* (gN/m2/day) daily leaf to litter nitrogen pool at cell level */
	double daily_leaf_to_litr1N;                                          /* (gN/m2/day) daily leaf to litter labile nitrogen pool at cell level */
	double daily_leaf_to_litr2N;                                          /* (gN/m2/day) daily leaf to litter unshelded cellulose nitrogen pool at cell level */
	double daily_leaf_to_litr3N;                                          /* (gN/m2/day) daily leaf to litter  shelded cellulose nitrogen pool at cell level */
	double daily_leaf_to_litr4N;                                          /* (gN/m2/day) daily leaf to litter lignin nitrogen pool at cell level */
	double daily_froot_to_litrN;                                          /* (gN/m2/day) daily fine root to litter nitrogen pool at cell level */
	double daily_froot_to_litr1N;                                         /* (gN/m2/day) daily fine root to litter labile nitrogen pool at cell level */
	double daily_froot_to_litr2N;                                         /* (gN/m2/day) daily fine root to litter unshelded cellulose nitrogen pool at cell level */
	double daily_froot_to_litr3N;                                         /* (gN/m2/day) daily fine root to litter shelded cellulose nitrogen pool at cell level */
	double daily_froot_to_litr4N;                                         /* (gN/m2/day) daily fine root to litter lignin nitrogen pool at cell level */
	double daily_litr_to_soilN;                                           /* (gN/m2/day) daily litter nitrogen pool to soil pool */
	double daily_litr1N_to_soil1N;                                        /* (gN/m2/day) daily decomposition fluxes litter to soil */
	double daily_litr2N_to_soil2N;                                        /* (gN/m2/day) daily decomposition fluxes litter to soil */
	double daily_litr3N_to_litr2N;                                        /* (gN/m2/day) daily decomposition fluxes litter to soil */
	double daily_litr4N_to_soil3N;                                        /* (gN/m2/day) daily decomposition fluxes litter to soil */
	/* soil */
	double daily_to_soilN;                                                /* (gN/m2/day) daily litter nitrogen to soil */
	double daily_to_soil1N;                                               /* (gN/m2/day) daily balance of microbial recycling pool nitrogen (fast) */
	double daily_to_soil2N;                                               /* (gN/m2/day) daily balance of microbial recycling pool nitrogen (medium) */
	double daily_to_soil3N;                                               /* (gN/m2/day) daily balance of microbial recycling pool nitrogen (slow) */
	double daily_to_soil4N;                                               /* (gN/m2/day) daily balance of carbon to litter lignin carbon pool at cell level */
	double daily_soil1N_to_soil2N;                                        /* (gN/m2/day) daily decomposition fluxes soil to soil */
	double daily_soil2N_to_soil3N;                                        /* (gN/m2/day) daily decomposition fluxes soil to soil */
	double daily_soil3N_to_soil4N;                                        /* (gN/m2/day) daily decomposition fluxes soil to soil */
	double daily_soil4N_to_soilMinN;                                      /* (gN/m2/day) daily decomposition fluxes soil to mineralized */

	/* litter and soil carbon pools */
	double cwd_C;                                                         /* (gC/m2) cwd cwd carbon */
	double cwd_2C;                                                        /* (gC/m2) cwd cwd unshielded carbon */
	double cwd_3C;                                                        /* (gC/m2) cwd cwd shielded carbon */
	double cwd_4C;                                                        /* (gC/m2) cwd cwd lignin carbon */
	double litrC;                                                         /* (gC/m2) litter carbon */
	double litr1C;                                                        /* (gC/m2) litter labile carbon */
	double litr2C;                                                        /* (gC/m2) litter unshielded cellulose carbon */
	double litr3C;                                                        /* (gC/m2) litter shielded cellulose carbon */
	double litr4C;                                                        /* (gC/m2) litter lignin carbon */
	double leaf_litrC;                                                    /* (gC/m2) leaf total litter carbon */
	double leaf_litr1C;                                                   /* (gC/m2) leaf litter labile carbon */
	double reserve_litrC;                                                    /* (gC/m2) reserve total litter carbon */
	double reserve_litr1C;                                                   /* (gC/m2) reserve litter labile carbon */
	double leaf_litr2C;                                                   /* (gC/m2) leaf litter unshielded cellulose carbon */
	double leaf_litr3C;                                                   /* (gC/m2) leaf litter shielded cellulose carbon */
	double leaf_litr4C;                                                   /* (gC/m2) leaf litter lignin carbon */
	double froot_litrC;                                                   /* (gC/m2) fine root total litter carbon */
	double froot_litr1C;                                                  /* (gC/m2) fine root litter labile carbon */
	double froot_litr2C;                                                  /* (gC/m2) fine root litter unshielded cellulose carbon */
	double froot_litr3C;                                                  /* (gC/m2) fine root litter shielded cellulose carbon */
	double froot_litr4C;                                                  /* (gC/m2) fine root litter lignin carbon */
	double soilC;                                                         /* (gC/m2) Soil at cell level */
	double soil1C;                                                        /* (gC/m2) microbial recycling pool carbon (fast) */
	double soil2C;                                                        /* (gC/m2) microbial recycling pool carbon (medium) */
	double soil3C;                                                        /* (gC/m2) microbial recycling pool carbon (slow) */
	double soil4C;                                                        /* (gC/m2) recalcitrant SOM carbon (humus, slowest) */

	/* litter and soil nitrogen pools */
	double litrN;                                                         /* (gN/m2) litter nitrogen */
	double litr1N;                                                        /* (gN/m2) litter labile nitrogen */
	double litr2N;                                                        /* (gN/m2) litter unshielded cellulose nitrogen */
	double litr3N;                                                        /* (gN/m2) litter shielded cellulose nitrogen */
	double litr4N;                                                        /* (gN/m2) litter lignin nitrogen */
	double leaf_litrN;                                                    /* (gN/m2) leaf total litter nitrogen */
	double leaf_litr1N;                                                   /* (gN/m2) leaf litter labile nitrogen */
	double leaf_litr2N;                                                   /* (gN/m2) leaf litter unshielded cellulose nitrogen */
	double leaf_litr3N;                                                   /* (gN/m2) leaf litter shielded cellulose nitrogen */
	double leaf_litr4N;                                                   /* (gN/m2) leaf litter lignin nitrogen */
	double froot_litrN;                                                   /* (gN/m2) fine root total litter nitrogen */
	double froot_litr1N;                                                  /* (gN/m2) fine root litter labile nitrogen */
	double froot_litr2N;                                                  /* (gN/m2) fine root litter unshielded cellulose nitrogen */
	double froot_litr3N;                                                  /* (gN/m2) fine root litter shielded cellulonitrogenN */
	double froot_litr4N;                                                  /* (gN/m2) fine root litter lignin nitrogen */
	double cwd_N;                                                         /* (gN/m2) cwd nitrogen */
	double cwd_2N;                                                        /* (gC/m2) cwd unshielded nitrogen */
	double cwd_3N;                                                        /* (gC/m2) cwd shielded nitrogen */
	double cwd_4N;                                                        /* (gC/m2) cwd lignin nitrogen */
	double soilN;                                                         /* (gN/m2) soil nitrogen */
	double soil1N;                                                        /* (gN/m2) microbial recycling pool nitrogen (fast) */
	double soil2N;                                                        /* (gN/m2) microbial recycling pool nitrogen (medium) */
	double soil3N;                                                        /* (gN/m2) microbial recycling pool nitrogen (slow) */
	double soil4N;                                                        /* (gN/m2) recalcitrant SOM nitrogen (humus, slowest) */
	double SoilMinN;                                                      /* (gN/m2) Soil mineral nitrogen */
	double Npool;                                                         /* (gN/m2) temporary plant nitrogen pool */
	double Nfix;                                                          /* (gN/m2) SUM of biological nitrogen fixation */
	double Nleach;                                                        /* (gN/m2) SUM of nitrogen lost for leacheing */
	double Nvol;                                                          /* (gN/m2) SUM of nitrogen lost for volatilization  */
	double NPP_gN_demand;


	/* temporary nitrogen variables for reconciliation of decomposition
	immobilization fluxes and plant growth N demands */
	//fixme to used to remove
	double mineralized;
	double potential_immob;
	double plitr1c_loss;                                                  /* litter carbon labile decomposition */
	double pmnf_l1s1;
	double plitr2c_loss;                                                  /* litter carbon cellulose decomposition */
	double pmnf_l2s2;
	double plitr4c_loss;                                                  /* litter carbon ligning decomposition */
	double pmnf_l4s3;
	double psoil1c_loss;                                                  /* soil carbon labile decomposition */
	double pmnf_s1s2;
	double psoil2c_loss;                                                  /* soil carbon unshielded cellulose decomposition */
	double pmnf_s2s3;
	double psoil3c_loss;                                                  /* soil carbon shielded cellulose decomposition */
	double pmnf_s3s4;
	double psoil4c_loss;                                                  /* soil carbon lignin decomposition */
	double kl4;                                                           /* decay rate of lignin litter */

	/* heterotrophic soil */
	double daily_litr_het_resp;                                           /* (gC/m2/day) heterotrophic respiration for total litter carbon */
	double daily_litr1_het_resp;                                          /* (gC/m2/day) heterotrophic respiration for litter labile carbon */
	double daily_litr2_het_resp;                                          /* (gC/m2/day) heterotrophic respiration for litter unshielded cellulose carbon */
	double daily_litr4_het_resp;                                          /* (gC/m2/day) heterotrophic respiration for litter lignin carbon */
	double daily_soil_het_resp;                                           /* (gC/m2/day) heterotrophic respiration for total soil carbon */
	double daily_soil1_het_resp;                                          /* (gC/m2/day) heterotrophic respiration for microbial recycling pool carbon (fast) */
	double daily_soil2_het_resp;                                          /* (gC/m2/day) heterotrophic respiration for microbial recycling pool carbon (medium) */
	double daily_soil3_het_resp;                                          /* (gC/m2/day) heterotrophic respiration for microbial recycling pool carbon (slow) */
	double daily_soil4_het_resp;                                          /* (gC/m2/day) heterotrophic respiration for recalcitrant SOM carbon (humus, slowest) */

	/* soil fraction */
	double leaf_litt_scel_frac;                                           /* (DIM) leaf litter shielded cellulose fraction */
	double leaf_litt_uscel_frac;                                          /* (DIM) leaf litter unshielded cellulose fraction */
	double froot_litt_scel_frac;                                          /* (DIM) fine root litter shielded cellulose fraction */
	double froot_litt_uscel_frac;                                         /* (DIM) fine root litter unshielded fraction */
	double deadwood_scel_frac;                                            /* (DIM) dead wood litter shielded cellulose fraction */
	double deadwood_uscel_frac;                                           /* (DIM) dead wood litter unshielded fraction */

	/* water */
	int days_since_rain;                                                  /* (days) consecutive days without rain */

	/* pools */
	double swc;                                                           /* (kgH2O/day) soil Water content */
	double asw;                                                           /* (mm/volume) current available soil water  */
	double cum_asw;
	double max_asw_fc;                                                    /* (mmKgH2O/m3) max available soil water at field capacity */
	double max_asw_sat;                                                   /* (mmKgH2O/m3) max available soil water at field capacity */
	double snow_pack;                                                     /* (Kg/m2)current amount of snow */
	double canopy_water_stored;                                           /* (mm/m2) canopy water stored at cell level */
	double canopy_snow_stored;                                            /* (mm/m2) canopy snow stored at cell level */
	double old_water_store;                                               /* (mm/m2) previous day soil water stored */

	/* fluxes */
	double daily_snow_melt;                                               /* (mm/m2/time) current amount of melted snow  */
	double daily_snow_subl;                                               /* (mm/m2/time) current amount of sublimated snow */
	double daily_out_flow, monthly_out_flow, annual_out_flow;             /* (mm/m2/time) current amount of water outflow (runoff) */
	double daily_canopy_rain_int, monthly_canopy_rain_int, annual_canopy_rain_int; /* (mm/m2/time) daily, monthly and canopy rain interception */
	double daily_canopy_snow_int, monthly_canopy_snow_int, annual_canopy_snow_int; /* (mm/m2/time) daily, monthly and canopy snow interception */
	double daily_canopy_transp, monthly_canopy_transp, annual_canopy_transp; /* (mm/m2/time) daily, monthly and canopy transpiration */
	double daily_canopy_evapo, monthly_canopy_evapo, annual_canopy_evapo; /* (mm/m2/time) daily, monthly and canopy evaporation */
	double daily_canopy_et, monthly_canopy_et, annual_canopy_et;          /* (mm/m2/time) daily, monthly and canopy evapotranspiration */
	double daily_soil_evapo, monthly_soil_evapo, annual_soil_evapo;       /* (mm/m2/time) daily, monthly and soil evaporation */
	double daily_et, monthly_et, annual_et;                               /* (mm/m2/time) daily, monthly and evapotranspiration */
	double daily_tot_w_flux, monthly_tot_w_flux, annual_tot_w_flux;       /* (mm/m2/time) daily, monthly and annual water fluxes */
	double daily_soil_evapo_watt;                                         /* (W/m2) current daily soil evaporation */
	double daily_soil_lh_flux;                                            /* (W/m2) current daily soil latent heat flux */
	double daily_soil_sh_flux;                                            /* (W/m2) current daily soil sensible heat flux */
	double daily_soil_resp, monthly_soil_resp, annual_soil_resp;          /* (gC/m2/time) Soil respiration flux */
	double daily_soil_respCO2, monthly_soil_respCO2,annual_soil_respCO2;  /* (gCO2/m2/time) Soil respiration flux */
	double daily_r_eco, monthly_r_eco, annual_r_eco;                      /* (gC/m2/time) daily, monthly and annual ecosystem resp */
	double daily_het_resp, monthly_het_resp, annual_het_resp;             /* (gC/m2/time) daily, monthly and annual heterotrophic resp */

	double daily_irrigation;

	double daily_nep, monthly_nep, annual_nep;                            /* (gC/m2/time) daily, monthly and annual net ecosystem production */
	double daily_nee, monthly_nee, annual_nee;                            /* (gC/m2/time) daily, monthly and annual net ecosystem exchange */

	/************************************************** BALANCES VARIABLES **************************************************/
	/* cell */
	double cell_carbon_in;                                                /* (gC/day/day) cell carbon balances in */
	double cell_carbon_out;                                               /* (gC/day/day) cell carbon balances out */
	double cell_carbon_store;                                             /* (gC/day) cell carbon balances store */
	double cell_carbon_old_store;                                         /* (gC/day) cell carbon balances old_store */
	double cell_carbon_balance;                                           /* (gC/day/day) cell carbon balances balance */

	double cell_carbontc_in;                                              /* (tC/cell/day) cell carbon balances in */
	double cell_carbontc_out;                                             /* (tC/cell/day) cell carbon balances out */
	double cell_carbontc_store;                                           /* (tC/cell) cell carbon balances store */
	double cell_carbontc_old_store;                                       /* (tC/cell) cell carbon balances old_store */
	double cell_carbontc_balance;                                         /* (tC/cell/day) cell carbon balances balance */

	double cell_water_in;                                                 /* (mm/day/day) cell water balances in */
	double cell_water_out;                                                /* (mm/day/day) cell water balances out */
	double cell_water_store;                                              /* (mm/day) cell water balances store */
	double cell_water_old_store;                                          /* (mm/day) cell water balances old_store */
	double cell_water_balance;                                            /* (mm/day/day) cell water balances balance */

	double cell_snow_in;                                                  /* (mm/day/day) cell snow balances in */
	double cell_snow_out;                                                 /* (mm/day/day) cell snow balances out */
	double cell_snow_store;                                               /* (mm/day) cell snow balances store */
	double cell_snow_old_store;                                           /* (mm/day) cell snow balances old_store */
	double cell_snow_balance;                                             /* (mm/day/day) cell snow balances balance */

	double cell_nitrogen_in;                                              /* (gN/day/day) cell nitrogen balances in */
	double cell_nitrogen_out;                                             /* (gN/day/day) cell nitrogen balances out */
	double cell_nitrogen_store;                                           /* (gN/day) cell nitrogen balances store */
	double cell_nitrogen_old_store;                                       /* (gN/day) cell nitrogen balances old_store */
	double cell_nitrogen_balance;                                         /* (gN/day/day) cell nitrogen balances balance */

	double cell_nitrogentc_in;                                            /* (tN/cell/day) cell nitrogen balances in */
	double cell_nitrogentc_out;                                           /* (tN/cell/day) cell nitrogen balances out */
	double cell_nitrogentc_store;                                         /* (tN/cell) cell nitrogen balances store */
	double cell_nitrogentc_old_store;                                     /* (tN/cell/day) cell nitrogen balances old_store */
	double cell_nitrogentc_balance;                                       /* (tN/cell/day) cell nitrogen balances balance */

	/* litter */
	double litr_carbon_in;                                                /* (gC/day/day) cell carbon balances in */
	double litr_carbon_out;                                               /* (gC/day/day) cell carbon balances out */
	double litr_carbon_store;                                             /* (gC/day) cell carbon balances store */
	double litr_carbon_old_store;                                         /* (gC/day) cell carbon balances old_store */
	double litr_carbon_balance;                                           /* (gC/day/day) cell carbon balances balance */

	double litr_nitrogen_in;                                              /* (gN/day/day) cell nitrogen balances in */
	double litr_nitrogen_out;                                             /* (gN/day/day) cell nitrogen balances out */
	double litr_nitrogen_store;                                           /* (gN/day) cell nitrogen balances store */
	double litr_nitrogen_old_store;                                       /* (gN/day) cell nitrogen balances old_store */
	double litr_nitrogen_balance;                                         /* (gN/day/day) cell nitrogen balances balance */

	/* soil */
	double soil_carbon_in;                                                /* (gC/day/day) cell carbon balances in */
	double soil_carbon_out;                                               /* (gC/day/day) cell carbon balances out */
	double soil_carbon_store;                                             /* (gC/day) cell carbon balances store */
	double soil_carbon_old_store;                                         /* (gC/day) cell carbon balances old_store */
	double soil_carbon_balance;                                           /* (gC/day/day) cell carbon balances balance */

	double soil_nitrogen_in;                                              /* (gN/day/day) cell nitrogen balances in */
	double soil_nitrogen_out;                                             /* (gN/day/day) cell nitrogen balances out */
	double soil_nitrogen_store;                                           /* (gN/day) cell nitrogen balances store */
	double soil_nitrogen_old_store;                                       /* (gN/day) cell nitrogen balances old_store */
	double soil_nitrogen_balance;                                         /* (gN/day/day) cell nitrogen balances balance */

	/********************************************** TOPOGRAPHIC VARIABLES *******************************************************/

	/* TOPOGRAPHIC CELL INITIALIZATION */
	/* from "topo.txt" file */
	double elev;                                                          /* (m) elevation */

} cell_t;

typedef struct {
	cell_t *cells;
	int cells_count;
	int cells_avail;
	int x_cells_count;
	int y_cells_count;
} matrix_t;

matrix_t* matrix_create(const soil_settings_t*const s, const int count, const char* const filename, dataset_t** dataset);
void matrix_free           ( matrix_t *m);
void dataset_free          ( dataset_t *p);
void simulation_summary    ( const matrix_t* const m);
void site_summary          ( const matrix_t* const m);
void topo_summary          ( const matrix_t* const m);
void soil_summary          ( const matrix_t* const m, const cell_t* const cell);
void forest_initialization ( const matrix_t* const m, const int day, const int month, const int year);
void litter_initialization ( const matrix_t* const m, const int day, const int month, const int year);
void soil_initialization   ( const matrix_t* const m, const int day, const int month, const int year);
int layer_add              ( cell_t* const c);
int fill_species_from_file ( species_t *const s);
#endif /* MATRIX_H_ */
