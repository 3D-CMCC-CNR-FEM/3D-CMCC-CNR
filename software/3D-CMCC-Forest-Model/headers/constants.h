/* constants.h */
#ifndef CONSTANTS_H_
#define CONSTANTS_H_

enum {
	JANUARY
	, FEBRUARY
	, MARCH
	, APRIL
	, MAY
	, JUNE
	, JULY
	, AUGUST
	, SEPTEMBER
	, OCTOBER
	, NOVEMBER
	, DECEMBER

	, MONTHS_COUNT
};

#define NO_DATA                -9999
/** constants **/
#define Q0                    1376.1      /* Solar constant (W/m^2) */
#define Q0_MJ                  0.082      /* Solar constant (MJ/m^2/day) */
#define QA                       -60      /* Intercept of Net vs Solar Radiation (W/m^2) - that is the net long wave radiation, here set as average value -60 instead of -90 */
#define QB                       0.8      /* Slope of Net vs Solar Radiation - that is more or less 1-ALBEDO*/ 
#define RAD2PAR                 0.45      /* (DIM) ratio PAR / SWtotal */
#define RAD2PPFD                0.45      /* (DIM) ratio PPFD / SWtotal */
#define GC_GDM                     2      /* convert grams C to grams DM (Campioli et al., 2013) */
#define GLUCOSE_C                0.4      /* fraction of Carbon in glucose reserve (Barbaroux et al., 2002) */
#define MOLPAR_MJ                2.3      /* Conversion of Solar Radiation to PAR - CONVERT MJ TO MOLPAR -- 1 MJ = 2.3 molPAR */
#define GDM_MOL                   24      /* Molecular weight of dry matter - CONVERT molC TO gDM (12 Molecular weight of C * 2 to have DM)(IPCC guidelines 1996) */
#define GC_MOL                12.011      /* Convert molC to grams of C */
#define GCO2_MOL                  44.     /* molecular weight of CO2, g mole-1 */
#define MJ_TO_W         11.574074074      /* Convert MJ/m2/day to W/m2 */
#define W_TO_MJ               0.0864      /* Convert W/m2 to MJ/m2/day */
#define W_MJ                 1000000      /* Convert Watt to MegaJoule/sec m^2 */
#define SBC_W              5.6696e-8      /* (W/(m^2 K^4)) Stefan-Boltzmann constant */
#define SBC_MJ              4.903e-9      /* (MJ/K^4/m^2/day) Stefan-Boltzmann constant in MJ day */
#define GCO2_MOL                  44.     /* molecular weight of CO2, g mole-1 */
#define KARM                    0.41      /* von Karman's constant (dimensionless) */
#define CP                    1010.0      /* (J/kg K) specific heat of air */
#define G_STD                9.80665      /* (m/s2) standard gravitational accel. */
#define P_STD               101325.0      /* (Pa) standard pressure at 0.0 m elevation */
#define T_STD                 288.15      /* (K) standard temp at 0.0 m elevation  */
#define AIRMASS           28.9644e-3      /* (kg/mol) molecular weight of air */
#define WATERMASS         18.0148e-3      /* (kg/mol) molecular weight of water */
#define LR_STD                0.0065      /* (-K/m) standard temperature lapse rate */
#define Rgas                  8.3144      /* (m3 Pa/ mol K) gas law constant */
#define EPS                   0.6219      /* (MW/MA) unitless ratio of molec weights */
#define ESTAR                0.61076      /* to compute saturation vapour pressure at the air temperature (KPa) */
#define PPFD50                  75.0      /* (umol/m2/s) PPFD for 1/2 stomatal closure */
#define PARD50                  6.48      /* (mol/m2/day) PARD for 1/2 stomatal closure */
#define EPAR                    4.55      /* (umol/J OR mol/MJ) PAR photon energy ratio */
#define SLAMKG                   0.1	  /* Convert SLA IN cm ^2/gC in m^2/KgC */
#define E20                      2.2      /* rate of change of saturated VPD with T at 20C */
#define GCtoGW                  1.57      /* conversion from conductance to CO2 to H2O */
#define VPDCONV             0.000622      /* convert VPD to saturation deficit = 18/29/1000 */
#define EVAPOCOEFF              1.32      /* Priestley Taylor Coefficient (Hobbins et al., 2001) */
#define MWratio                0.622      /* ratio molecular weight of water vapour air */
#define TempAbs               273.13      /* absolute temperature in Kelvin */
#define GDD_BASIS                  5      /* temperature C° basis for thermic sum */
#define EMLEAF                  0.95      /* Emissivity of thermal radiation by leaf (see MAESPA) */
#define EMCANOPY               0.987      /* Emissivity of thermal radiation by canopy (see Campbell & Norman, Environmental Biophysics, 2nd edition, page 273) */
#define EMSOIL                  0.95      /* Emissivity of thermal radiation by soil (see MAESPA) */
#define Pi               3.141592654      /* Pi greek */
#define LN2              0.693147181      /* logarithm constant */
#define eps3                    1e-3      /* differences in comparison */
#define eps                     1e-6      /* differences in comparison */
#define eps_pr                  1e-15      /* critical precision */
#define ZERO                  -1e-10      /* power limit of state variables before they are set to 0.0 to control rounding and overflow errors */
#define ZERO_C                 -1e-8      /* power limit of state variables of C cycle at cell level  TBD  if to force to 0*/
#define CRIT_PREC              1e-20      /* precision control */
#define e                 2.71828183      /* nepero/eulero constant */
#define DBH_ref                  1.3      /* reference DBH value */
#define GRPERC                   0.3      /* (DIM) percentage of growth resp per unit of C grown see BIOME, 0.25 for LPJ */
#define GRPERCMIN               0.25      /* (DIM) MINIMUM percentage of growth resp per unit of C grown see Ryan et al., */
#define GRPERCMAX               0.30      /* (DIM) MAXIMUM percentage of growth resp per unit of C grown see Ryan et al., */
#define O2CONC                20.946      /* % of [O2] see Verustraete 1994, 2002 */
#define DBH_MAX                  200      /* Maximum stem diameter at breast height (in cm) */
#define STEMPOWER_A	            2.08      /* Power in the Stem Mass v. diameter relationship for DBH < 9 cm */
#define STEMPOWER_B	            2.64      /* Power in the Stem Mass v. diameter relationship for 9 < DBH < 15 cm */
#define STEMPOWER_C	            2.51      /* Power in the Stem Mass v. diameter relationship for DBH > 15 cm */
#define MAX_N_TREE_LAYER           3      /* define the maximum number for tree layers (todo move into setting.txt) */
#define C_FRAC_TO_RETRANSL       0.1      /* fraction of C to re-translocate (see Bossell et al., 2006 and Campioli et al., 2013; 0.2 Hoffmann 1995 */
#define N_FRAC_TO_RETRANSL       0.1      /* fraction of N to re-translocate (see Bossell et al., 2006 and Campioli et al., 2013; 0.2 Hoffmann 1995 */
#define MAX_PEAK_LAI_PROJ          6      /* Maximum Peak Lai (m2/m2) */
#define WATER_STRESS_LIMIT       0.3      /* Minimum imposed water stress limit to soil water modifier */
#define SNOW_ABS                 0.6      /* Absorptivity of snow */


/* soil constants (BIOME-BGC method) */
#define BULK_DENITRIF_FRAC       0.5      /* Fraction of Nitrogen lost through Denitrification when N demand is lower than N available */
#define SOIL1_CN                12.0      /* C:N for fast microbial recycling pool */
#define SOIL2_CN                12.0      /* C:N for slow microbial recycling pool */
#define SOIL3_CN                10.0      /* C:N for recalcitrant SOM pool (humus) */
#define SOIL4_CN                10.0      /* C:N for recalcitrant SOM pool (humus) */
#define DENITRIF_PROPORTION     0.01      /* fraction of mineralization to volatile */
#define MOBILEN_PROPORTION       0.1      /* fraction mineral N avail for leaching */

/* use this block of constants to include the dynamics for slowest soil pool (s4) */
/* respiration fractions for fluxes between compartments (unitless) */
#define	RFL1S1                  0.39      /* transfer from litter 1 to soil 1 */
#define	RFL2S2                  0.55      /* transfer from litter 2 to soil 2 */
#define	RFL4S3                  0.29      /* transfer from litter 4 to soil 3 */
#define	RFS1S2                  0.28      /* transfer from soil 1 to soil 2 */
#define	RFS2S3                  0.46      /* transfer from soil 2 to soil 3 */
#define	RFS3S4                  0.55      /* transfer from soil 3 to soil 4 */

/* base (maximum) decomposition rate constants (1/day) */
#define KL1_BASE                0.7       /* labile litter pool */
#define KL2_BASE               0.07       /* cellulose litter pool */
#define KL4_BASE              0.014       /* lignin litter pool */
#define KS1_BASE               0.07       /* fast microbial recycling pool */
#define KS2_BASE              0.014       /* fast microbial recycling pool */
#define KS3_BASE             0.0014       /* slow microbial recycling pool */
#define KS4_BASE             0.0001       /* recalcitrant SOM (humus) pool */
#define KFRAG_BASE            0.001       /* physical fragmentation of coarse woody debris */

// thresholds/constant  5p7
#define thinning_threshold      35        /* [year] age limit to perform thinning when shelterwood is performed */
#define coppice_rf              0.10       /* reduction factor of stem density-replanting after each clear cut for coppice */
#define reg_threshold      30        /* [year] age threshold of the youngest forest class to set seedl_layer to 1 (regeneration layer exists) */

//MARCONI CONSTANTS::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//soil constants
/*
#define DRF                     0.04      //0.03
#define EFFRB                    .67
#define EFFAC                     .2
#define RBO                     0.02
#define SRB                      0.9
#define KRCVL                   0.25
#define KRCL                   0.074
#define KRCR                    0.02
#define KRB                     0.12
#define HRB                     0.04
#define EFFNO                    0.2
//#define um_no3   .67
//#define um_no2   .67
//#define um_no    .34
//#define um_n2o   .34

 */

#endif /* CONSTANTS_H */
