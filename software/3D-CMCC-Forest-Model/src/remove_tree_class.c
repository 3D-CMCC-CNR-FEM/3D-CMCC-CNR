/*
 * remove_tree_class.c
 *
 *  Created on: 02/set/2016
 *      Author: alessio-cmcc
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "new_forest_tree_class.h"
#include "constants.h"
#include "settings.h"
#include "common.h"
#include "tree_model.h"
#include "remove_tree_class.h"

extern settings_t* g_settings;

static void species_free(species_t* s) {
	assert(s);

	if ( s->name ) {
		free(s->name);
		s->name = NULL;
	}
}

static void age_free(age_t* a) {
	int species;

	assert(a);

	for ( species = 0; species < a->species_count; ++species ) {
		species_free(&a->species[species]);
	}
	free(a->species);
	a->species = NULL;
}

static void dbh_free(dbh_t* d) {
	int age;
	int species;

	assert(d);

	for ( age = 0; age < d->ages_count; ++age ) {
		for ( species = 0; species < d->ages[age].species_count; ++species ) {
			if ( d->ages[age].species[species].name ) {
				free(d->ages[age].species[species].name);
			}
		}
		free(d->ages[age].species);
	}
	free(d->ages);
	d->ages = NULL;
}

static void height_free(height_t* h) {
	int dbh;

	assert(h);

	for ( dbh = 0; dbh < h->dbhs_count; ++dbh ) {
		dbh_free(&h->dbhs[dbh]);
	}
	free(h->dbhs);
	h->dbhs = NULL;
}

static int height_remove(cell_t *c, const int height) {
	int i;
	int y;

	assert(c);

	if ( height > c->heights_count ) {
		return 0;
	}

	y = 0;
	for ( i = 0; i < c->heights_count; ++i ) {
		if ( height == i ) {
			height_free(&c->heights[i]);
			continue;
		}
		c->heights[y++] = c->heights[i];
	}
	--c->heights_count;
	++c->heights_avail;

	return 1;
}

static int dbh_remove(cell_t *c, const int height, const int dbh) {
	int i;
	int y;

	assert(c);

	if ( height > c->heights_count ) {
		return 0;
	}

	if ( dbh > c->heights[height].dbhs_count ) {
		return 0;
	}

	y = 0;
	for ( i = 0; i < c->heights[height].dbhs_count; ++i ) {
		if ( dbh == i ) {
			dbh_free(&c->heights[height].dbhs[i]);
			continue;
		}
		c->heights[height].dbhs[y++] = c->heights[height].dbhs[i];
	}
	--c->heights[height].dbhs_count;
	++c->heights[height].dbhs_avail;

	return 1;
}

static int age_remove(cell_t *c, const int height, const int dbh, const int age) {
	int i;
	int y;

	assert(c);

	if ( height > c->heights_count ) {
		return 0;
	}

	if ( dbh > c->heights[height].dbhs_count ) {
		return 0;
	}

	if ( age > c->heights[height].dbhs[dbh].ages_count ) {
		return 0;
	}

	y = 0;
	for ( i = 0; i < c->heights[height].dbhs[dbh].ages_count; ++i ) {
		if ( age == i ) {
			age_free(&c->heights[height].dbhs[dbh].ages[i]);
			continue;
		}
		c->heights[height].dbhs[dbh].ages[y++] = c->heights[height].dbhs[dbh].ages[i];
	}
	--c->heights[height].dbhs[dbh].ages_count;
	++c->heights[height].dbhs[dbh].ages_avail;

	return 1;
}

static int species_remove(cell_t *c, const int height, const int dbh, const int age, const int species) {
	int i;
	int y;

	assert(c);

	if ( height > c->heights_count ) {
		return 0;
	}

	if ( dbh > c->heights[height].dbhs_count ) {
		return 0;
	}

	if ( age > c->heights[height].dbhs[dbh].ages_count ) {
		return 0;
	}

	if ( species > c->heights[height].dbhs[dbh].ages[age].species_count ) {
		return 0;
	}

	y = 0;
	for ( i = 0; i < c->heights[height].dbhs[dbh].ages[age].species_count; ++i ) {
		if ( species == i ) {
			species_free(&c->heights[height].dbhs[dbh].ages[age].species[i]);
			continue;
		}
		c->heights[height].dbhs[dbh].ages[age].species[y++] = c->heights[height].dbhs[dbh].ages[age].species[i];
	}
	--c->heights[height].dbhs[dbh].ages[age].species_count;
	++c->heights[height].dbhs[dbh].ages[age].species_avail;

	return 1;
}

static int species_zeroed(cell_t *c, const int height, const int dbh, const int age, const int species) {
	int i;

	assert(c);

	if ( height > c->heights_count ) {
		return 0;
	}

	if ( dbh > c->heights[height].dbhs_count ) {
		return 0;
	}

	if ( age > c->heights[height].dbhs[dbh].ages_count ) {
		return 0;
	}

	if ( species > c->heights[height].dbhs[dbh].ages[age].species_count ) {
		return 0;
	}

	for ( i = 0; i < c->heights[height].dbhs[dbh].ages[age].species_count; ++i ) {
		if ( species == i ) {
			int y;

			for ( y = 0; y < VALUES; y++ ) {
				c->heights[height].dbhs[dbh].ages[age].species[i].value[y] = 0.;
			}
			for ( y = 0; y < COUNTERS; y++ ) {
				c->heights[height].dbhs[dbh].ages[age].species[i].counter[y] = 0;
			}
			//c->heights[height].dbhs[dbh].ages[age].species[i].management = 0;
			c->heights[height].dbhs[dbh].ages[age].species[i].phenology_phase = 0;
		}
	}

	return 1;
}

int tree_class_remove(cell_t *const c, const int height, const int dbh, const int age, const int species)
{
	assert(c);

	if ( c->heights_count ) {
		if ( height > c->heights_count ) return 0;
		if ( c->heights[height].dbhs_count ) {
			if (  dbh > c->heights[height].dbhs_count ) return 0;
			if ( c->heights[height].dbhs[dbh].ages_count ) {
				if ( age > c->heights[height].dbhs[dbh].ages_count ) return 0;
				if ( c->heights[height].dbhs[dbh].ages[age].species_count ) {
					if ( species > c->heights[height].dbhs[dbh].ages[age].species_count ) return 0;
				}
			}
		}
	}

#if 1
	/* remove class if N_TREE < 0 or if called by harvesting function */
	if ( ! c->heights[height].dbhs[dbh].ages[age].species[species].counter[N_TREE]
			|| g_settings->management ) {

		if ( ! species_remove(c, height, dbh, age, species) ) return 0;

		if ( 1 == c->heights[height].dbhs[dbh].ages_count ) {
			if ( c->heights[height].dbhs[dbh].ages[0].species_count <= 1 ) {
				if ( ! age_remove(c, height, dbh, 0) ) return 0;
			}
		}

		if ( 1 == c->heights[height].dbhs_count ) {
			if ( ! c->heights[height].dbhs[0].ages_count
				|| ( (1 == c->heights[height].dbhs[0].ages_count)
					&& (c->heights[height].dbhs[0].ages[0].species_count <= 1) ) ) {
					if ( ! dbh_remove(c, height, 0) ) return 0;
				}
		}
#if 0
		if ( 1 == c->heights_count ) {
			if ( ! c->heights[0].dbhs_count 
				|| ( (1 == c->heights[0].dbhs_count) 
					&& (c->heights[0].dbhs[0].ages_count <= 1)
					&& (c->heights[0].dbhs[0].ages[0].species_count <= 1) )
				) {
					if ( ! height_remove(c, 0) ) return 0;
			}
		}
#else
		if ( ! c->heights[height].dbhs_count 
			|| ( (1 == c->heights[height].dbhs_count) 
				&& (c->heights[height].dbhs[0].ages_count <= 1)
				&& (c->heights[height].dbhs[0].ages[0].species_count <= 1) )
			) {
				if ( ! height_remove(c, height) ) return 0;
		}
#endif
	}
#else
	/* zeroed species class if N_TREE < 0 or if called by harvesting function */
	if ( ! c->heights[height].dbhs[dbh].ages[age].species[species].counter[N_TREE]
			/*|| g_settings->management*/ ) {	
		if ( ! species_zeroed(c, height, dbh, age, species) ) return 0;
	}
#endif

	return 1;
}
