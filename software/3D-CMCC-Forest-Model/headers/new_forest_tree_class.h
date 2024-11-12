/*
 * forest_tree_class.h
 *
 *  Created on: 18/ago/2016
 *      Author: alessio-cmcc
 */

#ifndef FOREST_TREE_CLASS_H_
#define FOREST_TREE_CLASS_H_

#include "matrix.h"

//int add_tree_class_for_replanting (cell_t *const c, const int day, const int month, const int year, const int rsi);
int add_tree_class_for_replanting (cell_t *const c, const int day, const int month, const int year);
int add_tree_class_for_regeneration ( cell_t *const c );

//static int fill_cell(cell_t *const c);

#endif /* FOREST_TREE_CLASS_H_ */
