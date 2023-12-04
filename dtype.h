#ifndef PLS_CHAIN_DTYPE_H
#define PLS_CHAIN_DTYPE_H
#include <zlib.h>

#include "khash.h"

#define MIN(a,b)        (((a)<(b))?(a):(b))
#define MAX(a,b)        (((a)>(b))?(a):(b))
#define true 1
#define false 0

#define ARR_SIZE 5
#define STEP_SIZE 3

typedef struct _tuple_t {
	int comp_idx;
	int s_idx;
} tuple_t;

typedef struct _item_t {
	tuple_t value;
	int count;
} item_t;

void init_item(item_t * item, int comp_idx, int s_idx, int count){
	item->count = count;
	item->value.comp_idx = comp_idx;
	item->value.s_idx = s_idx;
}

int ieq_item(item_t * item, tuple_t value){
	if (item->value.comp_idx == value.comp_idx && item->value.s_idx == value.s_idx) {
		item->count ++;
		return true;
	} else {
		return false;
	}
}

typedef struct _tetra_t {
	int ia;
	int ja;
	int a_jdx;
	int count;
} tetra_t;

void init_tetra(tetra_t * tetra, int ia, int ja, int a_jdx, int count){
	tetra->ia = ia;
	tetra->ja = ja;
	tetra->a_jdx = a_jdx;
	tetra->count = count;
}

typedef struct _layer_t {
	size_t l, m;
	tetra_t * tetras;
} layer_t;

void append_layer(layer_t * layer, item_t item, int a_jdx){
	if (layer->m == 0) {
		layer->m = ARR_SIZE;
		layer->tetras = (tetra_t *) calloc(layer->m, sizeof(tetra_t));
	}
	if (layer->m == layer->l) {
		layer->m += STEP_SIZE;
		layer->tetras = (tetra_t *) realloc(layer->tetras, layer->m *sizeof(tetra_t));
	}
	layer->tetras[layer->l] = (tetra_t) {item.value.comp_idx, item.value.s_idx, a_jdx, item.count};
	layer->l ++;
}

void free_layer(layer_t * layers, int num_comp){
	for (int i = 0; i < num_comp; i ++){
		if (layers[i].tetras != NULL) {
			free(layers[i].tetras);
		}
	}
	free(layers);
}

khint64_t get_mask(size_t k){
    khint64_t mask = 0;
	for (int i = 0; i < k; i ++) mask = ((mask << 2) | 0b11);
    return mask;
}

#endif