#ifndef PLS_CHAIN_TREE_H
#define PLS_CHAIN_TREE_H
#include "dtype.h"
#include <stdio.h>

typedef struct _leaf_t {
    tetra_t tetra;
    int lb;
    int rb;
    struct _leaf_t * parent;
    struct _leaf_t ** children;
    size_t l, m;
} leaf_t;

void init_leaf(leaf_t * leaf, tetra_t tetra, int lb, int rb) {
    leaf->tetra = tetra;
    leaf->lb = lb;
    leaf->rb = rb;

    leaf->parent = NULL;
    leaf->children = NULL;
    leaf->l = 0;
    leaf->m = 0;
}

typedef struct _tree_t {
    leaf_t root;
    leaf_t ** leaves;
    size_t l, m;
} tree_t;

tetra_t null_tetra(){
    return (tetra_t) {-1, -1, -1, -1};
}

void init_tree(tree_t * tree) {
    init_leaf(&(tree->root), null_tetra(), -1, -1);
    tree->m = 1;
    tree->leaves = (leaf_t **) malloc(sizeof(leaf_t *));
    tree->l = 1;
    tree->leaves[0] = &tree->root;
}

leaf_t * append_leaves(leaf_t *** leaves, int * m, int * l, tetra_t tetra, int lb, int rb){
    if (*m == 0) {
        *m = ARR_SIZE;
        *leaves = (leaf_t **) calloc(*m, sizeof(leaf_t *));
    }
    if (*m == *l){
        *m += STEP_SIZE;
        *leaves = (leaf_t **) realloc(*leaves, *m * sizeof(leaf_t *));
    }
    leaf_t * appended = (leaf_t *) malloc(sizeof(leaf_t));
    init_leaf(appended, tetra, lb, rb);
    (*leaves)[*l] = appended;
    *l += 1;
    return appended;
}

void append_children(leaf_t * node, leaf_t * child){
    if (node->m == 0) {
        node->m = ARR_SIZE;
        node->children = (leaf_t **) calloc(node->m, sizeof(leaf_t *));
    }
    if (node->m == node->l){
        node->m += STEP_SIZE;
        node->children = (leaf_t **) realloc(node->children, node->m * sizeof(leaf_t *));
    }
    node->children[node->l] = child;
    node->l ++;
}

int append_tree(tree_t * tree, layer_t * layer) {
    // if (tree->l == 0) { // invalid tree
    //     return false;
    // }

    int m = 0, l = 0;
    leaf_t ** new_leaves = NULL;
    leaf_t * old_leaf;
    leaf_t * new_leaf;
    int prev_s, next_s;

    for (int i = 0; i < layer->l; i ++) {
        for (int j = 0; j < tree->l; j ++) {
            new_leaf = NULL;
            old_leaf = tree->leaves[j];
            prev_s = old_leaf->tetra.q_pos;
            next_s = layer->tetras[i].q_pos;
            if (prev_s == -1) {
                // previous layer is root only
                new_leaf = append_leaves(&new_leaves, &m, &l, layer->tetras[i], -1, next_s);
            } else {
                // previous layer is intermediate layer
                if ((next_s > old_leaf->lb) && (next_s < old_leaf->rb)){
                    new_leaf = append_leaves(&new_leaves, &m, &l, layer->tetras[i], next_s, old_leaf->rb);
                } else if ((old_leaf->lb == -1) && (next_s > prev_s)) {
                    new_leaf = append_leaves(&new_leaves, &m, &l, layer->tetras[i], -1, old_leaf->rb);
                }
            }
            if (new_leaf != NULL) {
                new_leaf->parent = old_leaf;
                append_children(old_leaf, new_leaf);
            }
        }
    }

    if (l == 0) {
        // no new leaves be found
        // append NULL leaves to previous layer to emulate the missing part, 
        // inhert info from previous layer such as qpos, lb, rb
        for (int j = 0; j < tree->l; j ++) {
            old_leaf = tree->leaves[j];
            new_leaf = append_leaves(&new_leaves, &m, &l, null_tetra(), old_leaf->lb, old_leaf->rb);
            new_leaf->tetra.q_pos = old_leaf->tetra.q_pos;
            new_leaf->parent = old_leaf;
            append_children(old_leaf, new_leaf);
        }
    }

    if (tree->m < m) {
        tree->leaves = (leaf_t **) realloc(tree->leaves, m*sizeof(leaf_t *));
        tree->m = m;
    }
    memset(tree->leaves, 0, tree->m * sizeof(leaf_t *));
    tree->l = l;
    memcpy(tree->leaves, new_leaves, tree->l * sizeof(leaf_t *));
    free(new_leaves);
    return true;
}

tetra_t * tree_traverse(leaf_t * leaf, int * m, int * l, int * w){
    *m = 0;
    *l = 0;
    *w = 0;
    tetra_t * path = NULL;

    leaf_t * node = leaf;
    while (node->parent != NULL) {
        if (*m == 0) {
            *m = ARR_SIZE;
            path = (tetra_t *) calloc(*m, sizeof(tetra_t));
        }
        if (*m == *l) {
            *m += STEP_SIZE;
            path = (tetra_t *) realloc(path, *m * sizeof(tetra_t));
        }
        path[*l] = node->tetra;
        *w += node->tetra.count;
        *l += 1;
        node = node->parent;
    }
    return path;
}

tetra_t * get_max_path(tree_t * tree, int *max_m, int * max_l, int * max_w, int is_asc) {
    *max_m = 0, *max_l = 0, *max_w = 0;
    tetra_t * max_path = NULL;

    int m, l, w;
    tetra_t * path = NULL;
    for (int i = 0; i < tree->l; i ++){
        path = tree_traverse(tree->leaves[i], &m, &l, &w);
        if (l > 0) {
            if (*max_l == 0) {
                *max_l = l;
                *max_w = w;
                *max_m = m;
                max_path = path;
            } else {
                if (w > *max_w) {
                    free(max_path);
                    max_path = NULL;
                    *max_l = l;
                    *max_w = w;
                    *max_m = m;
                    max_path = path;
                } else {
                    free(path);
                }
            }
        }
    }
    if (is_asc){
        // reverse max_path
        tetra_t tmp_tetra;
        int pivot = *max_l / 2;
        for (int i = 0; i < pivot; i ++){
            tmp_tetra = max_path[i];
            max_path[i] = max_path[*max_l - i - 1];
            max_path[*max_l - i - 1] = tmp_tetra;
        }
    }

    return max_path;
}

void free_leaf(leaf_t * node){
    if (node == NULL) return;
    for (int i = 0; i < node->l; i ++) {
        free_leaf(node->children[i]);
    }
    if (node->l != 0){
        free(node->children);
    }
    free(node);
}

void free_tree(tree_t * tree){
    if (tree == NULL) return;
    for (int i = 0; i < tree->root.l; i ++){
        free_leaf(tree->root.children[i]);
    }
    if (tree->root.l != 0){
        free(tree->root.children);
    }

    if (tree->l != 0) {
        free(tree->leaves);
        tree->leaves = NULL;
    }
    free(tree);
}



#endif