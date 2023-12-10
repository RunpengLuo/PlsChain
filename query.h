#ifndef PLS_CHAIN_QUERY_H
#define PLS_CHAIN_QUERY_H
#include "kseq.h"
#include "khash.h"
#include "kmer.h"
#include "dtype.h"
#include "tree.h"
#include "file_io.h"

#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <fcntl.h>

KHASH_MAP_INIT_INT64(1, tuple_t)

char *** read_idx_file(char * db_dir, int * num_comp, int ** sizes);
void free_idx_arrs(char *** idx_arrs, int * sizes, int num_comp);
int read_udb_file(char * db_dir, size_t * k, khash_t(1) * ukmer_table);
layer_t * proc_query(kseq_t * seq, size_t k, int num_comp, khash_t(1) * ukmer_table, char *** idx_arrs, int * status);

int query_file(char * db_dir, char * out_dir, char * query_file) {
    size_t k;

	// seq parsing
	gzFile fp;
	kseq_t *seq;
	int l;

	// hash table
	khash_t(1) *ukmer_table;

    // idx array
    char *** idx_arrs;
    int * sizes;
    int num_comp;

	// output
	FILE *fd_qry_total;
	int status;

	if ((fd_qry_total = open_file(out_dir, "qry_total.csv", "wb")) == NULL) {
		return 1;
	}

    if ((idx_arrs = read_idx_file(db_dir, &num_comp, &sizes)) == NULL) {
        fclose(fd_qry_total);
        return 1;
    }

    ukmer_table = kh_init(1);
    // error handling
    if ((status = read_udb_file(db_dir, &k, ukmer_table)) == 1){
        fclose(fd_qry_total);
        kh_destroy(1, ukmer_table);
        free_idx_arrs(idx_arrs, sizes, num_comp);
        return 1;
    }

    // process query
    if ((fp = gzopen(query_file, "r")) == Z_NULL) {
        fclose(fd_qry_total);
        fprintf(stderr, "Not able to open the component file %s\n", query_file);
        kh_destroy(1, ukmer_table);
        free_idx_arrs(idx_arrs, sizes, num_comp);
        return 1;
    }

    seq = kseq_init(fp);
    layer_t * res = NULL;
    int alpha, beta;
    int rcount = 0, rclassified = 0;
    while ((l = kseq_read(seq)) >= 0) {
        rcount ++;
	res = proc_query(seq, k, num_comp, ukmer_table, idx_arrs, &status);
        fprintf(fd_qry_total, "%s", seq->name.s);
        if (status != 1) {
            fprintf(fd_qry_total, ",fail\n");
        } else {
	        rclassified ++;
            for (int i = 0; i < res->l; i ++){
                alpha = (res->tetras)[i].alpha;
                beta = (res->tetras)[i].beta;

                fprintf(fd_qry_total, ",%s", idx_arrs[alpha][beta]);
            }
            fprintf(fd_qry_total, "\n");
            free(res->tetras);
            free(res);
        }
    }

    // clean up
    fclose(fd_qry_total);

    kh_destroy(1, ukmer_table);
    free_idx_arrs(idx_arrs, sizes, num_comp);

    kseq_destroy(seq);
    gzclose(fp);

    printf("Process %d reads, with %d reads been classified\n", rcount, rclassified);
    return 0;
}

void free_idx_arrs(char *** idx_arrs, int * sizes, int num_comp){
    int size, i, j;
    for (i = 0; i < num_comp; i ++){
        size = sizes[i];
        for (j = 0; j < size; j ++){
            free(idx_arrs[i][j]);
            idx_arrs[i][j] = NULL;
        }
        free(idx_arrs[i]);
        idx_arrs[i] = NULL;
    }
    free(idx_arrs);
    free(sizes);
}

char *** read_idx_file(char * db_dir, int * num_comp, int ** sizes){
    int i, j, s;
    char *** idx_arrs = NULL;

    FILE *fd_idx;
    char *line = NULL;
    size_t len = 0;
    ssize_t read;

	if ((fd_idx = open_file(db_dir, "comps_idx.txt", "rb")) == NULL) {
		return NULL;
	}

    i = 0;
    read = getline(&line, &len, fd_idx); // get first line
    line[read - 1] = '\0';
    *num_comp = atoi(line);
    *sizes = (int *) malloc(sizeof(int) * *num_comp);

    idx_arrs = (char ***) calloc(*num_comp, sizeof(char **));

    s = ARR_SIZE;
    idx_arrs[i] = (char **) calloc(s, sizeof(char *));
    (*sizes)[i] = 0;

    j = 0;
    while ((read = getline(&line, &len, fd_idx)) != -1) {
        if (read == 1){ //read a \n
            if (i < *num_comp - 1){
                i ++;
                s = ARR_SIZE;
                idx_arrs[i] = (char **) calloc(s, sizeof(char *));
                (*sizes)[i] = 0;
                j=0;
            }
            continue;
        }

        if (j == s) {
            s += STEP_SIZE;
            idx_arrs[i] = (char **) realloc(idx_arrs[i], s * sizeof(char *));
        }

        idx_arrs[i][j] = (char *) malloc(read * sizeof(char));
        memcpy(idx_arrs[i][j], line, (read - 1) * sizeof(char));
        idx_arrs[i][j][read - 1] = '\0';
        (*sizes)[i] ++;
        j++;
    }

    free(line);
    fclose(fd_idx);
    return idx_arrs;
}

int read_udb_file(char * db_dir, size_t * k, khash_t(1) * ukmer_table){
    FILE *fd_udb;
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    khint64_t mask;
    khiter_t kh;

    int absent;
    char *token;

	if ((fd_udb = open_file(db_dir, "comps_udb.txt", "rb")) == NULL) {
		return 1;
	}

    read = getline(&line, &len, fd_udb); // get first line
    line[read - 1] = '\0';
    *k = atoi(line);
    mask = get_mask(*k);

    while ((read = getline(&line, &len, fd_udb)) != -1) {
        line[read - 1] = '\0';

        token = strtok(line, "\t");
        khint64_t key = (khint64_t) strtoul(token, NULL, 10);
        token = strtok(NULL, "\t");

        tuple_t value;
        value.comp_idx = atoi(token);
        token = strtok(NULL, "\t");
        value.s_idx = atoi(token);

        kh = kh_put(1, ukmer_table, key, &absent);
        kh_val(ukmer_table, kh) = value;

        kh = kh_put(1, ukmer_table, revComp(key, mask, *k), &absent);
        kh_val(ukmer_table, kh) = value;
    }
    if (line != NULL) free(line);
    fclose(fd_udb);
    return 0;
}

// status = 2 if error, 1 if success, 0 if fail
layer_t * proc_query(kseq_t * seq, size_t k, int num_comp, khash_t(1) * ukmer_table, char *** idx_arrs, int * status){
    *status = 0;
    int i, len, q_pos;
	khint64_t key, mask;
	
	khiter_t kh;

	if (seq->seq.l < k) {
		// fprintf(stderr, "invalid k=%zu selection on sname %s\n", k, seq->name.s);
        *status = 2;
        return NULL;
	}

	len = 0;
	key = 0;
	mask = get_mask(k);

    int arr_s = ARR_SIZE;
    int a_idx = 0;
    item_t * arraylist = (item_t *) malloc(sizeof(item_t) * arr_s);
    init_item(&arraylist[a_idx], -1, -1, 0);
    tuple_t kmer_idd;

    // parse our kmers
	for (i = 0; i < (int) seq->seq.l; i ++){
		if (!(seq->seq.s[i] == 'A' || seq->seq.s[i] == 'C' || seq->seq.s[i] == 'G' || seq->seq.s[i] == 'T'))
        {
            key = 0;
            len = 0;
            continue;
        }
		key = (key << 0b10);
		key = (key | ((seq->seq.s[i] >> 0b1) & 0b11));
		key = mask & key;
        len++;
		if (len == k)
        {
            // key
            kh = kh_get(1, ukmer_table, key);
            if (kh != kh_end(ukmer_table)){
                kmer_idd = kh_val(ukmer_table, kh);
                if (!ieq_item(&arraylist[a_idx], kmer_idd)) {
                    if (arraylist[a_idx].count > 0) {
                        a_idx ++;
                        if (a_idx == arr_s) {
                            arr_s += STEP_SIZE;
                            arraylist = (item_t *) realloc(arraylist, arr_s * sizeof(item_t));
                        }
                    }
                    init_item(&arraylist[a_idx], kmer_idd.comp_idx, kmer_idd.s_idx, 1);
                }
            }
            len--;
        }
	}

    if (a_idx == 0) {
        // no unique k-mer be found
        // fprintf(stderr, "no unique-kmer for %s\n", seq->name.s);
        free(arraylist);
        *status = 0;
        return NULL;
    }

    item_t item;
    layer_t * layers = (layer_t *) calloc(num_comp, sizeof(layer_t));
    for (i = 0; i < num_comp; i ++){
        layers[i].l = 0;
        layers[i].m = 0;
        layers[i].tetras = NULL;
    }

    int alpha;
    for (q_pos = 0; q_pos < a_idx + 1; q_pos ++) {
        item = arraylist[q_pos];
        alpha = item.value.comp_idx;
        append_layer(&layers[alpha], item, q_pos);
    }
    free(arraylist);

    tree_t * chain_tree_asc = (tree_t *) malloc(sizeof(tree_t));
    tree_t * chain_tree_dsc = (tree_t *) malloc(sizeof(tree_t));
    init_tree(chain_tree_asc);
    init_tree(chain_tree_dsc);

    for (i = 0; i < num_comp; i ++) {
        append_tree(chain_tree_asc, &layers[i]);
        append_tree(chain_tree_dsc, &layers[num_comp - i - 1]);
    }
    free_layer(layers, num_comp);

    int w0, l0, m0;
    int w1, l1, m1;

    tetra_t * path0 = get_max_path(chain_tree_asc, &m0, &l0, &w0, true);
    tetra_t * path1 = get_max_path(chain_tree_dsc, &m1, &l1, &w1, false);
    free_tree(chain_tree_asc);
    free_tree(chain_tree_dsc);

    tetra_t * max_path;
    int max_l, max_m;
    if (w0 >= w1){
        if (path1 != NULL) free(path1);

        max_path = path0;
        max_l = l0;
        max_m = m0;
    } else {
        if (path0 != NULL) free(path0);

        max_path = path1;
        max_l = l1;
        max_m = m1;
    }

    layer_t * rtn = NULL;
    if (max_l != 0){
        *status = 1;
        rtn = (layer_t *) malloc(sizeof(layer_t));
        rtn->l = max_l;
        rtn->m = max_m;
        rtn->tetras = max_path;
    } else {
        *status = 0;
        // fprintf(stdout, "skip %s\n", seq->name.s);
    }
	return rtn;
}
#endif