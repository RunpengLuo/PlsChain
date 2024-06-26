#ifndef PLS_CHAIN_INDEX_H
#define PLS_CHAIN_INDEX_H
#include "kseq.h"
#include "khash.h"
#include "kmer.h"
#include "dtype.h"
#include "file_io.h"

#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <fcntl.h>

KHASH_MAP_INIT_INT64(0, tuple_t)

int kmer_profile(kseq_t * seq, size_t k, int comp_idx, int s_idx, khash_t(0) *ukmer_table){
	int i, len;
	khint64_t key, ckey, mask;
	
	khiter_t kh;
	int absent;

	if (seq->seq.l < k) {
		fprintf(stderr, "invalid k=%zu selection", k);
		return 1;
	}

	len = 0;
	key = 0;
	mask = get_mask(k);

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
			ckey = MIN(key, revComp(key, mask, k));
			kh = kh_put(0, ukmer_table, ckey, &absent);
			if (absent) {
				kh_val(ukmer_table, kh) = (tuple_t) {comp_idx, s_idx};
			} else {
				kh_val(ukmer_table, kh) = (tuple_t) {-1, -1}; // not unique anymore
			}
            len--;
        }
	}
	return 0;
}

int indexing(size_t k, int num_comp, char ** comp_files, char * out_dir) {
    int ret;
	int comp_idx, s_idx;

	// seq parsing
	gzFile fp;
	kseq_t *seq;
	int l;

	// hash table
	khash_t(0) *ukmer_table;
	khiter_t kh;

	// output
	FILE *fd_idx;
	FILE *fd_udb;
	FILE *fd_nme;
	FILE *fd_cfg;

    // build hash table
	if (k < 15 || k > 32) {
		fprintf(stderr, "k-mer size should be in [15,32]");
		return 1;
	}
	ukmer_table = kh_init(0);

	if ((fd_idx = open_file(out_dir, "comps_idx.txt", "wb")) == NULL) {
		return 1;
	}

	if ((fd_nme = open_file(out_dir, "comps_nme.txt", "wb")) == NULL) {
		return 1;
	}

	if ((fd_udb = open_file(out_dir, "comps_udb.txt", "wb")) == NULL) {
		return 1;
	}

	if ((fd_cfg = open_file(out_dir, "comps_cfg.txt", "wb")) == NULL) {
		return 1;
	}

	// record the max component length
	int per_comp = 0;
	int total_comp = 0;

	fprintf(fd_idx, "%d\n", num_comp);
	for (int i = 0; i < num_comp; i ++){
		fprintf(fd_nme, "%s\n", comp_files[i]);
		comp_idx = i;
		s_idx = 0;
		/* parsing the file */
		fp = gzopen(comp_files[i], "r");
		if (fp == Z_NULL) {
			fprintf(stderr, "Not able to open the component file %s\n", comp_files[i]);
        	kh_destroy(0, ukmer_table);
			return 1;
		}
		seq = kseq_init(fp);
		while ((l = kseq_read(seq)) >= 0) {
			per_comp = seq->seq.l > per_comp ? seq->seq.l : per_comp;
			fprintf(fd_idx, "%s\n", seq->name.s);
			ret = kmer_profile(seq, k, comp_idx, s_idx, ukmer_table);
			if (ret == 1){
				fprintf(stderr, "invalid component sequence, size less than k %zu\n", k);
				kh_destroy(0, ukmer_table);
				kseq_destroy(seq);
				gzclose(fp);
				return 1;
			}
			s_idx ++;
		}
		fprintf(fd_idx, "\n");
		total_comp += per_comp;
		per_comp = 0;

		kseq_destroy(seq);
		gzclose(fp);
	}

	fprintf(fd_cfg, "%d\n", total_comp);

    fprintf(fd_udb, "%zu\n", k);
	for (kh = 0; kh < kh_end(ukmer_table); kh ++) {
		tuple_t value;
		khint64_t key;
		if (!kh_exist(ukmer_table, kh)) continue;
		value = kh_val(ukmer_table, kh);
		if (value.comp_idx != -1) {
			// unique k-mer
			key = kh_key(ukmer_table, kh);
			fprintf(fd_udb, "%lu\t%d\t%d\n", key, value.comp_idx, value.s_idx);
		}
	}

	// clean up
	fclose(fd_idx);
	fclose(fd_nme);
	fclose(fd_udb);
	fclose(fd_cfg);
	kh_destroy(0, ukmer_table);

	return 0;
}
#endif