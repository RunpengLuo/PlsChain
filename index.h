#ifndef PLS_CLASS_INDEX_H
#define PLS_CLASS_INDEX_H
#include "kseq.h"
#include "khash.h"
#include "kmer.h"
#include "dtype.h"

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
	int slen;

    // build hash table
	if (k < 15 || k > 32) {
		fprintf(stderr, "k-mer size should be in [15,32]");
		return 1;
	}
	ukmer_table = kh_init(0);

	slen = strlen(out_dir) + strlen("comps_xyz.txt") + 3;
	char out_file[slen];

	snprintf(out_file, slen, "%s/%s", out_dir, "comps_idx.txt");
	out_file[slen-1] = '\0';
	fd_idx = fopen(out_file, "wb");
	if(fd_idx == NULL) {
        fprintf(stderr, "Not able to open the output file %s\n", out_file);
        return 1;
    }
	fprintf(fd_idx, "%d\n", num_comp);
	for (int i = 0; i < num_comp; i ++){
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

		kseq_destroy(seq);
		gzclose(fp);
	}

    out_file[0] = '\0';
	snprintf(out_file, slen, "%s/%s", out_dir, "comps_udb.txt");
	out_file[slen-1] = '\0';
	fd_udb = fopen(out_file, "wb");
	if(fd_udb == NULL) {
        fprintf(stderr, "Not able to open the output file %s\n", out_file);
        return 1;
    }

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
	fclose(fd_udb);
	kh_destroy(0, ukmer_table);

	return 0;
}
#endif