#ifndef PLS_CHAIN_FILE_IO_H
#define PLS_CHAIN_FILE_IO_H

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

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

// assert op_dir doesn't end with '/', and op_dir is created already
FILE * open_file(char * op_dir, char * op_file, const char * __restrict op_mode) {
    int slen = strlen(op_dir) + strlen(op_file) + 3;
    char filename[slen];
    snprintf(filename, slen, "%s/%s", op_dir, op_file);
    filename[slen-1] = '\0';
    FILE * op_fd = fopen(filename, op_mode);
    if(op_fd == NULL) {
        fprintf(stderr, "Not able to open the file %s with mode %s\n", filename, op_mode);
    }
    return op_fd;
}


#endif