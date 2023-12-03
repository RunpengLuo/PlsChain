#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#include <string.h>

#include "index.h"
#include "query.h"


int main(int argc, char *argv[]){
    int i_flag = 0;
    int q_flag = 0;
    int k_flag = 0;
    int o_flag = 0;
    int h_flag = 0;

    int err_flag = 0;

    size_t k = 0;
    char * i_dir = NULL; // index directory
    char * o_dir = NULL; // output directory

    char * qry_file = NULL;

    char ** comp_files = NULL; // plasmid components
    int num_comp = 0;

    int c;

    opterr = 0;
    while ((c = getopt (argc, argv, "ik:q:o:h")) != -1)
    switch (c)
    {
        case 'k':
            if (q_flag || h_flag) err_flag ++;
            k_flag = 1;
            k = atoi(optarg);
            break;
        case 'i':
            if (q_flag || h_flag) err_flag ++;
            i_flag = 1;
            break;
        case 'q':
            if (i_flag || k_flag || h_flag) err_flag ++;
            q_flag = 1;
            i_dir = optarg;
            break;
        case 'o':
            o_flag = 1;
            o_dir = optarg;
            break;
        case 'h':
            if (i_flag || k_flag || q_flag) err_flag ++;
            h_flag = 1;
            break;
        // case ':':       /* -k, -q or -o without operand */
        //     fprintf(stderr, "[ERROR] Option -%c requires an operand\n", optopt);
        //     return 1;
        case '?':
            if (isprint (optopt))
                fprintf(stderr, "[ERROR] Parsing error on option -%c, See `./plschain -h` for usage\n", optopt);
            else
                fprintf(stderr, "[ERROR] Unknown option character `\\x%x', See `./plschain -h` for usage\n", optopt);
            return 1;
        default:
            abort();
    }
    if (h_flag & !err_flag) {
        fprintf(stdout, "Usage: plschain -i -k INT -o DIRECTORY FILE1 FILE2 FILE3 ...\n");
        fprintf(stdout, "       plschain -q DIRECTORY -o DIRECTORY <query.fa>\n");
        fprintf(stdout, "Options:\n");
        fprintf(stdout, "    -i            Indexing mode\n");
        fprintf(stdout, "    -q DIRECTORY  Query mode, index directory\n");
        fprintf(stdout, "    -k INT        k-mer size [15,32]\n");
        fprintf(stdout, "    -o DIRECTORY  output directory\n");
        fprintf(stdout, "    -h            show this message\n");
        return 0;
    }
    if (i_flag == 0 && q_flag == 0) {
        fprintf(stderr, "[ERROR] No mode, See `./plschain -h` for usage\n");
        return 1;
    }
    if (o_flag == 0) {
        fprintf(stderr, "[ERROR] No output DIRECTORY specified\n");
        return 1;
    }

    if (err_flag){
        fprintf(stderr, "[ERROR] Parsing error, See `./plschain -h` for usage\n");
        return 1;
    }

    if (o_dir[strlen(o_dir)-1] == '/'){
            o_dir[strlen(o_dir)-1] = '\0';
    }
	if (mkdir(o_dir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1){
		fprintf(stderr, "%s\n", strerror(errno));
        return 1;
	}

    struct timeval start, end;
    double elapsed;
    gettimeofday(&start, NULL);
    int ret = 0;
    if (i_flag == 1){
        if (k_flag == 0) {
            fprintf(stderr, "[ERROR] No k-mer size specified for index (-i) mode\n");
            return 1;
        }
        if (optind == argc) {
            fprintf(stderr, "[ERROR] No Plasmid components are specified for index (-i) mode\n");
            return 1;
        }
        comp_files = &argv[optind];
        num_comp = argc - optind;
        // run index
        ret = indexing(k, num_comp, comp_files, o_dir);
    } else {
        /* q_flag == 1 */
        if (optind == argc) {
            fprintf(stderr, "[ERROR] No query file is specified for query (-q) mode\n");
            return 1;
        }
        if (optind != argc - 1) {
            // handle multiple query file later.
            fprintf(stderr, "[ERROR] Multiple query files are supplied\n");
            return 1;
        }

        if (i_dir[strlen(i_dir)-1] == '/'){
                i_dir[strlen(i_dir)-1] = '\0';
        }

        qry_file = argv[optind];
        ret = query_file(i_dir, o_dir, qry_file);
    }

    gettimeofday(&end, NULL);
    elapsed = end.tv_sec + end.tv_usec / 1e6 -
            start.tv_sec - start.tv_usec / 1e6; // in seconds

    fprintf(stdout, "[%s] Time elapsed: %fsec\n", ret == 0 ? "SUCCESS" : "FAILURE", elapsed);
    return ret;
}
