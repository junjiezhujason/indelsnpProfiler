#ifndef SIMPLEREPEAT_H
#define SIMPLEREPEAT_H
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

typedef struct{
    char *fn;
    FILE *fp;
    char **chrom;
    unsigned int nchrom;
} srFILE_t;

typedef struct{
    unsigned int bin;
    unsigned int chromID;
    unsigned int chromStart;
    unsigned int chromEnd;
    char name[256];
    unsigned int period;
    double copyNum;
    unsigned int consensusSize;
    unsigned int perMatch;
    unsigned int perIndel;
    unsigned int score;
    unsigned int A;
    unsigned int C;
    unsigned int G;
    unsigned int T;
    double entropy;
    char *sequence;
}sr1_t;

srFILE_t *sr_open(const char *fn, const char *mode); //read only support

sr1_t *sr_init1();

int sr_close(srFILE_t *fp);

int sr_read1(srFILE_t *fp, sr1_t *s);
int sr_write1(srFILE_t *fp, sr1_t *s);

#endif
