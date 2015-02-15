#include"simpleRepeat.h"

srFILE_t *sr_open(const char *fn, const char *mode){
    srFILE_t *fp = (srFILE_t*) malloc(sizeof(srFILE_t));
    fp->fp = fopen(fn, mode);
    if(fp == NULL) return NULL;
    fp->fn = (char*) malloc(sizeof(char)*strlen(fn));
    strcpy(fp->fn, fn);
    fp->nchrom = 0;
    return fp;
}

int sr_close(srFILE_t *fp){
    return fclose(fp->fp);
}

void nexttab(FILE *fp, char *buffer){
    int i;
    char c;

    i = 0;
    c = fgetc(fp);
    while((c != '\t') && (c != '\n') && (c != EOF)){
        buffer[i] = c;
        i++;
        c = fgetc(fp);
    }
    buffer[i] ='\0';
}

sr1_t *sr_init1(){
    return (sr1_t*) calloc(1, sizeof(sr1_t));
}

int sr_read1(srFILE_t *fp, sr1_t *s){
    static char buffer[256];
    
    nexttab(fp->fp, buffer); if(feof(fp->fp)) return EOF;   s->bin = atoi(buffer); 

    nexttab(fp->fp, buffer); if(feof(fp->fp)) return EOF;
    if(fp->nchrom && (strcmp(fp->chrom[fp->nchrom-1], buffer) == 0)){
        s->chromID = fp->nchrom-1;
    }else{
        fp->nchrom++;
        fp->chrom = (char**) realloc(fp->chrom, fp->nchrom*sizeof(char*));
        fp->chrom[fp->nchrom-1] = (char*) malloc((strlen(buffer)+1)*sizeof(char));
        strcpy(fp->chrom[fp->nchrom-1], buffer);
        s->chromID = fp->nchrom-1;
    }

    nexttab(fp->fp, buffer); if(feof(fp->fp)) return EOF;   s->chromStart = atoi(buffer) - 1; //change to 0 base
    nexttab(fp->fp, buffer); if(feof(fp->fp)) return EOF;   s->chromEnd = atoi(buffer) - 1; //change to 0 base
    nexttab(fp->fp, buffer); if(feof(fp->fp)) return EOF;   strcpy(s->name, buffer);
    nexttab(fp->fp, buffer); if(feof(fp->fp)) return EOF;   s->period = atoi(buffer);
    nexttab(fp->fp, buffer); if(feof(fp->fp)) return EOF;   s->copyNum = atof(buffer); 
    nexttab(fp->fp, buffer); if(feof(fp->fp)) return EOF;   s->consensusSize = atoi(buffer);
    nexttab(fp->fp, buffer); if(feof(fp->fp)) return EOF;   s->perMatch = atoi(buffer);
    nexttab(fp->fp, buffer); if(feof(fp->fp)) return EOF;   s->perIndel = atoi(buffer);
    nexttab(fp->fp, buffer); if(feof(fp->fp)) return EOF;   s->score = atoi(buffer);
    nexttab(fp->fp, buffer); if(feof(fp->fp)) return EOF;   s->A = atoi(buffer);
    nexttab(fp->fp, buffer); if(feof(fp->fp)) return EOF;   s->C = atoi(buffer);
    nexttab(fp->fp, buffer); if(feof(fp->fp)) return EOF;   s->G = atoi(buffer);
    nexttab(fp->fp, buffer); if(feof(fp->fp)) return EOF;   s->T = atoi(buffer);
    nexttab(fp->fp, buffer); if(feof(fp->fp)) return EOF;   s->entropy = atof(buffer);

    s->sequence = (char*) realloc(s->sequence, (s->consensusSize+1)*sizeof(char));
    nexttab(fp->fp,s->sequence); if(strlen(s->sequence) != s->consensusSize) return EOF;

    


    return 0;
}

int sr_write1(srFILE_t *fp, sr1_t *s){
    static char n[2] = {'\0', '\0'};
    fprintf(fp->fp, "%s%u\t%s\t%u\t%u\t%s\t%u\t%lf\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%lf\t%s",
      n, s->bin, fp->chrom[s->chromID], s->chromStart+1, s->chromEnd+1, s->name, s->period, s->copyNum, s->consensusSize,
      s->perMatch, s->perIndel, s->score, s->A, s->C, s->G, s->T, s->entropy, s->sequence);
    n[0] = '\n';
    return 0;
}
