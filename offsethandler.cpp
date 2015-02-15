#include"profiler.h"

FILE *foffset;
static char RNAME[300];
static int posn;
static int offset;
static int offsetn;
static int *totaltrueindel;

int handler_nextval();

int handler_init(FILE *_foffset, int *_totaltrueindel){
    foffset = _foffset;
    totaltrueindel = _totaltrueindel;
    handler_nextseq();
    return 0;
}

int handler_nextseq(){
    char temp[300];
    char *d;
    int i;

    do{
        handler_nextval();
    }while(posn != INT32_MAX);

    do{
        if(fgets(temp, 300, foffset) == NULL){
            return EOF;
        }
    }while(temp[0] != '>');

    d = strpbrk(temp, "\n");
    if(d == NULL){
        fprintf(stderr, "offsethandler: error in offset file\n"); 
        exit(1);
    }
    i = (int) (d-temp-1);
    strncpy(RNAME, temp+1, i);
    RNAME[i]=0;
   
    handler_nextval();
    offset = 0;

    return 0;
}

int handler_nextval(){
    char temp[40];
    char c;
    int i;
    c = fgetc(foffset);
    if((c == '>') || (c == EOF)){
        ungetc(c, foffset);
        posn = INT32_MAX;
        offsetn = 0;
    }else{
        ungetc(c, foffset);
        fgets(temp, 40, foffset);
        
        offset = offsetn;

        if(sscanf(temp, "%d %d", &posn, &offsetn) < 2){
            fprintf(stderr, "offsethandler: error in offset file\n"); 
            exit(1);
        }

        i = offsetn - offset; 
        if(i > 0){//insertion
            totaltrueindel[(i>10)?0:i]++;
        }else if(i < 0){//deletion
            totaltrueindel[(i<-10)?0:(10-i)]++;
        }else{
            fprintf(stderr, "offsethandler: error in offset file\n"); 
            exit(1);
        }
    }
    return 0;
}

int32_t posr2t(int32_t posr, const char *name){
    while(strcmp(name, RNAME)){
        handler_nextseq();
    }

    while(posr + offsetn >= posn){
        handler_nextval();
    }

    return MIN(posr + offset, posn);
}
