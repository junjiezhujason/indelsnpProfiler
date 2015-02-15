#include"profiler.h"


int seqreader::enlarge_buffer(){
    range(pos_start, pos_end);
    i = n; // i = old n
    n *= 2; // double n
    free(buffer);
    buffer = (char*) malloc(n*sizeof(char)); 
    if(buffer == NULL){
        fprintf(stderr, "seqreader::enlarge_buffer insufficient memory\n");
        exit(1);
    }
    memcpy(buffer, temp, i*sizeof(char));
    free(temp);
    temp = (char*) malloc((n+1)*sizeof(char));
    if(temp == NULL){
        fprintf(stderr, "seqreader::enlarge_buffer insufficient memory\n");
        exit(1);
    }

    int32_t j;
    char c;
    for(j = 0; j < i; j++){
        c = getbase();
        if(c == EOF) break;
        buffer[j+i] = c;
    }
    pos_end += j;

    i = 0; // new i
    return 0;
}

char* seqreader::range(int32_t start, int32_t end){
    if(start < pos_start){
        advance(pos_start-start);
    }
    if(start < pos_start){
        fprintf(stderr, "range som ting wong %d %d\n", start, pos_start);
        exit(1);
    }
    if(start <= end){
        while(!eos()){
            if(end - pos_start >= n){
                enlarge_buffer();
            }
            if(end - pos_start < n) break;
        }
        end = MIN(end, pos_end);
        int a = (int)(start - pos_start) + i;
        a = a % n;
        int b = (int)(end - pos_start) + i;
        b = b % n;
        if(a <= b){
            memcpy(temp, buffer+a, (b-a+1)*sizeof(char));
        }
        else{
            memcpy(temp, buffer+a, (n-a)*sizeof(char));
            memcpy(temp+n-a, buffer, (b+1)*sizeof(char));
        }
        temp[(b-a+n)%n+1] = '\0';
        return temp;
    }else{
//  start > end : nothing
        temp[0] = '\0';
        return temp;
    }
}

int seqreader::advance(int32_t d){
    int32_t j;
    char c;
    for(j = 0;j < d; j++){
        c = getbase();
        if(c == EOF) break;
        buffer[(j+i)%n] = c;
    }
    i = (i+j)%n;
    pos_start += j;
    pos_end += j;
    return 0;
}

int seqreader::advanceto(int32_t pos){
    return advance(pos - pos_start);
}

char seqreader::getbase(){
    char c;
    do{
        c = fgetc(fa);
        if(c == '\n'){
            c = fgetc(fa);
            if(c == '>'){
                ungetc('>', fa);
                _eos = true;
                return EOF;
            }
            ungetc(c, fa);
            c = '\n';
            continue;
        }else if(c == EOF){
            _eos = true;
            return EOF;
        }else{
            return c;
        }
    }while(c != EOF);
        fprintf(stderr, "getbase som ting wong\n");
    return 0;
}



int seqreader::nextseq(){
    char c;
    char *d;
    bool cont = true;
    int j = 0;
    while(cont){
        c = fgetc(fa);
        if(c == '>'){
            break;
        }else if(c == EOF){
            return EOF;
        }else{
            fprintf(stderr, "nextseq som ting wong\n");
        }
    }
    while(cont){
        c = fgetc(fa);
        if(c == '\n'){
            name[j] = 0;
            break;           
        }else if(c == EOF){
            return EOF;
        }else{
            name[j] = c;
            j++;
        }
    }
    d = strpbrk(name, " \t\n");
    if(d != NULL) *d = '\0';
    else name[299] = '\0';
    _eos = false;

    i = 0;
    pos_start = 0;

    for(pos_end = 0; pos_end < n; pos_end++){
        c = getbase();
        if(c == EOF) break;
        buffer[pos_end] = c;
    }
    pos_end--;
    return 0;
}

bool seqreader::eof(){
    return feof(fa);
}

char seqreader::at(int32_t pos){
    if(pos < pos_start){
        fprintf(stderr, "at som ting wong\n");
        return '\0';
    }
    while(pos > pos_end){
        if(eos()) return '\0';
        enlarge_buffer();
    }
    return buffer[(pos-pos_start+i)%n];
}
