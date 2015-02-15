#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>

#define CASE_INSENSITIVE 1

int simple_distance(const char *str1, const char *str2){
#if CASE_INSENSITIVE
    int i;
    for(i = 0; str1[i]; i++){
        if(toupper(str1[i]) != toupper(str2[i])) return 1;
    }
    if(str2[i]) return 1;
    return 0;
#else
    return strcmp(str1,str2) == 0 ? 0 : 1;
#endif
}

int length_distance(const char *str1, const char *str2){
    int x = strlen(str1) - strlen(str2);
    return x > 0 ? x : -x;
}

#define matrix(a,b) _matrix[((size_t)(a))*(len2+1)+(b)]
#define MIN(a,b,c) ((a)<(b))?(((a)<(c))?(a):(c)):(((b)<(c))?(b):(c))

int edit_distance(const char *word1, const char *word2){
    int len1, len2, i, j;
    int n0, n1;
    int k;
    int *_matrix = NULL;
    char *upperword1, *upperword2;
    char c1, c2;
    size_t nbytes;

    len1 = strlen(word1);
    len2 = strlen(word2);
    k = len1 + 1;

    nbytes = sizeof(int);
    nbytes *= k;
    nbytes *= len2+1;
    _matrix = (int*) malloc(nbytes);

    while(_matrix == NULL){
        k = k/2;
        nbytes = sizeof(int);
        nbytes *= k;
        nbytes *= len2+1;
        _matrix = (int*) malloc(nbytes);
    }
    if(k < 2){
        fprintf(stderr, "distance out of memory\n");
        exit(1);
    }

#if CASE_INSENSITIVE
    upperword1 = (char*) malloc(sizeof(char)*len1);
    upperword2 = (char*) malloc(sizeof(char)*len2);
#else
    upperword1 = word1;
    upperword2 = word2;
#endif

#if CASE_INSENSITIVE
    for (i = 0; i < len1; i++) {
        upperword1[i] = toupper(word1[i]);
    }
#endif

    for (i = 0; i < len2; i++) {
#if CASE_INSENSITIVE
        upperword2[i] = toupper(word2[i]);
#endif
        matrix(0, i) = i;
    }
    matrix(0, len2) = len2;

    matrix(0, 0) = 0;

    for (i = 1 ; i <= len1; i++) {
        n0 = (i)%k;
        matrix(n0, 0) = i;
        n1 = (i-1)%k;
        c1 = upperword1[i-1];
        for (j = 1; j <= len2; j++) {
            c2 = upperword2[j-1];
            if (c1 == c2) matrix(n0, j) = matrix(n1, j-1);
            else matrix(n0, j) = MIN(matrix(n1, j), matrix(n0, j-1), matrix(n1, j-1)) + 1;
            //else matrix(n0, j) = MIN(matrix(n1, j) + 1, matrix(n0, j-1) + 1, matrix(n1, j-1) + 1);
            //MIN(DEL,INS,SUB)
        }
    }
#if CASE_INSENSITIVE
    free(upperword1);
    free(upperword2);
#endif
    i = matrix((len1)%k, len2);
    free(_matrix);
    if(len1 == 499944) printf("distance end\n");
    return i;

}

