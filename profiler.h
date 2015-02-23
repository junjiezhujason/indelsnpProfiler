#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"vcf.h"
#include"hts.h"
#include"simpleRepeat.h"
#include<fstream>
using namespace std;

#define MULTI_SEQ 0

#define MAX(a,b) (a)>(b)?(a):(b)
#define MIN(a,b) (a)<(b)?(a):(b)

#define inext() {\
    ir = bcf_read(indelBCF, ih, ib);\
    if(ir == 0){\
        bcf_unpack(ib, BCF_UN_STR);\
    }else{\
        ib->rid = nseq;\
        ieof = true;\
    }\
}
#define tnext() {\
    tr = bcf_read(testBCF, th, tb);\
    if(tr == 0){\
        bcf_unpack(tb, BCF_UN_STR);\
    }else{\
        tb->rid = nseq;\
        tb->qual = qual_threshold;\
        teof = true;\
    }\
}

class VarCount
{
public:
    int snp;
    int indel[20];
    int lindel[9];
    ofstream file;
    string  filename;
    VarCount(string name="");  // constructor declaration
    ~VarCount();               // destructor declaration
    void update(int pos, char *allele0, char *allele1);
};

int simple_distance(const char *word1, const char *word2); // distance.cpp
int length_distance(const char *word1, const char *word2);
int edit_distance(const char *word1, const char *word2);

