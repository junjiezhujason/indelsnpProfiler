#ifndef VCOUNT_H
#define VCOUNT_H

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<fstream>
using namespace std;

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

#endif
