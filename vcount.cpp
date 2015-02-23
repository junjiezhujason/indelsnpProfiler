#include"vcount.h"

VarCount::VarCount(string name)
{
    int i;
    snp = 0;
    for(i = 0; i < 20; i++){
        indel[i] = 0;
    }
    for(i = 0; i < 9; i++){
        lindel[i] = 0;
    }
    filename = name;
    if (filename != ""){     
        file.open(filename.c_str());
        //cout << "File Open: " << name << endl;
    }
}

VarCount::~VarCount(void)
{ 
    if (filename != ""){
        file.close();
        //cout << "File Closed: " << filename << endl;
    }
}

void VarCount::update(int pos, char *allele0, char *allele1)
{ 
    int varlen;
    varlen = strlen(allele1) - strlen(allele0);
    if(varlen == 0){ //snps
        snp++;
    }else if(varlen > 0){ //insertion
        if(varlen <= 10){ 
            indel[varlen-1]++; // short 
        }else if(varlen <= 30){
            lindel[1]++;
        }else if(varlen <= 50){
            lindel[2]++;
        }else if(varlen <= 70){
            lindel[3]++;
        }else if(varlen <= 90){
            lindel[4]++;
        }else{
            lindel[0]++;
        }
    }else{ //deletion
        if(varlen >= -10){
            indel[9-varlen]++; // short
        }else if(varlen >= -30){
            lindel[5]++;
        }else if(varlen >= -50){
            lindel[6]++;
        }else if(varlen >= -70){
            lindel[7]++;
        }else if(varlen >= -90){
            lindel[8]++;
        }else{
            lindel[0]++;
        }
    }
    if (filename != ""){     
        file << pos << "\t" << allele0 << "\t" << allele1 << endl;
    }
}