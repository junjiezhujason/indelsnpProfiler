#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
using namespace std;

// WORKING PROPTYPE!

class VarCount
{
public:
    int snp;
    int indel[20];
    int lindel[9];
    ofstream file;
    string  filename;
    VarCount(string name="");  // constructor declaration
    ~VarCount(); // destructor declaration

    void update(int pos, char *allele0, char *allele1);

};

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
        file.open(filename);
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


int main(){
    VarCount totaltrue;                              // false negatives + true positives
    VarCount concord;                                // true positive, concordance (true calls)
    VarCount discord("catagory/discord.txt");        // discordance (false calls in new positions)
    VarCount pconcord("catagory/pconcord.txt");      // partial concordance (false calls in correct positions
    VarCount nonconcord("catagory/nonconcord.txt");  // false negatives (discordance + partial concordance)


    char a0[] = "A";
    char a1[] = "ATAG";
    char **allele ;
    allele[0] = a0;
    allele[1] = a1;

    discord.update(1, allele[0], allele[1]);
    discord.update(6, allele[0], allele[1]);
    discord.update(7, allele[0], allele[1]);

    
    return 0;
}