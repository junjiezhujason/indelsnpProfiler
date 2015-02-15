#include <iostream>



int * updatecounts (int varlen, int p_type[3]){
    // * options; totaltrue, concord, discord, pconcord, nonconcord; 
    int temp;
    static int *newcounts;
    newcounts = p_type;
    if(varlen == 0){ //snps
        newcounts[0]++;
            }else if(varlen > 0){ //insertion
                if(varlen <= 10){ // short
                    newcounts[1]++;
                }
            }else{ //deletion
                if(varlen >= -10){ // short
                    newcounts[2]++;
                }
            }
    return newcounts;
}

int main(){

    int concord[3];    // true positive, concordance (true calls)
    int discord[3];    // discordance (false calls in new positions)
    int pconcord[3];   // partial concordance (false calls in correct positions
    int nonconcord[3]; // false negatives (discordance + partial concordance)


    int totaltrue[3];  // false negatives + true positives

    int i;

    // initialize array
    for(i = 0; i < 3; i++){
        totaltrue[i] = 0;
    }

    for(i = 0; i < 3; i++){
        printf("%d\n",totaltrue[i]);
    }
    int *t;
    t = updatecounts(0, totaltrue);
    totaltrue = &t[0];
    for(i = 0; i < 3; i++){
        printf("%d\n",t[i]);
    }

    

    return 0;
}