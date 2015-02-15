#include"profiler.h"

//assume chromosome order is the same for both files!!


int main(int argc, char **argv){

    htsFile *indelBCF, *testBCF;

    bcf1_t *ib, *tb;
    bcf_hdr_t *ih, *th;

    int nseq;
    int ir, tr;
    int totaltruesnps, calledtruesnps, calledfalsesnps;
    bool ieof, teof;
    int indellen;
    int dALT;

    float qual_threshold;
    int approx_threshold;
    int i;
    int totaltrueindel[20]; //0-9 ins //10-19 del
    int calledtrueindel[20]; //0-9 ins //10-19 del
    int calledfalseindel[20]; //0-9 ins //10-19 del
    int ltotaltrueindel[9]; //1 11-30 ins //2 31-50 ins //3 51-70 ins //4 71-90 ins
    int lcalledtrueindel[9]; //5 11-30 del //6 31-50 del //7 51-70 del //8 71-90 del
    int lcalledfalseindel[9]; //0 91+ indel

    char *d;

    int (*distance)(const char *, const char *);

    if(argc != 5){
        fprintf(stderr, "Usage: profiler <indel.bcf> <test.bcf> <qual_threshold> <approx_threshold>\n");
        exit(1); 
    }

    indelBCF = bcf_open(argv[1], "r");
    if(indelBCF == NULL){
        fprintf(stderr, "Unable to open %s\n", argv[1]);
        fprintf(stderr, "Usage: profiler <indel.bcf> <test.bcf> <qual_threshold> <approx_threshold>\n");
        exit(1);
    }

    d = (char*) malloc((strlen(argv[2])+5)*sizeof(char));
        sprintf(d, "%s.bcf", argv[2]);
        testBCF = bcf_open(d, "r");
    if(testBCF == NULL){
        d[strlen(d)-3] = 'v';// sprintf(d, "%s.vcf", argv[2];
        testBCF = vcf_open(d, "r");
    }
    if(testBCF == NULL){
        fprintf(stderr, "Unable to open %s\n", d);
        fprintf(stderr, "Usage: profiler <indel.bcf> <test.bcf> <qual_threshold> <approx_threshold>\n");
        exit(1);
    }
    
    qual_threshold = (float) atof(argv[3]);
    approx_threshold = atoi(argv[4]);
    approx_threshold = approx_threshold > 0 ? approx_threshold : 0;
    distance = approx_threshold > 0 ? &length_distance : &simple_distance;


    ieof = teof = false;
    totaltruesnps = calledtruesnps = calledfalsesnps = 0;
    for(i = 0; i <= 19; i++){
        totaltrueindel[i] = calledtrueindel[i] = calledfalseindel[i] = 0;
    }
    for(i = 0; i <= 8; i++){
        ltotaltrueindel[i] = lcalledtrueindel[i] = lcalledfalseindel[i] = 0;
    }

    ib = bcf_init();
    tb = bcf_init();
    ih = bcf_hdr_read(indelBCF);
    th = bcf_hdr_read(testBCF);

    { const char **x = bcf_hdr_seqnames(ih, &nseq); free(x); }

    inext();
    tnext();

    while(!(ieof && teof)){
        if(ib->rid < tb->rid){//rid mismatch
            indellen = strlen(ib->d.allele[1]) - strlen(ib->d.allele[0]);
            if(indellen == 0){ //snps
                totaltruesnps++;
            }else if(indellen > 0){ //insertion
                if(indellen <= 10){
                    totaltrueindel[indellen-1]++;
                }else if(indellen <= 30){
                    ltotaltrueindel[1]++;
                }else if(indellen <= 50){
                    ltotaltrueindel[2]++;
                }else if(indellen <= 70){
                    ltotaltrueindel[3]++;
                }else if(indellen <= 90){
                    ltotaltrueindel[4]++;
                }else{
                    ltotaltrueindel[0]++;
                }
            }else{ //deletion
                if(indellen >= -10){
                    totaltrueindel[9-indellen]++;
                }else if(indellen >= -30){
                    ltotaltrueindel[5]++;
                }else if(indellen >= -50){
                    ltotaltrueindel[6]++;
                }else if(indellen >= -70){
                    ltotaltrueindel[7]++;
                }else if(indellen >= -90){
                    ltotaltrueindel[8]++;
                }else{
                    ltotaltrueindel[0]++;
                }
            }
            inext();
            continue;
        }
        if(ib->rid > tb->rid){//rid mismatch
          if(tb->qual >= qual_threshold){
            for(i = 1; i < tb->n_allele; i++){ //multi ALT on same test.bcf entry
                indellen = strlen(tb->d.allele[i]) - strlen(tb->d.allele[0]);
                if(indellen == 0){ //snps
                    calledfalsesnps++;
                }else if(indellen > 0){ //insertion
                    if(indellen <= 10){
                        calledfalseindel[indellen-1]++;
                    }else if(indellen <= 30){
                        lcalledfalseindel[1]++;
                    }else if(indellen <= 50){
                        lcalledfalseindel[2]++;
                    }else if(indellen <= 70){
                        lcalledfalseindel[3]++;
                    }else if(indellen <= 90){
                        lcalledfalseindel[4]++;
                    }else{
                        lcalledfalseindel[0]++;
                    }
                }else{ //deletion
                    if(indellen >= -10){
                        calledfalseindel[9-indellen]++;
                    }else if(indellen >= -30){
                        lcalledfalseindel[5]++;
                    }else if(indellen >= -50){
                        lcalledfalseindel[6]++;
                    }else if(indellen >= -70){
                        lcalledfalseindel[7]++;
                    }else if(indellen >= -90){
                        lcalledfalseindel[8]++;
                    }else{
                        lcalledfalseindel[0]++;
                    }
                }
            }
          }
          tnext();
          continue;
        }
        //rid match
        if(ib->pos + approx_threshold < tb->pos){//pos mismatch
            indellen = strlen(ib->d.allele[1]) - strlen(ib->d.allele[0]);
            if(indellen == 0){ //snps
                totaltruesnps++;
            }else if(indellen > 0){ //insertion
                if(indellen <= 10){
                    totaltrueindel[indellen-1]++;
                }else if(indellen <= 30){
                    ltotaltrueindel[1]++;
                }else if(indellen <= 50){
                    ltotaltrueindel[2]++;
                }else if(indellen <= 70){
                    ltotaltrueindel[3]++;
                }else if(indellen <= 90){
                    ltotaltrueindel[4]++;
                }else{
                    ltotaltrueindel[0]++;
                }
            }else{ //deletion
                if(indellen >= -10){
                    totaltrueindel[9-indellen]++;
                }else if(indellen >= -30){
                    ltotaltrueindel[5]++;
                }else if(indellen >= -50){
                    ltotaltrueindel[6]++;
                }else if(indellen >= -70){
                    ltotaltrueindel[7]++;
                }else if(indellen >= -90){
                    ltotaltrueindel[8]++;
                }else{
                    ltotaltrueindel[0]++;
                }
            }
            inext();
            continue;
        }
        if(ib->pos > tb->pos + approx_threshold){ //pos match up to +- approx_threshold
          if(tb->qual >= qual_threshold){
            for(i = 1; i < tb->n_allele; i++){ //multi ALT on same test.bcf entry
                indellen = strlen(tb->d.allele[i]) - strlen(tb->d.allele[0]);
                if(indellen == 0){ //snps
                    calledfalsesnps++;
                }else if(indellen > 0){ //insertion
                    if(indellen <= 10){
                        calledfalseindel[indellen-1]++;
                        if(indellen == 10){
                            printf("%d\n", tb->pos);
                            exit(0);
                        }
                    }else if(indellen <= 30){
                        lcalledfalseindel[1]++;
                    }else if(indellen <= 50){
                        lcalledfalseindel[2]++;
                    }else if(indellen <= 70){
                        lcalledfalseindel[3]++;
                    }else if(indellen <= 90){
                        lcalledfalseindel[4]++;
                    }else{
                        lcalledfalseindel[0]++;
                    }
                }else{ //deletion
                    if(indellen >= -10){
                        calledfalseindel[9-indellen]++;
                    }else if(indellen >= -30){
                        lcalledfalseindel[5]++;
                    }else if(indellen >= -50){
                        lcalledfalseindel[6]++;
                    }else if(indellen >= -70){
                        lcalledfalseindel[7]++;
                    }else if(indellen >= -90){
                        lcalledfalseindel[8]++;
                    }else{
                        lcalledfalseindel[0]++;
                    }
                }
            }
          }
          tnext();
          continue;
        }
        bool has_match = false;
        if(tb->qual > qual_threshold){
          for(i = 1; i < tb->n_allele; i++){ //multi ALT on same test.bcf entry
            dALT = distance(ib->d.allele[1], tb->d.allele[i]);
            indellen = strlen(tb->d.allele[i]) - strlen(tb->d.allele[0]);
            if(dALT > approx_threshold*2){
                if(indellen == 0){ //snps
                    calledfalsesnps++;
                }else if(indellen > 0){ //insertion
                    if(indellen <= 10){
                        calledfalseindel[indellen-1]++;
                    }else if(indellen <= 30){
                        lcalledfalseindel[1]++;
                    }else if(indellen <= 50){
                        lcalledfalseindel[2]++;
                    }else if(indellen <= 70){
                        lcalledfalseindel[3]++;
                    }else if(indellen <= 90){
                        lcalledfalseindel[4]++;
                    }else{
                        lcalledfalseindel[0]++;
                    }
                }else{ //deletion
                    if(indellen >= -10){
                        calledfalseindel[9-indellen]++;
                    }else if(indellen >= -30){
                        lcalledfalseindel[5]++;
                    }else if(indellen >= -50){
                        lcalledfalseindel[6]++;
                    }else if(indellen >= -70){
                        lcalledfalseindel[7]++;
                    }else if(indellen >= -90){
                        lcalledfalseindel[8]++;
                    }else{
                        lcalledfalseindel[0]++;
                    }
                }
            }else{
                has_match = true;
                indellen = strlen(ib->d.allele[1]) - strlen(ib->d.allele[0]);
                if(indellen == 0){ //snps
                    totaltruesnps++;
                    calledtruesnps++;
                }else if(indellen > 0){ //insertion
                    if(indellen <= 10){
                        totaltrueindel[indellen-1]++;
                        calledtrueindel[indellen-1]++;
                    }else if(indellen <= 30){
                        ltotaltrueindel[1]++;
                        lcalledtrueindel[1]++;
                    }else if(indellen <= 50){
                        ltotaltrueindel[2]++;
                        lcalledtrueindel[2]++;
                    }else if(indellen <= 70){
                        ltotaltrueindel[3]++;
                        lcalledtrueindel[3]++;
                    }else if(indellen <= 90){
                        ltotaltrueindel[4]++;
                        lcalledtrueindel[4]++;
                    }else{
                        ltotaltrueindel[0]++;
                        lcalledtrueindel[0]++;
                    }
                }else{ //deletion
                    if(indellen >= -10){
                        totaltrueindel[9-indellen]++;
                        calledtrueindel[9-indellen]++;
                    }else if(indellen >= -30){
                        ltotaltrueindel[5]++;
                        lcalledtrueindel[5]++;
                    }else if(indellen >= -50){
                        ltotaltrueindel[6]++;
                        lcalledtrueindel[6]++;
                    }else if(indellen >= -70){
                        ltotaltrueindel[7]++;
                        lcalledtrueindel[7]++;
                    }else if(indellen >= -90){
                        ltotaltrueindel[8]++;
                        lcalledtrueindel[8]++;
                    }else{
                        ltotaltrueindel[0]++;
                        lcalledtrueindel[0]++;
                    }
                }
            }
          }
        }
        if(has_match){
            inext();
            tnext();
        }else{
            tnext();
        }
    }




    //report
    int sum[6];
    sum[0] = sum[1] = sum[2] = 0;
    for(i = 0; i <= 19; i++){    
        sum[0] += totaltrueindel[i];
        sum[1] += calledtrueindel[i];
        sum[2] += calledfalseindel[i];
    }
    
    printf("qual threshold %f\n", qual_threshold);
    printf("approx threshold %d\n", approx_threshold);
    printf("total  calledtrue  calledfalse  missed(total-calledtrue)\n");
    printf("short indels\n");
    printf("insertion length 1-10\n");
    for(i = 0; i <= 9; i++)
        printf("%d %d %d %d\n", totaltrueindel[i], calledtrueindel[i], calledfalseindel[i], totaltrueindel[i] - calledtrueindel[i]);
    printf("deletion length 1-10\n");
    for(i = 10; i <= 19; i++)
        printf("%d %d %d %d\n", totaltrueindel[i], calledtrueindel[i], calledfalseindel[i], totaltrueindel[i] - calledtrueindel[i]);
    printf("overall short indels(sum)\n");
        printf("%d %d %d %d\n", sum[0], sum[1], sum[2], sum[0] - sum[1]);

    sum[0] = sum[1] = sum[2] = 0;
    for(i = 0; i <= 8; i++){    
        sum[0] += ltotaltrueindel[i];
        sum[1] += lcalledtrueindel[i];
        sum[2] += lcalledfalseindel[i];
    }

    printf("long indels\n");
    printf("insertion length 20, 40, 60, 80\n");
    for(i = 1; i <= 4; i++)
        printf("%d %d %d %d\n", ltotaltrueindel[i], lcalledtrueindel[i], lcalledfalseindel[i], ltotaltrueindel[i] - lcalledtrueindel[i]);
    printf("deletion length 20, 40, 60, 80\n");
    for(i = 5; i <= 8; i++)
        printf("%d %d %d %d\n", ltotaltrueindel[i], lcalledtrueindel[i], lcalledfalseindel[i], ltotaltrueindel[i] - lcalledtrueindel[i]);
    printf("indel length > 90\n");
        printf("%d %d %d %d\n", ltotaltrueindel[0], lcalledtrueindel[0], lcalledfalseindel[0], ltotaltrueindel[0] - lcalledtrueindel[0]);
    printf("overall long indels(sum)\n");
        printf("%d %d %d %d\n", sum[0], sum[1], sum[2], sum[0] - sum[1]);


    printf("SNPs\n");
        printf("%d %d %d %d\n", totaltruesnps, calledtruesnps, calledfalsesnps, totaltruesnps - calledtruesnps);
    bcf_hdr_destroy(ih);
    bcf_hdr_destroy(th);
    bcf_destroy(ib);
    bcf_destroy(tb);
    bcf_close(indelBCF);
    bcf_close(testBCF);
    return 0;
}


