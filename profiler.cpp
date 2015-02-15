#include"profiler.h"

//assume chromosome order is the same in each files!!


int main(int argc, char **argv){

    htsFile *indelBCF, *testBCF;
    srFILE_t *isr = NULL, *tsr = NULL;

    bcf1_t *ib, *tb;
    bcf_hdr_t *ih, *th;
    sr1_t *is, *ts;

    int nseq;
    int ir, tr;
    int32_t isr_rgnStart, tsr_rgnStart;

    float qual_threshold;
    int approx_threshold;
    int totaltruesnps, calledtruesnps, calledfalsesnps;
    bool ieof, teof;
    int indellen;
    int dALT;
    bool has_match; //duplicated match count as one

    int i;
    int totaltrueindel[20]; //0-9 ins //10-19 del
    int calledtrueindel[20]; //0-9 ins //10-19 del
    int calledfalseindel[20]; //0-9 ins //10-19 del
    int ltotaltrueindel[9]; //1 11-30 ins //2 31-50 ins //3 51-70 ins //4 71-90 ins
    int lcalledtrueindel[9]; //5 11-30 del //6 31-50 del //7 51-70 del //8 71-90 del
    int lcalledfalseindel[9]; //0 91+ indel
    int duplicate = 0;
    int iskip = 0;
    int tskip = 0;
    int qskip = 0;

    // ******  test info:
    int misssnp = 0;

    char *d;

    int (*distance)(const char *, const char *);

    if((argc != 5) && (argc != 6)){
        fprintf(stderr, "Usage: profiler <indel.bcf> <test.bcf> <qual_threshold> <approx_threshold> [simpleRepeat.txt] \n");
        exit(1); 
    }

    indelBCF = bcf_open(argv[1], "r");
    if(indelBCF == NULL){
        fprintf(stderr, "Unable to open %s\n", argv[1]);
        fprintf(stderr, "Usage: profiler <indel.bcf> <test.bcf> <qual_threshold> <approx_threshold> [simpleRepeat.txt] \n");
        exit(1);
    }

    // 
    d = (char*) malloc((strlen(argv[2])+5)*sizeof(char));
        sprintf(d, "%s.bcf", argv[2]);
        testBCF = bcf_open(d, "r");
    if(testBCF == NULL){
        d[strlen(d)-3] = 'v';// sprintf(d, "%s.vcf", argv[2];
        testBCF = vcf_open(d, "r");
    }
    if(testBCF == NULL){
        fprintf(stderr, "Unable to open %s\n", d);
        fprintf(stderr, "Usage: profiler <indel.bcf> <test.bcf> <qual_threshold> <approx_threshold> [simpleRepeat.txt] \n");
        exit(1);
    }
    
    qual_threshold = (float) atof(argv[3]);
    approx_threshold = atoi(argv[4]);
    approx_threshold = approx_threshold > 0 ? approx_threshold : 0;
    distance = approx_threshold > 0 ? &length_distance : &simple_distance;

    if(argc == 6){
      isr = sr_open(argv[5], "r");
      tsr = sr_open(argv[5], "r");
      if((isr == NULL) || (tsr == NULL)){
        fprintf(stderr, "Unable to open %s\n", argv[5]);
        fprintf(stderr, "Usage: profiler <indel.bcf> <test.bcf> <qual_threshold> <approx_threshold> [simpleRepeat.txt] \n");
        exit(1);
      }
    }

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
    if(isr){ 
        is = sr_init1(); ts = sr_init1();
        sr_read1(isr, is); sr_read1(tsr, ts);
    }
    { const char **x = bcf_hdr_seqnames(ih, &nseq); free(x); }

    inext(); has_match = false;
    tnext();

    while(!(ieof && teof)){
        //filter repeat region: skip ib within any [is->chromStart, is->chromEnd]
        //filter repeat region: skip tb within any [is->chromStart+approx, is->chromEnd-approx]
        if(isr){
            if(!ieof){
                if(ib->rid < (int32_t) is->chromID){
                    isr_rgnStart = INT32_MAX;
                }else if(ib->rid > (int32_t) is->chromID){
                    if(sr_read1(isr, is) == EOF){ //critical
                        fprintf(stderr, "error in simpleRepeat.txt\n");
                        exit(1);
                    }
                    continue;
                }else{
                    if(ib->pos > (int32_t) is->chromEnd){
                        if(sr_read1(isr, is) == EOF){ //non critical
                            is->chromID = isr->nchrom;    
                        }
                        continue;
                    }
                    isr_rgnStart = is->chromStart;
                }
                if(ib->pos >= (int32_t) isr_rgnStart){
                    inext(); has_match = false; iskip++;
                    continue;
                }
            }
            if(!teof){
                if(tb->rid < (int32_t) ts->chromID){
                    tsr_rgnStart = INT32_MAX;
                }else if(tb->rid > (int32_t) ts->chromID){
                    if(sr_read1(tsr, ts) == EOF){ //critical
                        fprintf(stderr, "error in simpleRepeat.txt\n");
                        exit(1);
                    }
                    continue;
                }else{
                    if(tb->pos > (int32_t) ts->chromEnd){
                        if(sr_read1(tsr, ts) == EOF){ //non critical
                            ts->chromID = tsr->nchrom;
                        }
                        continue;
                    }
                    tsr_rgnStart = ts->chromStart;
                }
                if(tb->pos >= (int32_t) tsr_rgnStart){
                    tnext(); has_match = false; tskip++;
                    continue;
                }
            }
        }
        if(tb->qual < qual_threshold){
            qskip++;
            tnext();
            continue;
        }


        // if all of the called variant entries have been looked up (EOF),
        // then the true variant entries left are true (and missed) calls
        if(ib->rid < tb->rid){//rid mismatch
            indellen = strlen(ib->d.allele[1]) - strlen(ib->d.allele[0]);
            if(indellen == 0){ //snps
                totaltruesnps++;


                misssnp++;


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
            inext(); has_match = false;
            continue;
        }

        // if all of the true variant entries have been looked up (EOF),
        // then the called variant entries left are false calls
        if(ib->rid > tb->rid){//rid mismatch
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
            tnext();
            continue;
        }


        // if the pos of true variant location has not reached the pos of 
        // the called variant, then move to the pos. of next true variant
        // this counts as true (and missed) calls

        if(ib->pos < tb->pos - approx_threshold){//pos mismatch
            indellen = strlen(ib->d.allele[1]) - strlen(ib->d.allele[0]);
            if(indellen == 0){ //snps
                totaltruesnps++;


                    
                    misssnp++;


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
            // **** CHECK IF THIS IS TRUE
//          if(indellen && !has_match) printf("%d\n", ib->pos); //missed indel //FIXME: total of these does not agree with totaltrue - calledtrue
            inext(); has_match = false;
            continue;
        }


        // if the pos of called variant location has not reached the pos of 
        // the true variant, then move to the pos. of next called variant
        // this counts as false calls

        if(ib->pos > tb->pos + approx_threshold){ //pos mismatch
//          printf("pos mismatch %d\n", tb->pos);
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
            tnext();
            continue;
        }


        // pos of true and called match within +- approx_threshold
        //

        for(i = 1; i < tb->n_allele; i++){ //multi ALT on same test.bcf entry
            dALT = distance(ib->d.allele[1], tb->d.allele[i]);
            indellen = strlen(tb->d.allele[i]) - strlen(tb->d.allele[0]);
            if(dALT > approx_threshold*2){
//              printf("pos mismatch %d\n", tb->pos);
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
              if(!has_match){
                has_match = true;
                indellen = strlen(ib->d.allele[1]) - strlen(ib->d.allele[0]);
                if(indellen == 0){ //snps
                    calledtruesnps++;
                }else if(indellen > 0){ //insertion
                    if(indellen <= 10){
                        calledtrueindel[indellen-1]++;
                    }else if(indellen <= 30){
                        lcalledtrueindel[1]++;
                    }else if(indellen <= 50){
                        lcalledtrueindel[2]++;
                    }else if(indellen <= 70){
                        lcalledtrueindel[3]++;
                    }else if(indellen <= 90){
                        lcalledtrueindel[4]++;
                    }else{
                        lcalledtrueindel[0]++;
                    }
                }else{ //deletion
                    if(indellen >= -10){
                        calledtrueindel[9-indellen]++;
                    }else if(indellen >= -30){
                        lcalledtrueindel[5]++;
                    }else if(indellen >= -50){
                        lcalledtrueindel[6]++;
                    }else if(indellen >= -70){
                        lcalledtrueindel[7]++;
                    }else if(indellen >= -90){
                        lcalledtrueindel[8]++;
                    }else{
                        lcalledtrueindel[0]++;
                    }
                }
              }else{
                duplicate++;
              }
            }
        }
        tnext();
    }


    // ============
    // BEGIN REPORT


    int sum[6];
    sum[0] = sum[1] = sum[2] = 0;
    for(i = 0; i <= 19; i++){    
        sum[0] += totaltrueindel[i];
        sum[1] += calledtrueindel[i];
        sum[2] += calledfalseindel[i];
    }


    printf("============ Filter Info ============\n");
    if(isr){
        printf("File: %s\n", isr->fn);
        printf("<indel.bcf>skip : %d\n", iskip);
        printf("<test.bcf>skip  : %d\n", tskip);
    }

    printf("============ Thresholds ============\n"); 
    printf("qual_thres : %f\n", qual_threshold);
    printf("qual_skip  : %d\n", qskip);
    printf("appr_thres : %d\n", approx_threshold);

    printf("============ Error Rates ============\n"); 
    printf("SNP_FA   : %lf\n", 1.0*calledfalsesnps/calledtruesnps);
    printf("SNP_MD   : %lf\n", 1.0*(totaltruesnps-calledtruesnps)/totaltruesnps);
    printf("INDEL_FA : %lf\n", 1.0*sum[2]/sum[1]);
    printf("INDEL_MD : %lf\n", 1.0*(sum[0]-sum[1])/sum[0]);

    printf("=============== Counts ==============\n"); 
    printf("total  calledtrue  calledfalse  missed(total-calledtrue)\n");
    printf("SNP: \n");
    printf("%d\t%d\t%d\t%d\n", totaltruesnps, calledtruesnps, calledfalsesnps, totaltruesnps - calledtruesnps);
    printf("INDEL: \n");
    printf("%d\t%d\t%d\t%d\n", sum[0], sum[1], sum[2], sum[0] - sum[1]);
    printf("INDEL breakdown: \n");
    printf("insertion length 1-10\n");
    for(i = 0; i <= 9; i++)
        printf("%d\t%d\t%d\t%d\n", totaltrueindel[i], calledtrueindel[i], calledfalseindel[i], totaltrueindel[i] - calledtrueindel[i]);
    printf("deletion length 1-10\n");
    for(i = 10; i <= 19; i++)
        printf("%d\t%d\t%d\t%d\n", totaltrueindel[i], calledtrueindel[i], calledfalseindel[i], totaltrueindel[i] - calledtrueindel[i]);
    
    printf("* Duplicates : %d\n", duplicate);

    // *** test
    printf("* missed SNPs : %d\n", misssnp);

    // END REPORT
    // ============


    bcf_hdr_destroy(ih);
    bcf_hdr_destroy(th);
    bcf_destroy(ib);
    bcf_destroy(tb);
    bcf_close(indelBCF);
    bcf_close(testBCF);
    if(isr){
        sr_close(isr);
        sr_close(tsr);
    }
    return 0;
}


