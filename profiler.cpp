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
    bool ieof, teof;
    int indellen;
    int dALT;
    bool has_match; //duplicated match count as one
    bool match_pos; 


    int i;

    VarCount totaltrue;                              // false negatives + true positives
    VarCount concord;                                // true positive, concordance (true calls)
    VarCount discord("catagory/discord.txt");        // discordance (false calls in new positions)
    VarCount pconcord("catagory/pconcord.txt");      // partial concordance (false calls in correct positions)
    VarCount nonconcord("catagory/nonconcord.txt");  // false negatives (discordance + partial concordance)
    VarCount missing("catagory/missing.txt");        // variants in the reference missed by callers
    
    
    int duplicate = 0;
    int iskip = 0;
    int tskip = 0;
    int qskip = 0;

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
        // ---- END of simple repeat filter stuff


        // if all of the called variant entries have been looked up (EOF),
        // then the true variant entries left are true (and missed) calls
        if(ib->rid < tb->rid){//rid mismatch
            totaltrue.update(ib->pos, ib->d.allele[0],ib->d.allele[1]);
            if (!has_match){ // if previously we had a match, then DO NOT count is as missed
                missing.update(ib->pos, ib->d.allele[0],ib->d.allele[1]);
            }
            inext(); has_match = false; match_pos = false;
            continue;
        }

        // if all of the true variant entries have been looked up (EOF),
        // then the called variant entries left are false calls
        if(ib->rid > tb->rid){//rid mismatch
            for(i = 1; i < tb->n_allele; i++){ //multi ALT on same test.bcf entry
                discord.update(tb->pos, tb->d.allele[0], tb->d.allele[i]);
            }
            tnext(); match_pos = false;
            continue;
        }

        // if the pos of true variant location has not reached the pos of 
        // the called variant, then move to the pos. of next true variant
        // this counts as true (and missed) calls
        if(ib->pos < tb->pos - approx_threshold){//pos mismatch
            totaltrue.update(ib->pos, ib->d.allele[0],ib->d.allele[1]);
            if (!match_pos){ // if previously we had a match, then DO NOT count it as missed
                missing.update(ib->pos, ib->d.allele[0],ib->d.allele[1]);
            }
            inext(); has_match = false; match_pos = false;
            continue;
        }

        // if the pos of called variant location has not reached the pos of 
        // the true variant, then move to the pos. of next called variant
        // this counts as false calls
        if(ib->pos > tb->pos + approx_threshold){ //pos mismatch
            for(i = 1; i < tb->n_allele; i++){ //multi ALT on same test.bcf entry
                discord.update(tb->pos, tb->d.allele[0], tb->d.allele[i]);
            }
            tnext(); match_pos = true;
            continue;
        }

        // pos of true and called match within +- approx_threshold
        // there might be multiple entries in <test.bcf> that match an entry in <indel.bcf>
        for(i = 1; i < tb->n_allele; i++){ //multi ALT on same test.bcf entry
            dALT = distance(ib->d.allele[1], tb->d.allele[i]);
            if(dALT > approx_threshold*2){
                pconcord.update(tb->pos, tb->d.allele[0], tb->d.allele[i]); 
            }else{
              if(!has_match){      
                has_match = true;
                concord.update(ib->pos, ib->d.allele[0], ib->d.allele[1]);
              }else{  // if we have already recorded the matched indel, we should not record again
                duplicate++;
              }
            }
        }
        tnext(); match_pos = true;
    }




    // ============
    // BEGIN REPORT

    int sum[6];
    sum[0] = sum[1] = sum[2] = 0;
    for(i = 0; i < 20; i++){    
        sum[0] += totaltrue.indel[i];
        sum[1] += concord.indel[i];
        sum[2] += pconcord.indel[i];
        sum[3] += discord.indel[i];
        sum[4] += missing.indel[i];
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
    // totaltrue   = FN+TP
    // calledtrue  = TP
    // calledfalse = FP
    // missed      = FN

    // FP/(TN+FP) = Type I Error, False Pos. Rate, 1-Specitivity.
    // TN/(TN+FP) = Specitivity 
    // FN/(FN+TP) = Type II Error  **** = miss detection rate ***
    // TP/(FN+TP) = Sensitiviy, Recall, True Pos. Rate 
    // TP/(TP+FP) = Precision, Pos. Pred. Value 
    // FP/(TP+FP) = False Discovery Proportion
    

    //printf("SNP_FA   : %lf\n", 1.0*calledfalsesnps/calledtruesnps);
    
    printf("SNP_Sensitivity  : %lf\n", 1.0*concord.snp/totaltrue.snp); 
    printf("INDEL_Sensitivity: %lf\n", 1.0*sum[1]/sum[0]);
    printf("SNP_FalDisRate   : %lf\n", 1.0*(discord.snp+pconcord.snp)/(discord.snp+pconcord.snp+concord.snp));
    printf("INDEL_FalDisRate : %lf\n", 1.0*(sum[3]+sum[2])/(sum[3]+sum[2]+sum[1]));


    printf("=============== Counts ==============\n"); 
    printf("totaltrue/concordance/p-concordance/discordance/missing\n");
    printf("-------------------------------------\n");
    printf("                 SNP                 \n");
    printf("-------------------------------------\n");
    printf("%d\t%d\t%d\t%d\t%d\n", totaltrue.snp, concord.snp, pconcord.snp, discord.snp, missing.snp);
    printf("-------------------------------------\n");
    printf("                INDEL                \n");
    printf("-------------------------------------\n");
    printf("%d\t%d\t%d\t%d\t%d\n", sum[0], sum[1], sum[2], sum[3], sum[4]);
    printf("-------------------------------------\n");
    printf("  Breakdown:  insertion length 1-10  \n");
    printf("-------------------------------------\n");
    for(i = 0; i <= 9; i++)
        printf("%d\t%d\t%d\t%d\t%d\n", totaltrue.indel[i], concord.indel[i], pconcord.indel[i], discord.indel[i], missing.indel[i]);
    printf("-------------------------------------\n");
    printf("  Breakdown:  deletion length 1-10   \n");
    printf("-------------------------------------\n");
    for(i = 10; i <= 19; i++)
        printf("%d\t%d\t%d\t%d\t%d\n", totaltrue.indel[i], concord.indel[i], pconcord.indel[i], discord.indel[i], missing.indel[i]);
    printf("-------------------------------------\n");
    printf("* Duplicates : %d\n", duplicate);

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





