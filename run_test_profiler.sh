#bin/bash:
# specify the result ID based on the makefile
DDIR=/data/junjie/chr15hg38/ex2       # where all the VCF files are read
DDIRMASK=/data/junjie/chr15hg38/masks # where the masks are read
DDIRRESULTS=$DDIR/profilerResults     # where all the results are saved
filter=repeat50Chr15.txt


##
aligner=novoalign
aID=1
caller=gatk
cID=3
caller=samtools
cID=2
##

make

mkdir catagory

rm profilerOutput.$aligner.$aID.$caller.$cID.$filter

./profiler $DDIR/indel.vcf $DDIR/$aligner.$aID.$caller.$cID 0 5 $DDIRMASK/$filter >> profilerOutput.$aligner.$aID.$caller.$cID.$filter



# aligner=bowtie2
# aID=2
# caller=gatk
# cID=3
# echo align command ================= > $DDIR/result.$aligner.$aID.$caller.$cID
# cat $aligner.sh >> $DDIR/result.$aligner.$aID.$caller.$cID
# 
# echo align start =================== >> $DDIR/result.$aligner.$aID.$caller.$cID
# (time . $aligner.sh) &>> $DDIR/result.$aligner.$aID.$caller.$cID
# 
# echo call command ================== >> $DDIR/result.$aligner.$aID.$caller.$cID
# cat $caller.sh >> $DDIR/result.$aligner.$aID.$caller.$cID
# 
# echo call start ==================== >> $DDIR/result.$aligner.$aID.$caller.$cID
# (time . $caller.sh) &>> $DDIR/result.$aligner.$aID.$caller.$cID

#echo >> result.$rID

#echo =============================== >> result.$rID

#./profiler $DDIR/indel.bcf $DDIR/$aligner.$aID.$caller.$cID 0 5 >> result.$rID

#echo =============================== >> result.$rID

#./profiler $DDIR/indel.bcf $DDIR/$aligner.$aID.$caller.$cID 0 5 $filter1 >> result.$rID

#echo =============================== >> result.$rID

#./profiler $DDIR/indel.bcf $DDIR/$aligner.$aID.$caller.$cID 0 5 $filter2 >> result.$rID

#echo =============================== >> result.$rID

#./profiler $DDIR/indel.bcf $DDIR/$aligner.$aID.$caller.$cID 0 5 $filter3 >> result.$rID

#echo =============================== >> result.$rID

#./profiler $DDIR/indel.bcf $DDIR/$aligner.$aID.$caller.$cID 0 5 $filter4 >> result.$rID

#echo =============================== >> result.$rID

#./profiler $DDIR/indel.bcf $DDIR/$aligner.$aID.$caller.$cID 0 5 $filter5 >> result.$rID


