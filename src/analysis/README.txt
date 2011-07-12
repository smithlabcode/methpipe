---------------MAPPING---------------------

1. map paired-end reads of ESC

   ./rmapbs -o test_ESC_1_mapped.mr -c test_ref_dir -m 3 -y test_ESC_1.fastq
   ./rmapbs -o test_ESC_2_mapped.mr -c test_ref_dir -m 3 -y test_ESC_2.fastq

2. map single-end reads of NHFF

   ./rmapbs -o test_NHFF_mapped.mr -c test_ref_dir -m 3 -y test_NHFF.fastq


--------------POSTMAPPING PAIRED-END MAPPED READS-----------------

 
1. sort mapped reads by read ID (name)

   sort -k 4,4 test_ESC_1_mapped.mr > test_ESC_1_sortname.mr
   sort -k 4,4 test_ESC_2_mapped.mr > test_ESC_2_sortname.mr

2. mask overlapping region of mates; for mates2: take reverse-complement and switch strands

   ./clipmates -T test_ESC_1_sortname.mr -A test_ESC_2_sortname.mr 
               -S test_ESC_clipmates_stat.txt -o test_ESC_clipmates.mr -L 500


---------------POSTMAPPING MAPPED READS---------------------------

1. sort by genome location: chrom, start, end, strand

   sort -k 1,1 -k 2,2g -k 3,3g -k 6,6 test_ESC_clipmates.mr > test_ESC_sorted_start.mr
   sort -k 1,1 -k 2,2g -k 3,3g -k 6,6 test_NHFF_mapped.mr > test_NHFF_sorted_start.mr

2. remove copy-duplicates, the reads mapped to the same start position 

   ./duplicate-remover -S test_ESC_duplremove_start_stat.txt 
                       -o test_ESC_duplremove_start.mr  test_ESC_sorted_start.mr

3. sort by genome location: chrom, end, start, strand

   sort -k 1,1 -k 3,3g -k 2,2g -k 6,6 test_ESC_duplremove_start.mr > test_ESC_sorted_end.mr
   sort -k 1,1 -k 3,3g -k 2,2g -k 6,6 test_NHFF_duplremove_start.mr > test_NHFF_sorted_end.mr
   
4. remove copy-duplicates, the reads mapped to the same end position

   ./duplicate-remover -S test_ESC_duplremove_end_stat.txt 
                       -o test_ESC_distinct.mr test_ESC_sorted_end.mr -B

   ./duplicate-remover -S test_NHFF_duplremove_end_stat.txt
                       -o test_NHFF_distinct.mr test_NHFF_sorted_end.mr -B

5. calculate methylation level

   ./methcounts -M 40 -c test_ref_dir -o test_ESC_methcounts.bed test_ESC_distinct.mr
   ./methcounts -N -M 40 -c test_ref_dir -o test_NHFF_methcounts.bed test_NHFF_distinct.mr

6. calculate allelic specific methylation scores

   ./allelicmeth -M 40 -c test_ref_dir -o test_ESC_allelicmeth.bed test_ESC_distinct.mr
   ./allelicmeth -N -M 40 -c test_ref_dir -o test_NHFF_allelicmeth.bed test_NHFF_distinct.mr

7. calculate bisulfite conversion rate

   ./bsrate -c test_ref_dir -M 40 -o test_ESC_bsrate.txt test_ESC_distinct.mr
   ./bsrate -c test_ref_dir -M 40 -o test_NHFF_bsrate.txt test_NHFF_distinct.mr

8. identify HMRs

   ./hmr_finder -i 20 -o test_ESC_HMR.bed test_ESC_methcounts.bed
   ./hmr_finder -i 20 -o test_NHFF_HMR.bed test_NHFF_methcounts.bed

9. calculate differential methylation scores

   ./methdiff -o test_NHFF_ESC_methdiff.bed test_NHFF_methcounts.bed test_ESC_methcounts.bed

10. identify DMRs

   ./dmr_finder -o test_NHFF_ESC_DMR.bed -d 500 -b 10 -C 10 -m 200 -c 0.7 test_NHFF_ESC_methdif.bed

