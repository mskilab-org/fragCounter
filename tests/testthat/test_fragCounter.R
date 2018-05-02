context("unit testing fragCounter operations")

library(fragCounter)

### Making small BAM
# cd /gpfs/commons/groups/imielinski_lab/projects/PurityPloidy/Flow/Remix_Bams/HCC1143_1_5_4 1020
# samtools view -H downsample.bam > header_downsample.sam
# samtools view downsample.bam  | grep "chr21" | cat  header_downsample.sam  - | samtools view -Sb - > chr21.unique.bam
## $ samtools view -S -b chr21.unique.sam > chr21.bam
## $ samtools index chr21.bam

### Sample 2bit file w/ chr21
# samtools faidx hg19.fa chr21 > hg19_chr21.fa
# /gpfs/commons/home/twalradt/software/Blat/faToTwoBit hg19_chr21.fa hg19_chr21.2bit

### Create sample bigWig file
# bw.path = '~/DB/UCSC/wgEncodeCrgMapabilityAlign100mer.bigWig'
# bw = rtracklayer::import(bw.path, selection = tiles)
# export(bw, "~/git/fragCounter/tests/testthat/chr21.bigWig", 'bigWig')

bw = 'chr21.bigWig'
twobit = 'hg19_chr21.2bit'

example_bam = 'chr21.bam'
example_bai = 'chr21.bai'

cov21 = 'cov21.rds'

cov = readRDS('samp.rds')

test_that("MAP.fun", {

  expect_equal(length(MAP.fun(twobitURL = twobit, bw.path = bw)), 240650)
  expect_equal(length(MAP.fun(twobitURL = twobit, bw.path = bw) %Q% (score > 0)), 171124)

})



test_that("GC.fun", {

  expect_equal(length(GC.fun(twobitURL = twobit)), 240650)
  expect_equal(length(GC.fun(twobitURL =twobit) %Q% (score > 0)), 175542)

})


test_that("PrepareCov", {

  expect_error(PrepareCov(bam = NULL, cov = NULL))
  expect_equal(length(PrepareCov(example_bam)), 15509063)
  expect_equal(max(width(PrepareCov(example_bam))), 200)
#  expect_equal(length(PrepareCov(example_bam, paired = FALSE)), 15509063)
#  expect_equal(max(width(PrepareCov(example_bam, paired = FALSE))), 200)

})


test_that("correctcov_stub", {

    expect_equal(length(correctcov_stub(cov.wig = cov21, gc.rds.dir = "gcMAP21", map.rds.dir = "gcMAP21")), 234722)

})


test_that("multicoco", {

    expect_equal(length(multicoco(cov)), 50001)
    expect_equal(sum(!is.na(multicoco(cov)$reads.corrected)), 48888)
    expect_equal(sum(!is.na(multicoco(cov, mono = TRUE)$reads.corrected)), 49571)

})









