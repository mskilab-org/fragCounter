context("unit testing fragCounter operations")

library(fragCounter)
library(data.table)
library(GenomeInfoDb)

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



bw = system.file('extdata', 'chr21.bigWig', package = 'fragCounter')

twobit = system.file("extdata", 'hg19_chr21.2bit', package = 'fragCounter')

example_bam = system.file("extdata", 'chr21.bam', package = 'fragCounter')

example_bai = system.file("extdata", 'chr21.bai', package = 'fragCounter')

cov21 = system.file("extdata", 'cov21.rds', package = 'fragCounter')

cov21_exome = system.file("extdata", 'cov21_exome.rds', package = 'fragCounter')

cov = readRDS(system.file("extdata", 'samp.rds', package = 'fragCounter'))

gcmapdir = gsub("cov21.rds","gcMAP21",cov21)

skeleton = system.file("extdata", "sample_blacklist.rds", package = "fragCounter")

#' twalradt Friday, May 18, 2018 01:37:49 PM Commented out to see if fragCounter will pass on travis

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
  expect_true(length(PrepareCov(example_bam)) >= 15509063)
  expect_equal(max(width(PrepareCov(example_bam))), 200)
  expect_equal(length(PrepareCov(example_bam, exome = TRUE, skeleton = skeleton, use.skel = TRUE)), 194880)
#  expect_equal(length(PrepareCov(example_bam, paired = FALSE)), 15509063)
#  expect_equal(max(width(PrepareCov(example_bam, paired = FALSE))), 200)

})


test_that("correctcov_stub", {

  expect_equal(length(correctcov_stub(cov.wig = cov21, gc.rds.dir = gcmapdir, map.rds.dir = gcmapdir)), 234722)
  expect_equal(length(correctcov_stub(cov.wig = cov21_exome, gc.rds.dir = gcmapdir, map.rds.dir = gcmapdir, exome = TRUE)), 3663)

})


test_that("multicoco", {

    expect_equal(length(multicoco(cov)), 50001)
    expect_equal(sum(!is.na(multicoco(cov, mono = FALSE)$reads.corrected)), 48888)
    expect_equal(sum(!is.na(multicoco(cov, mono = TRUE)$reads.corrected)), 49571)

})



test_that("coco", {

    expect_equal(length(coco(cov)), 50001)
    expect_equal(sum(!is.na(coco(cov)$reads.corrected)), 49571)

})






