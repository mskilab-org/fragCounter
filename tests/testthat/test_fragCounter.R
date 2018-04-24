context("unit testing fragCounter operations")

library(fragCounter)

### Making small BAM
# cd /gpfs/commons/groups/imielinski_lab/projects/PurityPloidy/Flow/Remix_Bams/HCC1143_1_5_4 1020
# samtools view -H downsample.bam > header_downsample.sam
# samtools view downsample.bam  | grep "chr21" | cat  header_downsample.sam  - | samtools view -Sb - > chr21.unique.bam
## $ samtools view -S -b chr21.unique.sam > chr21.bam
## $ samtools index chr21.bam

example_bam = 'chr21.bam'
example_bai = 'chr21.bai'

test_that("PrepareCov", {

  expect_error(PrepareCov(bam = NULL, cov = NULL))
  expect_equal(length(PrepareCov(bam = example_bam)),15509063)


})
