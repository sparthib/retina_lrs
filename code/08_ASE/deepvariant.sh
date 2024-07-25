
# Pull the image.
BIN_VERSION=1.6.1
ml load singularity
singularity pull docker://google/deepvariant:"${BIN_VERSION}"

# singularity pull docker://google/deepvariant:"${BIN_VERSION}"
# FATAL:   Image file already exists: "deepvariant_1.6.1.sif" - will not overwrite

PWD=/dcs04/hicks/data/sparthib/retina_lrs/09_deepVariant
INPUT_DIR="${PWD}/quickstart-testdata"
DATA_HTTP_DIR="https://storage.googleapis.com/deepvariant/quickstart-testdata"
mkdir -p ${INPUT_DIR}
wget -P ${INPUT_DIR} "${DATA_HTTP_DIR}"/NA12878_S1.chr20.10_10p1mb.bam
wget -P ${INPUT_DIR} "${DATA_HTTP_DIR}"/NA12878_S1.chr20.10_10p1mb.bam.bai
wget -P ${INPUT_DIR} "${DATA_HTTP_DIR}"/test_nist.b37_chr20_100kbp_at_10mb.bed
wget -P ${INPUT_DIR} "${DATA_HTTP_DIR}"/test_nist.b37_chr20_100kbp_at_10mb.vcf.gz
wget -P ${INPUT_DIR} "${DATA_HTTP_DIR}"/test_nist.b37_chr20_100kbp_at_10mb.vcf.gz.tbi
wget -P ${INPUT_DIR} "${DATA_HTTP_DIR}"/ucsc.hg19.chr20.unittest.fasta
wget -P ${INPUT_DIR} "${DATA_HTTP_DIR}"/ucsc.hg19.chr20.unittest.fasta.fai
wget -P ${INPUT_DIR} "${DATA_HTTP_DIR}"/ucsc.hg19.chr20.unittest.fasta.gz
wget -P ${INPUT_DIR} "${DATA_HTTP_DIR}"/ucsc.hg19.chr20.unittest.fasta.gz.fai
wget -P ${INPUT_DIR} "${DATA_HTTP_DIR}"/ucsc.hg19.chr20.unittest.fasta.gz.gzi

OUTPUT_DIR="${PWD}/test_output"


# Run DeepVariant.
singularity run -B /usr/lib/locale/:/usr/lib/locale/ \
  --bind ${INPUT_DIR} --bind ${OUTPUT_DIR} \
  docker://google/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=WGS \
  --ref="${INPUT_DIR}"/ucsc.hg19.chr20.unittest.fasta \
  --reads="${INPUT_DIR}"/NA12878_S1.chr20.10_10p1mb.bam \
  --output_vcf="${OUTPUT_DIR}"/output.vcf.gz \
  --output_gvcf="${OUTPUT_DIR}"/output.g.vcf.gz \
  --
  --num_shards=1

