#!/bin/bash

#===============================================================================
# SCRIPT: 2_align.sh
# PURPOSE: STAR Alignment for RiboTag RNA-seq Data
#
# DESCRIPTION:
# Performs splice-aware alignment of paired-end FASTQ files to the mouse
# reference genome (GRCm39/mm39) using STAR aligner.
#
# NOTE: This pipeline skips trimming since the data from Azenta is high quality
# (Q38+, >92% bases ≥Q30). If adapter contamination is detected in FastQC,
# add a trimming step.
#
# USAGE:
# sbatch 2_align.sh
#===============================================================================

#SBATCH --job-name=2_align_ribotag
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --array=0-3
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/2_align_%a.err"
#SBATCH --output="./logs/2_align_%a.out"

# Set up conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/rnaseq-quant

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/90-1265107649"
RAW_DIR="${BASE_DIR}/00_fastq"
OUTPUT_DIR="${BASE_DIR}/results/02_aligned"

# Reference genome and annotation (Mouse GRCm39)
GENOME_INDEX="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/genome/STAR_GRCm39"
GTF_FILE="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/refdata-gex-GRCm39-2024-A/genes/genes.gtf"

SAMPLES=($(cat ${BASE_DIR}/config/samples.txt))
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

echo "=== Starting alignment for sample: ${SAMPLE} ==="
echo "Timestamp: $(date)"
echo "Input files:"
echo "  R1: ${RAW_DIR}/${SAMPLE}_R1_001.fastq.gz"
echo "  R2: ${RAW_DIR}/${SAMPLE}_R2_001.fastq.gz"
echo "Output directory: ${OUTPUT_DIR}"
echo "Genome index: ${GENOME_INDEX}"
echo "GTF file: ${GTF_FILE}"
echo ""

# Check if input files exist
R1_FILE="${RAW_DIR}/${SAMPLE}_R1_001.fastq.gz"
R2_FILE="${RAW_DIR}/${SAMPLE}_R2_001.fastq.gz"

if [[ ! -f "${R1_FILE}" ]]; then
    echo "ERROR: R1 file not found: ${R1_FILE}"
    exit 1
fi
if [[ ! -f "${R2_FILE}" ]]; then
    echo "ERROR: R2 file not found: ${R2_FILE}"
    exit 1
fi

# Check if genome index exists
if [[ ! -d "${GENOME_INDEX}" ]]; then
    echo "ERROR: STAR genome index not found: ${GENOME_INDEX}"
    exit 1
fi

# Create sample-specific output directory
SAMPLE_OUTPUT="${OUTPUT_DIR}/${SAMPLE}"
mkdir -p ${SAMPLE_OUTPUT}

echo "=== Step 1: Running STAR alignment ==="
echo "Timestamp: $(date)"

# STAR alignment with RNA-seq optimized parameters
STAR \
    --runThreadN 32 \
    --genomeDir ${GENOME_INDEX} \
    --sjdbGTFfile ${GTF_FILE} \
    --readFilesIn ${R1_FILE} ${R2_FILE} \
    --readFilesCommand zcat \
    --outFileNamePrefix ${SAMPLE_OUTPUT}/${SAMPLE}_ \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMattributes Standard \
    --outFilterType BySJout \
    --outFilterMultimapNmax 20 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --sjdbScore 1 \
    --quantMode GeneCounts \
    --twopassMode Basic

echo "=== Step 2: Indexing BAM file ==="
echo "Timestamp: $(date)"

# Index the BAM file
samtools index -@ 32 ${SAMPLE_OUTPUT}/${SAMPLE}_Aligned.sortedByCoord.out.bam

echo "=== Step 3: Generating alignment statistics ==="
echo "Timestamp: $(date)"

# Generate alignment statistics
samtools flagstat ${SAMPLE_OUTPUT}/${SAMPLE}_Aligned.sortedByCoord.out.bam > ${SAMPLE_OUTPUT}/${SAMPLE}_flagstat.txt
samtools stats ${SAMPLE_OUTPUT}/${SAMPLE}_Aligned.sortedByCoord.out.bam > ${SAMPLE_OUTPUT}/${SAMPLE}_stats.txt

echo "=== Step 4: Checking output files ==="
BAM_FILE="${SAMPLE_OUTPUT}/${SAMPLE}_Aligned.sortedByCoord.out.bam"
COUNTS_FILE="${SAMPLE_OUTPUT}/${SAMPLE}_ReadsPerGene.out.tab"
LOG_FILE="${SAMPLE_OUTPUT}/${SAMPLE}_Log.final.out"

if [[ -f "${BAM_FILE}" ]]; then
    BAM_SIZE=$(ls -lh "${BAM_FILE}" | awk '{print $5}')
    echo "SUCCESS: BAM file created - ${BAM_FILE} (${BAM_SIZE})"

    # Quick read count check
    READ_COUNT=$(samtools view -c "${BAM_FILE}")
    echo "Total reads in BAM: ${READ_COUNT}"
else
    echo "ERROR: BAM file not created!"
    exit 1
fi

if [[ -f "${COUNTS_FILE}" ]]; then
    GENE_COUNT=$(tail -n +5 "${COUNTS_FILE}" | wc -l)
    echo "SUCCESS: Gene counts file created - ${COUNTS_FILE}"
    echo "Genes with counts: ${GENE_COUNT}"
else
    echo "WARNING: Gene counts file not created!"
fi

if [[ -f "${LOG_FILE}" ]]; then
    echo "SUCCESS: STAR log file created - ${LOG_FILE}"
    echo "Key alignment metrics:"
    grep -E "(Uniquely mapped reads|Number of reads mapped to multiple loci|% of reads mapped)" "${LOG_FILE}"
else
    echo "WARNING: STAR log file not found!"
fi

# Create symlinks in main output directory for easier access
ln -sf ${SAMPLE_OUTPUT}/${SAMPLE}_Aligned.sortedByCoord.out.bam ${OUTPUT_DIR}/${SAMPLE}_sorted.bam
ln -sf ${SAMPLE_OUTPUT}/${SAMPLE}_Aligned.sortedByCoord.out.bam.bai ${OUTPUT_DIR}/${SAMPLE}_sorted.bam.bai
ln -sf ${SAMPLE_OUTPUT}/${SAMPLE}_ReadsPerGene.out.tab ${OUTPUT_DIR}/${SAMPLE}_counts.tab

echo "=== Alignment complete for ${SAMPLE} ==="
echo "Timestamp: $(date)"
