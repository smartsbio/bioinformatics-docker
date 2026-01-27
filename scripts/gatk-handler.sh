#!/bin/bash
set -e

# GATK Handler Script
# Handles various GATK commands for variant calling and analysis

echo "🧬 GATK Handler Started"

# Get the input file (should be downloaded already)
INPUT_FILES=(/tmp/input/*)
if [[ ${#INPUT_FILES[@]} -eq 0 ]]; then
    echo "❌ No input files found in /tmp/input/"
    exit 1
fi

INPUT_FILE_PATH="${INPUT_FILES[0]}"
INPUT_FILENAME=$(basename "$INPUT_FILE_PATH")

echo "📁 Processing file: $INPUT_FILENAME"
echo "🎯 GATK command: $COMMAND"

# Parse additional parameters from environment
OUTPUT_FILE=${OUTPUT_FILE:-"output.vcf"}
REFERENCE_GENOME=${REFERENCE_GENOME:-""}
INTERVALS=${INTERVALS:-""}
THREADS=${THREADS:-"4"}
MEMORY=${MEMORY:-"16g"}

# Download reference genome if it uses @ notation
if [[ "$REFERENCE_GENOME" == @* ]]; then
    REF_S3_KEY="${REFERENCE_GENOME:1}"  # Remove @ prefix
    REF_FILENAME=$(basename "$REF_S3_KEY")
    REF_LOCAL_PATH="/tmp/reference/$REF_FILENAME"

    echo "📥 Downloading reference genome from S3: $REF_S3_KEY"
    mkdir -p /tmp/reference

    if ! aws s3 cp "s3://$S3_BUCKET/test-data/input/$REF_S3_KEY" "$REF_LOCAL_PATH" --no-progress; then
        echo "❌ Failed to download reference genome from S3"
        exit 1
    fi

    # Download reference index files (.fai and .dict)
    aws s3 cp "s3://$S3_BUCKET/test-data/input/${REF_S3_KEY}.fai" "${REF_LOCAL_PATH}.fai" --no-progress 2>/dev/null || echo "⚠️ No .fai index found"
    aws s3 cp "s3://$S3_BUCKET/test-data/input/${REF_S3_KEY%.fasta}.dict" "${REF_LOCAL_PATH%.fasta}.dict" --no-progress 2>/dev/null || echo "⚠️ No .dict file found"

    echo "✅ Reference genome downloaded: $REF_LOCAL_PATH"
    REFERENCE_GENOME="$REF_LOCAL_PATH"
fi

# GATK-specific parameters
MIN_BASE_QUALITY=${MIN_BASE_QUALITY:-"20"}
MIN_MAPPING_QUALITY=${MIN_MAPPING_QUALITY:-"20"}
STAND_CALL_CONF=${STAND_CALL_CONF:-"30"}
MAX_READS_PER_ALIGNMENT=${MAX_READS_PER_ALIGNMENT:-"10000"}

case "$COMMAND" in
    "HaplotypeCaller")
        echo "🔍 Running GATK HaplotypeCaller..."

        if [[ -z "$REFERENCE_GENOME" ]]; then
            echo "❌ Reference genome required for HaplotypeCaller"
            exit 1
        fi

        # Copy input BAM file
        cp "$INPUT_FILE_PATH" "/tmp/input.bam"

        GATK_CMD="gatk HaplotypeCaller"

        # Add input and reference
        GATK_CMD="$GATK_CMD -I /tmp/input.bam"
        GATK_CMD="$GATK_CMD -R $REFERENCE_GENOME"
        GATK_CMD="$GATK_CMD -O /tmp/output/$OUTPUT_FILE"

        # Add intervals if specified
        if [[ -n "$INTERVALS" ]]; then
            GATK_CMD="$GATK_CMD -L $INTERVALS"
        fi

        # Add quality parameters
        if [[ -n "$MIN_BASE_QUALITY" ]]; then
            GATK_CMD="$GATK_CMD --min-base-quality-score $MIN_BASE_QUALITY"
        fi

        if [[ -n "$STAND_CALL_CONF" ]]; then
            GATK_CMD="$GATK_CMD --standard-min-confidence-threshold-for-calling $STAND_CALL_CONF"
        fi

        # Add parallelization
        GATK_CMD="$GATK_CMD --native-pair-hmm-threads $THREADS"

        echo "🚀 Executing: $GATK_CMD"

        # Execute the command
        if eval "$GATK_CMD"; then
            echo "✅ GATK HaplotypeCaller completed successfully"

            # Display output file info
            OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILE" 2>/dev/null || stat -f%z "/tmp/output/$OUTPUT_FILE" 2>/dev/null || echo "unknown")
            echo "📊 Output file: $OUTPUT_FILE ($OUTPUT_SIZE bytes)"
        else
            echo "❌ GATK HaplotypeCaller failed"
            exit 1
        fi
        ;;

    "Mutect2")
        echo "🔍 Running GATK Mutect2 (somatic variant calling)..."

        if [[ -z "$REFERENCE_GENOME" ]]; then
            echo "❌ Reference genome required for Mutect2"
            exit 1
        fi

        # Copy input BAM file (tumor sample)
        cp "$INPUT_FILE_PATH" "/tmp/tumor.bam"

        GATK_CMD="gatk Mutect2"

        # Add input and reference
        GATK_CMD="$GATK_CMD -I /tmp/tumor.bam"
        GATK_CMD="$GATK_CMD -R $REFERENCE_GENOME"
        GATK_CMD="$GATK_CMD -O /tmp/output/$OUTPUT_FILE"

        # Add intervals if specified
        if [[ -n "$INTERVALS" ]]; then
            GATK_CMD="$GATK_CMD -L $INTERVALS"
        fi

        # Add parallelization
        GATK_CMD="$GATK_CMD --native-pair-hmm-threads $THREADS"

        echo "🚀 Executing: $GATK_CMD"

        # Execute the command
        if eval "$GATK_CMD"; then
            echo "✅ GATK Mutect2 completed successfully"

            # Display output file info
            OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILE" 2>/dev/null || stat -f%z "/tmp/output/$OUTPUT_FILE" 2>/dev/null || echo "unknown")
            echo "📊 Output file: $OUTPUT_FILE ($OUTPUT_SIZE bytes)"
        else
            echo "❌ GATK Mutect2 failed"
            exit 1
        fi
        ;;

    "GenotypeGVCFs")
        echo "🔍 Running GATK GenotypeGVCFs..."

        if [[ -z "$REFERENCE_GENOME" ]]; then
            echo "❌ Reference genome required for GenotypeGVCFs"
            exit 1
        fi

        # Copy input GVCF file
        cp "$INPUT_FILE_PATH" "/tmp/input.g.vcf"

        GATK_CMD="gatk GenotypeGVCFs"

        # Add input and reference
        GATK_CMD="$GATK_CMD -V /tmp/input.g.vcf"
        GATK_CMD="$GATK_CMD -R $REFERENCE_GENOME"
        GATK_CMD="$GATK_CMD -O /tmp/output/$OUTPUT_FILE"

        # Add intervals if specified
        if [[ -n "$INTERVALS" ]]; then
            GATK_CMD="$GATK_CMD -L $INTERVALS"
        fi

        echo "🚀 Executing: $GATK_CMD"

        # Execute the command
        if eval "$GATK_CMD"; then
            echo "✅ GATK GenotypeGVCFs completed successfully"

            # Display output file info
            OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILE" 2>/dev/null || stat -f%z "/tmp/output/$OUTPUT_FILE" 2>/dev/null || echo "unknown")
            echo "📊 Output file: $OUTPUT_FILE ($OUTPUT_SIZE bytes)"
        else
            echo "❌ GATK GenotypeGVCFs failed"
            exit 1
        fi
        ;;

    "VariantFiltration")
        echo "🔍 Running GATK VariantFiltration..."

        if [[ -z "$REFERENCE_GENOME" ]]; then
            echo "❌ Reference genome required for VariantFiltration"
            exit 1
        fi

        # Copy input VCF file
        cp "$INPUT_FILE_PATH" "/tmp/input.vcf"

        GATK_CMD="gatk VariantFiltration"

        # Add input and reference
        GATK_CMD="$GATK_CMD -V /tmp/input.vcf"
        GATK_CMD="$GATK_CMD -R $REFERENCE_GENOME"
        GATK_CMD="$GATK_CMD -O /tmp/output/$OUTPUT_FILE"

        # Add default filter expressions
        GATK_CMD="$GATK_CMD --filter-expression 'QD < 2.0' --filter-name 'QD2'"
        GATK_CMD="$GATK_CMD --filter-expression 'FS > 60.0' --filter-name 'FS60'"
        GATK_CMD="$GATK_CMD --filter-expression 'MQ < 40.0' --filter-name 'MQ40'"

        echo "🚀 Executing: $GATK_CMD"

        # Execute the command
        if eval "$GATK_CMD"; then
            echo "✅ GATK VariantFiltration completed successfully"

            # Display output file info
            OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILE" 2>/dev/null || stat -f%z "/tmp/output/$OUTPUT_FILE" 2>/dev/null || echo "unknown")
            echo "📊 Output file: $OUTPUT_FILE ($OUTPUT_SIZE bytes)"
        else
            echo "❌ GATK VariantFiltration failed"
            exit 1
        fi
        ;;

    "SelectVariants")
        echo "🔍 Running GATK SelectVariants..."

        if [[ -z "$REFERENCE_GENOME" ]]; then
            echo "❌ Reference genome required for SelectVariants"
            exit 1
        fi

        # Copy input VCF file
        cp "$INPUT_FILE_PATH" "/tmp/input.vcf"

        GATK_CMD="gatk SelectVariants"

        # Add input and reference
        GATK_CMD="$GATK_CMD -V /tmp/input.vcf"
        GATK_CMD="$GATK_CMD -R $REFERENCE_GENOME"
        GATK_CMD="$GATK_CMD -O /tmp/output/$OUTPUT_FILE"

        # Add intervals if specified
        if [[ -n "$INTERVALS" ]]; then
            GATK_CMD="$GATK_CMD -L $INTERVALS"
        fi

        # Select only SNPs by default
        GATK_CMD="$GATK_CMD --select-type-to-include SNP"

        echo "🚀 Executing: $GATK_CMD"

        # Execute the command
        if eval "$GATK_CMD"; then
            echo "✅ GATK SelectVariants completed successfully"

            # Display output file info
            OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILE" 2>/dev/null || stat -f%z "/tmp/output/$OUTPUT_FILE" 2>/dev/null || echo "unknown")
            echo "📊 Output file: $OUTPUT_FILE ($OUTPUT_SIZE bytes)"
        else
            echo "❌ GATK SelectVariants failed"
            exit 1
        fi
        ;;

    *)
        echo "❌ Unsupported GATK command: $COMMAND"
        echo "Supported commands: HaplotypeCaller, Mutect2, GenotypeGVCFs, VariantFiltration, SelectVariants"
        exit 1
        ;;
esac

echo "🎯 GATK handler completed successfully"
