#!/bin/bash
set -e

# GATK Handler Script
# Handles various GATK commands for variant calling and analysis

echo "🧬 GATK Handler Started"

# Get the input file (should be downloaded already)
# Use find to locate files in subdirectories (handles nested folder structures)
INPUT_FILE_PATH=$(find /tmp/input -type f | head -n 1)

if [[ -z "$INPUT_FILE_PATH" ]]; then
    echo "❌ No input files found in /tmp/input/"
    exit 1
fi

INPUT_FILENAME=$(basename "$INPUT_FILE_PATH")

echo "📁 Processing file: $INPUT_FILENAME"
echo "🎯 GATK command: $COMMAND"

# Parse additional parameters from environment
OUTPUT_FILE=${OUTPUT_FILE:-"output.vcf"}
REFERENCE_GENOME=${REFERENCE_GENOME:-""}
INTERVALS=${INTERVALS:-""}
THREADS=${THREADS:-"4"}
MEMORY=${MEMORY:-"16g"}

# Extract organization and workspace IDs from INPUT_S3_KEY
# Format: organizations/{orgId}/workspaces/{workspaceId}/files/{path}
if [[ -n "$INPUT_S3_KEY" ]]; then
    ORGANIZATION_ID=$(echo "$INPUT_S3_KEY" | cut -d'/' -f2)
    WORKSPACE_ID=$(echo "$INPUT_S3_KEY" | cut -d'/' -f4)
    echo "📂 Organization ID: $ORGANIZATION_ID"
    echo "📂 Workspace ID: $WORKSPACE_ID"
fi

# Find reference genome file
# If REFERENCE_GENOME uses @ notation, find or download it
if [[ "$REFERENCE_GENOME" == @* ]]; then
    REF_S3_KEY="${REFERENCE_GENOME:1}"  # Remove @ prefix
    REF_FILENAME=$(basename "$REF_S3_KEY")

    # Look for the file in /tmp/input (it may have been downloaded already)
    if [[ -f "/tmp/input/$REF_S3_KEY" ]]; then
        REFERENCE_GENOME="/tmp/input/$REF_S3_KEY"
        echo "✅ Found reference genome: $REFERENCE_GENOME"
    elif [[ -f "/tmp/input/$REF_FILENAME" ]]; then
        REFERENCE_GENOME="/tmp/input/$REF_FILENAME"
        echo "✅ Found reference genome: $REFERENCE_GENOME"
    else
        # Not found locally — download from S3 using the @ workspace notation
        echo "📥 Reference genome not in /tmp/input, downloading from S3..."
        if [[ -n "$ORGANIZATION_ID" && -n "$WORKSPACE_ID" ]]; then
            FULL_REF_S3_KEY="organizations/${ORGANIZATION_ID}/workspaces/${WORKSPACE_ID}/files/${REF_S3_KEY}"
            TARGET_DIR=$(dirname "/tmp/input/$REF_S3_KEY")
            mkdir -p "$TARGET_DIR"
            echo "📥 Downloading: s3://${S3_BUCKET}/${FULL_REF_S3_KEY} -> /tmp/input/${REF_S3_KEY}"
            if aws s3 cp "s3://${S3_BUCKET}/${FULL_REF_S3_KEY}" "/tmp/input/${REF_S3_KEY}" --no-progress; then
                REFERENCE_GENOME="/tmp/input/$REF_S3_KEY"
                echo "✅ Reference genome downloaded: $REFERENCE_GENOME"

                # Try to download pre-built .fai index (skips samtools faidx at runtime)
                if aws s3 cp "s3://${S3_BUCKET}/${FULL_REF_S3_KEY}.fai" "/tmp/input/${REF_S3_KEY}.fai" --no-progress 2>/dev/null; then
                    echo "✅ Reference .fai index downloaded"
                fi

                # Try to download pre-built .dict file (skips samtools dict at runtime)
                DICT_S3_KEY="${FULL_REF_S3_KEY%.*}.dict"
                DICT_LOCAL="${REF_S3_KEY%.*}.dict"
                if aws s3 cp "s3://${S3_BUCKET}/${DICT_S3_KEY}" "/tmp/input/${DICT_LOCAL}" --no-progress 2>/dev/null; then
                    echo "✅ Reference .dict file downloaded"
                fi
            else
                echo "❌ Failed to download reference genome from S3"
                exit 1
            fi
        else
            echo "❌ Reference genome not found and cannot determine S3 location"
            echo "Available files in /tmp/input:"
            ls -lah /tmp/input/
            exit 1
        fi
    fi
elif [[ -z "$REFERENCE_GENOME" ]]; then
    # If no reference specified, try to auto-detect .fasta/.fa files
    FASTA_FILES=(/tmp/input/*.fasta /tmp/input/*.fa /tmp/input/*/*.fasta /tmp/input/*/*.fa)
    for f in "${FASTA_FILES[@]}"; do
        if [[ -f "$f" ]]; then
            REFERENCE_GENOME="$f"
            echo "✅ Auto-detected reference genome: $REFERENCE_GENOME"
            break
        fi
    done
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

        # Create FASTA index if it doesn't exist
        if [[ ! -f "${REFERENCE_GENOME}.fai" ]]; then
            echo "🔧 Creating FASTA index for reference genome..."
            if ! samtools faidx "$REFERENCE_GENOME"; then
                echo "❌ Failed to create FASTA index"
                exit 1
            fi
            echo "✅ FASTA index created successfully"
        else
            echo "✅ FASTA index already exists"
        fi

        # Create sequence dictionary if it doesn't exist (.dict file)
        # Remove .fasta or .fa extension and add .dict
        DICT_FILE="${REFERENCE_GENOME%.fasta}"
        DICT_FILE="${DICT_FILE%.fa}.dict"
        if [[ ! -f "$DICT_FILE" ]]; then
            echo "🔧 Creating sequence dictionary for reference genome..."
            if ! samtools dict "$REFERENCE_GENOME" -o "$DICT_FILE"; then
                echo "⚠️ Failed to create sequence dictionary (continuing anyway)"
            else
                echo "✅ Sequence dictionary created successfully"
            fi
        else
            echo "✅ Sequence dictionary already exists"
        fi

        # Check if input is SAM or BAM, convert if necessary
        INPUT_BAM="/tmp/input.bam"
        if [[ "$INPUT_FILENAME" == *.sam ]]; then
            echo "📝 Converting SAM to BAM format..."

            # Check if SAM has read group information
            if samtools view -H "$INPUT_FILE_PATH" | grep -q "^@RG"; then
                echo "✅ SAM file has read group information"
                if ! samtools view -b -o "$INPUT_BAM" "$INPUT_FILE_PATH"; then
                    echo "❌ Failed to convert SAM to BAM"
                    exit 1
                fi
            else
                echo "⚠️ SAM file missing read groups, adding default read group..."
                # Add read group during conversion
                if ! samtools view -b "$INPUT_FILE_PATH" | \
                     samtools addreplacerg -r '@RG\tID:1\tSM:sample\tPL:ILLUMINA\tLB:lib1\tPU:unit1' -o "$INPUT_BAM" -; then
                    echo "❌ Failed to convert SAM to BAM with read groups"
                    exit 1
                fi
            fi
            echo "✅ SAM to BAM conversion completed"
        elif [[ "$INPUT_FILENAME" == *.bam ]]; then
            echo "📝 Input is already in BAM format, copying..."
            cp "$INPUT_FILE_PATH" "$INPUT_BAM"

            # Check if BAM has read group information (required by GATK HaplotypeCaller)
            if ! samtools view -H "$INPUT_BAM" | grep -q "^@RG"; then
                echo "⚠️ BAM file missing read groups, adding default read group..."
                TEMP_BAM="/tmp/input_with_rg.bam"
                if ! samtools addreplacerg -r '@RG\tID:1\tSM:sample\tPL:ILLUMINA\tLB:lib1\tPU:unit1' -o "$TEMP_BAM" "$INPUT_BAM"; then
                    echo "❌ Failed to add read groups to BAM"
                    exit 1
                fi
                mv "$TEMP_BAM" "$INPUT_BAM"
                echo "✅ Read groups added successfully"
            else
                echo "✅ BAM file has read group information"
            fi
        else
            echo "⚠️ Unknown input format (expected .sam or .bam), attempting to use as-is..."
            cp "$INPUT_FILE_PATH" "$INPUT_BAM"
        fi

        # Create BAM index
        if [[ -f "${INPUT_BAM}.bai" ]]; then
            echo "✅ BAM index already exists"
        else
            echo "🔧 Creating BAM index with samtools..."
            if ! samtools index "$INPUT_BAM"; then
                echo "❌ Failed to create BAM index"
                exit 1
            fi
            echo "✅ BAM index created successfully"
        fi

        GATK_CMD="gatk HaplotypeCaller"

        # Add input and reference
        GATK_CMD="$GATK_CMD -I $INPUT_BAM"
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
