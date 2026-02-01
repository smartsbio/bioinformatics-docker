#!/bin/bash
set -e

# HISAT2 Handler Script
# Handles RNA-seq read alignment using HISAT2

echo "🧬 HISAT2 Handler Started"

# Get the input file (should be downloaded already)
# Use find to locate files in subdirectories (handles nested folder structures)
INPUT_FILE_PATH=$(find /tmp/input -type f | head -n 1)

if [[ -z "$INPUT_FILE_PATH" ]]; then
    echo "❌ No input files found in /tmp/input/"
    exit 1
fi

INPUT_FILENAME=$(basename "$INPUT_FILE_PATH")

echo "📁 Processing file: $INPUT_FILENAME"
echo "🎯 HISAT2 command: $COMMAND"

# Extract organization and workspace IDs from INPUT_S3_KEY
# Format: organizations/{orgId}/workspaces/{workspaceId}/files/{path}
if [[ -n "$INPUT_S3_KEY" ]]; then
    ORGANIZATION_ID=$(echo "$INPUT_S3_KEY" | cut -d'/' -f2)
    WORKSPACE_ID=$(echo "$INPUT_S3_KEY" | cut -d'/' -f4)
    echo "📂 Organization ID: $ORGANIZATION_ID"
    echo "📂 Workspace ID: $WORKSPACE_ID"
fi

# Parse additional parameters from environment
INDEX_PREFIX=${INDEX_PREFIX:-"reference_index"}
THREADS=${THREADS:-"4"}
PHRED_QUALITY=${PHRED_QUALITY:-"phred33"}
STRANDNESS=${STRANDNESS:-""}
MAX_INTRON_LENGTH=${MAX_INTRON_LENGTH:-"500000"}
MIN_INTRON_LENGTH=${MIN_INTRON_LENGTH:-"20"}
REFERENCE_GENOME=${REFERENCE_GENOME:-""}

# Set command-specific default output filenames
if [[ -z "$OUTPUT_FILE" ]]; then
    case "$COMMAND" in
        "extract-splicesites")
            OUTPUT_FILE="splicesites.txt"
            ;;
        "extract-exons")
            OUTPUT_FILE="exons.txt"
            ;;
        "build")
            OUTPUT_FILE="index"  # Not used, but set for consistency
            ;;
        *)
            OUTPUT_FILE="aligned_reads.sam"
            ;;
    esac
fi

# Strip @ notation from output file path
OUTPUT_FILE="${OUTPUT_FILE#@}"

case "$COMMAND" in
    "build")
        echo "📇 Running HISAT2 index building..."

        # Find reference genome file (prefer .fasta or .fa files)
        REFERENCE_FILE=$(find /tmp/input -type f \( -name "*.fasta" -o -name "*.fa" \) | head -n 1)

        if [[ -z "$REFERENCE_FILE" ]]; then
            echo "⚠️  No .fasta or .fa file found, using first available file"
            REFERENCE_FILE="$INPUT_FILE_PATH"
        fi

        echo "📚 Using reference file: $(basename $REFERENCE_FILE)"

        # Copy input file (reference genome)
        cp "$REFERENCE_FILE" "/tmp/reference.fasta"

        HISAT2_BUILD_CMD="hisat2-build"
        HISAT2_BUILD_CMD="$HISAT2_BUILD_CMD -p $THREADS"
        HISAT2_BUILD_CMD="$HISAT2_BUILD_CMD /tmp/reference.fasta"
        HISAT2_BUILD_CMD="$HISAT2_BUILD_CMD /tmp/output/$INDEX_PREFIX"

        echo "🚀 Executing: $HISAT2_BUILD_CMD"

        if eval "$HISAT2_BUILD_CMD"; then
            echo "✅ HISAT2 index build completed successfully"

            # Display generated index files
            echo "📁 Generated index files:"
            ls -la /tmp/output/${INDEX_PREFIX}* | awk '{print $9, $5 " bytes"}'
        else
            echo "❌ HISAT2 index build failed"
            exit 1
        fi
        ;;

    "align")
        echo "🧬 Running HISAT2 read alignment..."

        if [[ -z "$REFERENCE_GENOME" ]]; then
            echo "❌ Reference genome required for alignment"
            exit 1
        fi

        # Strip @ notation and file extension from reference genome
        # Example: @exp1/phix174.fasta -> exp1/phix174
        CLEAN_REFERENCE="${REFERENCE_GENOME#@}"  # Remove @ prefix
        CLEAN_REFERENCE="${CLEAN_REFERENCE%.*}"   # Remove file extension
        echo "📚 Using reference index: $CLEAN_REFERENCE"

        # Check if HISAT2 index exists, if not build it automatically
        INDEX_BASE="/tmp/index/$(basename $CLEAN_REFERENCE)"
        if [[ ! -f "${INDEX_BASE}.1.ht2" ]]; then
            echo "⚠️  HISAT2 index not found, building it automatically..."

            # Download reference genome FASTA file
            REFERENCE_FASTA="${REFERENCE_GENOME#@}"  # Remove @ prefix (keep extension)
            REFERENCE_S3_KEY="organizations/${ORGANIZATION_ID}/workspaces/${WORKSPACE_ID}/files/${REFERENCE_FASTA}"

            echo "📥 Downloading reference genome: s3://${S3_BUCKET}/${REFERENCE_S3_KEY}"
            mkdir -p /tmp/index

            if aws s3 cp "s3://${S3_BUCKET}/${REFERENCE_S3_KEY}" "/tmp/index/reference.fasta" --no-progress; then
                echo "✅ Reference genome downloaded"

                # Build the index
                echo "🔨 Building HISAT2 index..."
                BUILD_CMD="hisat2-build -p $THREADS /tmp/index/reference.fasta ${INDEX_BASE}"

                echo "🚀 Executing: $BUILD_CMD"
                if eval "$BUILD_CMD"; then
                    echo "✅ HISAT2 index built successfully"
                else
                    echo "❌ Failed to build HISAT2 index"
                    exit 1
                fi
            else
                echo "❌ Failed to download reference genome from S3"
                exit 1
            fi
        else
            echo "✅ HISAT2 index already exists"
        fi

        # Copy input file (reads)
        cp "$INPUT_FILE_PATH" "/tmp/reads.fastq"

        HISAT2_CMD="hisat2"
        HISAT2_CMD="$HISAT2_CMD -p $THREADS"

        # Add quality encoding
        if [[ "$PHRED_QUALITY" == "phred33" ]]; then
            HISAT2_CMD="$HISAT2_CMD --phred33"
        elif [[ "$PHRED_QUALITY" == "phred64" ]]; then
            HISAT2_CMD="$HISAT2_CMD --phred64"
        fi

        # Add strandness if specified
        if [[ -n "$STRANDNESS" ]]; then
            HISAT2_CMD="$HISAT2_CMD --rna-strandness $STRANDNESS"
        fi

        # Add intron length parameters
        HISAT2_CMD="$HISAT2_CMD --min-intronlen $MIN_INTRON_LENGTH"
        HISAT2_CMD="$HISAT2_CMD --max-intronlen $MAX_INTRON_LENGTH"

        # Specify index and reads (use the built index path)
        HISAT2_CMD="$HISAT2_CMD -x ${INDEX_BASE}"
        HISAT2_CMD="$HISAT2_CMD -U /tmp/reads.fastq"

        # Create subdirectories in output path if needed
        OUTPUT_DIR=$(dirname "/tmp/output/$OUTPUT_FILE")
        if [[ "$OUTPUT_DIR" != "/tmp/output" ]]; then
            mkdir -p "$OUTPUT_DIR"
            echo "📁 Created output directory: $OUTPUT_DIR"
        fi

        HISAT2_CMD="$HISAT2_CMD -S /tmp/output/$OUTPUT_FILE"

        echo "🚀 Executing: $HISAT2_CMD"

        if eval "$HISAT2_CMD"; then
            echo "✅ HISAT2 read alignment completed successfully"

            # Display output file info
            OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILE" 2>/dev/null || stat -f%z "/tmp/output/$OUTPUT_FILE" 2>/dev/null || echo "unknown")
            echo "📊 Aligned reads: $OUTPUT_FILE ($OUTPUT_SIZE bytes)"
        else
            echo "❌ HISAT2 alignment failed"
            exit 1
        fi
        ;;

    "align-paired")
        echo "🧬 Running HISAT2 paired-end alignment..."

        # For paired-end alignment, we need two input files
        # This is simplified - in real scenario would handle multiple input files
        cp "$INPUT_FILE_PATH" "/tmp/reads_1.fastq"

        # Second file would come from additional input
        READ2_FILE=${READ2_FILE:-""}
        if [[ -z "$READ2_FILE" ]]; then
            echo "⚠️ Note: Second read file not specified, using single-end mode"
            # Fall back to single-end alignment
            HISAT2_CMD="hisat2 -p $THREADS --phred33"
            HISAT2_CMD="$HISAT2_CMD -x /tmp/output/$INDEX_PREFIX"
            HISAT2_CMD="$HISAT2_CMD -U /tmp/reads_1.fastq"
            HISAT2_CMD="$HISAT2_CMD -S /tmp/output/$OUTPUT_FILE"
        else
            HISAT2_CMD="hisat2 -p $THREADS --phred33"

            # Add strandness if specified
            if [[ -n "$STRANDNESS" ]]; then
                HISAT2_CMD="$HISAT2_CMD --rna-strandness $STRANDNESS"
            fi

            HISAT2_CMD="$HISAT2_CMD -x /tmp/output/$INDEX_PREFIX"
            HISAT2_CMD="$HISAT2_CMD -1 /tmp/reads_1.fastq"
            HISAT2_CMD="$HISAT2_CMD -2 $READ2_FILE"
            HISAT2_CMD="$HISAT2_CMD -S /tmp/output/$OUTPUT_FILE"
        fi

        echo "🚀 Executing: $HISAT2_CMD"

        if eval "$HISAT2_CMD"; then
            echo "✅ HISAT2 paired-end alignment completed successfully"

            OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILE" 2>/dev/null || stat -f%z "/tmp/output/$OUTPUT_FILE" 2>/dev/null || echo "unknown")
            echo "📊 Aligned reads: $OUTPUT_FILE ($OUTPUT_SIZE bytes)"
        else
            echo "❌ HISAT2 paired-end alignment failed"
            exit 1
        fi
        ;;

    "extract-splicesites")
        echo "🧬 Running HISAT2 splice site extraction..."

        # Copy input file (GTF/GFF annotation)
        cp "$INPUT_FILE_PATH" "/tmp/annotation.gtf"

        EXTRACT_CMD="hisat2_extract_splice_sites.py"
        EXTRACT_CMD="$EXTRACT_CMD /tmp/annotation.gtf"
        EXTRACT_CMD="$EXTRACT_CMD > /tmp/output/$OUTPUT_FILE"

        echo "🚀 Executing: $EXTRACT_CMD"

        if eval "$EXTRACT_CMD"; then
            echo "✅ Splice site extraction completed successfully"

            OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILE" 2>/dev/null || stat -f%z "/tmp/output/$OUTPUT_FILE" 2>/dev/null || echo "unknown")
            echo "📊 Splice sites: $OUTPUT_FILE ($OUTPUT_SIZE bytes)"
        else
            echo "❌ Splice site extraction failed"
            exit 1
        fi
        ;;

    "extract-exons")
        echo "🧬 Running HISAT2 exon extraction..."

        # Copy input file (GTF/GFF annotation)
        cp "$INPUT_FILE_PATH" "/tmp/annotation.gtf"

        EXTRACT_CMD="hisat2_extract_exons.py"
        EXTRACT_CMD="$EXTRACT_CMD /tmp/annotation.gtf"
        EXTRACT_CMD="$EXTRACT_CMD > /tmp/output/$OUTPUT_FILE"

        echo "🚀 Executing: $EXTRACT_CMD"

        if eval "$EXTRACT_CMD"; then
            echo "✅ Exon extraction completed successfully"

            OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILE" 2>/dev/null || stat -f%z "/tmp/output/$OUTPUT_FILE" 2>/dev/null || echo "unknown")
            echo "📊 Exons: $OUTPUT_FILE ($OUTPUT_SIZE bytes)"
        else
            echo "❌ Exon extraction failed"
            exit 1
        fi
        ;;

    *)
        echo "❌ Unsupported HISAT2 command: $COMMAND"
        echo "Supported commands: build, align, align-paired, extract-splicesites, extract-exons"
        exit 1
        ;;
esac

echo "🎯 HISAT2 handler completed successfully"
