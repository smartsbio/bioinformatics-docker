#!/bin/bash
set -e

# HISAT2 Handler Script
# Handles RNA-seq read alignment using HISAT2

echo "🧬 HISAT2 Handler Started"

# Get the input file (should be downloaded already)
INPUT_FILES=(/tmp/input/*)
if [[ ${#INPUT_FILES[@]} -eq 0 ]]; then
    echo "❌ No input files found in /tmp/input/"
    exit 1
fi

INPUT_FILE_PATH="${INPUT_FILES[0]}"
INPUT_FILENAME=$(basename "$INPUT_FILE_PATH")

echo "📁 Processing file: $INPUT_FILENAME"
echo "🎯 HISAT2 command: $COMMAND"

# Parse additional parameters from environment
OUTPUT_FILE=${OUTPUT_FILE:-"aligned_reads.sam"}
INDEX_PREFIX=${INDEX_PREFIX:-"reference_index"}
THREADS=${THREADS:-"4"}
PHRED_QUALITY=${PHRED_QUALITY:-"phred33"}
STRANDNESS=${STRANDNESS:-""}
MAX_INTRON_LENGTH=${MAX_INTRON_LENGTH:-"500000"}
MIN_INTRON_LENGTH=${MIN_INTRON_LENGTH:-"20"}

case "$COMMAND" in
    "build")
        echo "📇 Running HISAT2 index building..."

        # Copy input file (reference genome)
        cp "$INPUT_FILE_PATH" "/tmp/reference.fasta"

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

        # Copy input file (reads)
        cp "$INPUT_FILE_PATH" "/tmp/reads.fastq"

        # Check if index files are provided (in real scenario)
        if [[ -z "$INDEX_PREFIX" ]]; then
            echo "❌ Index prefix required for alignment"
            exit 1
        fi

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

        # Specify index and reads
        HISAT2_CMD="$HISAT2_CMD -x /tmp/output/$INDEX_PREFIX"
        HISAT2_CMD="$HISAT2_CMD -U /tmp/reads.fastq"
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
        EXTRACT_CMD="$EXTRACT_CMD > /tmp/output/splicesites.txt"

        echo "🚀 Executing: $EXTRACT_CMD"

        if eval "$EXTRACT_CMD"; then
            echo "✅ Splice site extraction completed successfully"

            OUTPUT_SIZE=$(stat -c%s "/tmp/output/splicesites.txt" 2>/dev/null || stat -f%z "/tmp/output/splicesites.txt" 2>/dev/null || echo "unknown")
            echo "📊 Splice sites: splicesites.txt ($OUTPUT_SIZE bytes)"
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
        EXTRACT_CMD="$EXTRACT_CMD > /tmp/output/exons.txt"

        echo "🚀 Executing: $EXTRACT_CMD"

        if eval "$EXTRACT_CMD"; then
            echo "✅ Exon extraction completed successfully"

            OUTPUT_SIZE=$(stat -c%s "/tmp/output/exons.txt" 2>/dev/null || stat -f%z "/tmp/output/exons.txt" 2>/dev/null || echo "unknown")
            echo "📊 Exons: exons.txt ($OUTPUT_SIZE bytes)"
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
