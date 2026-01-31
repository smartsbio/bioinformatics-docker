#!/bin/bash
set -e

# FMLRC Handler Script
# Handles long-read error correction using FM-index based methods

echo "🧬 FMLRC Handler Started"

# Get the input file (should be downloaded already)
# Use find to locate files in subdirectories (handles nested folder structures)
INPUT_FILE_PATH=$(find /tmp/input -type f | head -n 1)

if [[ -z "$INPUT_FILE_PATH" ]]; then
    echo "❌ No input files found in /tmp/input/"
    exit 1
fi

INPUT_FILENAME=$(basename "$INPUT_FILE_PATH")

echo "📁 Processing file: $INPUT_FILENAME"
echo "🎯 FMLRC command: $COMMAND"

# Extract organization and workspace IDs from INPUT_S3_KEY
# Format: organizations/{orgId}/workspaces/{workspaceId}/files/{path}
if [[ -n "$INPUT_S3_KEY" ]]; then
    ORGANIZATION_ID=$(echo "$INPUT_S3_KEY" | cut -d'/' -f2)
    WORKSPACE_ID=$(echo "$INPUT_S3_KEY" | cut -d'/' -f4)
    echo "📂 Organization ID: $ORGANIZATION_ID"
    echo "📂 Workspace ID: $WORKSPACE_ID"
fi

# Parse additional parameters from environment
OUTPUT_FILE=${OUTPUT_FILE:-"corrected_reads.fastq"}
# Strip @ notation from output file path
OUTPUT_FILE="${OUTPUT_FILE#@}"
KMER_SIZE=${KMER_SIZE:-"21"}
MIN_COUNT=${MIN_COUNT:-"2"}
CACHE_SIZE=${CACHE_SIZE:-"8"}
THREADS=${THREADS:-"4"}
SHORT_READS_FILE=${SHORT_READS_FILE:-""}

# Create subdirectories in output path if needed
OUTPUT_DIR=$(dirname "/tmp/output/$OUTPUT_FILE")
if [[ "$OUTPUT_DIR" != "/tmp/output" ]]; then
    mkdir -p "$OUTPUT_DIR"
    echo "📁 Created output directory: $OUTPUT_DIR"
fi

case "$COMMAND" in
    "correct")
        echo "🔧 Running FMLRC long-read correction..."
        
        if [[ -z "$SHORT_READS_FILE" ]]; then
            echo "❌ Short reads file required for error correction"
            exit 1
        fi
        
        # Copy input files
        cp "$INPUT_FILE_PATH" "/tmp/long_reads.fastq"
        # Note: SHORT_READS_FILE would be another input file in real implementation
        
        # Step 1: Convert short reads to FMLRC format
        echo "📊 Converting short reads to FMLRC format..."

        # Extract raw sequences (fmlrc-convert requires raw DNA bases only, no headers)
        if [[ "$SHORT_READS_FILE" == *.fastq ]] || [[ "$SHORT_READS_FILE" == *.fq ]]; then
            echo "📊 Extracting raw sequences from FASTQ..."
            # Use awk to extract only sequence lines (every 2nd line starting from line 2)
            if ! awk 'NR%4==2' "$SHORT_READS_FILE" > /tmp/short_reads_raw.txt; then
                echo "❌ Raw sequence extraction failed"
                exit 1
            fi
            SHORT_READS_INPUT="/tmp/short_reads_raw.txt"
        elif [[ "$SHORT_READS_FILE" == *.fa ]] || [[ "$SHORT_READS_FILE" == *.fasta ]]; then
            echo "📊 Extracting raw sequences from FASTA..."
            # Use grep to get only sequence lines (not header lines starting with >)
            if ! grep -v '^>' "$SHORT_READS_FILE" > /tmp/short_reads_raw.txt; then
                echo "❌ Raw sequence extraction failed"
                exit 1
            fi
            SHORT_READS_INPUT="/tmp/short_reads_raw.txt"
        else
            SHORT_READS_INPUT="$SHORT_READS_FILE"
        fi

        # fmlrc-convert syntax: fmlrc-convert -i input_file output_file
        CONVERT_CMD="fmlrc-convert -i $SHORT_READS_INPUT /tmp/short_reads.fmlrc"

        echo "🚀 Executing: $CONVERT_CMD"

        if ! eval "$CONVERT_CMD"; then
            echo "❌ FMLRC convert failed"
            exit 1
        fi

        echo "✅ Short reads converted successfully"
        
        # Step 2: Perform error correction
        echo "🔧 Performing long-read error correction..."
        # fmlrc syntax: fmlrc [options] <comp_msbwt.npy> <long_reads.fa> <corrected_reads.fa>
        # -k: k-mer size, -p: threads (not -t), -m: minimum count
        FMLRC_CMD="fmlrc -k $KMER_SIZE -p $THREADS"

        if [[ -n "$MIN_COUNT" ]]; then
            FMLRC_CMD="$FMLRC_CMD -m $MIN_COUNT"
        fi

        FMLRC_CMD="$FMLRC_CMD /tmp/short_reads.fmlrc /tmp/long_reads.fastq /tmp/output/$OUTPUT_FILE"

        echo "🚀 Executing: $FMLRC_CMD"

        if eval "$FMLRC_CMD"; then
            echo "✅ FMLRC error correction completed successfully"

            # Display output file info
            OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILE" 2>/dev/null || stat -f%z "/tmp/output/$OUTPUT_FILE" 2>/dev/null || echo "unknown")
            echo "📊 Corrected reads: $OUTPUT_FILE ($OUTPUT_SIZE bytes)"
        else
            echo "❌ FMLRC error correction failed"
            exit 1
        fi
        ;;
        
    "convert")
        echo "🔄 Running FMLRC convert (short reads preprocessing)..."

        # Copy input file
        cp "$INPUT_FILE_PATH" "/tmp/short_reads.fastq"

        # Extract raw sequences (fmlrc-convert requires raw DNA bases only, no headers)
        echo "📊 Extracting raw sequences from FASTQ..."
        # Use awk to extract only sequence lines (every 2nd line starting from line 2)
        if ! awk 'NR%4==2' /tmp/short_reads.fastq > /tmp/short_reads_raw.txt; then
            echo "❌ Raw sequence extraction failed"
            exit 1
        fi
        echo "✅ Extracted raw sequences"

        # fmlrc-convert syntax: fmlrc-convert -i input_file output_file
        CONVERT_CMD="fmlrc-convert -i /tmp/short_reads_raw.txt /tmp/output/$OUTPUT_FILE"

        echo "🚀 Executing: $CONVERT_CMD"

        if eval "$CONVERT_CMD"; then
            echo "✅ FMLRC convert completed successfully"

            # Display output file info
            OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILE" 2>/dev/null || stat -f%z "/tmp/output/$OUTPUT_FILE" 2>/dev/null || echo "unknown")
            echo "📊 FMLRC index: $OUTPUT_FILE ($OUTPUT_SIZE bytes)"
        else
            echo "❌ FMLRC convert failed"
            exit 1
        fi
        ;;
        
    "batch-correct")
        echo "🔧 Running FMLRC batch correction..."
        
        # For batch processing multiple long-read files
        cp "$INPUT_FILE_PATH" "/tmp/long_reads.fastq"
        
        if [[ -z "$SHORT_READS_FILE" ]]; then
            echo "❌ Short reads file required for batch correction"
            exit 1
        fi
        
        # Convert short reads first
        # Extract raw sequences (fmlrc-convert requires raw DNA bases only, no headers)
        if [[ "$SHORT_READS_FILE" == *.fastq ]] || [[ "$SHORT_READS_FILE" == *.fq ]]; then
            echo "📊 Extracting raw sequences from FASTQ..."
            # Use awk to extract only sequence lines (every 2nd line starting from line 2)
            if ! awk 'NR%4==2' "$SHORT_READS_FILE" > /tmp/short_reads_raw.txt; then
                echo "❌ Raw sequence extraction failed"
                exit 1
            fi
            SHORT_READS_INPUT="/tmp/short_reads_raw.txt"
        elif [[ "$SHORT_READS_FILE" == *.fa ]] || [[ "$SHORT_READS_FILE" == *.fasta ]]; then
            echo "📊 Extracting raw sequences from FASTA..."
            # Use grep to get only sequence lines (not header lines starting with >)
            if ! grep -v '^>' "$SHORT_READS_FILE" > /tmp/short_reads_raw.txt; then
                echo "❌ Raw sequence extraction failed"
                exit 1
            fi
            SHORT_READS_INPUT="/tmp/short_reads_raw.txt"
        else
            SHORT_READS_INPUT="$SHORT_READS_FILE"
        fi

        # fmlrc-convert syntax: fmlrc-convert -i input_file output_file
        CONVERT_CMD="fmlrc-convert -i $SHORT_READS_INPUT /tmp/short_reads.fmlrc"

        if eval "$CONVERT_CMD"; then
            echo "✅ Short reads converted for batch processing"
        else
            echo "❌ Batch convert failed"
            exit 1
        fi
        
        # Correct each long-read file
        BATCH_OUTPUT_DIR="/tmp/output/batch_corrected"
        mkdir -p "$BATCH_OUTPUT_DIR"

        # fmlrc syntax: fmlrc [options] <comp_msbwt.npy> <long_reads.fa> <corrected_reads.fa>
        FMLRC_CMD="fmlrc -k $KMER_SIZE -p $THREADS -m $MIN_COUNT"
        FMLRC_CMD="$FMLRC_CMD /tmp/short_reads.fmlrc /tmp/long_reads.fastq $BATCH_OUTPUT_DIR/$OUTPUT_FILE"

        if eval "$FMLRC_CMD"; then
            echo "✅ FMLRC batch correction completed successfully"

            echo "📁 Batch corrected files:"
            ls -la "$BATCH_OUTPUT_DIR/"
        else
            echo "❌ FMLRC batch correction failed"
            exit 1
        fi
        ;;
        
    *)
        echo "❌ Unsupported FMLRC command: $COMMAND"
        echo "Supported commands: correct, convert, batch-correct"
        exit 1
        ;;
esac

echo "🎯 FMLRC handler completed successfully"