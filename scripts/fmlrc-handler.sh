#!/bin/bash
set -e

# FMLRC Handler Script
# Handles long-read error correction using FM-index based methods

echo "🧬 FMLRC Handler Started"

# Get the input file (should be downloaded already)
INPUT_FILES=(/tmp/input/*)
if [[ ${#INPUT_FILES[@]} -eq 0 ]]; then
    echo "❌ No input files found in /tmp/input/"
    exit 1
fi

INPUT_FILE_PATH="${INPUT_FILES[0]}"
INPUT_FILENAME=$(basename "$INPUT_FILE_PATH")

echo "📁 Processing file: $INPUT_FILENAME"
echo "🎯 FMLRC command: $COMMAND"

# Parse additional parameters from environment
OUTPUT_FILE=${OUTPUT_FILE:-"corrected_reads.fastq"}
KMER_SIZE=${KMER_SIZE:-"21"}
MIN_COUNT=${MIN_COUNT:-"2"}
CACHE_SIZE=${CACHE_SIZE:-"8"}
THREADS=${THREADS:-"4"}
SHORT_READS_FILE=${SHORT_READS_FILE:-""}

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
        CONVERT_CMD="fmlrc-convert"
        CONVERT_CMD="$CONVERT_CMD -k $KMER_SIZE"
        CONVERT_CMD="$CONVERT_CMD -t $THREADS"
        CONVERT_CMD="$CONVERT_CMD $SHORT_READS_FILE /tmp/short_reads.fmlrc"
        
        echo "🚀 Executing: $CONVERT_CMD"
        
        if ! eval "$CONVERT_CMD"; then
            echo "❌ FMLRC convert failed"
            exit 1
        fi
        
        echo "✅ Short reads converted successfully"
        
        # Step 2: Perform error correction
        echo "🔧 Performing long-read error correction..."
        FMLRC_CMD="fmlrc"
        FMLRC_CMD="$FMLRC_CMD -k $KMER_SIZE"
        FMLRC_CMD="$FMLRC_CMD -t $THREADS"
        FMLRC_CMD="$FMLRC_CMD -C $CACHE_SIZE"
        
        if [[ -n "$MIN_COUNT" ]]; then
            FMLRC_CMD="$FMLRC_CMD -m $MIN_COUNT"
        fi
        
        FMLRC_CMD="$FMLRC_CMD /tmp/short_reads.fmlrc /tmp/long_reads.fastq /tmp/output/$OUTPUT_FILE"
        
        echo "🚀 Executing: $FMLRC_CMD"
        
        if eval "$FMLRC_CMD"; then
            echo "✅ FMLRC error correction completed successfully"
            
            # Display output file info
            OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILE" 2>/dev/null || echo "unknown")
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
        
        CONVERT_CMD="fmlrc-convert"
        CONVERT_CMD="$CONVERT_CMD -k $KMER_SIZE"
        CONVERT_CMD="$CONVERT_CMD -t $THREADS"
        CONVERT_CMD="$CONVERT_CMD /tmp/short_reads.fastq /tmp/output/$OUTPUT_FILE"
        
        echo "🚀 Executing: $CONVERT_CMD"
        
        if eval "$CONVERT_CMD"; then
            echo "✅ FMLRC convert completed successfully"
            
            # Display output file info
            OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILE" 2>/dev/null || echo "unknown")
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
        CONVERT_CMD="fmlrc-convert -k $KMER_SIZE -t $THREADS $SHORT_READS_FILE /tmp/short_reads.fmlrc"
        
        if eval "$CONVERT_CMD"; then
            echo "✅ Short reads converted for batch processing"
        else
            echo "❌ Batch convert failed"
            exit 1
        fi
        
        # Correct each long-read file
        BATCH_OUTPUT_DIR="/tmp/output/batch_corrected"
        mkdir -p "$BATCH_OUTPUT_DIR"
        
        FMLRC_CMD="fmlrc -k $KMER_SIZE -t $THREADS -C $CACHE_SIZE"
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