#!/bin/bash
set -e

# Trimmomatic Handler Script
# Handles adapter trimming and quality filtering for sequencing data

echo "🧬 Trimmomatic Handler Started"

# Get the input file (should be downloaded already)
# Use find to locate files in subdirectories (handles nested folder structures)
INPUT_FILE_PATH=$(find /tmp/input -type f | head -n 1)

if [[ -z "$INPUT_FILE_PATH" ]]; then
    echo "❌ No input files found in /tmp/input/"
    exit 1
fi

INPUT_FILENAME=$(basename "$INPUT_FILE_PATH")

echo "📁 Processing file: $INPUT_FILENAME"
echo "🎯 Trimmomatic command: $COMMAND"

# Parse additional parameters from environment
OUTPUT_FILE=${OUTPUT_FILE:-"trimmed_reads.fastq"}
THREADS=${THREADS:-"4"}
QUALITY_THRESHOLD=${QUALITY_THRESHOLD:-"20"}
MIN_LENGTH=${MIN_LENGTH:-"36"}
LEADING=${LEADING:-"3"}
TRAILING=${TRAILING:-"3"}
SLIDING_WINDOW=${SLIDING_WINDOW:-"4:15"}
ADAPTER_FILE=${ADAPTER_FILE:-""}

case "$COMMAND" in
    "single-end")
        echo "✂️ Running Trimmomatic single-end trimming..."

        # Use the input file directly (trimmomatic handles .fastq.gz natively)
        OUTPUT_FILE="${OUTPUT_FILE:-"trimmed_reads.fastq.gz"}"

        TRIMMOMATIC_CMD="trimmomatic SE -phred33"

        # Add threads parameter
        if [[ -n "$THREADS" ]]; then
            TRIMMOMATIC_CMD="$TRIMMOMATIC_CMD -threads $THREADS"
        fi

        # Add input and output files
        TRIMMOMATIC_CMD="$TRIMMOMATIC_CMD $INPUT_FILE_PATH /tmp/output/$OUTPUT_FILE"
        
        # Add adapter trimming if specified
        if [[ -n "$ADAPTER_FILE" ]]; then
            TRIMMOMATIC_CMD="$TRIMMOMATIC_CMD ILLUMINACLIP:$ADAPTER_FILE:2:30:10"
        fi
        
        # Add quality trimming parameters
        if [[ -n "$LEADING" ]]; then
            TRIMMOMATIC_CMD="$TRIMMOMATIC_CMD LEADING:$LEADING"
        fi
        
        if [[ -n "$TRAILING" ]]; then
            TRIMMOMATIC_CMD="$TRIMMOMATIC_CMD TRAILING:$TRAILING"
        fi
        
        if [[ -n "$SLIDING_WINDOW" ]]; then
            TRIMMOMATIC_CMD="$TRIMMOMATIC_CMD SLIDINGWINDOW:$SLIDING_WINDOW"
        fi
        
        if [[ -n "$MIN_LENGTH" ]]; then
            TRIMMOMATIC_CMD="$TRIMMOMATIC_CMD MINLEN:$MIN_LENGTH"
        fi
        
        echo "🚀 Executing: $TRIMMOMATIC_CMD"
        
        # Execute the command
        if eval "$TRIMMOMATIC_CMD"; then
            echo "✅ Trimmomatic single-end trimming completed successfully"
            
            # Display output file info
            OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILE" 2>/dev/null || echo "unknown")
            echo "📊 Output file: $OUTPUT_FILE ($OUTPUT_SIZE bytes)"
        else
            echo "❌ Trimmomatic single-end trimming failed"
            exit 1
        fi
        ;;
        
    "paired-end")
        echo "✂️ Running Trimmomatic paired-end trimming..."

        # R1: use the primary input file directly (preserves .gz extension for trimmomatic)
        R1_FILE="$INPUT_FILE_PATH"

        # R2: locate using INPUT_S3_KEY_2 filename, fall back to second file in /tmp/input
        R2_FILE=""
        if [[ -n "$INPUT_S3_KEY_2" ]]; then
            R2_FILENAME=$(basename "$INPUT_S3_KEY_2")
            R2_FILE=$(find /tmp/input -name "$R2_FILENAME" | head -n 1)
        fi
        if [[ -z "$R2_FILE" ]]; then
            R2_FILE=$(find /tmp/input -type f ! -name "$(basename "$R1_FILE")" | head -n 1)
        fi
        if [[ -z "$R2_FILE" ]]; then
            echo "❌ No R2 input file found for paired-end processing"
            exit 1
        fi

        echo "📁 R1: $(basename "$R1_FILE")"
        echo "📁 R2: $(basename "$R2_FILE")"

        # Output file names — use OUTPUT_FILE / OUTPUT_FILE_2 env vars from ECS
        OUTPUT_R1="${OUTPUT_FILE:-"trimmed_reads_R1.fastq.gz"}"
        OUTPUT_R2="${OUTPUT_FILE_2:-"trimmed_reads_R2.fastq.gz"}"
        OUTPUT_UNPAIRED_R1="unpaired_${OUTPUT_R1}"
        OUTPUT_UNPAIRED_R2="unpaired_${OUTPUT_R2}"

        TRIMMOMATIC_CMD="trimmomatic PE -phred33"

        # Add threads parameter
        if [[ -n "$THREADS" ]]; then
            TRIMMOMATIC_CMD="$TRIMMOMATIC_CMD -threads $THREADS"
        fi

        # Input and output files (trimmomatic handles .fastq.gz natively)
        TRIMMOMATIC_CMD="$TRIMMOMATIC_CMD $R1_FILE $R2_FILE"
        TRIMMOMATIC_CMD="$TRIMMOMATIC_CMD /tmp/output/$OUTPUT_R1 /tmp/output/$OUTPUT_UNPAIRED_R1"
        TRIMMOMATIC_CMD="$TRIMMOMATIC_CMD /tmp/output/$OUTPUT_R2 /tmp/output/$OUTPUT_UNPAIRED_R2"

        # Add adapter trimming if specified
        if [[ -n "$ADAPTER_FILE" ]]; then
            TRIMMOMATIC_CMD="$TRIMMOMATIC_CMD ILLUMINACLIP:$ADAPTER_FILE:2:30:10"
        fi

        # Add quality trimming parameters
        if [[ -n "$LEADING" ]]; then
            TRIMMOMATIC_CMD="$TRIMMOMATIC_CMD LEADING:$LEADING"
        fi

        if [[ -n "$TRAILING" ]]; then
            TRIMMOMATIC_CMD="$TRIMMOMATIC_CMD TRAILING:$TRAILING"
        fi

        if [[ -n "$SLIDING_WINDOW" ]]; then
            TRIMMOMATIC_CMD="$TRIMMOMATIC_CMD SLIDINGWINDOW:$SLIDING_WINDOW"
        fi

        if [[ -n "$MIN_LENGTH" && "$MIN_LENGTH" != "0" ]]; then
            TRIMMOMATIC_CMD="$TRIMMOMATIC_CMD MINLEN:$MIN_LENGTH"
        fi

        echo "🚀 Executing: $TRIMMOMATIC_CMD"

        # Execute the command
        if eval "$TRIMMOMATIC_CMD"; then
            echo "✅ Trimmomatic paired-end trimming completed successfully"

            # Display output file info
            echo "📁 Generated files:"
            ls -la /tmp/output/
        else
            echo "❌ Trimmomatic paired-end trimming failed"
            exit 1
        fi
        ;;
        
    "quality-trim")
        echo "✂️ Running Trimmomatic quality-only trimming..."

        # Use the input file directly (trimmomatic handles .fastq.gz natively)
        OUTPUT_FILE="${OUTPUT_FILE:-"trimmed_reads.fastq.gz"}"

        TRIMMOMATIC_CMD="trimmomatic SE -phred33"

        # Add threads parameter
        if [[ -n "$THREADS" ]]; then
            TRIMMOMATIC_CMD="$TRIMMOMATIC_CMD -threads $THREADS"
        fi

        # Add input and output files
        TRIMMOMATIC_CMD="$TRIMMOMATIC_CMD $INPUT_FILE_PATH /tmp/output/$OUTPUT_FILE"
        
        # Quality-only trimming (no adapter removal)
        if [[ -n "$QUALITY_THRESHOLD" ]]; then
            TRIMMOMATIC_CMD="$TRIMMOMATIC_CMD AVGQUAL:$QUALITY_THRESHOLD"
        fi
        
        if [[ -n "$LEADING" ]]; then
            TRIMMOMATIC_CMD="$TRIMMOMATIC_CMD LEADING:$LEADING"
        fi
        
        if [[ -n "$TRAILING" ]]; then
            TRIMMOMATIC_CMD="$TRIMMOMATIC_CMD TRAILING:$TRAILING"
        fi
        
        if [[ -n "$MIN_LENGTH" ]]; then
            TRIMMOMATIC_CMD="$TRIMMOMATIC_CMD MINLEN:$MIN_LENGTH"
        fi
        
        echo "🚀 Executing: $TRIMMOMATIC_CMD"
        
        # Execute the command
        if eval "$TRIMMOMATIC_CMD"; then
            echo "✅ Trimmomatic quality trimming completed successfully"
            
            # Display output file info
            OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILE" 2>/dev/null || echo "unknown")
            echo "📊 Output file: $OUTPUT_FILE ($OUTPUT_SIZE bytes)"
        else
            echo "❌ Trimmomatic quality trimming failed"
            exit 1
        fi
        ;;
        
    *)
        echo "❌ Unsupported Trimmomatic command: $COMMAND"
        echo "Supported commands: single-end, paired-end, quality-trim"
        exit 1
        ;;
esac

echo "🎯 Trimmomatic handler completed successfully"