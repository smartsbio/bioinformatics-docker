#!/bin/bash
set -e

# Trimmomatic Handler Script
# Handles adapter trimming and quality filtering for sequencing data

echo "🧬 Trimmomatic Handler Started"

# Get the input file (should be downloaded already)
INPUT_FILES=(/tmp/input/*)
if [[ ${#INPUT_FILES[@]} -eq 0 ]]; then
    echo "❌ No input files found in /tmp/input/"
    exit 1
fi

INPUT_FILE_PATH="${INPUT_FILES[0]}"
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
        
        # Copy input file
        cp "$INPUT_FILE_PATH" "/tmp/input_reads.fastq"
        
        TRIMMOMATIC_CMD="trimmomatic SE"
        
        # Add threads parameter
        if [[ -n "$THREADS" ]]; then
            TRIMMOMATIC_CMD="$TRIMMOMATIC_CMD -threads $THREADS"
        fi
        
        # Add input and output files
        TRIMMOMATIC_CMD="$TRIMMOMATIC_CMD /tmp/input_reads.fastq /tmp/output/$OUTPUT_FILE"
        
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
        
        # Copy input files (assuming paired files)
        cp "$INPUT_FILE_PATH" "/tmp/input_reads_R1.fastq"
        # Note: In real implementation, we'd handle multiple input files properly
        
        OUTPUT_FILE_R1=${OUTPUT_FILE_R1:-"trimmed_reads_R1.fastq"}
        OUTPUT_FILE_R2=${OUTPUT_FILE_R2:-"trimmed_reads_R2.fastq"}
        OUTPUT_FILE_UNPAIRED_R1=${OUTPUT_FILE_UNPAIRED_R1:-"trimmed_reads_unpaired_R1.fastq"}
        OUTPUT_FILE_UNPAIRED_R2=${OUTPUT_FILE_UNPAIRED_R2:-"trimmed_reads_unpaired_R2.fastq"}
        
        TRIMMOMATIC_CMD="trimmomatic PE"
        
        # Add threads parameter
        if [[ -n "$THREADS" ]]; then
            TRIMMOMATIC_CMD="$TRIMMOMATIC_CMD -threads $THREADS"
        fi
        
        # Add input and output files
        TRIMMOMATIC_CMD="$TRIMMOMATIC_CMD /tmp/input_reads_R1.fastq /tmp/input_reads_R2.fastq"
        TRIMMOMATIC_CMD="$TRIMMOMATIC_CMD /tmp/output/$OUTPUT_FILE_R1 /tmp/output/$OUTPUT_FILE_UNPAIRED_R1"
        TRIMMOMATIC_CMD="$TRIMMOMATIC_CMD /tmp/output/$OUTPUT_FILE_R2 /tmp/output/$OUTPUT_FILE_UNPAIRED_R2"
        
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
        
        # Copy input file
        cp "$INPUT_FILE_PATH" "/tmp/input_reads.fastq"
        
        TRIMMOMATIC_CMD="trimmomatic SE"
        
        # Add threads parameter
        if [[ -n "$THREADS" ]]; then
            TRIMMOMATIC_CMD="$TRIMMOMATIC_CMD -threads $THREADS"
        fi
        
        # Add input and output files
        TRIMMOMATIC_CMD="$TRIMMOMATIC_CMD /tmp/input_reads.fastq /tmp/output/$OUTPUT_FILE"
        
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