#!/bin/bash
set -e

# FastQC Handler Script
# Handles quality control analysis for sequencing data

echo "🧬 FastQC Handler Started"

# Get the input file (should be downloaded already)
# Use find to locate files in subdirectories (handles nested folder structures)
INPUT_FILE_PATH=$(find /tmp/input -type f | head -n 1)

if [[ -z "$INPUT_FILE_PATH" ]]; then
    echo "❌ No input files found in /tmp/input/"
    exit 1
fi

INPUT_FILENAME=$(basename "$INPUT_FILE_PATH")

echo "📁 Processing file: $INPUT_FILENAME"
echo "🎯 FastQC command: $COMMAND"

# Parse additional parameters from environment
OUTPUT_DIR=${OUTPUT_DIR:-"/tmp/output"}
THREADS=${THREADS:-"1"}
KMER_SIZE=${KMER_SIZE:-"7"}
FORMAT=${FORMAT:-"auto"}
EXTRACT=${EXTRACT:-"true"}

case "$COMMAND" in
    "quality-report")
        echo "📊 Running FastQC quality analysis..."

        FASTQC_CMD="fastqc"

        # Add threads parameter
        if [[ -n "$THREADS" ]]; then
            FASTQC_CMD="$FASTQC_CMD --threads $THREADS"
        fi

        # Add output directory
        FASTQC_CMD="$FASTQC_CMD --outdir $OUTPUT_DIR"

        # Add format if specified
        if [[ "$FORMAT" != "auto" ]]; then
            FASTQC_CMD="$FASTQC_CMD --format $FORMAT"
        fi

        # Add kmer size
        if [[ "$KMER_SIZE" != "7" ]]; then
            FASTQC_CMD="$FASTQC_CMD --kmers $KMER_SIZE"
        fi

        # Extract results
        if [[ "$EXTRACT" == "true" ]]; then
            FASTQC_CMD="$FASTQC_CMD --extract"
        fi

        # Add input file - use actual input file path
        FASTQC_CMD="$FASTQC_CMD \"$INPUT_FILE_PATH\""
        
        echo "🚀 Executing: $FASTQC_CMD"
        
        # Execute the command
        if eval "$FASTQC_CMD"; then
            echo "✅ FastQC quality analysis completed successfully"
            
            # List output files
            echo "📁 Generated files:"
            ls -la $OUTPUT_DIR/
            
            # Show file sizes
            for file in $OUTPUT_DIR/*; do
                if [[ -f "$file" ]]; then
                    filename=$(basename "$file")
                    filesize=$(stat -c%s "$file" 2>/dev/null || echo "unknown")
                    echo "📊 Output file: $filename ($filesize bytes)"
                fi
            done
        else
            echo "❌ FastQC quality analysis failed"
            exit 1
        fi
        ;;
        
    "batch-report")
        echo "📊 Running FastQC batch quality analysis..."
        
        # For batch processing, process all input files
        FASTQC_CMD="fastqc"
        
        # Add threads parameter
        if [[ -n "$THREADS" ]]; then
            FASTQC_CMD="$FASTQC_CMD --threads $THREADS"
        fi
        
        # Add output directory
        FASTQC_CMD="$FASTQC_CMD --outdir $OUTPUT_DIR"
        
        # Add format if specified
        if [[ "$FORMAT" != "auto" ]]; then
            FASTQC_CMD="$FASTQC_CMD --format $FORMAT"
        fi
        
        # Extract results
        if [[ "$EXTRACT" == "true" ]]; then
            FASTQC_CMD="$FASTQC_CMD --extract"
        fi
        
        # Add all input files
        for input_file in "${INPUT_FILES[@]}"; do
            FASTQC_CMD="$FASTQC_CMD $input_file"
        done
        
        echo "🚀 Executing: $FASTQC_CMD"
        
        # Execute the command
        if eval "$FASTQC_CMD"; then
            echo "✅ FastQC batch analysis completed successfully"
            
            # List output files
            echo "📁 Generated files:"
            ls -la $OUTPUT_DIR/
        else
            echo "❌ FastQC batch analysis failed"
            exit 1
        fi
        ;;
        
    "custom-limits")
        echo "📊 Running FastQC with custom limits..."

        FASTQC_CMD="fastqc"

        # Add threads parameter
        if [[ -n "$THREADS" ]]; then
            FASTQC_CMD="$FASTQC_CMD --threads $THREADS"
        fi

        # Add output directory
        FASTQC_CMD="$FASTQC_CMD --outdir $OUTPUT_DIR"

        # Add custom limits file if specified
        if [[ -n "$LIMITS_FILE" ]]; then
            FASTQC_CMD="$FASTQC_CMD --limits $LIMITS_FILE"
        fi

        # Add contaminants file if specified
        if [[ -n "$CONTAMINANTS_FILE" ]]; then
            FASTQC_CMD="$FASTQC_CMD --contaminants $CONTAMINANTS_FILE"
        fi

        # Add adapters file if specified
        if [[ -n "$ADAPTERS_FILE" ]]; then
            FASTQC_CMD="$FASTQC_CMD --adapters $ADAPTERS_FILE"
        fi

        # Add input file - use actual input file path
        FASTQC_CMD="$FASTQC_CMD \"$INPUT_FILE_PATH\""
        
        echo "🚀 Executing: $FASTQC_CMD"
        
        # Execute the command
        if eval "$FASTQC_CMD"; then
            echo "✅ FastQC custom analysis completed successfully"
            
            # List output files
            echo "📁 Generated files:"
            ls -la $OUTPUT_DIR/
        else
            echo "❌ FastQC custom analysis failed"
            exit 1
        fi
        ;;
        
    *)
        echo "❌ Unsupported FastQC command: $COMMAND"
        echo "Supported commands: quality-report, batch-report, custom-limits"
        exit 1
        ;;
esac

echo "🎯 FastQC handler completed successfully"