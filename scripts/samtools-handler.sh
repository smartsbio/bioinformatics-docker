#!/bin/bash
set -e

# SAMtools Handler Script
# Handles various SAMtools commands and operations

echo "🧬 SAMtools Handler Started"

# Get the input file (should be downloaded already)
INPUT_FILES=(/tmp/input/*)
if [[ ${#INPUT_FILES[@]} -eq 0 ]]; then
    echo "❌ No input files found in /tmp/input/"
    exit 1
fi

INPUT_FILE_PATH="${INPUT_FILES[0]}"
INPUT_FILENAME=$(basename "$INPUT_FILE_PATH")

echo "📁 Processing file: $INPUT_FILENAME"
echo "🎯 SAMtools command: $COMMAND"

# Parse additional parameters from environment
OUTPUT_FORMAT=${OUTPUT_FORMAT:-"bam"}
OUTPUT_FILE=${OUTPUT_FILE:-"output.${OUTPUT_FORMAT}"}
THREADS=${THREADS:-"4"}
COMPRESSION_LEVEL=${COMPRESSION_LEVEL:-"6"}
REGION=${REGION:-""}
FLAGS=${FLAGS:-""}
REFERENCE_GENOME=${REFERENCE_GENOME:-""}

case "$COMMAND" in
    "view")
        echo "🔍 Running SAMtools view command..."
        
        # Build SAMtools view command based on output format
        SAMTOOLS_CMD="samtools view"
        
        # Add threads parameter
        if [[ -n "$THREADS" ]]; then
            SAMTOOLS_CMD="$SAMTOOLS_CMD --threads $THREADS"
        fi
        
        # Add compression level for BAM/CRAM
        case "$OUTPUT_FORMAT" in
            "bam")
                SAMTOOLS_CMD="$SAMTOOLS_CMD -b"
                if [[ -n "$COMPRESSION_LEVEL" ]]; then
                    SAMTOOLS_CMD="$SAMTOOLS_CMD -l $COMPRESSION_LEVEL"
                fi
                ;;
            "sam")
                SAMTOOLS_CMD="$SAMTOOLS_CMD -h"
                ;;
            "cram")
                SAMTOOLS_CMD="$SAMTOOLS_CMD -C"
                if [[ -n "$REFERENCE_GENOME" ]]; then
                    SAMTOOLS_CMD="$SAMTOOLS_CMD -T $REFERENCE_GENOME"
                fi
                ;;
            *)
                echo "⚠️ Unknown output format: $OUTPUT_FORMAT, defaulting to SAM"
                SAMTOOLS_CMD="$SAMTOOLS_CMD -h"
                ;;
        esac
        
        # Add region if specified
        if [[ -n "$REGION" ]]; then
            SAMTOOLS_CMD="$SAMTOOLS_CMD $INPUT_FILE_PATH $REGION"
        else
            SAMTOOLS_CMD="$SAMTOOLS_CMD $INPUT_FILE_PATH"
        fi
        
        # Add any additional flags
        if [[ -n "$FLAGS" ]]; then
            SAMTOOLS_CMD="$SAMTOOLS_CMD $FLAGS"
        fi
        
        echo "🚀 Executing: $SAMTOOLS_CMD > /tmp/output/$OUTPUT_FILE"
        
        # Execute the command
        if eval "$SAMTOOLS_CMD > /tmp/output/$OUTPUT_FILE"; then
            echo "✅ SAMtools view completed successfully"
            
            # Display output file info
            OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILE" 2>/dev/null || echo "unknown")
            echo "📊 Output file: $OUTPUT_FILE ($OUTPUT_SIZE bytes)"
        else
            echo "❌ SAMtools view failed"
            exit 1
        fi
        ;;
        
    "sort")
        echo "🔄 Running SAMtools sort command..."
        
        SAMTOOLS_CMD="samtools sort"
        
        # Add threads parameter
        if [[ -n "$THREADS" ]]; then
            SAMTOOLS_CMD="$SAMTOOLS_CMD --threads $THREADS"
        fi
        
        # Add compression level
        if [[ -n "$COMPRESSION_LEVEL" ]]; then
            SAMTOOLS_CMD="$SAMTOOLS_CMD -l $COMPRESSION_LEVEL"
        fi
        
        SAMTOOLS_CMD="$SAMTOOLS_CMD $INPUT_FILE_PATH -o /tmp/output/$OUTPUT_FILE"
        
        echo "🚀 Executing: $SAMTOOLS_CMD"
        
        if eval "$SAMTOOLS_CMD"; then
            echo "✅ SAMtools sort completed successfully"
        else
            echo "❌ SAMtools sort failed"
            exit 1
        fi
        ;;
        
    "index")
        echo "📇 Running SAMtools index command..."
        
        # For indexing, we need to work with the file in place
        cp "$INPUT_FILE_PATH" "/tmp/output/$INPUT_FILENAME"
        
        SAMTOOLS_CMD="samtools index /tmp/output/$INPUT_FILENAME"
        
        echo "🚀 Executing: $SAMTOOLS_CMD"
        
        if eval "$SAMTOOLS_CMD"; then
            echo "✅ SAMtools index completed successfully"
            
            # List all generated files
            echo "📁 Generated files:"
            ls -la /tmp/output/
        else
            echo "❌ SAMtools index failed"
            exit 1
        fi
        ;;
        
    "stats")
        echo "📈 Running SAMtools stats command..."
        
        SAMTOOLS_CMD="samtools stats $INPUT_FILE_PATH"
        
        echo "🚀 Executing: $SAMTOOLS_CMD > /tmp/output/stats.txt"
        
        if eval "$SAMTOOLS_CMD > /tmp/output/stats.txt"; then
            echo "✅ SAMtools stats completed successfully"
        else
            echo "❌ SAMtools stats failed"
            exit 1
        fi
        ;;
        
    *)
        echo "❌ Unsupported SAMtools command: $COMMAND"
        echo "Supported commands: view, sort, index, stats"
        exit 1
        ;;
esac

echo "🎯 SAMtools handler completed successfully"