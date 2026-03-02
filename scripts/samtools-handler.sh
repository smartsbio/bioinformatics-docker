#!/bin/bash
set -e

# SAMtools Handler Script
# Handles various SAMtools commands and operations

echo "🧬 SAMtools Handler Started"

# Get the input file (should be downloaded already)
# Use find to locate files in subdirectories (handles nested folder structures)
INPUT_FILE_PATH=$(find /tmp/input -type f | head -n 1)

if [[ -z "$INPUT_FILE_PATH" ]]; then
    echo "❌ No input files found in /tmp/input/"
    exit 1
fi

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
        # Convert OUTPUT_FORMAT to lowercase for case-insensitive matching
        OUTPUT_FORMAT_LOWER=$(echo "$OUTPUT_FORMAT" | tr '[:upper:]' '[:lower:]')
        case "$OUTPUT_FORMAT_LOWER" in
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

        # Create subdirectories in output path if needed
        OUTPUT_DIR=$(dirname "/tmp/output/$OUTPUT_FILE")
        if [[ "$OUTPUT_DIR" != "/tmp/output" ]]; then
            mkdir -p "$OUTPUT_DIR"
            echo "📁 Created output directory: $OUTPUT_DIR"
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

        # Limit memory per thread to prevent OOM on constrained Fargate instances.
        # Default: 768M/thread (samtools default) × 9 threads = 6+ GB — far exceeds 2 GB Fargate limit.
        # 128M × (8+1) threads = 1152 MB — safe within 2 GB.
        SORT_MEMORY=${SORT_MEMORY:-"128M"}

        SAMTOOLS_CMD="samtools sort"

        # Add threads parameter
        if [[ -n "$THREADS" ]]; then
            SAMTOOLS_CMD="$SAMTOOLS_CMD --threads $THREADS"
        fi

        # Cap memory per thread
        SAMTOOLS_CMD="$SAMTOOLS_CMD -m $SORT_MEMORY"

        # Add compression level
        if [[ -n "$COMPRESSION_LEVEL" ]]; then
            SAMTOOLS_CMD="$SAMTOOLS_CMD -l $COMPRESSION_LEVEL"
        fi

        # Create subdirectories in output path if needed
        OUTPUT_DIR=$(dirname "/tmp/output/$OUTPUT_FILE")
        if [[ "$OUTPUT_DIR" != "/tmp/output" ]]; then
            mkdir -p "$OUTPUT_DIR"
            echo "📁 Created output directory: $OUTPUT_DIR"
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

        # Create subdirectories in output path if needed
        OUTPUT_DIR=$(dirname "/tmp/output/$INPUT_FILENAME")
        if [[ "$OUTPUT_DIR" != "/tmp/output" ]]; then
            mkdir -p "$OUTPUT_DIR"
            echo "📁 Created output directory: $OUTPUT_DIR"
        fi

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
        
    "flagstat")
        echo "📊 Running SAMtools flagstat command..."
        
        SAMTOOLS_CMD="samtools flagstat $INPUT_FILE_PATH"
        
        echo "🚀 Executing: $SAMTOOLS_CMD > /tmp/output/flagstat.txt"
        
        if eval "$SAMTOOLS_CMD > /tmp/output/flagstat.txt"; then
            echo "✅ SAMtools flagstat completed successfully"
            
            # Show summary
            echo "📈 Alignment statistics:"
            cat "/tmp/output/flagstat.txt"
        else
            echo "❌ SAMtools flagstat failed"
            exit 1
        fi
        ;;
        
    "depth")
        echo "📏 Running SAMtools depth command..."
        
        SAMTOOLS_CMD="samtools depth"
        
        # Add region if specified
        if [[ -n "$REGION" ]]; then
            SAMTOOLS_CMD="$SAMTOOLS_CMD -r $REGION"
        fi
        
        SAMTOOLS_CMD="$SAMTOOLS_CMD $INPUT_FILE_PATH"
        
        echo "🚀 Executing: $SAMTOOLS_CMD > /tmp/output/depth.txt"
        
        if eval "$SAMTOOLS_CMD > /tmp/output/depth.txt"; then
            echo "✅ SAMtools depth completed successfully"
            
            # Show sample of depth data
            echo "📊 Depth data sample:"
            head -20 "/tmp/output/depth.txt" | tail -10
        else
            echo "❌ SAMtools depth failed"
            exit 1
        fi
        ;;
        
    "merge")
        echo "🔗 Running SAMtools merge command..."
        
        # For merge, we'd typically need multiple input files
        # This is a simplified version
        SAMTOOLS_CMD="samtools merge"
        
        # Add threads parameter
        if [[ -n "$THREADS" ]]; then
            SAMTOOLS_CMD="$SAMTOOLS_CMD --threads $THREADS"
        fi
        
        # Add compression level
        if [[ -n "$COMPRESSION_LEVEL" ]]; then
            SAMTOOLS_CMD="$SAMTOOLS_CMD -l $COMPRESSION_LEVEL"
        fi
        
        SAMTOOLS_CMD="$SAMTOOLS_CMD /tmp/output/$OUTPUT_FILE $INPUT_FILE_PATH"
        
        echo "🚀 Executing: $SAMTOOLS_CMD"
        
        if eval "$SAMTOOLS_CMD"; then
            echo "✅ SAMtools merge completed successfully"
        else
            echo "❌ SAMtools merge failed"
            exit 1
        fi
        ;;
        
    "idxstats")
        echo "📊 Running SAMtools idxstats command..."
        
        SAMTOOLS_CMD="samtools idxstats $INPUT_FILE_PATH"
        
        echo "🚀 Executing: $SAMTOOLS_CMD > /tmp/output/idxstats.txt"
        
        if eval "$SAMTOOLS_CMD > /tmp/output/idxstats.txt"; then
            echo "✅ SAMtools idxstats completed successfully"
            
            # Show index statistics
            echo "📈 Index statistics:"
            cat "/tmp/output/idxstats.txt"
        else
            echo "❌ SAMtools idxstats failed"
            exit 1
        fi
        ;;
        
    "coverage")
        echo "📊 Running SAMtools coverage command..."
        
        SAMTOOLS_CMD="samtools coverage"
        
        # Add region if specified
        if [[ -n "$REGION" ]]; then
            SAMTOOLS_CMD="$SAMTOOLS_CMD -r $REGION"
        fi
        
        SAMTOOLS_CMD="$SAMTOOLS_CMD $INPUT_FILE_PATH"
        
        echo "🚀 Executing: $SAMTOOLS_CMD > /tmp/output/coverage.txt"
        
        if eval "$SAMTOOLS_CMD > /tmp/output/coverage.txt"; then
            echo "✅ SAMtools coverage completed successfully"
            
            # Show coverage summary
            echo "📊 Coverage summary:"
            head -10 "/tmp/output/coverage.txt"
        else
            echo "❌ SAMtools coverage failed"
            exit 1
        fi
        ;;
        
    "tview")
        echo "👁️ Running SAMtools tview command (text viewer)..."
        
        SAMTOOLS_CMD="samtools tview"
        
        # Add reference genome if available
        if [[ -n "$REFERENCE_GENOME" ]]; then
            SAMTOOLS_CMD="$SAMTOOLS_CMD $INPUT_FILE_PATH $REFERENCE_GENOME"
        else
            SAMTOOLS_CMD="$SAMTOOLS_CMD $INPUT_FILE_PATH"
        fi
        
        # For non-interactive use, capture a region
        if [[ -n "$REGION" ]]; then
            SAMTOOLS_CMD="$SAMTOOLS_CMD -p $REGION"
        fi
        
        echo "🚀 Executing: $SAMTOOLS_CMD > /tmp/output/tview.txt"
        
        # Note: tview is interactive, this captures limited output
        if timeout 5 eval "$SAMTOOLS_CMD > /tmp/output/tview.txt 2>&1" || [[ $? -eq 124 ]]; then
            echo "✅ SAMtools tview completed (non-interactive mode)"
        else
            echo "❌ SAMtools tview failed"
            exit 1
        fi
        ;;
        
    *)
        echo "❌ Unsupported SAMtools command: $COMMAND"
        echo "Supported commands: view, sort, index, stats, flagstat, depth, merge, idxstats, coverage, tview"
        exit 1
        ;;
esac

echo "🎯 SAMtools handler completed successfully"