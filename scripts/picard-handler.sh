#!/bin/bash
set -e

# Picard Handler Script
# Handles various Picard tools for SAM/BAM/VCF file manipulation

echo "🧬 Picard Handler Started"

# Get the input file (should be downloaded already)
INPUT_FILES=(/tmp/input/*)
if [[ ${#INPUT_FILES[@]} -eq 0 ]]; then
    echo "❌ No input files found in /tmp/input/"
    exit 1
fi

INPUT_FILE_PATH="${INPUT_FILES[0]}"
INPUT_FILENAME=$(basename "$INPUT_FILE_PATH")

echo "📁 Processing file: $INPUT_FILENAME"
echo "🎯 Picard command: $COMMAND"

# Parse additional parameters from environment
OUTPUT_FILE=${OUTPUT_FILE:-"output.bam"}
SORT_ORDER=${SORT_ORDER:-"coordinate"}
VALIDATION_STRINGENCY=${VALIDATION_STRINGENCY:-"LENIENT"}
MAX_RECORDS_IN_RAM=${MAX_RECORDS_IN_RAM:-"500000"}
COMPRESSION_LEVEL=${COMPRESSION_LEVEL:-"5"}
CREATE_INDEX=${CREATE_INDEX:-"true"}

case "$COMMAND" in
    "sort")
        echo "🔄 Running Picard SortSam..."
        
        # Copy input file
        cp "$INPUT_FILE_PATH" "/tmp/input.bam"
        
        PICARD_CMD="picard SortSam"
        PICARD_CMD="$PICARD_CMD INPUT=/tmp/input.bam"
        PICARD_CMD="$PICARD_CMD OUTPUT=/tmp/output/$OUTPUT_FILE"
        
        # Add sort order
        if [[ -n "$SORT_ORDER" ]]; then
            PICARD_CMD="$PICARD_CMD SORT_ORDER=$SORT_ORDER"
        fi
        
        # Add validation stringency
        if [[ -n "$VALIDATION_STRINGENCY" ]]; then
            PICARD_CMD="$PICARD_CMD VALIDATION_STRINGENCY=$VALIDATION_STRINGENCY"
        fi
        
        # Add memory parameters
        if [[ -n "$MAX_RECORDS_IN_RAM" ]]; then
            PICARD_CMD="$PICARD_CMD MAX_RECORDS_IN_RAM=$MAX_RECORDS_IN_RAM"
        fi
        
        # Add compression
        if [[ -n "$COMPRESSION_LEVEL" ]]; then
            PICARD_CMD="$PICARD_CMD COMPRESSION_LEVEL=$COMPRESSION_LEVEL"
        fi
        
        # Create index
        if [[ "$CREATE_INDEX" == "true" ]]; then
            PICARD_CMD="$PICARD_CMD CREATE_INDEX=true"
        fi
        
        echo "🚀 Executing: $PICARD_CMD"
        
        # Execute the command
        if eval "$PICARD_CMD"; then
            echo "✅ Picard SortSam completed successfully"
            
            # Display output file info
            OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILE" 2>/dev/null || echo "unknown")
            echo "📊 Output file: $OUTPUT_FILE ($OUTPUT_SIZE bytes)"
            
            # List all generated files
            echo "📁 Generated files:"
            ls -la /tmp/output/
        else
            echo "❌ Picard SortSam failed"
            exit 1
        fi
        ;;
        
    "mark-duplicates")
        echo "🔍 Running Picard MarkDuplicates..."
        
        # Copy input file
        cp "$INPUT_FILE_PATH" "/tmp/input.bam"
        
        METRICS_FILE=${METRICS_FILE:-"duplicate_metrics.txt"}
        
        PICARD_CMD="picard MarkDuplicates"
        PICARD_CMD="$PICARD_CMD INPUT=/tmp/input.bam"
        PICARD_CMD="$PICARD_CMD OUTPUT=/tmp/output/$OUTPUT_FILE"
        PICARD_CMD="$PICARD_CMD METRICS_FILE=/tmp/output/$METRICS_FILE"
        
        # Add validation stringency
        if [[ -n "$VALIDATION_STRINGENCY" ]]; then
            PICARD_CMD="$PICARD_CMD VALIDATION_STRINGENCY=$VALIDATION_STRINGENCY"
        fi
        
        # Create index
        if [[ "$CREATE_INDEX" == "true" ]]; then
            PICARD_CMD="$PICARD_CMD CREATE_INDEX=true"
        fi
        
        # Remove duplicates option
        if [[ "$REMOVE_DUPLICATES" == "true" ]]; then
            PICARD_CMD="$PICARD_CMD REMOVE_DUPLICATES=true"
        fi
        
        echo "🚀 Executing: $PICARD_CMD"
        
        # Execute the command
        if eval "$PICARD_CMD"; then
            echo "✅ Picard MarkDuplicates completed successfully"
            
            # Display output files info
            echo "📁 Generated files:"
            ls -la /tmp/output/
        else
            echo "❌ Picard MarkDuplicates failed"
            exit 1
        fi
        ;;
        
    "merge")
        echo "🔗 Running Picard MergeSamFiles..."
        
        # Copy input files (for merging multiple files)
        cp "$INPUT_FILE_PATH" "/tmp/input1.bam"
        # Note: In real implementation, we'd handle multiple input files
        
        PICARD_CMD="picard MergeSamFiles"
        PICARD_CMD="$PICARD_CMD INPUT=/tmp/input1.bam"
        # Add more inputs as needed: INPUT=/tmp/input2.bam INPUT=/tmp/input3.bam
        PICARD_CMD="$PICARD_CMD OUTPUT=/tmp/output/$OUTPUT_FILE"
        
        # Add validation stringency
        if [[ -n "$VALIDATION_STRINGENCY" ]]; then
            PICARD_CMD="$PICARD_CMD VALIDATION_STRINGENCY=$VALIDATION_STRINGENCY"
        fi
        
        # Create index
        if [[ "$CREATE_INDEX" == "true" ]]; then
            PICARD_CMD="$PICARD_CMD CREATE_INDEX=true"
        fi
        
        # Merge sequence dictionaries
        PICARD_CMD="$PICARD_CMD MERGE_SEQUENCE_DICTIONARIES=true"
        
        echo "🚀 Executing: $PICARD_CMD"
        
        # Execute the command
        if eval "$PICARD_CMD"; then
            echo "✅ Picard MergeSamFiles completed successfully"
            
            # Display output file info
            OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILE" 2>/dev/null || echo "unknown")
            echo "📊 Output file: $OUTPUT_FILE ($OUTPUT_SIZE bytes)"
        else
            echo "❌ Picard MergeSamFiles failed"
            exit 1
        fi
        ;;
        
    "validate")
        echo "✅ Running Picard ValidateSamFile..."
        
        # Copy input file
        cp "$INPUT_FILE_PATH" "/tmp/input.bam"
        
        VALIDATION_REPORT=${OUTPUT_FILE:-"validation_report.txt"}
        
        PICARD_CMD="picard ValidateSamFile"
        PICARD_CMD="$PICARD_CMD INPUT=/tmp/input.bam"
        PICARD_CMD="$PICARD_CMD OUTPUT=/tmp/output/$VALIDATION_REPORT"
        
        # Add validation mode
        PICARD_CMD="$PICARD_CMD MODE=VERBOSE"
        
        # Add validation stringency
        if [[ -n "$VALIDATION_STRINGENCY" ]]; then
            PICARD_CMD="$PICARD_CMD VALIDATION_STRINGENCY=$VALIDATION_STRINGENCY"
        fi
        
        echo "🚀 Executing: $PICARD_CMD"
        
        # Execute the command
        if eval "$PICARD_CMD"; then
            echo "✅ Picard ValidateSamFile completed successfully"
            
            # Display output file info
            OUTPUT_SIZE=$(stat -c%s "/tmp/output/$VALIDATION_REPORT" 2>/dev/null || echo "unknown")
            echo "📊 Validation report: $VALIDATION_REPORT ($OUTPUT_SIZE bytes)"
        else
            echo "❌ Picard ValidateSamFile failed"
            exit 1
        fi
        ;;
        
    "collect-metrics")
        echo "📊 Running Picard CollectMultipleMetrics..."
        
        # Copy input file
        cp "$INPUT_FILE_PATH" "/tmp/input.bam"
        
        BASE_NAME=${BASE_NAME:-"metrics"}
        
        PICARD_CMD="picard CollectMultipleMetrics"
        PICARD_CMD="$PICARD_CMD INPUT=/tmp/input.bam"
        PICARD_CMD="$PICARD_CMD OUTPUT=/tmp/output/$BASE_NAME"
        
        # Add validation stringency
        if [[ -n "$VALIDATION_STRINGENCY" ]]; then
            PICARD_CMD="$PICARD_CMD VALIDATION_STRINGENCY=$VALIDATION_STRINGENCY"
        fi
        
        # Add reference if provided
        if [[ -n "$REFERENCE_SEQUENCE" ]]; then
            PICARD_CMD="$PICARD_CMD REFERENCE_SEQUENCE=$REFERENCE_SEQUENCE"
        fi
        
        echo "🚀 Executing: $PICARD_CMD"
        
        # Execute the command
        if eval "$PICARD_CMD"; then
            echo "✅ Picard CollectMultipleMetrics completed successfully"
            
            # Display all generated metric files
            echo "📁 Generated metric files:"
            ls -la /tmp/output/
        else
            echo "❌ Picard CollectMultipleMetrics failed"
            exit 1
        fi
        ;;
        
    *)
        echo "❌ Unsupported Picard command: $COMMAND"
        echo "Supported commands: sort, mark-duplicates, merge, validate, collect-metrics"
        exit 1
        ;;
esac

echo "🎯 Picard handler completed successfully"