#!/bin/bash
set -e

# Picard Handler Script
# Handles various Picard tools for SAM/BAM/VCF file manipulation

echo "🧬 Picard Handler Started"

# Get the input file (should be downloaded already)
# Use find to locate files in subdirectories (handles nested folder structures)
INPUT_FILE_PATH=$(find /tmp/input -type f | head -n 1)

if [[ -z "$INPUT_FILE_PATH" ]]; then
    echo "❌ No input files found in /tmp/input/"
    exit 1
fi

INPUT_FILENAME=$(basename "$INPUT_FILE_PATH")

echo "📁 Processing file: $INPUT_FILENAME"
echo "🎯 Picard command: $COMMAND"

# Extract organization and workspace IDs from INPUT_S3_KEY
# Format: organizations/{orgId}/workspaces/{workspaceId}/files/{path}
if [[ -n "$INPUT_S3_KEY" ]]; then
    ORGANIZATION_ID=$(echo "$INPUT_S3_KEY" | cut -d'/' -f2)
    WORKSPACE_ID=$(echo "$INPUT_S3_KEY" | cut -d'/' -f4)
    echo "📂 Organization ID: $ORGANIZATION_ID"
    echo "📂 Workspace ID: $WORKSPACE_ID"
fi

# Parse additional parameters from environment
OUTPUT_FILE=${OUTPUT_FILE:-"output.bam"}
# Strip @ notation from output file path
OUTPUT_FILE="${OUTPUT_FILE#@}"
SORT_ORDER=${SORT_ORDER:-"coordinate"}
VALIDATION_STRINGENCY=${VALIDATION_STRINGENCY:-"LENIENT"}
MAX_RECORDS_IN_RAM=${MAX_RECORDS_IN_RAM:-"500000"}
COMPRESSION_LEVEL=${COMPRESSION_LEVEL:-"5"}
CREATE_INDEX=${CREATE_INDEX:-"true"}

# Create subdirectories in output path if needed
OUTPUT_DIR=$(dirname "/tmp/output/$OUTPUT_FILE")
if [[ "$OUTPUT_DIR" != "/tmp/output" ]]; then
    mkdir -p "$OUTPUT_DIR"
    echo "📁 Created output directory: $OUTPUT_DIR"
fi

case "$COMMAND" in
    "SortSam"|"sort")
        echo "🔄 Running Picard SortSam..."

        PICARD_CMD="picard SortSam"
        PICARD_CMD="$PICARD_CMD INPUT=$INPUT_FILE_PATH"
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

    "MarkDuplicates"|"mark-duplicates")
        echo "🔍 Running Picard MarkDuplicates..."

        METRICS_FILE=${METRICS_FILE:-"duplicate_metrics.txt"}

        PICARD_CMD="picard MarkDuplicates"
        PICARD_CMD="$PICARD_CMD INPUT=$INPUT_FILE_PATH"
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

    "MergeSamFiles"|"merge")
        echo "🔗 Running Picard MergeSamFiles..."

        PICARD_CMD="picard MergeSamFiles"
        PICARD_CMD="$PICARD_CMD INPUT=$INPUT_FILE_PATH"
        # Note: For merging multiple files, add more INPUT parameters
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

    "ValidateSamFile"|"validate")
        echo "✅ Running Picard ValidateSamFile..."

        VALIDATION_REPORT=${OUTPUT_FILE:-"validation_report.txt"}

        PICARD_CMD="picard ValidateSamFile"
        PICARD_CMD="$PICARD_CMD INPUT=$INPUT_FILE_PATH"

        # Add validation mode
        PICARD_CMD="$PICARD_CMD MODE=VERBOSE"

        # Add validation stringency
        if [[ -n "$VALIDATION_STRINGENCY" ]]; then
            PICARD_CMD="$PICARD_CMD VALIDATION_STRINGENCY=$VALIDATION_STRINGENCY"
        fi

        echo "🚀 Executing: $PICARD_CMD"

        # Execute the command (allow non-zero exit for validation warnings)
        if eval "$PICARD_CMD" 2>&1 | tee "/tmp/output/$VALIDATION_REPORT"; then
            echo "✅ Picard ValidateSamFile completed"
        else
            # ValidateSamFile may return non-zero for warnings, which is okay
            echo "⚠️  Picard ValidateSamFile found validation issues (see report)"
        fi

        # Always save the report
        echo "📊 Validation report saved: $VALIDATION_REPORT"
        ls -la /tmp/output/
        ;;

    "CollectAlignmentSummaryMetrics"|"CollectInsertSizeMetrics"|"CollectGcBiasMetrics"|"CollectSequencingArtifactMetrics"|"CollectQualityYieldMetrics"|"CollectRnaSeqMetrics"|"CollectTargetedPcrMetrics"|"CollectWgsMetrics"|"CollectHsMetrics"|"collect-metrics")
        echo "📊 Running Picard $COMMAND..."

        METRICS_FILE=${METRICS_FILE:-"metrics.txt"}

        PICARD_CMD="picard $COMMAND"
        PICARD_CMD="$PICARD_CMD INPUT=$INPUT_FILE_PATH"
        PICARD_CMD="$PICARD_CMD OUTPUT=/tmp/output/$METRICS_FILE"

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
            echo "✅ Picard $COMMAND completed successfully"

            # Display all generated metric files
            echo "📁 Generated metric files:"
            ls -la /tmp/output/
        else
            echo "❌ Picard $COMMAND failed"
            exit 1
        fi
        ;;

    *)
        echo "❌ Unsupported Picard command: $COMMAND"
        echo "Supported commands: ValidateSamFile, SortSam, MarkDuplicates, MergeSamFiles, CollectAlignmentSummaryMetrics, CollectInsertSizeMetrics, and more"
        exit 1
        ;;
esac

echo "🎯 Picard handler completed successfully"