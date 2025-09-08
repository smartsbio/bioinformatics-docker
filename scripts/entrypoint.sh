#!/bin/bash
set -e

# Bioinformatics Process Entrypoint
# This script handles the main processing workflow for bioinformatics tools

echo "🧬 Starting bioinformatics process: $(date)"

# Environment Variables (set by ECS task)
TOOL_NAME=${TOOL_NAME:-"samtools"}
COMMAND=${COMMAND:-"view"}
INPUT_FILE=${INPUT_FILE:-""}
OUTPUT_FILE=${OUTPUT_FILE:-""}
S3_BUCKET=${S3_BUCKET:-""}
INPUT_S3_KEY=${INPUT_S3_KEY:-""}
OUTPUT_PATH=${OUTPUT_PATH:-""}

# Process ID for logging
PROCESS_ID=${PROCESS_ID:-"unknown"}

echo "🔧 Tool: $TOOL_NAME"
echo "📝 Command: $COMMAND"
echo "📂 Process ID: $PROCESS_ID"

# Validate required environment variables
if [[ -z "$S3_BUCKET" || -z "$INPUT_S3_KEY" || -z "$OUTPUT_PATH" ]]; then
    echo "❌ Error: Missing required environment variables"
    echo "Required: S3_BUCKET, INPUT_S3_KEY, OUTPUT_PATH"
    exit 1
fi

# Create working directories
mkdir -p /tmp/input /tmp/output

echo "📥 Downloading input files..."

# Download input file from S3
if [[ -n "$INPUT_S3_KEY" ]]; then
    INPUT_FILENAME=$(basename "$INPUT_S3_KEY")
    echo "📥 Downloading: s3://$S3_BUCKET/$INPUT_S3_KEY -> /tmp/input/$INPUT_FILENAME"
    
    if ! aws s3 cp "s3://$S3_BUCKET/$INPUT_S3_KEY" "/tmp/input/$INPUT_FILENAME"; then
        echo "❌ Failed to download input file from S3"
        exit 1
    fi
    
    echo "✅ Input file downloaded successfully"
else
    echo "❌ No input file specified"
    exit 1
fi

# Change to input directory for processing
cd /tmp/input

echo "🔬 Running bioinformatics tool: $TOOL_NAME $COMMAND"

# Route to appropriate tool handler
case "$TOOL_NAME" in
    "samtools")
        /usr/local/bin/scripts/samtools-handler.sh
        ;;
    *)
        echo "❌ Unsupported tool: $TOOL_NAME"
        exit 1
        ;;
esac

echo "📤 Uploading output files to S3..."

# Upload all files from output directory to S3
if [[ -d "/tmp/output" ]] && [[ -n "$(ls -A /tmp/output)" ]]; then
    for file in /tmp/output/*; do
        if [[ -f "$file" ]]; then
            filename=$(basename "$file")
            s3_output_key="${OUTPUT_PATH}${filename}"
            
            echo "📤 Uploading: $file -> s3://$S3_BUCKET/$s3_output_key"
            
            if ! aws s3 cp "$file" "s3://$S3_BUCKET/$s3_output_key"; then
                echo "❌ Failed to upload $filename to S3"
                exit 1
            fi
            
            echo "✅ Uploaded: $filename"
        fi
    done
else
    echo "❌ No output files found to upload"
    exit 1
fi

echo "🎉 Process completed successfully at $(date)"
echo "📊 Output files uploaded to: s3://$S3_BUCKET/$OUTPUT_PATH"