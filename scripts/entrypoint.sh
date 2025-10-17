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
    
    if ! aws s3 cp "s3://$S3_BUCKET/$INPUT_S3_KEY" "/tmp/input/$INPUT_FILENAME" --no-progress; then
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
    "bwa")
        /usr/local/bin/scripts/bwa-handler.sh
        ;;
    "bowtie2")
        /usr/local/bin/scripts/bowtie2-handler.sh
        ;;
    "bedtools")
        /usr/local/bin/scripts/bedtools-handler.sh
        ;;
    "fastqc")
        /usr/local/bin/scripts/fastqc-handler.sh
        ;;
    "trimmomatic")
        /usr/local/bin/scripts/trimmomatic-handler.sh
        ;;
    "picard")
        /usr/local/bin/scripts/picard-handler.sh
        ;;
    "vcftools")
        /usr/local/bin/scripts/vcftools-handler.sh
        ;;
    "fmlrc")
        /usr/local/bin/scripts/fmlrc-handler.sh
        ;;
    "homer")
        /usr/local/bin/scripts/homer-handler.sh
        ;;
    "annovar")
        /usr/local/bin/scripts/annovar-handler.sh
        ;;
    "genomicranges")
        /usr/local/bin/scripts/genomicranges-handler.sh
        ;;
    *)
        echo "❌ Unsupported tool: $TOOL_NAME"
        echo "Supported tools: samtools, bwa, bowtie2, bedtools, fastqc, trimmomatic, picard, vcftools, fmlrc, homer, annovar, genomicranges"
        exit 1
        ;;
esac

echo "📤 Uploading output files to S3..."

# Array to collect file metadata for JSON output
declare -a uploaded_files=()

# Upload all files from output directory to S3
if [[ -d "/tmp/output" ]] && [[ -n "$(ls -A /tmp/output)" ]]; then
    for file in /tmp/output/*; do
        if [[ -f "$file" ]]; then
            filename=$(basename "$file")
            # Ensure OUTPUT_PATH ends with / before concatenating filename
            OUTPUT_PATH_NORMALIZED="${OUTPUT_PATH%/}/"  # Remove trailing / if exists, then add it back
            s3_output_key="${OUTPUT_PATH_NORMALIZED}${filename}"

            echo "📤 Uploading: $file -> s3://$S3_BUCKET/$s3_output_key"

            # Get file size and checksum before upload
            file_size=$(stat -c%s "$file" 2>/dev/null || stat -f%z "$file" 2>/dev/null || echo "0")
            file_md5=$(md5sum "$file" 2>/dev/null | awk '{print $1}' || md5 -q "$file" 2>/dev/null || echo "unknown")

            if ! aws s3 cp "$file" "s3://$S3_BUCKET/$s3_output_key" --no-progress; then
                echo "❌ Failed to upload $filename to S3"
                exit 1
            fi

            echo "✅ Uploaded: $filename (${file_size} bytes)"

            # Collect file metadata for JSON output
            uploaded_files+=("{\"filename\":\"$filename\",\"s3Key\":\"$s3_output_key\",\"s3Bucket\":\"$S3_BUCKET\",\"size\":$file_size,\"checksum\":\"$file_md5\"}")
        fi
    done
else
    echo "❌ No output files found to upload"
    exit 1
fi

echo "🎉 Process completed successfully at $(date)"
echo "📊 Output files uploaded to: s3://$S3_BUCKET/$OUTPUT_PATH"

# Output structured JSON for Process Manager to parse
# This special marker helps Process Manager identify the JSON output
echo "===SMARTS_BIO_OUTPUT_FILES_JSON_START==="
echo -n "["
for i in "${!uploaded_files[@]}"; do
    echo -n "${uploaded_files[$i]}"
    if [ $i -lt $((${#uploaded_files[@]} - 1)) ]; then
        echo -n ","
    fi
done
echo "]"
echo "===SMARTS_BIO_OUTPUT_FILES_JSON_END==="