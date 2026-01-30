#!/bin/bash
set -e

# Bowtie2 Handler Script
# Handles various Bowtie2 commands and operations for sequence alignment

echo "🧬 Bowtie2 Handler Started"

# Get the input file (should be downloaded already)
INPUT_FILES=(/tmp/input/*)
if [[ ${#INPUT_FILES[@]} -eq 0 ]]; then
    echo "❌ No input files found in /tmp/input/"
    exit 1
fi

INPUT_FILE_PATH="${INPUT_FILES[0]}"
INPUT_FILENAME=$(basename "$INPUT_FILE_PATH")

echo "📁 Processing file: $INPUT_FILENAME"
echo "🎯 Bowtie2 command: $COMMAND"

# Extract organization and workspace IDs from INPUT_S3_KEY
# Format: organizations/{orgId}/workspaces/{workspaceId}/files/{path}
if [[ -n "$INPUT_S3_KEY" ]]; then
    ORGANIZATION_ID=$(echo "$INPUT_S3_KEY" | cut -d'/' -f2)
    WORKSPACE_ID=$(echo "$INPUT_S3_KEY" | cut -d'/' -f4)
    echo "📂 Organization ID: $ORGANIZATION_ID"
    echo "📂 Workspace ID: $WORKSPACE_ID"
fi

# Parse additional parameters from environment
OUTPUT_FILE=${OUTPUT_FILE:-"output.sam"}
THREADS=${THREADS:-"4"}
REFERENCE_GENOME=${REFERENCE_GENOME:-""}
PRESET=${PRESET:-"sensitive"}
ALIGNMENT_MODE=${ALIGNMENT_MODE:-"end-to-end"}

# Bowtie2-specific parameters
MAX_INSERT_SIZE=${MAX_INSERT_SIZE:-"500"}
MIN_INSERT_SIZE=${MIN_INSERT_SIZE:-"0"}
MAX_VALID_ALIGNMENTS=${MAX_VALID_ALIGNMENTS:-"5"}
SEED_MISMATCHES=${SEED_MISMATCHES:-"0"}
SEED_LENGTH=${SEED_LENGTH:-"22"}
SEED_INTERVAL=${SEED_INTERVAL:-"S,1,1.15"}
MAX_RESEED_ATTEMPTS=${MAX_RESEED_ATTEMPTS:-"2"}
MAX_MISMATCH_QUAL_SUM=${MAX_MISMATCH_QUAL_SUM:-"70"}
MAX_N_CEILING=${MAX_N_CEILING:-"L,0,0.15"}

case "$COMMAND" in
    "build")
        echo "📇 Running Bowtie2 build command..."
        
        if [[ -z "$REFERENCE_GENOME" ]]; then
            echo "❌ Reference genome file required for building index"
            exit 1
        fi
        
        # Copy reference genome to working directory
        cp "$INPUT_FILE_PATH" "/tmp/reference.fasta"
        
        BOWTIE2_CMD="bowtie2-build"
        
        # Add threads parameter
        if [[ -n "$THREADS" ]]; then
            BOWTIE2_CMD="$BOWTIE2_CMD --threads $THREADS"
        fi
        
        # Add reference and index base name
        BOWTIE2_CMD="$BOWTIE2_CMD /tmp/reference.fasta /tmp/output/reference_index"
        
        echo "🚀 Executing: $BOWTIE2_CMD"
        
        if eval "$BOWTIE2_CMD"; then
            echo "✅ Bowtie2 index build completed successfully"
            
            # List generated index files
            echo "📁 Generated index files:"
            ls -la /tmp/output/
        else
            echo "❌ Bowtie2 index build failed"
            exit 1
        fi
        ;;
        
    "align")
        echo "🔍 Running Bowtie2 alignment command..."

        if [[ -z "$REFERENCE_GENOME" ]]; then
            echo "❌ Reference genome index required for alignment"
            exit 1
        fi

        # Strip @ notation and file extension from reference genome
        # Example: @exp1/phix174.fasta -> exp1/phix174
        CLEAN_REFERENCE="${REFERENCE_GENOME#@}"  # Remove @ prefix
        CLEAN_REFERENCE="${CLEAN_REFERENCE%.*}"   # Remove file extension
        echo "📚 Using reference index: $CLEAN_REFERENCE"

        # Check if bowtie2 index exists, if not build it
        INDEX_BASE="/tmp/index/$(basename $CLEAN_REFERENCE)"
        if [[ ! -f "${INDEX_BASE}.1.bt2" ]]; then
            echo "⚠️  Bowtie2 index not found, building it automatically..."

            # Download reference genome FASTA file
            REFERENCE_FASTA="${REFERENCE_GENOME#@}"  # Remove @ prefix (keep extension)
            REFERENCE_S3_KEY="organizations/${ORGANIZATION_ID}/workspaces/${WORKSPACE_ID}/files/${REFERENCE_FASTA}"

            echo "📥 Downloading reference genome: s3://${S3_BUCKET}/${REFERENCE_S3_KEY}"
            mkdir -p /tmp/index

            if aws s3 cp "s3://${S3_BUCKET}/${REFERENCE_S3_KEY}" "/tmp/index/reference.fasta" --no-progress; then
                echo "✅ Reference genome downloaded"

                # Build the index
                echo "🔨 Building bowtie2 index..."
                BUILD_CMD="bowtie2-build"
                if [[ -n "$THREADS" ]]; then
                    BUILD_CMD="$BUILD_CMD --threads $THREADS"
                fi
                BUILD_CMD="$BUILD_CMD /tmp/index/reference.fasta ${INDEX_BASE}"

                echo "🚀 Executing: $BUILD_CMD"
                if eval "$BUILD_CMD"; then
                    echo "✅ Bowtie2 index built successfully"
                else
                    echo "❌ Failed to build bowtie2 index"
                    exit 1
                fi
            else
                echo "❌ Failed to download reference genome from S3"
                exit 1
            fi
        else
            echo "✅ Bowtie2 index already exists"
        fi

        # Copy reads
        cp "$INPUT_FILE_PATH" "/tmp/reads.fastq"

        BOWTIE2_CMD="bowtie2"

        # Add threads parameter
        if [[ -n "$THREADS" ]]; then
            BOWTIE2_CMD="$BOWTIE2_CMD --threads $THREADS"
        fi

        # Add preset parameter
        if [[ -n "$PRESET" ]]; then
            case "$PRESET" in
                "very-fast"|"fast"|"sensitive"|"very-sensitive"|"very-fast-local"|"fast-local"|"sensitive-local"|"very-sensitive-local")
                    BOWTIE2_CMD="$BOWTIE2_CMD --$PRESET"
                    ;;
                *)
                    echo "⚠️ Unknown preset: $PRESET, using default"
                    ;;
            esac
        fi

        # Add alignment mode
        if [[ "$ALIGNMENT_MODE" == "local" ]]; then
            BOWTIE2_CMD="$BOWTIE2_CMD --local"
        fi

        # Add insert size parameters for paired-end
        if [[ -n "$MIN_INSERT_SIZE" ]]; then
            BOWTIE2_CMD="$BOWTIE2_CMD -I $MIN_INSERT_SIZE"
        fi

        if [[ -n "$MAX_INSERT_SIZE" ]]; then
            BOWTIE2_CMD="$BOWTIE2_CMD -X $MAX_INSERT_SIZE"
        fi

        # Add reporting parameters
        if [[ -n "$MAX_VALID_ALIGNMENTS" ]]; then
            BOWTIE2_CMD="$BOWTIE2_CMD -k $MAX_VALID_ALIGNMENTS"
        fi

        # Add seed parameters
        if [[ -n "$SEED_MISMATCHES" ]]; then
            BOWTIE2_CMD="$BOWTIE2_CMD -N $SEED_MISMATCHES"
        fi

        if [[ -n "$SEED_LENGTH" ]]; then
            BOWTIE2_CMD="$BOWTIE2_CMD -L $SEED_LENGTH"
        fi

        # Add alignment parameters
        if [[ -n "$MAX_MISMATCH_QUAL_SUM" ]]; then
            BOWTIE2_CMD="$BOWTIE2_CMD --mp $MAX_MISMATCH_QUAL_SUM"
        fi

        # Add index and reads (use the built index path)
        BOWTIE2_CMD="$BOWTIE2_CMD -x ${INDEX_BASE} -U /tmp/reads.fastq"

        # Create subdirectories in output path if needed
        OUTPUT_DIR=$(dirname "/tmp/output/$OUTPUT_FILE")
        if [[ "$OUTPUT_DIR" != "/tmp/output" ]]; then
            mkdir -p "$OUTPUT_DIR"
            echo "📁 Created output directory: $OUTPUT_DIR"
        fi

        echo "🚀 Executing: $BOWTIE2_CMD -S /tmp/output/$OUTPUT_FILE"

        # Execute the command
        if eval "$BOWTIE2_CMD -S /tmp/output/$OUTPUT_FILE"; then
            echo "✅ Bowtie2 alignment completed successfully"
            
            # Display output file info
            OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILE" 2>/dev/null || echo "unknown")
            echo "📊 Output file: $OUTPUT_FILE ($OUTPUT_SIZE bytes)"
        else
            echo "❌ Bowtie2 alignment failed"
            exit 1
        fi
        ;;
        
    "align-paired")
        echo "🔍 Running Bowtie2 paired-end alignment command..."

        if [[ -z "$REFERENCE_GENOME" ]]; then
            echo "❌ Reference genome index required for alignment"
            exit 1
        fi

        # Strip @ notation and file extension from reference genome
        # Example: @exp1/phix174.fasta -> exp1/phix174
        CLEAN_REFERENCE="${REFERENCE_GENOME#@}"  # Remove @ prefix
        CLEAN_REFERENCE="${CLEAN_REFERENCE%.*}"   # Remove file extension
        echo "📚 Using reference index: $CLEAN_REFERENCE"

        # Check if bowtie2 index exists, if not build it
        INDEX_BASE="/tmp/index/$(basename $CLEAN_REFERENCE)"
        if [[ ! -f "${INDEX_BASE}.1.bt2" ]]; then
            echo "⚠️  Bowtie2 index not found, building it automatically..."

            # Download reference genome FASTA file
            REFERENCE_FASTA="${REFERENCE_GENOME#@}"  # Remove @ prefix (keep extension)
            REFERENCE_S3_KEY="organizations/${ORGANIZATION_ID}/workspaces/${WORKSPACE_ID}/files/${REFERENCE_FASTA}"

            echo "📥 Downloading reference genome: s3://${S3_BUCKET}/${REFERENCE_S3_KEY}"
            mkdir -p /tmp/index

            if aws s3 cp "s3://${S3_BUCKET}/${REFERENCE_S3_KEY}" "/tmp/index/reference.fasta" --no-progress; then
                echo "✅ Reference genome downloaded"

                # Build the index
                echo "🔨 Building bowtie2 index..."
                BUILD_CMD="bowtie2-build"
                if [[ -n "$THREADS" ]]; then
                    BUILD_CMD="$BUILD_CMD --threads $THREADS"
                fi
                BUILD_CMD="$BUILD_CMD /tmp/index/reference.fasta ${INDEX_BASE}"

                echo "🚀 Executing: $BUILD_CMD"
                if eval "$BUILD_CMD"; then
                    echo "✅ Bowtie2 index built successfully"
                else
                    echo "❌ Failed to build bowtie2 index"
                    exit 1
                fi
            else
                echo "❌ Failed to download reference genome from S3"
                exit 1
            fi
        else
            echo "✅ Bowtie2 index already exists"
        fi

        # For paired-end, we expect two input files
        # Copy reads (assuming paired files)
        cp "$INPUT_FILE_PATH" "/tmp/reads1.fastq"
        # Note: In real implementation, we'd handle multiple input files

        BOWTIE2_CMD="bowtie2"

        # Add threads parameter
        if [[ -n "$THREADS" ]]; then
            BOWTIE2_CMD="$BOWTIE2_CMD --threads $THREADS"
        fi

        # Add preset parameter
        if [[ -n "$PRESET" ]]; then
            case "$PRESET" in
                "very-fast"|"fast"|"sensitive"|"very-sensitive"|"very-fast-local"|"fast-local"|"sensitive-local"|"very-sensitive-local")
                    BOWTIE2_CMD="$BOWTIE2_CMD --$PRESET"
                    ;;
                *)
                    echo "⚠️ Unknown preset: $PRESET, using default"
                    ;;
            esac
        fi

        # Add insert size parameters
        if [[ -n "$MIN_INSERT_SIZE" ]]; then
            BOWTIE2_CMD="$BOWTIE2_CMD -I $MIN_INSERT_SIZE"
        fi

        if [[ -n "$MAX_INSERT_SIZE" ]]; then
            BOWTIE2_CMD="$BOWTIE2_CMD -X $MAX_INSERT_SIZE"
        fi

        # Add index and paired reads (use the built index path)
        BOWTIE2_CMD="$BOWTIE2_CMD -x ${INDEX_BASE} -1 /tmp/reads1.fastq -2 /tmp/reads2.fastq"

        # Create subdirectories in output path if needed
        OUTPUT_DIR=$(dirname "/tmp/output/$OUTPUT_FILE")
        if [[ "$OUTPUT_DIR" != "/tmp/output" ]]; then
            mkdir -p "$OUTPUT_DIR"
            echo "📁 Created output directory: $OUTPUT_DIR"
        fi

        echo "🚀 Executing: $BOWTIE2_CMD -S /tmp/output/$OUTPUT_FILE"

        # Execute the command
        if eval "$BOWTIE2_CMD -S /tmp/output/$OUTPUT_FILE"; then
            echo "✅ Bowtie2 paired-end alignment completed successfully"
            
            # Display output file info
            OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILE" 2>/dev/null || echo "unknown")
            echo "📊 Output file: $OUTPUT_FILE ($OUTPUT_SIZE bytes)"
        else
            echo "❌ Bowtie2 paired-end alignment failed"
            exit 1
        fi
        ;;
        
    "inspect")
        echo "🔍 Running Bowtie2 inspect command..."

        if [[ -z "$REFERENCE_GENOME" ]]; then
            echo "❌ Reference genome index required for inspection"
            exit 1
        fi

        # Strip @ notation and file extension from reference genome
        # Example: @exp1/phix174.fasta -> exp1/phix174
        CLEAN_REFERENCE="${REFERENCE_GENOME#@}"  # Remove @ prefix
        CLEAN_REFERENCE="${CLEAN_REFERENCE%.*}"   # Remove file extension
        echo "📚 Using reference index: $CLEAN_REFERENCE"

        BOWTIE2_CMD="bowtie2-inspect"

        # Add index base name (use cleaned reference without @ and extension)
        BOWTIE2_CMD="$BOWTIE2_CMD $CLEAN_REFERENCE"
        
        echo "🚀 Executing: $BOWTIE2_CMD > /tmp/output/$OUTPUT_FILE"
        
        # Execute the command
        if eval "$BOWTIE2_CMD > /tmp/output/$OUTPUT_FILE"; then
            echo "✅ Bowtie2 inspect completed successfully"
            
            # Display output file info
            OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILE" 2>/dev/null || echo "unknown")
            echo "📊 Output file: $OUTPUT_FILE ($OUTPUT_SIZE bytes)"
        else
            echo "❌ Bowtie2 inspect failed"
            exit 1
        fi
        ;;
        
    *)
        echo "❌ Unsupported Bowtie2 command: $COMMAND"
        echo "Supported commands: build, align, align-paired, inspect"
        exit 1
        ;;
esac

echo "🎯 Bowtie2 handler completed successfully"