#!/bin/bash
set -e

# BWA Handler Script
# Handles various BWA commands and operations for sequence alignment

echo "🧬 BWA Handler Started"

# Get the input file (should be downloaded already)
INPUT_FILES=(/tmp/input/*)
if [[ ${#INPUT_FILES[@]} -eq 0 ]]; then
    echo "❌ No input files found in /tmp/input/"
    exit 1
fi

INPUT_FILE_PATH="${INPUT_FILES[0]}"
INPUT_FILENAME=$(basename "$INPUT_FILE_PATH")

echo "📁 Processing file: $INPUT_FILENAME"
echo "🎯 BWA command: $COMMAND"

# Parse additional parameters from environment
OUTPUT_FILE=${OUTPUT_FILE:-"output.sam"}
THREADS=${THREADS:-"4"}
REFERENCE_GENOME=${REFERENCE_GENOME:-""}
ALGORITHM=${ALGORITHM:-"mem"}
INDEX_ALGORITHM=${INDEX_ALGORITHM:-"bwtsw"}

# BWA-specific parameters
MIN_SEED_LENGTH=${MIN_SEED_LENGTH:-"19"}
BAND_WIDTH=${BAND_WIDTH:-"100"}
OFF_DIAGONAL_X_DROPOFF=${OFF_DIAGONAL_X_DROPOFF:-"100"}
INTERNAL_SEEDS_LENGTH_CUTOFF=${INTERNAL_SEEDS_LENGTH_CUTOFF:-"1.5"}
MINIMUM_SCORE_TO_OUTPUT=${MINIMUM_SCORE_TO_OUTPUT:-"30"}
MATCHING_SCORE=${MATCHING_SCORE:-"1"}
MISMATCH_PENALTY=${MISMATCH_PENALTY:-"4"}
GAP_OPEN_PENALTIES=${GAP_OPEN_PENALTIES:-"6,6"}
GAP_EXTENSION_PENALTY=${GAP_EXTENSION_PENALTY:-"1,1"}
CLIPPING_PENALTY=${CLIPPING_PENALTY:-"5,5"}
UNPAIRED_READ_PAIR_PENALTY=${UNPAIRED_READ_PAIR_PENALTY:-"17"}

case "$COMMAND" in
    "index")
        echo "📇 Running BWA index command..."
        
        if [[ -z "$REFERENCE_GENOME" ]]; then
            echo "❌ Reference genome file required for indexing"
            exit 1
        fi
        
        # Copy reference genome to working directory
        cp "$INPUT_FILE_PATH" "/tmp/reference.fasta"
        
        BWA_CMD="bwa index"
        
        # Add algorithm parameter
        if [[ -n "$INDEX_ALGORITHM" ]]; then
            BWA_CMD="$BWA_CMD -a $INDEX_ALGORITHM"
        fi
        
        BWA_CMD="$BWA_CMD /tmp/reference.fasta"
        
        echo "🚀 Executing: $BWA_CMD"
        
        if eval "$BWA_CMD"; then
            echo "✅ BWA index completed successfully"
            
            # Copy all index files to output
            cp /tmp/reference.fasta* /tmp/output/
            
            # List generated files
            echo "📁 Generated index files:"
            ls -la /tmp/output/
        else
            echo "❌ BWA index failed"
            exit 1
        fi
        ;;
        
    "mem")
        echo "🔍 Running BWA mem command..."
        
        if [[ -z "$REFERENCE_GENOME" ]]; then
            echo "❌ Reference genome required for alignment"
            exit 1
        fi
        
        # For mem, we need a pre-built index
        # Copy reference and reads
        cp "$INPUT_FILE_PATH" "/tmp/reads.fastq"
        
        BWA_CMD="bwa mem"
        
        # Add threads parameter
        if [[ -n "$THREADS" ]]; then
            BWA_CMD="$BWA_CMD -t $THREADS"
        fi
        
        # Add BWA mem specific parameters
        if [[ -n "$MIN_SEED_LENGTH" ]]; then
            BWA_CMD="$BWA_CMD -k $MIN_SEED_LENGTH"
        fi
        
        if [[ -n "$BAND_WIDTH" ]]; then
            BWA_CMD="$BWA_CMD -w $BAND_WIDTH"
        fi
        
        if [[ -n "$OFF_DIAGONAL_X_DROPOFF" ]]; then
            BWA_CMD="$BWA_CMD -d $OFF_DIAGONAL_X_DROPOFF"
        fi
        
        if [[ -n "$MINIMUM_SCORE_TO_OUTPUT" ]]; then
            BWA_CMD="$BWA_CMD -T $MINIMUM_SCORE_TO_OUTPUT"
        fi
        
        if [[ -n "$MATCHING_SCORE" ]]; then
            BWA_CMD="$BWA_CMD -A $MATCHING_SCORE"
        fi
        
        if [[ -n "$MISMATCH_PENALTY" ]]; then
            BWA_CMD="$BWA_CMD -B $MISMATCH_PENALTY"
        fi
        
        # Add reference and reads
        BWA_CMD="$BWA_CMD $REFERENCE_GENOME /tmp/reads.fastq"
        
        echo "🚀 Executing: $BWA_CMD > /tmp/output/$OUTPUT_FILE"
        
        # Execute the command
        if eval "$BWA_CMD > /tmp/output/$OUTPUT_FILE"; then
            echo "✅ BWA mem alignment completed successfully"
            
            # Display output file info
            OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILE" 2>/dev/null || echo "unknown")
            echo "📊 Output file: $OUTPUT_FILE ($OUTPUT_SIZE bytes)"
        else
            echo "❌ BWA mem alignment failed"
            exit 1
        fi
        ;;
        
    "aln")
        echo "🔍 Running BWA aln command..."
        
        if [[ -z "$REFERENCE_GENOME" ]]; then
            echo "❌ Reference genome required for alignment"
            exit 1
        fi
        
        # Copy reads
        cp "$INPUT_FILE_PATH" "/tmp/reads.fastq"
        
        BWA_CMD="bwa aln"
        
        # Add threads parameter
        if [[ -n "$THREADS" ]]; then
            BWA_CMD="$BWA_CMD -t $THREADS"
        fi
        
        # Add reference and reads
        BWA_CMD="$BWA_CMD $REFERENCE_GENOME /tmp/reads.fastq"
        
        echo "🚀 Executing: $BWA_CMD > /tmp/output/$OUTPUT_FILE"
        
        # Execute the command
        if eval "$BWA_CMD > /tmp/output/$OUTPUT_FILE"; then
            echo "✅ BWA aln alignment completed successfully"
            
            # Display output file info
            OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILE" 2>/dev/null || echo "unknown")
            echo "📊 Output file: $OUTPUT_FILE ($OUTPUT_SIZE bytes)"
        else
            echo "❌ BWA aln alignment failed"
            exit 1
        fi
        ;;
        
    "samse")
        echo "🔄 Running BWA samse command..."
        
        if [[ -z "$REFERENCE_GENOME" ]]; then
            echo "❌ Reference genome required for samse"
            exit 1
        fi
        
        # Copy alignment file (from aln step)
        cp "$INPUT_FILE_PATH" "/tmp/alignment.sai"
        
        BWA_CMD="bwa samse $REFERENCE_GENOME /tmp/alignment.sai /tmp/reads.fastq"
        
        echo "🚀 Executing: $BWA_CMD > /tmp/output/$OUTPUT_FILE"
        
        # Execute the command
        if eval "$BWA_CMD > /tmp/output/$OUTPUT_FILE"; then
            echo "✅ BWA samse completed successfully"
            
            # Display output file info
            OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILE" 2>/dev/null || echo "unknown")
            echo "📊 Output file: $OUTPUT_FILE ($OUTPUT_SIZE bytes)"
        else
            echo "❌ BWA samse failed"
            exit 1
        fi
        ;;
        
    "sampe")
        echo "🔄 Running BWA sampe command (paired-end)..."
        
        if [[ -z "$REFERENCE_GENOME" ]]; then
            echo "❌ Reference genome required for sampe"
            exit 1
        fi
        
        # Copy alignment files (from aln step for paired reads)
        cp "$INPUT_FILE_PATH" "/tmp/alignment1.sai"
        
        BWA_CMD="bwa sampe $REFERENCE_GENOME /tmp/alignment1.sai /tmp/alignment2.sai /tmp/reads1.fastq /tmp/reads2.fastq"
        
        echo "🚀 Executing: $BWA_CMD > /tmp/output/$OUTPUT_FILE"
        
        # Execute the command
        if eval "$BWA_CMD > /tmp/output/$OUTPUT_FILE"; then
            echo "✅ BWA sampe completed successfully"
            
            # Display output file info
            OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILE" 2>/dev/null || echo "unknown")
            echo "📊 Output file: $OUTPUT_FILE ($OUTPUT_SIZE bytes)"
        else
            echo "❌ BWA sampe failed"
            exit 1
        fi
        ;;
        
    "bwasw")
        echo "🔍 Running BWA bwasw command (long sequences)..."
        
        if [[ -z "$REFERENCE_GENOME" ]]; then
            echo "❌ Reference genome required for bwasw"
            exit 1
        fi
        
        # Copy reads
        cp "$INPUT_FILE_PATH" "/tmp/reads.fastq"
        
        BWA_CMD="bwa bwasw"
        
        # Add threads parameter
        if [[ -n "$THREADS" ]]; then
            BWA_CMD="$BWA_CMD -t $THREADS"
        fi
        
        # Add reference and reads
        BWA_CMD="$BWA_CMD $REFERENCE_GENOME /tmp/reads.fastq"
        
        echo "🚀 Executing: $BWA_CMD > /tmp/output/$OUTPUT_FILE"
        
        # Execute the command
        if eval "$BWA_CMD > /tmp/output/$OUTPUT_FILE"; then
            echo "✅ BWA bwasw alignment completed successfully"
            
            # Display output file info
            OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILE" 2>/dev/null || echo "unknown")
            echo "📊 Output file: $OUTPUT_FILE ($OUTPUT_SIZE bytes)"
        else
            echo "❌ BWA bwasw alignment failed"
            exit 1
        fi
        ;;
        
    "mem-paired")
        echo "🔍 Running BWA mem command (paired-end)..."
        
        if [[ -z "$REFERENCE_GENOME" ]]; then
            echo "❌ Reference genome required for alignment"
            exit 1
        fi
        
        # Handle paired-end reads
        cp "$INPUT_FILE_PATH" "/tmp/reads1.fastq"
        # Note: In real implementation, would handle second file from environment
        
        BWA_CMD="bwa mem"
        
        # Add threads parameter
        if [[ -n "$THREADS" ]]; then
            BWA_CMD="$BWA_CMD -t $THREADS"
        fi
        
        # Add BWA mem specific parameters
        if [[ -n "$MIN_SEED_LENGTH" ]]; then
            BWA_CMD="$BWA_CMD -k $MIN_SEED_LENGTH"
        fi
        
        if [[ -n "$BAND_WIDTH" ]]; then
            BWA_CMD="$BWA_CMD -w $BAND_WIDTH"
        fi
        
        # Mark shorter split hits as secondary (recommended)
        BWA_CMD="$BWA_CMD -M"
        
        # Add read group if specified
        if [[ -n "$READ_GROUP" ]]; then
            BWA_CMD="$BWA_CMD -R '$READ_GROUP'"
        fi
        
        # Add reference and paired reads
        BWA_CMD="$BWA_CMD $REFERENCE_GENOME /tmp/reads1.fastq /tmp/reads2.fastq"
        
        echo "🚀 Executing: $BWA_CMD > /tmp/output/$OUTPUT_FILE"
        
        # Execute the command
        if eval "$BWA_CMD > /tmp/output/$OUTPUT_FILE"; then
            echo "✅ BWA mem paired-end alignment completed successfully"
            
            # Display output file info
            OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILE" 2>/dev/null || echo "unknown")
            echo "📊 Output file: $OUTPUT_FILE ($OUTPUT_SIZE bytes)"
        else
            echo "❌ BWA mem paired-end alignment failed"
            exit 1
        fi
        ;;
        
    "mem-pacbio")
        echo "🔍 Running BWA mem command (PacBio preset)..."
        
        if [[ -z "$REFERENCE_GENOME" ]]; then
            echo "❌ Reference genome required for alignment"
            exit 1
        fi
        
        cp "$INPUT_FILE_PATH" "/tmp/reads.fastq"
        
        BWA_CMD="bwa mem"
        
        # PacBio specific parameters (adapted for long, error-prone reads)
        BWA_CMD="$BWA_CMD -x pacbio"
        
        # Add threads parameter
        if [[ -n "$THREADS" ]]; then
            BWA_CMD="$BWA_CMD -t $THREADS"
        fi
        
        # Add reference and reads
        BWA_CMD="$BWA_CMD $REFERENCE_GENOME /tmp/reads.fastq"
        
        echo "🚀 Executing: $BWA_CMD > /tmp/output/$OUTPUT_FILE"
        
        # Execute the command
        if eval "$BWA_CMD > /tmp/output/$OUTPUT_FILE"; then
            echo "✅ BWA mem PacBio alignment completed successfully"
            
            # Display output file info
            OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILE" 2>/dev/null || echo "unknown")
            echo "📊 Output file: $OUTPUT_FILE ($OUTPUT_SIZE bytes)"
        else
            echo "❌ BWA mem PacBio alignment failed"
            exit 1
        fi
        ;;
        
    "mem-ont2d")
        echo "🔍 Running BWA mem command (Oxford Nanopore 2D preset)..."
        
        if [[ -z "$REFERENCE_GENOME" ]]; then
            echo "❌ Reference genome required for alignment"
            exit 1
        fi
        
        cp "$INPUT_FILE_PATH" "/tmp/reads.fastq"
        
        BWA_CMD="bwa mem"
        
        # Oxford Nanopore 2D specific parameters
        BWA_CMD="$BWA_CMD -x ont2d"
        
        # Add threads parameter
        if [[ -n "$THREADS" ]]; then
            BWA_CMD="$BWA_CMD -t $THREADS"
        fi
        
        # Add reference and reads
        BWA_CMD="$BWA_CMD $REFERENCE_GENOME /tmp/reads.fastq"
        
        echo "🚀 Executing: $BWA_CMD > /tmp/output/$OUTPUT_FILE"
        
        # Execute the command
        if eval "$BWA_CMD > /tmp/output/$OUTPUT_FILE"; then
            echo "✅ BWA mem Oxford Nanopore alignment completed successfully"
            
            # Display output file info
            OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILE" 2>/dev/null || echo "unknown")
            echo "📊 Output file: $OUTPUT_FILE ($OUTPUT_SIZE bytes)"
        else
            echo "❌ BWA mem Oxford Nanopore alignment failed"
            exit 1
        fi
        ;;
        
    *)
        echo "❌ Unsupported BWA command: $COMMAND"
        echo "Supported commands: index, mem, aln, samse, sampe, bwasw, mem-paired, mem-pacbio, mem-ont2d"
        exit 1
        ;;
esac

echo "🎯 BWA handler completed successfully"