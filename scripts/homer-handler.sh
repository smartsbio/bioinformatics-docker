#!/bin/bash
set -e

# HOMER Handler Script
# Handles motif discovery and ChIP-Seq analysis

echo "🧬 HOMER Handler Started"

# Get the input file (should be downloaded already)
INPUT_FILES=(/tmp/input/*)
if [[ ${#INPUT_FILES[@]} -eq 0 ]]; then
    echo "❌ No input files found in /tmp/input/"
    exit 1
fi

INPUT_FILE_PATH="${INPUT_FILES[0]}"
INPUT_FILENAME=$(basename "$INPUT_FILE_PATH")

echo "📁 Processing file: $INPUT_FILENAME"
echo "🎯 HOMER command: $COMMAND"

# Extract organization and workspace IDs from INPUT_S3_KEY
# Format: organizations/{orgId}/workspaces/{workspaceId}/files/{path}
if [[ -n "$INPUT_S3_KEY" ]]; then
    ORGANIZATION_ID=$(echo "$INPUT_S3_KEY" | cut -d'/' -f2)
    WORKSPACE_ID=$(echo "$INPUT_S3_KEY" | cut -d'/' -f4)
    echo "📂 Organization ID: $ORGANIZATION_ID"
    echo "📂 Workspace ID: $WORKSPACE_ID"
fi

# Parse additional parameters from environment
OUTPUT_DIR=${OUTPUT_DIR:-"/tmp/output/homer_results"}
# Strip @ notation from output directory path
OUTPUT_DIR="${OUTPUT_DIR#@}"
# Ensure OUTPUT_DIR starts with /tmp/output/
if [[ "$OUTPUT_DIR" != /tmp/output/* ]]; then
    OUTPUT_DIR="/tmp/output/$OUTPUT_DIR"
fi
GENOME=${GENOME:-"hg38"}
MOTIF_LENGTH=${MOTIF_LENGTH:-"8,10,12"}
NUM_MOTIFS=${NUM_MOTIFS:-"25"}
SIZE=${SIZE:-"200"}
MASK=${MASK:-"true"}
CPUS=${CPUS:-"4"}

# Create output directory
mkdir -p "$OUTPUT_DIR"

case "$COMMAND" in
    "find-motifs-genome")
        echo "🔍 Running HOMER motif discovery from genome positions..."
        
        # Copy input file (should be peak file in BED format)
        cp "$INPUT_FILE_PATH" "/tmp/peaks.bed"
        
        HOMER_CMD="findMotifsGenome.pl"
        HOMER_CMD="$HOMER_CMD /tmp/peaks.bed"
        HOMER_CMD="$HOMER_CMD $GENOME"
        HOMER_CMD="$HOMER_CMD $OUTPUT_DIR"
        
        # Add motif length parameter
        if [[ -n "$MOTIF_LENGTH" ]]; then
            HOMER_CMD="$HOMER_CMD -len $MOTIF_LENGTH"
        fi
        
        # Add size parameter
        if [[ -n "$SIZE" ]]; then
            HOMER_CMD="$HOMER_CMD -size $SIZE"
        fi
        
        # Add mask parameter
        if [[ "$MASK" == "true" ]]; then
            HOMER_CMD="$HOMER_CMD -mask"
        fi
        
        # Add CPU parameter
        if [[ -n "$CPUS" ]]; then
            HOMER_CMD="$HOMER_CMD -p $CPUS"
        fi
        
        # Add number of motifs to find
        if [[ -n "$NUM_MOTIFS" ]]; then
            HOMER_CMD="$HOMER_CMD -S $NUM_MOTIFS"
        fi
        
        echo "🚀 Executing: $HOMER_CMD"
        
        if eval "$HOMER_CMD"; then
            echo "✅ HOMER motif discovery completed successfully"
            
            # Display results
            echo "📁 Generated files:"
            ls -la "$OUTPUT_DIR/"
            
            # Show motif results if available
            if [[ -f "$OUTPUT_DIR/homerResults.html" ]]; then
                echo "📊 Results available in: homerResults.html"
            fi
        else
            echo "❌ HOMER motif discovery failed"
            exit 1
        fi
        ;;
        
    "find-motifs-sequences")
        echo "🔍 Running HOMER motif discovery from sequences..."
        
        # Copy input file (should be FASTA sequences)
        cp "$INPUT_FILE_PATH" "/tmp/sequences.fasta"
        
        HOMER_CMD="findMotifs.pl"
        HOMER_CMD="$HOMER_CMD /tmp/sequences.fasta"
        HOMER_CMD="$HOMER_CMD fasta"
        HOMER_CMD="$HOMER_CMD $OUTPUT_DIR"
        
        # Add motif length parameter
        if [[ -n "$MOTIF_LENGTH" ]]; then
            HOMER_CMD="$HOMER_CMD -len $MOTIF_LENGTH"
        fi
        
        # Add CPU parameter
        if [[ -n "$CPUS" ]]; then
            HOMER_CMD="$HOMER_CMD -p $CPUS"
        fi
        
        # Add number of motifs to find
        if [[ -n "$NUM_MOTIFS" ]]; then
            HOMER_CMD="$HOMER_CMD -S $NUM_MOTIFS"
        fi
        
        echo "🚀 Executing: $HOMER_CMD"
        
        if eval "$HOMER_CMD"; then
            echo "✅ HOMER sequence motif discovery completed successfully"
            
            echo "📁 Generated files:"
            ls -la "$OUTPUT_DIR/"
        else
            echo "❌ HOMER sequence motif discovery failed"
            exit 1
        fi
        ;;
        
    "annotate-peaks")
        echo "📝 Running HOMER peak annotation..."
        
        # Copy input file (peak file)
        cp "$INPUT_FILE_PATH" "/tmp/peaks.bed"
        
        ANNOTATION_FILE="/tmp/output/annotated_peaks.txt"
        
        HOMER_CMD="annotatePeaks.pl"
        HOMER_CMD="$HOMER_CMD /tmp/peaks.bed"
        HOMER_CMD="$HOMER_CMD $GENOME"
        
        # Add additional parameters
        if [[ -n "$SIZE" ]]; then
            HOMER_CMD="$HOMER_CMD -size $SIZE"
        fi
        
        echo "🚀 Executing: $HOMER_CMD > $ANNOTATION_FILE"
        
        if eval "$HOMER_CMD > $ANNOTATION_FILE"; then
            echo "✅ HOMER peak annotation completed successfully"
            
            # Display output file info
            OUTPUT_SIZE=$(stat -c%s "$ANNOTATION_FILE" 2>/dev/null || echo "unknown")
            echo "📊 Annotated peaks: $(basename $ANNOTATION_FILE) ($OUTPUT_SIZE bytes)"
        else
            echo "❌ HOMER peak annotation failed"
            exit 1
        fi
        ;;
        
    "make-tag-directory")
        echo "📚 Running HOMER tag directory creation..."
        
        # Copy input file (alignment file)
        cp "$INPUT_FILE_PATH" "/tmp/alignments.bam"
        
        TAG_DIR="/tmp/output/tag_directory"
        
        HOMER_CMD="makeTagDirectory"
        HOMER_CMD="$HOMER_CMD $TAG_DIR"
        HOMER_CMD="$HOMER_CMD /tmp/alignments.bam"
        
        # Add genome parameter
        if [[ -n "$GENOME" ]]; then
            HOMER_CMD="$HOMER_CMD -genome $GENOME"
        fi
        
        # Add checkGC parameter
        HOMER_CMD="$HOMER_CMD -checkGC"
        
        echo "🚀 Executing: $HOMER_CMD"
        
        if eval "$HOMER_CMD"; then
            echo "✅ HOMER tag directory created successfully"
            
            echo "📁 Tag directory contents:"
            ls -la "$TAG_DIR/"
        else
            echo "❌ HOMER tag directory creation failed"
            exit 1
        fi
        ;;
        
    "find-peaks")
        echo "🏔️ Running HOMER peak finding..."
        
        # For peak finding, we need tag directories
        # This is a simplified version
        TAG_DIR="/tmp/input_tags"
        CONTROL_DIR="/tmp/control_tags"
        PEAKS_FILE="/tmp/output/peaks.txt"
        
        # Copy input as tag directory (simplified)
        mkdir -p "$TAG_DIR"
        cp "$INPUT_FILE_PATH" "$TAG_DIR/"
        
        HOMER_CMD="findPeaks"
        HOMER_CMD="$HOMER_CMD $TAG_DIR"
        HOMER_CMD="$HOMER_CMD -style factor"
        HOMER_CMD="$HOMER_CMD -o $PEAKS_FILE"
        
        # Add size parameter
        if [[ -n "$SIZE" ]]; then
            HOMER_CMD="$HOMER_CMD -size $SIZE"
        fi
        
        # Add minimum distance
        HOMER_CMD="$HOMER_CMD -minDist 150"
        
        echo "🚀 Executing: $HOMER_CMD"
        
        if eval "$HOMER_CMD"; then
            echo "✅ HOMER peak finding completed successfully"
            
            # Display output file info
            OUTPUT_SIZE=$(stat -c%s "$PEAKS_FILE" 2>/dev/null || echo "unknown")
            echo "📊 Peaks found: $(basename $PEAKS_FILE) ($OUTPUT_SIZE bytes)"
        else
            echo "❌ HOMER peak finding failed"
            exit 1
        fi
        ;;
        
    *)
        echo "❌ Unsupported HOMER command: $COMMAND"
        echo "Supported commands: find-motifs-genome, find-motifs-sequences, annotate-peaks, make-tag-directory, find-peaks"
        exit 1
        ;;
esac

echo "🎯 HOMER handler completed successfully"