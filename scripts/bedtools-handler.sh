#!/bin/bash
set -e

# BEDTools Handler Script
# Handles various BEDTools commands for genomic intervals analysis

echo "🧬 BEDTools Handler Started"

# Get the input file (should be downloaded already)
INPUT_FILES=(/tmp/input/*)
if [[ ${#INPUT_FILES[@]} -eq 0 ]]; then
    echo "❌ No input files found in /tmp/input/"
    exit 1
fi

INPUT_FILE_PATH="${INPUT_FILES[0]}"
INPUT_FILENAME=$(basename "$INPUT_FILE_PATH")

echo "📁 Processing file: $INPUT_FILENAME"
echo "🎯 BEDTools command: $COMMAND"

# Extract organization and workspace IDs from INPUT_S3_KEY
# Format: organizations/{orgId}/workspaces/{workspaceId}/files/{path}
if [[ -n "$INPUT_S3_KEY" ]]; then
    ORGANIZATION_ID=$(echo "$INPUT_S3_KEY" | cut -d'/' -f2)
    WORKSPACE_ID=$(echo "$INPUT_S3_KEY" | cut -d'/' -f4)
    echo "📂 Organization ID: $ORGANIZATION_ID"
    echo "📂 Workspace ID: $WORKSPACE_ID"
fi

# Parse additional parameters from environment
OUTPUT_FILE=${OUTPUT_FILE:-"output.bed"}
# Strip @ notation from output file path
OUTPUT_FILE="${OUTPUT_FILE#@}"
OVERLAP_FRACTION=${OVERLAP_FRACTION:-"0.0"}
MINIMUM_OVERLAP=${MINIMUM_OVERLAP:-"1"}
STRAND_SPECIFIC=${STRAND_SPECIFIC:-"false"}
RECIPROCAL_OVERLAP=${RECIPROCAL_OVERLAP:-"false"}

case "$COMMAND" in
    "intersect")
        echo "🔍 Running BEDTools intersect command..."
        
        # Copy input files
        cp "$INPUT_FILE_PATH" "/tmp/input_intervals.bed"
        
        BEDTOOLS_CMD="bedtools intersect"
        
        # Add input file
        BEDTOOLS_CMD="$BEDTOOLS_CMD -a /tmp/input_intervals.bed"
        
        # Add second file (for intersection)
        # Note: In real implementation, we'd handle multiple input files
        BEDTOOLS_CMD="$BEDTOOLS_CMD -b /tmp/input_intervals.bed"
        
        # Add overlap parameters
        if [[ "$OVERLAP_FRACTION" != "0.0" ]]; then
            BEDTOOLS_CMD="$BEDTOOLS_CMD -f $OVERLAP_FRACTION"
        fi
        
        if [[ "$MINIMUM_OVERLAP" != "1" ]]; then
            BEDTOOLS_CMD="$BEDTOOLS_CMD -F $MINIMUM_OVERLAP"
        fi
        
        if [[ "$RECIPROCAL_OVERLAP" == "true" ]]; then
            BEDTOOLS_CMD="$BEDTOOLS_CMD -r"
        fi
        
        if [[ "$STRAND_SPECIFIC" == "true" ]]; then
            BEDTOOLS_CMD="$BEDTOOLS_CMD -s"
        fi
        
        # Create subdirectories in output path if needed
        OUTPUT_DIR=$(dirname "/tmp/output/$OUTPUT_FILE")
        if [[ "$OUTPUT_DIR" != "/tmp/output" ]]; then
            mkdir -p "$OUTPUT_DIR"
            echo "📁 Created output directory: $OUTPUT_DIR"
        fi

        echo "🚀 Executing: $BEDTOOLS_CMD > /tmp/output/$OUTPUT_FILE"

        # Execute the command
        if eval "$BEDTOOLS_CMD > /tmp/output/$OUTPUT_FILE"; then
            echo "✅ BEDTools intersect completed successfully"
            
            # Display output file info
            OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILE" 2>/dev/null || echo "unknown")
            echo "📊 Output file: $OUTPUT_FILE ($OUTPUT_SIZE bytes)"
        else
            echo "❌ BEDTools intersect failed"
            exit 1
        fi
        ;;
        
    "merge")
        echo "🔄 Running BEDTools merge command..."
        
        # Copy input file
        cp "$INPUT_FILE_PATH" "/tmp/input_intervals.bed"
        
        BEDTOOLS_CMD="bedtools merge -i /tmp/input_intervals.bed"
        
        if [[ "$STRAND_SPECIFIC" == "true" ]]; then
            BEDTOOLS_CMD="$BEDTOOLS_CMD -s"
        fi
        
        # Create subdirectories in output path if needed
        OUTPUT_DIR=$(dirname "/tmp/output/$OUTPUT_FILE")
        if [[ "$OUTPUT_DIR" != "/tmp/output" ]]; then
            mkdir -p "$OUTPUT_DIR"
            echo "📁 Created output directory: $OUTPUT_DIR"
        fi

        echo "🚀 Executing: $BEDTOOLS_CMD > /tmp/output/$OUTPUT_FILE"

        # Execute the command
        if eval "$BEDTOOLS_CMD > /tmp/output/$OUTPUT_FILE"; then
            echo "✅ BEDTools merge completed successfully"
            
            # Display output file info
            OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILE" 2>/dev/null || echo "unknown")
            echo "📊 Output file: $OUTPUT_FILE ($OUTPUT_SIZE bytes)"
        else
            echo "❌ BEDTools merge failed"
            exit 1
        fi
        ;;
        
    "sort")
        echo "🔄 Running BEDTools sort command..."
        
        # Copy input file
        cp "$INPUT_FILE_PATH" "/tmp/input_intervals.bed"
        
        BEDTOOLS_CMD="bedtools sort -i /tmp/input_intervals.bed"
        
        # Create subdirectories in output path if needed
        OUTPUT_DIR=$(dirname "/tmp/output/$OUTPUT_FILE")
        if [[ "$OUTPUT_DIR" != "/tmp/output" ]]; then
            mkdir -p "$OUTPUT_DIR"
            echo "📁 Created output directory: $OUTPUT_DIR"
        fi

        echo "🚀 Executing: $BEDTOOLS_CMD > /tmp/output/$OUTPUT_FILE"

        # Execute the command
        if eval "$BEDTOOLS_CMD > /tmp/output/$OUTPUT_FILE"; then
            echo "✅ BEDTools sort completed successfully"
            
            # Display output file info
            OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILE" 2>/dev/null || echo "unknown")
            echo "📊 Output file: $OUTPUT_FILE ($OUTPUT_SIZE bytes)"
        else
            echo "❌ BEDTools sort failed"
            exit 1
        fi
        ;;
        
    "coverage")
        echo "📊 Running BEDTools coverage command..."
        
        # Copy input files
        cp "$INPUT_FILE_PATH" "/tmp/input_intervals.bed"
        
        BEDTOOLS_CMD="bedtools coverage"
        
        # Add input file
        BEDTOOLS_CMD="$BEDTOOLS_CMD -a /tmp/input_intervals.bed"
        
        # Add second file (for coverage calculation)
        # Note: In real implementation, we'd handle multiple input files
        BEDTOOLS_CMD="$BEDTOOLS_CMD -b /tmp/input_intervals.bed"
        
        if [[ "$STRAND_SPECIFIC" == "true" ]]; then
            BEDTOOLS_CMD="$BEDTOOLS_CMD -s"
        fi
        
        # Create subdirectories in output path if needed
        OUTPUT_DIR=$(dirname "/tmp/output/$OUTPUT_FILE")
        if [[ "$OUTPUT_DIR" != "/tmp/output" ]]; then
            mkdir -p "$OUTPUT_DIR"
            echo "📁 Created output directory: $OUTPUT_DIR"
        fi

        echo "🚀 Executing: $BEDTOOLS_CMD > /tmp/output/$OUTPUT_FILE"

        # Execute the command
        if eval "$BEDTOOLS_CMD > /tmp/output/$OUTPUT_FILE"; then
            echo "✅ BEDTools coverage completed successfully"
            
            # Display output file info
            OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILE" 2>/dev/null || echo "unknown")
            echo "📊 Output file: $OUTPUT_FILE ($OUTPUT_SIZE bytes)"
        else
            echo "❌ BEDTools coverage failed"
            exit 1
        fi
        ;;
        
    "subtract")
        echo "➖ Running BEDTools subtract command..."
        
        # Copy input files
        cp "$INPUT_FILE_PATH" "/tmp/input_intervals.bed"
        
        BEDTOOLS_CMD="bedtools subtract"
        
        # Add input file
        BEDTOOLS_CMD="$BEDTOOLS_CMD -a /tmp/input_intervals.bed"
        
        # Add second file (for subtraction)
        # Note: In real implementation, we'd handle multiple input files
        BEDTOOLS_CMD="$BEDTOOLS_CMD -b /tmp/input_intervals.bed"
        
        if [[ "$STRAND_SPECIFIC" == "true" ]]; then
            BEDTOOLS_CMD="$BEDTOOLS_CMD -s"
        fi
        
        # Create subdirectories in output path if needed
        OUTPUT_DIR=$(dirname "/tmp/output/$OUTPUT_FILE")
        if [[ "$OUTPUT_DIR" != "/tmp/output" ]]; then
            mkdir -p "$OUTPUT_DIR"
            echo "📁 Created output directory: $OUTPUT_DIR"
        fi

        echo "🚀 Executing: $BEDTOOLS_CMD > /tmp/output/$OUTPUT_FILE"

        # Execute the command
        if eval "$BEDTOOLS_CMD > /tmp/output/$OUTPUT_FILE"; then
            echo "✅ BEDTools subtract completed successfully"
            
            # Display output file info
            OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILE" 2>/dev/null || echo "unknown")
            echo "📊 Output file: $OUTPUT_FILE ($OUTPUT_SIZE bytes)"
        else
            echo "❌ BEDTools subtract failed"
            exit 1
        fi
        ;;
        
    "closest")
        echo "📍 Running BEDTools closest command..."
        
        # Copy input files
        cp "$INPUT_FILE_PATH" "/tmp/input_intervals.bed"
        
        BEDTOOLS_CMD="bedtools closest"
        
        # Add input file
        BEDTOOLS_CMD="$BEDTOOLS_CMD -a /tmp/input_intervals.bed"
        
        # Add second file (for finding closest)
        BEDTOOLS_CMD="$BEDTOOLS_CMD -b /tmp/input_intervals.bed"
        
        if [[ "$STRAND_SPECIFIC" == "true" ]]; then
            BEDTOOLS_CMD="$BEDTOOLS_CMD -s"
        fi
        
        # Create subdirectories in output path if needed
        OUTPUT_DIR=$(dirname "/tmp/output/$OUTPUT_FILE")
        if [[ "$OUTPUT_DIR" != "/tmp/output" ]]; then
            mkdir -p "$OUTPUT_DIR"
            echo "📁 Created output directory: $OUTPUT_DIR"
        fi

        echo "🚀 Executing: $BEDTOOLS_CMD > /tmp/output/$OUTPUT_FILE"

        # Execute the command
        if eval "$BEDTOOLS_CMD > /tmp/output/$OUTPUT_FILE"; then
            echo "✅ BEDTools closest completed successfully"
            
            # Display output file info
            OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILE" 2>/dev/null || echo "unknown")
            echo "📊 Output file: $OUTPUT_FILE ($OUTPUT_SIZE bytes)"
        else
            echo "❌ BEDTools closest failed"
            exit 1
        fi
        ;;
        
    "complement")
        echo "🔄 Running BEDTools complement command..."
        
        # Copy input file
        cp "$INPUT_FILE_PATH" "/tmp/input_intervals.bed"
        
        # Complement requires a genome file
        BEDTOOLS_CMD="bedtools complement -i /tmp/input_intervals.bed"
        
        # Add genome parameter (simplified)
        if [[ -n "$GENOME_FILE" ]]; then
            BEDTOOLS_CMD="$BEDTOOLS_CMD -g $GENOME_FILE"
        else
            # Create a simple genome file for demonstration
            echo -e "chr1\t250000000\nchr2\t250000000\nchr3\t200000000" > /tmp/genome.txt
            BEDTOOLS_CMD="$BEDTOOLS_CMD -g /tmp/genome.txt"
        fi
        
        # Create subdirectories in output path if needed
        OUTPUT_DIR=$(dirname "/tmp/output/$OUTPUT_FILE")
        if [[ "$OUTPUT_DIR" != "/tmp/output" ]]; then
            mkdir -p "$OUTPUT_DIR"
            echo "📁 Created output directory: $OUTPUT_DIR"
        fi

        echo "🚀 Executing: $BEDTOOLS_CMD > /tmp/output/$OUTPUT_FILE"

        # Execute the command
        if eval "$BEDTOOLS_CMD > /tmp/output/$OUTPUT_FILE"; then
            echo "✅ BEDTools complement completed successfully"
            
            # Display output file info
            OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILE" 2>/dev/null || echo "unknown")
            echo "📊 Output file: $OUTPUT_FILE ($OUTPUT_SIZE bytes)"
        else
            echo "❌ BEDTools complement failed"
            exit 1
        fi
        ;;
        
    "flank")
        echo "↔️ Running BEDTools flank command..."
        
        # Copy input file
        cp "$INPUT_FILE_PATH" "/tmp/input_intervals.bed"
        
        FLANK_SIZE=${FLANK_SIZE:-"1000"}
        
        BEDTOOLS_CMD="bedtools flank -i /tmp/input_intervals.bed -b $FLANK_SIZE"
        
        # Add genome parameter
        if [[ -n "$GENOME_FILE" ]]; then
            BEDTOOLS_CMD="$BEDTOOLS_CMD -g $GENOME_FILE"
        else
            # Create a simple genome file for demonstration
            echo -e "chr1\t250000000\nchr2\t250000000\nchr3\t200000000" > /tmp/genome.txt
            BEDTOOLS_CMD="$BEDTOOLS_CMD -g /tmp/genome.txt"
        fi
        
        if [[ "$STRAND_SPECIFIC" == "true" ]]; then
            BEDTOOLS_CMD="$BEDTOOLS_CMD -s"
        fi
        
        # Create subdirectories in output path if needed
        OUTPUT_DIR=$(dirname "/tmp/output/$OUTPUT_FILE")
        if [[ "$OUTPUT_DIR" != "/tmp/output" ]]; then
            mkdir -p "$OUTPUT_DIR"
            echo "📁 Created output directory: $OUTPUT_DIR"
        fi

        echo "🚀 Executing: $BEDTOOLS_CMD > /tmp/output/$OUTPUT_FILE"

        # Execute the command
        if eval "$BEDTOOLS_CMD > /tmp/output/$OUTPUT_FILE"; then
            echo "✅ BEDTools flank completed successfully"
            
            # Display output file info
            OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILE" 2>/dev/null || echo "unknown")
            echo "📊 Output file: $OUTPUT_FILE ($OUTPUT_SIZE bytes)"
        else
            echo "❌ BEDTools flank failed"
            exit 1
        fi
        ;;
        
    "slop")
        echo "📏 Running BEDTools slop command..."
        
        # Copy input file
        cp "$INPUT_FILE_PATH" "/tmp/input_intervals.bed"
        
        SLOP_SIZE=${SLOP_SIZE:-"500"}
        
        BEDTOOLS_CMD="bedtools slop -i /tmp/input_intervals.bed -b $SLOP_SIZE"
        
        # Add genome parameter
        if [[ -n "$GENOME_FILE" ]]; then
            BEDTOOLS_CMD="$BEDTOOLS_CMD -g $GENOME_FILE"
        else
            # Create a simple genome file for demonstration
            echo -e "chr1\t250000000\nchr2\t250000000\nchr3\t200000000" > /tmp/genome.txt
            BEDTOOLS_CMD="$BEDTOOLS_CMD -g /tmp/genome.txt"
        fi
        
        if [[ "$STRAND_SPECIFIC" == "true" ]]; then
            BEDTOOLS_CMD="$BEDTOOLS_CMD -s"
        fi
        
        # Create subdirectories in output path if needed
        OUTPUT_DIR=$(dirname "/tmp/output/$OUTPUT_FILE")
        if [[ "$OUTPUT_DIR" != "/tmp/output" ]]; then
            mkdir -p "$OUTPUT_DIR"
            echo "📁 Created output directory: $OUTPUT_DIR"
        fi

        echo "🚀 Executing: $BEDTOOLS_CMD > /tmp/output/$OUTPUT_FILE"

        # Execute the command
        if eval "$BEDTOOLS_CMD > /tmp/output/$OUTPUT_FILE"; then
            echo "✅ BEDTools slop completed successfully"
            
            # Display output file info
            OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILE" 2>/dev/null || echo "unknown")
            echo "📊 Output file: $OUTPUT_FILE ($OUTPUT_SIZE bytes)"
        else
            echo "❌ BEDTools slop failed"
            exit 1
        fi
        ;;
        
    "window")
        echo "🪟 Running BEDTools window command..."
        
        # Copy input files
        cp "$INPUT_FILE_PATH" "/tmp/input_intervals.bed"
        
        WINDOW_SIZE=${WINDOW_SIZE:-"1000"}
        
        BEDTOOLS_CMD="bedtools window"
        
        # Add input file
        BEDTOOLS_CMD="$BEDTOOLS_CMD -a /tmp/input_intervals.bed"
        
        # Add second file (for window search)
        BEDTOOLS_CMD="$BEDTOOLS_CMD -b /tmp/input_intervals.bed"
        
        # Add window size
        BEDTOOLS_CMD="$BEDTOOLS_CMD -w $WINDOW_SIZE"
        
        if [[ "$STRAND_SPECIFIC" == "true" ]]; then
            BEDTOOLS_CMD="$BEDTOOLS_CMD -s"
        fi
        
        # Create subdirectories in output path if needed
        OUTPUT_DIR=$(dirname "/tmp/output/$OUTPUT_FILE")
        if [[ "$OUTPUT_DIR" != "/tmp/output" ]]; then
            mkdir -p "$OUTPUT_DIR"
            echo "📁 Created output directory: $OUTPUT_DIR"
        fi

        echo "🚀 Executing: $BEDTOOLS_CMD > /tmp/output/$OUTPUT_FILE"

        # Execute the command
        if eval "$BEDTOOLS_CMD > /tmp/output/$OUTPUT_FILE"; then
            echo "✅ BEDTools window completed successfully"
            
            # Display output file info
            OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILE" 2>/dev/null || echo "unknown")
            echo "📊 Output file: $OUTPUT_FILE ($OUTPUT_SIZE bytes)"
        else
            echo "❌ BEDTools window failed"
            exit 1
        fi
        ;;
        
    "cluster")
        echo "🔗 Running BEDTools cluster command..."
        
        # Copy input file
        cp "$INPUT_FILE_PATH" "/tmp/input_intervals.bed"
        
        MAX_DISTANCE=${MAX_DISTANCE:-"0"}
        
        BEDTOOLS_CMD="bedtools cluster -i /tmp/input_intervals.bed"
        
        if [[ "$MAX_DISTANCE" != "0" ]]; then
            BEDTOOLS_CMD="$BEDTOOLS_CMD -d $MAX_DISTANCE"
        fi
        
        if [[ "$STRAND_SPECIFIC" == "true" ]]; then
            BEDTOOLS_CMD="$BEDTOOLS_CMD -s"
        fi
        
        # Create subdirectories in output path if needed
        OUTPUT_DIR=$(dirname "/tmp/output/$OUTPUT_FILE")
        if [[ "$OUTPUT_DIR" != "/tmp/output" ]]; then
            mkdir -p "$OUTPUT_DIR"
            echo "📁 Created output directory: $OUTPUT_DIR"
        fi

        echo "🚀 Executing: $BEDTOOLS_CMD > /tmp/output/$OUTPUT_FILE"

        # Execute the command
        if eval "$BEDTOOLS_CMD > /tmp/output/$OUTPUT_FILE"; then
            echo "✅ BEDTools cluster completed successfully"
            
            # Display output file info
            OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILE" 2>/dev/null || echo "unknown")
            echo "📊 Output file: $OUTPUT_FILE ($OUTPUT_SIZE bytes)"
        else
            echo "❌ BEDTools cluster failed"
            exit 1
        fi
        ;;
        
    *)
        echo "❌ Unsupported BEDTools command: $COMMAND"
        echo "Supported commands: intersect, merge, sort, coverage, subtract, closest, complement, flank, slop, window, cluster"
        exit 1
        ;;
esac

echo "🎯 BEDTools handler completed successfully"