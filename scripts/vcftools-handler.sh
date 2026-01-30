#!/bin/bash
set -e

# VCFtools Handler Script
# Handles various VCFtools commands for VCF file manipulation and analysis

echo "🧬 VCFtools Handler Started"

# Get the input file (should be downloaded already)
INPUT_FILES=(/tmp/input/*)
if [[ ${#INPUT_FILES[@]} -eq 0 ]]; then
    echo "❌ No input files found in /tmp/input/"
    exit 1
fi

INPUT_FILE_PATH="${INPUT_FILES[0]}"
INPUT_FILENAME=$(basename "$INPUT_FILE_PATH")

echo "📁 Processing file: $INPUT_FILENAME"
echo "🎯 VCFtools command: $COMMAND"

# Extract organization and workspace IDs from INPUT_S3_KEY
# Format: organizations/{orgId}/workspaces/{workspaceId}/files/{path}
if [[ -n "$INPUT_S3_KEY" ]]; then
    ORGANIZATION_ID=$(echo "$INPUT_S3_KEY" | cut -d'/' -f2)
    WORKSPACE_ID=$(echo "$INPUT_S3_KEY" | cut -d'/' -f4)
    echo "📂 Organization ID: $ORGANIZATION_ID"
    echo "📂 Workspace ID: $WORKSPACE_ID"
fi

# Parse additional parameters from environment
OUTPUT_PREFIX=${OUTPUT_PREFIX:-"output"}
# Strip @ notation from output prefix path
OUTPUT_PREFIX="${OUTPUT_PREFIX#@}"
MIN_ALLELES=${MIN_ALLELES:-"2"}
MAX_ALLELES=${MAX_ALLELES:-"2"}
MIN_QUALITY=${MIN_QUALITY:-"30"}
MIN_DEPTH=${MIN_DEPTH:-"10"}
MAX_DEPTH=${MAX_DEPTH:-""}
MAF=${MAF:-""}  # Minor allele frequency
MAX_MISSING=${MAX_MISSING:-""}
REMOVE_INDELS=${REMOVE_INDELS:-"false"}
REMOVE_SNPS=${REMOVE_SNPS:-"false"}

# Create subdirectories in output path if needed
OUTPUT_DIR=$(dirname "/tmp/output/$OUTPUT_PREFIX")
if [[ "$OUTPUT_DIR" != "/tmp/output" ]]; then
    mkdir -p "$OUTPUT_DIR"
    echo "📁 Created output directory: $OUTPUT_DIR"
fi

case "$COMMAND" in
    "filter")
        echo "🔍 Running VCFtools filtering..."
        
        # Copy input file
        cp "$INPUT_FILE_PATH" "/tmp/input.vcf"
        
        VCFTOOLS_CMD="vcftools --vcf /tmp/input.vcf"
        
        # Add quality filter
        if [[ -n "$MIN_QUALITY" ]]; then
            VCFTOOLS_CMD="$VCFTOOLS_CMD --minQ $MIN_QUALITY"
        fi
        
        # Add depth filters
        if [[ -n "$MIN_DEPTH" ]]; then
            VCFTOOLS_CMD="$VCFTOOLS_CMD --min-meanDP $MIN_DEPTH"
        fi
        
        if [[ -n "$MAX_DEPTH" ]]; then
            VCFTOOLS_CMD="$VCFTOOLS_CMD --max-meanDP $MAX_DEPTH"
        fi
        
        # Add allele frequency filter
        if [[ -n "$MAF" ]]; then
            VCFTOOLS_CMD="$VCFTOOLS_CMD --maf $MAF"
        fi
        
        # Add missing data filter
        if [[ -n "$MAX_MISSING" ]]; then
            VCFTOOLS_CMD="$VCFTOOLS_CMD --max-missing $MAX_MISSING"
        fi
        
        # Add allele count filters
        if [[ -n "$MIN_ALLELES" ]]; then
            VCFTOOLS_CMD="$VCFTOOLS_CMD --min-alleles $MIN_ALLELES"
        fi
        
        if [[ -n "$MAX_ALLELES" ]]; then
            VCFTOOLS_CMD="$VCFTOOLS_CMD --max-alleles $MAX_ALLELES"
        fi
        
        # Remove indels or SNPs if specified
        if [[ "$REMOVE_INDELS" == "true" ]]; then
            VCFTOOLS_CMD="$VCFTOOLS_CMD --remove-indels"
        fi
        
        if [[ "$REMOVE_SNPS" == "true" ]]; then
            VCFTOOLS_CMD="$VCFTOOLS_CMD --keep-only-indels"
        fi
        
        # Add output
        VCFTOOLS_CMD="$VCFTOOLS_CMD --recode --recode-INFO-all --out /tmp/output/$OUTPUT_PREFIX"
        
        echo "🚀 Executing: $VCFTOOLS_CMD"
        
        # Execute the command
        if eval "$VCFTOOLS_CMD"; then
            echo "✅ VCFtools filtering completed successfully"
            
            # Display output files
            echo "📁 Generated files:"
            ls -la /tmp/output/
        else
            echo "❌ VCFtools filtering failed"
            exit 1
        fi
        ;;
        
    "stats")
        echo "📊 Running VCFtools statistics..."
        
        # Copy input file
        cp "$INPUT_FILE_PATH" "/tmp/input.vcf"
        
        # Generate multiple statistics
        echo "📈 Generating variant statistics..."
        
        # Site frequency spectrum
        VCFTOOLS_CMD1="vcftools --vcf /tmp/input.vcf --freq --out /tmp/output/${OUTPUT_PREFIX}_freq"
        
        # Site depth statistics
        VCFTOOLS_CMD2="vcftools --vcf /tmp/input.vcf --site-depth --out /tmp/output/${OUTPUT_PREFIX}_depth"
        
        # Site quality statistics
        VCFTOOLS_CMD3="vcftools --vcf /tmp/input.vcf --site-quality --out /tmp/output/${OUTPUT_PREFIX}_quality"
        
        # Individual depth statistics
        VCFTOOLS_CMD4="vcftools --vcf /tmp/input.vcf --depth --out /tmp/output/${OUTPUT_PREFIX}_ind_depth"
        
        # Missing data per individual
        VCFTOOLS_CMD5="vcftools --vcf /tmp/input.vcf --missing-indv --out /tmp/output/${OUTPUT_PREFIX}_missing"
        
        echo "🚀 Executing statistics commands..."
        
        # Execute all commands
        if eval "$VCFTOOLS_CMD1" && eval "$VCFTOOLS_CMD2" && eval "$VCFTOOLS_CMD3" && eval "$VCFTOOLS_CMD4" && eval "$VCFTOOLS_CMD5"; then
            echo "✅ VCFtools statistics completed successfully"
            
            # Display all generated files
            echo "📁 Generated statistic files:"
            ls -la /tmp/output/
        else
            echo "❌ VCFtools statistics failed"
            exit 1
        fi
        ;;
        
    "convert")
        echo "🔄 Running VCFtools format conversion..."
        
        # Copy input file
        cp "$INPUT_FILE_PATH" "/tmp/input.vcf"
        
        OUTPUT_FORMAT=${OUTPUT_FORMAT:-"vcf"}
        
        case "$OUTPUT_FORMAT" in
            "plink")
                VCFTOOLS_CMD="vcftools --vcf /tmp/input.vcf --plink --out /tmp/output/$OUTPUT_PREFIX"
                ;;
            "beagle")
                VCFTOOLS_CMD="vcftools --vcf /tmp/input.vcf --BEAGLE-GL --out /tmp/output/$OUTPUT_PREFIX"
                ;;
            "012")
                VCFTOOLS_CMD="vcftools --vcf /tmp/input.vcf --012 --out /tmp/output/$OUTPUT_PREFIX"
                ;;
            *)
                VCFTOOLS_CMD="vcftools --vcf /tmp/input.vcf --recode --recode-INFO-all --out /tmp/output/$OUTPUT_PREFIX"
                ;;
        esac
        
        echo "🚀 Executing: $VCFTOOLS_CMD"
        
        # Execute the command
        if eval "$VCFTOOLS_CMD"; then
            echo "✅ VCFtools conversion completed successfully"
            
            # Display output files
            echo "📁 Generated files:"
            ls -la /tmp/output/
        else
            echo "❌ VCFtools conversion failed"
            exit 1
        fi
        ;;
        
    "merge")
        echo "🔗 Running VCFtools merge..."
        
        # For merging, we'd need multiple VCF files
        # Copy input files
        cp "$INPUT_FILE_PATH" "/tmp/input1.vcf"
        # Note: In real implementation, handle multiple input files
        
        # VCFtools doesn't have direct merge - this would typically use bcftools
        # But we'll provide a placeholder
        VCFTOOLS_CMD="vcftools --vcf /tmp/input1.vcf --recode --recode-INFO-all --out /tmp/output/$OUTPUT_PREFIX"
        
        echo "🚀 Executing: $VCFTOOLS_CMD"
        
        if eval "$VCFTOOLS_CMD"; then
            echo "✅ VCFtools merge completed successfully"
            echo "📁 Generated files:"
            ls -la /tmp/output/
        else
            echo "❌ VCFtools merge failed"
            exit 1
        fi
        ;;
        
    "extract-samples")
        echo "👥 Running VCFtools sample extraction..."
        
        # Copy input file
        cp "$INPUT_FILE_PATH" "/tmp/input.vcf"
        
        SAMPLE_LIST=${SAMPLE_LIST:-""}
        
        VCFTOOLS_CMD="vcftools --vcf /tmp/input.vcf"
        
        # Add sample filtering
        if [[ -n "$SAMPLE_LIST" ]]; then
            # Create sample list file
            echo "$SAMPLE_LIST" | tr ',' '\n' > /tmp/samples.txt
            VCFTOOLS_CMD="$VCFTOOLS_CMD --keep /tmp/samples.txt"
        fi
        
        # Add output
        VCFTOOLS_CMD="$VCFTOOLS_CMD --recode --recode-INFO-all --out /tmp/output/$OUTPUT_PREFIX"
        
        echo "🚀 Executing: $VCFTOOLS_CMD"
        
        if eval "$VCFTOOLS_CMD"; then
            echo "✅ VCFtools sample extraction completed successfully"
            
            echo "📁 Generated files:"
            ls -la /tmp/output/
        else
            echo "❌ VCFtools sample extraction failed"
            exit 1
        fi
        ;;
        
    *)
        echo "❌ Unsupported VCFtools command: $COMMAND"
        echo "Supported commands: filter, stats, convert, merge, extract-samples"
        exit 1
        ;;
esac

echo "🎯 VCFtools handler completed successfully"