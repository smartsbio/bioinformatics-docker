#!/bin/bash
set -e

# GenomicRanges Handler Script
# Handles R-based genomic analysis using Bioconductor GenomicRanges

echo "🧬 GenomicRanges Handler Started"

# Get the input file (should be downloaded already)
INPUT_FILES=(/tmp/input/*)
if [[ ${#INPUT_FILES[@]} -eq 0 ]]; then
    echo "❌ No input files found in /tmp/input/"
    exit 1
fi

INPUT_FILE_PATH="${INPUT_FILES[0]}"
INPUT_FILENAME=$(basename "$INPUT_FILE_PATH")

echo "📁 Processing file: $INPUT_FILENAME"
echo "🎯 GenomicRanges command: $COMMAND"

# Parse additional parameters from environment
OUTPUT_FILE=${OUTPUT_FILE:-"genomic_analysis.csv"}
GENOME_BUILD=${GENOME_BUILD:-"hg38"}
CHR_PREFIX=${CHR_PREFIX:-"chr"}
OUTPUT_FORMAT=${OUTPUT_FORMAT:-"csv"}
WINDOW_SIZE=${WINDOW_SIZE:-"1000"}
OVERLAP_TYPE=${OVERLAP_TYPE:-"any"}

case "$COMMAND" in
    "interval-operations")
        echo "🔍 Running GenomicRanges interval operations..."
        
        # Copy input file (should be BED or GRanges format)
        cp "$INPUT_FILE_PATH" "/tmp/intervals.bed"
        
        # Create R script for interval operations
        cat > "/tmp/genomicranges_intervals.R" << 'EOF'
library(GenomicRanges)
library(rtracklayer)
library(IRanges)

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]
genome_build <- args[3]
operation <- ifelse(length(args) >= 4, args[4], "merge")

# Read intervals from BED file
cat("📊 Reading intervals from", input_file, "\n")
intervals <- import(input_file, format = "BED")

# Perform interval operations
cat("🔧 Performing", operation, "operation...\n")
result <- switch(operation,
    "merge" = reduce(intervals),
    "gaps" = gaps(intervals),
    "flank" = flank(intervals, width = 1000),
    "resize" = resize(intervals, width = 2000),
    "shift" = shift(intervals, shift = 500),
    intervals  # default: return original
)

# Add metadata
mcols(result)$operation <- operation
mcols(result)$genome_build <- genome_build

# Export results
cat("💾 Writing results to", output_file, "\n")
export(result, output_file, format = "BED")

# Summary statistics
cat("✅ Interval operations completed\n")
cat("📊 Original intervals:", length(intervals), "\n")
cat("📊 Result intervals:", length(result), "\n")
cat("📊 Total width original:", sum(width(intervals)), "bp\n")
cat("📊 Total width result:", sum(width(result)), "bp\n")
EOF

        echo "🚀 Executing R script for interval operations..."
        
        if Rscript "/tmp/genomicranges_intervals.R" "/tmp/intervals.bed" "/tmp/output/$OUTPUT_FILE" "$GENOME_BUILD" "merge"; then
            echo "✅ GenomicRanges interval operations completed successfully"
            
            # Display output file info
            OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILE" 2>/dev/null || echo "unknown")
            echo "📊 Genomic intervals: $OUTPUT_FILE ($OUTPUT_SIZE bytes)"
        else
            echo "❌ GenomicRanges interval operations failed"
            exit 1
        fi
        ;;
        
    "overlap-analysis")
        echo "🔍 Running GenomicRanges overlap analysis..."
        
        # For overlap analysis, we need two input files
        # This is simplified - in practice, you'd have two separate input files
        cp "$INPUT_FILE_PATH" "/tmp/query_intervals.bed"
        cp "$INPUT_FILE_PATH" "/tmp/subject_intervals.bed"  # Placeholder
        
        # Create R script for overlap analysis
        cat > "/tmp/genomicranges_overlap.R" << 'EOF'
library(GenomicRanges)
library(rtracklayer)
library(IRanges)

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
query_file <- args[1]
subject_file <- args[2]
output_file <- args[3]
overlap_type <- ifelse(length(args) >= 4, args[4], "any")

# Read intervals
cat("📊 Reading query intervals from", query_file, "\n")
query <- import(query_file, format = "BED")

cat("📊 Reading subject intervals from", subject_file, "\n")
subject <- import(subject_file, format = "BED")

# Perform overlap analysis
cat("🔧 Performing overlap analysis (type:", overlap_type, ")...\n")

# Find overlaps
overlaps <- findOverlaps(query, subject, type = overlap_type)
query_hits <- queryHits(overlaps)
subject_hits <- subjectHits(overlaps)

# Create result data frame
overlap_results <- data.frame(
    query_chr = as.character(seqnames(query[query_hits])),
    query_start = start(query[query_hits]),
    query_end = end(query[query_hits]),
    subject_chr = as.character(seqnames(subject[subject_hits])),
    subject_start = start(subject[subject_hits]),
    subject_end = end(subject[subject_hits]),
    overlap_width = width(pintersect(query[query_hits], subject[subject_hits])),
    stringsAsFactors = FALSE
)

# Write results
cat("💾 Writing overlap results to", output_file, "\n")
write.csv(overlap_results, output_file, row.names = FALSE)

# Summary statistics
cat("✅ Overlap analysis completed\n")
cat("📊 Query intervals:", length(query), "\n")
cat("📊 Subject intervals:", length(subject), "\n")
cat("📊 Overlaps found:", length(overlaps), "\n")
cat("📊 Query intervals with overlaps:", length(unique(query_hits)), "\n")
EOF

        echo "🚀 Executing R script for overlap analysis..."
        
        if Rscript "/tmp/genomicranges_overlap.R" "/tmp/query_intervals.bed" "/tmp/subject_intervals.bed" "/tmp/output/$OUTPUT_FILE" "$OVERLAP_TYPE"; then
            echo "✅ GenomicRanges overlap analysis completed successfully"
            
            OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILE" 2>/dev/null || echo "unknown")
            echo "📊 Overlap analysis: $OUTPUT_FILE ($OUTPUT_SIZE bytes)"
        else
            echo "❌ GenomicRanges overlap analysis failed"
            exit 1
        fi
        ;;
        
    "windowed-analysis")
        echo "📊 Running GenomicRanges windowed analysis..."
        
        cp "$INPUT_FILE_PATH" "/tmp/genomic_data.bed"
        
        # Create R script for windowed analysis
        cat > "/tmp/genomicranges_windows.R" << 'EOF'
library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
library(BSgenome)

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]
window_size <- as.numeric(args[3])
genome_build <- args[4]

# Read genomic data
cat("📊 Reading genomic data from", input_file, "\n")
genomic_data <- import(input_file, format = "BED")

# Create genome windows
cat("🔧 Creating", window_size, "bp windows for", genome_build, "...\n")

# Simplified approach - create windows for main chromosomes
chromosomes <- paste0("chr", c(1:22, "X", "Y"))
chr_lengths <- c(rep(200000000, 24))  # Simplified chromosome lengths

# Create seqinfo
seqinfo_obj <- Seqinfo(seqnames = chromosomes, seqlengths = chr_lengths, genome = genome_build)

# Create tiling windows
windows <- tileGenome(seqinfo_obj, tilewidth = window_size, cut.last.tile.in.chrom = TRUE)

# Count overlaps between data and windows
cat("🔧 Counting overlaps in windows...\n")
overlap_counts <- countOverlaps(windows, genomic_data)

# Create result data frame
window_results <- data.frame(
    chr = as.character(seqnames(windows)),
    start = start(windows),
    end = end(windows),
    width = width(windows),
    count = overlap_counts,
    density = overlap_counts / width(windows) * 1000,  # per kb
    stringsAsFactors = FALSE
)

# Write results
cat("💾 Writing windowed analysis to", output_file, "\n")
write.csv(window_results, output_file, row.names = FALSE)

# Summary statistics
cat("✅ Windowed analysis completed\n")
cat("📊 Total windows:", nrow(window_results), "\n")
cat("📊 Windows with data:", sum(window_results$count > 0), "\n")
cat("📊 Mean density:", mean(window_results$density), "per kb\n")
EOF

        echo "🚀 Executing R script for windowed analysis..."
        
        if Rscript "/tmp/genomicranges_windows.R" "/tmp/genomic_data.bed" "/tmp/output/$OUTPUT_FILE" "$WINDOW_SIZE" "$GENOME_BUILD"; then
            echo "✅ GenomicRanges windowed analysis completed successfully"
            
            OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILE" 2>/dev/null || echo "unknown")
            echo "📊 Windowed analysis: $OUTPUT_FILE ($OUTPUT_SIZE bytes)"
        else
            echo "❌ GenomicRanges windowed analysis failed"
            exit 1
        fi
        ;;
        
    "annotation")
        echo "📝 Running GenomicRanges annotation..."
        
        cp "$INPUT_FILE_PATH" "/tmp/regions.bed"
        
        # Create R script for annotation
        cat > "/tmp/genomicranges_annotation.R" << 'EOF'
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]
genome_build <- args[3]

# Read regions to annotate
cat("📊 Reading regions from", input_file, "\n")
regions <- import(input_file, format = "BED")

# Load transcript database
cat("🔧 Loading transcript annotations for", genome_build, "...\n")
if (genome_build == "hg38") {
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
} else {
    cat("⚠️ Using hg38 annotations as default\n")
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
}

# Get gene annotations
genes <- genes(txdb)
exons <- exons(txdb)
promoters_regions <- promoters(genes, upstream = 2000, downstream = 500)

# Annotate regions
cat("🔧 Annotating regions...\n")

# Find overlaps with different features
gene_overlaps <- countOverlaps(regions, genes)
exon_overlaps <- countOverlaps(regions, exons)
promoter_overlaps <- countOverlaps(regions, promoters_regions)

# Create annotation results
annotation_results <- data.frame(
    chr = as.character(seqnames(regions)),
    start = start(regions),
    end = end(regions),
    width = width(regions),
    gene_overlaps = gene_overlaps,
    exon_overlaps = exon_overlaps,
    promoter_overlaps = promoter_overlaps,
    annotation = ifelse(promoter_overlaps > 0, "promoter",
                       ifelse(exon_overlaps > 0, "exonic",
                             ifelse(gene_overlaps > 0, "intronic", "intergenic"))),
    stringsAsFactors = FALSE
)

# Write results
cat("💾 Writing annotation results to", output_file, "\n")
write.csv(annotation_results, output_file, row.names = FALSE)

# Summary statistics
cat("✅ Annotation completed\n")
cat("📊 Total regions:", nrow(annotation_results), "\n")
cat("📊 Promoter regions:", sum(annotation_results$annotation == "promoter"), "\n")
cat("📊 Exonic regions:", sum(annotation_results$annotation == "exonic"), "\n")
cat("📊 Intronic regions:", sum(annotation_results$annotation == "intronic"), "\n")
cat("📊 Intergenic regions:", sum(annotation_results$annotation == "intergenic"), "\n")
EOF

        echo "🚀 Executing R script for annotation..."
        
        if Rscript "/tmp/genomicranges_annotation.R" "/tmp/regions.bed" "/tmp/output/$OUTPUT_FILE" "$GENOME_BUILD"; then
            echo "✅ GenomicRanges annotation completed successfully"
            
            OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILE" 2>/dev/null || echo "unknown")
            echo "📊 Annotation results: $OUTPUT_FILE ($OUTPUT_SIZE bytes)"
        else
            echo "❌ GenomicRanges annotation failed"
            exit 1
        fi
        ;;
        
    *)
        echo "❌ Unsupported GenomicRanges command: $COMMAND"
        echo "Supported commands: interval-operations, overlap-analysis, windowed-analysis, annotation"
        exit 1
        ;;
esac

echo "🎯 GenomicRanges handler completed successfully"