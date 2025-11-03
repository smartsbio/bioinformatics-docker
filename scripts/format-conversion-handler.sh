#!/bin/bash
# Format Conversion Handler
# Handles all bioinformatics file format conversions on AWS ECS

set -e  # Exit on error
set -o pipefail

echo "========================================="
echo "🔄 Format Conversion Handler Started"
echo "========================================="

# ==========================================
# Environment Variables
# ==========================================
CONVERSION_COMMAND="${CONVERSION_COMMAND:-}"
SOURCE_FORMAT="${SOURCE_FORMAT:-auto}"
TARGET_FORMAT="${TARGET_FORMAT:-}"
COMPRESSION_LEVEL="${COMPRESSION_LEVEL:-6}"
DEFAULT_QUALITY="${DEFAULT_QUALITY:-30}"
PRESERVE_METADATA="${PRESERVE_METADATA:-true}"
EXTRACT_ARCHIVE="${EXTRACT_ARCHIVE:-true}"
PRESERVE_STRUCTURE="${PRESERVE_STRUCTURE:-false}"
REFERENCE_GENOME="${REFERENCE_GENOME:-}"
IS_ARCHIVE="${IS_ARCHIVE:-false}"
THREADS="${THREADS:-4}"

echo "📋 Configuration:"
echo "   Conversion: $CONVERSION_COMMAND"
echo "   Source Format: $SOURCE_FORMAT"
echo "   Target Format: $TARGET_FORMAT"
echo "   Input: /tmp/input/"
echo "   Output: /tmp/output/"
echo ""

# ==========================================
# SEQUENCE CONVERSION FUNCTIONS
# ==========================================

convert_fasta_to_fastq() {
    local input_file="$1"
    local output_file="$2"
    local quality_score="${DEFAULT_QUALITY:-30}"

    # Convert quality score to ASCII character (Phred+33 encoding)
    local quality_char=$(printf "\\$(printf '%03o' $((quality_score + 33)))")

    echo "Converting FASTA to FASTQ: $(basename "$input_file") → $(basename "$output_file")"
    echo "Using quality score: $quality_score (ASCII: $quality_char)"

    if command -v seqtk &>/dev/null; then
        seqtk seq -F "$quality_char" "$input_file" > "$output_file"
    else
        awk -v q="$quality_char" '
            /^>/ {
                if (seq) {
                    print "@" substr(header, 2);
                    print seq;
                    print "+";
                    for (i=1; i<=length(seq); i++) printf q;
                    printf "\n";
                }
                header = $0;
                seq = "";
                next;
            }
            { seq = seq $0 }
            END {
                if (seq) {
                    print "@" substr(header, 2);
                    print seq;
                    print "+";
                    for (i=1; i<=length(seq); i++) printf q;
                    printf "\n";
                }
            }
        ' "$input_file" > "$output_file"
    fi

    echo "✅ Conversion complete"
}

convert_fastq_to_fasta() {
    local input_file="$1"
    local output_file="$2"

    echo "Converting FASTQ to FASTA: $(basename "$input_file") → $(basename "$output_file")"

    if command -v seqtk &>/dev/null; then
        seqtk seq -A "$input_file" > "$output_file"
    else
        awk 'NR%4==1 { print ">" substr($0, 2) } NR%4==2 { print }' "$input_file" > "$output_file"
    fi

    echo "✅ Conversion complete"
}

convert_genbank_to_fasta() {
    local input_file="$1"
    local output_file="$2"

    echo "Converting GenBank to FASTA: $(basename "$input_file") → $(basename "$output_file")"

    python3 <<EOF
from Bio import SeqIO
try:
    count = SeqIO.convert("$input_file", "genbank", "$output_file", "fasta")
    print(f"✅ Converted {count} sequence(s)")
except Exception as e:
    print(f"❌ Error: {e}")
    exit(1)
EOF
}

convert_embl_to_fasta() {
    local input_file="$1"
    local output_file="$2"

    echo "Converting EMBL to FASTA: $(basename "$input_file") → $(basename "$output_file")"

    python3 <<EOF
from Bio import SeqIO
try:
    count = SeqIO.convert("$input_file", "embl", "$output_file", "fasta")
    print(f"✅ Converted {count} sequence(s)")
except Exception as e:
    print(f"❌ Error: {e}")
    exit(1)
EOF
}

convert_ab1_to_fasta() {
    local input_file="$1"
    local output_file="$2"

    echo "Converting AB1 to FASTA: $(basename "$input_file") → $(basename "$output_file")"

    python3 <<EOF
from Bio import SeqIO
import os
try:
    record = SeqIO.read("$input_file", "abi")
    if not record.id or record.id == "<unknown id>":
        record.id = os.path.splitext(os.path.basename("$input_file"))[0]
        record.description = "Sanger sequencing read"
    SeqIO.write(record, "$output_file", "fasta")
    print(f"✅ Converted sequence: {record.id} ({len(record.seq)} bp)")
except Exception as e:
    print(f"❌ Error: {e}")
    exit(1)
EOF
}

convert_ab1_to_fastq() {
    local input_file="$1"
    local output_file="$2"

    echo "Converting AB1 to FASTQ: $(basename "$input_file") → $(basename "$output_file")"

    python3 <<EOF
from Bio import SeqIO
import os
try:
    record = SeqIO.read("$input_file", "abi")
    if not record.id or record.id == "<unknown id>":
        record.id = os.path.splitext(os.path.basename("$input_file"))[0]
        record.description = "Sanger sequencing read"
    SeqIO.write(record, "$output_file", "fastq")
    print(f"✅ Converted sequence: {record.id} ({len(record.seq)} bp)")
except Exception as e:
    print(f"❌ Error: {e}")
    exit(1)
EOF
}

convert_seq_to_fasta() {
    local input_file="$1"
    local output_file="$2"

    echo "Converting SEQ to FASTA: $(basename "$input_file") → $(basename "$output_file")"

    local filename=$(basename "$input_file" .seq)
    echo ">$filename" > "$output_file"
    tr -d '[:space:]' < "$input_file" | tr -d '[:digit:]' >> "$output_file"

    echo "✅ Conversion complete"
}

convert_fasta_to_genbank() {
    local input_file="$1"
    local output_file="$2"

    echo "Converting FASTA to GenBank: $(basename "$input_file") → $(basename "$output_file")"

    python3 <<EOF
from Bio import SeqIO
from datetime import datetime
try:
    records = []
    for record in SeqIO.parse("$input_file", "fasta"):
        record.annotations["molecule_type"] = "DNA"
        record.annotations["date"] = datetime.now().strftime("%d-%b-%Y").upper()
        records.append(record)
    count = SeqIO.write(records, "$output_file", "genbank")
    print(f"✅ Converted {count} sequence(s)")
except Exception as e:
    print(f"❌ Error: {e}")
    exit(1)
EOF
}

# ==========================================
# ALIGNMENT CONVERSION FUNCTIONS
# ==========================================

convert_sam_to_bam() {
    local input_file="$1"
    local output_file="$2"

    echo "Converting SAM to BAM: $(basename "$input_file") → $(basename "$output_file")"

    samtools view -b -@ "$THREADS" -l "$COMPRESSION_LEVEL" -o "$output_file" "$input_file"

    if [ $? -eq 0 ]; then
        echo "✅ BAM file created"
        echo "Creating BAM index..."
        samtools index "$output_file"
        echo "✅ Index created"
    else
        echo "❌ Conversion failed"
        exit 1
    fi
}

convert_bam_to_sam() {
    local input_file="$1"
    local output_file="$2"

    echo "Converting BAM to SAM: $(basename "$input_file") → $(basename "$output_file")"

    samtools view -h -@ "$THREADS" -o "$output_file" "$input_file"

    [ $? -eq 0 ] && echo "✅ SAM file created" || (echo "❌ Conversion failed" && exit 1)
}

convert_bam_to_cram() {
    local input_file="$1"
    local output_file="$2"

    if [ -z "$REFERENCE_GENOME" ]; then
        echo "❌ Error: Reference genome required for CRAM conversion"
        exit 1
    fi

    echo "Converting BAM to CRAM: $(basename "$input_file") → $(basename "$output_file")"
    echo "Reference: $(basename "$REFERENCE_GENOME")"

    samtools view -C -T "$REFERENCE_GENOME" -@ "$THREADS" -o "$output_file" "$input_file"

    if [ $? -eq 0 ]; then
        echo "✅ CRAM file created"
        echo "Creating CRAM index..."
        samtools index "$output_file"
        echo "✅ Index created"
    else
        echo "❌ Conversion failed"
        exit 1
    fi
}

convert_cram_to_bam() {
    local input_file="$1"
    local output_file="$2"

    echo "Converting CRAM to BAM: $(basename "$input_file") → $(basename "$output_file")"

    local ref_arg=""
    if [ -n "$REFERENCE_GENOME" ]; then
        ref_arg="-T $REFERENCE_GENOME"
    fi

    samtools view -b $ref_arg -@ "$THREADS" -l "$COMPRESSION_LEVEL" -o "$output_file" "$input_file"

    if [ $? -eq 0 ]; then
        echo "✅ BAM file created"
        samtools index "$output_file"
        echo "✅ Index created"
    else
        echo "❌ Conversion failed"
        exit 1
    fi
}

convert_bam_to_bed() {
    local input_file="$1"
    local output_file="$2"

    echo "Converting BAM to BED: $(basename "$input_file") → $(basename "$output_file")"

    if command -v bedtools &>/dev/null; then
        bedtools bamtobed -i "$input_file" > "$output_file"
    else
        samtools view "$input_file" | \
            awk 'BEGIN {OFS="\t"} !and($2, 4) {
                print $3, $4-1, $4+length($10)-1, $1, $5, (and($2, 16) ? "-" : "+")
            }' > "$output_file"
    fi

    if [ $? -eq 0 ]; then
        local count=$(wc -l < "$output_file")
        echo "✅ BED file created ($count alignments)"
    else
        echo "❌ Conversion failed"
        exit 1
    fi
}

convert_sam_to_bed() {
    local input_file="$1"
    local output_file="$2"

    echo "Converting SAM to BED: $(basename "$input_file") → $(basename "$output_file")"
    convert_bam_to_bed <(samtools view -b "$input_file") "$output_file"
}

convert_bam_to_fastq() {
    local input_file="$1"
    local output_file="$2"

    echo "Converting BAM to FASTQ: $(basename "$input_file") → $(basename "$output_file")"

    local paired=$(samtools view -c -f 1 "$input_file" 2>/dev/null | head -1)

    if [ "$paired" -gt 0 ]; then
        echo "Detected paired-end reads"
        local output_dir=$(dirname "$output_file")
        local output_base=$(basename "$output_file" .fastq)
        local r1_file="$output_dir/${output_base}_R1.fastq"
        local r2_file="$output_dir/${output_base}_R2.fastq"
        local singleton_file="$output_dir/${output_base}_singleton.fastq"

        samtools fastq -@ "$THREADS" -1 "$r1_file" -2 "$r2_file" -s "$singleton_file" -0 /dev/null -n "$input_file"

        echo "✅ Created paired-end FASTQ files"
    else
        echo "Detected single-end reads"
        samtools fastq -@ "$THREADS" -0 "$output_file" "$input_file"
        [ $? -eq 0 ] && echo "✅ FASTQ file created" || (echo "❌ Conversion failed" && exit 1)
    fi
}

convert_sam_to_fastq() {
    local input_file="$1"
    local output_file="$2"

    echo "Converting SAM to FASTQ: $(basename "$input_file") → $(basename "$output_file")"
    convert_bam_to_fastq <(samtools view -b "$input_file") "$output_file"
}

# ==========================================
# VARIANT CONVERSION FUNCTIONS
# ==========================================

convert_vcf_to_bcf() {
    local input_file="$1"
    local output_file="$2"

    echo "Converting VCF to BCF: $(basename "$input_file") → $(basename "$output_file")"

    bcftools view -O b --threads "$THREADS" -o "$output_file" "$input_file"

    if [ $? -eq 0 ]; then
        echo "✅ BCF file created"
        bcftools index "$output_file"
        echo "✅ Index created"
    else
        echo "❌ Conversion failed"
        exit 1
    fi
}

convert_bcf_to_vcf() {
    local input_file="$1"
    local output_file="$2"

    echo "Converting BCF to VCF: $(basename "$input_file") → $(basename "$output_file")"

    bcftools view -O v --threads "$THREADS" -o "$output_file" "$input_file"

    if [ $? -eq 0 ]; then
        local count=$(bcftools view -H "$output_file" | wc -l)
        echo "✅ VCF file created ($count variants)"
    else
        echo "❌ Conversion failed"
        exit 1
    fi
}

convert_vcf_to_bed() {
    local input_file="$1"
    local output_file="$2"

    echo "Converting VCF to BED: $(basename "$input_file") → $(basename "$output_file")"

    bcftools query -f '%CHROM\t%POS0\t%END\t%ID\t%QUAL\n' "$input_file" > "$output_file" 2>/dev/null || \
        bcftools query -f '%CHROM\t%POS0\t%POS1\t%ID\t%QUAL\n' "$input_file" | \
        awk 'BEGIN {OFS="\t"} {
            id = ($4 == "." ? "variant_"NR : $4)
            qual = ($5 == "." ? "0" : $5)
            print $1, $2, $3, id, qual
        }' > "$output_file"

    if [ $? -eq 0 ]; then
        local count=$(wc -l < "$output_file")
        echo "✅ BED file created ($count variants)"
    else
        echo "❌ Conversion failed"
        exit 1
    fi
}

convert_bcf_to_bed() {
    convert_vcf_to_bed "$1" "$2"
}

convert_vcf_to_tsv() {
    local input_file="$1"
    local output_file="$2"

    echo "Converting VCF to TSV: $(basename "$input_file") → $(basename "$output_file")"

    echo -e "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" > "$output_file"
    bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\n' "$input_file" >> "$output_file"

    if [ $? -eq 0 ]; then
        local count=$(expr $(wc -l < "$output_file") - 1)
        echo "✅ TSV file created ($count variants)"
    else
        echo "❌ Conversion failed"
        exit 1
    fi
}

convert_vcf_to_csv() {
    local input_file="$1"
    local output_file="$2"

    echo "Converting VCF to CSV: $(basename "$input_file") → $(basename "$output_file")"

    echo "CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO" > "$output_file"
    bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\n' "$input_file" | \
        awk 'BEGIN {FS="\t"; OFS=","} {
            for (i=1; i<=NF; i++) {
                if ($i ~ /,/) {
                    gsub(/"/, "\"\"", $i)
                    $i = "\"" $i "\""
                }
            }
            print
        }' >> "$output_file"

    if [ $? -eq 0 ]; then
        local count=$(expr $(wc -l < "$output_file") - 1)
        echo "✅ CSV file created ($count variants)"
    else
        echo "❌ Conversion failed"
        exit 1
    fi
}

# ==========================================
# ANNOTATION CONVERSION FUNCTIONS
# ==========================================

convert_gff3_to_gtf() {
    local input_file="$1"
    local output_file="$2"

    echo "Converting GFF3 to GTF: $(basename "$input_file") → $(basename "$output_file")"

    if command -v gffread &>/dev/null; then
        gffread "$input_file" -T -o "$output_file"
    else
        awk 'BEGIN {OFS="\t"}
        !/^#/ {
            attrs = $9
            gsub(/=/, " \"", attrs)
            gsub(/;/, "\"; ", attrs)
            attrs = attrs "\""
            print $1, $2, $3, $4, $5, $6, $7, $8, attrs
        }' "$input_file" > "$output_file"
    fi

    if [ $? -eq 0 ]; then
        local count=$(grep -v "^#" "$output_file" | wc -l)
        echo "✅ GTF file created ($count features)"
    else
        echo "❌ Conversion failed"
        exit 1
    fi
}

convert_gtf_to_gff3() {
    local input_file="$1"
    local output_file="$2"

    echo "Converting GTF to GFF3: $(basename "$input_file") → $(basename "$output_file")"

    if command -v gffread &>/dev/null; then
        gffread "$input_file" -o "$output_file"
    else
        echo "##gff-version 3" > "$output_file"
        awk 'BEGIN {OFS="\t"}
        !/^#/ {
            attrs = $9
            gene_id = ""; transcript_id = ""; gene_name = ""
            if (match(attrs, /gene_id "([^"]+)"/, arr)) gene_id = arr[1]
            if (match(attrs, /transcript_id "([^"]+)"/, arr)) transcript_id = arr[1]
            if (match(attrs, /gene_name "([^"]+)"/, arr)) gene_name = arr[1]
            gff_attrs = ""
            if (gene_id) gff_attrs = gff_attrs "ID=" gene_id ";"
            if (transcript_id) gff_attrs = gff_attrs "Parent=" transcript_id ";"
            if (gene_name) gff_attrs = gff_attrs "Name=" gene_name ";"
            sub(/;$/, "", gff_attrs)
            print $1, $2, $3, $4, $5, $6, $7, $8, gff_attrs
        }' "$input_file" >> "$output_file"
    fi

    if [ $? -eq 0 ]; then
        local count=$(grep -v "^#" "$output_file" | wc -l)
        echo "✅ GFF3 file created ($count features)"
    else
        echo "❌ Conversion failed"
        exit 1
    fi
}

convert_gff_to_gtf() {
    convert_gff3_to_gtf "$1" "$2"
}

convert_gtf_to_gff() {
    convert_gtf_to_gff3 "$1" "$2"
}

convert_bed_to_gff3() {
    local input_file="$1"
    local output_file="$2"

    echo "Converting BED to GFF3: $(basename "$input_file") → $(basename "$output_file")"

    echo "##gff-version 3" > "$output_file"
    awk 'BEGIN {OFS="\t"}
    !/^#/ && !/^track/ && !/^browser/ {
        chrom = $1; start = $2 + 1; end = $3
        name = (NF >= 4 && $4 != "" ? $4 : "feature_" NR)
        score = (NF >= 5 && $5 != "" ? $5 : ".")
        strand = (NF >= 6 && $6 != "" ? $6 : ".")
        print chrom, "BED", "feature", start, end, score, strand, ".", "ID=" name ";Name=" name
    }' "$input_file" >> "$output_file"

    if [ $? -eq 0 ]; then
        local count=$(grep -v "^#" "$output_file") | wc -l)
        echo "✅ GFF3 file created ($count features)"
    else
        echo "❌ Conversion failed"
        exit 1
    fi
}

convert_gff3_to_bed() {
    local input_file="$1"
    local output_file="$2"

    echo "Converting GFF3 to BED: $(basename "$input_file") → $(basename "$output_file")"

    awk 'BEGIN {OFS="\t"}
    !/^#/ {
        chrom = $1; start = $4 - 1; end = $5
        score = ($6 == "." ? "0" : $6)
        strand = ($7 == "." ? "+" : $7)
        attrs = $9; name = "feature"
        if (match(attrs, /ID=([^;]+)/, arr)) name = arr[1]
        else if (match(attrs, /Name=([^;]+)/, arr)) name = arr[1]
        print chrom, start, end, name, score, strand
    }' "$input_file" > "$output_file"

    if [ $? -eq 0 ]; then
        local count=$(wc -l < "$output_file")
        echo "✅ BED file created ($count features)"
    else
        echo "❌ Conversion failed"
        exit 1
    fi
}

convert_gtf_to_bed() {
    local input_file="$1"
    local output_file="$2"

    echo "Converting GTF to BED: $(basename "$input_file") → $(basename "$output_file")"

    awk 'BEGIN {OFS="\t"}
    !/^#/ {
        chrom = $1; start = $4 - 1; end = $5
        score = ($6 == "." ? "0" : $6)
        strand = ($7 == "." ? "+" : $7)
        attrs = $9; name = "feature"
        if (match(attrs, /gene_id "([^"]+)"/, arr)) name = arr[1]
        else if (match(attrs, /transcript_id "([^"]+)"/, arr)) name = arr[1]
        print chrom, start, end, name, score, strand
    }' "$input_file" > "$output_file"

    if [ $? -eq 0 ]; then
        local count=$(wc -l < "$output_file")
        echo "✅ BED file created ($count features)"
    else
        echo "❌ Conversion failed"
        exit 1
    fi
}

# ==========================================
# ARCHIVE HANDLING
# ==========================================

extract_archive() {
    local archive_file="$1"
    local extract_dir="$2"

    echo "📦 Extracting archive: $(basename "$archive_file")"

    mkdir -p "$extract_dir"

    case "$archive_file" in
        *.zip)
            unzip -q -o "$archive_file" -d "$extract_dir"
            ;;
        *.tar.gz|*.tgz)
            tar -xzf "$archive_file" -C "$extract_dir"
            ;;
        *.tar.bz2)
            tar -xjf "$archive_file" -C "$extract_dir"
            ;;
        *.tar.xz)
            tar -xJf "$archive_file" -C "$extract_dir"
            ;;
        *.7z)
            7z x "$archive_file" -o"$extract_dir" >/dev/null
            ;;
        *)
            echo "❌ Unsupported archive format"
            exit 1
            ;;
    esac

    if [ $? -eq 0 ]; then
        local file_count=$(find "$extract_dir" -type f | wc -l)
        echo "✅ Extracted $file_count file(s)"
    else
        echo "❌ Extraction failed"
        exit 1
    fi
}

# ==========================================
# MAIN PROCESSING LOGIC
# ==========================================

# Get input files
INPUT_FILES=(/tmp/input/*)

if [ ${#INPUT_FILES[@]} -eq 0 ] || [ ! -e "${INPUT_FILES[0]}" ]; then
    echo "❌ Error: No input files found in /tmp/input/"
    exit 1
fi

echo "📁 Found ${#INPUT_FILES[@]} input file(s)"

# Check if archive extraction is needed
if [ "$IS_ARCHIVE" = "true" ] && [ "$EXTRACT_ARCHIVE" = "true" ]; then
    echo "📦 Archive mode: Extracting files..."
    mkdir -p /tmp/extract

    for archive in "${INPUT_FILES[@]}"; do
        extract_archive "$archive" /tmp/extract
    done

    # Update input files list to extracted files
    INPUT_FILES=(/tmp/extract/*)
    echo "📁 Processing ${#INPUT_FILES[@]} extracted file(s)"
fi

# Process each file
CONVERTED_COUNT=0
FAILED_COUNT=0

for input_file in "${INPUT_FILES[@]}"; do
    [ ! -f "$input_file" ] && continue

    filename=$(basename "$input_file")
    filename_noext="${filename%.*}"
    target_ext=".$TARGET_FORMAT"

    # Adjust extension for special cases
    case "$TARGET_FORMAT" in
        fastq) target_ext=".fq" ;;
        fasta) target_ext=".fa" ;;
        genbank) target_ext=".gb" ;;
    esac

    output_file="/tmp/output/${filename_noext}${target_ext}"

    echo ""
    echo "🔄 Processing: $filename"

    # Call appropriate conversion function
    if type "convert_$CONVERSION_COMMAND" &>/dev/null; then
        "convert_$CONVERSION_COMMAND" "$input_file" "$output_file"

        if [ $? -eq 0 ]; then
            ((CONVERTED_COUNT++))
        else
            ((FAILED_COUNT++))
        fi
    else
        echo "❌ Conversion function not found: convert_$CONVERSION_COMMAND"
        ((FAILED_COUNT++))
    fi
done

echo ""
echo "========================================="
echo "📊 Conversion Summary"
echo "========================================="
echo "✅ Successfully converted: $CONVERTED_COUNT file(s)"
echo "❌ Failed conversions: $FAILED_COUNT file(s)"
echo ""

if [ $CONVERTED_COUNT -eq 0 ]; then
    echo "❌ No files were successfully converted"
    exit 1
fi

echo "✅ Format conversion completed successfully"
exit 0
