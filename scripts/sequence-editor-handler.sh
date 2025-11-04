#!/bin/bash
# Sequence Editor Handler
# Handles DNA/RNA/protein sequence editing operations on AWS ECS
# Supports FASTA and FASTQ formats with BioPython and seqtk

set -e  # Exit on error
set -o pipefail

echo "========================================="
echo "✂️  Sequence Editor Handler Started"
echo "========================================="

# ==========================================
# Environment Variables
# ==========================================
OPERATION="${OPERATION:-}"
START_POSITION="${START_POSITION:-}"
END_POSITION="${END_POSITION:-}"
TRIM_START="${TRIM_START:-0}"
TRIM_END="${TRIM_END:-0}"
READING_FRAME="${READING_FRAME:-1}"
GENETIC_CODE="${GENETIC_CODE:-1}"
QUALITY_THRESHOLD="${QUALITY_THRESHOLD:-20}"
MIN_QUALITY_SCORE="${MIN_QUALITY_SCORE:-20}"
WINDOW_SIZE="${WINDOW_SIZE:-4}"
MIN_LENGTH="${MIN_LENGTH:-}"
MAX_LENGTH="${MAX_LENGTH:-}"
CASE_TYPE="${CASE_TYPE:-upper}"
MASK_START="${MASK_START:-}"
MASK_END="${MASK_END:-}"
MASK_CHAR="${MASK_CHAR:-N}"
SEARCH_PATTERN="${SEARCH_PATTERN:-}"
REPLACE_WITH="${REPLACE_WITH:-}"
USE_REGEX="${USE_REGEX:-false}"
SEQUENCE_IDS="${SEQUENCE_IDS:-}"
SAMPLE_SIZE="${SAMPLE_SIZE:-}"
RANDOM_SEED="${RANDOM_SEED:-}"

echo "📋 Configuration:"
echo "   Operation: $OPERATION"
echo "   Input: /tmp/input/"
echo "   Output: /tmp/output/"
echo ""

# ==========================================
# Utility Functions
# ==========================================

# Detect file format
detect_format() {
    local file="$1"
    local first_char=$(head -c 1 "$file")

    if [ "$first_char" = ">" ]; then
        echo "fasta"
    elif [ "$first_char" = "@" ]; then
        echo "fastq"
    else
        echo "unknown"
    fi
}

# Count sequences
count_sequences() {
    local file="$1"
    local format="$2"

    if [ "$format" = "fasta" ]; then
        grep -c "^>" "$file" || echo "0"
    elif [ "$format" = "fastq" ]; then
        echo $(( $(wc -l < "$file") / 4 ))
    else
        echo "0"
    fi
}

# ==========================================
# TRIM OPERATION
# ==========================================
operation_trim() {
    local input="$1"
    local output="$2"
    local format="$3"

    echo "🔪 Trimming sequences..."
    echo "   Trim start: $TRIM_START bases"
    echo "   Trim end: $TRIM_END bases"

    if [ -n "$START_POSITION" ] && [ -n "$END_POSITION" ]; then
        echo "   Extract range: $START_POSITION-$END_POSITION"
        operation_extract "$input" "$output" "$format"
        return $?
    fi

    if command -v seqtk &>/dev/null && [ "$format" = "fastq" ]; then
        # Use seqtk for FASTQ trimming
        seqtk trimfq -b "$TRIM_START" -e "$TRIM_END" "$input" > "$output"
    else
        # Use BioPython
        python3 <<EOF
from Bio import SeqIO
import sys

input_file = "$input"
output_file = "$output"
trim_start = int("$TRIM_START")
trim_end = int("$TRIM_END")
format_type = "$format"

try:
    with open(output_file, 'w') as out_handle:
        for record in SeqIO.parse(input_file, format_type):
            seq_len = len(record.seq)
            new_start = trim_start
            new_end = seq_len - trim_end if trim_end > 0 else seq_len

            if new_start < new_end:
                record.seq = record.seq[new_start:new_end]
                if format_type == "fastq" and hasattr(record, 'letter_annotations'):
                    record.letter_annotations["phred_quality"] = record.letter_annotations["phred_quality"][new_start:new_end]
                SeqIO.write(record, out_handle, format_type)

    print("✅ Trim complete", file=sys.stderr)
    sys.exit(0)
except Exception as e:
    print(f"❌ Trim failed: {e}", file=sys.stderr)
    sys.exit(1)
EOF
    fi

    echo "✅ Trim operation complete"
}

# ==========================================
# EXTRACT OPERATION
# ==========================================
operation_extract() {
    local input="$1"
    local output="$2"
    local format="$3"

    echo "📍 Extracting subsequence: $START_POSITION-$END_POSITION"

    python3 <<EOF
from Bio import SeqIO
import sys

input_file = "$input"
output_file = "$output"
start_pos = int("$START_POSITION") - 1  # Convert to 0-indexed
end_pos = int("$END_POSITION")
format_type = "$format"

try:
    with open(output_file, 'w') as out_handle:
        for record in SeqIO.parse(input_file, format_type):
            actual_start = min(start_pos, len(record.seq))
            actual_end = min(end_pos, len(record.seq))

            if actual_start < actual_end:
                record.seq = record.seq[actual_start:actual_end]
                if format_type == "fastq" and hasattr(record, 'letter_annotations'):
                    record.letter_annotations["phred_quality"] = record.letter_annotations["phred_quality"][actual_start:actual_end]
                SeqIO.write(record, out_handle, format_type)

    print("✅ Extract complete", file=sys.stderr)
    sys.exit(0)
except Exception as e:
    print(f"❌ Extract failed: {e}", file=sys.stderr)
    sys.exit(1)
EOF

    echo "✅ Extract operation complete"
}

# ==========================================
# REVERSE COMPLEMENT OPERATION
# ==========================================
operation_reverse_complement() {
    local input="$1"
    local output="$2"
    local format="$3"

    echo "🔄 Generating reverse complement..."

    if command -v seqtk &>/dev/null && [ "$format" = "fasta" ]; then
        seqtk seq -r "$input" > "$output"
    else
        python3 <<EOF
from Bio import SeqIO
import sys

input_file = "$input"
output_file = "$output"
format_type = "$format"

try:
    with open(output_file, 'w') as out_handle:
        for record in SeqIO.parse(input_file, format_type):
            rc_record = record.reverse_complement(id=True, description=True)
            SeqIO.write(rc_record, out_handle, format_type)

    print("✅ Reverse complement complete", file=sys.stderr)
    sys.exit(0)
except Exception as e:
    print(f"❌ Reverse complement failed: {e}", file=sys.stderr)
    sys.exit(1)
EOF
    fi

    echo "✅ Reverse complement operation complete"
}

# ==========================================
# TRANSLATION OPERATION
# ==========================================
operation_translate() {
    local input="$1"
    local output="$2"

    echo "🧬 Translating DNA/RNA to protein..."
    echo "   Reading frame: $READING_FRAME"
    echo "   Genetic code: $GENETIC_CODE"

    python3 <<EOF
from Bio import SeqIO
from Bio.Seq import Seq
import sys

input_file = "$input"
output_file = "$output"
reading_frame = int("$READING_FRAME")
genetic_code = int("$GENETIC_CODE")

try:
    with open(output_file, 'w') as out_handle:
        for record in SeqIO.parse(input_file, "fasta"):
            # Handle reading frame
            if reading_frame < 0:
                # Reverse complement for negative frames
                seq = record.seq.reverse_complement()
                frame_offset = abs(reading_frame) - 1
            else:
                seq = record.seq
                frame_offset = reading_frame - 1

            # Adjust for frame
            seq = seq[frame_offset:]

            # Translate
            protein_seq = seq.translate(table=genetic_code, to_stop=False)

            # Create new record
            record.seq = protein_seq
            record.id = f"{record.id}_protein"
            record.description = f"{record.description} (translated frame {reading_frame:+d})"

            SeqIO.write(record, out_handle, "fasta")

    print("✅ Translation complete", file=sys.stderr)
    sys.exit(0)
except Exception as e:
    print(f"❌ Translation failed: {e}", file=sys.stderr)
    sys.exit(1)
EOF

    echo "✅ Translation operation complete"
}

# ==========================================
# TRANSCRIPTION OPERATION
# ==========================================
operation_transcribe() {
    local input="$1"
    local output="$2"

    echo "🧬 Transcribing DNA to RNA (T → U)..."

    python3 <<EOF
from Bio import SeqIO
import sys

input_file = "$input"
output_file = "$output"

try:
    with open(output_file, 'w') as out_handle:
        for record in SeqIO.parse(input_file, "fasta"):
            record.seq = record.seq.transcribe()
            record.id = f"{record.id}_RNA"
            record.description = f"{record.description} (transcribed)"
            SeqIO.write(record, out_handle, "fasta")

    print("✅ Transcription complete", file=sys.stderr)
    sys.exit(0)
except Exception as e:
    print(f"❌ Transcription failed: {e}", file=sys.stderr)
    sys.exit(1)
EOF

    echo "✅ Transcription operation complete"
}

# ==========================================
# QUALITY FILTERING OPERATION (FASTQ)
# ==========================================
operation_quality_filter() {
    local input="$1"
    local output="$2"

    echo "📊 Filtering by quality score >= $MIN_QUALITY_SCORE..."

    python3 <<EOF
from Bio import SeqIO
import sys

input_file = "$input"
output_file = "$output"
min_quality = float("$MIN_QUALITY_SCORE")

try:
    filtered_count = 0
    total_count = 0

    with open(output_file, 'w') as out_handle:
        for record in SeqIO.parse(input_file, "fastq"):
            total_count += 1
            # Calculate average quality
            qualities = record.letter_annotations["phred_quality"]
            avg_quality = sum(qualities) / len(qualities) if qualities else 0

            if avg_quality >= min_quality:
                SeqIO.write(record, out_handle, "fastq")
                filtered_count += 1

    print(f"✅ Quality filter complete: {filtered_count}/{total_count} sequences passed", file=sys.stderr)
    sys.exit(0)
except Exception as e:
    print(f"❌ Quality filter failed: {e}", file=sys.stderr)
    sys.exit(1)
EOF

    echo "✅ Quality filtering operation complete"
}

# ==========================================
# QUALITY TRIMMING OPERATION (FASTQ)
# ==========================================
operation_quality_trim() {
    local input="$1"
    local output="$2"

    echo "✂️  Trimming by quality (threshold=$QUALITY_THRESHOLD, window=$WINDOW_SIZE)..."

    python3 <<EOF
from Bio import SeqIO
import sys

input_file = "$input"
output_file = "$output"
threshold = float("$QUALITY_THRESHOLD")
window_size = int("$WINDOW_SIZE")

try:
    trimmed_count = 0

    with open(output_file, 'w') as out_handle:
        for record in SeqIO.parse(input_file, "fastq"):
            qualities = record.letter_annotations["phred_quality"]

            # Find trim positions using sliding window
            trim_start = 0
            trim_end = len(qualities)

            # Trim from start
            for i in range(len(qualities) - window_size + 1):
                window_avg = sum(qualities[i:i+window_size]) / window_size
                if window_avg >= threshold:
                    trim_start = i
                    break

            # Trim from end
            for i in range(len(qualities) - window_size, -1, -1):
                window_avg = sum(qualities[i:i+window_size]) / window_size
                if window_avg >= threshold:
                    trim_end = i + window_size
                    break

            # Apply trim if valid
            if trim_start < trim_end and (trim_end - trim_start) >= window_size:
                record.seq = record.seq[trim_start:trim_end]
                record.letter_annotations["phred_quality"] = qualities[trim_start:trim_end]
                SeqIO.write(record, out_handle, "fastq")
                trimmed_count += 1

    print(f"✅ Quality trim complete: {trimmed_count} sequences trimmed", file=sys.stderr)
    sys.exit(0)
except Exception as e:
    print(f"❌ Quality trim failed: {e}", file=sys.stderr)
    sys.exit(1)
EOF

    echo "✅ Quality trimming operation complete"
}

# ==========================================
# LENGTH FILTERING OPERATION
# ==========================================
operation_length_filter() {
    local input="$1"
    local output="$2"
    local format="$3"

    echo "📏 Filtering by length (min=$MIN_LENGTH, max=$MAX_LENGTH)..."

    python3 <<EOF
from Bio import SeqIO
import sys

input_file = "$input"
output_file = "$output"
format_type = "$format"
min_len = int("$MIN_LENGTH") if "$MIN_LENGTH" else 0
max_len = int("$MAX_LENGTH") if "$MAX_LENGTH" else float('inf')

try:
    filtered_count = 0
    total_count = 0

    with open(output_file, 'w') as out_handle:
        for record in SeqIO.parse(input_file, format_type):
            total_count += 1
            seq_len = len(record.seq)

            if min_len <= seq_len <= max_len:
                SeqIO.write(record, out_handle, format_type)
                filtered_count += 1

    print(f"✅ Length filter complete: {filtered_count}/{total_count} sequences passed", file=sys.stderr)
    sys.exit(0)
except Exception as e:
    print(f"❌ Length filter failed: {e}", file=sys.stderr)
    sys.exit(1)
EOF

    echo "✅ Length filtering operation complete"
}

# ==========================================
# CASE CONVERSION OPERATION
# ==========================================
operation_case_convert() {
    local input="$1"
    local output="$2"
    local format="$3"

    echo "🔤 Converting case to: $CASE_TYPE"

    python3 <<EOF
from Bio import SeqIO
import sys

input_file = "$input"
output_file = "$output"
format_type = "$format"
case_type = "$CASE_TYPE"

try:
    with open(output_file, 'w') as out_handle:
        for record in SeqIO.parse(input_file, format_type):
            if case_type == "upper":
                record.seq = record.seq.upper()
            elif case_type == "lower":
                record.seq = record.seq.lower()

            SeqIO.write(record, out_handle, format_type)

    print("✅ Case conversion complete", file=sys.stderr)
    sys.exit(0)
except Exception as e:
    print(f"❌ Case conversion failed: {e}", file=sys.stderr)
    sys.exit(1)
EOF

    echo "✅ Case conversion operation complete"
}

# ==========================================
# REMOVE AMBIGUOUS BASES OPERATION
# ==========================================
operation_remove_ambiguous() {
    local input="$1"
    local output="$2"
    local format="$3"

    echo "🧹 Removing ambiguous bases (N)..."

    python3 <<EOF
from Bio import SeqIO
import sys

input_file = "$input"
output_file = "$output"
format_type = "$format"

try:
    with open(output_file, 'w') as out_handle:
        for record in SeqIO.parse(input_file, format_type):
            # Find positions to keep (not N)
            keep_positions = [i for i, base in enumerate(str(record.seq)) if base.upper() != 'N']

            if keep_positions:
                new_seq = "".join([str(record.seq)[i] for i in keep_positions])
                record.seq = type(record.seq)(new_seq)

                if format_type == "fastq" and hasattr(record, 'letter_annotations'):
                    record.letter_annotations["phred_quality"] = [record.letter_annotations["phred_quality"][i] for i in keep_positions]

                SeqIO.write(record, out_handle, format_type)

    print("✅ Remove ambiguous complete", file=sys.stderr)
    sys.exit(0)
except Exception as e:
    print(f"❌ Remove ambiguous failed: {e}", file=sys.stderr)
    sys.exit(1)
EOF

    echo "✅ Remove ambiguous operation complete"
}

# ==========================================
# MASK REGION OPERATION
# ==========================================
operation_mask_region() {
    local input="$1"
    local output="$2"
    local format="$3"

    echo "🎭 Masking region $MASK_START-$MASK_END with '$MASK_CHAR'..."

    python3 <<EOF
from Bio import SeqIO
import sys

input_file = "$input"
output_file = "$output"
format_type = "$format"
mask_start = int("$MASK_START") - 1  # Convert to 0-indexed
mask_end = int("$MASK_END")
mask_char = "$MASK_CHAR"

try:
    with open(output_file, 'w') as out_handle:
        for record in SeqIO.parse(input_file, format_type):
            seq_str = str(record.seq)
            seq_len = len(seq_str)

            actual_start = min(mask_start, seq_len)
            actual_end = min(mask_end, seq_len)

            if actual_start < actual_end:
                masked_seq = seq_str[:actual_start] + mask_char * (actual_end - actual_start) + seq_str[actual_end:]
                record.seq = type(record.seq)(masked_seq)

            SeqIO.write(record, out_handle, format_type)

    print("✅ Mask region complete", file=sys.stderr)
    sys.exit(0)
except Exception as e:
    print(f"❌ Mask region failed: {e}", file=sys.stderr)
    sys.exit(1)
EOF

    echo "✅ Mask region operation complete"
}

# ==========================================
# DEDUPLICATE OPERATION
# ==========================================
operation_deduplicate() {
    local input="$1"
    local output="$2"
    local format="$3"

    echo "🔍 Removing duplicate sequences..."

    python3 <<EOF
from Bio import SeqIO
import sys

input_file = "$input"
output_file = "$output"
format_type = "$format"

try:
    seen_sequences = set()
    unique_count = 0
    total_count = 0

    with open(output_file, 'w') as out_handle:
        for record in SeqIO.parse(input_file, format_type):
            total_count += 1
            seq_str = str(record.seq)

            if seq_str not in seen_sequences:
                seen_sequences.add(seq_str)
                SeqIO.write(record, out_handle, format_type)
                unique_count += 1

    print(f"✅ Deduplicate complete: {unique_count}/{total_count} unique sequences", file=sys.stderr)
    sys.exit(0)
except Exception as e:
    print(f"❌ Deduplicate failed: {e}", file=sys.stderr)
    sys.exit(1)
EOF

    echo "✅ Deduplicate operation complete"
}

# ==========================================
# SUBSAMPLE OPERATION
# ==========================================
operation_subsample() {
    local input="$1"
    local output="$2"
    local format="$3"

    echo "🎲 Randomly sampling $SAMPLE_SIZE sequences..."
    [ -n "$RANDOM_SEED" ] && echo "   Using random seed: $RANDOM_SEED"

    if command -v seqtk &>/dev/null; then
        if [ -n "$RANDOM_SEED" ]; then
            seqtk sample -s "$RANDOM_SEED" "$input" "$SAMPLE_SIZE" > "$output"
        else
            seqtk sample "$input" "$SAMPLE_SIZE" > "$output"
        fi
    else
        python3 <<EOF
from Bio import SeqIO
import sys
import random

input_file = "$input"
output_file = "$output"
format_type = "$format"
sample_size = int("$SAMPLE_SIZE")
random_seed = "$RANDOM_SEED"

try:
    if random_seed:
        random.seed(int(random_seed))

    # Read all sequences
    records = list(SeqIO.parse(input_file, format_type))
    total_count = len(records)

    # Sample
    sample_count = min(sample_size, total_count)
    sampled_records = random.sample(records, sample_count)

    # Write
    with open(output_file, 'w') as out_handle:
        SeqIO.write(sampled_records, out_handle, format_type)

    print(f"✅ Subsample complete: {sample_count}/{total_count} sequences sampled", file=sys.stderr)
    sys.exit(0)
except Exception as e:
    print(f"❌ Subsample failed: {e}", file=sys.stderr)
    sys.exit(1)
EOF
    fi

    echo "✅ Subsample operation complete"
}

# ==========================================
# EXTRACT BY ID OPERATION
# ==========================================
operation_extract_by_id() {
    local input="$1"
    local output="$2"
    local format="$3"

    echo "🔖 Extracting sequences by ID..."
    echo "   IDs: $SEQUENCE_IDS"

    python3 <<EOF
from Bio import SeqIO
import sys

input_file = "$input"
output_file = "$output"
format_type = "$format"
sequence_ids = "$SEQUENCE_IDS".split(',')
id_set = set(id.strip() for id in sequence_ids)

try:
    extracted_count = 0

    with open(output_file, 'w') as out_handle:
        for record in SeqIO.parse(input_file, format_type):
            if record.id in id_set:
                SeqIO.write(record, out_handle, format_type)
                extracted_count += 1

    print(f"✅ Extract by ID complete: {extracted_count} sequences extracted", file=sys.stderr)
    sys.exit(0)
except Exception as e:
    print(f"❌ Extract by ID failed: {e}", file=sys.stderr)
    sys.exit(1)
EOF

    echo "✅ Extract by ID operation complete"
}

# ==========================================
# FIND AND REPLACE OPERATION
# ==========================================
operation_find_replace() {
    local input="$1"
    local output="$2"
    local format="$3"

    echo "🔍 Find and replace operation..."
    echo "   Search: $SEARCH_PATTERN"
    echo "   Replace with: $REPLACE_WITH"
    [ "$USE_REGEX" = "true" ] && echo "   Using regex: yes"

    python3 <<EOF
from Bio import SeqIO
import sys
import re

input_file = "$input"
output_file = "$output"
format_type = "$format"
search_pattern = "$SEARCH_PATTERN"
replace_with = "$REPLACE_WITH"
use_regex = "$USE_REGEX" == "true"

try:
    replaced_count = 0
    total_count = 0

    with open(output_file, 'w') as out_handle:
        for record in SeqIO.parse(input_file, format_type):
            total_count += 1
            seq_str = str(record.seq)
            original_len = len(seq_str)

            # Perform replacement
            if use_regex:
                try:
                    new_seq = re.sub(search_pattern, replace_with, seq_str, flags=re.IGNORECASE)
                except re.error as e:
                    print(f"❌ Invalid regex pattern: {e}", file=sys.stderr)
                    sys.exit(1)
            else:
                # Case-insensitive simple replacement
                new_seq = seq_str.upper().replace(search_pattern.upper(), replace_with.upper())

            # Update record
            record.seq = type(record.seq)(new_seq)

            # Warning if length changed for FASTQ (quality scores may be invalid)
            if format_type == "fastq" and len(new_seq) != original_len:
                print(f"⚠️  Warning: Sequence length changed for {record.id} - quality scores may be invalid", file=sys.stderr)

            if new_seq != seq_str:
                replaced_count += 1

            SeqIO.write(record, out_handle, format_type)

    print(f"✅ Find and replace complete: {replaced_count}/{total_count} sequences modified", file=sys.stderr)
    sys.exit(0)
except Exception as e:
    print(f"❌ Find and replace failed: {e}", file=sys.stderr)
    sys.exit(1)
EOF

    echo "✅ Find and replace operation complete"
}

# ==========================================
# CONCATENATE OPERATION
# ==========================================
operation_concatenate() {
    local input="$1"
    local output="$2"
    local format="$3"

    echo "🔗 Concatenating all sequences into one..."

    python3 <<EOF
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import sys

input_file = "$input"
output_file = "$output"
format_type = "$format"

try:
    # Read all sequences
    records = list(SeqIO.parse(input_file, format_type))
    total_count = len(records)

    if total_count == 0:
        print("❌ No sequences found to concatenate", file=sys.stderr)
        sys.exit(1)

    # Concatenate sequences
    concatenated_seq = ""
    concatenated_quality = [] if format_type == "fastq" else None

    for record in records:
        concatenated_seq += str(record.seq)
        if format_type == "fastq" and hasattr(record, 'letter_annotations'):
            concatenated_quality.extend(record.letter_annotations["phred_quality"])

    # Create new concatenated record
    new_record = SeqRecord(
        Seq(concatenated_seq),
        id="concatenated_sequence",
        description=f"Concatenated from {total_count} sequences"
    )

    # Add quality scores for FASTQ
    if format_type == "fastq" and concatenated_quality:
        new_record.letter_annotations["phred_quality"] = concatenated_quality

    # Write output
    with open(output_file, 'w') as out_handle:
        SeqIO.write(new_record, out_handle, format_type)

    print(f"✅ Concatenate complete: {total_count} sequences → 1 sequence ({len(concatenated_seq)} bases)", file=sys.stderr)
    sys.exit(0)
except Exception as e:
    print(f"❌ Concatenate failed: {e}", file=sys.stderr)
    sys.exit(1)
EOF

    echo "✅ Concatenate operation complete"
}

# ==========================================
# Main Processing Logic
# ==========================================

# Find all input files (supports batch processing)
INPUT_FILES=($(find /tmp/input -type f \( -name "*.fasta" -o -name "*.fa" -o -name "*.fastq" -o -name "*.fq" \) | sort))

if [ ${#INPUT_FILES[@]} -eq 0 ]; then
    echo "❌ No input files found in /tmp/input/"
    exit 1
fi

echo "📁 Found ${#INPUT_FILES[@]} input file(s)"

# Check if batch mode (multiple files)
BATCH_MODE=false
if [ ${#INPUT_FILES[@]} -gt 1 ]; then
    BATCH_MODE=true
    echo "🔄 Batch mode: Processing ${#INPUT_FILES[@]} files"
fi

# Handle concatenate operation specially - process all files together
if [ "$OPERATION" = "concatenate" ] && [ "$BATCH_MODE" = "true" ]; then
    echo "🔗 Concatenating ${#INPUT_FILES[@]} files together..."

    # Detect format from first file
    FIRST_FILE="${INPUT_FILES[0]}"
    FORMAT=$(detect_format "$FIRST_FILE")
    echo "📋 Format: $FORMAT"

    # Create temporary file list for Python
    TEMP_FILE_LIST="/tmp/file_list.txt"
    printf '%s\n' "${INPUT_FILES[@]}" > "$TEMP_FILE_LIST"

    # Output file
    OUTPUT_FILE="/tmp/output/concatenated_${#INPUT_FILES[@]}_files.${FORMAT}"

    # Run concatenate with Python
    python3 <<EOF
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import sys

file_list_path = "$TEMP_FILE_LIST"
output_file = "$OUTPUT_FILE"
format_type = "$FORMAT"

try:
    # Read file list
    with open(file_list_path, 'r') as f:
        input_files = [line.strip() for line in f if line.strip()]

    print(f"📖 Reading {len(input_files)} files...", file=sys.stderr)

    # Read all sequences
    all_sequences = []
    for input_file in input_files:
        sequences = list(SeqIO.parse(input_file, format_type))
        all_sequences.extend(sequences)
        print(f"  Loaded {len(sequences)} sequences from {input_file}", file=sys.stderr)

    print(f"📊 Total sequences loaded: {len(all_sequences)}", file=sys.stderr)

    # Concatenate all sequences
    concatenated_seq = ""
    concatenated_quality = [] if format_type == "fastq" else None

    for record in all_sequences:
        concatenated_seq += str(record.seq)
        if format_type == "fastq" and hasattr(record, 'letter_annotations'):
            concatenated_quality.extend(record.letter_annotations["phred_quality"])

    # Create concatenated record
    new_record = SeqRecord(
        Seq(concatenated_seq),
        id="concatenated_sequence",
        description=f"Concatenated from {len(input_files)} files ({len(all_sequences)} sequences)"
    )

    # Add quality for FASTQ
    if format_type == "fastq" and concatenated_quality:
        new_record.letter_annotations["phred_quality"] = concatenated_quality

    # Write output
    with open(output_file, 'w') as out_handle:
        SeqIO.write(new_record, out_handle, format_type)

    print(f"✅ Concatenate complete: {len(input_files)} files → 1 sequence ({len(concatenated_seq)} bases)", file=sys.stderr)
    sys.exit(0)
except Exception as e:
    print(f"❌ Concatenate failed: {e}", file=sys.stderr)
    sys.exit(1)
EOF

    # Clean up
    rm -f "$TEMP_FILE_LIST"

    # Verify output
    if [ ! -f "$OUTPUT_FILE" ]; then
        echo "❌ Output file was not created"
        exit 1
    fi

    echo "📁 Output file: $(basename "$OUTPUT_FILE")"
    echo ""
    echo "========================================="
    echo "✅ Sequence Editor Handler Complete"
    echo "========================================="
    exit 0
fi

# Batch mode (non-concatenate): process each file separately
if [ "$BATCH_MODE" = "true" ]; then
    echo "🔄 Processing each file separately..."
    SUCCESS_COUNT=0

    for INPUT_FILE in "${INPUT_FILES[@]}"; do
        echo ""
        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        echo "📄 Processing: $(basename "$INPUT_FILE")"
        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

        # Detect format
        FORMAT=$(detect_format "$INPUT_FILE")
        echo "📋 Format: $FORMAT"

        if [ "$FORMAT" = "unknown" ]; then
            echo "⚠️  Skipping unknown format"
            continue
        fi

        # Count input sequences
        INPUT_COUNT=$(count_sequences "$INPUT_FILE" "$FORMAT")
        echo "📊 Input sequences: $INPUT_COUNT"

        # Generate output filename
        INPUT_BASENAME=$(basename "$INPUT_FILE")
        OUTPUT_EXT="${INPUT_BASENAME##*.}"
        OUTPUT_BASENAME="${INPUT_BASENAME%.*}_${OPERATION}.${OUTPUT_EXT}"
        OUTPUT_FILE="/tmp/output/$OUTPUT_BASENAME"

        # Execute operation
        case "$OPERATION" in
            trim)
                operation_trim "$INPUT_FILE" "$OUTPUT_FILE" "$FORMAT" && SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
                ;;
            extract)
                operation_extract "$INPUT_FILE" "$OUTPUT_FILE" "$FORMAT" && SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
                ;;
            reverse_complement)
                operation_reverse_complement "$INPUT_FILE" "$OUTPUT_FILE" "$FORMAT" && SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
                ;;
            translate)
                operation_translate "$INPUT_FILE" "$OUTPUT_FILE" && SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
                ;;
            transcribe)
                operation_transcribe "$INPUT_FILE" "$OUTPUT_FILE" && SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
                ;;
            quality_filter)
                operation_quality_filter "$INPUT_FILE" "$OUTPUT_FILE" && SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
                ;;
            quality_trim)
                operation_quality_trim "$INPUT_FILE" "$OUTPUT_FILE" && SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
                ;;
            length_filter)
                operation_length_filter "$INPUT_FILE" "$OUTPUT_FILE" "$FORMAT" && SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
                ;;
            case_convert)
                operation_case_convert "$INPUT_FILE" "$OUTPUT_FILE" "$FORMAT" && SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
                ;;
            remove_ambiguous)
                operation_remove_ambiguous "$INPUT_FILE" "$OUTPUT_FILE" "$FORMAT" && SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
                ;;
            mask_region)
                operation_mask_region "$INPUT_FILE" "$OUTPUT_FILE" "$FORMAT" && SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
                ;;
            deduplicate)
                operation_deduplicate "$INPUT_FILE" "$OUTPUT_FILE" "$FORMAT" && SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
                ;;
            subsample)
                operation_subsample "$INPUT_FILE" "$OUTPUT_FILE" "$FORMAT" && SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
                ;;
            extract_by_id)
                operation_extract_by_id "$INPUT_FILE" "$OUTPUT_FILE" "$FORMAT" && SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
                ;;
            find_replace)
                operation_find_replace "$INPUT_FILE" "$OUTPUT_FILE" "$FORMAT" && SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
                ;;
            *)
                echo "❌ Unknown operation: $OPERATION"
                ;;
        esac
    done

    echo ""
    echo "========================================="
    echo "✅ Batch Processing Complete"
    echo "   Processed: $SUCCESS_COUNT/${#INPUT_FILES[@]} files"
    echo "========================================="
    exit 0
fi

# Single file mode
INPUT_FILE="${INPUT_FILES[0]}"
echo "📁 Input file: $(basename "$INPUT_FILE")"

# Detect format
FORMAT=$(detect_format "$INPUT_FILE")
echo "📋 Format: $FORMAT"

if [ "$FORMAT" = "unknown" ]; then
    echo "❌ Unknown file format"
    exit 1
fi

# Count input sequences
INPUT_COUNT=$(count_sequences "$INPUT_FILE" "$FORMAT")
echo "📊 Input sequences: $INPUT_COUNT"
echo ""

# Generate output filename
INPUT_BASENAME=$(basename "$INPUT_FILE")
OUTPUT_EXT="${INPUT_BASENAME##*.}"
OUTPUT_BASENAME="${INPUT_BASENAME%.*}_${OPERATION}.${OUTPUT_EXT}"
OUTPUT_FILE="/tmp/output/$OUTPUT_BASENAME"

echo "🎯 Executing operation: $OPERATION"
echo ""

# Route to appropriate operation
case "$OPERATION" in
    trim)
        operation_trim "$INPUT_FILE" "$OUTPUT_FILE" "$FORMAT"
        ;;
    extract)
        operation_extract "$INPUT_FILE" "$OUTPUT_FILE" "$FORMAT"
        ;;
    reverse_complement)
        operation_reverse_complement "$INPUT_FILE" "$OUTPUT_FILE" "$FORMAT"
        ;;
    translate)
        operation_translate "$INPUT_FILE" "$OUTPUT_FILE"
        ;;
    transcribe)
        operation_transcribe "$INPUT_FILE" "$OUTPUT_FILE"
        ;;
    quality_filter)
        operation_quality_filter "$INPUT_FILE" "$OUTPUT_FILE"
        ;;
    quality_trim)
        operation_quality_trim "$INPUT_FILE" "$OUTPUT_FILE"
        ;;
    length_filter)
        operation_length_filter "$INPUT_FILE" "$OUTPUT_FILE" "$FORMAT"
        ;;
    case_convert)
        operation_case_convert "$INPUT_FILE" "$OUTPUT_FILE" "$FORMAT"
        ;;
    remove_ambiguous)
        operation_remove_ambiguous "$INPUT_FILE" "$OUTPUT_FILE" "$FORMAT"
        ;;
    mask_region)
        operation_mask_region "$INPUT_FILE" "$OUTPUT_FILE" "$FORMAT"
        ;;
    deduplicate)
        operation_deduplicate "$INPUT_FILE" "$OUTPUT_FILE" "$FORMAT"
        ;;
    subsample)
        operation_subsample "$INPUT_FILE" "$OUTPUT_FILE" "$FORMAT"
        ;;
    extract_by_id)
        operation_extract_by_id "$INPUT_FILE" "$OUTPUT_FILE" "$FORMAT"
        ;;
    find_replace)
        operation_find_replace "$INPUT_FILE" "$OUTPUT_FILE" "$FORMAT"
        ;;
    concatenate)
        operation_concatenate "$INPUT_FILE" "$OUTPUT_FILE" "$FORMAT"
        ;;
    *)
        echo "❌ Unknown operation: $OPERATION"
        exit 1
        ;;
esac

# Verify output
if [ ! -f "$OUTPUT_FILE" ]; then
    echo "❌ Output file was not created"
    exit 1
fi

# Count output sequences
OUTPUT_COUNT=$(count_sequences "$OUTPUT_FILE" "$FORMAT")
echo ""
echo "📊 Output sequences: $OUTPUT_COUNT"
echo "📁 Output file: $(basename "$OUTPUT_FILE")"
echo ""
echo "========================================="
echo "✅ Sequence Editor Handler Complete"
echo "========================================="

exit 0
