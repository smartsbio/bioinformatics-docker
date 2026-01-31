#!/bin/bash
set -e

# Zip Handler Script
# Handles compression, decompression, and listing of various archive formats

echo "📦 Zip Handler Started"

# Parse additional parameters from environment
COMPRESSION_TYPE=${COMPRESSION_TYPE:-"zip"}
COMPRESSION_LEVEL=${COMPRESSION_LEVEL:-"6"}
PASSWORD=${PASSWORD:-""}
PRESERVE_PATHS=${PRESERVE_PATHS:-"true"}
RECURSIVE=${RECURSIVE:-"true"}
OUTPUT_FILENAME=${OUTPUT_FILENAME:-"archive"}

echo "🎯 Operation: $COMMAND"
echo "📁 Compression type: $COMPRESSION_TYPE"

case "$COMMAND" in
    "compress")
        echo "🗜️ Running compression operation..."

        # Get all input files (recursively to handle subdirectories)
        mapfile -t INPUT_FILES < <(find /tmp/input -type f -o -type d -mindepth 1 -maxdepth 1)
        if [[ ${#INPUT_FILES[@]} -eq 0 ]]; then
            echo "❌ No input files found in /tmp/input/"
            exit 1
        fi

        echo "📁 Found ${#INPUT_FILES[@]} file(s) to compress"

        # Convert compression type to lowercase for case-insensitive matching
        COMPRESSION_TYPE_LOWER=$(echo "$COMPRESSION_TYPE" | tr '[:upper:]' '[:lower:]')

        case "$COMPRESSION_TYPE_LOWER" in
            "zip")
                OUTPUT_FILE="${OUTPUT_FILENAME}.zip"
                echo "🚀 Creating ZIP archive: $OUTPUT_FILE"

                ZIP_CMD="zip"

                # Add compression level (-0 to -9)
                if [[ -n "$COMPRESSION_LEVEL" ]]; then
                    ZIP_CMD="$ZIP_CMD -$COMPRESSION_LEVEL"
                fi

                # Add password protection if specified
                if [[ -n "$PASSWORD" ]]; then
                    ZIP_CMD="$ZIP_CMD -P \"$PASSWORD\""
                fi

                # Add recursive flag
                if [[ "$RECURSIVE" == "true" ]]; then
                    ZIP_CMD="$ZIP_CMD -r"
                fi

                ZIP_CMD="$ZIP_CMD /tmp/output/$OUTPUT_FILE"

                # Add all input files
                for file in "${INPUT_FILES[@]}"; do
                    ZIP_CMD="$ZIP_CMD \"$file\""
                done

                echo "🚀 Executing: $ZIP_CMD"

                if eval "$ZIP_CMD"; then
                    echo "✅ ZIP compression completed successfully"
                    OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILE" 2>/dev/null || stat -f%z "/tmp/output/$OUTPUT_FILE")
                    echo "📊 Output file: $OUTPUT_FILE ($OUTPUT_SIZE bytes)"
                else
                    echo "❌ ZIP compression failed"
                    exit 1
                fi
                ;;

            "gzip"|"gz")
                echo "🚀 Creating GZIP compressed file(s)"

                # GZIP compresses files individually
                for file in "${INPUT_FILES[@]}"; do
                    filename=$(basename "$file")
                    echo "📄 Compressing: $filename"

                    GZIP_CMD="gzip"

                    # Add compression level (-1 to -9)
                    if [[ -n "$COMPRESSION_LEVEL" ]]; then
                        GZIP_CMD="$GZIP_CMD -$COMPRESSION_LEVEL"
                    fi

                    # Keep original file and output to stdout
                    GZIP_CMD="$GZIP_CMD -c \"$file\""

                    if eval "$GZIP_CMD > /tmp/output/${filename}.gz"; then
                        OUTPUT_SIZE=$(stat -c%s "/tmp/output/${filename}.gz" 2>/dev/null || stat -f%z "/tmp/output/${filename}.gz")
                        echo "✅ Compressed: ${filename}.gz ($OUTPUT_SIZE bytes)"
                    else
                        echo "❌ GZIP compression failed for $filename"
                        exit 1
                    fi
                done

                echo "✅ GZIP compression completed successfully"
                ;;

            "bzip2"|"bz2")
                echo "🚀 Creating BZIP2 compressed file(s)"

                # BZIP2 compresses files individually
                for file in "${INPUT_FILES[@]}"; do
                    filename=$(basename "$file")
                    echo "📄 Compressing: $filename"

                    BZIP2_CMD="bzip2"

                    # Add compression level (-1 to -9)
                    if [[ -n "$COMPRESSION_LEVEL" ]]; then
                        BZIP2_CMD="$BZIP2_CMD -$COMPRESSION_LEVEL"
                    fi

                    # Keep original file and output to stdout
                    BZIP2_CMD="$BZIP2_CMD -c \"$file\""

                    if eval "$BZIP2_CMD > /tmp/output/${filename}.bz2"; then
                        OUTPUT_SIZE=$(stat -c%s "/tmp/output/${filename}.bz2" 2>/dev/null || stat -f%z "/tmp/output/${filename}.bz2")
                        echo "✅ Compressed: ${filename}.bz2 ($OUTPUT_SIZE bytes)"
                    else
                        echo "❌ BZIP2 compression failed for $filename"
                        exit 1
                    fi
                done

                echo "✅ BZIP2 compression completed successfully"
                ;;

            "xz")
                echo "🚀 Creating XZ compressed file(s)"

                # XZ compresses files individually
                for file in "${INPUT_FILES[@]}"; do
                    filename=$(basename "$file")
                    echo "📄 Compressing: $filename"

                    XZ_CMD="xz"

                    # Add compression level (-0 to -9)
                    if [[ -n "$COMPRESSION_LEVEL" ]]; then
                        XZ_CMD="$XZ_CMD -$COMPRESSION_LEVEL"
                    fi

                    # Keep original file and output to stdout
                    XZ_CMD="$XZ_CMD -c \"$file\""

                    if eval "$XZ_CMD > /tmp/output/${filename}.xz"; then
                        OUTPUT_SIZE=$(stat -c%s "/tmp/output/${filename}.xz" 2>/dev/null || stat -f%z "/tmp/output/${filename}.xz")
                        echo "✅ Compressed: ${filename}.xz ($OUTPUT_SIZE bytes)"
                    else
                        echo "❌ XZ compression failed for $filename"
                        exit 1
                    fi
                done

                echo "✅ XZ compression completed successfully"
                ;;

            "7z"|"7zip")
                OUTPUT_FILE="${OUTPUT_FILENAME}.7z"
                echo "🚀 Creating 7-Zip archive: $OUTPUT_FILE"

                SEVENZ_CMD="7z a"

                # Add compression level (0-9, where 0=copy, 9=ultra)
                # 7z uses -mx parameter for compression level
                if [[ -n "$COMPRESSION_LEVEL" ]]; then
                    SEVENZ_CMD="$SEVENZ_CMD -mx=$COMPRESSION_LEVEL"
                fi

                # Add password protection if specified
                if [[ -n "$PASSWORD" ]]; then
                    SEVENZ_CMD="$SEVENZ_CMD -p\"$PASSWORD\""
                fi

                SEVENZ_CMD="$SEVENZ_CMD /tmp/output/$OUTPUT_FILE"

                # Add all input files
                for file in "${INPUT_FILES[@]}"; do
                    SEVENZ_CMD="$SEVENZ_CMD \"$file\""
                done

                echo "🚀 Executing: 7z a [options] /tmp/output/$OUTPUT_FILE [files...]"

                if eval "$SEVENZ_CMD"; then
                    echo "✅ 7-Zip compression completed successfully"
                    OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILE" 2>/dev/null || stat -f%z "/tmp/output/$OUTPUT_FILE")
                    echo "📊 Output file: $OUTPUT_FILE ($OUTPUT_SIZE bytes)"
                else
                    echo "❌ 7-Zip compression failed"
                    exit 1
                fi
                ;;

            "tar.gz"|"tgz")
                OUTPUT_FILE="${OUTPUT_FILENAME}.tar.gz"
                echo "🚀 Creating TAR.GZ archive: $OUTPUT_FILE"

                TAR_CMD="tar -czf /tmp/output/$OUTPUT_FILE"

                # Add all input files
                for file in "${INPUT_FILES[@]}"; do
                    TAR_CMD="$TAR_CMD \"$file\""
                done

                echo "🚀 Executing: $TAR_CMD"

                if eval "$TAR_CMD"; then
                    echo "✅ TAR.GZ compression completed successfully"
                    OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILE" 2>/dev/null || stat -f%z "/tmp/output/$OUTPUT_FILE")
                    echo "📊 Output file: $OUTPUT_FILE ($OUTPUT_SIZE bytes)"
                else
                    echo "❌ TAR.GZ compression failed"
                    exit 1
                fi
                ;;

            "tar.bz2"|"tbz2")
                OUTPUT_FILE="${OUTPUT_FILENAME}.tar.bz2"
                echo "🚀 Creating TAR.BZ2 archive: $OUTPUT_FILE"

                TAR_CMD="tar -cjf /tmp/output/$OUTPUT_FILE"

                # Add all input files
                for file in "${INPUT_FILES[@]}"; do
                    TAR_CMD="$TAR_CMD \"$file\""
                done

                echo "🚀 Executing: $TAR_CMD"

                if eval "$TAR_CMD"; then
                    echo "✅ TAR.BZ2 compression completed successfully"
                    OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILE" 2>/dev/null || stat -f%z "/tmp/output/$OUTPUT_FILE")
                    echo "📊 Output file: $OUTPUT_FILE ($OUTPUT_SIZE bytes)"
                else
                    echo "❌ TAR.BZ2 compression failed"
                    exit 1
                fi
                ;;

            "tar.xz"|"txz")
                OUTPUT_FILE="${OUTPUT_FILENAME}.tar.xz"
                echo "🚀 Creating TAR.XZ archive: $OUTPUT_FILE"

                TAR_CMD="tar -cJf /tmp/output/$OUTPUT_FILE"

                # Add all input files
                for file in "${INPUT_FILES[@]}"; do
                    TAR_CMD="$TAR_CMD \"$file\""
                done

                echo "🚀 Executing: $TAR_CMD"

                if eval "$TAR_CMD"; then
                    echo "✅ TAR.XZ compression completed successfully"
                    OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILE" 2>/dev/null || stat -f%z "/tmp/output/$OUTPUT_FILE")
                    echo "📊 Output file: $OUTPUT_FILE ($OUTPUT_SIZE bytes)"
                else
                    echo "❌ TAR.XZ compression failed"
                    exit 1
                fi
                ;;

            *)
                echo "❌ Unsupported compression type: $COMPRESSION_TYPE"
                echo "Supported types: zip, gzip/gz, bzip2/bz2, xz, 7z/7zip, tar.gz/tgz, tar.bz2/tbz2, tar.xz/txz"
                exit 1
                ;;
        esac
        ;;

    "decompress")
        echo "📂 Running decompression operation..."

        # Get input file (find first file in subdirectories)
        INPUT_FILE=$(find /tmp/input -type f | head -n 1)
        if [[ -z "$INPUT_FILE" ]]; then
            echo "❌ No input files found in /tmp/input/"
            exit 1
        fi
        INPUT_FILENAME=$(basename "$INPUT_FILE")

        echo "📁 Decompressing file: $INPUT_FILENAME"

        # Auto-detect compression type from extension if not specified
        if [[ "$COMPRESSION_TYPE" == "auto" ]] || [[ -z "$COMPRESSION_TYPE" ]]; then
            case "$INPUT_FILENAME" in
                *.zip)
                    COMPRESSION_TYPE="zip"
                    ;;
                *.gz)
                    if [[ "$INPUT_FILENAME" == *.tar.gz ]] || [[ "$INPUT_FILENAME" == *.tgz ]]; then
                        COMPRESSION_TYPE="tar.gz"
                    else
                        COMPRESSION_TYPE="gzip"
                    fi
                    ;;
                *.bz2)
                    if [[ "$INPUT_FILENAME" == *.tar.bz2 ]] || [[ "$INPUT_FILENAME" == *.tbz2 ]]; then
                        COMPRESSION_TYPE="tar.bz2"
                    else
                        COMPRESSION_TYPE="bzip2"
                    fi
                    ;;
                *.xz)
                    if [[ "$INPUT_FILENAME" == *.tar.xz ]] || [[ "$INPUT_FILENAME" == *.txz ]]; then
                        COMPRESSION_TYPE="tar.xz"
                    else
                        COMPRESSION_TYPE="xz"
                    fi
                    ;;
                *.7z)
                    COMPRESSION_TYPE="7z"
                    ;;
                *)
                    echo "❌ Unable to auto-detect compression type from extension"
                    exit 1
                    ;;
            esac
            echo "🔍 Auto-detected compression type: $COMPRESSION_TYPE"
        fi

        COMPRESSION_TYPE_LOWER=$(echo "$COMPRESSION_TYPE" | tr '[:upper:]' '[:lower:]')

        # Helper function to organize extracted files
        # If archive contains single file/folder at root -> extract directly
        # If archive contains multiple files/folders at root -> create containing folder
        organize_extracted_files() {
            local archive_name="$1"

            # Get base name without extension
            local base_name="${archive_name%%.*}"

            # Count items at root level of /tmp/extract
            local item_count=$(ls -A /tmp/extract | wc -l)

            echo "📊 Archive contains $item_count top-level item(s)"

            if [[ $item_count -eq 1 ]]; then
                # Single item - move directly to /tmp/output
                echo "📦 Single item detected - extracting to root level"
                mv /tmp/extract/* /tmp/output/
            else
                # Multiple items - create containing folder with archive name
                echo "📦 Multiple items detected - creating containing folder: $base_name"
                mkdir -p "/tmp/output/$base_name"
                mv /tmp/extract/* "/tmp/output/$base_name/"
            fi

            # Clean up temp extraction directory
            rm -rf /tmp/extract
        }

        case "$COMPRESSION_TYPE_LOWER" in
            "zip")
                echo "🚀 Extracting ZIP archive"

                # Extract to temporary directory first
                mkdir -p /tmp/extract

                UNZIP_CMD="unzip"

                # Add password if specified
                if [[ -n "$PASSWORD" ]]; then
                    UNZIP_CMD="$UNZIP_CMD -P \"$PASSWORD\""
                fi

                # Overwrite files without prompting
                UNZIP_CMD="$UNZIP_CMD -o"

                UNZIP_CMD="$UNZIP_CMD \"$INPUT_FILE\" -d /tmp/extract"

                echo "🚀 Executing: unzip [options] \"$INPUT_FILENAME\" -d /tmp/extract"

                if eval "$UNZIP_CMD"; then
                    echo "✅ ZIP extraction completed successfully"
                    organize_extracted_files "$INPUT_FILENAME"
                    echo "📁 Final output:"
                    ls -lah /tmp/output/
                else
                    echo "❌ ZIP extraction failed"
                    exit 1
                fi
                ;;

            "gzip"|"gz")
                echo "🚀 Decompressing GZIP file"

                # Remove .gz extension for output filename
                OUTPUT_FILENAME="${INPUT_FILENAME%.gz}"

                GUNZIP_CMD="gunzip -c \"$INPUT_FILE\""

                if eval "$GUNZIP_CMD > /tmp/output/$OUTPUT_FILENAME"; then
                    echo "✅ GZIP decompression completed successfully"
                    OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILENAME" 2>/dev/null || stat -f%z "/tmp/output/$OUTPUT_FILENAME")
                    echo "📊 Output file: $OUTPUT_FILENAME ($OUTPUT_SIZE bytes)"
                else
                    echo "❌ GZIP decompression failed"
                    exit 1
                fi
                ;;

            "bzip2"|"bz2")
                echo "🚀 Decompressing BZIP2 file"

                # Remove .bz2 extension for output filename
                OUTPUT_FILENAME="${INPUT_FILENAME%.bz2}"

                BUNZIP2_CMD="bunzip2 -c \"$INPUT_FILE\""

                if eval "$BUNZIP2_CMD > /tmp/output/$OUTPUT_FILENAME"; then
                    echo "✅ BZIP2 decompression completed successfully"
                    OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILENAME" 2>/dev/null || stat -f%z "/tmp/output/$OUTPUT_FILENAME")
                    echo "📊 Output file: $OUTPUT_FILENAME ($OUTPUT_SIZE bytes)"
                else
                    echo "❌ BZIP2 decompression failed"
                    exit 1
                fi
                ;;

            "xz")
                echo "🚀 Decompressing XZ file"

                # Remove .xz extension for output filename
                OUTPUT_FILENAME="${INPUT_FILENAME%.xz}"

                UNXZ_CMD="unxz -c \"$INPUT_FILE\""

                if eval "$UNXZ_CMD > /tmp/output/$OUTPUT_FILENAME"; then
                    echo "✅ XZ decompression completed successfully"
                    OUTPUT_SIZE=$(stat -c%s "/tmp/output/$OUTPUT_FILENAME" 2>/dev/null || stat -f%z "/tmp/output/$OUTPUT_FILENAME")
                    echo "📊 Output file: $OUTPUT_FILENAME ($OUTPUT_SIZE bytes)"
                else
                    echo "❌ XZ decompression failed"
                    exit 1
                fi
                ;;

            "7z"|"7zip")
                echo "🚀 Extracting 7-Zip archive"

                # Extract to temporary directory first
                mkdir -p /tmp/extract

                SEVENZ_CMD="7z x"

                # Add password if specified
                if [[ -n "$PASSWORD" ]]; then
                    SEVENZ_CMD="$SEVENZ_CMD -p\"$PASSWORD\""
                fi

                # Output to directory
                SEVENZ_CMD="$SEVENZ_CMD -o/tmp/extract"

                # Assume yes on all queries
                SEVENZ_CMD="$SEVENZ_CMD -y"

                SEVENZ_CMD="$SEVENZ_CMD \"$INPUT_FILE\""

                echo "🚀 Executing: 7z x [options] \"$INPUT_FILENAME\""

                if eval "$SEVENZ_CMD"; then
                    echo "✅ 7-Zip extraction completed successfully"
                    organize_extracted_files "$INPUT_FILENAME"
                    echo "📁 Final output:"
                    ls -lah /tmp/output/
                else
                    echo "❌ 7-Zip extraction failed"
                    exit 1
                fi
                ;;

            "tar.gz"|"tgz")
                echo "🚀 Extracting TAR.GZ archive"

                # Extract to temporary directory first
                mkdir -p /tmp/extract

                TAR_CMD="tar -xzf \"$INPUT_FILE\" -C /tmp/extract"

                if eval "$TAR_CMD"; then
                    echo "✅ TAR.GZ extraction completed successfully"
                    organize_extracted_files "$INPUT_FILENAME"
                    echo "📁 Final output:"
                    ls -lah /tmp/output/
                else
                    echo "❌ TAR.GZ extraction failed"
                    exit 1
                fi
                ;;

            "tar.bz2"|"tbz2")
                echo "🚀 Extracting TAR.BZ2 archive"

                # Extract to temporary directory first
                mkdir -p /tmp/extract

                TAR_CMD="tar -xjf \"$INPUT_FILE\" -C /tmp/extract"

                if eval "$TAR_CMD"; then
                    echo "✅ TAR.BZ2 extraction completed successfully"
                    organize_extracted_files "$INPUT_FILENAME"
                    echo "📁 Final output:"
                    ls -lah /tmp/output/
                else
                    echo "❌ TAR.BZ2 extraction failed"
                    exit 1
                fi
                ;;

            "tar.xz"|"txz")
                echo "🚀 Extracting TAR.XZ archive"

                # Extract to temporary directory first
                mkdir -p /tmp/extract

                TAR_CMD="tar -xJf \"$INPUT_FILE\" -C /tmp/extract"

                if eval "$TAR_CMD"; then
                    echo "✅ TAR.XZ extraction completed successfully"
                    organize_extracted_files "$INPUT_FILENAME"
                    echo "📁 Final output:"
                    ls -lah /tmp/output/
                else
                    echo "❌ TAR.XZ extraction failed"
                    exit 1
                fi
                ;;

            *)
                echo "❌ Unsupported compression type: $COMPRESSION_TYPE"
                echo "Supported types: zip, gzip/gz, bzip2/bz2, xz, 7z/7zip, tar.gz/tgz, tar.bz2/tbz2, tar.xz/txz"
                exit 1
                ;;
        esac
        ;;

    "list")
        echo "📋 Listing archive contents..."

        # Get input file (find first file in subdirectories)
        INPUT_FILE=$(find /tmp/input -type f | head -n 1)
        if [[ -z "$INPUT_FILE" ]]; then
            echo "❌ No input files found in /tmp/input/"
            exit 1
        fi
        INPUT_FILENAME=$(basename "$INPUT_FILE")

        echo "📁 Archive: $INPUT_FILENAME"

        # Auto-detect compression type from extension
        case "$INPUT_FILENAME" in
            *.zip)
                COMPRESSION_TYPE="zip"
                ;;
            *.7z)
                COMPRESSION_TYPE="7z"
                ;;
            *.tar.gz|*.tgz)
                COMPRESSION_TYPE="tar.gz"
                ;;
            *.tar.bz2|*.tbz2)
                COMPRESSION_TYPE="tar.bz2"
                ;;
            *.tar.xz|*.txz)
                COMPRESSION_TYPE="tar.xz"
                ;;
            *)
                echo "❌ Unable to detect archive type from extension"
                exit 1
                ;;
        esac

        COMPRESSION_TYPE_LOWER=$(echo "$COMPRESSION_TYPE" | tr '[:upper:]' '[:lower:]')

        case "$COMPRESSION_TYPE_LOWER" in
            "zip")
                LIST_CMD="unzip -l"

                if [[ -n "$PASSWORD" ]]; then
                    LIST_CMD="$LIST_CMD -P \"$PASSWORD\""
                fi

                LIST_CMD="$LIST_CMD \"$INPUT_FILE\""

                if eval "$LIST_CMD > /tmp/output/contents.txt"; then
                    echo "✅ ZIP listing completed"
                    cat /tmp/output/contents.txt
                else
                    echo "❌ ZIP listing failed"
                    exit 1
                fi
                ;;

            "7z"|"7zip")
                LIST_CMD="7z l"

                if [[ -n "$PASSWORD" ]]; then
                    LIST_CMD="$LIST_CMD -p\"$PASSWORD\""
                fi

                LIST_CMD="$LIST_CMD \"$INPUT_FILE\""

                if eval "$LIST_CMD > /tmp/output/contents.txt"; then
                    echo "✅ 7-Zip listing completed"
                    cat /tmp/output/contents.txt
                else
                    echo "❌ 7-Zip listing failed"
                    exit 1
                fi
                ;;

            "tar.gz"|"tgz")
                if tar -tzf "$INPUT_FILE" > /tmp/output/contents.txt; then
                    echo "✅ TAR.GZ listing completed"
                    cat /tmp/output/contents.txt
                else
                    echo "❌ TAR.GZ listing failed"
                    exit 1
                fi
                ;;

            "tar.bz2"|"tbz2")
                if tar -tjf "$INPUT_FILE" > /tmp/output/contents.txt; then
                    echo "✅ TAR.BZ2 listing completed"
                    cat /tmp/output/contents.txt
                else
                    echo "❌ TAR.BZ2 listing failed"
                    exit 1
                fi
                ;;

            "tar.xz"|"txz")
                if tar -tJf "$INPUT_FILE" > /tmp/output/contents.txt; then
                    echo "✅ TAR.XZ listing completed"
                    cat /tmp/output/contents.txt
                else
                    echo "❌ TAR.XZ listing failed"
                    exit 1
                fi
                ;;

            *)
                echo "❌ Unsupported archive type: $COMPRESSION_TYPE"
                exit 1
                ;;
        esac
        ;;

    *)
        echo "❌ Unsupported command: $COMMAND"
        echo "Supported commands: compress, decompress, list"
        exit 1
        ;;
esac

echo "🎯 Zip handler completed successfully"
