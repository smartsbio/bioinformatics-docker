# Bioinformatics Docker Images

Custom Docker images for bioinformatics tools with AWS integration, designed for the smarts.bio platform.

## 🎯 Overview

This repository contains Docker images that combine popular bioinformatics tools with AWS CLI and cloud-native processing capabilities. Each image is built for scalable, containerized bioinformatics workflows.

**Currently includes 13+ bioinformatics tools:**
- SAMtools, BWA, Bowtie2, BEDTools, FastQC, Trimmomatic, Picard, VCFtools, FMLRC, HOMER, ANNOVAR (placeholder), GenomicRanges (R), and **Zip Toolkit** (NEW!)
- **Compression utilities**: gzip, bzip2, xz, zip, 7z, tar

## 🧬 Currently Supported Tools

### SAMtools
- **Version**: 1.19.2
- **Supported Commands**:
  - `view` - Format conversion (SAM ↔ BAM ↔ CRAM)
  - `sort` - Sort alignment files
  - `index` - Create index files
  - `stats` - Generate alignment statistics
  - `flagstat` - Flag statistics
  - `depth` - Per-position depth
  - `coverage` - Coverage statistics
  - And more...

### Zip Toolkit (NEW!)
- **Compression Tools**: gzip, bzip2, xz, zip, 7z, p7zip
- **Supported Commands**:
  - `compress` - Create compressed archives
  - `decompress` - Extract from archives
  - `list` - List archive contents
- **Supported Formats**:
  - ZIP - Standard zip archives with optional password protection
  - GZIP - Individual file compression (.gz)
  - BZIP2 - High compression (.bz2)
  - XZ - Maximum compression (.xz)
  - 7-Zip - High compression archives with password support
  - TAR.GZ - Tar archives with gzip compression
  - TAR.BZ2 - Tar archives with bzip2 compression
  - TAR.XZ - Tar archives with xz compression

## 🚀 Usage

### Environment Variables

| Variable | Required | Description | Example |
|----------|----------|-------------|---------|
| `TOOL_NAME` | ✅ | Bioinformatics tool to use | `samtools`, `zip` |
| `COMMAND` | ✅ | Tool-specific command | `view`, `compress` |
| `S3_BUCKET` | ✅ | S3 bucket for input/output files | `my-bucket` |
| `INPUT_S3_KEY` | ✅ | S3 key for input file | `data/input.sam` |
| `OUTPUT_PATH` | ✅ | S3 path for output files | `results/` |
| `OUTPUT_FORMAT` | ❌ | Output format (for applicable tools) | `bam` |
| `OUTPUT_FILE` | ❌ | Output filename | `output.bam` |
| `COMPRESSION_TYPE` | ❌ | Compression format (zip tool) | `zip`, `gzip`, `tar.gz` |
| `COMPRESSION_LEVEL` | ❌ | Compression level 0-9 (zip tool) | `6` |
| `PASSWORD` | ❌ | Archive password (zip/7z only) | `secret123` |

### Docker Run Examples

#### SAMtools Example

```bash
docker run --rm \
  -e TOOL_NAME=samtools \
  -e COMMAND=view \
  -e OUTPUT_FORMAT=bam \
  -e OUTPUT_FILE=converted.bam \
  -e S3_BUCKET=my-bucket \
  -e INPUT_S3_KEY=data/sample.sam \
  -e OUTPUT_PATH=results/ \
  -e AWS_ACCESS_KEY_ID=$AWS_ACCESS_KEY_ID \
  -e AWS_SECRET_ACCESS_KEY=$AWS_SECRET_ACCESS_KEY \
  -e AWS_DEFAULT_REGION=us-east-1 \
  smartsbio/bioinformatics-docker:latest
```

#### Zip Toolkit Examples

**Compress files to ZIP:**
```bash
docker run --rm \
  -e TOOL_NAME=zip \
  -e COMMAND=compress \
  -e COMPRESSION_TYPE=zip \
  -e COMPRESSION_LEVEL=9 \
  -e OUTPUT_FILENAME=data-archive \
  -e S3_BUCKET=my-bucket \
  -e INPUT_S3_KEY=data/file1.fastq \
  -e OUTPUT_PATH=archives/ \
  -e AWS_ACCESS_KEY_ID=$AWS_ACCESS_KEY_ID \
  -e AWS_SECRET_ACCESS_KEY=$AWS_SECRET_ACCESS_KEY \
  -e AWS_DEFAULT_REGION=us-east-1 \
  smartsbio/bioinformatics-docker:latest
```

**Decompress GZIP file:**
```bash
docker run --rm \
  -e TOOL_NAME=zip \
  -e COMMAND=decompress \
  -e COMPRESSION_TYPE=gzip \
  -e S3_BUCKET=my-bucket \
  -e INPUT_S3_KEY=data/sequences.fastq.gz \
  -e OUTPUT_PATH=extracted/ \
  -e AWS_ACCESS_KEY_ID=$AWS_ACCESS_KEY_ID \
  -e AWS_SECRET_ACCESS_KEY=$AWS_SECRET_ACCESS_KEY \
  -e AWS_DEFAULT_REGION=us-east-1 \
  smartsbio/bioinformatics-docker:latest
```

**List archive contents:**
```bash
docker run --rm \
  -e TOOL_NAME=zip \
  -e COMMAND=list \
  -e COMPRESSION_TYPE=zip \
  -e S3_BUCKET=my-bucket \
  -e INPUT_S3_KEY=archives/data.zip \
  -e OUTPUT_PATH=temp/ \
  -e AWS_ACCESS_KEY_ID=$AWS_ACCESS_KEY_ID \
  -e AWS_SECRET_ACCESS_KEY=$AWS_SECRET_ACCESS_KEY \
  -e AWS_DEFAULT_REGION=us-east-1 \
  smartsbio/bioinformatics-docker:latest
```

### AWS ECS Integration

This image is designed to work seamlessly with AWS ECS Fargate:

```json
{
  "family": "bioinformatics-task",
  "taskRoleArn": "arn:aws:iam::account:role/ecsTaskRole",
  "executionRoleArn": "arn:aws:iam::account:role/ecsTaskExecutionRole",
  "networkMode": "awsvpc",
  "requiresCompatibilities": ["FARGATE"],
  "cpu": "1024",
  "memory": "2048",
  "containerDefinitions": [
    {
      "name": "bioinformatics-container",
      "image": "smartsbio/bioinformatics-docker:latest",
      "environment": [
        {"name": "TOOL_NAME", "value": "samtools"},
        {"name": "COMMAND", "value": "view"},
        {"name": "OUTPUT_FORMAT", "value": "bam"}
      ]
    }
  ]
}
```

## 🛠 Development

### Building Locally

```bash
# Build the Docker image
docker build -t bioinformatics-docker:latest .

# Test locally
docker-compose up
```

### Adding New Tools

1. Update the `Dockerfile` to install the new tool
2. Create a handler script in `scripts/[tool]-handler.sh`
3. Update the `entrypoint.sh` to route to the new handler
4. Add documentation and examples

## 🔄 CI/CD

This repository uses GitHub Actions to automatically build and push Docker images:

- **Push to `main`**: Builds and tags as `latest`
- **Git tags**: Builds and tags with version number
- **Pull requests**: Builds for testing only

## 📦 Image Tags

| Tag | Description | Built From |
|-----|-------------|------------|
| `latest` | Latest stable version | `main` branch |
| `v1.0.0` | Specific version | Git tag |
| `dev` | Development version | `develop` branch |

## 🧪 Testing

Local testing with sample data:

```bash
# Download test data
wget -O test-data/sample.sam "https://github.com/samtools/samtools/raw/develop/test/dat/test.sam"

# Run test
docker-compose -f docker-compose.test.yml up
```

## 📁 Repository Structure

```
bioinformatics-docker/
├── Dockerfile                 # Main Docker image
├── docker-compose.yml         # Local development
├── docker-compose.test.yml    # Testing setup
├── scripts/
│   ├── entrypoint.sh          # Main entry point (routes to handlers)
│   ├── samtools-handler.sh    # SAMtools operations (view, sort, index, etc.)
│   ├── zip-handler.sh         # Compression/decompression operations
│   ├── bwa-handler.sh         # BWA alignment
│   ├── bowtie2-handler.sh     # Bowtie2 alignment
│   ├── bedtools-handler.sh    # BEDTools operations
│   ├── fastqc-handler.sh      # FastQC quality control
│   ├── trimmomatic-handler.sh # Read trimming
│   ├── picard-handler.sh      # Picard tools
│   ├── vcftools-handler.sh    # VCF manipulation
│   └── [other tools...]       # Additional tool handlers
├── test-data/                 # Sample files for testing
├── .github/workflows/         # CI/CD pipelines
└── README.md                  # This file
```

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/new-tool`)
3. Add your bioinformatics tool with proper handler
4. Test locally with `docker-compose`
5. Submit a pull request

## 📝 License

This project is part of the smarts.bio platform. See the main repository for licensing information.

## 🆘 Support

For issues and questions:
- Create an issue in this repository
- Contact the smarts.bio development team
- Check the main documentation at [smarts.bio](https://smarts.bio)