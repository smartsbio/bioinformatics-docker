# Bioinformatics Docker Images

Custom Docker images for bioinformatics tools with AWS integration, designed for the smarts.bio platform.

## 🎯 Overview

This repository contains Docker images that combine popular bioinformatics tools with AWS CLI and cloud-native processing capabilities. Each image is built for scalable, containerized bioinformatics workflows.

## 🧬 Currently Supported Tools

### SAMtools
- **Version**: 1.19.2
- **Supported Commands**: 
  - `view` - Format conversion (SAM ↔ BAM ↔ CRAM)
  - `sort` - Sort alignment files
  - `index` - Create index files
  - `stats` - Generate alignment statistics

## 🚀 Usage

### Environment Variables

| Variable | Required | Description | Example |
|----------|----------|-------------|---------|
| `TOOL_NAME` | ✅ | Bioinformatics tool to use | `samtools` |
| `COMMAND` | ✅ | Tool-specific command | `view` |
| `S3_BUCKET` | ✅ | S3 bucket for input/output files | `my-bucket` |
| `INPUT_S3_KEY` | ✅ | S3 key for input file | `data/input.sam` |
| `OUTPUT_PATH` | ✅ | S3 path for output files | `results/` |
| `OUTPUT_FORMAT` | ❌ | Output format (for applicable tools) | `bam` |
| `OUTPUT_FILE` | ❌ | Output filename | `output.bam` |

### Docker Run Example

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
│   ├── entrypoint.sh          # Main entry point
│   ├── samtools-handler.sh    # SAMtools operations
│   └── [future tools...]      # Additional tool handlers
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