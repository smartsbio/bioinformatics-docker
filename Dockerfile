# Bioinformatics Docker Image with AWS Integration
# Starting with SAMtools - more tools to be added incrementally

FROM ubuntu:22.04

# Avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies and AWS CLI
RUN apt-get update && apt-get install -y \
    # Basic system tools
    curl \
    wget \
    unzip \
    jq \
    # Build tools for compiling SAMtools
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libcurl4-openssl-dev \
    # AWS CLI
    awscli \
    # Python for additional scripting
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

# Install SAMtools from source for latest version
ENV SAMTOOLS_VERSION=1.19.2
RUN cd /tmp && \
    wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    tar -xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VERSION} && \
    make && \
    make install && \
    cd / && \
    rm -rf /tmp/samtools-*

# Install htslib (dependency for SAMtools)
ENV HTSLIB_VERSION=1.19.1
RUN cd /tmp && \
    wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 && \
    tar -xjf htslib-${HTSLIB_VERSION}.tar.bz2 && \
    cd htslib-${HTSLIB_VERSION} && \
    make && \
    make install && \
    cd / && \
    rm -rf /tmp/htslib-*

# Create directories for our scripts and data
RUN mkdir -p /usr/local/bin/scripts /tmp/input /tmp/output

# Copy our custom scripts
COPY scripts/ /usr/local/bin/scripts/
RUN chmod +x /usr/local/bin/scripts/*.sh

# Set environment variables
ENV PATH="/usr/local/bin/scripts:${PATH}"

# Verify installations
RUN samtools --version && \
    aws --version && \
    echo "✅ SAMtools and AWS CLI installed successfully"

# Set working directory
WORKDIR /tmp

# Default entrypoint
ENTRYPOINT ["/usr/local/bin/scripts/entrypoint.sh"]