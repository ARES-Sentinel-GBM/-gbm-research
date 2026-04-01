# ARES-GBM Research Tools — Docker Image
# Containerized environment for GBM metabolic modeling

FROM python:3.11-slim

LABEL maintainer="ARES-Sentinel-GBM Team <admin@ares-bio.com>"
LABEL description="Docker image for ARES-GBM research pipeline"
LABEL version="1.0.0"

# Set working directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements first for better caching
COPY requirements.txt .

# Install Python dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Copy application files
COPY scripts/ ./scripts/
COPY data/ ./data/
COPY notebooks/ ./notebooks/

# Create results directory
RUN mkdir -p /app/results

# Set environment variables
ENV PYTHONUNBUFFERED=1
ENV PYTHONDONTWRITEBYTECODE=1

# Default command
CMD ["python", "scripts/flux_analysis.py", "--help"]
