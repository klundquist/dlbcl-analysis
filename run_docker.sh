#!/bin/bash
# Script to build and run the Docker container

# Set the image name
IMAGE_NAME="dlbcl-analysis"
CONTAINER_NAME="dlbcl-analysis"

# Check if Docker is installed
if ! command -v docker &> /dev/null; then
    echo "Error: Docker is not installed or not in your PATH"
    echo "Please install Docker from https://docs.docker.com/get-docker/"
    exit 1
fi

# Build the Docker image (this will take some time the first run)
echo "Building Docker image with all required packages..."
docker build -t $IMAGE_NAME .

# Remove any old containers
echo "Cleaning up any existing containers..."
docker rm -f $CONTAINER_NAME 2>/dev/null || true

# Run the container
echo "Starting Docker container..."
echo "
Welcome to the DLBCL Analysis Container!

This container has all required packages pre-installed.
To run the analysis pipeline, type:
> source(\"scripts/01_get_process_data.R\")
> run_pipeline()
"

docker run -it --name $CONTAINER_NAME \
    -v "$(pwd):/project" \
    $IMAGE_NAME --quiet --no-save

# Note: If the container exits, you can restart it with:
# docker start -i dlbcl-analysis