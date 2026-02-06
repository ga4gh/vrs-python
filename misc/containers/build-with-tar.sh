#!/bin/bash

# Build script that creates a tar.gz with only necessary files for container build
# Usage: ./build-with-tar.sh [ASSEMBLY]

set -e

ASSEMBLY=${1:-GRCh38}
TAR_NAME="build-context.tar.gz"
BUILD_DIR="build-context"

echo "Building container with assembly: $ASSEMBLY"
echo "Creating build context tar.gz..."

# Clean up any existing build context
rm -rf "$BUILD_DIR" "$TAR_NAME"

# Create build directory
mkdir -p "$BUILD_DIR/misc/containers"

# Copy necessary files for the container build
echo "Copying files to build context..."

# Container-specific files
cp misc/containers/Dockerfile "$BUILD_DIR/misc/containers/"
cp misc/containers/entrypoint.sh "$BUILD_DIR/misc/containers/"
cp misc/containers/build-${ASSEMBLY}.bash "$BUILD_DIR/misc/containers/"

# Create the tar.gz
echo "Creating tar.gz..."
rm -rf "$TAR_NAME"
tar -czf "$TAR_NAME" -C "$BUILD_DIR" .

# Clean up build directory
rm -rf "$BUILD_DIR"

echo "Build context created: $TAR_NAME"

# Detect container runtime
if command -v docker >/dev/null 2>&1; then
    CONTAINER_CMD="docker"
    echo "Using Docker for build..."
elif command -v podman >/dev/null 2>&1; then
    CONTAINER_CMD="podman"
    echo "Using Podman for build..."
else
    echo "Error: Neither docker nor podman found in PATH"
    exit 1
fi

# Run container build with the tar.gz as context
cat "$TAR_NAME" | $CONTAINER_CMD build \
    --arch linux/arm64,linux/amd64 \
    --build-arg ASSEMBLY="$ASSEMBLY" \
    --target data \
    -t ghcr.io/theferrit32/vrs-python:${ASSEMBLY}-data \
    -f ./misc/containers/Dockerfile

cat "$TAR_NAME" | $CONTAINER_CMD build \
    --arch linux/arm64,linux/amd64 \
    --build-arg ASSEMBLY="$ASSEMBLY" \
    --target build \
    -t ghcr.io/theferrit32/vrs-python:${ASSEMBLY}-build \
    -f ./misc/containers/Dockerfile

cat "$TAR_NAME" | $CONTAINER_CMD build \
    --arch linux/arm64,linux/amd64 \
    --build-arg ASSEMBLY="$ASSEMBLY" \
    --target build \
    -t ghcr.io/theferrit32/vrs-python:${ASSEMBLY} \
    -f ./misc/containers/Dockerfile

# Clean up tar file
# rm -f "$TAR_NAME"

echo "Build completed successfully!"
