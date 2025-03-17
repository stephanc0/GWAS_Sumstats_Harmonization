#!/bin/bash

# Define the installation directory inside the current working directory
INSTALL_DIR="$(pwd)/metal"

# Ensure the directory exists
mkdir -p "$INSTALL_DIR"

# METAL download URL
METAL_URL="https://csg.sph.umich.edu/abecasis/Metal/download/Linux-metal.tar.gz"

# Download METAL
echo "Downloading METAL..."
wget --no-check-certificate -O "$INSTALL_DIR/METAL.tar.gz" "$METAL_URL"

# Verify the download was successful
if [[ ! -f "$INSTALL_DIR/METAL.tar.gz" || $(wc -c <"$INSTALL_DIR/METAL.tar.gz") -lt 1000 ]]; then
    echo "Error: METAL download failed."
    rm -f "$INSTALL_DIR/METAL.tar.gz"
    exit 1
fi

# Extract METAL into the installation directory
echo "Extracting METAL..."
tar -xzf "$INSTALL_DIR/METAL.tar.gz" -C "$INSTALL_DIR" --strip-components=1

# Ensure the binary is executable
if [[ -f "$INSTALL_DIR/metal" ]]; then
    chmod +x "$INSTALL_DIR/metal"
    echo "METAL installed successfully in: $INSTALL_DIR"
else
    echo "Error: METAL binary not found after extraction. Exiting."
    exit 1
fi

# Cleanup
rm "$INSTALL_DIR/METAL.tar.gz"
