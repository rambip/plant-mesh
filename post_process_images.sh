#!/bin/bash

# Check if ImageMagick is installed
if ! command -v convert &> /dev/null; then
    echo "Error: ImageMagick is not installed. Please install it first:"
    echo "sudo apt-get install imagemagick  # For Debian/Ubuntu"
    echo "sudo yum install ImageMagick      # For CentOS/RHEL"
    exit 1
fi

# Check if all required arguments are provided
if [ $# -ne 3 ]; then
    echo "Usage: $0 source_directory destination_directory brightness_value"
    echo "Example: $0 ~/images ~/adjusted 120"
    echo "         brightness: 100 = normal, >100 = brighter, <100 = darker"
    exit 1
fi

source_dir="$1"
dest_dir="$2"
brightness="$3"

# Create destination directory if it doesn't exist
mkdir -p "$dest_dir"

# Process each PNG file in the source directory
found_files=false
for img in "$source_dir"/*.png; do
    if [ -f "$img" ]; then
        found_files=true
        filename=$(basename "$img")
        echo "Processing: $filename"
        
        # Use modulate with all three parameters explicitly set
        convert "$img" \
            -fill white  -tint 120 \
            "$dest_dir/${filename%.*}.png"
    fi
done

if [ "$found_files" = false ]; then
    echo "No PNG files found in $source_dir"
    exit 1
fi

echo "Done! Processed images are in $dest_dir"
