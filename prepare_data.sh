#!/bin/bash

DATA_PATH="data"

# Install unzip if not present
if ! command -v unzip &> /dev/null
then
    echo "unzip is not installed. Installing now..."
    sudo apt-get update
    sudo apt-get install -y unzip
fi

# Check if the data is a zip file
if [ -f "$DATA_PATH.zip" ]; then
  echo "Unzipping data.zip..."
  unzip -o "$DATA_PATH.zip" -d "$DATA_PATH"
else
  echo "No zip file found at $DATA_PATH/data.zip"
fi

# Check if the data directory exists
if [ ! -d "$DATA_PATH" ]; then
  echo "Error: The directory '$DATA_PATH' does not exist."
  exit 1
fi

# Go trough all subdirectories and decompress .gz files
cd "$DATA_PATH" || exit

echo "Decompressing all .gz files in subdirectories..."
find . -type f -name "*.gz" -print0 | while IFS= read -r -d $'\0' file
do
    gunzip "$file"
done

echo "All .gz files have been decompressed."
