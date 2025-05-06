#!/bin/bash


# (using this before uploading to zenodo)

# Check if required tools are available
for cmd in samtools zstd; do
    if ! command -v $cmd &> /dev/null; then
        echo "Error: $cmd is required but not installed. Please install it first."
        exit 1
    fi
done

# Function to compress SAM files
compress_sam() {
    local file="$1"
    if [[ -f "$file" && ! "$file" =~ \.gz$ ]]; then
        echo "Compressing $file..."
        bgzip -@ 10 -f "$file"
    fi
}

# Function to compress TSV files
compress_tsv() {
    local file="$1"
    if [[ -f "$file" && ! "$file" =~ \.zst$ ]]; then
        echo "Compressing $file..."
        zstd -f -T0 "$file"
    fi
}

# Function to compress FASTA files
compress_fasta() {
    local file="$1"
    if [[ -f "$file" && ! "$file" =~ \.gz$ ]]; then
        echo "Compressing $file..."
        bgzip -@ 10 -f "$file"
    fi
}

# Process simulated data
echo "Processing simulated data..."
cd results/simulated

# Compress FASTA file s in simulated_data directories
for fastafile in $(find . -type f -name "simulated_*.fa"); do
    compress_fasta "$fastafile"
done

# Compress SAM files in all run_parameters directories
for samfile in $(find . -type f -name "*output.sam"); do
    compress_sam "$samfile"
done

# Compress TSV files in all run_parameters directories
for tsvfile in $(find . -type f -name "*output.tsv" -o -name "*_results.tsv"); do
    compress_tsv "$tsvfile"
done

# Process real data
echo "Processing real data..."
cd ../real_data

# Compress SAM files in all run_parameters directories
for samfile in $(find . -type f -name "*output.sam"); do
    compress_sam "$samfile"
done

# Compress TSV files in all run_parameters directories
for tsvfile in $(find . -type f -name "*output.tsv" -o -name "*_results.tsv"); do
    compress_tsv "$tsvfile"
done

echo "Compression complete!"

# Print summary of compressed files
echo -e "\nCompression Summary:"
echo "FASTA files compressed: $(find . -type f -name "*.fa.gz" | wc -l)"
echo "SAM files compressed: $(find . -type f -name "*.sam.gz" | wc -l)"
echo "TSV files compressed: $(find . -type f -name "*.tsv.zst" | wc -l)" 



# Function to create tar archive
create_archive() {
    local dir="$1"
    local archive_name="$2"
    
    echo "Creating archive for $dir..."
    cd "$dir"
    tar --use-compress-program=pigz -cf "$archive_name.tar.gz"  .
    
    # Print archive size
    echo "Archive size: $(du -h "$archive_name.tar.gz" | cut -f1)"
    cd - > /dev/null
}

# Create archives for both directories
create_archive "simulated" "simulated_data"
create_archive "real_data" "real_data"

echo "Archive creation complete!" 