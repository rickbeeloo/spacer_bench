#!/bin/bash


# Function to create environment and install tool
create_env() {
    local env_name=$1
    local package=$2
    local channel=${3:-bioconda}  # default to bioconda if no channel specified

    echo "Creating environment: $env_name with package: $package"
    
    # Remove environment if it exists
    micromamba env remove -n $env_name # ignore errors if it doesn't exist yolo
    
    # Create new environment with only the specific tool
    micromamba create -n $env_name -c $channel $package -y
}


tools=("bowtie1" "bowtie2" "minimap2" "bbmap" "strobealign" "blast" "mmseqs" "spacepharer" "spacer-containment" "mummer4", "lexicmap","vsearch", "bwa", "hisat2", "sassy")

# Create individual environments for each tool
for tool in "${tools[@]}"; do
    if [ "$tool" != "sassy" ]; then
        create_env "${tool}_env" "$tool"
    fi
done

# Create individual environments for each tool
# create_env "bowtie1_env" "bowtie"
# create_env "bowtie2_env" "bowtie2"
# create_env "minimap2_env" "minimap2"
# create_env "bbmap_env" "bbmap"
# create_env "strobealign_env" "strobealign"
# create_env "blast_env" "blast"
# create_env "mmseqs_env" "mmseqs2"
# create_env "spacepharers_env" "spacepharer"



# For spacer-containment, assuming it's available through a specific channel
# Adjust the channel and package name as needed
create_env "spacer_containment_env" "spacer-containment" 

echo "All environments created successfully!"

# function to check if a tool is installed and print the version
check_tool_version() {
    local tool_name=$1
    echo "#####" >> versions.txt
    micromamba run -n $tool_name"_env" $tool_name --version  >> versions.txt
}   

# Print versions of installed tools
echo "Checking installed versions:"
echo "#####" > versions.txt
# check_tool_version "bowtie"
# check_tool_version "bowtie2"
# check_tool_version "minimap2"
# check_tool_version "bbmap"
# check_tool_version "strobealign"
# check_tool_version "blast"
# check_tool_version "mmseqs"

micromamba run -n bowtie2_env bowtie2 --version >> versions.txt
echo "#####" >> versions.txt
micromamba run -n minimap2_env minimap2 --version >> versions.txt
echo "#####" >> versions.txt
micromamba run -n bbmap_env bbmap.sh --version >> versions.txt
echo "#####" >> versions.txt
micromamba run -n strobealign_env strobealign --version >> versions.txt
echo "#####" >> versions.txt
micromamba run -n blast_env blastn -version >> versions.txt
echo "#####" >> versions.txt
micromamba run -n mmseqs_env mmseqs version >> versions.txt
echo "#####" >> versions.txt
micromamba run -n spacepharer_env spacepharer --version >> versions.txt
echo "#####" >> versions.txt


###### 
# mmseqs 16.747c6 (current bioconda version) results in frequent crashes. as it is borken, we will install the latest version from github
micromamba activate mmseqs_env
wget https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz 
tar -xvzf mmseqs-linux-avx2.tar.gz
mv ./mmseqs/bin/mmseqs $CONDA_PREFIX/bin/mmseqs
chmod +x $CONDA_PREFIX/bin/mmseqs
rm -rf mmseqs-linux-avx2.tar.gz mmseqs
echo "#####\n MMseqs version: \n" >> tool_versions.txt
mmseqs version >> tool_versions.txt

###### 
# mummer 4 from main branch has bugs in SAM format output (see https://github.com/mummer4/mummer/issues/24)
# we will install the latest version of the develop branch (https://github.com/mummer4/mummer/tree/develop)
micromamba activate mummer4_env
# install dependencies
micromamba install gcc make yaggo -y
wget https://github.com/mummer4/mummer/archive/refs/heads/develop.zip # https://github.com/mummer4/mummer/releases/download/v4.0.1/mummer-4.0.1.tar.gz testing now
rm -rf mummer-develop
unzip develop.zip
cd mummer-develop
rm $CONDA_PREFIX/bin/mummer -rf
autoreconf -fi
./configure --prefix=$CONDA_PREFIX/bin/mummer
make  -j 6
make install
cd ..
rm -rf mummer-develop develop.zip
cp $CONDA_PREFIX/bin/mummer/bin/nucmer $CONDA_PREFIX/bin/nucmer
chmod +x $CONDA_PREFIX/bin/nucmer
echo "#####\n Mummer version: \n" >> tool_versions.txt
nucmer --version >> tool_versions.txt


######
# Sassy - Rust tool that needs to be built from source
micromamba activate sassy_env
# Install Rust if not already installed
if ! command -v cargo &> /dev/null; then
    echo "Installing Rust..."
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
    source ~/.cargo/env
fi

# Clone and build Sassy
echo "Building Sassy from source..."
git clone https://github.com/RagnarGrootKoerkamp/sassy.git
cd sassy
cargo build --release
# Copy the binary to the conda environment
cp target/release/sassy $CONDA_PREFIX/bin/sassy
chmod +x $CONDA_PREFIX/bin/sassy
cd ..
rm -rf sassy
echo "#####\n Sassy version: \n" >> tool_versions.txt
sassy --version >> tool_versions.txt 2>&1 || echo "Version info not available" >> tool_versions.txt

######
# hisat2
# create_env "hisat2_env" "hisat2"
micromamba activate hisat2_env
echo "#####\n Hisat2 version: \n" >> tool_versions.txt
hisat2 --version >> tool_versions.txt
# wget https://github.com/DaehwanKimLab/hisat2/archive/refs/tags/v2.2.1.tar.gz
# tar -xvzf v2.2.1.tar.gz
# cd hisat2-2.2.1
# ./install_hisat2.sh
# cd ..
# rm -rf v2.2.1.tar.gz hisat2-2.2.1

######
# Bwa   
# create_env "bwa_env" "bwa"
micromamba activate bwa_env
echo "#####\n Bwa version: \n" >> tool_versions.txt
bwa 2>&1 |  head -n 3 |tail -n 1 >> tool_versions.txt

######