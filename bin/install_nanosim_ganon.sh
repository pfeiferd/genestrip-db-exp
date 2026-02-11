# Homebrew
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
# Rest of Homebrew
echo >>  ~/.bashrc
echo 'eval "$(/home/linuxbrew/.linuxbrew/bin/brew shellenv bash)"' >>  ~/.bashrc
eval "$(/home/linuxbrew/.linuxbrew/bin/brew shellenv bash)"
# Miniconda - works on Mac - but on Linux?
brew install miniconda


## THIS WORKED for me on LINUX

# Conda
# On Linux, this worked for me
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh

### IMPORTANT EXIT SHELL now and reconnect - so that changes take effect. ###

cd  ~

# Bioconda
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
# Nanosim
git clone https://github.com/bcgsc/NanoSim.git
conda create --name nanosim python=3.7
conda activate nanosim
cd NanoSim
conda install --file requirements.txt -c conda-forge -c bioconda

# Get additional regex dependency (this worked for me using "pip")
pip install regex

cd  ~/NanoSim/src
# Execute a test command:
./read_analysis.py

# DON'T FORGET: Every time you login again, you have to reset the conda environment:
conda activate nanosim


# Ganon:
# Doesn't seem to work on a MAC (I can't compile it, I can't install it via conda).
# Therefor just on Linux:
conda install -c bioconda -c conda-forge ganon

# Requires FTP to work.
ganon build --db-prefix viral_cg_rs --source refseq --organism-group viral --complete-genomes --threads 24


# Run ganon on simulated virus fastq from paper:
ganon classify --db-prefix viral_cg_rs data/projects/viral/fastq/viral_fasta2fastq_fasta1.fastq.gz --output-all results --threads 32
