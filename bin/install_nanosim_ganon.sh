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
# Therefore just on Linux:
conda install -c bioconda -c conda-forge ganon
