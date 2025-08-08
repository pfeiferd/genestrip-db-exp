**Genestrip-DB-Exp** - Experiments concerning [Genestrip](https://github.com/pfeiferd/genestrip)
===============================================
  
This project exists to run accuracy, quality and performance experiments with respect to [Genestrip]([Genestrip](https://github.com/pfeiferd/genestrip)).

## License

[Genestrip-DB-Exp is free for any kind of use.](./LICENSE.txt) 
However, the associated software, [Genestrip](https://github.com/pfeiferd/genestrip), has a [more restrictive License](https://github.com/pfeiferd/genestrip#license). 

## Building and installing

Genestrip-DB-Exp requires [Maven 2 or 3](https://maven.apache.org/) and is compatible with the JDK 11 or higher.

To build it, `cd` to the installation directory `genestrip-db-exp`. Given a matching Maven and JDK installation, `mvn install` will compile and install the Java program library.

The rest of the of system is a bunch of shell scripts that can be executed on a Linux x64 or on a Mac Arm architecture.
You need a large disk with at least 2 TB free space and a high network bandwidth (with Internet access) to download and store the necessary data.
You need at least 32 GB of RAM for the result data to be computed. The computations may take several days...

This version of Genestrip-DB-Exp is tied to Genestrip [Version 2.0](https://github.com/pfeiferd/genestrip/releases/tag/v2.0).
The original experiments were based on the [RefSeq Release 230](https://ftp.ncbi.nlm.nih.gov/refseq/release/RELEASE_NUMBER).
Using later version of the RefSeq might bring (slightly) different results than the original runs.

## Experiment preparation

You need at least 2,5 TB of disk space (ideally on SSD drive) for the experiments.
You need excellent internet bandwidth in order to download over 1,5 TB of data.

Please `cd` to `genestrip-db-exp/bin` and execute the shell scripts there
in this given order:

1) `sh ./install_krakenuniq.sh` installs KrakenUniq in a respective folder. Please consult the [KrakenUniq README](https://github.com/fbreitwieser/krakenuniq/blob/master/README.md#installation) if this step fails in order to provide potential fixes.
2) `sh ./install_sra_sdk_linux64.sh` or `sh ./install_sra_sdk_mac_arm.sh` installs the [SRA tools](https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit) on Linux or a Mac, respectively.
3) `sh ./download_hs_fastqs.sh` downloads fastq files based on human saliva samples from the [SRA](https://www.ncbi.nlm.nih.gov/sra/).
4) `sh ./download_tick_fastqs.sh` downloads fastq files related to tick analysis from the [SRA](https://www.ncbi.nlm.nih.gov/sra/).
5) `sh ./download_ku_microbial_db.sh` downloads and installs the [MicrobialDB](https://benlangmead.github.io/aws-indexes/k2) for KrakenUniq.
6) `sh ./make_ku_viraldb.sh` creates a complete viral DB for KrakenUniq. This incurs a download of viral genomes from the RefSeq.
7) `sh ./make_ku_human_virusdb.sh` creates a human virus DB for KrakenUniq. This incurs a download of viral genomes from the RefSeq.

## Running the experiments (without performance experiments)

`sh ./runexps.sh` runs all the non-performance-related experiments and produces related result files right under `genestrip-db-exp/results/`. 
Beware: This incurs a download of all viral and bacterial genomes from the RefSeq and triggers the generation of three Genestrip databases.

## Performance experiments

### DB generation performance

Please first `cd` to `genestrip-db-exp/bin`.

**Important:** For the script under 2. you must use `bash` under Linux.

1) Under *Linux*, please install [`cgmemtime`](https://github.com/gsauthof/cgmemtime) which should be used
instead of `/usr/bin/time`. To do so, execute `sh ./install_cgmemtime_linux.sh`. (However, `/usr/bin/time` is used under *macOS*.)
2) To run the database generation-experiments, execute
`sh ./run_gendb_perf_exps_macos.sh` or `bash ./run_gendb_perf_exps_linux.sh`, respectively.

Afterward, you will find log files like `db_gen_human_virus.log` etc. under `genestrip-db-exp/results/logs`.

### Classification performance

Please first `cd` to `genestrip-db-exp/bin`.

**To run the classification performance experiments with Genestrip:**

**Important:** `sh ./run_gendb_perf_exps_macos.sh` or `bash ./run_gendb_perf_exps_linux.sh` most have been
run beforehand so that the necessary databases exist.
Under Linux, `cgmemtime` must have been installed via `sh ./install_cgmemtime_linux.sh` (see above).
For the following scripts you must use `bash` under Linux.

1) Execute `sh ./run_match_ticks_perf_exps_macos.sh` or `bash ./run_match_ticks_perf_exps_linux.sh`.
2) Execute `sh ./run_match_saliva_perf_exps_macos.sh` or `bash ./run_match_saliva_perf_exps_linux.sh`.

Afterward, you will find log files like `match_tick-borne_tick1.log` etc. under `genestrip-db-exp/results/logs`.

**To run the classification performance experiments with KrakenUniq:**

**Important:** KrakenUniq does not support multithreading under macOS,
so for fair results, corresponding scripts should be run under Linux (although there exist corresponding scripts for macOS).
For the following scripts you must use `bash` under Linux.

1) Execute `sh ./run_ku_ticks_perf_exps_macos.sh` or `bash ./run_ku_ticks_perf_exps_linux.sh`.
2) Execute `sh ./run_ku_saliva_perf_exps_macos.sh` or `bash ./run_ku_saliva_perf_exps_linux.sh`.

Afterward, you will find log files like `match_ku_mb_tick1.log` etc. under `genestrip-db-exp/results/logs`.

**Note:** There are additional scripts, files and folders that are currently not used to produce the experiments' results.
