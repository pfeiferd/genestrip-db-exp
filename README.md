**Genestrip-DB-Exp** - Experiments concerning [Genestrip](https://github.com/pfeiferd/genestrip)
===============================================
  
This project exists to run accuracy, quality and performance experiments with respect to [Genestrip]([Genestrip](https://github.com/pfeiferd/genestrip)).

## License

[Genestrip-DB-Exp is free for any kind of use.](./LICENSE.txt) 
However, the associated software, [Genestrip](https://github.com/pfeiferd/genestrip), has a [more restrictive License](https://github.com/pfeiferd/genestrip#license). 

## Building and installing

Genestrip-DB-Exp requires [Maven 2 or 3](https://maven.apache.org/) and is compatible with the JDK 11 or higher.

To build it, `cd` to the installation directory `genestrip-db-exp`. Given a matching Maven and JDK installation, `mvn install` will compile and install the Java program library.

The rest of the of system is a bunch of shell scripts that can be executed on an Ubuntu Linux x64 architecture.
You need a large disk with at least 2.5 TB free space and a high network bandwidth (with Internet access) to download and store the necessary data.
You need at least 128 GB of RAM for the result data to be computed. The computations may take several days...

This version of Genestrip-DB-Exp is tied to Genestrip [Version 2.4](https://github.com/pfeiferd/genestrip/releases/tag/v2.5).
The original experiments were based on the [RefSeq Release 233](https://ftp.ncbi.nlm.nih.gov/refseq/release/RELEASE_NUMBER).
Using later version of the RefSeq might bring (slightly) different results than the original runs.

## Experiment preparation

You need at least 2,5 TB of disk space (ideally on SSD drive) for the experiments.
You need excellent internet bandwidth in order to download over 1,5 TB of data.

Please first install [Ganon](https://pirovc.github.io/ganon/start/#install),
[ISS](https://insilicoseq.readthedocs.io/en/latest/iss/install.html) and 
[NanoSim](https://github.com/bcgsc/NanoSim?tab=readme-ov-file#installation) manually. 
(We were not able to support a compact installation script for these tools.)
Ganon should then be executable anywhere from the command line via `ganon`.
Similarly, ISS then be executable anywhere from the command line via `iss`.
Regarding NanoSim, we have used conda for its installation and
*always* run `conda activate nanosim` to enable it execution.
Moreover, `nanosimdir` should be an exported shell variable pointing to the `src` folder of NanoSim, i.e,
where its cors script `read_analysis.py` and `simulator.py` are stored.

Then, please `cd` to `genestrip-db-exp/bin` and execute the shell scripts there
in this given order:

1) `sh ./install_tools.sh` installs most necessary tools including KrakenUniq, Kraken 2, SRA tools and cgmemtime in respective folders. Please consult the [KrakenUniq README](https://github.com/fbreitwieser/krakenuniq/blob/master/README.md#installation) if this step partially fails in order to provide potential fixes.
2) `sh ./make_db.sh` prepares the *k*-mer databases necessary for the experiments. 
Beware: This incurs a download of all viral and bacterial genomes from the RefSeq, triggers the generation of three Genestrip databases, downloads large KrakenUniq and Kraken 2 databases and generates ganon databases.
3) `sh ./make_fastqs.sh` downloads or generates fastq files necessary for the experiments.

## Running the experiments (without performance experiments)

Please first `cd` to `genestrip-db-exp/bin`.

1) `sh ./classify.sh` performs classifications on viral data.
2) `sh ./classifyTicks.sh` performs classifications on tick data.
3) `sh ./runexps.sh` runs all the rest of non-performance-related experiments and produces related result files right under `genestrip-db-exp/results/`. 

## Performance experiments

### DB generation performance

Please first `cd` to `genestrip-db-exp/bin`. To run the database generation-experiments execute `sh ./run_gendb_perf_exps_linux.sh`.

Afterward, you will find log files like `db_gen_human_virus.log` etc. under `genestrip-db-exp/results/logs`.

### Classification performance

Please first `cd` to `genestrip-db-exp/bin`.

**To run the classification performance experiments with Genestrip:**

**Important:** `sh ./run_gendb_perf_exps_linux.sh` most have been
run beforehand so that the necessary databases exist.

1) Execute `sh ./run_match_ticks_perf_exps_linux.sh`.
2) Execute `sh ./run_match_saliva_perf_exps_linux.sh`.

Afterward, you will find log files like `match_tick-borne_tick1.log` etc. under `genestrip-db-exp/results/logs`.

**To run the classification performance experiments with KrakenUniq:**

1) Execute `sh ./run_ku_ticks_perf_exps_linux.sh`.
2) Execute `sh ./run_ku_saliva_perf_exps_linux.sh`.

Afterward, you will find log files like `match_ku_mb_tick1.log` etc. under `genestrip-db-exp/results/logs`.

**Note:** There are additional scripts, files and folders that are currently not needed to produce the experiments' results.
