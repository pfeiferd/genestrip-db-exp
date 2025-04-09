package org.metagene.genestrip.exp;

import org.metagene.genestrip.kucomp.KrakenMatchComparator;

import java.io.File;
import java.io.IOException;

public class ExperimentRunner {
    public static final String[] DB_NAMES = new String[] { "babesia",
            "borrelia",
            "borrelia_plasmid",
            "chronicb",
            "human_virus2",
    "parasites",
    "plasmopara",
    "protozoa",
    "vineyard"};
    public static final String[] OFF_NAMES = DB_NAMES;

    private final File baseDir;
    private final String[] offNames;
    private final String[] dbBaseNames;

    public ExperimentRunner(File baseDir, String[] offNames, String[] dbBaseNames) {
        this.baseDir = baseDir;
        this.dbBaseNames = dbBaseNames;
        this.offNames = offNames;
    }

    public void runMatchComparison(String dbName, String csvFile) throws IOException {
        String dbName2 = getMinUpdateDBName(dbName);
        File ticksCSVFile = new File(baseDir, csvFile);
        new MatchComparator(baseDir).writeCompleteReport(dbName2, dbName, ticksCSVFile.getAbsolutePath());
    }

    public void runDatabaseComparisons() throws IOException {
        DatabaseComparator comparator = new DatabaseComparator(baseDir);
        for (String dbBaseName : dbBaseNames) {
            comparator.reportRankCounts(dbBaseName);
        }
        String[] minUpdateDBNames = new String[dbBaseNames.length];
        int i = 0;
        for (String dbBaseName : dbBaseNames) {
            minUpdateDBNames[i++] = getMinUpdateDBName(dbBaseName);
        }
        comparator.reportComparisons(offNames, dbBaseNames, minUpdateDBNames);
    }

    protected String getMinUpdateDBName(String dbBaseName) {
        return dbBaseName + "-minupdate";
    }

    public static void main(String[] args) throws IOException {
        /*
        ExperimentRunner runner = new ExperimentRunner(new File("./data"), OFF_NAMES, DB_NAMES);
        runner.runDatabaseComparisons();
        runner.runMatchComparison("chronicb", "fastq/ticks.txt");
         */
        KrakenMatchComparator c = new KrakenMatchComparator(new File("./data"));
        c.compareWithKUResults("viral", "viral_ku_comp_fasta.txt", "viral_ku_comp.txt", new String[] { "fastq1" });
    }
}
