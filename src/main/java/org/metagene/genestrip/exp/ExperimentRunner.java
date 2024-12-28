package org.metagene.genestrip.exp;

import java.io.File;
import java.io.IOException;

public class ExperimentRunner {
    public static final String[] DB_NAMES = new String[] { "chronicb" };

    private final File baseDir;
    private final String[] dbBaseNames;

    public ExperimentRunner(File baseDir, String[] dbBaseNames) {
        this.baseDir = baseDir;
        this.dbBaseNames = dbBaseNames;
    }

    public void runMatchComparison(String dbName, String csvFile) throws IOException {
        String dbName2 = getMinUpdatDBName(dbName);
        File ticksCSVFile = new File(baseDir, csvFile);
        new MatchComparator(baseDir).writeCompleteReport(dbName, dbName2, ticksCSVFile.getAbsolutePath());
    }

    public void runDatabaseComparisons() throws IOException {
        DatabaseComparator comparator = new DatabaseComparator(baseDir);
        for (String dbBaseName : dbBaseNames) {
            comparator.reportRankCounts(dbBaseName);
        }
        String[] minUpdatDBNames = new String[dbBaseNames.length];
        int i = 0;
        for (String dbBaseName : dbBaseNames) {
            minUpdatDBNames[i++] = getMinUpdatDBName(dbBaseName);
            i++;
        }
        comparator.reportComparisons(dbBaseNames, minUpdatDBNames);
    }

    protected String getMinUpdatDBName(String dbBaseName) {
        return dbBaseName + "-minupdate";
    }

    public static void main(String[] args) throws IOException {
        ExperimentRunner runner = new ExperimentRunner(new File("."), DB_NAMES);
        runner.runDatabaseComparisons();
        runner.runMatchComparison("chronicb", "./data/fastq/ticks.txt");
    }
}
