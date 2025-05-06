package org.metagene.genestrip.exp;

import org.metagene.genestrip.kucomp.KrakenDBComparator;
import org.metagene.genestrip.kucomp.KrakenMatchComparator;

import java.io.File;
import java.io.IOException;

public class ExperimentRunner {
    /*
    public static final String[] DB_NAMES = new String[]{"babesia",
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

    @Deprecated
    public ExperimentRunner(File baseDir, String[] offNames, String[] dbBaseNames) {
        this.baseDir = baseDir;
        this.dbBaseNames = dbBaseNames;
        this.offNames = offNames;
    }

    @Deprecated
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

    @Deprecated
    protected String getMinUpdateDBName(String dbBaseName) {
        return dbBaseName + "-minupdate";
    }
    */

    public static void main(String[] args) throws IOException {
        GenestripComparator c1 = new GenestripComparator(new File("./data"));
        KrakenMatchComparator c2 = new KrakenMatchComparator(new File("./data"));
        KrakenDBComparator c3 = new KrakenDBComparator(new File("./data"));

        c1.writeUnfoldedTaxids("human_virus");

        // Figure 2
        c3.reportKMerComparisons("viral", "viral_db");
        // Figure 3
        c2.compareWithKUResults("viral", "viral_ku_comp_fasta.txt", "viral_ku_comp.txt");
        // Text in context with Figure 3
        c2.accuracyCheckForSimulatedViralReads("viral", "viral_ku_comp.txt");
        // Figure 4
        c2.compareWithKUResults("viral", null, "saliva.txt");
        // Figure 5
        c1.compareCommonDBEntries("viral", "human_virus");
        // Use taxids from Genestrip's human virus database as basis for tax ids.
        c3.reportKrakenDBComparison("human_virus", "viral_db", "human_virus_db");
        // Figure 6
        c1.compareResults("viral", "human_virus", "saliva.txt");
        c2.compareKUWithKUResults("viral", "human_virus", "saliva.txt");

        // Not needed: ?
        //c2.compareWithKUResults("human_virus", null, "saliva.txt");
    }
}
