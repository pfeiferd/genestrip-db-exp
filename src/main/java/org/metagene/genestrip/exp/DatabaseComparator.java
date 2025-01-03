package org.metagene.genestrip.exp;

import it.unimi.dsi.fastutil.objects.Object2LongMap;
import org.metagene.genestrip.GSCommon;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSMaker;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.store.KMerSortedArray;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

public class DatabaseComparator {
    private static final DecimalFormat DF = new DecimalFormat("0.00", new DecimalFormatSymbols(Locale.US));

    private final File baseDir;

    public DatabaseComparator(File baseDir) {
        this.baseDir = baseDir;
    }

    public void reportRankCounts(String dbName) throws IOException {
        File file = new File(getOutDir(dbName),  dbName + "_ranks.csv");
        try (PrintStream out = new PrintStream(new FileOutputStream(file))) {
            int[] counts = getRankCounts(dbName);
            out.println("rank; nodes;");

            for (Rank rank : Rank.values()) {
                out.print(rank);
                out.print(";");
                out.print(counts[rank.ordinal()]);
                out.println(";");
            }
        }
    }

    public int[] getRankCounts(String dbName) throws IOException {
        Database database = getDatabase(dbName, false);
        Object2LongMap<String> stats = database.getStats();
        SmallTaxTree taxTree = database.getTaxTree();

        List<String> sortedTaxIds = new ArrayList<String>(stats.keySet());
        taxTree.sortTaxidsViaTree(sortedTaxIds);

        int[] countsPerRank = new int[Rank.values().length];

        for (String taxId : sortedTaxIds) {
            if (taxId != null) {
                SmallTaxTree.SmallTaxIdNode taxNode = taxTree.getNodeByTaxId(taxId);
                if (taxNode != null && taxNode.getRank() != null) {
                    countsPerRank[taxNode.getRank().ordinal()]++;
                }
            }
        }
        return countsPerRank;
    }

    public void reportComparisons(String[] offNames, String[] db1Name1, String[] db2Names) throws IOException {
        File file = new File(baseDir, "dbs_comp.csv");
        try (PrintStream out = new PrintStream(new FileOutputStream(file))) {
            out.println("name; db 1; temp; db 2; total kmers; moved kmers; moved kmers percent;");
            int i = 0;
            for (String dbName : db1Name1) {
                reportComparison(offNames[i], dbName, false, db2Names[i], out);
                reportComparison(offNames[i], dbName, true, db2Names[i], out);
                i++;
            }
        }
    }

    public void reportComparison(String offName, String db1Name, boolean temp, String db2Name, PrintStream out)  throws IOException {
        KMerSortedArray<String> store1 = getDatabase(db1Name, temp).getKmerStore();
        KMerSortedArray<String> store2 = getDatabase(db2Name, false).getKmerStore();

        long movedKMers = compareDBs(store1, store2);

        out.print(offName);
        out.print(";");
        out.print(db1Name);
        out.print(";");
        out.print(temp);
        out.print(";");
        out.print(db2Name);
        out.print(";");
        out.print(store1.getEntries());
        out.print(";");
        out.print(movedKMers);
        out.print(";");
        out.print(DF.format((100d * movedKMers) / store1.getEntries()));
        out.println(";");
    }

    public long compareDBs(final KMerSortedArray<String> store1, final KMerSortedArray<String> store2)  {
        if (store1.getEntries() != store2.getEntries()) {
            throw new RuntimeException("Store's entries do not match");
        }
        final long[] movedKmers = new long[] { 0 };
        store1.visit(new KMerSortedArray.KMerSortedArrayVisitor<String>() {
            @Override
            public void nextValue(KMerSortedArray<String> trie, long kmer, short index, long i) {
                String taxId1 = store1.getValueForIndex(index);
                String taxId2 = store2.getLong(kmer, null);
                if (taxId1 == null || taxId2 == null) {
                    throw new RuntimeException("Inconsistent databases");
                }

                if (!taxId1.equals(taxId2)) {
                    movedKmers[0]++;
                }
            }
        });
        return movedKmers[0];
    }

    protected Database getDatabase(String dbName, boolean temp) throws IOException {
        GSCommon config = new GSCommon(baseDir);

        GSProject project = new GSProject(config, dbName, null, null, null, null, null, false, null,
                null, null, null, false);

        GSMaker maker = new GSMaker(project);

        ObjectGoal<Database, GSProject> storeGoal = (ObjectGoal<Database, GSProject>) maker.getGoal(temp ? GSGoalKey.LOAD_TEMPDB : GSGoalKey.LOAD_DB);
        return storeGoal.get();
    }

    protected File getOutDir(String dbName) throws IOException {
        GSCommon config = new GSCommon(baseDir);
        GSProject project = new GSProject(config, dbName, null, null, null, null, null, false, null,
                null, null, null, false);
        return project.getResultsDir();
    }
}
