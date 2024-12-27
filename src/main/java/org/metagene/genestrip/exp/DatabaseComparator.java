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
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class DatabaseComparator {
    private final File baseDir;

    public DatabaseComparator(File baseDir) {
        this.baseDir = baseDir;
    }

    public void reportRankCounts(String dbName) throws IOException {
        int[] counts = getRankCounts(dbName);

        for (Rank rank : Rank.values()) {
            System.out.println("Rank: " + rank + " Nodes: " + counts[rank.ordinal()]);
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
                if (taxNode != null) {
                    countsPerRank[taxNode.getRank().ordinal()]++;
                }
            }
        }
        return countsPerRank;
    }

    public void reportComparison(String db1Name, boolean temp, String db2Name)  throws IOException {
        File baseDir = null; // APITest.getBaseDir();

        KMerSortedArray<String> store1 = getDatabase(db1Name, temp).getKmerStore();
        KMerSortedArray<String> store2 = getDatabase(db2Name, false).getKmerStore();

        long movedKMers = compareDBs(store1, store2);

        // TODO: Write this to a file
        System.out.println("Database 1: " + db1Name);
        System.out.println("Database 2: " + db2Name);
        System.out.println("Total KMers: " + store1.getEntries());
        System.out.println("Moved KMers: " + movedKMers);
        System.out.println("Moved KMers %: " + (double)(movedKMers / store1.getEntries()) * 100);
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
                if (taxId1 != null || taxId2 != null) {
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

        GSProject project = new GSProject(config, dbName, null, null, null, null, null, false, "64320,12637+",
                null, null, null, false);

        GSMaker maker = new GSMaker(project);

        ObjectGoal<Database, GSProject> storeGoal = (ObjectGoal<Database, GSProject>) maker.getGoal(temp ? GSGoalKey.TEMPDB : GSGoalKey.LOAD_DB);
        return storeGoal.get();
    }

    public static void main(String[] args) throws Exception {
        DatabaseComparator comparator = new DatabaseComparator(new File("."));
        comparator.reportRankCounts(args[0]);
        comparator.reportComparison(args[0], false, args[1]);
    }
}
