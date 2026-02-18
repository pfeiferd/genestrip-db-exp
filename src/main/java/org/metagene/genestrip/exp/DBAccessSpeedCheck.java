package org.metagene.genestrip.exp;

import org.metagene.genestrip.bloom.AbstractKMerBloomFilter;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.store.KMerSortedArray;
import org.metagene.genestrip.tax.SmallTaxTree;

import java.io.File;
import java.io.IOException;
import java.util.Random;

public class DBAccessSpeedCheck extends GenestripComparator {
    private final Random rand = new Random(42);

    public DBAccessSpeedCheck(File baseDir) {
        super(baseDir, null);
    }

    public long[] testDBAccessSpeed(String db, int tries) throws IOException {
        Database database = getDatabase(db, false);
        KMerSortedArray<SmallTaxTree.SmallTaxIdNode> kMerSortedArray = database.convertKMerStore();
        AbstractKMerBloomFilter filter = kMerSortedArray.getFilter();
        long[] data = new long[tries];
        for (int i = 0; i < tries; i++) {
            data[i] = rand.nextLong();
        }
        long[] res = new long[3];
        long start = System.currentTimeMillis();
        for (int i = 0; i < tries; i++) {
            filter.containsLong(data[i]);
        }
        res[0] = System.currentTimeMillis() - start;
        kMerSortedArray.setUseFilter(false);
        start = System.currentTimeMillis();
        for (int i = 0; i < tries; i++) {
            kMerSortedArray.getLong(data[i], null);
        }
        res[1] = System.currentTimeMillis() - start;
        kMerSortedArray.setUseFilter(true);
        start = System.currentTimeMillis();
        for (int i = 0; i < tries; i++) {
            kMerSortedArray.getLong(data[i], null);
        }
        res[2] = System.currentTimeMillis() - start;

        return res;
    }


    public static void main(String[] args) throws IOException {
        DBAccessSpeedCheck check = new DBAccessSpeedCheck(new File("./data"));

        for (String db : new String[] { "viral", "human_virus", "tick-borne" }) {
            long[] res = check.testDBAccessSpeed(db, 1000000000);

            System.out.println(res[0]);
            System.out.println(res[1]);
            System.out.println(res[2]);

            System.out.println("+++ " + db + " +++");
            System.out.println("Filter versus L: " + res[1] / (double) res[0]);
            System.out.println("Filter + L versus L: " +  res[1] / (double) res[2]);
        }
    }
}
