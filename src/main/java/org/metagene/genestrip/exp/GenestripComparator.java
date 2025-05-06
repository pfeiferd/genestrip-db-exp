package org.metagene.genestrip.exp;

import it.unimi.dsi.fastutil.objects.Object2LongMap;
import org.metagene.genestrip.*;
import org.metagene.genestrip.goals.MatchResultGoal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.match.CountsPerTaxid;
import org.metagene.genestrip.match.MatchingResult;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.TaxTree;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

public class GenestripComparator {
    public static final String SPECIES_OR_BELOW = "species or below";
    public static final String GENUS = "genus";
    public static final String ABOVE_GENUS = "null";

    protected final File baseDir;

    public GenestripComparator(File baseDir) {
        this.baseDir = baseDir;
    }

    public void compareCommonDBEntries(String dbName1, String dbName2) throws IOException {
        Database db1 = getDatabase(dbName1, false);
        SmallTaxTree tree1 = db1.getTaxTree();
        Object2LongMap<String> stats1 = db1.getStats();
        db1 = null;

        Database db2 = getDatabase(dbName2, false);
        Object2LongMap<String> stats2 = db2.getStats();
        SmallTaxTree tree2 = db2.getTaxTree();
        db2 = null;

        File out = new File(baseDir, dbName1 + "_" + dbName2 + "_gs_gs_dbcomp.csv");
        try (PrintStream ps = new PrintStream(new FileOutputStream(out))) {
            ps.println("taxid; rank; kmers 1; kmers 2;");
            for (SmallTaxTree.SmallTaxIdNode node1 : tree1) {
                String taxId = node1.getTaxId();
                SmallTaxTree.SmallTaxIdNode node2 = tree2.getNodeByTaxId(taxId);
                // We only report on tax ids which are in both (Genestrip) dbs:
                if (node1 != null && node2 != null) {
                    long count1 = stats1.getOrDefault(taxId, 0);
                    long count2 = stats2.getOrDefault(taxId, 0);
                    ps.print(taxId);
                    ps.print(";");
                    ps.print(getRankString(node1));
                    ps.print(";");
                    ps.print(correctDBValue(count1));
                    ps.print(";");
                    ps.print(correctDBValue(count2));
                    ps.println(";");
                }
            }
        }
    }

    protected String getRankString(SmallTaxTree.SmallTaxIdNode node) {
        boolean norank = false;
        while (node != null && Rank.NO_RANK.equals(node.getRank())) {
            norank = true;
            node = node.getParent();
        }
        Rank r = node == null ? null : node.getRank();
        if (r != null && (r.isBelow(Rank.SPECIES) ||
                (norank && r.equals(Rank.SPECIES)))) {
            /*
            if (count1 != count2) {
                System.out.println("Diff below species for: " + key);
            }
            */
            return SPECIES_OR_BELOW;
        } else if (r != null && r.equals(Rank.SPECIES)) {
            /*
            if (count1 != count2) {
                System.out.println("Diff. species for: " + key);
            }
             */
            return SPECIES_OR_BELOW;
        } else if (r != null && r.equals(Rank.GENUS)) {
           return GENUS;
        } else {
           return ABOVE_GENUS;
        }
    }

    protected long correctDBValue(long v) {
        return v + 1;
    }

    protected void writeUnfoldedTaxids(String dbName) throws IOException {
        GSCommon config = new GSCommon(baseDir);
        GSProject project = new GSProject(config, dbName, null, null, null, null, null, null,
                null, null, null, false);
        GSMaker maker = new GSMaker(project);
        Set<TaxTree.TaxIdNode> nodes = ((ObjectGoal<Set<TaxTree.TaxIdNode>, GSProject>) maker.getGoal(GSGoalKey.TAXNODES)).get();
        File out = new File(project.getProjectDir(), "unfolded_taxids.txt");
        try (PrintStream ps = new PrintStream(new FileOutputStream(out))) {
            for (TaxTree.TaxIdNode node : nodes) {
                ps.println(node.getTaxId());
            }
        }
    }

    protected Database getDatabase(String dbName, boolean temp) throws IOException {
        GSCommon config = new GSCommon(baseDir);
        GSProject project = new GSProject(config, dbName, null, null, null, null, null, null,
                null, null, null, false);
        GSMaker maker = new GSMaker(project);
        ObjectGoal<Database, GSProject> storeGoal = (ObjectGoal<Database, GSProject>) maker.getGoal(temp ? GSGoalKey.LOAD_TEMPDB : GSGoalKey.LOAD_DB);
        maker.getGoal(GSGoalKey.DBINFO).make(); // Ensure DB info is around too...

        Database db = storeGoal.get();
        maker.dumpAll();
        return db;
    }

    private int errs;
    private long kMersErrSum;
    private long kMersErrSquareSum;
    private long readsErrSum;
    private long readsErrSquareSum;

    public void compareResults(String dbName1, String dbName2, String csvFile) throws IOException {
        errs = 0;
        kMersErrSum = 0;
        kMersErrSquareSum = 0;
        readsErrSum = 0;
        readsErrSquareSum = 0;

        SmallTaxTree[] taxTreeRef1 = new SmallTaxTree[1];
        Map<String, MatchingResult> matches1 = match(dbName1, csvFile, taxTreeRef1);
        SmallTaxTree[] taxTreeRef2 = new SmallTaxTree[1];
        Map<String, MatchingResult> matches2 = match(dbName2, csvFile, taxTreeRef2);

        Map<String, Map<String, long[]>> allResults = new HashMap<>();

        for (String key : matches1.keySet()) {
            MatchingResult res1 = matches1.get(key);
            MatchingResult res2 = matches2.get(key);
            if (res1.getGlobalStats().getKMers() != res2.getGlobalStats().getKMers()) {
                throw new RuntimeException("Match results do not match");
            }

            Map<String, CountsPerTaxid> stats1 = res1.getTaxid2Stats();
            Map<String, CountsPerTaxid> stats2 = res2.getTaxid2Stats();

            File out = new File(baseDir, dbName1 + "_" + dbName2 + "_" + key + "_gs_gs_comp.csv");
            try (PrintStream ps = new PrintStream(new FileOutputStream(out))) {
                ps.println("taxid; rank; kmers 1; kmers 2; ukmers 1; ukmers 2; reads 1; reads 2");
                for (SmallTaxTree.SmallTaxIdNode node1 : taxTreeRef1[0]) {
                    String taxId = node1.getTaxId();
                    // We only report on tax ids which are in both dbs, so "taxTreeRef2":
                    SmallTaxTree.SmallTaxIdNode node2 = taxTreeRef2[0].getNodeByTaxId(taxId);
                    if (node2 == null) {
                        continue;
                    }
                    CountsPerTaxid c1 = stats1.get(taxId);
                    CountsPerTaxid c2 = stats2.get(taxId);
                    ps.print(taxId);
                    ps.print(';');
                    String rs = getRankString(node1);
                    ps.print(rs);
                    ps.print(';');
                    ps.print(correctDBValue(c1 == null ? 0 : c1.getKMers()));
                    ps.print(';');
                    ps.print(correctDBValue(c2 == null ? 0 : c2.getKMers()));
                    ps.print(';');
                    ps.print(correctDBValue(c1 == null ? 0 : c1.getUniqueKMers()));
                    ps.print(';');
                    ps.print(correctDBValue(c2 == null ? 0 : c2.getUniqueKMers()));
                    ps.print(';');
                    ps.print(correctDBValue(c1 == null ? 0 : c1.getReads()));
                    ps.print(';');
                    ps.print(correctDBValue(c2 == null ? 0 : c2.getReads()));
                    ps.println(';');
                    if (SPECIES_OR_BELOW.equals(rs)) {
                        sumErrorStats(c1, c2);
                    }
                }
            }
        }
        System.out.println("Error counts: " + errs);
        System.out.println("Mean kmers error: " + ((double) kMersErrSum) / errs);
        System.out.println("Mean kmers error std dev: " + getKMersErrStdDev());
        System.out.println("Mean reads error: " + ((double) readsErrSum) / errs);
        System.out.println("Mean reads error std dev: " + getReadsErrStdDev());
    }

    protected double getKMersErrStdDev() {
        return Math.sqrt((kMersErrSquareSum  - ((double) kMersErrSum * kMersErrSum) / errs) / (errs - 1));
    }

    protected double getReadsErrStdDev() {
        return Math.sqrt((readsErrSquareSum  - ((double) readsErrSum * readsErrSum) / errs) / (errs - 1));
    }

    protected void sumErrorStats(CountsPerTaxid c1, CountsPerTaxid c2) {
            long k1 = c1 == null ? 0 : c1.getKMers();
            long k2 = c2 == null ? 0 : c2.getKMers();
            long err = Math.abs(k1 - k2);
            kMersErrSum += err;
            kMersErrSquareSum += err * err;

            long r1 = c1 == null ? 0 : c1.getReads();
            long r2 = c2 == null ? 0 : c2.getReads();
            err = Math.abs(r1 - r2);
            readsErrSum += err;
            readsErrSquareSum += err * err;
    }

    public Map<String, MatchingResult> match(String dbName, String csvFile, SmallTaxTree[] taxTreeRef) throws IOException {
        GSCommon config = new GSCommon(baseDir);
        GSProject project = new GSProject(config, dbName, null, null, csvFile, null, null, null,
                null, null, null, false);
        project.initConfigParam(GSConfigKey.KRAKEN_STYLE_MATCH, true);
        project.initConfigParam(GSConfigKey.THREADS, -1);
        GSMaker maker = new GSMaker(project);

        ObjectGoal<Database, GSProject> dbGoal = (ObjectGoal<Database, GSProject>) maker.getGoal(GSGoalKey.LOAD_DB);
        if (taxTreeRef != null) {
            taxTreeRef[0] = dbGoal.get().getTaxTree();
        }
        MatchResultGoal matchGoal = (MatchResultGoal) maker.getGoal(GSGoalKey.MATCHRES);
        Map<String, MatchingResult> matches = matchGoal.get();
        maker.dumpAll();
        return matches;
    }
}
