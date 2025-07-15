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
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.*;

public class GenestripComparator {
    protected static final DecimalFormat DF = new DecimalFormat("#,###.00", new DecimalFormatSymbols(Locale.US));

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
                if (node2 != null) {
                    long count1 = stats1.getOrDefault(taxId, 0);
                    long count2 = stats2.getOrDefault(taxId, 0);
                    ps.print(taxId);
                    ps.print(";");
                    ps.print(getRankString(node1));
                    ps.print(";");
                    ps.print(correctDBValue(count1, false));
                    ps.print(";");
                    ps.print(correctDBValue(count2, false));
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
        } else if (r != null && (r.equals(Rank.GENUS) || r.isBelow(Rank.GENUS))) {
            return GENUS;
        } else {
            return ABOVE_GENUS;
        }
    }

    protected long correctDBValue(long v, boolean full) {
        return full ? v : v + 1;
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

    public Map<String, ErrCompInfo> compareResults(String dbName1, String dbName2, String csvFile) throws IOException {
        SmallTaxTree[] taxTreeRef1 = new SmallTaxTree[1];
        Map<String, MatchingResult> matches1 = match(dbName1, csvFile, taxTreeRef1);
        SmallTaxTree[] taxTreeRef2 = new SmallTaxTree[1];
        Map<String, MatchingResult> matches2 = match(dbName2, csvFile, taxTreeRef2);

        Map<String, ErrCompInfo> result = new LinkedHashMap<String, ErrCompInfo>(); // Maintains order of keys.
        for (String key : matches1.keySet()) {
            MatchingResult res1 = matches1.get(key);
            MatchingResult res2 = matches2.get(key);
            if (res1.getGlobalStats().getKMers() != res2.getGlobalStats().getKMers()) {
                throw new RuntimeException("Match results do not match");
            }

            ErrCompInfo errCompInfo = new ErrCompInfo(res1.getGlobalStats().getKMers(),
                    res1.getGlobalStats().getReads());
            result.put(key, errCompInfo);

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
                    long k1 = c1 == null ? 0 : c1.getKMers();
                    ps.print(correctDBValue(k1, false));
                    ps.print(';');
                    long k2 = c2 == null ? 0 : c2.getKMers();
                    ps.print(correctDBValue(k2, false));
                    ps.print(';');
                    ps.print(correctDBValue(c1 == null ? 0 : c1.getUniqueKMers(), false));
                    ps.print(';');
                    ps.print(correctDBValue(c2 == null ? 0 : c2.getUniqueKMers(), false));
                    ps.print(';');
                    long r1 = c1 == null ? 0 : c1.getReads();
                    ps.print(correctDBValue(r1, false));
                    ps.print(';');
                    long r2 = c2 == null ? 0 : c2.getReads();
                    ps.print(correctDBValue(r2, false));
                    ps.println(';');
                    if (SPECIES_OR_BELOW.equals(rs)) {
                        errCompInfo.sumErrorStats(k1, r1, k2, r2);
                    }
                }
            }
        }
        return result;
    }

    public void writeErrInfos(String dbName1, String dbName2, Map<String, ErrCompInfo> map1) throws IOException {
        File errOut = new File(baseDir, dbName1 + "_" + dbName2 + "_errors_gs_ku_comp.csv");
        try (PrintStream errPs = new PrintStream(new FileOutputStream(errOut))) {
            errPs.println("no; key; reads; gs read len; " +
                    "errs; kmer err; kmer err std dev; read err; read err std dev; kmer diffs; read diffs; kmer diffs %; read diffs %;");
            int counter = 0;
            for (String key : map1.keySet()) {
                ErrCompInfo errCompInfo1 = map1.get(key);
                errPs.print(counter);
                errPs.print(';');
                errPs.print(key);
                errPs.print(';');
                errPs.print(errCompInfo1.getReads());
                errPs.print(';');
                errPs.print(DF.format(((double) errCompInfo1.getKMers()) / errCompInfo1.getReads()));
                errPs.print(';');
                errPs.print(errCompInfo1.getErrs());
                errPs.print(';');
                errPs.print(DF.format(100 * errCompInfo1.getMeanKMersErr()));
                errPs.print(';');
                errPs.print(DF.format(100 * errCompInfo1.getKMersErrStdDev()));
                errPs.print(';');
                errPs.print(DF.format(100 * errCompInfo1.getMeanReadsErr()));
                errPs.print(';');
                errPs.print(DF.format(100 * errCompInfo1.getReadsErrStdDev()));
                errPs.print(';');
                errPs.print(errCompInfo1.getKmerDiffs());
                errPs.print(';');
                errPs.print(errCompInfo1.getReadDiffs());
                errPs.println(';');
                errPs.print(';');
                errPs.print(DF.format(100 * errCompInfo1.getKmerDiffsRatio()));
                errPs.print(';');
                errPs.print(DF.format(100 * errCompInfo1.getReadDiffsRatio()));
                counter++;
            }
        }
    }

    public void combineErrInfos(String dbName1, String dbName2, Map<String, ErrCompInfo> map1, Map<String, ErrCompInfo> map2) throws IOException {
        File errOut = new File(baseDir, dbName1 + "_" + dbName2 + "_errors_gs_ku_comp.csv");
        try (PrintStream errPs = new PrintStream(new FileOutputStream(errOut))) {
            errPs.println("no; key; reads; gs read len; " +
                    "gs errs; gs kmer err; gs kmer err std dev; gs read err; gs read err std dev; " +
                    "ku errs; ku kmer err; ku kmer err std dev; ku read err; ku read err std dev;");
            int counter = 0;
            for (String key : map1.keySet()) {
                ErrCompInfo errCompInfo1 = map1.get(key);
                ErrCompInfo errCompInfo2 = map2 != null ? map2.get(key) : null;
                if (errCompInfo2 == null) {
                    errCompInfo2 = new ErrCompInfo(0,0);
                }
                errPs.print(counter);
                errPs.print(';');
                errPs.print(key);
                errPs.print(';');
                errPs.print(errCompInfo1.getReads());
                errPs.print(';');
                errPs.print(DF.format(((double) errCompInfo1.getKMers()) / errCompInfo1.getReads()));
                errPs.print(';');
                errPs.print(errCompInfo1.getErrs());
                errPs.print(';');
                errPs.print(DF.format(100 * errCompInfo1.getMeanKMersErr()));
                errPs.print(';');
                errPs.print(DF.format(100 * errCompInfo1.getKMersErrStdDev()));
                errPs.print(';');
                errPs.print(DF.format(100 * errCompInfo1.getMeanReadsErr()));
                errPs.print(';');
                errPs.print(DF.format(100 * errCompInfo1.getReadsErrStdDev()));
                errPs.print(';');
                errPs.print(errCompInfo2.getErrs());
                errPs.print(';');
                errPs.print(DF.format(100 * errCompInfo2.getMeanKMersErr()));
                errPs.print(';');
                errPs.print(DF.format(100 * errCompInfo2.getKMersErrStdDev()));
                errPs.print(';');
                errPs.print(DF.format(100 * errCompInfo2.getMeanReadsErr()));
                errPs.print(';');
                errPs.print(DF.format(100 * errCompInfo2.getReadsErrStdDev()));
                errPs.println(';');
                counter++;
            }
        }
    }

    public Map<String, MatchingResult> match(String dbName, String csvFile, SmallTaxTree[] taxTreeRef) throws IOException {
        GSCommon config = new GSCommon(baseDir);
        GSProject project = new GSProject(config, dbName, null, null, csvFile, null, null, null,
                null, null, null, false);
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

    public static class ErrCompInfo {
        private long kmers;
        private long reads;

        private int errs;
        private double kMersErrSum;
        private double kMersErrSquareSum;
        private double readsErrSum;
        private double readsErrSquareSum;

        private int kmerDiffs;
        private int readDiffs;

        public ErrCompInfo(long kmers, long reads) {
            this.kmers = kmers;
            this.reads = reads;
        }

        public int getKmerDiffs() {
            return kmerDiffs;
        }

        public int getReadDiffs() {
            return readDiffs;
        }

        public double getKmerDiffsRatio() {
            return ((double) kmerDiffs) / errs;
        }

        public double getReadDiffsRatio() {
            return ((double) readDiffs) / errs;
        }

        public long getKMers() {
            return kmers;
        }

        public long getReads() {
            return reads;
        }

        public int getErrs() {
            return errs;
        }

        public double getKMersErrStdDev() {
            return Math.sqrt((kMersErrSquareSum - ((double) kMersErrSum * kMersErrSum) / errs) / (errs - 1));
        }

        public double getReadsErrStdDev() {
            return Math.sqrt((readsErrSquareSum - ((double) readsErrSum * readsErrSum) / errs) / (errs - 1));
        }

        public double getMeanKMersErr() {
            return kMersErrSum / errs;
        }

        public double getMeanReadsErr() {
            return readsErrSum / errs;
        }

        public void sumErrorStats(long k1, long r1, long k2, long r2) {
            errs++;
            if (k1 != k2) {
                kmerDiffs++;
            }
            if (r1 != r2) {
                readDiffs++;
            }
            double err = ((double) Math.abs(k1 - k2)) / (k1 + 1);
            kMersErrSum += err;
            kMersErrSquareSum += err * err;

            err = ((double) Math.abs(r1 - r2)) / (r1 + 1);
            readsErrSum += err;
            readsErrSquareSum += err * err;
        }
    }
}
