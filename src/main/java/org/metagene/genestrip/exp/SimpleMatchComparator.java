package org.metagene.genestrip.exp;

import it.unimi.dsi.fastutil.objects.Object2LongMap;
import it.unimi.dsi.fastutil.objects.Object2LongOpenHashMap;
import org.metagene.genestrip.*;
import org.metagene.genestrip.goals.MatchResultGoal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.match.CountsPerTaxid;
import org.metagene.genestrip.match.MatchingResult;
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
import java.util.HashMap;
import java.util.Locale;
import java.util.Map;

public class SimpleMatchComparator {
    protected final File baseDir;

    public SimpleMatchComparator(File baseDir) {
        this.baseDir = baseDir;
    }

    public void compareCommonDBEntries(String dbName1, String dbName2) throws IOException {
        Database db1 = getDatabase(dbName1, false);
        SmallTaxTree tree1 = db1.getTaxTree();
        Object2LongMap<String> stats1 = db1.getStats();
        db1 = null;
        Object2LongMap<String> stats2 = getDatabase(dbName2, false).getStats();
        File out = new File(baseDir, dbName1 + "_" + dbName2 + "_dbcomp.csv");
        try (PrintStream ps = new PrintStream(new FileOutputStream(out))) {
            ps.println("taxid; rank; " + dbName1 + " kmers; " + dbName2 + " kmers;");
            for (String key : stats1.keySet()) {
                if (key != null) {
                    long count1 = stats1.get(key);
                    SmallTaxTree.SmallTaxIdNode node = tree1.getNodeByTaxId(key);
                    if (node != null) {
                        long count2 = stats2.getOrDefault(key, -1);
                        if (count2 != -1) {
                            ps.print(key);
                            ps.print(";");
                            boolean norank = false;
                            while (node != null && Rank.NO_RANK.equals(node.getRank())) {
                                norank = true;
                                node = node.getParent();
                            }
                            Rank r = node == null ? null : node.getRank();
                            if (r != null && (r.isBelow(Rank.SPECIES) || (norank && r.equals(Rank.SPECIES)))) {
                                ps.print("below species");
                            }
                            else if (r != null && r.equals(Rank.SPECIES)) {
                                ps.print("species");
                            }
                            else if (r != null && r.equals(Rank.GENUS)) {
                                ps.print("genus");
                            }
                            else {
                                ps.print("null");
                            }
                            ps.print(";");
                            ps.print(correctDBValue(count1));
                            ps.print(";");
                            ps.print(correctDBValue(count2));
                            ps.println(";");
                        }
                    }
                }
            }
        }
    }

    protected long correctDBValue(long v) {
        return v + 1;
    }

    protected Database getDatabase(String dbName, boolean temp) throws IOException {
        GSCommon config = new GSCommon(baseDir);
        GSProject project = new GSProject(config, dbName, null, null, null, null, null, null,
                null, null, null, false);
        GSMaker maker = new GSMaker(project);
        ObjectGoal<Database, GSProject> storeGoal = (ObjectGoal<Database, GSProject>) maker.getGoal(temp ? GSGoalKey.LOAD_TEMPDB : GSGoalKey.LOAD_DB);
        Database db = storeGoal.get();
        maker.dumpAll();
        return db;
    }

    public void compareResults(String dbName1, String dbName2, String csvFile) throws IOException {
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

            File out = new File(baseDir, key + "_comp.csv");
            try (PrintStream ps = new PrintStream(new FileOutputStream(out))) {
                ps.println("taxid; rank; genestrip kmers 1; genestrip kmers 2; genestrip ukmers 1; genestrip ukmers 2; genestrip reads 1; genestrip reads 2");
                for (String taxid : stats1.keySet()) {
                    if (taxTreeRef2[0].getNodeByTaxId(taxid) == null) {
                        continue;
                    }
                    CountsPerTaxid c1 = stats1.get(taxid);
                    CountsPerTaxid c2 = stats2.get(taxid);
                    ps.print(c1.getTaxid());
                    ps.print(';');
                    ps.print(c1.getRank());
                    ps.print(';');
                    ps.print(c1.getKMers());
                    ps.print(';');
                    ps.print(c2 == null ? 0 : c2.getKMers());
                    ps.print(';');
                    ps.print(c1.getUniqueKMers());
                    ps.print(';');
                    ps.print(c2 == null ? 0 : c2.getUniqueKMers());
                    ps.print(';');
                    ps.print(c1.getReads());
                    ps.print(';');
                    ps.print(c2 == null ? 0 : c2.getReads());
                    ps.println(';');
                }
                for (String taxid : stats2.keySet()) {
                    if (taxTreeRef1[0].getNodeByTaxId(taxid) == null) {
                        continue;
                    }
                    CountsPerTaxid c1 = stats1.get(taxid);
                    if (c1 == null) {
                        CountsPerTaxid c2 = stats2.get(taxid);
                        ps.print(c2.getTaxid());
                        ps.print(';');
                        ps.print(c2.getRank());
                        ps.print(';');
                        ps.print(0);
                        ps.print(';');
                        ps.print(c2.getKMers());
                        ps.print(';');
                        ps.print(0);
                        ps.print(';');
                        ps.print(c2.getUniqueKMers());
                        ps.print(';');
                        ps.print(0);
                        ps.print(';');
                        ps.print(c2.getReads());
                        ps.println(';');
                    }
                }
            }
        }
    }

    public Map<String, MatchingResult> match(String dbName, String csvFile, SmallTaxTree[] taxTreeRef) throws IOException {
        GSCommon config = new GSCommon(baseDir);
        GSProject project = new GSProject(config, dbName, null, null, csvFile, null, null, null,
                null, null, null, false);
        project.initConfigParam(GSConfigKey.KRAKEN_STYLE_MATCH, true);
        project.initConfigParam(GSConfigKey.THREADS, -1);
        GSMaker maker = new GSMaker(project);

        ObjectGoal<Database, GSProject> dbGoal = (ObjectGoal<Database, GSProject>) maker.getGoal(GSGoalKey.LOAD_DB);
        taxTreeRef[0] = dbGoal.get().getTaxTree();
        MatchResultGoal matchGoal = (MatchResultGoal) maker.getGoal(GSGoalKey.MATCHRES);
        Map<String, MatchingResult> matches = matchGoal.get();
        maker.dumpAll();
        return matches;
    }
}
