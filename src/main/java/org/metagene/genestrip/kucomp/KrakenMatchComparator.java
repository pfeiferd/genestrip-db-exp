package org.metagene.genestrip.kucomp;

import org.metagene.genestrip.*;
import org.metagene.genestrip.exp.GenestripComparator;
import org.metagene.genestrip.goals.MatchResultGoal;
import org.metagene.genestrip.goals.kraken.KrakenResCountGoal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.match.CountsPerTaxid;
import org.metagene.genestrip.match.FastqKMerMatcher;
import org.metagene.genestrip.match.MatchingResult;
import org.metagene.genestrip.refseq.AccessionMap;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.util.ByteArrayUtil;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class KrakenMatchComparator extends GenestripComparator {
    public KrakenMatchComparator(File baseDir) {
        super(baseDir);
    }

    public void compareKUWithKUResults(String dbName1, String dbName2, String csvFile) throws IOException {
        Database db1 = getDatabase(dbName1, false);
        SmallTaxTree tree1 = db1.getTaxTree();
        Database db2 = getDatabase(dbName2, false);
        SmallTaxTree tree2 = db2.getTaxTree();

        GSCommon config = new GSCommon(baseDir);

        GSProject project1 = new GSProject(config, dbName1, null, null, csvFile, null, null, null,
                null, null, null, false);
        GSMaker maker1 = new GSMaker(project1);
        ObjectGoal<Map<String, List<KrakenResCountGoal.KrakenResStats>>, GSProject> countGoal1 =
                (ObjectGoal<Map<String, List<KrakenResCountGoal.KrakenResStats>>, GSProject>) maker1.getGoal(GSGoalKey.KRAKENCOUNT);
        Map<String, List<KrakenResCountGoal.KrakenResStats>> stats1 = countGoal1.get();
        maker1.dumpAll();

        GSProject project2 = new GSProject(config, dbName2, null, null, csvFile, null, null, null,
                null, null, null, false);
        GSMaker maker2 = new GSMaker(project2);
        ObjectGoal<Map<String, List<KrakenResCountGoal.KrakenResStats>>, GSProject> countGoal2 =
                (ObjectGoal<Map<String, List<KrakenResCountGoal.KrakenResStats>>, GSProject>) maker2.getGoal(GSGoalKey.KRAKENCOUNT);
        Map<String, List<KrakenResCountGoal.KrakenResStats>> stats2 = countGoal2.get();
        maker2.dumpAll();

        for (String key : stats1.keySet()) {
            List<KrakenResCountGoal.KrakenResStats> list1 = stats1.get(key);
            Map<String, KrakenResCountGoal.KrakenResStats> map1 = new HashMap<>();
            for (KrakenResCountGoal.KrakenResStats stats : list1) {
                map1.put(stats.getTaxid(), stats);
            }
            List<KrakenResCountGoal.KrakenResStats> list2 = stats2.get(key);
            Map<String, KrakenResCountGoal.KrakenResStats> map2 = new HashMap<>();
            for (KrakenResCountGoal.KrakenResStats stats : list2) {
                map2.put(stats.getTaxid(), stats);
            }

            File out = new File(baseDir, dbName1 + "_" + dbName2 + "_" + key + "_ku_ku_comp.csv");
            try (PrintStream ps = new PrintStream(new FileOutputStream(out))) {
                ps.println("taxid; rank; kmers 1; kmers 2; reads 1; reads 2;");
                for (KrakenResCountGoal.KrakenResStats kustats1 : list1) {
                    SmallTaxTree.SmallTaxIdNode node1 = tree1.getNodeByTaxId(kustats1.getTaxid());
                    SmallTaxTree.SmallTaxIdNode node2 = tree2.getNodeByTaxId(kustats1.getTaxid());
                    // We only report on tax ids which are in both (Genestrip) dbs:
                    if (node1 != null && node2 != null) {
                        KrakenResCountGoal.KrakenResStats kustats2 = map2.get(kustats1.getTaxid());
                        ps.print(kustats1.getTaxid());
                        ps.print(';');
                        ps.print(getRankString(node1));
                        ps.print(';');
                        ps.print(correctDBValue(kustats1.getKmers()));
                        ps.print(';');
                        ps.print(correctDBValue(kustats2 == null ? 0 : kustats2.getKmers()));
                        ps.print(';');
                        ps.print(correctDBValue(kustats1.getReads()));
                        ps.print(';');
                        ps.print(correctDBValue(kustats2 == null ? 0 : kustats2.getReads()));
                        ps.println(';');
                    }
                }

                for (KrakenResCountGoal.KrakenResStats kustats2 : list2) {
                    if (!map1.containsKey(kustats2.getTaxid())) {
                        SmallTaxTree.SmallTaxIdNode node1 = tree1.getNodeByTaxId(kustats2.getTaxid());
                        SmallTaxTree.SmallTaxIdNode node2 = tree2.getNodeByTaxId(kustats2.getTaxid());
                        // We only report on tax ids which are in both (Genestrip) dbs:
                        if (node1 != null && node2 != null) {
                            ps.print(kustats2.getTaxid());
                            ps.print(';');
                            ps.print(getRankString(node1));
                            ps.print(';');
                            ps.print(correctDBValue(0));
                            ps.print(';');
                            ps.print(correctDBValue(kustats2.getKmers()));
                            ps.print(';');
                            ps.print(correctDBValue(0));
                            ps.print(';');
                            ps.print(correctDBValue(kustats2.getReads()));
                            ps.println(';');
                        }
                    }
                }
            }
        }
    }

    public void compareWithKUResults(String db, String csvFile1, String csvFile2) throws IOException {
        GSCommon config = new GSCommon(baseDir);
        if (csvFile1 != null) {
            GSProject project = new GSProject(config, db, null, null, csvFile1, null, null, null,
                    null, null, null, false);
            GSMaker maker = new GSMaker(project);
            maker.getGoal(GSGoalKey.FASTA2FASTQ).make();
            maker.dumpAll();
        }

        GSProject project = new GSProject(config, db, null, null, csvFile2, null, null, null,
                null, null, null, false);
        project.initConfigParam(GSConfigKey.KRAKEN_STYLE_MATCH, true);

        GSMaker maker2 = new GSMaker(project);
        ObjectGoal<Map<String, MatchingResult>, GSProject> matchResGoal = (ObjectGoal<Map<String, MatchingResult>, GSProject>) maker2.getGoal(GSGoalKey.MATCHRES);
        Map<String, MatchingResult> matchResult = matchResGoal.get();

        ObjectGoal<Map<String, List<KrakenResCountGoal.KrakenResStats>>, GSProject> countGoal =
                (ObjectGoal<Map<String, List<KrakenResCountGoal.KrakenResStats>>, GSProject>) maker2.getGoal(GSGoalKey.KRAKENCOUNT);
        Map<String, List<KrakenResCountGoal.KrakenResStats>> stats = countGoal.get();
        maker2.dumpAll();

        for (String key : matchResult.keySet()) {
            List<KrakenResCountGoal.KrakenResStats> list = stats.get(key);
            Map<String, CountsPerTaxid> gstats = new HashMap<>(matchResult.get(key).getTaxid2Stats());

            int differentKMerValues = 0;
            int differentReadValues = 0;

            File out = new File(baseDir, db + "_" + key + "_gs_ku_comp.csv");
            try (PrintStream ps = new PrintStream(new FileOutputStream(out))) {
                ps.println("taxid; rank; kmers 1; kmers 2; reads 1; reads 2");
                for (KrakenResCountGoal.KrakenResStats kustats : list) {
                    CountsPerTaxid gcounts = gstats.get(kustats.getTaxid());
                    ps.print(kustats.getTaxid());
                    ps.print(';');
                    ps.print(gcounts == null ? "" : gcounts.getRank());
                    ps.print(';');
                    long gkmers = gcounts == null ? 0 : gcounts.getKMers();
                    ps.print(correctDBValue(gkmers));
                    ps.print(';');
                    ps.print(correctDBValue(kustats.getKmers()));
                    ps.print(';');
                    long greads = gcounts == null ? 0 : gcounts.getReads();
                    ps.print(correctDBValue(greads));
                    ps.print(';');
                    ps.print(correctDBValue(kustats.getReads()));
                    ps.println(';');
                    gstats.remove(kustats.getTaxid());
                    if (kustats.getKmers() != gkmers) {
                        differentKMerValues++;
                    }
                    if (kustats.getReads() != greads) {
                        differentReadValues++;
                    }
                }

                for (CountsPerTaxid gcounts : gstats.values()) {
                    if (gcounts.getTaxid() != null) {
                        if (gcounts.getKMers() != 0 || gcounts.getReads() != 0) {
                            ps.print(gcounts.getTaxid());
                            ps.print(';');
                            ps.print(gcounts.getRank());
                            ps.print(';');
                            ps.print(correctDBValue(gcounts.getKMers()));
                            ps.print(';');
                            ps.print(correctDBValue(0));
                            ps.print(';');
                            ps.print(correctDBValue(gcounts.getReads()));
                            ps.print(';');
                            ps.print(correctDBValue(0));
                            ps.println(';');
                            if (gcounts.getKMers() != 0) {
                                differentKMerValues++;
                            }
                            if (gcounts.getReads() != 0) {
                                differentReadValues++;
                            }
                        }
                    }
                }
            }

            System.out.println("Different kmer values: " + differentKMerValues);
            System.out.println("Different read values: " + differentReadValues);
        }
    }

    public void accuracyCheckForSimulatedViralReads(String db, String csvFile2) throws IOException {
        GSCommon config = new GSCommon(baseDir);
        GSProject project2 = new GSProject(config, db, null, null, csvFile2, null, null, null,
                null, null, null, false);
        project2.initConfigParam(GSConfigKey.KRAKEN_STYLE_MATCH, true);
        project2.initConfigParam(GSConfigKey.THREADS, -1);

        GSMaker maker2 = new GSMaker(project2);

        ObjectGoal<AccessionMap, GSProject> accessCollGoal = (ObjectGoal<AccessionMap, GSProject>) maker2.getGoal(GSGoalKey.ACCMAP);
        AccessionMap map = accessCollGoal.get();
        ObjectGoal<Database, GSProject> dbGoal = (ObjectGoal<Database, GSProject>) maker2.getGoal(GSGoalKey.LOAD_DB);
        SmallTaxTree smallTaxTree = dbGoal.get().getTaxTree();
        MatchResultGoal matchResGoal = (MatchResultGoal) maker2.getGoal(GSGoalKey.MATCHRES);
        int[] counters = new int[4];
        matchResGoal.setAfterMatchCallback(new FastqKMerMatcher.AfterMatchCallback() {
            @Override
            public void afterMatch(FastqKMerMatcher.MyReadEntry myReadEntry, boolean b) {
                if (myReadEntry.classNode != null) {
                    byte[] desc = myReadEntry.readDescriptor;
                    int pos = ByteArrayUtil.indexOf(desc, 5, desc.length, '_');
                    TaxTree.TaxIdNode node = map.get(desc, 2, pos);
                    if (node != null) {
                        SmallTaxTree.SmallTaxIdNode snode = smallTaxTree.getNodeByTaxId(node.getTaxId());
                        SmallTaxTree.SmallTaxIdNode lca = smallTaxTree.getLeastCommonAncestor(snode, myReadEntry.classNode);
                        while (lca != null && Rank.NO_RANK.equals(lca.getRank())) {
                            lca = lca.getParent();
                        }
                        if (lca != null) {
                            Rank r = lca.getRank();
                            if (r != null) {
                                if (Rank.GENUS.equals(r) || r.isBelow(Rank.GENUS)) {
                                    counters[1]++;
                                }
                                if (Rank.SPECIES.equals(r) || r.isBelow(Rank.SPECIES)) {
                                    counters[2]++;
                                }
                                if (Rank.STRAIN.equals(r) || r.isBelow(Rank.STRAIN)) {
                                    counters[3]++;
                                }
                            }
                        }
                    }
                }
                counters[0]++;
            }
        });
        Map<String, MatchingResult> matchResult = matchResGoal.get();
        System.out.println("Correct classifications GENUS: " + counters[1]);
        System.out.println("Correct classifications SPECIES: " + counters[2]);
        System.out.println("Correct classifications STRAIN: " + counters[3]);
        System.out.println("Total count: " + counters[0]);
        maker2.dumpAll();
    }

    public static void main(String[] args) throws IOException {
        KrakenMatchComparator c = new KrakenMatchComparator(new File("./data"));
        c.compareWithKUResults("viral", "viral_ku_comp_fasta.txt", "viral_ku_comp.txt");
    }
}
