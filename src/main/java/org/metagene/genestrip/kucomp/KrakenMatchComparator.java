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
        SmallTaxTree tree1 = getDatabase(dbName1, false).getTaxTree();
        SmallTaxTree tree2 = getDatabase(dbName2, false).getTaxTree();

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
                for (SmallTaxTree.SmallTaxIdNode node1 : tree1) {
                    String taxId = node1.getTaxId();
                    SmallTaxTree.SmallTaxIdNode node2 = tree2.getNodeByTaxId(taxId);
                    // We only report on tax ids which are in both (Genestrip) dbs:
                    if (node2 != null) {
                        KrakenResCountGoal.KrakenResStats kustats1 = map1.get(taxId);
                        KrakenResCountGoal.KrakenResStats kustats2 = map2.get(taxId);
                        ps.print(taxId);
                        ps.print(';');
                        ps.print(getRankString(node1));
                        ps.print(';');
                        ps.print(correctDBValue(kustats1 == null ? 0 : kustats1.getKmers()));
                        ps.print(';');
                        ps.print(correctDBValue(kustats2 == null ? 0 : kustats2.getKmers()));
                        ps.print(';');
                        ps.print(correctDBValue(kustats1 == null ? 0 : kustats1.getReads()));
                        ps.print(';');
                        ps.print(correctDBValue(kustats2 == null ? 0 : kustats2.getReads()));
                        ps.println(';');
                    }
                }
            }
        }
    }

    public void compareWithKUResults(String dbName, String csvFile1, String csvFile2) throws IOException {
        GSCommon config = new GSCommon(baseDir);
        if (csvFile1 != null) {
            GSProject project = new GSProject(config, dbName, null, null, csvFile1, null, null, null,
                    null, null, null, false);
            GSMaker maker = new GSMaker(project);
            maker.getGoal(GSGoalKey.FASTA2FASTQ).make();
            maker.dumpAll();
        }

        GSProject project = new GSProject(config, dbName, null, null, csvFile2, null, null, null,
                null, null, null, false);
        project.initConfigParam(GSConfigKey.KRAKEN_STYLE_MATCH, true);

        GSMaker maker2 = new GSMaker(project);
        ObjectGoal<Database, GSProject> storeGoal = (ObjectGoal<Database, GSProject>) maker2.getGoal(GSGoalKey.LOAD_DB);
        SmallTaxTree tree = storeGoal.get().getTaxTree();
        ObjectGoal<Map<String, MatchingResult>, GSProject> matchResGoal = (ObjectGoal<Map<String, MatchingResult>, GSProject>) maker2.getGoal(GSGoalKey.MATCHRES);
        Map<String, MatchingResult> matchResult = matchResGoal.get();

        ObjectGoal<Map<String, List<KrakenResCountGoal.KrakenResStats>>, GSProject> countGoal =
                (ObjectGoal<Map<String, List<KrakenResCountGoal.KrakenResStats>>, GSProject>) maker2.getGoal(GSGoalKey.KRAKENCOUNT);
        Map<String, List<KrakenResCountGoal.KrakenResStats>> stats = countGoal.get();
        maker2.dumpAll();

        for (String key : matchResult.keySet()) {
            List<KrakenResCountGoal.KrakenResStats> list = stats.get(key);
            Map<String, KrakenResCountGoal.KrakenResStats> map = new HashMap<>();
            for (KrakenResCountGoal.KrakenResStats stat : list) {
                map.put(stat.getTaxid(), stat);
            }
            Map<String, CountsPerTaxid> gstats = new HashMap<>(matchResult.get(key).getTaxid2Stats());

            int differentKMerValues = 0;
            int differentReadValues = 0;

            File out = new File(baseDir, dbName + "_" + key + "_gs_ku_comp.csv");
            try (PrintStream ps = new PrintStream(new FileOutputStream(out))) {
                ps.println("taxid; rank; kmers 1; kmers 2; reads 1; reads 2");
                for (SmallTaxTree.SmallTaxIdNode node : tree) {
                    String taxId = node.getTaxId();
                    CountsPerTaxid gcounts = gstats.get(taxId);
                    KrakenResCountGoal.KrakenResStats kustats = map.get(taxId);
                    ps.print(taxId);
                    ps.print(';');
                    ps.print(getRankString(node));
                    ps.print(';');
                    long gkmers = gcounts == null ? 0 : gcounts.getKMers();
                    ps.print(correctDBValue(gkmers));
                    ps.print(';');
                    long kukmers = kustats == null ? 0 : kustats.getKmers();
                    ps.print(correctDBValue(kukmers));
                    ps.print(';');
                    long greads = gcounts == null ? 0 : gcounts.getReads();
                    ps.print(correctDBValue(greads));
                    ps.print(';');
                    long kureads = kustats == null ? 0 : kustats.getReads();
                    ps.print(correctDBValue(kureads));
                    ps.println(';');
                    if (kukmers != gkmers) {
                        differentKMerValues++;
                    }
                    if (kureads != greads) {
                        differentReadValues++;
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
