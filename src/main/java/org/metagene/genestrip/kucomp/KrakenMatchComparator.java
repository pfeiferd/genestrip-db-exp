package org.metagene.genestrip.kucomp;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.metagene.genestrip.*;
import org.metagene.genestrip.exp.GenestripComparator;
import org.metagene.genestrip.goals.MatchResultGoal;
import org.metagene.genestrip.goals.kraken.KrakenResCountGoal;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.io.StreamingResourceStream;
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

import java.io.*;
import java.util.*;

import static org.metagene.genestrip.GSGoalKey.FASTQ_MAP;

public class KrakenMatchComparator extends GenestripComparator {
    public KrakenMatchComparator(File baseDir) {
        super(baseDir);
    }

    public Map<String, ErrCompInfo> compareKUWithKUResults(String dbName1, String dbName2, String csvFile) throws IOException {
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

        Map<String, ErrCompInfo> result = new LinkedHashMap<>(); // Maintains order of keys
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
            ErrCompInfo errCompInfo = new ErrCompInfo(0, 0);
            result.put(key, errCompInfo);

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
                        String rs = getRankString(node1);
                        ps.print(rs);
                        ps.print(';');
                        long k1 = kustats1 == null ? 0 : kustats1.getKmers();
                        ps.print(LF.format(correctDBValue(k1, false)));
                        ps.print(';');
                        long k2 = kustats2 == null ? 0 : kustats2.getKmers();
                        ps.print(LF.format(correctDBValue(k2, false)));
                        ps.print(';');
                        long r1 = kustats1 == null ? 0 : kustats1.getReads();
                        ps.print(LF.format(correctDBValue(r1, false)));
                        ps.print(';');
                        long r2 = kustats2 == null ? 0 : kustats2.getReads();
                        ps.print(LF.format(correctDBValue(r2, false)));
                        ps.println(';');
                        if (SPECIES_OR_BELOW.equals(rs)) {
                            errCompInfo.sumErrorStats(k1, r1, k2, r2);
                        }
                    }
                }
            }
        }
        return result;
    }

    private static final CSVFormat FORMAT = CSVFormat.DEFAULT.builder().setQuote(null).setCommentMarker('#')
            .setDelimiter(';').setRecordSeparator('\n').setHeader().build();

    public void aggregateCompareWithKUResults(String dbName, String csvFile1) throws IOException {
        GSCommon config = new GSCommon(baseDir);
        GSProject project = new GSProject(config, dbName, null, null, csvFile1, null, null, null,
                null, null, null, false);
        GSMaker maker2 = new GSMaker(project);

        ObjectGoal<Map<String, StreamingResourceStream>, GSProject> mapGoal = (ObjectGoal<Map<String, StreamingResourceStream>, GSProject>) maker2.getGoal(FASTQ_MAP);
        ObjectGoal<Database, GSProject> storeGoal = (ObjectGoal<Database, GSProject>) maker2.getGoal(GSGoalKey.LOAD_DB);
        SmallTaxTree tree = storeGoal.get().getTaxTree();

        for (String key : mapGoal.get().keySet()) {
            File in = new File(baseDir, dbName + "_" + key + "_gs_ku_comp.csv");
            if (in.exists()) {
                Map<String, long[]> sumsMap = new HashMap<>();
                try (CSVParser parser = FORMAT
                        .parse(new InputStreamReader(StreamProvider.getInputStreamForFile(in)))) {

                    for (CSVRecord record : parser) {
                        String taxid = record.get("taxid");
                        SmallTaxTree.SmallTaxIdNode node = tree.getNodeByTaxId(taxid);
                        while (node != null) {
                            if (Rank.GENUS.equals(node.getRank())) {
                                break;
                            }
                            node = node.getParent();
                        }
                        if (node != null) {
                            long[] sums = sumsMap.get(node.getTaxId());
                            if (sums == null) {
                                sums = new long[4];
                                sumsMap.put(node.getTaxId(), sums);
                            }
                            sums[0] += Long.valueOf(record.get("kmers 1"));
                            sums[1] += Long.valueOf(record.get("kmers 2"));
                            sums[2] += Long.valueOf(record.get("reads 1"));
                            sums[3] += Long.valueOf(record.get("reads 2"));
                        } else {
                            System.err.println("Warning missing taxid node for: " + taxid);
                        }
                    }
                }
                File out = new File(baseDir, dbName + "_" + key + "_genus_agg_gs_ku_comp.csv");
                try (PrintStream ps = new PrintStream(new FileOutputStream(out))) {
                    ps.println("taxid;kmers 1;kmers 2;reads 1;reads 2");
                    for (String taxid : sumsMap.keySet()) {
                        long[] sums = sumsMap.get(taxid);
                        ps.print(taxid);
                        ps.print(';');
                        for (int i = 0; i < sums.length; i++) {
                            ps.print(LF.format(sums[i]));
                            ps.print(';');
                        }
                        ps.println();
                    }
                }
            }
        }
    }

    public Map<String, ErrCompInfo> compareWithKUResults(String dbName, String kuDBName, String csvFile1, String csvFile2, boolean full) throws IOException {
        GSCommon config = new GSCommon(baseDir);
        if (csvFile1 != null) {
            GSProject project = new GSProject(config, dbName, null, null, csvFile1, null, null, null,
                    null, null, null, false);
            project.initConfigParam(GSConfigKey.ALWAYS_ASSUME_GZIP, false);
            GSMaker maker = new GSMaker(project);
            maker.getGoal(GSGoalKey.FASTA2FASTQ).make();
            maker.dumpAll();
        }

        GSProject project = new GSProject(config, dbName, null, null, csvFile2, null, null, null,
                null, null, null, false);

        GSMaker maker2 = new GSMaker(project);
        ObjectGoal<Map<String, List<KrakenResCountGoal.KrakenResStats>>, GSProject> countGoal =
                (ObjectGoal<Map<String, List<KrakenResCountGoal.KrakenResStats>>, GSProject>) maker2.getGoal(GSGoalKey.KRAKENCOUNT);

        Map<String, List<KrakenResCountGoal.KrakenResStats>> stats = countGoal.get();
        ObjectGoal<Database, GSProject> storeGoal = (ObjectGoal<Database, GSProject>) maker2.getGoal(GSGoalKey.LOAD_DB);
        SmallTaxTree tree = storeGoal.get().getTaxTree();

        ObjectGoal<Map<String, MatchingResult>, GSProject> matchResGoal = (ObjectGoal<Map<String, MatchingResult>, GSProject>) maker2.getGoal(GSGoalKey.MATCHRES);
        Map<String, MatchingResult> matchResult = matchResGoal.get();
        maker2.dumpAll();

        KrakenDBComparator krakenDBComparator = new KrakenDBComparator(baseDir);
        Map<String, Long> kuTaxid2KMer = krakenDBComparator.getKrakenDBCounts(krakenDBComparator.getKrakenCountsFile(kuDBName));

        Map<String, ErrCompInfo> result = new LinkedHashMap<>(); // Maintains order of keys...
        for (String key : matchResult.keySet()) {
            MatchingResult res1 = matchResult.get(key);

            ErrCompInfo errCompInfo = new ErrCompInfo(res1.getGlobalStats().getKMers(),
                    res1.getGlobalStats().getReads());
            result.put(key, errCompInfo);

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
                if (full) {
                    ps.print("name;");
                }
                ps.println("taxid;rank;kmers 1;kmers 2;reads 1;reads 2");
                for (SmallTaxTree.SmallTaxIdNode node : tree) {
                    String taxId = node.getTaxId();
                    if (full || kuTaxid2KMer.get(taxId) != null) {
                        CountsPerTaxid gcounts = gstats.get(taxId);
                        KrakenResCountGoal.KrakenResStats kustats = map.get(taxId);
                        map.remove(taxId);
                        long gkmers = gcounts == null ? 0 : gcounts.getKMers();
                        long kukmers = kustats == null ? 0 : kustats.getKmers();
                        long greads = gcounts == null ? 0 : gcounts.getReads();
                        long kureads = kustats == null ? 0 : kustats.getReads();
                        // No need to report if all is zero.
                        if (gkmers != 0 || kukmers != 0 || greads != 0 || kureads != 0) {
                            if (full) {
                                ps.print(node.getName());
                                ps.print(';');
                            }
                            ps.print(taxId);
                            ps.print(';');
                            String rs = full ? (node.getRank() == null ? null : node.getRank().getName()) : getRankString(node);
                            ps.print(rs);
                            ps.print(';');
                            ps.print(LF.format(correctDBValue(gkmers, full)));
                            ps.print(';');
                            ps.print(LF.format(correctDBValue(kukmers, full)));
                            ps.print(';');
                            ps.print(LF.format(correctDBValue(greads, full)));
                            ps.print(';');
                            ps.print(LF.format(correctDBValue(kureads, full)));
                            ps.println(';');
                            if (kukmers != gkmers) {
                                differentKMerValues++;
                            }
                            if (kureads != greads) {
                                differentReadValues++;
                            }
                            errCompInfo.sumErrorStats(gkmers, greads, kukmers, kureads);
                        }
                    }
                }
                if (full) {
                    for (String taxId : map.keySet()) {
                        KrakenResCountGoal.KrakenResStats kustats = map.get(taxId);
                        long kukmers = kustats == null ? 0 : kustats.getKmers();
                        long kureads = kustats == null ? 0 : kustats.getReads();
                        ps.print("??");
                        ps.print(';');
                        ps.print(taxId);
                        ps.print(';');
                        String rs = null;
                        ps.print(rs);
                        ps.print(';');
                        ps.print(LF.format(correctDBValue(0, full)));
                        ps.print(';');
                        ps.print(LF.format(correctDBValue(kukmers, full)));
                        ps.print(';');
                        ps.print(LF.format(correctDBValue(0, full)));
                        ps.print(';');
                        ps.print(LF.format(correctDBValue(kureads, full)));
                        ps.println(';');
                        if (kukmers != 0) {
                            differentKMerValues++;
                        }
                        if (kureads != 0) {
                            differentReadValues++;
                        }
                        errCompInfo.sumErrorStats(0, 0, kukmers, kureads);
                    }
                }
            }

            System.out.println("Different kmer values: " + differentKMerValues);
            System.out.println("Different read values: " + differentReadValues);
        }
        return result;
    }

    public void accuracyCheckForSimulatedViralReads(String db, String csvFile2) throws IOException {
        GSCommon config = new GSCommon(baseDir);
        GSProject project2 = new GSProject(config, db, null, null, csvFile2, null, null, null,
                null, null, null, false);
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
                synchronized (counters) {
                    if (myReadEntry.classNode != null) {
                        byte[] desc = myReadEntry.readDescriptor;
                        int pos = ByteArrayUtil.indexOf(desc, 5, desc.length, '_');
                        TaxTree.TaxIdNode node = map.get(desc, 2, pos, false);
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
            }
        });
        matchResGoal.make();
        System.out.println("Correct classifications GENUS: " + counters[1]);
        System.out.println("Correct classifications SPECIES: " + counters[2]);
        System.out.println("Correct classifications STRAIN: " + counters[3]);
        System.out.println("Total count: " + counters[0]);
        maker2.dumpAll();
    }
}
