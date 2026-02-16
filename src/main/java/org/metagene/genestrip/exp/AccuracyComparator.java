package org.metagene.genestrip.exp;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.metagene.genestrip.*;
import org.metagene.genestrip.goals.MatchResultGoal;
import org.metagene.genestrip.goals.kraken.KrakenResCountGoal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.match.FastqKMerMatcher;
import org.metagene.genestrip.match.MatchingResult;
import org.metagene.genestrip.refseq.AccessionMap;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.util.ByteArrayUtil;

import java.io.*;
import java.util.*;

public class AccuracyComparator extends GenestripComparator {
    private static final CSVFormat GANON_ALL_FORMAT = CSVFormat.DEFAULT.builder().setQuote(null).setCommentMarker('#')
            .setDelimiter('\t').setRecordSeparator('\n').build();

    private final boolean format;
    private final TaxTree taxTree;
    private final AccessionMap accessionMap;

    public AccuracyComparator(File baseDir, boolean format) throws IOException {
        super(baseDir, null);
        this.format = format;

        GSCommon config = new GSCommon(baseDir);
        // Project name does not matter for as long as it exits.
        GSProject project = new GSProject(config, "viral", null, null, null, null, null, null,
                null, null, null, false);
        project.initConfigParam(GSConfigKey.THREADS, -1);

        GSMaker maker = new GSMaker(project);

        ObjectGoal<TaxTree, GSProject> taxTreeGoal = (ObjectGoal<TaxTree, GSProject>) maker.getGoal(GSGoalKey.TAXTREE);

        taxTree = taxTreeGoal.get();
        ObjectGoal<AccessionMap, GSProject> accessCollGoal = (ObjectGoal<AccessionMap, GSProject>) maker.getGoal(GSGoalKey.ACCMAP);
        accessionMap = accessCollGoal.get();
        accessCollGoal.cleanThis();

        maker.dumpAll();
    }

    public void writeReportFile(String db, String checkDB, String mapFile, boolean nanoSim) throws IOException {
        GSCommon config = new GSCommon(baseDir);
        GSProject project = new GSProject(config, db, null, null, null, null, null, null,
                null, null, null, false);
        SmallTaxTree checkTree = null;
        if (checkDB != null) {
            checkTree = getDatabase(checkDB, false).getTaxTree();
        }

        try (PrintStream ps = new PrintStream(new FileOutputStream(new File(project.getResultsDir(), db + (checkDB == null ? "" : "_" + checkDB) + "_accuracy.csv")))) {
            Map<String, int[]> resGenestrip = accuracyForSimulatedReadsGenestrip(db, mapFile, checkTree, nanoSim);
            Map<String, int[]> resKU = accuracyForSimulatedReadsKU(db, mapFile, checkTree, false, nanoSim);
            Map<String, int[]> resK2 = accuracyForSimulatedReadsKU(db, mapFile, checkTree, true, nanoSim);

            ps.println("fastq key; system; classified; correct genus; correct species; total; precision genus; recall genus; f1 genus; precision species; recall species; f1 species;");
            for (String fastqKey : resGenestrip.keySet()) {
                int[] counts = resGenestrip.get(fastqKey);
                int total = counts[5]; // No correct results without ground truth available.
                printCounts(ps, fastqKey, "genestrip", counts, total);

                counts = resKU.get(fastqKey);
                printCounts(ps, fastqKey, "krakenUniq", counts, total);

                counts = resK2.get(fastqKey);
                printCounts(ps, fastqKey, "k2", counts, total);

                counts = accuracyForSimulatedReadsGanon(db, "ganon/" + db + "_" + fastqKey + ".all", checkTree, nanoSim);
                printCounts(ps, fastqKey, "ganon", counts, total);

                counts = accuracyForSimulatedReadsGanon(db, "ganon/" + db + "_lowfp_" + fastqKey + ".all", checkTree, nanoSim);
                printCounts(ps, fastqKey, "ganon_lowfp", counts, total);
            }
        }
    }

    public void writeTickBorneSimReport() throws IOException {
        String db = "tick-borne";
        String checkDB = db;
        String mapFile = "ticks_sim.txt";
        GSCommon config = new GSCommon(baseDir);
        GSProject project = new GSProject(config, "tick-borne", null, null, null, null, null, null,
                null, null, null, false);
        SmallTaxTree checkTree = getDatabase(checkDB, false).getTaxTree();

        try (PrintStream ps = new PrintStream(new FileOutputStream(new File(project.getResultsDir(), db + (checkDB == null ? "" : "_" + checkDB) + "_accuracy.csv")))) {
            Map<String, int[]> resGenestrip = accuracyForSimulatedReadsGenestrip(db, mapFile, checkTree, true);
            Map<String, int[]> resKU = accuracyForSimulatedReadsKU("microbial", mapFile, checkTree, false, true);
            Map<String, int[]> resK2 = accuracyForSimulatedReadsKU("standard", mapFile, checkTree, true, true);

            ps.println("fastq key; system; classified; correct genus; correct species; total; precision genus; recall genus; f1 genus; precision species; recall species; f1 species;");
            for (String fastqKey : resGenestrip.keySet()) {
                int[] counts = resGenestrip.get(fastqKey);
                int total = counts[5]; // No correct results without ground truth available.
                printCounts(ps, fastqKey, "genestrip", counts, total);

                counts = resKU.get(fastqKey);
                printCounts(ps, fastqKey, "krakenUniq", counts, total);

                counts = resK2.get(fastqKey);
                printCounts(ps, fastqKey, "k2", counts, total);

                counts = accuracyForSimulatedReadsGanon("standard", "ganon/" + db + "_" + fastqKey + ".all", checkTree, true);
                printCounts(ps, fastqKey, "ganon", counts, total);
            }
        }
    }

    private void printCounts(PrintStream ps, String key, String system, int[] counts, int total) {
        ps.print(key);
        ps.print(';');
        ps.print(system);
        ps.print(';');
        ps.print(format(counts[0]));
        ps.print(';');
        ps.print(format(counts[1]));
        ps.print(';');
        ps.print(format(counts[2]));
        ps.print(';');
        ps.print(format(total));
        ps.print(';');
        for (int i = 1; i <= 2; i++) {
            double precision = ((double) counts[i]) / counts[0];
            double recall = ((double) counts[i]) / total;
            double f1 = 2 / ((1 / precision) + (1 / recall));
            ps.print(format(precision));
            ps.print(';');
            ps.print(format(recall));
            ps.print(';');
            ps.print(format(f1));
            ps.print(';');
        }
        ps.println();
    }

    private String format(int i) {
        return format ? LF.format(i) : String.valueOf(i);
    }

    private String format(double d) {
        return format ? DF.format(d) : String.valueOf(d);
    }

    public int[] accuracyForSimulatedReadsGanon(String db, String ganonReportFile, SmallTaxTree checkTree, boolean nanosim) throws IOException {
        Map<String, TaxTree.TaxIdNode> matchesMap = getUniqueClassMap(new File(ganonReportFile), false);
        int[] counters = new int[5];
        for (String descr : matchesMap.keySet()) {
            TaxTree.TaxIdNode node = nodeFromDesc(descr.getBytes(), nanosim);
            if (node != null) {
                TaxTree.TaxIdNode classNode = matchesMap.get(descr);
                if (isAsRequestedOrBelowInCheckTree(node, checkTree)) {
                    /*
                    if (node != classNode) {
                        System.out.println("stop");
                    }
                     */
                    updateMatchCounts(classNode, node, counters);
                }
                else if (classNode != null && isAsRequestedOrBelowInCheckTree(classNode, checkTree)) {
                    counters[0]++; // Count as classified.
                }
            } else {
                counters[4]++;
            }
        }
        System.out.println("+++ Ganon +++");
        System.out.println("Report file: " + ganonReportFile);
        printConsoleRes(counters);
        return counters;
    }

    protected TaxTree.TaxIdNode nodeFromDesc(byte[] desc, boolean nanosim) {
        if (!nanosim) {
            int startPos;
            if (desc[0] == '@') {
                startPos = desc[1] == '>' ? 2 : 1;
            }
            else {
                startPos = desc[0] == '>' ? 1 : 0;
            }
            int endPos = ByteArrayUtil.indexOf(desc, 5, desc.length, '_');
            // ByteArrayUtil.println(desc, System.out);
            return accessionMap.get(desc, startPos, endPos, false);
        }
        else {
            int startPos = ByteArrayUtil.indexOf(desc, 1, desc.length, '-');
            int endPos = ByteArrayUtil.indexOf(desc, 2, desc.length, '_');
            int nextDash = ByteArrayUtil.indexOf(desc, startPos + 1, desc.length, '-');
            desc[nextDash] = '_';
            return accessionMap.get(desc, startPos + 1, endPos, false, true);
        }
    }

    public Map<String, int[]> accuracyForSimulatedReadsGenestrip(String db, String csvFile2, SmallTaxTree checkTree, boolean nanosim) throws IOException {
        GSCommon config = new GSCommon(baseDir);
        GSProject project = new GSProject(config, db, null, null, csvFile2, null, null, null,
                null, null, null, false);
        project.initConfigParam(GSConfigKey.THREADS, -1);

        GSMaker maker = new GSMaker(project);
        MatchResultGoal matchResGoal = (MatchResultGoal) maker.getGoal(GSGoalKey.MATCHRES);
        Map<String, int[]> result = new HashMap<>();
        int[] counters = new int[6];
        matchResGoal.setAfterMatchCallback(new MatchResultGoal.AfterMatchCallback() {
            @Override
            public void afterMatch(FastqKMerMatcher.MatcherReadEntry myReadEntry, boolean b) {
                synchronized (counters) {
                    handleMatch(myReadEntry.classNode == null ? null : myReadEntry.classNode.getTaxId(),  myReadEntry.readDescriptor, counters, checkTree, nanosim);
                }
            }

            @Override
            public void afterKey(String key, MatchingResult res) {
                System.out.println("+++ Genestrip +++");
                System.out.println("Key: " + key);
                printConsoleRes(counters);
                System.out.println("Total count: " + counters[5]);
                result.put(key, Arrays.copyOf(counters, counters.length));
                Arrays.fill(counters, 0);
            }
        });
        matchResGoal.make();
        maker.dumpAll();
        return result;
    }

    protected void handleMatch(String classTaxId, byte[] desc, int[] counters, SmallTaxTree checkTree, boolean nanosim) {
        TaxTree.TaxIdNode node = nodeFromDesc(desc, nanosim);
        if (node != null) {
            // ByteArrayUtil.println(desc, System.out);
            TaxTree.TaxIdNode classNode = classTaxId == null ? null : taxTree.getNodeByTaxId(classTaxId);
            if (isAsRequestedOrBelowInCheckTree(node, checkTree)) {
                counters[5]++;
                updateMatchCounts(classNode, node, counters);
            }
            else if (classNode != null && isAsRequestedOrBelowInCheckTree(classNode, checkTree)) {
                counters[0]++; // Count as classified.
            }
        } else {
            counters[4]++;
        }
    }

    private boolean isAsRequestedOrBelowInCheckTree(TaxTree.TaxIdNode node, SmallTaxTree checkTree) {
        if (checkTree != null) {
            SmallTaxTree.SmallTaxIdNode snode = checkTree.getNodeByTaxId(node.getTaxId());
            while (snode != null) {
                if (snode.isRequested()) {
                    return true;
                }
                snode = snode.getParent();
            }
            return false;
        }
        else {
            return true;
        }
    }

    private void updateMatchCounts(TaxTree.TaxIdNode classNode, TaxTree.TaxIdNode node, int[] counters) {
        if (classNode != null) {
            counters[0]++;
            TaxTree.TaxIdNode lca = taxTree.getLowestCommonAncestor(node, classNode);
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

    public Map<String, int[]> accuracyForSimulatedReadsKU(String db, String csvFile2, SmallTaxTree checkTree, boolean k2, boolean nanosim) throws IOException {
        GSCommon config = new GSCommon(baseDir);
        GSProject project = new GSProject(config, db, null, null, csvFile2, null, null, null,
                null, null, null, false);
        project.initConfigParam(GSConfigKey.THREADS, -1);
        if (k2) {
            project.initConfigParam(GSConfigKey.KRAKEN_BIN, "./k2/kraken2/k2");
            project.initConfigParam(GSConfigKey.KRAKEN_DB, "./k2/" + db + "_db");
            project.initConfigParam(GSConfigKey.KRAKEN_EXEC_EXPR, "{0} classify --threads 10 --db {1} {2}");
        }

        GSMaker maker = new GSMaker(project);

        Map<String, int[]> result = new HashMap<>();
        int[] counters = new int[6];
        KrakenResCountGoal krakenResCountGoal = (KrakenResCountGoal) maker.getGoal(GSGoalKey.KRAKENCOUNT);
        krakenResCountGoal.setAfterMatchCallback(new KrakenResCountGoal.AfterMatchCallback() {
            @Override
            public void afterMatch(String krakenTaxid, byte[] readDescriptor) {
                synchronized (counters) {
                    handleMatch(krakenTaxid, readDescriptor, counters, checkTree, nanosim);
                }
            }

            @Override
            public void afterKey(String key, List<KrakenResCountGoal.KrakenResStats> res) {
                System.out.println("+++ Kraken 2 / KrakenUniq +++");
                System.out.println("Key: " + key);
                System.out.println("Total count: " + counters[5]);
                printConsoleRes(counters);
                result.put(key, Arrays.copyOf(counters, counters.length));
                Arrays.fill(counters, 0);
            }
        });
        krakenResCountGoal.make();
        maker.dumpAll();
        return result;
    }

    private void printConsoleRes(int[] counters) {
        System.out.println("Correct classifications GENUS: " + counters[1]);
        System.out.println("Correct classifications SPECIES: " + counters[2]);
        System.out.println("Correct classifications STRAIN: " + counters[3]);
        System.out.println("Bad Ground Truth Descriptor: " + counters[4]);
        System.out.println("Total classified: " + counters[0]);
    }

    public void writeReportFile2(String name, String db, String checkDB, String... fastqKeys) throws IOException {
        GSCommon config = new GSCommon(baseDir);
        GSProject project = new GSProject(config, db, null, null, null, null, null, null,
                null, null, null, false);
        SmallTaxTree checkTree = null;
        if (checkDB != null) {
            checkTree = getDatabase(checkDB, false).getTaxTree();
        }

        try (PrintStream ps = new PrintStream(new FileOutputStream(new File(project.getResultsDir(), db + "_" + name + "_rel_accuracy.csv")))) {
            ps.println("fastq key; system; classified; correct genus; correct species; total; precision genus; recall genus; f1 genus; precision species; recall species; f1 species;");
            for (String fastqKey : fastqKeys) {
                int[] counts = accuracyVia2ReportFiles(
                        "data/projects/" + db + "/krakenout/" + db + "_matchres_" + fastqKey + ".out",
                        "data/projects/" + checkDB + "/krakenout/" + checkDB + "_matchres_" + fastqKey + ".out", checkTree, true);
                printCounts(ps, fastqKey, "genestrip", counts, counts[5]);

                counts = accuracyVia2ReportFiles("ganon/" + db + "_" + fastqKey + ".all", "ganon/" + checkDB + "_" + fastqKey + ".all", checkTree, false);
                printCounts(ps, fastqKey, "ganon", counts, counts[5]);

                counts = accuracyVia2ReportFiles("ganon/" + db + "_lowfp_" + fastqKey + ".all", "ganon/" + checkDB + "_lowfp_" + fastqKey + ".all", checkTree, false);
                printCounts(ps, fastqKey, "ganon", counts, counts[5]);

                counts = accuracyVia2ReportFiles("ku/" + db + "_" + fastqKey + ".tsv", "ku/" + checkDB + "_" + fastqKey + ".tsv", checkTree, true);
                printCounts(ps, fastqKey, "krakenUniq", counts, counts[5]);

                counts = accuracyVia2ReportFiles("k2/" + db + "_" + fastqKey + ".tsv", "ku/" + checkDB + "_" + fastqKey + ".tsv", checkTree, true);
                printCounts(ps, fastqKey, "k2", counts, counts[5]);
            }
        }
    }

    protected int[] accuracyVia2ReportFiles(String groundTruthReportFile, String resFile, SmallTaxTree checkTree, boolean kraken) throws IOException {
        Map<String, TaxTree.TaxIdNode> gtMap = getUniqueClassMap(new File(groundTruthReportFile), kraken);
        int[] counters = new int[6];
        // Count the positives.
        for (String descr : gtMap.keySet()) {
            TaxTree.TaxIdNode match = gtMap.get(descr);
            if (isAsRequestedOrBelowInCheckTree(match, checkTree)) {
                counters[5]++;
            }
        }
        Map<String, TaxTree.TaxIdNode> resMap = getUniqueClassMap(new File(resFile), kraken);
        /*
        HashMap<TaxTree.TaxIdNode, Integer> truePositiveTaxons = new HashMap<>();
        HashMap<TaxTree.TaxIdNode, Integer> falsePositiveTaxons = new HashMap<>();
        HashMap<TaxTree.TaxIdNode, TaxTree.TaxIdNode> falsePositiveMismatchSamples = new HashMap<>();
        int total = 0;
         */
        for (String key : resMap.keySet()) {
            TaxTree.TaxIdNode classNode = resMap.get(key);
            TaxTree.TaxIdNode node = gtMap.get(key);
            if (node != null && isAsRequestedOrBelowInCheckTree(node, checkTree)) {
                /*
                Integer c = truePositiveTaxons.get(classNode);
                if (c == null) {
                    c = 0;
                }
                truePositiveTaxons.put(classNode, c + 1);
                 */
                updateMatchCounts(classNode, node, counters);
            }
            else if (isAsRequestedOrBelowInCheckTree(classNode, checkTree)) {
                /*
                Integer c = falsePositiveTaxons.get(classNode);
                if (c == null) {
                    c = 0;
                    falsePositiveMismatchSamples.put(classNode, node);
                }
                falsePositiveTaxons.put(classNode, c + 1);
                total++;
                 */
                counters[0]++; // Count as classified.
            }
        }
        /*
        System.out.println(falsePositiveMismatchSamples);
        System.out.println(falsePositiveTaxons);
        System.out.println(falsePositiveTaxons.size());
        System.out.println(total);
        System.out.println("True positive taxons: " + truePositiveTaxons.size());
        System.out.println(truePositiveTaxons);
         */
        return counters;
    }

    // Ganon does multiple matches.
    // We reduce it to a single match by computing the LCA in case of the multiple matches.
    // The method works also for Kraken-based output files.
    protected Map<String, TaxTree.TaxIdNode> getUniqueClassMap(File ganonFile, boolean kraken) throws IOException {
        Map<String, TaxTree.TaxIdNode> gtMap = new HashMap<>();
        try (CSVParser parser = GANON_ALL_FORMAT
                .parse(new InputStreamReader(new FileInputStream(ganonFile)))) {
            for (CSVRecord record : parser.getRecords()) {
                String descr = record.get(kraken ? 1 : 0);
                String taxid = record.get(kraken ? 2 : 1);
                TaxTree.TaxIdNode node = taxTree.getNodeByTaxId(taxid);
                if (node != null) {
                    TaxTree.TaxIdNode match = gtMap.get(descr);
                    if (match != null) {
                        match = taxTree.getLowestCommonAncestor(match, node);
                    } else {
                        match = node;
                    }
                    gtMap.put(descr, match);
                }
                else if (!"0".equals(taxid)) {
                    System.err.println("Invalid taxid in report file: " + taxid);
                }
            }
        }
        return gtMap;
    }

    public static void main(String[] args) throws IOException {
        AccuracyComparator comp = new AccuracyComparator(new File("./data"), false);
        /*
        comp.writeReportFile("viral", null, "viral_acc_comp.txt", false);
        comp.writeReportFile("human_virus", "human_virus", "viral_acc_comp.txt", false);

        comp.writeReportFile2("sim", "viral", "human_virus", "fastq1", "iss_hiseq", "iss_miseq");
         */

        comp.writeReportFile2("saliva", "viral", "human_virus", "ERR1395613");

        // Simulated tick files:
        //comp.writeTickBorneSimReport();
    }
}
