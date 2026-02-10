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
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.util.ByteArrayUtil;

import java.io.*;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
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

    public void writeReportFile(String db, String checkDB, String... fastqKeys) throws IOException {
        GSCommon config = new GSCommon(baseDir);
        GSProject project = new GSProject(config, db, null, null, null, null, null, null,
                null, null, null, false);
        SmallTaxTree checkTree = null;
        if (checkDB != null) {
            checkTree = getDatabase(checkDB, false).getTaxTree();
        }

        try (PrintStream ps = new PrintStream(new FileOutputStream(new File(project.getResultsDir(), db + "_accuracyReport.csv")))) {
            Map<String, int[]> resGenestrip = accuracyForSimulatedReadsGenestrip(db, "viral_acc_comp.txt", checkTree);
            Map<String, int[]> resKU = accuracyForSimulatedReadsKU(db, "viral_acc_comp.txt", checkTree);

            ps.println("fastq key; system; correct genus; correct species; total; precision genus; recall genus; f1 genus; precision species; recall species; f1 species;");
            for (String fastqKey : fastqKeys) {
                int[] counts = resGenestrip.get(fastqKey);
                int total = counts[5]; // No correct results without ground truth available.
                printCounts(ps, fastqKey, "genestrip", counts, total);

                counts = resKU.get(fastqKey);
                printCounts(ps, fastqKey, "krakenUniq", counts, total);

                counts = accuracyForSimulatedReadsGanon(db, "ganon/" + db + "_" + fastqKey + ".all", checkTree);
                printCounts(ps, fastqKey, "ganon", counts, total);

                counts = accuracyForSimulatedReadsGanon(db, "ganon/" + db + "_lowfp_" + fastqKey + ".all", checkTree);
                printCounts(ps, fastqKey, "ganon_lowfp", counts, total);
            }
        }
    }

    private void printCounts(PrintStream ps, String key, String system, int[] counts, int total) {
        ps.print(key);
        ps.print(';');
        ps.print(system);
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

    public int[] accuracyForSimulatedReadsGanon(String db, String ganonReportFile, SmallTaxTree checkTree) throws IOException {
        Map<String, TaxTree.TaxIdNode> matchesMap = new HashMap<String, TaxTree.TaxIdNode>();
        try (CSVParser parser = GANON_ALL_FORMAT
                .parse(new InputStreamReader(new FileInputStream(new File(ganonReportFile))))) {
            // Ganon does multiple matches.
            // We reduce it to a single match by computing the LCA in case of the multiple matches.
            for (CSVRecord record : parser.getRecords()) {
                String descr = record.get(0);
                String ganonTaxid = record.get(1);
                if (ganonTaxid != null) {
                    TaxTree.TaxIdNode classNode = taxTree.getNodeByTaxId(ganonTaxid);
                    if (classNode != null) {
                        TaxTree.TaxIdNode match = matchesMap.get(descr);
                        if (match != null) {
                            match = taxTree.getLeastCommonAncestor(match, classNode);
                        } else {
                            match = classNode;
                        }
                        matchesMap.put(descr, match);
                    } else {
                        System.out.println("Invalid taxid from Ganon: " + ganonTaxid);
                    }
                }
            }
        }

        int[] counters = new int[5];
        for (String descr : matchesMap.keySet()) {
            byte[] desc = descr.getBytes();
            int pos = ByteArrayUtil.indexOf(desc, 5, desc.length, '_');
            TaxTree.TaxIdNode node = accessionMap.get(desc, desc[0] == '>' ? 1 : 0, pos, false);
            if (node != null) {
                if (isAsSpeciesOrBelowInCheckTree(node, checkTree)) {
                    TaxTree.TaxIdNode classNode = matchesMap.get(descr);
                    updateMatchCounts(classNode, node, counters);
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

    public Map<String, int[]> accuracyForSimulatedReadsGenestrip(String db, String csvFile2, SmallTaxTree checkTree) throws IOException {
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
                    handleMatch(myReadEntry.classNode == null ? null : myReadEntry.classNode.getTaxId(),  myReadEntry.readDescriptor, counters, checkTree);
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

    protected void handleMatch(String classTaxId, byte[] desc, int[] counters, SmallTaxTree checkTree) {
        int pos = ByteArrayUtil.indexOf(desc, 5, desc.length, '_');
        TaxTree.TaxIdNode node = accessionMap.get(desc, desc[1] == '>' ? 2 : 1, pos, false);
        if (node != null) {
            if (isAsSpeciesOrBelowInCheckTree(node, checkTree)) {
                TaxTree.TaxIdNode classNode = classTaxId == null ? null : taxTree.getNodeByTaxId(classTaxId);
                updateMatchCounts(classNode, node, counters);
                counters[5]++;
            }
        } else {
            counters[4]++;
        }
    }

    private boolean isAsSpeciesOrBelowInCheckTree(TaxTree.TaxIdNode node, SmallTaxTree checkTree) {
        if (checkTree != null) {
            SmallTaxTree.SmallTaxIdNode snode = checkTree.getNodeByTaxId(node.getTaxId());
            while (snode != null && !snode.getRank().equals(Rank.SPECIES)) {
                snode = snode.getParent();
            }
            return snode != null;
        }
        else {
            return true;
        }
    }

    private void updateMatchCounts(TaxTree.TaxIdNode classNode, TaxTree.TaxIdNode node, int[] counters) {
        if (classNode != null) {
            counters[0]++;
            TaxTree.TaxIdNode lca = taxTree.getLeastCommonAncestor(node, classNode);
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

    public Map<String, int[]> accuracyForSimulatedReadsKU(String db, String csvFile2, SmallTaxTree checkTree) throws IOException {
        GSCommon config = new GSCommon(baseDir);
        GSProject project = new GSProject(config, db, null, null, csvFile2, null, null, null,
                null, null, null, false);
        project.initConfigParam(GSConfigKey.THREADS, -1);

        GSMaker maker = new GSMaker(project);

        Map<String, int[]> result = new HashMap<>();
        int[] counters = new int[6];
        KrakenResCountGoal krakenResCountGoal = (KrakenResCountGoal) maker.getGoal(GSGoalKey.KRAKENCOUNT);
        krakenResCountGoal.setAfterMatchCallback(new KrakenResCountGoal.AfterMatchCallback() {
            @Override
            public void afterMatch(String krakenTaxid, byte[] readDescriptor) {
                synchronized (counters) {
                    handleMatch(krakenTaxid, readDescriptor, counters, checkTree);
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

    public static void main(String[] args) throws IOException {
        AccuracyComparator comp = new AccuracyComparator(new File("./data"), false);
        comp.writeReportFile("viral", null, "fastq1", "iss_hiseq", "iss_miseq");
        comp.writeReportFile("human_virus", "human_virus", "fastq1", "iss_hiseq", "iss_miseq");
    }
}
