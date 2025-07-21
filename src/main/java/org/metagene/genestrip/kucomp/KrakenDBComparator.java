package org.metagene.genestrip.kucomp;

import it.unimi.dsi.fastutil.objects.Object2LongMap;
import org.metagene.genestrip.GSCommon;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSMaker;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.exp.GenestripComparator;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.TaxTree;

import java.io.*;
import java.util.*;

public class KrakenDBComparator extends GenestripComparator {
    public KrakenDBComparator(File baseDir) {
        super(baseDir);
    }

    public void reportKrakenDBComparison(String genestripDB, String krakenDB1, String krakenDB2) throws IOException {
        SmallTaxTree tree = getDatabase(genestripDB, false).getTaxTree();
        Map<String, Long> kuTaxid2KMer1 = getKrakenDBCounts(getKrakenCountsFile(krakenDB1));
        Map<String, Long> kuTaxid2KMer2 = getKrakenDBCounts(getKrakenCountsFile(krakenDB2));

        int diff = 0;
        File out = new File(baseDir, krakenDB1 + "_" + krakenDB2 + "_ku_ku_dbcomp.csv");
        try (PrintStream ps = new PrintStream(new FileOutputStream(out))) {
            ps.println("taxid; rank; kmers 1; kmers 2;");
            for (SmallTaxTree.SmallTaxIdNode node : tree) {
                // We only report on tax ids which are in respective Genestrip db:
                if (node != null) {
                    String key = node.getTaxId();
                    Long count2 = kuTaxid2KMer2.get(key);
                    if (count2 == null) {
                        count2 = 0L;
                    }
                    Long count1 = kuTaxid2KMer1.get(key);
                    if (count1 == null) {
                        count1 = 0L;
                    }
                    ps.print(key);
                    ps.print(";");
                    String rs = getRankString(node);
                    ps.print(rs);
                    if (SPECIES_OR_BELOW.equals(rs)) {
                        if (!count1.equals(count2)) {
                            diff++;
                        }
                    }
                    ps.print(";");
                    ps.print(LF.format(correctDBValue(count1, false)));
                    ps.print(";");
                    ps.print(LF.format(correctDBValue(count2, false)));
                    ps.println(";");
                }
            }
        }
        System.out.println(diff);
    }

    public void reportKMerComparisons(String genestripDB, String krakenDB, String filterTaxid, boolean full) throws IOException {
        Map<String, Long> kuTaxid2KMer = getKrakenDBCounts(getKrakenCountsFile(krakenDB));

        File countsFile = new File(baseDir, genestripDB + (filterTaxid == null ? "" : "_" + filterTaxid) + "_gs_ku_dbcomp.csv");
        try (PrintStream out = new PrintStream(new FileOutputStream(countsFile))) {
            printJointStoreInfo(getDatabase(genestripDB, false), out, kuTaxid2KMer, filterTaxid, full);
        }
    }

    protected Map<String, Long> getKrakenDBCounts(File file) throws IOException {
        Map<String, Long> kuTaxid2KMer = new HashMap<>();
        try (BufferedReader br = new BufferedReader(new FileReader(file))) {
            String line = br.readLine();
            while (line != null) {
                int tab = line.indexOf('\t');
                if (tab != -1) {
                    String taxid = line.substring(0, tab);
                    Long kmers = Long.valueOf(line.substring(tab + 1));
                    kuTaxid2KMer.put(taxid, kmers);
                }
                line = br.readLine();
            }
        }
        return kuTaxid2KMer;
    }

    protected File getKrakenCountsFile(String krakenDB) {
        return new File(baseDir, "../ku/" + krakenDB + "/database.kdb.counts");
    }

    public void printJointStoreInfo(Database database, PrintStream out, Map<String, Long> kuTaxid2KMer, String filterTaxid, boolean full) {
        out.println("taxid; rank; kmers 1; kmers 2;");

        long inDataErr = 0;
        long inDataEntries = 0;
        long err = 0;
        long entries = 0;
        List<SmallTaxTree.SmallTaxIdNode> missingNodesInKu = new ArrayList<>();
        Map<SmallTaxTree.SmallTaxIdNode, Long> differences = new HashMap<>();

        Map<String, Long> kmerSumsPerGenusG = new HashMap<>();
        Map<String, Long> kmerSumsPerGenusKU = new HashMap<>();
        Object2LongMap<String> stats = database.getStats();
        for (SmallTaxTree.SmallTaxIdNode taxNode : database.getTaxTree()) {
            if (taxNode == null || taxNode.getRank() == null) {
                continue;
            }
            if (filterTaxid != null) {
                boolean found = false;
                for (SmallTaxTree.SmallTaxIdNode n = taxNode; n != null; n = n.getParent()) {
                    if (filterTaxid.equals(n.getTaxId())) {
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    continue;
                }
            }
            String taxId = taxNode.getTaxId();
            Long l = kuTaxid2KMer.get(taxId);
            long h = l == null ? 0 : l;
            long g = stats.getLong(taxId);
            long diff = Math.abs(h - g);
            if (l == null && g > 0) {
                missingNodesInKu.add(taxNode);
            } else {
                inDataErr += diff;
                inDataEntries++;
            }
            if (l != null && diff > 0) {
                differences.put(taxNode, diff);
            }
            err += diff;
            entries++;

            if (l != null || full) {
    //            if (taxNode.getRank().equals(Rank.GENUS) || taxNode.getRank().isBelow(Rank.GENUS)) {
                    out.print(taxNode.getTaxId());
                    out.print(';');
                    out.print(getRankString(taxNode));
                    out.print(';');
                    out.print(LF.format(correctDBValue(g, full)));
                    out.print(';');
                    out.print(LF.format(correctDBValue(h, full)));
                    out.println(';');
    //            }
            }
            kuTaxid2KMer.remove(taxId);
            SmallTaxTree.SmallTaxIdNode genusNode = taxNode;
            for (; genusNode != null && genusNode.getRank() != null && genusNode.getRank().isBelow(Rank.GENUS); genusNode = genusNode.getParent()) {
            }
            if (genusNode != null && Rank.GENUS.equals(genusNode.getRank())) {
                Long v1 = kmerSumsPerGenusG.get(genusNode.getTaxId());
                kmerSumsPerGenusG.put(genusNode.getTaxId(), v1 == null ? g : v1 + g);
                Long t1 = kmerSumsPerGenusG.get("TOTAL");
                kmerSumsPerGenusG.put("TOTAL", t1 == null ? g : t1 + g);
                Long v2 = kmerSumsPerGenusKU.get(genusNode.getTaxId());
                kmerSumsPerGenusKU.put(genusNode.getTaxId(), v2 == null ? h : v2 + h);
                Long t2 = kmerSumsPerGenusKU.get("TOTAL");
                kmerSumsPerGenusKU.put("TOTAL", t2 == null ? h : t2 + h);
            }
        }

        System.out.println("Genestrip:");
        System.out.println(kmerSumsPerGenusG);
        System.out.println("KU:");
        System.out.println(kmerSumsPerGenusKU);
        double macroAverageSum = 0;
        int count = 0;
        for (String key : kmerSumsPerGenusKU.keySet()) {
            Long h = kmerSumsPerGenusKU.get(key);
            Long g = kmerSumsPerGenusG.get(key);
            if (g != null && h != null) {
                System.out.println(key + ": " + ((double) g) / h);
            }
            if (!"TOTAL".equals(key)) {
                macroAverageSum += ((double) g) / h;
                count++;
            }
        }
        System.out.println("MacroAverage: " + macroAverageSum / count);

        System.out.println("Absolute error: " + err);
        System.out.println("Entries: " + entries);
        System.out.println("Absolute in data error: " + inDataErr);
        System.out.println("In data entries: " + inDataEntries);
        System.out.println("Mean absolute error: " + ((double) err) / entries);
        System.out.println("Mean absolute in data error: " + ((double) inDataErr) / inDataEntries);
        System.out.println("Tax ids of genomes in KrakenUniq but not in Genestrip:");
        System.out.println(kuTaxid2KMer.size());
        System.out.println(kuTaxid2KMer);
        System.out.println("Tax ids of genomes in Genestrip but not in KrakenUniq:");
        System.out.println(missingNodesInKu.size());
        System.out.println(missingNodesInKu);
        System.out.println("Differences:");
        System.out.println(differences.size());
        System.out.println(differences);

        for (SmallTaxTree.SmallTaxIdNode missingInKU : missingNodesInKu) {
            if (Rank.SPECIES.equals(missingInKU.getRank()) || missingInKU.getRank().isBelow(Rank.SPECIES) || missingInKU.getSubNodes() == null || missingInKU.getSubNodes().length == 0) {
                System.out.println(missingInKU.getTaxId());
            }
        }
    }

    public void writeTaxidsTxtFromKUDB(Rank maxRank, String krakenDB, String genestripDB) throws IOException {
        Map<String, Long> kuTaxid2KMer = getKrakenDBCounts(getKrakenCountsFile(krakenDB));

        GSCommon config = new GSCommon(baseDir);
        GSProject project = new GSProject(config, genestripDB, null, null, null, null, null, null,
                null, null, null, false);
        File taxidsTxt = new File(project.getProjectDir(), "taxids.txt");

        GSMaker maker = new GSMaker(project);
        TaxTree taxTree = ((ObjectGoal<TaxTree, GSProject>) maker.getGoal(GSGoalKey.TAXTREE)).get();

        try (PrintStream ps = new PrintStream(new FileOutputStream(taxidsTxt))) {
            for (String taxid : kuTaxid2KMer.keySet()) {
                TaxTree.TaxIdNode node = taxTree.getNodeByTaxId(taxid);
                if (node != null && (maxRank.equals(node.getRank()) || node.getRank().isBelow(maxRank))) {
                    ps.println(taxid);
                }
            }
        }

        maker.dumpAll();
    }
}
