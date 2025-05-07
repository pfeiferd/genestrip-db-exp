package org.metagene.genestrip.kucomp;

import it.unimi.dsi.fastutil.objects.Object2LongMap;
import org.metagene.genestrip.exp.GenestripComparator;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.tax.SmallTaxTree;

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
                String key = node.getTaxId();
                // We only report on tax ids which are in respective Genestrip db:
                if (node != null) {
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
                    ps.print(correctDBValue(count1));
                    ps.print(";");
                    ps.print(correctDBValue(count2));
                    ps.println(";");
                }
            }
        }
        System.out.println(diff);
    }

    public void reportKMerComparisons(String genestripDB, String krakenDB) throws IOException {
        Map<String, Long> kuTaxid2KMer = getKrakenDBCounts(getKrakenCountsFile(krakenDB));

        File countsFile = new File(baseDir, genestripDB + "_gs_ku_dbcomp.csv");
        try (PrintStream out = new PrintStream(new FileOutputStream(countsFile))) {
            printJointStoreInfo(getDatabase(genestripDB, false), out, kuTaxid2KMer);
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

    public void printJointStoreInfo(Database database, PrintStream out, Map<String, Long> kuTaxid2KMer) {
        out.println("taxid; rank; kmers 1; kmers 2;");

        long inDataErr = 0;
        long inDataEntries = 0;
        long err = 0;
        long entries = 0;
        List<SmallTaxTree.SmallTaxIdNode> missingNodesInKu = new ArrayList<>();
        Map<SmallTaxTree.SmallTaxIdNode, Long> differences = new HashMap<>();

        Object2LongMap<String> stats = database.getStats();
        for (SmallTaxTree.SmallTaxIdNode taxNode : database.getTaxTree()) {
            if (taxNode == null || taxNode.getRank() == null) {
                continue;
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

            if (l != null) {
                out.print(taxNode.getTaxId());
                out.print(';');
                out.print(getRankString(taxNode));
                out.print(';');
                out.print(correctDBValue(g));
                out.print(';');
                out.print(correctDBValue(h));
                out.println(';');
            }
            kuTaxid2KMer.remove(taxId);
        }

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
    }
}
