package org.metagene.genestrip.kucomp;

import it.unimi.dsi.fastutil.objects.Object2LongMap;
import org.metagene.genestrip.exp.DatabaseComparator;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree;

import java.io.*;
import java.util.*;

public class KrakenComparator extends DatabaseComparator {
    public KrakenComparator(File baseDir) {
        super(baseDir);
    }

    public void reportKMerComparisons(String genestripDB, boolean temp, String krakenDB) throws IOException {
        File file = getKrakenCountsFile(krakenDB);
        String file1 = "gku_kmer_counts.txt";
        Map<String, Long> kuTaxid2KMer = new HashMap<String, Long>();
        try (BufferedReader br = new BufferedReader(new FileReader(file))) {
            String line = br.readLine();
            while (line != null) {
                int tab = line.indexOf('\t');
                if (tab != -1) {
                    String taxid = line.substring(0, tab);
                    Long kmers = Long.valueOf(line.substring(tab + 1, line.length()));
                    kuTaxid2KMer.put(taxid, kmers);
                }
                line = br.readLine();
            }
        }
        File countsFile = new File(getOutDir(genestripDB), "gku_kmer_counts.csv");

        try (PrintStream out = new PrintStream(new FileOutputStream(countsFile))) {
            printJointStoreInfo(getDatabase(genestripDB, temp), out, kuTaxid2KMer);
        }
    }

    protected File getKrakenCountsFile(String krakenDB) {
        return new File(baseDir, "../ku/" + krakenDB + "/database.kdb.counts");
    }

    public void printJointStoreInfo(Database database, PrintStream out, Map<String, Long> kuTaxid2KMer) {

        out.println("name;rank;taxid;genestrip stored kmers;ku stored kmers;");

        visit(database.getTaxTree().getNodeByTaxId("1"), database.getStats(), out, kuTaxid2KMer);

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

    protected void visit(SmallTaxTree.SmallTaxIdNode node, Object2LongMap<String> stats, PrintStream out, Map<String, Long> kuTaxid2KMer) {
        handleNode(node, stats, out, kuTaxid2KMer);
        if (node.getSubNodes() != null) {
            for (SmallTaxTree.SmallTaxIdNode subNode : node.getSubNodes()) {
                visit(subNode, stats, out, kuTaxid2KMer);
            }
        }
    }

    private long inDataErr;
    private long inDataEntries;
    private long err;
    private long entries;
    private List<SmallTaxTree.SmallTaxIdNode> missingNodesInKu = new ArrayList<>();
    private Map<SmallTaxTree.SmallTaxIdNode, Long> differences = new HashMap<>();

    protected void handleNode(SmallTaxTree.SmallTaxIdNode taxNode, Object2LongMap<String> stats, PrintStream out, Map<String, Long> kuTaxid2KMer) {
        if (taxNode == null || taxNode.getRank() == null || !taxNode.getRank().isBelow(Rank.GENUS)) {
            return;
        }
        String taxId = taxNode.getTaxId();
        Long l = kuTaxid2KMer.get(taxId);
        long h = l == null ? 0 : l;
        long g = stats.getLong(taxId);
        long diff = Math.abs(h - g);
        if (l == null && g > 0) {
            missingNodesInKu.add(taxNode);
        }
        else {
            inDataErr += diff;
            inDataEntries ++;
        }
        if (l != null && diff > 0) {
            differences.put(taxNode, diff);
        }
        err += diff;
        entries ++;

        out.print(taxNode.getName());
        out.print(';');
        out.print(taxNode.getRank());
        out.print(';');
        out.print(taxNode.getTaxId());
        out.print(';');
        out.print(g);
        out.print(';');
        out.print(h);
        out.println(';');
        kuTaxid2KMer.remove(taxId);
    }

    public static void main(String[] args) throws IOException {
        KrakenComparator c = new KrakenComparator(new File("./data"));
        c.reportKMerComparisons("viral", false, "viral_db");
        //c.reportKMerComparisons("chronicb_std_ku", false, "standard_db");
    }
}
