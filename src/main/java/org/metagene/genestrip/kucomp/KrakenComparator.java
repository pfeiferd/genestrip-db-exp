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
        Object2LongMap<String> stats = database.getStats();

        out.println("name;rank;taxid;genestrip stored kmers; ku stored kmers;");

        /*
        out.print("TOTAL;");
        out.print(Rank.NO_RANK);
        out.print(';');
        out.print("1;");
        out.print(stats.getLong(null));
        out.println(';');
        out.println('0;');
         */

        List<String> sortedTaxIds = new ArrayList<>(stats.keySet());
        SmallTaxTree taxTree = database.getTaxTree();
        taxTree.sortTaxidsViaTree(sortedTaxIds);

        for (String taxId : sortedTaxIds) {
            if (taxId != null) {
                SmallTaxTree.SmallTaxIdNode taxNode = taxTree.getNodeByTaxId(taxId);
                if (taxNode != null) {
                    if (taxNode.getRank().isBelow(Rank.GENUS)) {
                        out.print(taxNode.getName());
                        out.print(';');
                        out.print(taxNode.getRank());
                        out.print(';');
                        out.print(taxNode.getTaxId());
                        out.print(';');
                        out.print(stats.getLong(taxId));
                        out.print(';');
                        Long l = kuTaxid2KMer.get(taxId);
                        out.print(l == null ? "0" : l);
                        out.println(';');
                    }
                }
            }
        }
    }

    public static void main(String[] args) throws IOException {
        new KrakenComparator(new File("./data")).reportKMerComparisons("human_virus", false, "viral_db");
    }
}
