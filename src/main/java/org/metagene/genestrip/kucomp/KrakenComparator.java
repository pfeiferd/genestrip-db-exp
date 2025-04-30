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

    public void reportKrakenDBComparison(String genestripDB, boolean temp, String krakenDB1, String krakenDB2) throws IOException {
        SmallTaxTree tree = getDatabase(genestripDB, false).getTaxTree();
        Map<String, Long> kuTaxid2KMer1 = getKrakenDBCounts(getKrakenCountsFile(krakenDB1));
        Map<String, Long> kuTaxid2KMer2 = getKrakenDBCounts(getKrakenCountsFile(krakenDB1));

        File out = new File(baseDir, krakenDB1 + "_" + krakenDB2 + "_kudbcomp.csv");
        try (PrintStream ps = new PrintStream(new FileOutputStream(out))) {
            ps.println("taxid; rank; " + krakenDB1 + " kmers; " + krakenDB2 + " kmers;");
            for (String key : kuTaxid2KMer1.keySet()) {
                Long count1 = kuTaxid2KMer2.get(key);
                if (count1 != null) {
                    SmallTaxTree.SmallTaxIdNode node = tree.getNodeByTaxId(key);
                    if (node != null) {
                        Long count2 = kuTaxid2KMer2.get(key);
                        if (count2 != null) {

                        }
                        ps.print(key);
                        ps.print(";");
                        boolean norank = false;
                        while (node != null && Rank.NO_RANK.equals(node.getRank())) {
                            norank = true;
                            node = node.getParent();
                        }
                        Rank r = node == null ? null : node.getRank();
                        if (r != null && (r.isBelow(Rank.SPECIES) ||
                                (norank && r.equals(Rank.SPECIES)))) {
                            ps.print("species or below");
                        } else if (r != null && r.equals(Rank.SPECIES)) {
                            ps.print("species or below");
                        } else if (r != null && r.equals(Rank.GENUS)) {
                            ps.print("genus");
                        } else {
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

    protected long correctDBValue(long v) {
        return v + 1;
    }

    public void reportKMerComparisons(String genestripDB, boolean temp, String krakenDB) throws IOException {
        Map<String, Long> kuTaxid2KMer = getKrakenDBCounts(getKrakenCountsFile(krakenDB));

        String file1 = "gku_kmer_counts.txt";
        File countsFile = new File(getOutDir(genestripDB), "gku_kmer_counts.csv");
        try (PrintStream out = new PrintStream(new FileOutputStream(countsFile))) {
            printJointStoreInfo(getDatabase(genestripDB, temp), out, kuTaxid2KMer);
        }
    }

    protected Map<String, Long> getKrakenDBCounts(File file) throws IOException {
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
        return kuTaxid2KMer;
    }

    protected File getKrakenCountsFile(String krakenDB) {
        return new File(baseDir, "../ku/" + krakenDB + "/database.kdb.counts");
    }

    public void printJointStoreInfo(Database database, PrintStream out, Map<String, Long> kuTaxid2KMer) {

        out.println("name;rank;taxid;genestrip stored kmers;ku stored kmers;");

        visit(database.getTaxTree().getRoot(), database.getStats(), out, kuTaxid2KMer);

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
        if (taxNode == null || taxNode.getRank() == null /* || !taxNode.getRank().isBelow(Rank.GENUS)*/) {
            return;
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
