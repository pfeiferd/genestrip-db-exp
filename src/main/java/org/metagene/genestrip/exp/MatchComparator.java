package org.metagene.genestrip.exp;

import org.metagene.genestrip.GSCommon;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSMaker;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.goals.MatchGoal;
import org.metagene.genestrip.goals.MatchResultGoal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.match.CountsPerTaxid;
import org.metagene.genestrip.match.MatchingResult;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.tax.SmallTaxTree;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.HashMap;
import java.util.Locale;
import java.util.Map;
import java.util.SortedMap;

public class MatchComparator {
    private static final DecimalFormat DF = new DecimalFormat("0.00", new DecimalFormatSymbols(Locale.US));

    protected final File baseDir;

    public MatchComparator(File baseDir) {
        this.baseDir = baseDir;
    }

    public void reportComparisonForScatterPlot(File out, boolean unique, Map<String, long[]> results) throws IOException {
        int base = unique ? 1 : 0;
        try (PrintStream ps = new PrintStream(new FileOutputStream(out))) {
            ps.println("kmers 1; kmers 2;");
            for (String taxid : results.keySet()) {
                long[] counts = results.get(taxid);
                ps.print(counts[base + 0]);
                ps.print(';');
                ps.print(counts[base + 2]);
                ps.println(';');
            }
        }
    }

    public String[] getMaxResults(Map<String, long[]> results, int index, int max) {
        String[] maxResults = new String[max];
        long[] maxCounts = new long[max];

        for (String taxid : results.keySet()) {
            long[] counts = results.get(taxid);
            for (int i = 0; i < maxCounts.length; i++) {
                if (counts[index] > maxCounts[i]) {
                    for (int j = maxCounts.length - 1; j > i; j--) {
                        maxResults[j] = maxResults[j - 1];
                        maxCounts[j] = maxCounts[j - 1];
                    }
                    maxCounts[i] = counts[index];
                    maxResults[i] = taxid;
                    break;
                }
            }
        }
        return maxResults;
    }

    public int getRank(String taxid, Map<String, long[]> results, int index) {
        int rank = 0;
        long[] counts = results.get(taxid);
        for (String taxid2 : results.keySet()) {
            long[] counts2 = results.get(taxid2);
            if (counts2[index] > counts[index]) {
                rank++;
            }
        }
        return rank;
    }

    public Map<String, Map<String, long[]>> compareResults(String dbName1, String dbName2, String csvFile, Map<String,String> nameMap) throws IOException {
        Map<String, MatchingResult> matches1 = match(dbName1, csvFile, nameMap);
        Map<String, MatchingResult> matches2 = match(dbName2, csvFile, nameMap);

        Map<String, Map<String, long[]>> allResults = new HashMap<>();

        for (String key : matches1.keySet()) {
            MatchingResult res1 = matches1.get(key);
            MatchingResult res2 = matches2.get(key);
            if (res1.getGlobalStats().getKMers() != res2.getGlobalStats().getKMers()) {
                throw new RuntimeException("Match results do not match");
            }

            Map<String, CountsPerTaxid> stats1 = res1.getTaxid2Stats();
            Map<String, CountsPerTaxid> stats2 = res2.getTaxid2Stats();

            Map<String, long[]> combinedRes = new HashMap<>();
            fillMap(combinedRes, stats1, 0);
            fillMap(combinedRes, stats2, 2);

            allResults.put(key, combinedRes);
        }

        return allResults;
    }

    private void fillMap(Map<String, long[]> combinedRes, Map<String, CountsPerTaxid> stats, int base) {
        for (String taxid : stats.keySet()) {
            CountsPerTaxid c1 = stats.get(taxid);
            long[] point = combinedRes.get(taxid);
            if (point == null) {
                point = new long[4];
                combinedRes.put(taxid, point);
            }
            point[base + 0] = c1.getKMers();
            point[base + 1] = c1.getUniqueKMers();
        }
    }

    public Map<String, MatchingResult> match(String dbName, String csvFile, Map<String,String> nameMap) throws IOException {
        GSCommon config = new GSCommon(baseDir);

        GSProject project = new GSProject(config, dbName, null, null, csvFile, null, null, null,
                null, null, null, false);

        GSMaker maker = new GSMaker(project);

        MatchResultGoal matchGoal = (MatchResultGoal) maker.getGoal(GSGoalKey.MATCHRES);
        matchGoal.cleanThis();
        matchGoal.make();
        Map<String, MatchingResult> matches = matchGoal.get();
        if (nameMap != null) {
            SmallTaxTree tree = ((ObjectGoal<Database, GSProject>) maker.getGoal(GSGoalKey.LOAD_DB)).get().getTaxTree();
            for (String res : matches.keySet()) {
                MatchingResult res1 = matches.get(res);
                for (String taxid : res1.getTaxid2Stats().keySet()) {
                    nameMap.put(taxid, tree.getNodeByTaxId(taxid).getName());
                }
            }
        }
        maker.dumpAll();
        return matches;
    }

    public void writeCompleteReport(String dbName1, String dbName2, String csvFile) throws IOException {
        File outDir = getOutDir(dbName1);
        Map<String,String> nameMap = new HashMap<>();
        Map<String, Map<String, long[]>> allResults = compareResults(dbName1, dbName2, csvFile, nameMap);
        for (String key : allResults.keySet()) {
            Map<String, long[]> results = allResults.get(key);
            reportComparisonForScatterPlot(new File(outDir,  key + "_kmer_scatter.csv"), false, results);
            reportComparisonForScatterPlot(new File(outDir,  key + "_ukmer_scatter.csv"), true, results);

            String[] maxResults = getMaxResults(results, 3, 10);

            File file = new File(outDir,  key + "_kmer_changes.csv");
            try (PrintStream out = new PrintStream(new FileOutputStream(file))) {
                out.println("key; taxid; name; old rank; old unique kmers; new unique kmers; new rank; unique kmer change percent;");
                for (int i = 0; i < maxResults.length; i++) {
                    String taxid = maxResults[i];
                    long[] counts = results.get(taxid);
                    out.print(key);
                    out.print(';');
                    out.print(taxid);
                    out.print(';');
                    String name = nameMap.get(taxid);
                    out.print(name == null ? "" : name);
                    out.print(';');
                    out.print(getRank(taxid, results, 1) + 1);
                    out.print(';');
                    out.print(counts[1]);
                    out.print(';');
                    out.print(counts[3]);
                    out.print(';');
                    out.print(i + 1);
                    out.print(';');
                    out.print(DF.format(100d * (counts[3] - counts[1]) / counts[1]));
                    out.println(';');
                }
            }
        }
    }

    protected File getOutDir(String dbName) throws IOException {
        GSCommon config = new GSCommon(baseDir);
        GSProject project = new GSProject(config, dbName, null, null, null, null, null, null,
                null, null, null, false);
        return project.getResultsDir();
    }
}
