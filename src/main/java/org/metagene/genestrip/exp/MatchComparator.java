package org.metagene.genestrip.exp;

import org.metagene.genestrip.GSCommon;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSMaker;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.goals.MatchGoal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.match.CountsPerTaxid;
import org.metagene.genestrip.match.MatchingResult;
import org.metagene.genestrip.store.Database;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Map;
import java.util.SortedMap;

public class MatchComparator {
    private final File baseDir;

    public MatchComparator(File baseDir) {
        this.baseDir = baseDir;
    }

    public void reportComparisonForScatterPlot(File out, boolean unique, Map<String, long[]> results) throws IOException {
        int base = unique ? 0 : 2;
        try (PrintStream ps = new PrintStream(new FileOutputStream(out))) {
            for (String taxid : results.keySet()) {
                long[] counts = results.get(taxid);
                ps.print(base + counts[0]);
                ps.print('\t');
                ps.print(base + counts[1]);
                ps.println();
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
                    for (int j = i; j < maxCounts.length - 1; j++) {
                        maxResults[j + 1] = maxResults[j];
                        maxCounts[j + 1] = maxCounts[j];
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

    public Map<String, Map<String, long[]>> compareResults(String dbName1, String dbName2, String csvFile) throws IOException {
        Map<String, MatchingResult> matches1 = match(dbName1, csvFile);
        Map<String, MatchingResult> matches2 = match(dbName1, csvFile);

        Map<String, Map<String, long[]>> allResults = new HashMap<>();

        for (String key : matches1.keySet()) {
            MatchingResult res1 = matches1.get(key);
            MatchingResult res2 = matches2.get(key);
            if (res1.getTotalKMers() != res2.getTotalKMers()) {
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

    public Map<String, MatchingResult> match(String dbName, String csvFile) throws IOException {
        GSCommon config = new GSCommon(baseDir);

        GSProject project = new GSProject(config, dbName, null, null, csvFile, null, null, false, "64320,12637+",
                null, null, null, false);

        GSMaker maker = new GSMaker(project);

        MatchGoal matchGoal = (MatchGoal) maker.getGoal(GSGoalKey.MATCHLR);
        matchGoal.make();
        return matchGoal.getMatchResults();
    }

    public void writeCompleteReport(File baseDir, String dbName1, String dbName2, String csvFile) throws IOException {
        MatchComparator comparator = new MatchComparator(baseDir);

        GSCommon config = new GSCommon(baseDir);

        GSProject project = new GSProject(config, dbName1, null, null, csvFile, null, null, false, "64320,12637+",
                null, null, null, false);
        File outDir = project.getKrakenOutDir();

        Map<String, Map<String, long[]>> allResults = comparator.compareResults(dbName1, dbName2, csvFile);
        for (String key : allResults.keySet()) {
            Map<String, long[]> results = allResults.get(key);
            comparator.reportComparisonForScatterPlot(new File(outDir,  key + "_kmer_scatter.out"), false, results);
            comparator.reportComparisonForScatterPlot(new File(outDir,  key + "_ukmer_scatter.out"), true, results);

            String[] maxResults = comparator.getMaxResults(results, 1, 10);

            // TODO: Write this to a file...
            for (int i = 0; i < maxResults.length; i++) {
                String taxid = maxResults[i];
                long[] counts = results.get(taxid);
                System.out.println("Key: " + key);
                System.out.println("Taxid: " + taxid);
                System.out.println("Old Rank: " + i);
                System.out.println("Old unique k-mers: " + counts[1]);
                System.out.println("New unique k-mers: " + counts[3]);
                System.out.println("New Rank: " + comparator.getRank(taxid, results, 3));
                System.out.println("Unique k-mers change %: " + 100d * counts[3] / counts[1]);
                System.out.println();
            }
        }
    }
}
