package org.metagene.genestrip.kucomp;

import org.metagene.genestrip.GSCommon;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSMaker;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.exp.MatchComparator;
import org.metagene.genestrip.goals.MatchResultGoal;
import org.metagene.genestrip.goals.kraken.KrakenResCountGoal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.match.CountsPerTaxid;
import org.metagene.genestrip.match.MatchingResult;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.List;
import java.util.Map;

public class KrakenMatchComparator {
    private final File baseDir;

    public KrakenMatchComparator(File baseDir) {
        this.baseDir = baseDir;
    }

    public void compareWithKUResults(String db, String csvFile1, String csvFile2, String[] keys) throws IOException {
        GSCommon config = new GSCommon(baseDir);
        GSProject project = new GSProject(config, db, null, null, csvFile1, null, null, false, null,
                null, null, null, false);
        GSMaker maker = new GSMaker(project);
        maker.getGoal(GSGoalKey.FASTA2FASTQ).make();

        GSProject project2 = new GSProject(config, db, null, null, csvFile2, null, null, false, null,
                null, null, null, false);

        GSMaker maker2 = new GSMaker(project2);
        ObjectGoal<Map<String, MatchingResult>, GSProject> matchResGoal = (ObjectGoal<Map<String, MatchingResult>, GSProject>) maker2.getGoal(GSGoalKey.MATCHRES);
        Map<String, MatchingResult> matchResult = matchResGoal.get();

        ObjectGoal<Map<String, List<KrakenResCountGoal.KrakenResStats>>, GSProject> countGoal =
                (ObjectGoal<Map<String, List<KrakenResCountGoal.KrakenResStats>>, GSProject>) maker2.getGoal(GSGoalKey.KRAKENCOUNT);
        Map<String, List<KrakenResCountGoal.KrakenResStats>> stats = countGoal.get();

        List<KrakenResCountGoal.KrakenResStats> list = stats.get(keys[0]);
        MatchingResult match = matchResult.get(keys[0]);

        File out = new File(project.getResultsDir(), keys[0] + ".csv");
        try (PrintStream ps = new PrintStream(new FileOutputStream(out))) {
            ps.println("taxid; genestrip kmers; ku kmers;");
            for (KrakenResCountGoal.KrakenResStats kustats : list) {
                CountsPerTaxid gcounts = match.getTaxid2Stats().get(kustats.getTaxid());
                ps.print(kustats.getTaxid());
                ps.print(';');
                ps.print(gcounts.getKMers());
                ps.print(';');
                ps.print(kustats.getKmers());
                ps.println(';');
            }
        }
    }

    public static void main(String[] args) throws IOException {
        KrakenMatchComparator c = new KrakenMatchComparator(new File("./data"));
        c.compareWithKUResults("viral", "viral_ku_comp_fasta.txt", "viral_ku_comp.txt", new String[] { "fastq1" });
    }
}
