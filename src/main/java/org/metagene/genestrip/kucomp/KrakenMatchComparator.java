package org.metagene.genestrip.kucomp;

import org.metagene.genestrip.GSCommon;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSMaker;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.exp.MatchComparator;
import org.metagene.genestrip.goals.MatchResultGoal;

import java.io.File;
import java.io.IOException;

public class KrakenMatchComparator  {
    private final File baseDir;

    public KrakenMatchComparator(File baseDir) {
        this.baseDir = baseDir;
    }

    public void compareWithKUResults(String db, String csvFile1, String csvFile2) throws IOException {
        GSCommon config = new GSCommon(baseDir);
        GSProject project = new GSProject(config, db, null, null, csvFile1, null, null, false, null,
                null, null, null, false);
        GSMaker maker = new GSMaker(project);
        maker.getGoal(GSGoalKey.FASTA2FASTQ).make();

        GSProject project2 = new GSProject(config, db, null, null, csvFile2, null, null, false, null,
                null, null, null, false);

        GSMaker maker2 = new GSMaker(project2);
        maker2.getGoal(GSGoalKey.MATCH).make();
        maker2.getGoal(GSGoalKey.KRAKENRES).make();
    }

    public static void main(String[] args) throws IOException {
        KrakenMatchComparator c = new KrakenMatchComparator(new File("./data"));
        c.compareWithKUResults("viral", "viral_ku_comp_fasta.txt", "viral_ku_comp.txt");
    }
}
