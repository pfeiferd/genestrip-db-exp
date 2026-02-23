package org.metagene.genestrip.exp;

import org.metagene.genestrip.*;

import java.io.File;
import java.io.IOException;

public class GetViralFastq {
    public static void main(String[] args) throws IOException {
        GSCommon config = new GSCommon(new File("./data"));
        GSProject project = new GSProject(config, "viral", null, null, "viral_ku_comp_fasta.txt", null, null, null,
                null, null, null, false);
        project.initConfigParam(GSConfigKey.ALWAYS_ASSUME_GZIP, false);
        GSMaker maker = new GSMaker(project);
        maker.getGoal(GSGoalKey.FASTA2FASTQ).make();
        maker.dumpAll();
    }
}
