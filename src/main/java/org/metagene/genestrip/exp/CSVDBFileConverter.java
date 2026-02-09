package org.metagene.genestrip.exp;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.metagene.genestrip.*;
import org.metagene.genestrip.goals.refseq.ExtractRefSeqCSVGoal;

import java.io.*;

public class CSVDBFileConverter {
    protected static final CSVFormat CSV_FORMAT = CSVFormat.DEFAULT.builder().setQuote(null).setCommentMarker('#')
            .setDelimiter(';').setRecordSeparator('\n').build();

    private final File baseDir;

    public CSVDBFileConverter(File file) {
        this.baseDir = file;
    }

    public void csv2GanonInputFileFormat(String db, String pathPrefix) throws IOException {
        GSCommon config = new GSCommon(baseDir);
        GSProject project = new GSProject(config, db, null, null, null, null, null, null,
                null, null, null, false);
        project.initConfigParam(GSConfigKey.THREADS, -1);

        GSMaker maker = new GSMaker(project);

        ExtractRefSeqCSVGoal extractRefSeqCSVGoal = (ExtractRefSeqCSVGoal) maker.getGoal(GSGoalKey.EXTRACT_REFSEQ_CSV);
        extractRefSeqCSVGoal.make();

        try (CSVParser parser = CSV_FORMAT
                .parse(new InputStreamReader(new FileInputStream(extractRefSeqCSVGoal.getFile())))) {
            File ganonInputFile = new File(project.getResultsDir(), db + "_ganon.tsv");
            try (PrintStream out = new PrintStream(new FileOutputStream(ganonInputFile))) {
                int i = 0;
                for (CSVRecord record : parser.getRecords()) {
                    if (i > 0) {
                        String descr = record.get(0);
                        String taxid = record.get(1);
                        out.print(pathPrefix);
                        out.print(descr);
                        out.print(".fa.gz");
                        out.print('\t');
                        out.print(descr);
                        out.print('\t');
                        out.print(taxid);
                        out.println();
                    }
                    i++;
                }
            }
        }
    }


    // This format also works for KrakenUniq to build the custom database
    public void csv2KuMapFileFormat(String db) throws IOException {
        GSCommon config = new GSCommon(baseDir);
        GSProject project = new GSProject(config, db, null, null, null, null, null, null,
                null, null, null, false);
        project.initConfigParam(GSConfigKey.THREADS, -1);

        GSMaker maker = new GSMaker(project);

        ExtractRefSeqCSVGoal extractRefSeqCSVGoal = (ExtractRefSeqCSVGoal) maker.getGoal(GSGoalKey.EXTRACT_REFSEQ_CSV);
        extractRefSeqCSVGoal.make();

        try (CSVParser parser = CSV_FORMAT
                .parse(new InputStreamReader(new FileInputStream(extractRefSeqCSVGoal.getFile())))) {
            File ganonInputFile = new File(project.getResultsDir(), db + "_ku.map");
            try (PrintStream out = new PrintStream(new FileOutputStream(ganonInputFile))) {
                int i = 0;
                for (CSVRecord record : parser.getRecords()) {
                    if (i > 0) {
                        String descr = record.get(0);
                        String taxid = record.get(1);
                        out.print(descr);
                        out.print('\t');
                        out.print(taxid);
                        out.println();
                    }
                    i++;
                }
            }
        }
    }

    public static void main(String[] args) throws IOException {
        CSVDBFileConverter converter = new CSVDBFileConverter(new File("./data"));

        String[] dbs = new String[] { "human_virus", "viral" };
        for (int i = 0; i < dbs.length; i++) {
            converter.csv2GanonInputFileFormat(dbs[i], "../data/projects/" + dbs[i] + "/fasta/");
            converter.csv2KuMapFileFormat(dbs[i]);
        }
    }

}
