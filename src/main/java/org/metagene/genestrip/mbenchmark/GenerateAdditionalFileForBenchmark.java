package org.metagene.genestrip.mbenchmark;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.metagene.genestrip.GSCommon;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSMaker;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.genbank.AssemblySummaryReader;
import org.metagene.genestrip.io.BufferedLineReader;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.tax.TaxTree;

import java.io.*;
import java.util.*;

public class GenerateAdditionalFileForBenchmark {
    private static final CSVFormat FORMAT = CSVFormat.DEFAULT.builder().setQuote(null).setCommentMarker('#')
            .setDelimiter('\t').setRecordSeparator('\n').build();

    private static final String INPUT_RESOURCE = "nucleotide_database_filelist.txt";

    private final File baseDir;

    public GenerateAdditionalFileForBenchmark(File baseDir) {
        this.baseDir = baseDir;
    }

    public Set<String> generate(String dbName) throws IOException {
        GSCommon config = new GSCommon(baseDir);
        GSProject project = new GSProject(config, dbName, null, null, null, null, null, false, null,
                null, null, null, false);
        GSMaker maker = new GSMaker(project);
        maker.getGoal(GSGoalKey.ASSEMBLYDOWNLOAD).make();

        Set<String> files = new HashSet<String>();
        try (BufferedReader lr = new BufferedReader(new InputStreamReader(getClass().getClassLoader().getResourceAsStream(INPUT_RESOURCE)))) {
            String line;
            while ((line = lr.readLine()) != null) {
                String fileName = line.trim().substring(line.lastIndexOf('/') + 1);
                String name = fileName.substring(0, fileName.indexOf('_', 4));
                files.add(name);
            }
        }
        ;
        Set<String> found = new HashSet<>();
        Set<String> taxids = new HashSet<>();
        try (PrintStream ps = new PrintStream(new File(project.getProjectDir(), "additional.txt"))) {
            try (CSVParser parser = FORMAT
                    .parse(new InputStreamReader(StreamProvider.getInputStreamForFile(new File(project.getCommon().getGenbankDir(), AssemblySummaryReader.ASSEMLY_SUM_GENBANK))))) {
                for (CSVRecord record : parser) {
                    String taxid = record.get(5);
                    String speciesTaxid = record.get(6);
                    String latest = record.get(10);
                    String complete = record.get(11);
                    String name = record.get(17);
                    String ftp = record.get(19);
                    AssemblySummaryReader.FTPEntryQuality quality = AssemblySummaryReader.FTPEntryQuality.fromString(complete, latest);
                    AssemblySummaryReader.FTPEntryWithQuality ewq = new AssemblySummaryReader.FTPEntryWithQuality(ftp, quality, null, false);
                    if (files.contains(name)) {
                        found.add(name);
                        ps.print(taxid);
                        ps.print(' ');
                        ps.print(ewq.getFileName());
                        ps.print(' ');
                        ps.print(ewq.getFtpURL());
                        ps.print('/');
                        ps.println(ewq.getFileName());
                        taxids.add(taxid);
                    }
                }
            }
        }
        try (PrintStream ps = new PrintStream(new File(project.getProjectDir(), "taxids.txt"))) {
            for (String taxid : taxids) {
                ps.println(taxid);
            }
        }
        files.removeAll(found);
        return files;
    }

    public static void main(String[] args) throws IOException {
        Set<String> missing = new GenerateAdditionalFileForBenchmark(new File("./data")).generate("mbenchmark");
        System.out.println("Missing Accessions (" + missing.size() + " Values):");
        System.out.println(missing);
    }
}
