package org.metagene.genestrip.kucomp;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.metagene.genestrip.*;
import org.metagene.genestrip.exp.GenestripComparator;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.refseq.AccessionMap;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.util.ByteArrayUtil;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class GanonMatchComparator extends GenestripComparator {
    private static final CSVFormat GANON_ALL_FORMAT = CSVFormat.DEFAULT.builder().setQuote(null).setCommentMarker('#')
            .setDelimiter('\t').setRecordSeparator('\n').build();

    public GanonMatchComparator(File baseDir, File resultsDir) {
        super(baseDir, resultsDir);
    }

    public void accuracyCheckForSimulatedViralReadsGanon(String db, String ganonReportFile) throws IOException {
        GSCommon config = new GSCommon(baseDir);
        GSProject project = new GSProject(config, db, null, null, null, null, null, null,
                null, null, null, false);
        project.initConfigParam(GSConfigKey.THREADS, -1);

        GSMaker maker = new GSMaker(project);

        ObjectGoal<TaxTree, GSProject> taxTreeGoal = (ObjectGoal<TaxTree, GSProject>) maker.getGoal(GSGoalKey.TAXTREE);
        TaxTree taxTree = taxTreeGoal.get();

        Map<String, TaxTree.TaxIdNode> matchesMap = new HashMap<String, TaxTree.TaxIdNode>();
        try (CSVParser parser = GANON_ALL_FORMAT
                .parse(new InputStreamReader(new FileInputStream(new File(ganonReportFile))))) {
            // Ganon does multiple matches.
            // We reduce it to a single match by computing the LCA in case of the multiple matches.
            for (CSVRecord record : parser.getRecords()) {
                String descr = record.get(0);
                String ganonTaxid = record.get(1);
                if (ganonTaxid != null) {
                    TaxTree.TaxIdNode classNode = taxTree.getNodeByTaxId(ganonTaxid);
                    if (classNode != null) {
                        TaxTree.TaxIdNode match = matchesMap.get(descr);
                        if (match != null) {
                            match = taxTree.getLeastCommonAncestor(match, classNode);
                        } else {
                            match = classNode;
                        }
                        matchesMap.put(descr, match);
                    }
                    else {
                        System.out.println("Invalid taxid from Ganon: " + ganonTaxid);
                    }
                }
            }
        }

        ObjectGoal<AccessionMap, GSProject> accessCollGoal = (ObjectGoal<AccessionMap, GSProject>) maker.getGoal(GSGoalKey.ACCMAP);
        AccessionMap map = accessCollGoal.get();

        int[] counters = new int[5];
        for (String descr : matchesMap.keySet()) {
            byte[] desc = descr.getBytes();
            int pos = ByteArrayUtil.indexOf(desc, 5, desc.length, '_');
            TaxTree.TaxIdNode node = map.get(desc, 1, pos, false);
            if (node != null) {
                TaxTree.TaxIdNode classNode = matchesMap.get(descr);
                TaxTree.TaxIdNode lca = taxTree.getLeastCommonAncestor(node, classNode);
                while (lca != null && Rank.NO_RANK.equals(lca.getRank())) {
                    lca = lca.getParent();
                }
                if (lca != null) {
                    Rank r = lca.getRank();
                    if (r != null) {
                        if (Rank.GENUS.equals(r) || r.isBelow(Rank.GENUS)) {
                            counters[1]++;
                        }
                        if (Rank.SPECIES.equals(r) || r.isBelow(Rank.SPECIES)) {
                            counters[2]++;
                        }
                        if (Rank.STRAIN.equals(r) || r.isBelow(Rank.STRAIN)) {
                            counters[3]++;
                        }
                    }
                }
            }
            else {
                counters[4]++;
            }
            counters[0]++;
        }
        System.out.println("+++ Ganon +++");
        System.out.println("Correct classifications GENUS: " + counters[1]);
        System.out.println("Correct classifications SPECIES: " + counters[2]);
        System.out.println("Correct classifications STRAIN: " + counters[3]);
        System.out.println("Bad Ground Truth Descriptor: " + counters[4]);
        System.out.println("Total classified: " + counters[0]);
        maker.dumpAll();
    }
}
