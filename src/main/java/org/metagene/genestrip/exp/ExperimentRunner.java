package org.metagene.genestrip.exp;

import org.metagene.genestrip.kucomp.KrakenDBComparator;
import org.metagene.genestrip.kucomp.KrakenMatchComparator;

import java.io.File;
import java.io.IOException;
import java.util.Map;

public class ExperimentRunner {
    public static void main(String[] args) throws IOException {
        KrakenDBComparator c2 = new KrakenDBComparator(new File("./data"));

        KrakenMatchComparator c1 = new KrakenMatchComparator(new File("./data"));
        c1.writeUnfoldedTaxids("human_virus");

        // Figure 2
        c2.reportKMerComparisons("viral", "viral_db", null, false);
        // Figure 3
        c1.compareWithKUResults("viral", "viral_db", "viral_ku_comp_fasta.txt", "viral_ku_comp.txt", false);
        // Text in context with Figure 3
        c1.accuracyCheckForSimulatedViralReads("viral", "viral_ku_comp.txt");


        // Figure 4
        Map<String, GenestripComparator.ErrCompInfo> res3 = c1.compareWithKUResults("viral", "viral_db", null, "saliva.txt", false);
        // Table 4
        c1.writeErrInfos("viral", "viral_db", res3);

        // Figure 5
        c1.compareCommonDBEntries("viral", "human_virus");
        // Use taxids from Genestrip's human virus database as basis for tax ids.
        c2.reportKrakenDBComparison("human_virus", "viral_db", "human_virus_db");

        // Figure 6
        Map<String, GenestripComparator.ErrCompInfo> res1 = c1.compareResults("viral", "human_virus", "saliva.txt");
        Map<String, GenestripComparator.ErrCompInfo> res2 = c1.compareKUWithKUResults("viral", "human_virus", "saliva.txt");
        // Table 5
        c1.combineErrInfos("viral", "human_virus", res1, res2);

        // Not needed: ?
        //c1.compareWithKUResults("human_virus", null, "saliva.txt");

        // Section "The Databases \texttt{MicrobialDB} and \texttt{tb}"
        c2.reportKMerComparisons("tick-borne", "microbial_db", null, false);
        c2.reportKMerComparisons("tick-borne", "microbial_db", "943", true);

        // Figure 7 and Figure 8
        c1.compareWithKUResults("tick-borne", "microbial_db", null, "ticks.txt", true);
        c1.aggregateCompareWithKUResults("tick-borne", "ticks.txt");
    }
}
