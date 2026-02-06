package org.metagene.genestrip.exp;

import org.metagene.genestrip.kucomp.GanonMatchComparator;
import org.metagene.genestrip.kucomp.KrakenDBComparator;
import org.metagene.genestrip.kucomp.KrakenMatchComparator;
import org.metagene.genestrip.tax.Rank;

import java.io.File;
import java.io.IOException;
import java.util.Map;

public class ExperimentRunner {
    public static void main(String[] args) throws IOException {
        GanonMatchComparator gmc = new GanonMatchComparator(new File("./data"), new File("./results"));
        gmc.accuracyCheckForSimulatedViralReadsGanon("viral", "ganon/virals.all");

        KrakenMatchComparator c1 = new KrakenMatchComparator(new File("./data"), new File("./results"));
        c1.accuracyCheckForSimulatedViralReads("viral", "viral_ku_comp.txt");
        /*
        KrakenDBComparator c2 = new KrakenDBComparator(new File("./data"), new File("./results"));

        KrakenMatchComparator c1 = new KrakenMatchComparator(new File("./data"), new File("./results"));

        c1.writeUnfoldedTaxids("human_virus");

        System.out.println("** Figure 2 **");
        c2.reportKMerComparisons("viral", "viral_db", null, false);

        System.out.println("** Figure 3 **");
        c1.compareWithKUResults("viral", "viral_db", "viral_ku_comp_fasta.txt", "viral_ku_comp.txt", false);

        System.out.println("** Text in context of Figure 3 **");
        c1.accuracyCheckForSimulatedViralReads("viral", "viral_ku_comp.txt");
        c1.accuracyCheckForSimulatedViralReadsKU("viral", "viral_ku_comp.txt");
        System.out.println("** Figure 4 **");
        Map<String, GenestripComparator.ErrCompInfo> res3 = c1.compareWithKUResults("viral", "viral_db", null, "saliva.txt", false);
        System.out.println("** Table 3 **");
        c1.writeErrInfos("viral", "viral_db", res3);

        System.out.println("** Figure 5 **");
        c1.compareCommonDBEntries("viral", "human_virus");
        // Use taxids from Genestrip's human virus database as basis for tax ids.
        c2.reportKrakenDBComparison("human_virus", "viral_db", "human_virus_db");

        System.out.println("** Figure 6 **");
        Map<String, GenestripComparator.ErrCompInfo> res1 = c1.compareResults("viral", "human_virus", "saliva.txt");
        Map<String, GenestripComparator.ErrCompInfo> res2 = c1.compareKUWithKUResults("viral", "human_virus", "saliva.txt");

        System.out.println("** Table 4 **");
        c1.combineErrInfos("viral", "human_virus", res1, res2);

        // Not needed: ?
        //c1.compareWithKUResults("human_virus", null, "saliva.txt");

        System.out.println("** Section \"The Databases \\texttt{MicrobialDB} and \\texttt{tb}\" **");
        c2.reportKMerComparisons("tick-borne", "microbial_db", null, false);
        c2.reportKMerComparisons("tick-borne", "microbial_db", "943", true);

        System.out.println("** Figure 7 and Figure 8 **");
        c1.compareWithKUResults("tick-borne", "microbial_db", null, "seventicks.txt", true);
        c1.aggregateCompareWithKUResults("tick-borne", "ticks.txt");

         */
    }
}
