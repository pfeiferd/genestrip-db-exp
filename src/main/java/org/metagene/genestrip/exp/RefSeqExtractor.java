package org.metagene.genestrip.exp;

import org.metagene.genestrip.kucomp.GanonMatchComparator;

import java.io.File;
import java.io.IOException;

public class RefSeqExtractor {
    public static void main(String[] args) throws IOException {
        GenestripComparator gc = new GenestripComparator(new File("./data"), new File("./results"));
        gc.extractRefSeqFastas(args[0]);
    }
}
