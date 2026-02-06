package org.metagene.genestrip.exp;

import org.junit.Test;

import java.util.Random;

import static org.junit.Assert.assertEquals;

public class ErrCompInfoTest {
    @Test
    public void testErrorStdDev() {
        for (int k = 0; k < 100; k++) {
            GenestripComparator.ErrCompInfo errCompInfo = new GenestripComparator.ErrCompInfo(0,0,0);

            int n = 10000;
            double rdiffrelmean = 0;
            double kdiffrelmean = 0;

            int bound = 10000;

            Random r = new Random(k);
            for (int i = 0; i < n; i++) {
                int k1 = r.nextInt(bound);
                int r1 = r.nextInt(bound);
                int k2 = r.nextInt(bound);
                int r2 = r.nextInt(bound);
                errCompInfo.sumErrorStats(k1, r1, k2, r2);
                rdiffrelmean += ((double) Math.abs(r1 - r2)) / (r1 + 1);
                kdiffrelmean += ((double) Math.abs(k1 - k2)) / (k1 + 1);
            }
            rdiffrelmean /= n;
            kdiffrelmean /= n;

            assertEquals(rdiffrelmean, errCompInfo.getMeanReadsErr(), 0.00001);
            assertEquals(kdiffrelmean, errCompInfo.getMeanKMersErr(), 0.00001);

            double rdiffrelstddev = 0;
            double kdiffrelstddev = 0;
            r = new Random(k);
            for (int i = 0; i < n; i++) {
                int k1 = r.nextInt(bound);
                int r1 = r.nextInt(bound);
                int k2 = r.nextInt(bound);
                int r2 = r.nextInt(bound);
                double rdiffrel = ((double) Math.abs(r1 - r2)) / (r1 + 1);
                double kdiffrel = ((double) Math.abs(k1 - k2)) / (k1 + 1);

                rdiffrelstddev += rdiffrel * rdiffrel;
                kdiffrelstddev += kdiffrel * kdiffrel;
            }
            rdiffrelstddev = Math.sqrt(rdiffrelstddev / (n - 1));
            kdiffrelstddev = Math.sqrt(kdiffrelstddev / (n - 1));

            System.out.println(rdiffrelstddev + " " + errCompInfo.getReadsErrStdDev());
            System.out.println(kdiffrelstddev + " " + errCompInfo.getKMersErrStdDev());
            assertEquals(rdiffrelstddev, errCompInfo.getReadsErrStdDev(), 0.3);
            assertEquals(kdiffrelstddev, errCompInfo.getKMersErrStdDev(), 0.3);
        }
    }
}
