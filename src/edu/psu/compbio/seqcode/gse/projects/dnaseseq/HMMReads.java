package edu.psu.compbio.seqcode.gse.projects.dnaseseq;

import java.io.*;
import java.util.*;
import cern.jet.random.Binomial;
import cern.jet.random.Normal;
import cern.jet.random.engine.DRand;
import edu.psu.compbio.seqcode.gse.datasets.general.*;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.*;
import edu.psu.compbio.seqcode.gse.projects.readdb.Client;
import edu.psu.compbio.seqcode.gse.projects.readdb.ClientException;
import edu.psu.compbio.seqcode.gse.utils.Closeable;

public class HMMReads implements Closeable {

    private Client client;
    private double subsample;
    private int smooth;

    public HMMReads() throws ClientException, IOException {
        client = new Client();
        subsample = 2;
        smooth = 0;
    }
    public void subSample(double d) {
        subsample = d;
    }
    public double subSample() {return subsample;}
    public int smooth() {return smooth;}
    public void smooth(int s) {
        smooth = s;
    }
    public ReadCounts getReadCounts(Region region,
                                    Collection<SeqAlignment> fgalignments,
                                    Collection<SeqAlignment> bgalignments) throws IOException, ClientException {
        ReadCounts fg = getReadCounts(region, fgalignments);
        ReadCounts bg = getReadCounts(region, bgalignments);
        int fgc[] = fg.getCounts();
        int bgc[] = bg.getCounts();
        for (int i = 0; i < fgc.length; i++) {
            fgc[i] = (int)Math.max(0, fgc[i] - bgc[i]);
        }
        return new ReadCounts(fgc, region);
    }


    public ReadCounts getReadCounts(Region region,
                                    Collection<SeqAlignment> alignments) throws IOException, ClientException {
        int output[] = new int[region.getWidth()];
        int regionstart = region.getStart();
        for (int i = 0; i < output.length; i++) {
            output[i] = 0;
        }
        for (SeqAlignment a : alignments) {
            TreeMap<Integer,Integer> m = client.getHistogram(Integer.toString(a.getDBID()),
                                                             region.getGenome().getChromID(region.getChrom()),
                                                             false,
                                                             false,
                                                             1, // binsize
                                                             10, // dedup
                                                             region.getStart(),
                                                             region.getEnd(),
                                                             null,
                                                             null,
                                                             true);
            for (int pos : m.keySet()) {
                output[pos - regionstart] += m.get(pos);
            }
        }
        if (subsample < 1) {
            Binomial binomial = new Binomial(5, .5, new DRand());
            int newout[] = new int[output.length];
            for (int i = 0; i < newout.length; i++) {
                if (output[i] > 0) {
                    newout[i] = binomial.nextInt(output[i], subsample);
                } else {
                    newout[i] = 0;
                }



            }
            output = newout;
        }
        if (smooth > 0) {
            int newout[] = new int[output.length];
            for (int i = 0; i < output.length; i++) {
                for (int j = 0; j < output[i]; j++) {
                    int pos = (int)(i + Normal.staticNextDouble(0, smooth));
                    if (pos >=0 && pos < newout.length) {
                        newout[pos]++;
                    }
                }
            }
            output = newout;
        }
        return new ReadCounts(output,region);
    }

    public void close() {
        try {
            client.close();
            client = null;
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    public boolean isClosed() {
        return client == null;
    }


}