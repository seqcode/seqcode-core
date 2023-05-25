package org.seqcode.deepseq.experiments;

import org.junit.Test;
import org.junit.Assert;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.seqcode.deepseq.StrandedBaseCount;
import org.seqcode.deepseq.StrandedPair;
import org.seqcode.deepseq.hitloaders.HDF5HitLoader;
import org.seqcode.deepseq.hitloaders.HitLoader;
import org.seqcode.genome.Genome;
import org.seqcode.genome.location.Region;

public class HDF5HitCacheTest {
	
	ClassLoader classLoader = getClass().getClassLoader();
	
	@Test
	public void testHDF5Cache() throws Exception {
		String bam1File = classLoader.getResource("org/seqcode/deepseq/resources/test_300reads_shuffle.bam").getFile();
		String bam2File = classLoader.getResource("org/seqcode/deepseq/resources/test_300reads_shuffle2.bam").getFile();
		String genFile = classLoader.getResource("org/seqcode/deepseq/resources/mm10.fa.fai").getFile();
		
		Genome gen = new Genome("Genome", new File(genFile), true);
		ExptConfig ec = new ExptConfig(gen, new String[] {"--loadpairs"});
		ec.keepHDF5 = true;
		
		HDF5HitLoader hl1 = new HDF5HitLoader(gen, new File(bam1File), true, true, true, true, false);
		HDF5HitLoader hl2 = new HDF5HitLoader(gen, new File(bam2File), true, true, true, true, false);
		
		List<HitLoader> hList = new ArrayList<HitLoader>() {
			{
				add(hl1);
				add(hl2);
			}
		};
		
		long start = System.currentTimeMillis();
		HDF5HitCache hc = new HDF5HitCache(ec, hList, "test");
		long end = System.currentTimeMillis();
		System.err.println((end - start) + "ms");
		
		List<StrandedPair> pairList = hc.getPairsByMid(new Region(gen, "1", 100000000, 150000000));
		int[] r1PosLists = new int[8];
		int index=0;
		for (StrandedPair sp: pairList) {
			r1PosLists[index] = sp.getR1Coordinate();
			index++;
		}
		Assert.assertArrayEquals(r1PosLists, new int[] {134078910, 134485871, 134986990, 135744005, 128874656, 131897632, 133366982, 134421875});
		
		hc.close();
	}
	
}
