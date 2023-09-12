package org.seqcode.deepseq.hitloaders;
/*
import java.io.File;
import org.junit.Test;
import static org.junit.jupiter.api.Assertions.assertEquals;
import org.seqcode.genome.Genome;

public class HDF5HitLoaderTest {
	
	ClassLoader classLoader = getClass().getClassLoader();
	
	@Test
	public void testHDF5() throws Exception {
		String bamFile = classLoader.getResource("org/seqcode/deepseq/resources/test_300reads_shuffle.bam").getFile();
		String genFile = classLoader.getResource("org/seqcode/deepseq/resources/mm10.fa.fai").getFile();
		
		Genome gen = new Genome("Genome", new File(genFile), true);
		
		HDF5HitLoader hl = new HDF5HitLoader(gen, new File(bamFile), false, false, false, true, false);
		
		hl.sourceAllHits();
		
		assertEquals(hl.getHitPairInfo().getElement("1", 0, "r1Pos", 1), 134986990);
		assertEquals(hl.getHitPairInfo().getElement("1", 0, "r2Strand", 1), 1);
		
		hl.getHitPairInfo().sortByReference();
		
		assertEquals(hl.getHitPairInfo().getElement("1", 0, "r1Pos", 1), 49974129);
		assertEquals(hl.getHitPairInfo().getElement("1", 0, "r2Strand", 1), 1);
		
		assertEquals(hl.getHitPairInfo().getLength("10", 0), 3);
		
		hl.close();
	}
	
}
*/
