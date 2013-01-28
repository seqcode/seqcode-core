package edu.psu.compbio.seqcode.gse.utils.stats;

public class TestStat {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		double v = StatUtil.hyperGeometricCDF(5,3596,139,219);
		System.out.println(v);
	}

}
