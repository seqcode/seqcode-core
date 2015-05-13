package edu.psu.compbio.seqcode.projects.kunz.chromeSOM;

public class test {

	public static void main(String[] args) 
	{
		for(double sig = .1; sig<1; sig+=.1)
			for(int o = 0; o<=3; o++)
			{
				for(int i = 1000; i>=0; i-=100)
				{
					double r = Math.exp(-1*((o*o)/(2*sig*sig)));
					System.out.println(o+ " - " + " - "+ r);
				}
			}
	}
}
