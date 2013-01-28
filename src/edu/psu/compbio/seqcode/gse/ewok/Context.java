package edu.psu.compbio.seqcode.gse.ewok;

import edu.psu.compbio.seqcode.gse.datasets.species.Genome;

/** Context is a set global variables that we don't want to pass around
 * to lots of method calls.  The static variable defaultContext should be
 * set so that methods have something to grab if no other context is provided 
 */
public class Context {
    private Genome genome;
    private static Context defaultcontext;

    public Context (Genome g) {
        genome = g;
    }
    public static void setDefaultContext(Context c) {defaultcontext = c;}
    public static Context defaultContext() {return defaultcontext;}
    public Genome getGenome() {return genome;}
}
