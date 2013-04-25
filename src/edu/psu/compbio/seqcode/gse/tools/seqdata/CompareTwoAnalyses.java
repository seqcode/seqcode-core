package edu.psu.compbio.seqcode.gse.tools.seqdata;

import java.sql.SQLException;
import java.util.*;

import edu.psu.compbio.seqcode.gse.datasets.general.*;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.*;
import edu.psu.compbio.seqcode.gse.datasets.species.*;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.*;


/**
 * Base class for comparing two chipseq analyses.
 *
 * main() should instantiate the object, call parseArgs(),
 * and then call printOutputEvents().
 *
 * Subclasses must override getOutputEvents() to provide the
 * list of output events.
 *
 * parseArgs() here looks for 
 * --one 'analysisname;version'
 * --two 'analysisname;version'
<<<<<<< HEAD
 * [--maxdist 100] how far apart can events be to be considered the same by containsMatch
 * [--topevents 1000] only include the top n events in analysis
 * [--sortbypval] sort events by pvalue rather than fold enrichment
=======
 * [--maxdist 50]  maximum distance to consider events in one and two the same
 * 
>>>>>>> 9bab5e357267fff52fd7c2d9cd8512a94e3d6a2b
 */

public abstract class CompareTwoAnalyses {

    private SeqAnalysis one, two;
    private int maxDistance = 20;
    private Genome genome;
    private List<Region> analysisRegions;
    private int topEvents = -1, firstCheck = 100;
    private boolean sortEventsByPval = false;

    public CompareTwoAnalyses() {}

    public void parseArgs(String args[]) throws NotFoundException, SQLException {
        one = Args.parseSeqAnalysis(args,"one");
        two = Args.parseSeqAnalysis(args,"two");
        maxDistance = Args.parseInteger(args,"maxdist",maxDistance);
        topEvents = Args.parseInteger(args,"topevents",-1);
        sortEventsByPval = Args.parseFlags(args).contains("sortbypval");
        firstCheck = Args.parseInteger(args,"firstcheck",firstCheck);
        genome = Args.parseGenome(args).cdr();
        analysisRegions = Args.parseRegionsOrDefault(args);
    }
    public SeqAnalysis getAnalysisOne() {return one;}
    public SeqAnalysis getAnalysisTwo() {return two;}
    public int getMaxDistance() { return maxDistance;}
    public Genome getGenome() {return genome;}
    public List<SeqAnalysisResult> getResultsOne() throws SQLException {
        List<SeqAnalysisResult> output = new ArrayList<SeqAnalysisResult>();
        for (Region r : getAnalysisRegions()) {
            output.addAll(getResultsOne(r));
        }
        return filterEvents(output);
    }
    public List<SeqAnalysisResult> filterEvents(List<SeqAnalysisResult> input) {
        if (topEvents < 1) {
            return input;
        } else {
            if (sortEventsByPval) {
                Collections.sort(input,new SeqAnalysisResultPvalueComparator());
            } else {
                Collections.sort(input,new SeqAnalysisResultEnrichmentComparator());
            }
        }
        return input.subList(0,Math.min(topEvents, input.size()));

    }
    public List<SeqAnalysisResult> getResultsTwo() throws SQLException {
        List<SeqAnalysisResult> output = new ArrayList<SeqAnalysisResult>();
        for (Region r : getAnalysisRegions()) {
            output.addAll(getResultsTwo(r));
        }
        return filterEvents(output);
    }
    public List<SeqAnalysisResult> getResultsOne(Region region) throws SQLException {
        return one.getResults(region);
    }
    public List<SeqAnalysisResult> getResultsTwo(Region region) throws SQLException {
        return two.getResults(region);
    }
    public List<Region> getAnalysisRegions() {
        return analysisRegions;
    }
    /** returns true if the sorted list contains r, according
        to whatever overlap criteria we're using
     */
    public boolean containsMatch(List<SeqAnalysisResult> list,
                                 SeqAnalysisResult r) {
        int i = Collections.binarySearch(list,r);
        if (i >= 0) {
            return true;
        }
        int inspt = -1 - i;
        try {
            if (inspt < list.size() && list.get(inspt).distance(r) <= maxDistance) {
                return true;
            }
        } catch (IllegalArgumentException e) {
            // ignore it. happens if the events aren't on the same chromosome
        }
        try {
            if (inspt > 0 && list.get(inspt-1).distance(r) <= maxDistance) {
                return true;
            }
        } catch (IllegalArgumentException e) {
            // ignore it.
        }

        return false;
    }
    public abstract List<SeqAnalysisResult> getOutputEvents() throws SQLException ;
    public void printOutputEvents() throws SQLException {
        for (SeqAnalysisResult r : getOutputEvents()) {
            System.out.println(String.format("%s:%d-%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2e\t%.1f",
                                             r.getChrom(), r.getStart(), r.getEnd(),
                                             r.getPosition(),
                                             r.getFG(), r.getBG(), r.getStrength(), r.getShape(),
                                             r.getPValue(), r.getFoldEnrichment()));
        }
    }
}