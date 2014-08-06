package edu.psu.compbio.seqcode.gse.projects.gps.utilities;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.location.Gene;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.Expander;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.NamedGeneratorFactory;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.NamedStrandedGeneratorFactory;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.NamedTypedGeneratorFactory;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.RefGeneGeneratorFactory;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.RepeatMaskedGenerator;

public class AnnotationLoader {

	protected Expander<Region, ? extends Region> annotExpander;
	protected String sourceName;
	protected int maxAnnotDist=50000;
	protected boolean overlapOnly=false;
	
	public AnnotationLoader(Genome gen, String name, String type, int maxAnnotDist, boolean overlapOnly){
		sourceName = name;
		this.maxAnnotDist=maxAnnotDist;
		this.overlapOnly=overlapOnly;
		
		if(type.equals("refGene")){
			annotExpander = (new RefGeneGeneratorFactory()).getExpander(gen);
		}else if(type.equals("namedRegion")){
			annotExpander = (new NamedGeneratorFactory(sourceName)).getExpander(gen);
		}else if (type.equals("namedStrandedRegions")){
			annotExpander = (new NamedStrandedGeneratorFactory(sourceName)).getExpander(gen);
		}else if (type.equals("namedTypedRegions")){
			annotExpander = (new NamedTypedGeneratorFactory(sourceName)).getExpander(gen);
		}else if (type.equals("repeatMasker")){
			annotExpander = new RepeatMaskedGenerator(gen);
		}else if (type.equals("file")){
			annotExpander = new TranscriptFileExpander(gen, name);
		}
	}

	public Collection<Region> getAnnotations(Region coords){
		ArrayList<Region> results = new ArrayList<Region>();
		Region query;
		if (overlapOnly) {
            query = coords;
        } else {
            query = coords.expand(maxAnnotDist, maxAnnotDist);
        }
        Iterator<? extends Region> iter = annotExpander.execute(query);
        while (iter.hasNext()) {
            results.add(iter.next());
        }
        return(results);
    }
	
	/** Return genes if the expander is over genes. 
	 * There is probably a better way of doing this that doesn't rely on instanceof 
	 */ 
	public Collection<Gene> getGenes(Region coords){
		ArrayList<Gene> results = new ArrayList<Gene>();
		Region query;
		if (overlapOnly) {
            query = coords;
        } else {
            query = coords.expand(maxAnnotDist, maxAnnotDist);
        }
        Iterator<? extends Region> iter = annotExpander.execute(query);
        while (iter.hasNext()) {
        	Region r = iter.next();
        	if(r instanceof Gene)
        		results.add((Gene)r);
        }
        return(results);
    }
}
