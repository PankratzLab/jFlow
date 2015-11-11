package seq.pathway;

import java.util.HashSet;
import java.util.Set;

import cnv.var.LocusSet;
import common.Logger;
import filesys.GeneData;

public class Pathway extends LocusSet<GeneData> {
	private static final long serialVersionUID = 1L;
	private String pathwayName;
	private Set<String> geneNames;
	private boolean complete;

	public Pathway(String pathwayName, GeneData[] loci, boolean sort, boolean complete, Logger log) {

		super(loci, sort, log);
		this.pathwayName = pathwayName;
		this.geneNames = new HashSet<String>();
		this.complete = complete;
		for (int i = 0; i < loci.length; i++) {
			geneNames.add(loci[i].getGeneName());
		}
	}

	public Pathway(GeneData[] loci, boolean sort, Logger log, String pathwayName, Set<String> geneNames, boolean complete) {
		super(loci, sort, log);
		this.pathwayName = pathwayName;
		this.geneNames = geneNames;
		this.complete = complete;
	}

	

	public boolean isComplete() {
		return complete;
	}

	public String getPathwayName() {
		return pathwayName;
	}

	public boolean containsGene(String gene) {
		return geneNames.contains(gene);
	}

}
