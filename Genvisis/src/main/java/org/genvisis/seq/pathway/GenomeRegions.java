package seq.pathway;

import java.io.Serializable;

import common.Files;
import common.Logger;
import filesys.GeneTrack;

public class GenomeRegions implements Serializable {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private GeneTrack geneTrack;
	private Pathways pathways;
	private Logger log;

	public GenomeRegions(GeneTrack geneTrack, Pathways pathways, Logger log) {
		super();
		this.geneTrack = geneTrack;
		this.pathways = pathways;
		this.log = log;
	}

	public Logger getLog() {
		return log;
	}

	public GeneTrack getGeneTrack() {
		return geneTrack;
	}

	public Pathways getPathways() {
		return pathways;
	}

	public void serialize(String filename) {
		Files.writeSerial(this, filename);
	}

	public static GenomeRegions load(String filename) {
		return (GenomeRegions) Files.readSerial(filename, false, false);
	}

}
