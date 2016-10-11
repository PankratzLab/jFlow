package org.genvisis.seq.pathway;

import java.io.Serializable;

import org.genvisis.common.Logger;
import org.genvisis.common.SerializedFiles;
import org.genvisis.filesys.GeneTrack;

public class GenomeRegions implements Serializable {
	/**
	 *
	 */
	private static final long serialVersionUID = 1L;
	private final GeneTrack geneTrack;
	private final Pathways pathways;
	private final Logger log;

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
		SerializedFiles.writeSerial(this, filename);
	}

	public static GenomeRegions load(String filename) {
		return (GenomeRegions) SerializedFiles.readSerial(filename, false, false);
	}

}
