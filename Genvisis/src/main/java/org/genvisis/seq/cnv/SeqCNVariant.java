package org.genvisis.seq.cnv;


import org.genvisis.filesys.CNVariant;

public class SeqCNVariant extends CNVariant {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private CNVExtraInfo[] cExtraInfos;

	public SeqCNVariant(String familyID, String individualID, byte chr, int start, int stop, int cn, double score, int numMarkers, int source, CNVExtraInfo[] cExtraInfos) {
		super(familyID, individualID, chr, start, stop, cn, score, numMarkers, source);
		this.cExtraInfos = cExtraInfos;
	}

	public CNVExtraInfo[] getcExtraInfos() {
		return cExtraInfos;
	}

	public void setcExtraInfos(CNVExtraInfo[] cExtraInfos) {
		this.cExtraInfos = cExtraInfos;
	}

	

}
