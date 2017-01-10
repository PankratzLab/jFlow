package org.genvisis.seq.cnv;


import java.util.ArrayList;

import org.genvisis.common.Array;
import org.genvisis.filesys.CNVariant;

public class SeqCNVariant extends CNVariant {
	/**
	 *
	 */
	private static final long serialVersionUID = 1L;
	private CNVExtraInfo[] cExtraInfos;

	public SeqCNVariant(String familyID, String individualID, byte chr, int start, int stop, int cn,
											double score, int numMarkers, int source, CNVExtraInfo[] cExtraInfos) {
		super(familyID, individualID, chr, start, stop, cn, score, numMarkers, source);
		this.cExtraInfos = cExtraInfos;
	}

	public SeqCNVariant(CNVariant cnv, CNVExtraInfo[] cExtraInfos) {
		super(cnv);
		this.cExtraInfos = cExtraInfos;
	}

	public CNVExtraInfo[] getcExtraInfos() {
		return cExtraInfos;
	}

	public void setcExtraInfos(CNVExtraInfo[] cExtraInfos) {
		this.cExtraInfos = cExtraInfos;
	}


	@Override
	public String toAnalysisString() {
		ArrayList<String> extraInfo = new ArrayList<String>();
		if (cExtraInfos != null) {
			for (CNVExtraInfo cnvExtraInfo : cExtraInfos) {
				extraInfo.add(cnvExtraInfo.getdExtra());
			}
		}
		return toPlinkFormat() + "\t" + Array.toStr(extraInfo);
	}

	@Override
	public String[] getHeader() {
		ArrayList<String> extraHeaders = new ArrayList<String>();
		if (cExtraInfos != null) {
			for (CNVExtraInfo cnvExtraInfo : cExtraInfos) {
				extraHeaders.add(cnvExtraInfo.getsExtra());
			}
		}
		return Array.concatAll(PLINK_CNV_HEADER, Array.toStringArray(extraHeaders));
	}



}
