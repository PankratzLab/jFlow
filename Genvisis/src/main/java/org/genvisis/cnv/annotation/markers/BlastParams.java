package org.genvisis.cnv.annotation.markers;

import java.util.ArrayList;

import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Logger;

import htsjdk.variant.vcf.VCFHeaderLine;

public class BlastParams implements AnalysisParams {

	/**
	 *
	 */
	private static final String KEY = "GENVISIS_BLAST_PARAMETERS";
	private static final String DATA_DELIMITER = ",";
	private static final String[] parseKeys =
																					new String[] {"fileSeq=", "ref=",
																												"maxAlignmentsReported=", "reportWordSize=",
																												"blastWordSize=", "evalueCutoff=", "date=",
																												"markerFingerPrint=", "notes="};

	private String fileSeq;
	private String fastaDb;
	private int maxAlignmentsReported;
	private int reportWordSize;
	private int blastWordSize;
	private String dateStamp;
	private double evalueCutoff;
	private long markerFingerPrint;
	private String notes;
	private boolean sawValidHeaderLine;
	private final Logger log;

	public BlastParams(Logger log) {
		this.log = log;
		sawValidHeaderLine = false;
	}

	public BlastParams(	String fileSeq, String fastaDb, int maxAlignmentsReported, int reportWordSize,
											int blastWordSize, String dateStamp, double evalueCutoff,
											long markerFingerPrint, String notes, Logger log) {
		super();
		this.fileSeq = fileSeq;
		this.fastaDb = fastaDb;
		this.maxAlignmentsReported = maxAlignmentsReported;
		this.reportWordSize = reportWordSize;
		this.blastWordSize = blastWordSize;
		this.dateStamp = dateStamp;
		this.evalueCutoff = evalueCutoff;
		this.markerFingerPrint = markerFingerPrint;
		this.notes = notes;
		sawValidHeaderLine = false;
		this.log = log;
	}

	public String getFileSeq() {
		return fileSeq;
	}

	public String getFastaDb() {
		return fastaDb;
	}

	public String getNotes() {
		return notes;
	}

	public void setNotes(String notes) {
		this.notes = notes;
	}

	public int getMaxAlignmentsReported() {
		return maxAlignmentsReported;
	}

	public int getReportWordSize() {
		return reportWordSize;
	}

	public int getBlastWordSize() {
		return blastWordSize;
	}

	public String getDateStamp() {
		return dateStamp;
	}

	public double getEvalueCutoff() {
		return evalueCutoff;
	}

	public long getMarkerFingerPrint() {
		return markerFingerPrint;
	}

	public Logger getLog() {
		return log;
	}

	@Override
	public VCFHeaderLine developHeaderLine() {
		ArrayList<String> valueString = new ArrayList<String>();
		valueString.add(parseKeys[0] + fileSeq);
		valueString.add(parseKeys[1] + fastaDb);
		valueString.add(parseKeys[2] + maxAlignmentsReported);
		valueString.add(parseKeys[3] + reportWordSize);
		valueString.add(parseKeys[4] + blastWordSize);
		valueString.add(parseKeys[5] + evalueCutoff);
		valueString.add(parseKeys[6] + dateStamp);
		valueString.add(parseKeys[7] + markerFingerPrint);
		valueString.add(parseKeys[8] + notes);

		String value = ArrayUtils.toStr(ArrayUtils.toStringArray(valueString), DATA_DELIMITER);
		return new VCFHeaderLine(KEY, value);
	}

	public boolean isSawValidHeaderLine() {
		return sawValidHeaderLine;
	}

	@Override
	public String getKey() {
		return KEY;
	}

	@Override
	public void parseHeaderLine(VCFHeaderLine vcfHeaderLine) {
		if (vcfHeaderLine.getKey().equals(KEY)) {
			sawValidHeaderLine = true;
			String[] tmp = vcfHeaderLine.getValue().split(DATA_DELIMITER);
			for (String element : tmp) {
				if (element.startsWith(parseKeys[0])) {
					fileSeq = element.split("=")[1];
				} else if (element.startsWith(parseKeys[1])) {
					fastaDb = element.split("=")[1];
				} else if (element.startsWith(parseKeys[2])) {
					try {
						maxAlignmentsReported = Integer.parseInt(element.split("=")[1]);
					} catch (NumberFormatException nfe) {
						log.reportError("Could not parse max alignments reported " + element);
					}
				} else if (element.startsWith(parseKeys[3])) {
					try {
						reportWordSize = Integer.parseInt(element.split("=")[1]);
					} catch (NumberFormatException nfe) {
						log.reportError("Could not parse reportWordSize " + element);
					}
				} else if (element.startsWith(parseKeys[4])) {
					try {
						blastWordSize = Integer.parseInt(element.split("=")[1]);
					} catch (NumberFormatException nfe) {
						log.reportError("Could not parse blastWordSize " + element);
					}
				} else if (element.startsWith(parseKeys[5])) {
					try {
						evalueCutoff = Double.parseDouble(element.split("=")[1]);
					} catch (NumberFormatException nfe) {
						log.reportError("Could not parse evalueCutoff " + element);
					}
				} else if (element.startsWith(parseKeys[6])) {
					dateStamp = element.split("=")[1];
				} else if (element.startsWith(parseKeys[7])) {
					try {
						markerFingerPrint = Long.parseLong(element.split("=")[1]);
					} catch (NumberFormatException nfe) {
						log.reportError("Could not parse marker fingerprint " + element);
					}
				} else if (element.startsWith(parseKeys[8])) {
					String[] tmpS = element.split("=");
					if (tmpS.length > 1) {
						notes = element.split("=")[1];
					}
				}
			}
		}
	}
}
