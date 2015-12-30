package cnv.annotation;

import java.util.ArrayList;

import common.Array;
import common.Logger;

import htsjdk.variant.vcf.VCFHeaderLine;

public class BlastParams implements AnalysisParams {

	/**
	 * 
	 */
	private static final String KEY = "GENVISIS_BLAST_PARAMETERS";
	private static final String DATA_DELIMITER = ",";
	private static final String[] parseKeys = new String[] { "fileSeq=", "ref=", "maxAlignmentsReported=", "reportWordSize=", "blastWordSize=", "evalueCutoff=", "date=", "markerFingerPrint=","notes=" };

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
	private Logger log;

	public BlastParams(Logger log) {
		this.log = log;
		this.sawValidHeaderLine = false;
	}

	public BlastParams(String fileSeq, String fastaDb, int maxAlignmentsReported, int reportWordSize, int blastWordSize, String dateStamp, double evalueCutoff, long markerFingerPrint,String notes, Logger log) {
		super();
		this.fileSeq = fileSeq;
		this.fastaDb = fastaDb;
		this.maxAlignmentsReported = maxAlignmentsReported;
		this.reportWordSize = reportWordSize;
		this.blastWordSize = blastWordSize;
		this.dateStamp = dateStamp;
		this.evalueCutoff = evalueCutoff;
		this.markerFingerPrint =markerFingerPrint;
		this.notes= notes;
		this.sawValidHeaderLine = false;
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

		String value = Array.toStr(Array.toStringArray(valueString), DATA_DELIMITER);
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
			for (int i = 0; i < tmp.length; i++) {
				if (tmp[i].startsWith(parseKeys[0])) {
					fileSeq = tmp[i].split("=")[1];
				} else if (tmp[i].startsWith(parseKeys[1])) {
					fastaDb = tmp[i].split("=")[1];
				} else if (tmp[i].startsWith(parseKeys[2])) {
					try {
						maxAlignmentsReported = Integer.parseInt(tmp[i].split("=")[1]);
					} catch (NumberFormatException nfe) {
						log.reportTimeError("Could not parse max alignments reported " + tmp[i]);
					}
				} else if (tmp[i].startsWith(parseKeys[3])) {
					try {
						reportWordSize = Integer.parseInt(tmp[i].split("=")[1]);
					} catch (NumberFormatException nfe) {
						log.reportTimeError("Could not parse reportWordSize " + tmp[i]);
					}
				} else if (tmp[i].startsWith(parseKeys[4])) {
					try {
						blastWordSize = Integer.parseInt(tmp[i].split("=")[1]);
					} catch (NumberFormatException nfe) {
						log.reportTimeError("Could not parse blastWordSize " + tmp[i]);
					}
				} else if (tmp[i].startsWith(parseKeys[5])) {
					try {
						evalueCutoff = Double.parseDouble(tmp[i].split("=")[1]);
					} catch (NumberFormatException nfe) {
						log.reportTimeError("Could not parse evalueCutoff " + tmp[i]);
					}
				} else if (tmp[i].startsWith(parseKeys[6])) {
					dateStamp = tmp[i].split("=")[1];
				} else if (tmp[i].startsWith(parseKeys[7])) {
					try {
						markerFingerPrint = Long.parseLong(tmp[i].split("=")[1]);
					} catch (NumberFormatException nfe) {
						log.reportTimeError("Could not parse marker fingerprint " + tmp[i]);
					}
				} else if (tmp[i].startsWith(parseKeys[8])) {
					notes = tmp[i].split("=")[1];
				}
			}
		}
	}
}
