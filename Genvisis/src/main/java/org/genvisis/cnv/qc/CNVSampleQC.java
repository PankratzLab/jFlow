package org.genvisis.cnv.qc;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

public class CNVSampleQC {
	public static final String[] QC_HEADS = { "LRR_mean", "LRR_median", "LRR_SD", "BAF_mean", "BAF_median", "BAF_SD", "BAF_DRIFT", "WF", "GCWF" };
	public static final String[] ID_HEADS = { "Sample", "FID", "IID" };
	private String sample;
	private String FID;
	private String IID;
	private double lrrMean;
	private double lrrMedian;
	private double lrrSDev;
	private double bafMean;
	private double bafMedian;
	private double bafSD;
	private double bafDrift;
	private double WF;
	private double GCWF;

	public CNVSampleQC(String sample, String fID, String iID, double lrrMean, double lrrMedian, double lrrSDev, double bafMean, double bafMedian, double bafSD, double bafDrift, double wF, double gCWF) {
		super();
		this.sample = sample;
		FID = fID;
		IID = iID;
		this.lrrMean = lrrMean;
		this.lrrMedian = lrrMedian;
		this.lrrSDev = lrrSDev;
		this.bafMean = bafMean;
		this.bafMedian = bafMedian;
		this.bafSD = bafSD;
		this.bafDrift = bafDrift;
		WF = wF;
		GCWF = gCWF;
	}

	public static CNVSampleQC[] getCNVSampleQCFromFile(Project proj, String QCFile) {
		ArrayList<CNVSampleQC> cnvSampleQCs = new ArrayList<CNVSampleQC>();
		Logger log = proj.getLog();
		
		try {
			BufferedReader reader = Files.getAppropriateReader(proj.PROJECT_DIRECTORY.getValue()+QCFile);
			String[] line;
			do {
				line = reader.readLine().trim().split("\t", -1);
			} while (reader.ready() && (ext.indexFactors(QC_HEADS, line, false, true)[0] == -1) && (ext.indexFactors(ID_HEADS, line, false, true)[0] == -1));
			if (!reader.ready()) {
				log.reportError("Error - reached the end of the file without finding a line with the following tokens: " + Array.toStr(QC_HEADS));
				log.reportError("      - perhaps the delimiter is set incorrectly? Determing most stable delimiter...");
				reader.close();
			}
			int[] QCIndices = ext.indexFactors(QC_HEADS, line, false, true);
			int[] IDIndices = ext.indexFactors(ID_HEADS, line, false, true);
			while (reader.ready()) {
				line = reader.readLine().split("\t");
				cnvSampleQCs.add(new CNVSampleQC(line[IDIndices[0]], line[IDIndices[1]], line[IDIndices[2]], Double.parseDouble(line[QCIndices[0]]), Double.parseDouble(line[QCIndices[1]]), Double.parseDouble(line[QCIndices[2]]), Double.parseDouble(line[QCIndices[3]]), Double.parseDouble(line[QCIndices[4]]), Double.parseDouble(line[QCIndices[5]]), Double.parseDouble(line[QCIndices[6]]), Double.parseDouble(line[QCIndices[7]]), Double.parseDouble(line[QCIndices[8]])));
			}
			

		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + proj.PROJECT_DIRECTORY.getValue() + QCFile + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + proj.PROJECT_DIRECTORY.getValue() + QCFile + "\"");
			System.exit(1);
		}
		return cnvSampleQCs.toArray(new CNVSampleQC[cnvSampleQCs.size()]);
	}

	public static Hashtable<String, CNVSampleQC> getSampleQCs(Project proj, String QCFile) {
		Hashtable<String, CNVSampleQC> cnvSampleQCHash = new Hashtable<String, CNVSampleQC>();
		if (QCFile == null) {
			return null;
		}
		CNVSampleQC[] cnvSampleQCs = getCNVSampleQCFromFile(proj, QCFile);
		for (int i = 0; i < cnvSampleQCs.length; i++) {
			cnvSampleQCHash.put(cnvSampleQCs[i].getFID() + "\t" + cnvSampleQCs[i].getIID(), cnvSampleQCs[i]);
		}
		return cnvSampleQCHash;

	}

	public String getDisplayString() {
		return this.sample + "\t" + this.FID + "\t" + this.IID + "\t" + this.lrrMean + "\t" + this.lrrMedian + "\t" + this.lrrSDev + "\t" + this.bafMean + "\t" + this.bafMedian + "\t" + this.bafSD + "\t" + this.bafDrift + "\t" + this.WF + "\t" + this.GCWF;
	}

	public String getSample() {
		return sample;
	}
	public String getFID() {
		return FID;
	}

	public String getIID() {
		return IID;
	}

	public double getLrrMean() {
		return lrrMean;
	}

	public double getLrrMedian() {
		return lrrMedian;
	}

	public double getLrrSDev() {
		return lrrSDev;
	}

	public double getBafMean() {
		return bafMean;
	}

	public double getBafMedian() {
		return bafMedian;
	}

	public double getBafSD() {
		return bafSD;
	}

	public double getBafDrift() {
		return bafDrift;
	}


	public double getWF() {
		return WF;
	}

	public double getGCWF() {
		return GCWF;
	}

	public static String[] getQcHeads() {
		return QC_HEADS;
	}

	public static String[] getIdHeads() {
		return ID_HEADS;
	}

}
