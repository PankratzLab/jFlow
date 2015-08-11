package cnv.annotation;

import htsjdk.samtools.Cigar;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Hashtable;

import seq.analysis.Blast;
import seq.analysis.Blast.BlastResults;
import seq.manage.CigarOps;
import stats.Histogram.DynamicHistogram;
import common.Array;
import common.ArraySpecialList.*;
import common.Files;
import cnv.annotation.BlastAnnotationTypes.BLAST_ANNOTATION_TYPES;
import cnv.annotation.BlastAnnotationTypes.BlastAnnotation;
import cnv.annotation.LocusAnnotation.Builder;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import cnv.filesys.Project.ARRAY;
import cnv.qc.MarkerBlast.MarkerFastaEntry;
import filesys.Segment;

/**
 * Handles
 */
public class BlastAnnotationWriter extends AnnotationFileWriter {

	private Project proj;
	private String[] blastResultFiles;
	private MarkerFastaEntry[] markerFastaEntries;
	private int minAlignmentLength;
	private int maxGaps;
	private int maxMismatches;
	private int maxAlignmentsReported;

	public BlastAnnotationWriter(Project proj, String outputFile, String[] blastResultFiles, int minAlignmentLength, int maxGaps, int maxMismatches, int maxAlignmentsReported) {
		super(proj, Array.concatAll(BlastAnnotationTypes.getBaseAnnotations(), new Annotation[] { MarkerBlastHistogramAnnotation.getDefaultBlastAnnotation() }, new Annotation[] { MarkerSeqAnnotation.getDefault() }), outputFile, false);
		this.proj = proj;
		this.blastResultFiles = blastResultFiles;
		this.minAlignmentLength = minAlignmentLength;
		this.maxGaps = maxGaps;
		this.maxMismatches = maxMismatches;
		this.maxAlignmentsReported = maxAlignmentsReported;
	}

	public void setMarkerFastaEntries(MarkerFastaEntry[] markerFastaEntries) {
		this.markerFastaEntries = markerFastaEntries;
	}

	public void summarizeResultFiles(boolean skipDefaultValue) {
		proj.getLog().reportTimeInfo(" Free memory " + proj.getLog().memoryPercentFree());
		proj.getLog().reportTimeInfo("Max number of alignments reported set to " + maxAlignmentsReported);
		LocusAnnotation[] annotations = summarizeResultFile(proj, blastResultFiles, markerFastaEntries, minAlignmentLength, maxGaps, maxMismatches, maxAlignmentsReported);
		for (int i = 0; i < annotations.length; i++) {
			if ((i + 1) % 100000 == 0) {
				proj.getLog().reportTimeInfo("Written " + i + " markers ");
				proj.getLog().reportTimeInfo(" Free memory " + proj.getLog().memoryPercentFree());

			}

			write(annotations[i], skipDefaultValue);
		}
	}

	/**
	 * @return the {@link LocusAnnotation } array for the markers in the project. This could also be passed to another {@link AnnotationFileWriter} instead of writing here
	 */
	private static LocusAnnotation[] summarizeResultFile(Project proj, String[] blastResultFile, MarkerFastaEntry[] markerFastaEntries, int minAlignmentLength, int maxGaps, int maxMismatches, int maxAlignmentsReported) {
		int seqLength = proj.getArrayType().getProbeLength();
		LocusAnnotation[] anDatas = initializeSummaries(proj, minAlignmentLength);
		MarkerBlastHistogramAnnotation[] mHistogramAnnotations = new MarkerBlastHistogramAnnotation[proj.getMarkerNames().length];
		ArrayBlastAnnotationList[][] intLists = new ArrayBlastAnnotationList[anDatas.length][BLAST_ANNOTATION_TYPES.values().length];
		for (int i = 0; i < intLists.length; i++) {
			for (int j = 0; j < intLists[i].length; j++) {
				intLists[i][j] = new ArrayBlastAnnotationList(maxAlignmentsReported);
			}
		}

		Hashtable<String, Integer> markerIndices = proj.getMarkerIndices();

		for (int fileIndex = 0; fileIndex < blastResultFile.length; fileIndex++) {
			proj.getLog().reportTimeInfo("Loading results from " + blastResultFile[fileIndex]);
			try {
				BufferedReader reader = Files.getAppropriateReader(blastResultFile[fileIndex]);
				int numEntries = 0;
				while (reader.ready()) {
					String[] line = reader.readLine().trim().split("\t");
					if ((line.length == Blast.BLAST_HEADER.length - 1 || line.length == Blast.BLAST_HEADER.length) && !line[0].startsWith(Blast.BLAST_HEADER[0])) {
						numEntries++;
						if (numEntries % 1000000 == 0) {
							proj.getLog().reportTimeInfo("Processed " + numEntries + " blast results");
							// reader.close();
							// break;
						}
						BlastResults blastResults = new BlastResults(line, proj.getLog());

						if (blastResults.getAlignmentLength() >= minAlignmentLength && blastResults.getGapOpens() <= maxGaps && blastResults.getMismatches() <= maxMismatches) {
							String marker = blastResults.getQueryID();
							if (proj.getArrayType() == ARRAY.AFFY_GW6 || proj.getArrayType() == ARRAY.AFFY_GW6_CN) {
								if (marker.endsWith("_A") || marker.endsWith("_B")) {
									marker = marker.substring(0, marker.length() - 2);
								} else {
									proj.getLog().reportTimeError("Query id did not end with _A or _B which is required for an AFFY array");
								}

							}
							int markerIndex = markerIndices.get(marker);
							
							Segment markerSeg = anDatas[markerIndex].getSeg().getBufferedSegment(1);
							if (mHistogramAnnotations[markerIndex] == null) {
								DynamicHistogram histogram = new DynamicHistogram(minAlignmentLength, proj.getArrayType().getProbeLength(), 0);
								mHistogramAnnotations[markerIndex] = new MarkerBlastHistogramAnnotation(MarkerBlastHistogramAnnotation.DEFAULT_NAME, MarkerBlastHistogramAnnotation.DEFAULT_DESCRIPTION, histogram);
							}
							mHistogramAnnotations[markerIndex].getDynamicHistogram().addDataPointToHistogram(blastResults.getAlignmentLength());

							for (int i = 0; i < BLAST_ANNOTATION_TYPES.values().length; i++) {
								if (BlastAnnotationTypes.shouldBeAnnotatedAs(proj, blastResults, BLAST_ANNOTATION_TYPES.values()[i], markerSeg, proj.getLog())) {
									BlastAnnotation blastAnnotation = new BlastAnnotation(Blast.convertBtopToCigar(blastResults, seqLength, proj.getLog()), blastResults.getSegment(),blastResults.determineStrand());
									intLists[markerIndex][i].add(blastAnnotation);
								}
							}
						}
					}
				}
				reader.close();
			} catch (FileNotFoundException fnfe) {
				proj.getLog().reportError("Error: file \"" + blastResultFile[fileIndex] + "\" not found in current directory");
				return null;
			} catch (IOException ioe) {
				proj.getLog().reportError("Error reading file \"" + blastResultFile[fileIndex] + "\"");
				return null;
			}
		}
		summarize(proj, maxAlignmentsReported, anDatas, mHistogramAnnotations, markerFastaEntries, intLists);
		return anDatas;
	}

	private static void summarize(Project proj, int maxAlignmentsReported, LocusAnnotation[] anDatas, MarkerBlastHistogramAnnotation[] mHistogramAnnotations, MarkerFastaEntry[] markerFastaEntries, ArrayBlastAnnotationList[][] intLists) {
		Hashtable<String, Integer> markerSeqIndices = new Hashtable<String, Integer>();
		if (markerFastaEntries != null) {
			for (int i = 0; i < markerFastaEntries.length; i++) {
				markerSeqIndices.put(markerFastaEntries[i].getName(), i);
			}
		}
		for (int i = 0; i < anDatas.length; i++) {
			if (i % 100000 == 0) {
				proj.getLog().reportTimeInfo("Summarized " + i + " markers ");
				proj.getLog().reportTimeInfo(" Free memory " + proj.getLog().memoryPercentFree());
			}
			for (int j = 0; j < BLAST_ANNOTATION_TYPES.values().length; j++) {
				BlastAnnotation[] tmp = intLists[i][j].toArray(new BlastAnnotation[intLists[i][j].size()]);
				if (tmp.length > 0) {
					Cigar[] cigarstmp = new Cigar[tmp.length];
					for (int k = 0; k < tmp.length; k++) {
						cigarstmp[k] = tmp[k].getCigar();
					}
					int[] order = CigarOps.sortByRefMatchLength(cigarstmp);

					BlastAnnotation[] toReport = new BlastAnnotation[Math.min(tmp.length, maxAlignmentsReported)];
					for (int k = 0; k < toReport.length; k++) {
						toReport[k] = tmp[order[k]];
					}

					anDatas[i].getAnnotations()[j].setData(Array.toStr(BlastAnnotation.toAnnotationString(toReport), BLAST_ANNOTATION_TYPES.values()[j].getSep()));
				} else {
					anDatas[i].getAnnotations()[j].setData(BLAST_ANNOTATION_TYPES.values()[j].getDefaultValue());
				}
			}
			if (mHistogramAnnotations[i] != null) {
				mHistogramAnnotations[i].setDataToHistogram(true);
				anDatas[i].addAnnotation(mHistogramAnnotations[i]);
			}
			if (markerSeqIndices.containsKey(anDatas[i].getLocusName())) {
				int seqIndex = markerSeqIndices.get(anDatas[i].getLocusName());
				MarkerFastaEntry markerFastaEntry = markerFastaEntries[seqIndex];
				MarkerSeqAnnotation markerSeqAnnotation = new MarkerSeqAnnotation();
				markerSeqAnnotation.setDesignData(markerFastaEntry.getSequence(), markerFastaEntry.getInterrogationPosition(),markerFastaEntry.getStrand());
				anDatas[i].addAnnotation(markerSeqAnnotation);
			}
		}
	}

	/**
	 * @param proj
	 * @return initialized blast summaries for all markers
	 */
	private static LocusAnnotation[] initializeSummaries(Project proj, int minAlignmentLength) {
		MarkerSet markerSet = proj.getMarkerSet();
		// TODO ablookup for annovar etc...
		byte[] chrs = markerSet.getChrs();
		int[] pos = markerSet.getPositions();
		String[] markerNames = proj.getMarkerNames();
		LocusAnnotation[] anDatas = new LocusAnnotation[markerNames.length];
		for (int i = 0; i < anDatas.length; i++) {
			Builder builder = new Builder();
			DynamicHistogram histogram = new DynamicHistogram(minAlignmentLength, proj.getArrayType().getProbeLength(), 0);
			MarkerBlastHistogramAnnotation markerBlastHistogramAnnotation = new MarkerBlastHistogramAnnotation(MarkerBlastHistogramAnnotation.DEFAULT_NAME, MarkerBlastHistogramAnnotation.DEFAULT_DESCRIPTION, histogram);
			builder.annotations(Array.concatAll(BlastAnnotationTypes.getAnnotationDatas(), new AnnotationData[] { markerBlastHistogramAnnotation }, new AnnotationData[] { MarkerSeqAnnotation.getDefault() }));
			Segment markerSeg = new Segment(chrs[i], pos[i], pos[i]);
			anDatas[i] = builder.build(markerNames[i], markerSeg);
		}
		return anDatas;
	}

}
