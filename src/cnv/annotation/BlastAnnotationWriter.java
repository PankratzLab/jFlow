package cnv.annotation;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarOperator;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Hashtable;

import seq.analysis.Blast;
import seq.analysis.Blast.BlastResults;
import seq.manage.CigarOps;
import common.Array;
import common.ArraySpecialList.*;
import common.Files;
import cnv.annotation.BlastAnnotationTypes.BLAST_ANNOTATION_TYPES;
import cnv.annotation.BlastAnnotationTypes.BlastAnnotation;
import cnv.annotation.LocusAnnotation.Builder;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import filesys.Segment;

/**
 * Handles
 */
public class BlastAnnotationWriter extends AnnotationFileWriter {

	private Project proj;
	private String[] blastResultFiles;
	private int minAlignmentLength;
	private int maxGaps;
	private int maxMismatches;
	private int maxAlignmentsReported;

	public BlastAnnotationWriter(Project proj, String outputFile, String[] blastResultFiles, int minAlignmentLength, int maxGaps, int maxMismatches, int maxAlignmentsReported) {
		super(proj, BlastAnnotationTypes.getBaseAnnotations(), outputFile, true);
		this.proj = proj;
		this.blastResultFiles = blastResultFiles;
		this.minAlignmentLength = minAlignmentLength;
		this.maxGaps = maxGaps;
		this.maxMismatches = maxMismatches;
		this.maxAlignmentsReported = maxAlignmentsReported;
	}

	public void summarizeResultFiles() {
		LocusAnnotation[] annotations = summarizeResultFile(proj, blastResultFiles, minAlignmentLength, maxGaps, maxMismatches, maxAlignmentsReported);
		for (int i = 0; i < annotations.length; i++) {
			if ((i + 1) % 100000 == 0) {
				proj.getLog().reportTimeInfo("Written " + i + " markers ");
				proj.getLog().reportTimeInfo(" Free memory " + proj.getLog().memoryPercentFree());

			}
			write(annotations[i]);
		}
	}

	/**
	 * @return the {@link LocusAnnotation } array for the markers in the project. This could also be passed to another {@link AnnotationFileWriter} instead of writing here
	 */
	private static LocusAnnotation[] summarizeResultFile(Project proj, String[] blastResultFile, int minAlignmentLength, int maxGaps, int maxMismatches, int maxAlignmentsReported) {
		int seqLength = proj.getArrayType().getProbeLength();
		LocusAnnotation[] anDatas = initializeSummaries(proj);
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
							int markerIndex = markerIndices.get(marker);
							Segment markerSeg = anDatas[markerIndex].getSeg().getBufferedSegment(1);

							for (int i = 0; i < BLAST_ANNOTATION_TYPES.values().length; i++) {
								if (BlastAnnotationTypes.shouldBeAnnotatedAs(proj, blastResults, BLAST_ANNOTATION_TYPES.values()[i], markerSeg, proj.getLog())) {
									BlastAnnotation blastAnnotation = new BlastAnnotation(Blast.convertBtopToCigar(blastResults, seqLength, proj.getLog()), blastResults.getSegment());
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
		for (int i = 0; i < anDatas.length; i++) {
			if (i % 100000 == 0) {
				proj.getLog().reportTimeInfo("Summarized " + i + " markers ");
				proj.getLog().reportTimeInfo(" Free memory " + proj.getLog().memoryPercentFree());
			}
			for (int j = 0; j < anDatas[i].getAnnotations().length; j++) {
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
		}
		return anDatas;
	}

	/**
	 * @param proj
	 * @return initialized blast summaries for all markers
	 */
	private static LocusAnnotation[] initializeSummaries(Project proj) {
		MarkerSet markerSet = proj.getMarkerSet();
		byte[] chrs = markerSet.getChrs();
		int[] pos = markerSet.getPositions();
		String[] markerNames = proj.getMarkerNames();
		LocusAnnotation[] anDatas = new LocusAnnotation[markerNames.length];
		for (int i = 0; i < anDatas.length; i++) {
			Builder builder = new Builder();
			builder.annotations(BlastAnnotationTypes.getAnnotationDatas());
			Segment markerSeg = new Segment(chrs[i], pos[i], pos[i]);
			anDatas[i] = builder.build(markerNames[i], markerSeg);
		}
		return anDatas;
	}

}
