package cnv.annotation;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarOperator;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Arrays;
import java.util.Hashtable;

import seq.analysis.Blast;
import seq.analysis.Blast.BlastResults;
import seq.manage.CigarOps;
import common.Array;
import common.ArraySpecialList.*;
import common.Files;
import cnv.annotation.BlastAnnotationTypes.BLAST_ANNOTATION_TYPES;
import cnv.annotation.LocusAnnotation.Builder;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import filesys.Segment;

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

	/**
	 * @param proj
	 * @return initilized blast summaries for all markers
	 */
	private static LocusAnnotation[] initializeSummaries(Project proj) {
		MarkerSet markerSet = proj.getMarkerSet();
		byte[] chrs = markerSet.getChrs();
		int[] pos = markerSet.getPositions();
		String[] markerNames = proj.getMarkerNames();

		LocusAnnotation[] anDatas = new LocusAnnotation[markerNames.length];
		//AnnotationData[] annotationDatas = BlastAnnotationTypes.getAnnotationDatas();
		for (int i = 0; i < anDatas.length; i++) {
			Builder builder = new Builder();
			builder.annotations(BlastAnnotationTypes.getAnnotationDatas());
			Segment markerSeg = new Segment(chrs[i], pos[i], pos[i]);
			anDatas[i] = builder.build(markerNames[i], markerSeg);
		}
		return anDatas;
	}

	public void summarizeResultFiles() {
		LocusAnnotation[] annotations = summarizeResultFile(proj, blastResultFiles[0], minAlignmentLength, maxGaps, maxMismatches, maxAlignmentsReported);
		for (int i = 0; i < annotations.length; i++) {
			if (i % 10000 == 0) {
				proj.getLog().reportTimeInfo("Written " + i + " markers ");
				proj.getLog().reportTimeInfo(" Free memory " + proj.getLog().memoryPercentFree());
				// reader.close();
				// break;
			}
			//System.out.println("HI"+annotations[i].getAnnotations()[0].getData());
			write(annotations[i]);
		}
	}

	private static LocusAnnotation[] summarizeResultFile(Project proj, String blastResultFile, int minAlignmentLength, int maxGaps, int maxMismatches, int maxAlignmentsReported) {
		int seqLength = proj.getArrayType().getProbeLength();
		LocusAnnotation[] anDatas = initializeSummaries(proj);
		ArrayCigarList[][] intLists = new ArrayCigarList[anDatas.length][BLAST_ANNOTATION_TYPES.values().length];
		for (int i = 0; i < intLists.length; i++) {
			for (int j = 0; j < intLists[i].length; j++) {
				intLists[i][j] = new ArrayCigarList(maxAlignmentsReported);
			}
		}
		Hashtable<String, Integer> markerIndices = proj.getMarkerIndices();
		try {
			BufferedReader reader = Files.getAppropriateReader(blastResultFile);
			int numEntries = 0;
			while (reader.ready()) {
				String[] line = reader.readLine().trim().split("\t");
				if ((line.length == Blast.BLAST_HEADER.length - 1 || line.length == Blast.BLAST_HEADER.length) && !line[0].startsWith(Blast.BLAST_HEADER[0])) {
					numEntries++;
					if (numEntries % 100000 == 0) {
						proj.getLog().reportTimeInfo("Processed " + numEntries + " blast results");
						// reader.close();
						// break;
					}
					BlastResults blastResults = new BlastResults(line, proj.getLog());
					String marker = blastResults.getQueryID();
					int markerIndex = markerIndices.get(marker);
					Segment markerSeg = anDatas[markerIndex].getSeg().getBufferedSegment(1);

					for (int i = 0; i < BLAST_ANNOTATION_TYPES.values().length; i++) {
						if (BlastAnnotationTypes.shouldBeAnnotatedAs(proj, blastResults, BLAST_ANNOTATION_TYPES.values()[i], markerSeg, proj.getLog())) {
							intLists[markerIndex][i].add(Blast.convertBtopToCigar(blastResults, seqLength, proj.getLog()));
						}
					}
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			proj.getLog().reportError("Error: file \"" + blastResultFile + "\" not found in current directory");
			return null;
		} catch (IOException ioe) {
			proj.getLog().reportError("Error reading file \"" + blastResultFile + "\"");
			return null;
		}
		for (int i = 0; i < anDatas.length; i++) {
			if (i % 10000 == 0) {
				proj.getLog().reportTimeInfo("Summarized " + i + " markers ");
				proj.getLog().reportTimeInfo(" Free memory " + proj.getLog().memoryPercentFree());
				// reader.close();
				// break;
			}
			for (int j = 0; j < anDatas[i].getAnnotations().length; j++) {
				Cigar[] tmp = intLists[i][j].toArray(new Cigar[intLists[i][j].size()]);
				if (tmp.length > 0) {
					tmp = CigarOps.sortByRefMatchLength(tmp);
					tmp = Array.subArray(tmp, Array.intArray(Math.min(tmp.length, maxAlignmentsReported)));
					anDatas[i].getAnnotations()[j].setData(Array.toStr(CigarOps.toStringArray(tmp), BLAST_ANNOTATION_TYPES.values()[j].getSep()));
					
					System.out.println(anDatas[i].getLocusName()+"\n"+Array.toStr(CigarOps.toStringArray(tmp), BLAST_ANNOTATION_TYPES.values()[j].getSep()));
					
					//System.out.println(anDatas[i].getAnnotations()[j].getData()+"\n"+Array.toStr(CigarOps.toStringArray(tmp), BLAST_ANNOTATION_TYPES.values()[j].getSep()));

				} 
				//else {
//					anDatas[i].getAnnotations()[j].setData(CigarOps.getConstantCigar(seqLength, CigarOperator.X).toString());
//				}
			}
		}
		return anDatas;
	}

	public static void test() {
		Project proj = new Project("/home/pankrat2/lanej/projects/aric_exome.properties", false);
		String outfile = proj.PROJECT_DIRECTORY.getValue() + "Blasts/blast.anno.gz";
		String[] blastResultFiles = new String[] { proj.PROJECT_DIRECTORY.getValue() + "Blasts/GPL18544_humanexome-12v1_a.csv.blasted.ws.30.rep.0.tmp6" };
		int minAlignmentLength = proj.getArrayType().getProbeLength() - 10;
		int maxGaps = 0;
		int maxMismatches = 0;
		BlastAnnotationWriter blastAnnotation = new BlastAnnotationWriter(proj, outfile, blastResultFiles, minAlignmentLength, maxGaps, maxMismatches, 15);
		blastAnnotation.summarizeResultFiles();
		blastAnnotation.close();
	}

	public static void main(String[] args) {
		test();
	}

}
