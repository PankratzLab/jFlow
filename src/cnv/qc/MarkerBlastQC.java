package cnv.qc;

import java.util.ArrayList;

import cnv.annotation.AnnotationParser;
import cnv.annotation.MarkerAnnotationLoader;
import cnv.annotation.MarkerBlastAnnotation;
import cnv.annotation.MarkerGCAnnotation;
import cnv.annotation.AnnotationFileLoader.QUERY_ORDER;
import cnv.annotation.BlastAnnotationTypes.BLAST_ANNOTATION_TYPES;
import cnv.annotation.BlastAnnotationTypes.BlastAnnotation;
import cnv.annotation.BlastAnnotationTypes.PROBE_TAG;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import common.Array;
import common.Files;
import common.Logger;
import common.ext;

public class MarkerBlastQC {
	
	public static final double DEFAULT_CROSS_HYBE_THRESHOLD = 0.8;
	
	public static String defaultOneHitWondersFilename(String blastVCF) {
		return ext.rootRootOf(blastVCF) + "_OneHitWonders.txt";
	}

	public static void getOneHitWonders(Project proj, String blastVCF, String outFile, double crossHybePercent, Logger log) {
		if (blastVCF == null) {
			blastVCF = proj.BLAST_ANNOTATION_FILENAME.getValue();
		}
		if (outFile == null) {
			outFile = defaultOneHitWondersFilename(blastVCF);
		}
		MarkerSet markerSet = proj.getMarkerSet();
		String[] markerNames = markerSet.getMarkerNames();
		MarkerAnnotationLoader markerAnnotationLoader = new MarkerAnnotationLoader(proj, null, proj.BLAST_ANNOTATION_FILENAME.getValue(), proj.getMarkerSet(), true);
		markerAnnotationLoader.setReportEvery(500000);
		MarkerGCAnnotation[] gcAnnotations = MarkerGCAnnotation.initForMarkers(proj, markerNames, markerAnnotationLoader.getMarkerSet(), markerAnnotationLoader.getIndices());
		MarkerBlastAnnotation[] blastResults = MarkerBlastAnnotation.initForMarkers(markerNames);
		ArrayList<AnnotationParser[]> parsers = new ArrayList<AnnotationParser[]>();
		parsers.add(gcAnnotations);
		parsers.add(blastResults);
		ArrayList<String> oneHitters = new ArrayList<String>();

		//TODO: Update to more efficient Annotation Loading for entire file
		markerAnnotationLoader.fillAnnotations(null, parsers, QUERY_ORDER.ONE_PER_IN_ORDER);

		for (int i = 0; i < blastResults.length; i++) {
			MarkerBlastAnnotation current = blastResults[i];
			ArrayList<BlastAnnotation> perfectMatches = current.getAnnotationsFor(BLAST_ANNOTATION_TYPES.PERFECT_MATCH, log);
			if (perfectMatches.size() == 1) {
				int[] alignmentHistogram = current.getAlignmentHistogram(proj);
				int sub = (int) Math.round((double)crossHybePercent*alignmentHistogram.length);
				int numHits = Array.sum(Array.subArray(alignmentHistogram, sub));
				if (numHits == 1){
					// Perfect match is the only hit within crossHybePercent
					oneHitters.add(markerNames[i]);
				} else if (numHits == 2) {
					PROBE_TAG perfectTag = perfectMatches.get(0).getTag();
					if (perfectTag != PROBE_TAG.BOTH) {
						for (BlastAnnotation annotation : current.getAnnotationsFor(BLAST_ANNOTATION_TYPES.ON_T_ALIGNMENTS_NON_PERFECT, log)) {
							if (annotation.getTag() != perfectTag && annotation.getRefLoc().getSize() == proj.getArrayType().getProbeLength() - 1) {
								// Second match is an on target read to the second probe, missing only the target base
								oneHitters.add(markerNames[i]);
								break;
							}
						}
					}
				}
			}
		}
		
		log.reportTime(oneHitters.size() + " one hit wonder markers identified out of " + markerNames.length + " total markers");
		log.report("Writing results to " + outFile);
		Files.writeArrayList(oneHitters, outFile);
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;
		String blastVCF = null;
		double crossHybridizationThreshold = DEFAULT_CROSS_HYBE_THRESHOLD;
		String usage = "\n" + "cnv.qc.MarkerBlast requires 3 arguments\n";
		usage += "   (1) Project file name (i.e. proj=" + filename + " (default))\n" + "";
		usage += "   (2) Blast.vcf filename  (i.e. blastVCF=" + "./blast.vcf.gz " + " (default based on project properties))\n" + "";
		usage += "   (3) Cross hybridization threshold  (i.e. crossHybridizationThreshold=" + crossHybridizationThreshold + " (default))\n" + "";


		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("blastVCF=")) {
				blastVCF = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("crossHybridizationThreshold=")) {
				crossHybridizationThreshold = ext.parseDoubleArg(args[i]);
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			Project proj = null;

			if (filename == null) {
				System.err.println(usage);
				System.exit(1);
			} else {
				proj = new Project(filename, false);
			}
			if (proj == null) {
				System.err.println("Invalid project");
				System.exit(1);
			}
			if (blastVCF == null) {
				blastVCF = proj.BLAST_ANNOTATION_FILENAME.getValue();
			}
			getOneHitWonders(proj, blastVCF, defaultOneHitWondersFilename(blastVCF), crossHybridizationThreshold, proj.getLog());
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
