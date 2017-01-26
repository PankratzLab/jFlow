package org.genvisis.seq.qc;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.filesys.GeneData;
import org.genvisis.filesys.GeneTrack;
import org.genvisis.filesys.LocusSet;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.manage.ReferenceGenome;
import org.genvisis.stats.Rscript.GeomText;
import org.genvisis.stats.Rscript.RScatter;
import org.genvisis.stats.Rscript.SCATTER_TYPE;

/**
 * @author lane0212 Used for determining capture effeciency of specific targets
 */
public class CaptureQC {

	private static final String UCSC = "UCSC";
	private static final String[] PLOT_BY_POS_PERCENT = new String[] {"Percent_Covered_at_depth_1",
																																		"Percent_Covered_at_depth_10",
																																		"Percent_Covered_at_depth_20",
																																		"Percent_Covered_at_depth_30",
																																		"Percent_Covered_at_depth_40"};
	private static final String[] PLOT_BY_POS_ACTUAL = new String[] {"averageCoverage"};

	public static void captureQC(	String referenceGenomeFasta, String bamQCSummary,
																String extraPostionFile, String geneTrackFile, String[] geneNames,
																String outputDir, String root, boolean allInOne, Logger log) {
		if (extraPostionFile != null) {
			GeomText.fromFile(extraPostionFile, log);
		}

		GeneTrack geneTrack = GeneTrack.load(geneTrackFile, false);
		ReferenceGenome referenceGenome = new ReferenceGenome(referenceGenomeFasta, log);
		String output = outputDir + root + "fullSummary.capture.txt";

		if (!allInOne || !Files.exists(output)) {

			ArrayList<GeneData> gds = new ArrayList<GeneData>();
			for (String geneName : geneNames) {
				GeneData[] tmp = geneTrack.lookupAllGeneData(geneName);
				for (GeneData element : tmp) {
					gds.add(element);
				}
			}
			LocusSet<GeneData> geneSet = new LocusSet<GeneData>(gds.toArray(new GeneData[gds.size()]),
																													true, log) {

				/**
				* 
				*/
				private static final long serialVersionUID = 1L;

			};

			// for (int i = 0; i < geneNames.length; i++) {

			try {
				// ArrayList<RScatter> rsScatters = new ArrayList<RScatter>();
				PrintWriter writer = new PrintWriter(new FileWriter(output, false));
				String[] header = Files.getHeaderOfFile(bamQCSummary, log);
				writer.println("GENE_NAME\tExon\tPosition\tGC_REF\tHIDE"	+ ArrayUtils.toStr(header)
												+ "\tInternalKey");

				int UCSCIndex = ext.indexOfStr(UCSC, header);
				BufferedReader reader = Files.getAppropriateReader(bamQCSummary);
				reader.readLine();
				while (reader.ready()) {

					String[] line = reader.readLine().trim().split("\t");
					Segment seg = new Segment(line[UCSCIndex]);
					GeneData[] overlapping = geneSet.getOverLappingLoci(seg);
					if (overlapping != null) {
						double gcRegion = referenceGenome.getGCContentFor(seg);
						for (GeneData currentData : overlapping) {
							for (int j = 0; j < currentData.getExonBoundaries().length; j++) {
								Segment exon = new Segment(	currentData.getChr(),
																						currentData.getExonBoundaries()[j][0],
																						currentData.getExonBoundaries()[j][1]);
								if (seg.overlaps(exon)) {
									if (allInOne) {
										writer.println(currentData.getGeneName()	+ "\t" + (j + 1) + "\t"
																		+ seg.getStart() + "\t" + gcRegion + "\t" + ArrayUtils.toStr(line)
																		+ "\t0");
									} else {// print bp resolution
										for (int k = seg.getStart(); k <= seg.getStop(); k++) {
											if (exon.getStart() <= k && exon.getStop() >= k) {
												writer.println(currentData.getGeneName()	+ "\t" + (j + 1) + "\t" + k + "\t"
																				+ gcRegion + "\t" + ArrayUtils.toStr(line) + "\t" + j + "_" + k);
											}
										}
									}
								}
							}
						}
					}
				}
				reader.close();
				writer.close();
			} catch (FileNotFoundException fnfe) {
				log.reportError("Error: file \"" + bamQCSummary + "\" not found in current directory");
				return;
			} catch (IOException ioe) {
				log.reportError("Error reading file \"" + bamQCSummary + "\"");
				return;
			}

		}

		plot(new String[] {"Full Summary"}, log, null, 0, output, true);

	}

	private static void plot(	String[] geneNames, Logger log, GeomText[] geomTexts, int i,
														String output, boolean allinone) {

		if (allinone) {
			for (int j = 0; j < PLOT_BY_POS_PERCENT.length; j++) {

				double[] data = ArrayUtils.toDoubleArray(HashVec.loadFileToStringArray(output, true,
																																					new int[] {ext.indexOfStr(PLOT_BY_POS_PERCENT[j],
																																																		Files.getHeaderOfFile(output,
																																																													log))},
																																					false));
				String[] genes = HashVec.loadFileToStringArray(	output, true,
																												new int[] {ext.indexOfStr("GENE_NAME",
																																									Files.getHeaderOfFile(output,
																																																				log))},
																												false);

				double average = ArrayUtils.mean(data, true);
				String rootExon = ext.rootOf(output, false) + "_coverageHist_" + j;
				RScatter rsScatterPos = new RScatter(	output, rootExon + ".rscript",
																							ext.removeDirectoryInfo(rootExon), rootExon + ".jpeg",
																							"InternalKey", new String[] {PLOT_BY_POS_PERCENT[j]},
																							SCATTER_TYPE.HIST, log);
				String format = PLOT_BY_POS_PERCENT[j].replaceAll("_", " ")
																							.replaceAll("Percent", "Proportion of target ");
				rsScatterPos.setxRange(new double[] {0, 1});
				rsScatterPos.setTitle("Average "	+ format + "  (" + ArrayUtils.unique(genes).length + " genes, "
															+ data.length + " targeted regions)" + " = "
															+ ext.formDeci(average, 3));
				rsScatterPos.setxLabel(format);
				rsScatterPos.setyLabel("Counts (targeted regions) labelled with percent of total");
				rsScatterPos.setOverWriteExisting(true);
				rsScatterPos.execute();
			}
			// rsScatterPosos.setVertLines(new VertLine[]{new VertLine()});

		} else {
			String root = ext.rootOf(output, false);
			RScatter rsScatterPos = new RScatter(	output, root + ".rscript", ext.removeDirectoryInfo(root),
																						root + ".pdf", "Position", PLOT_BY_POS_PERCENT,
																						SCATTER_TYPE.POINT, log);
			rsScatterPos.setTitle(geneNames[i] + " capture");
			rsScatterPos.setxLabel("Position");
			rsScatterPos.setyLabel("Proportion Covered at Depth");
			if (geomTexts != null) {
				rsScatterPos.setgTexts(geomTexts);
				rsScatterPos.setyRange(new double[] {0, 2});

			}
			rsScatterPos.setOverWriteExisting(true);
			rsScatterPos.setWidth(25);
			rsScatterPos.execute();

			String rootExon = ext.rootOf(output, false) + "_coverage";

			RScatter rsScatterActual = new RScatter(output, rootExon + ".rscript",
																							ext.removeDirectoryInfo(rootExon), rootExon + ".pdf",
																							"Position", PLOT_BY_POS_ACTUAL, SCATTER_TYPE.POINT,
																							log);
			rsScatterActual.setTitle(geneNames[i] + " capture");
			rsScatterActual.setxLabel("Position");
			rsScatterActual.setyLabel("Average Coverage");
			if (geomTexts != null) {
				for (GeomText geomText : geomTexts) {
					geomText.setY(geomText.getY() * 10 + 80);
				}
				rsScatterActual.setyRange(new double[] {0, 100});
				rsScatterActual.setgTexts(geomTexts);
			}
			rsScatterActual.setOverWriteExisting(true);
			rsScatterActual.setWidth(25);
			rsScatterActual.execute();
		}
	}

	// private static class ExtraPositions {
	// private static final String[] HEADER = new String[] { "UCSC", "GROUP", "VALUE" };
	// private String file;
	//
	// }

	public static void main(String[] args) {
		int numArgs = args.length;
		String bamQCSummary = "rrd_bamQC.targets.libraryResults.summary.txt";
		String geneTrackFile = "/panfs/roc/groups/5/pankrat2/public/bin/NCBI/RefSeq_hg19.gtrack";
		String referenceGenomeFasta = "/panfs/roc/groups/5/pankrat2/public/bin/ref/hg19_canonical.fa";
		String[] geneNames = new String[] {"TP53"};
		String outputDir = null;
		String root = "GeneSummary";
		// String logfile = null;
		String extraPostionFile = null;
		boolean allInOne = true;
		String usage = "\n" + "seq.qc.CaptureQC requires 0-1 arguments\n";
		usage += "   (1) bamQC file (i.e. bamQCSummary=" + bamQCSummary + " (default))\n" + "";
		usage += "   (2) gene track file(i.e. geneTrackFile=" + geneTrackFile + " (default))\n" + "";
		usage += "   (3) comma delimited gene names(i.e. geneNames="	+ ArrayUtils.toStr(geneNames, ",")
							+ " (default))\n" + "";
		usage += "   (4) output directory (i.e. outputDir= ( no default))\n" + "";
		usage += "   (5) reference genome fasta file (i.e. referenceGenomeFasta="	+ referenceGenomeFasta
							+ "  (default))\n" + "";
		usage += "   (6) extra positions file  (i.e. extrap= ( no default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("bamQCSummary=")) {
				bamQCSummary = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("geneTrackFile=")) {
				geneTrackFile = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("geneNames=")) {
				geneNames = arg.split("=")[1].split(",");
				numArgs--;
			} else if (arg.startsWith("outputDir=")) {
				outputDir = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("extrap=")) {
				extraPostionFile = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("root=")) {
				root = arg.split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + arg);
			}
		}
		if (numArgs != 0) {

			System.err.println(usage);
			System.exit(1);
		}
		try {
			if (outputDir == null) {
				outputDir = ext.parseDirectoryOfFile(bamQCSummary);
			}
			Logger log = new Logger(outputDir + "captureQC.log");
			if (Files.exists(geneNames[0])) {
				log.reportTimeInfo("Detected file input for gene names, loading " + geneNames[0]);
				geneNames = HashVec.loadFileToStringArray(geneNames[0], false, new int[] {0}, true);
				geneNames = ArrayUtils.unique(geneNames);
				log.reportTimeInfo("Loaded " + geneNames.length + " genes");
			}
			captureQC(referenceGenomeFasta, bamQCSummary, extraPostionFile, geneTrackFile, geneNames,
								outputDir, root, allInOne, log);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
