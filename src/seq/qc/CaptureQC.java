package seq.qc;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import seq.manage.ReferenceGenome;
import stats.Rscript.GeomText;
import stats.Rscript.RScatter;
import stats.Rscript.SCATTER_TYPE;
import cnv.var.LocusSet;
import common.Array;
import common.Files;
import common.Logger;
import common.ext;
import filesys.GeneData;
import filesys.GeneTrack;
import filesys.Segment;

/**
 * @author lane0212 Used for determining capture effeciency of specific targets
 */
public class CaptureQC {

	private static final String UCSC = "UCSC";
	private static final String[] PLOT_BY_POS_PERCENT = new String[] { "Percent_Covered_at_depth_1", "Percent_Covered_at_depth_10", "Percent_Covered_at_depth_20" };
	private static final String[] PLOT_BY_POS_ACTUAL = new String[] { "averageCoverage" };

	public static void captureQC(String referenceGenomeFasta, String bamQCSummary, String extraPostionFile, String geneTrackFile, String[] geneNames, String outputDir, Logger log) {
		GeomText[] geomTexts = null;
		if (extraPostionFile != null) {
			geomTexts = GeomText.fromFile(extraPostionFile, log);
		}
	
		GeneTrack geneTrack = GeneTrack.load(geneTrackFile, false);
		ReferenceGenome referenceGenome = new ReferenceGenome(referenceGenomeFasta, log);
		for (int i = 0; i < geneNames.length; i++) {
			String output = outputDir + geneNames[i] + ".capture.txt";
			GeneData[] geneData = geneTrack.lookupAllGeneData(geneNames[i]);
			LocusSet<GeneData> geneSet = new LocusSet<GeneData>(geneData, true, log) {

				/**
				 * 
				 */
				private static final long serialVersionUID = 1L;

			};

			try {
				//ArrayList<RScatter> rsScatters = new ArrayList<RScatter>();
				PrintWriter writer = new PrintWriter(new FileWriter(output));
				String[] header = Files.getHeaderOfFile(bamQCSummary, log);
				writer.println("GENE_NAME\tExon\tPosition\tGC_REF\tHIDE" + Array.toStr(header));
				int UCSCIndex = ext.indexOfStr(UCSC, header);
				BufferedReader reader = Files.getAppropriateReader(bamQCSummary);
				reader.readLine();
				while (reader.ready()) {
					String[] line = reader.readLine().trim().split("\t");
					Segment seg = new Segment(line[UCSCIndex]);
					GeneData[] overlapping = geneSet.getOverLappingLoci(seg);
					if (overlapping != null) {
						double gcRegion = referenceGenome.getGCContentFor(seg);
						for (int geneIndex = 0; geneIndex < overlapping.length; geneIndex++) {
							GeneData currentData = overlapping[geneIndex];
							for (int j = 0; j < currentData.getExonBoundaries().length; j++) {
								Segment exon = new Segment(currentData.getChr(), currentData.getExonBoundaries()[j][0], currentData.getExonBoundaries()[j][1]);
								if (seg.overlaps(exon)) {
									for (int k = seg.getStart(); k <= seg.getStop(); k++) {
										if (exon.getStart() <= k && exon.getStop() >= k) {
											writer.println(currentData.getGeneName() + "\t" + (j + 1) + "\t" + k + "\t" + gcRegion + "\t" + Array.toStr(line));
										}
									}
								}
							}
						}
					}
				}
				reader.close();
				writer.close();
				String root = ext.rootOf(output, false);
				RScatter rsScatterPos = new RScatter(output, root + ".rscript", ext.removeDirectoryInfo(root), root + ".pdf", "Position", PLOT_BY_POS_PERCENT, SCATTER_TYPE.POINT, log);
				rsScatterPos.setTitle(geneNames[i] + " capture");
				rsScatterPos.setxLabel("Position");
				rsScatterPos.setyLabel("Proportion Covered at Depth");
				if (geomTexts != null) {
					rsScatterPos.setgTexts(geomTexts);
					rsScatterPos.setyRange(new double[] { 0, 2 });

				}
				rsScatterPos.setOverWriteExisting(true);
				rsScatterPos.setWidth(25);
				rsScatterPos.execute();

				String rootExon = ext.rootOf(output, false) + "_coverage";

				RScatter rsScatterActual = new RScatter(output, rootExon + ".rscript", ext.removeDirectoryInfo(rootExon), rootExon + ".pdf", "Position", PLOT_BY_POS_ACTUAL, SCATTER_TYPE.POINT, log);
				rsScatterActual.setTitle(geneNames[i] + " capture");
				rsScatterActual.setxLabel("Position");
				rsScatterActual.setyLabel("Average Coverage");
				if (geomTexts != null) {
					for (int j = 0; j < geomTexts.length; j++) {
						geomTexts[j].setY(geomTexts[j].getY()*10 + 80);
					}
					rsScatterActual.setyRange(new double[] { 0, 100 });
					rsScatterActual.setgTexts(geomTexts);
				}
				rsScatterActual.setOverWriteExisting(true);
				rsScatterActual.setWidth(25);
				rsScatterActual.execute();

//				
//				rsScatters.add(rsScatterPos);
//				rsScatters.add(rsScatterActual);
//				RScatters rScatters = new RScatters(rsScatters.toArray(new RScatter[rsScatters.size()]), output + ".rscript", output + ".pdf", COLUMNS_MULTIPLOT.COLUMNS_MULTIPLOT_1, PLOT_DEVICE.PDF, log);
//				rScatters.execute();

			} catch (FileNotFoundException fnfe) {
				log.reportError("Error: file \"" + bamQCSummary + "\" not found in current directory");
				return;
			} catch (IOException ioe) {
				log.reportError("Error reading file \"" + bamQCSummary + "\"");
				return;
			}
		}

	}

//	private static class ExtraPositions {
//		private static final String[] HEADER = new String[] { "UCSC", "GROUP", "VALUE" };
//		private String file;
//
//	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String bamQCSummary = "rrd_bamQC.targets.libraryResults.summary.txt";
		String geneTrackFile = "/panfs/roc/groups/5/pankrat2/public/bin/NCBI/RefSeq_hg19.gtrack";
		String referenceGenomeFasta = "/panfs/roc/groups/5/pankrat2/public/bin/ref/hg19_canonical.fa";
		String[] geneNames = new String[] { "TP53" };
		String outputDir = null;
		//String logfile = null;
		String extraPostionFile = null;

		String usage = "\n" + "seq.qc.CaptureQC requires 0-1 arguments\n";
		usage += "   (1) bamQC file (i.e. bamQCSummary=" + bamQCSummary + " (default))\n" + "";
		usage += "   (2) gene track file(i.e. geneTrackFile=" + geneTrackFile + " (default))\n" + "";
		usage += "   (3) comma delimited gene names(i.e. geneNames=" + Array.toStr(geneNames, ",") + " (default))\n" + "";
		usage += "   (4) output directory (i.e. outputDir= ( no default))\n" + "";
		usage += "   (5) reference genome fasta file (i.e. referenceGenomeFasta=" + referenceGenomeFasta + "  (default))\n" + "";
		usage += "   (6) extra positions file  (i.e. extrap= ( no default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("bamQCSummary=")) {
				bamQCSummary = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("geneTrackFile=")) {
				geneTrackFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("geneNames=")) {
				geneNames = args[i].split("=")[1].split(",");
				numArgs--;
			} else if (args[i].startsWith("outputDir=")) {
				outputDir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("extrap=")) {
				extraPostionFile = args[i].split("=")[1];
				numArgs--;
			}  else {
				System.err.println("Error - invalid argument: " + args[i]);
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
			captureQC(referenceGenomeFasta, bamQCSummary, extraPostionFile, geneTrackFile, geneNames, outputDir, log);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
