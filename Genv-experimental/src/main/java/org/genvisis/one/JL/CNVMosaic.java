package org.genvisis.one.JL;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Set;
import java.util.concurrent.Callable;

import org.genvisis.cnv.analysis.BeastScore;
import org.genvisis.cnv.analysis.MosaicismDetect;
import org.genvisis.cnv.analysis.MosaicismDetect.MosaicBuilder;
import org.genvisis.cnv.analysis.MosaicismDetect.MosaicProducer;
import org.genvisis.cnv.filesys.MarkerSet;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.var.MosaicRegion;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.WorkerTrain;
import org.genvisis.common.WorkerTrain.AbstractProducer;
import org.genvisis.common.ext;
import org.genvisis.filesys.CNVariant;
import org.genvisis.filesys.LocusSet;
import org.genvisis.filesys.Segment;
import org.genvisis.stats.Rscript.COLUMNS_MULTIPLOT;
import org.genvisis.stats.Rscript.PLOT_DEVICE;
import org.genvisis.stats.Rscript.RScatter;
import org.genvisis.stats.Rscript.RScatters;
import org.genvisis.stats.Rscript.Restrictions;
import org.genvisis.stats.Rscript.SCATTER_TYPE;

public class CNVMosaic {

	private static LocusSet<MosaicRegion> callMosaic(Project proj, String output, int numThreads) {

		LocusSet<MosaicRegion> results = null;
		SampleData sampleData = proj.getSampleData(0, false);
		ArrayList<String> samples = new ArrayList<String>();
		String[] samps = proj.getSamples();
		for (int i = 0; i < samps.length; i++) {
			if (!sampleData.individualShouldBeExcluded(samps[i])) {
				samples.add(samps[i]);
			}
		}
		MarkerSet markerSet = proj.getMarkerSet();
		int[][] indicesByChr = markerSet.getIndicesByChr();
		MosaicBuilder builder = new MosaicBuilder();
		builder.indicesByChr(indicesByChr);

		ArrayList<Segment> callSegs = new ArrayList<Segment>();

		for (int i = 0; i < indicesByChr.length; i++) {
			if (i > 0 && i < 23 && indicesByChr[i].length > 0) {
				callSegs.add(new Segment((byte) i, 0, Integer.MAX_VALUE));
			}
		}
		LocusSet<Segment> segs = new LocusSet<Segment>(	callSegs.toArray(new Segment[callSegs.size()]),
																										true, proj.getLog()) {

			/**
			 *
			 */
			private static final long serialVersionUID = 1L;

		};
		String ser = output + ".ser";

		if (!Files.exists(ser) || !Files.exists(output)) {
			ArrayList<MosaicRegion> all = new ArrayList<MosaicRegion>();

			try {
				PrintWriter writer = new PrintWriter(new FileWriter(output));
				writer.println(Array.toStr(CNVariant.PLINK_CNV_HEADER)	+ "\t"
												+ Array.toStr(MosaicRegion.ADD_HEADER) + "\tEXCLUDED\tDNA\tUCSC");

				MosaicProducer producer = new MosaicProducer(	proj, builder, Array.toStringArray(samples),
																											markerSet, segs);
				WorkerTrain<LocusSet<MosaicRegion>> train =
																									new WorkerTrain<LocusSet<MosaicRegion>>(producer,
																																													numThreads,
																																													2,
																																													proj.getLog());
				int index = 0;
				while (train.hasNext()) {
					index++;
					LocusSet<MosaicRegion> tmp = train.next();
					tmp.addAll(all);
					for (int i = 0; i < tmp.getLoci().length; i++) {
						String sampKey = sampleData.lookup(tmp.getLoci()[i].getFamilyID()	+ "\t"
																								+ tmp.getLoci()[i].getIndividualID())[0];
						writer.println(tmp.getLoci()[i].toAnalysisString()	+ "\t"
														+ sampleData.individualShouldBeExcluded(sampKey) + "\t" + sampKey + "\t"
														+ tmp.getLoci()[i].getUCSClocation());
					}
					proj.getLog().reportTimeInfo("Calling mos for sample " + index + " of " + samples.size());
				}
				writer.close();
			} catch (Exception e) {
				proj.getLog().reportError("Error writing to " + output);
				proj.getLog().reportException(e);
			}
			results = new LocusSet<MosaicRegion>(	all.toArray(new MosaicRegion[all.size()]), true,
																						proj.getLog()) {

				/**
				 *
				 */
				private static final long serialVersionUID = 1L;

			};
			results.writeSerial(ser);
		} else {
			results = LocusSet.readSerialMRSet(ser, proj.getLog());
		}
		return null;

	}

	private static void computeCNVMosaic(Project proj) {

		String mosaicOut = proj.PROJECT_DIRECTORY.getValue() + "mosaicCalls.mos";
		callMosaic(proj, mosaicOut, proj.NUM_THREADS.getValue());
		int[] cns = new int[] {0};
		Restrictions[] restrictions = new Restrictions[] {new Restrictions(	new String[] {"SITES",},
																																				new double[] {75},
																																				new String[] {">",}, null)};
		plot(proj, mosaicOut, cns, restrictions, true);
		// plot(proj, mosaicOut);
		String[] cnvFiles = proj.CNV_FILENAMES.getValue();
		SampleData sampleData = proj.getSampleData(0, false);
		for (String cnvFile : cnvFiles) {
			String output = ext.addToRoot(cnvFile, ".mosaicMetrics");
			if (!Files.exists(output)) {
				LocusSet<CNVariant> cnvs = CNVariant.loadLocSet(cnvFile, proj.getLog());
				Hashtable<String, LocusSet<CNVariant>> indSets = CNVariant.breakIntoInds(	cnvs,
																																									proj.getLog());
				int numSampsRemoved = 0;
				int numCNVsRemoved = 0;
				ArrayList<String> remove = new ArrayList<String>();
				for (String ind : indSets.keySet()) {
					String sampKey = sampleData.lookup(ind)[0];
					if (sampleData.individualShouldBeExcluded(sampKey)) {
						numSampsRemoved++;
						numCNVsRemoved += indSets.get(ind).getLoci().length;
						remove.add(ind);
					}
				}
				proj.getLog().reportTimeInfo(indSets.size() + " current samples");
				for (String remover : remove) {
					indSets.remove(remover);
				}
				proj.getLog().reportTimeInfo(indSets.size() + " current samples after excludes");

				proj.getLog().reportTimeInfo("Removed "	+ numSampsRemoved + "  samples and "
																			+ numCNVsRemoved + " cnvs");
				MosaicForceProducer producer = new MosaicForceProducer(proj, indSets);
				WorkerTrain<MosaicRegion[]> train = new WorkerTrain<MosaicRegion[]>(producer,
																																						proj.NUM_THREADS.getValue(),
																																						2, proj.getLog());

				try {
					PrintWriter writer = new PrintWriter(new FileWriter(output));
					writer.println(Array.toStr(CNVariant.PLINK_CNV_HEADER)	+ "\t"
													+ Array.toStr(MosaicRegion.ADD_HEADER) + "\tEXCLUDED\tDNA\tUCSC");
					int index = 0;
					while (train.hasNext()) {
						MosaicRegion[] regions = train.next();
						for (MosaicRegion region : regions) {
							String sampKey = sampleData.lookup(region.getFamilyID()	+ "\t"
																									+ region.getIndividualID())[0];
							writer.println(region.toAnalysisString()	+ "\t"
															+ sampleData.individualShouldBeExcluded(sampKey) + "\t" + sampKey
															+ "\t" + region.getUCSClocation());
						}
						index++;
						proj.getLog().reportTimeInfo(index + " of " + indSets.size());
					}
					writer.close();
				} catch (Exception e) {
					proj.getLog().reportError("Error writing to " + output);
					proj.getLog().reportException(e);
				}
			}

		}

	}

	private static void plot(	Project proj, String output, int[] cns, Restrictions[] restrictions,
														boolean mos) {

		ArrayList<RScatter> rscatters = new ArrayList<RScatter>();

		for (int cn : cns) {
			String outPlot8 = ext.addToRoot(output, ".markMarkbox_" + cn);
			RScatter rScatterMarkBOX = new RScatter(output, outPlot8 + ".rscript",
																							ext.removeDirectoryInfo(outPlot8), outPlot8 + ".jpeg",
																							"numCustomFMarkers", new String[] {"customF"}, null,
																							SCATTER_TYPE.BOX, proj.getLog());
			rScatterMarkBOX.setyLabel("customF");
			rScatterMarkBOX.setxLabel("NumCustomF Markers");
			rScatterMarkBOX.setTitle("customF v NumCustomF Markers, copy number " + cn);
			Restrictions[] restrictionss = new Restrictions[] {new Restrictions(new String[] {"TYPE"},
																																					new double[] {cn},
																																					new String[] {"=="},
																																					null)};
			rScatterMarkBOX.setRestrictions(restrictionss);

			String outPlot9 = ext.addToRoot(output, ".markMarkbox_limit.d" + cn);
			RScatter rScatterMarkBOXLimit = new RScatter(	output, outPlot9 + ".rscript",
																										ext.removeDirectoryInfo(outPlot9),
																										outPlot9 + ".jpeg", "numCustomFMarkers",
																										new String[] {"customF"}, null,
																										SCATTER_TYPE.BOX, proj.getLog());
			rScatterMarkBOXLimit.setyLabel("customF");
			rScatterMarkBOXLimit.setxLabel("NumCustomF Markers");
			rScatterMarkBOXLimit.setTitle("customF v NumCustomF Markers, copy number " + cn);
			Restrictions[] restrictions2 = new Restrictions[] {new Restrictions(new String[] {"TYPE",
																																												"SITES"},
																																					new double[] {cn, 20},
																																					new String[] {"==", "<"},
																																					"&")};

			rScatterMarkBOXLimit.setRestrictions(restrictions2);

			// rScatterMarkMarkType.setScaleDensity(true);
			// rScatterMarkBOX.setyRange(new double[] { 0, 100 });
			// rScatterMarkBOX.setxRange(new double[] { 0, 200 });
			// rScatterMarkBOX.setOverWriteExisting(true);
			// rScatterMarkBOX.execute();
			rscatters.add(rScatterMarkBOX);
			rscatters.add(rScatterMarkBOXLimit);

		}

		String outPlot1 = ext.addToRoot(output, ".height");
		RScatter rScatterHeight = new RScatter(	output, outPlot1 + ".rscript",
																						ext.removeDirectoryInfo(outPlot1), outPlot1 + ".jpeg",
																						"customF", new String[] {"beastHeight"}, "TYPE",
																						SCATTER_TYPE.POINT, proj.getLog());
		rScatterHeight.setyLabel("BEAST_HEIGHT");
		rScatterHeight.setxLabel("CustomF");
		rScatterHeight.setTitle("CUSTOM_F vs BEAST Height");
		rScatterHeight.setRestrictions(restrictions);
		rScatterHeight.setOverWriteExisting(true);
		rScatterHeight.execute();

		if (mos) {
			String outPlotMos = ext.addToRoot(output, ".heightMos");
			RScatter rScatterHeightMos = new RScatter(output, outPlotMos + ".rscript",
																								ext.removeDirectoryInfo(outPlotMos),
																								outPlotMos + ".jpeg", "customF",
																								new String[] {"beastHeight"}, "bpWeightedScore",
																								SCATTER_TYPE.POINT, proj.getLog());
			rScatterHeightMos.setyLabel("BEAST_HEIGHT");
			rScatterHeightMos.setxLabel("CustomF");
			rScatterHeightMos.setTitle("CUSTOM_F vs BEAST Height v bpWeightedScore");
			rScatterHeightMos.setScaleDensity(true);
			rScatterHeightMos.setRestrictions(restrictions);
			rScatterHeightMos.setMidScaleDensity(10);
			rScatterHeightMos.setOverWriteExisting(true);
			rScatterHeightMos.execute();
			rscatters.add(rScatterHeightMos);
		}

		rscatters.add(rScatterHeight);
		rScatterHeight.setOverWriteExisting(true);
		rScatterHeight.execute();
		String outPlot3 = ext.addToRoot(output, ".score");
		RScatter rScatterScore = new RScatter(output, outPlot3 + ".rscript",
																					ext.removeDirectoryInfo(outPlot3), outPlot3 + ".jpeg",
																					"customF", new String[] {"beastScore"}, "TYPE",
																					SCATTER_TYPE.POINT, proj.getLog());
		rScatterScore.setRestrictions(restrictions);

		rScatterScore.setyLabel("BEAST_SCORE");
		rScatterScore.setxLabel("CustomF");
		rScatterScore.setTitle("CUSTOM_F vs BEAST SCORE");
		rscatters.add(rScatterScore);

		String outPlot4 = ext.addToRoot(output, ".penn");
		RScatter rScatterPennScore = new RScatter(output, outPlot4 + ".rscript",
																							ext.removeDirectoryInfo(outPlot4), outPlot4 + ".jpeg",
																							"customF", new String[] {"SCORE"}, "TYPE",
																							SCATTER_TYPE.POINT, proj.getLog());
		rScatterPennScore.setRestrictions(restrictions);

		rScatterPennScore.setyLabel("PENNCNV SCORE");
		rScatterPennScore.setxLabel("CustomF");
		rScatterPennScore.setTitle("CUSTOM_F vs PENNCNV SCORE");
		rscatters.add(rScatterPennScore);

		String outPlot2 = ext.addToRoot(output, ".markers");
		RScatter rScatterMarkers = new RScatter(output, outPlot2 + ".rscript",
																						ext.removeDirectoryInfo(outPlot2), outPlot2 + ".jpeg",
																						"customF", new String[] {"numCustomFMarkers"}, "TYPE",
																						SCATTER_TYPE.POINT, proj.getLog());
		rScatterMarkers.setRestrictions(restrictions);

		rScatterMarkers.setyLabel("NumCustomF Markers");
		rScatterMarkers.setxLabel("CustomF");
		rScatterMarkers.setTitle("CUSTOM_F vs Markers");
		rscatters.add(rScatterMarkers);

		String outPlot7 = ext.addToRoot(output, ".markersPenn");
		RScatter rScatterMarkersPenn = new RScatter(output, outPlot7 + ".rscript",
																								ext.removeDirectoryInfo(outPlot7),
																								outPlot7 + ".jpeg", "customF",
																								new String[] {"SITES"}, "TYPE", SCATTER_TYPE.POINT,
																								proj.getLog());
		rScatterMarkersPenn.setRestrictions(restrictions);

		rScatterMarkersPenn.setyLabel("SITES Markers");
		rScatterMarkersPenn.setxLabel("CustomF");
		rScatterMarkersPenn.setTitle("CUSTOM_F vs SITES");
		rscatters.add(rScatterMarkersPenn);

		String outPlot5 = ext.addToRoot(output, ".markMark");
		RScatter rScatterMarkMark = new RScatter(	output, outPlot5 + ".rscript",
																							ext.removeDirectoryInfo(outPlot5), outPlot5 + ".jpeg",
																							"SITES", new String[] {"numCustomFMarkers"},
																							"customF", SCATTER_TYPE.POINT, proj.getLog());
		rScatterMarkMark.setRestrictions(restrictions);

		rScatterMarkMark.setyLabel("NumCustomF Markers");
		rScatterMarkMark.setxLabel("PennMarkers");
		rScatterMarkMark.setTitle("CUSTOM_F vs PennMarkers v NumCustomF Markers");
		rScatterMarkMark.setScaleDensity(true);

		rScatterMarkMark.setyRange(new double[] {0, 100});
		rScatterMarkMark.setxRange(new double[] {0, 200});
		rscatters.add(rScatterMarkMark);

		String outPlot6 = ext.addToRoot(output, ".markMarkType");
		RScatter rScatterMarkMarkType = new RScatter(	output, outPlot6 + ".rscript",
																									ext.removeDirectoryInfo(outPlot6),
																									outPlot6 + ".jpeg", "SITES",
																									new String[] {"numCustomFMarkers"}, "TYPE",
																									SCATTER_TYPE.POINT, proj.getLog());
		rScatterMarkMarkType.setRestrictions(restrictions);

		rScatterMarkMarkType.setyLabel("NumCustomF Markers");
		rScatterMarkMarkType.setxLabel("PennMarkers");
		rScatterMarkMarkType.setTitle("TYPE vs PennMarkers v NumCustomF Markers");
		// rScatterMarkMarkType.setScaleDensity(true);
		rScatterMarkMarkType.setyRange(new double[] {0, 100});
		rScatterMarkMarkType.setxRange(new double[] {0, 200});

		RScatters rs = new RScatters(	rscatters.toArray(new RScatter[rscatters.size()]),
																	output + ".rscript", output + "final.pdf",
																	COLUMNS_MULTIPLOT.COLUMNS_MULTIPLOT_1, PLOT_DEVICE.PDF,
																	proj.getLog());
		rs.execute();
	}

	private static class MosaicForceProducer extends AbstractProducer<MosaicRegion[]> {
		private final Project proj;
		private final Hashtable<String, LocusSet<CNVariant>> indSets;
		private int index;
		private final String[] fidsIids;
		private final SampleData sampleData;
		// private int[][] indicesByChr;

		public MosaicForceProducer(Project proj, Hashtable<String, LocusSet<CNVariant>> indSets) {
			super();
			this.proj = proj;
			this.indSets = indSets;
			index = 0;
			Set<String> tmp = indSets.keySet();
			fidsIids = new String[tmp.size()];
			int i = 0;
			for (String fidIid : tmp) {
				fidsIids[i] = fidIid;
				i++;
			}
			sampleData = proj.getSampleData(0, false);
		}

		@Override
		public boolean hasNext() {

			return index < fidsIids.length;
		}

		@Override
		public Callable<MosaicRegion[]> next() {
			final String sample = sampleData.lookup(fidsIids[index])[0];
			final LocusSet<CNVariant> sampCNVs = indSets.get(fidsIids[index]);

			Callable<MosaicRegion[]> callable = new Callable<MosaicRegion[]>() {

				@Override
				public MosaicRegion[] call() throws Exception {
					MarkerSet markerSet = proj.getMarkerSet();
					int[][] indicesByChr = markerSet.getIndicesByChr();
					Sample samp = proj.getFullSampleFromRandomAccessFile(sample);
					ArrayList<MosaicRegion> all = new ArrayList<MosaicRegion>();
					MosaicBuilder builderMosaic = new MosaicBuilder();
					builderMosaic.verbose(true);
					MosaicismDetect md = builderMosaic.build(	proj, sample, markerSet,
																										Array.toDoubleArray(samp.getBAFs()));
					// MosaicQuantWorker worker = new MosaicQuantWorker(sampCNVs.getLoci(), proj, sample,
					// MOSAIC_TYPE.values(), 5);
					// WorkerHive<MosaicQuantResults[]> hive = new
					// WorkerHive<MosaicismQuant.MosaicQuantResults[]>(2, 10, proj.getLog());
					// hive.addCallable(worker);
					// hive.execute(true);
					// ArrayList<MosaicQuantResults[]> results = hive.getResults();

					for (int i = 0; i < sampCNVs.getLoci().length; i++) {
						System.out.println(sampCNVs.getLoci()[i].toPlinkFormat());
						LocusSet<MosaicRegion> mosSet = md.callMosaic(sampCNVs.getLoci()[i], true);

						// double f = results.get(0)[0].getFs()[i];
						// int numFMarkers = results.get(0)[0].getNumMarkers()[i];
						double f = Double.NaN;
						int numFMarkers = -1;
						if (mosSet.getLoci().length > 0) {

							for (int j = 0; j < mosSet.getLoci().length; j++) {

								MosaicRegion recaptured = new MosaicRegion(	sampCNVs.getLoci()[i],
																														mosSet.getLoci()[j]);
								recaptured.setF(f);
								recaptured.setNumCustomFMarkers(mosSet.getLoci()[j].getNumMarkers());
								recaptured.setNumFMarkers(numFMarkers);
								all.add(recaptured);
							}
						} else {
							// public MosaicRegion(CNVariant cnv, double bpWeightedScore, double
							// nearestStateScore, double pdfScore, double delta, double f, double customF) {

							MosaicRegion blank = new MosaicRegion(sampCNVs.getLoci()[i], 0, -1, 0, 0, f, 0);
							blank.setNumFMarkers(numFMarkers);
							all.add(blank);
						}
						System.out.println("DFSDF\t" + all.get(all.size() - 1).toPlinkFormat());

					}
					if (all.size() != sampCNVs.getLoci().length) {
						throw new IllegalStateException("Did not call all cnvs for sample " + sample);
					}
					MosaicRegion[] finals = all.toArray(new MosaicRegion[all.size()]);

					BeastScore beastScore = BeastScore.beastInd(proj, sampleData, samp.getLRRs(), finals,
																											markerSet.getChrs(), markerSet.getPositions(),
																											indicesByChr);
					beastScore.computeBeastScores();
					for (int i = 0; i < beastScore.getBeastScores().length; i++) {
						finals[i].setBeastScore(beastScore.getBeastScores()[i]);
						finals[i].setBeastHeight(beastScore.getBeastHeights()[i]);
						finals[i].setBeastLength(beastScore.getBeastLengths()[i]);
					}

					return finals;
				}
			};
			index++;
			return callable;
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;

		String usage = "\n"	+ "one.JL.CNVMosaic requires 0-1 arguments\n"
										+ "   (1) project  (i.e. proj=" + filename + " (default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("proj=")) {
				filename = arg.split("=")[1];
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
			Project proj = new Project(filename, false);
			computeCNVMosaic(proj);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
