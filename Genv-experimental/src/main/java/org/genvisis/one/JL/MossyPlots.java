package org.genvisis.one.JL;

import java.io.File;
import java.util.ArrayList;
import java.util.Set;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.ext;
import org.genvisis.seq.manage.VCFOps.VcfPopulation;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import org.genvisis.stats.Rscript;
import org.genvisis.stats.Rscript.COLUMNS_MULTIPLOT;
import org.genvisis.stats.Rscript.GeomText;
import org.genvisis.stats.Rscript.LEGEND_POSITION;
import org.genvisis.stats.Rscript.PLOT_DEVICE;
import org.genvisis.stats.Rscript.RScatter;
import org.genvisis.stats.Rscript.RScatters;
import org.genvisis.stats.Rscript.Restrictions;
import org.genvisis.stats.Rscript.SCATTER_TYPE;

public class MossyPlots {

	private static void plot(Project proj, String mosBaseFile) {
		Set<String> tns = VcfPopulation
																		.load(proj.PROJECT_DIRECTORY.getValue()	+ "TN.vpop",
																					POPULATION_TYPE.TUMOR_NORMAL, proj.getLog())
																		.getTumorSamples();
		String outDir = proj.PROJECT_DIRECTORY.getValue() + ext.rootOf(mosBaseFile) + "mosPlots/";
		new File(outDir).mkdirs();
		proj.getLog().reportTimeInfo(mosBaseFile);
		SampleData sampleData = proj.getSampleData(0, false);
		String[][] mos = HashVec.loadFileToStringMatrix(mosBaseFile, false, null, false);
		String mosFile = outDir + "mosResults_excludesTagged.txt";

		ArrayList<String> filtered = new ArrayList<String>();
		filtered.add(Array.toStr(mos[0]) + "\tExclude\tTUMOR\tChr");
		ArrayList<GeomText> geomTexts = new ArrayList<Rscript.GeomText>();
		// public static final String[] MOSAICISM_HEADER = { "Sample", "Band", "LRR N", "mean LRR", "BAF
		// N", "SD of BAF (0.15-0.85)", "IQR of BAF (0.15-0.85)", "%Homo", "BandPercentMosaicism",
		// "BpWeightedAverage", "NumberRegionsDetected" };

		for (int i = 1; i < mos.length; i++) {
			try {
				Double.parseDouble(mos[i][5]);
			} catch (Exception e) {
				mos[i][5] = "NaN";
			}
			try {
				Double.parseDouble(mos[i][6]);
			} catch (Exception e) {
				mos[i][6] = "NaN";
			}
			if (Double.parseDouble(mos[i][9]) >= 0) {
				String data = Array.toStr(mos[i])	+ "\t"
											+ (sampleData.individualShouldBeExcluded(mos[i][0]) ? 1 : 0) + "\t"
											+ (tns.contains(mos[i][0]) ? "TUMOR" : "GERMLINE") + "\t"
											+ mos[i][1].replaceAll("p", "").replaceAll("q", "");
				if (Double.parseDouble(mos[i][9]) >= .25) {
					geomTexts.add(new GeomText(	Double.parseDouble(mos[i][3]), Double.parseDouble(mos[i][9]),
																			.25, mos[i][1], .1));
				}
				filtered.add(data);
			}
		}
		Files.writeArray(Array.toStringArray(filtered), mosFile);
		ArrayList<RScatter> rscScatters = new ArrayList<RScatter>();

		String tNAll = outDir + "tumorNormalAll.plot";

		RScatter rsScatterAll = new RScatter(	mosFile, tNAll + ".rscript", ext.rootOf(tNAll),
																					tNAll + ".jpeg", "mean LRR",
																					new String[] {"BpWeightedAverage"}, "TUMOR",
																					SCATTER_TYPE.BASIC_POINT, proj.getLog());
		rsScatterAll.setOverWriteExisting(true);
		rsScatterAll.setyRange(new double[] {0, 1});
		rsScatterAll.setxRange(new double[] {-.5, .5});
		rsScatterAll.setxLabel("Average LRR of chromosomal arm");
		rsScatterAll.setyLabel("Mosaic Metric of chromosomal arm");
		rsScatterAll.setTitle("All Samples");
		rsScatterAll.setFontsize(14);

		// Restrictions[] restrictionss = new Restrictions[] { new Restrictions(new String[] { "Exclude"
		// }, new String[] { "true" }, new String[] { "==" }, null) };
		// rsScatter.setRestrictions(restrictionss);
		rsScatterAll.execute();
		rscScatters.add(rsScatterAll);

		String tumorNormalNoExclude = outDir + "tumorNormal.plot";

		RScatter rsScatterTN = new RScatter(mosFile, tumorNormalNoExclude + ".rscript",
																				ext.rootOf(tumorNormalNoExclude),
																				tumorNormalNoExclude + ".jpeg", "mean LRR",
																				new String[] {"BpWeightedAverage"}, "TUMOR",
																				SCATTER_TYPE.BASIC_POINT, proj.getLog());
		rsScatterTN.setOverWriteExisting(true);
		rsScatterTN.setyRange(new double[] {0, 1});
		rsScatterTN.setxRange(new double[] {-.5, .5});
		Restrictions[] restrictionTN = new Restrictions[] {new Restrictions(new String[] {"Exclude"},
																																				new double[] {0},
																																				new String[] {"=="}, null)};
		rsScatterTN.setRestrictions(restrictionTN);
		rsScatterTN.setxLabel("Average LRR of chromosomal arm");
		rsScatterTN.setyLabel("Mosaic Metric of chromosomal arm");
		rsScatterTN.setTitle("Non-excluded Samples");
		rsScatterTN.execute();
		rsScatterTN.setFontsize(14);

		rscScatters.add(rsScatterTN);

		String tumorNormalNoExcludeSamp = outDir + "tumorNormalSamp.plot";

		RScatter rsScatterTNSamp = new RScatter(mosFile, tumorNormalNoExcludeSamp + ".rscript",
																						ext.rootOf(tumorNormalNoExcludeSamp),
																						tumorNormalNoExcludeSamp + ".jpeg", "mean LRR",
																						new String[] {"BpWeightedAverage"}, "Sample",
																						SCATTER_TYPE.BASIC_POINT, proj.getLog());
		rsScatterTNSamp.setOverWriteExisting(true);

		rsScatterTNSamp.setyRange(new double[] {0, 1});
		rsScatterTNSamp.setxRange(new double[] {-.5, .5});
		Restrictions[] restrictionTNSamp = new Restrictions[] {new Restrictions(new String[] {"Exclude"},
																																						new double[] {0},
																																						new String[] {"=="},
																																						null)};
		rsScatterTNSamp.setlPosition(LEGEND_POSITION.NONE);
		rsScatterTNSamp.setRestrictions(restrictionTNSamp);
		rsScatterTNSamp.setxLabel("Average LRR of chromosomal arm");
		rsScatterTNSamp.setyLabel("Mosaic Metric of chromosomal arm");
		rsScatterTNSamp.setTitle("Non-excluded Samples");
		rsScatterTNSamp.execute();
		rsScatterTNSamp.setFontsize(14);

		rscScatters.add(rsScatterTNSamp);

		String tumorNormalNoExcludeFocus = outDir + "tumorNormalFocus.plot";
		double filter = 0.10;
		RScatter rsScatterTNFocus = new RScatter(	mosFile, tumorNormalNoExcludeFocus + ".rscript",
																							ext.rootOf(tumorNormalNoExcludeFocus),
																							tumorNormalNoExcludeFocus + ".jpeg", "mean LRR",
																							new String[] {"BpWeightedAverage"}, "TUMOR",
																							SCATTER_TYPE.BASIC_POINT, proj.getLog());
		rsScatterTNFocus.setOverWriteExisting(true);
		rsScatterTNFocus.setyRange(new double[] {0, 1});
		rsScatterTNFocus.setxRange(new double[] {-.5, .5});
		rsScatterTNFocus.setTitle("Samples with Mosaic Metric >=" + filter);
		rsScatterTNFocus.setFontsize(14);
		// rsScatterTNFocus.setgTexts(geomTexts.toArray(new GeomText[geomTexts.size()]));
		// rsScatterTNFocus.setHorizontalLines(new HorizontalLine[] { horizontalLine });
		Restrictions mm = new Restrictions(	new String[] {"BpWeightedAverage"}, new double[] {filter},
																				new String[] {">="}, null);
		Restrictions[] restrictionTNFocus = new Restrictions[] {mm,
																														new Restrictions(	new String[] {"Exclude"},
																																							new double[] {0},
																																							new String[] {"=="},
																																							null)};
		rsScatterTNFocus.setRestrictions(restrictionTNFocus);
		rsScatterTNFocus.setxLabel("Average LRR of chromosomal arm");
		rsScatterTNFocus.setyLabel("Mosaic Metric of chromosomal arm");

		rsScatterTNFocus.execute();
		rscScatters.add(rsScatterTNFocus);

		String tumorNormalNoExcludeFocus2 = outDir + "tumorNormalFocus2.plot";
		double filter2 = 0.25;
		RScatter rsScatterTNFocus2 = new RScatter(mosFile, tumorNormalNoExcludeFocus2 + ".rscript",
																							ext.rootOf(tumorNormalNoExcludeFocus2),
																							tumorNormalNoExcludeFocus2 + ".jpeg", "mean LRR",
																							new String[] {"BpWeightedAverage"}, "Sample",
																							SCATTER_TYPE.BASIC_POINT, proj.getLog());
		rsScatterTNFocus2.setOverWriteExisting(true);
		rsScatterTNFocus2.setyRange(new double[] {0, 1});
		rsScatterTNFocus2.setxRange(new double[] {-.5, .5});
		rsScatterTNFocus2.setTitle("Samples with Mosaic Metric >=" + filter2);
		rsScatterTNFocus2.setFontsize(14);
		// rsScatterTNFocus.setgTexts(geomTexts.toArray(new GeomText[geomTexts.size()]));
		// rsScatterTNFocus.setHorizontalLines(new HorizontalLine[] { horizontalLine });
		Restrictions mm2 = new Restrictions(new String[] {"BpWeightedAverage"}, new double[] {filter2},
																				new String[] {">="}, null);
		Restrictions[] restrictionTNFocus2 = new Restrictions[] {	mm2,
																															new Restrictions(	new String[] {"Exclude"},
																																								new double[] {0},
																																								new String[] {"=="},
																																								null)};
		rsScatterTNFocus2.setRestrictions(restrictionTNFocus2);
		rsScatterTNFocus2.setxLabel("Average LRR of chromosomal arm");
		rsScatterTNFocus2.setyLabel("Mosaic Metric of chromosomal arm");

		rsScatterTNFocus2.execute();
		rscScatters.add(rsScatterTNFocus2);

		String tumorNormalNoExcludeCompallColors = outDir + "tumorNormalCompMethodallColors.plot";

		RScatter rsScatterTNBandCompallColors = new RScatter(	mosFile,
																													tumorNormalNoExcludeCompallColors
																																		+ ".rscript",
																													ext.rootOf(tumorNormalNoExcludeCompallColors),
																													tumorNormalNoExcludeCompallColors + ".jpeg",
																													"SD of BAF (0.15-0.85)",
																													new String[] {"IQR of BAF (0.15-0.85)"},
																													"TUMOR", SCATTER_TYPE.BASIC_POINT,
																													proj.getLog());
		rsScatterTNBandCompallColors.setOverWriteExisting(true);
		rsScatterTNBandCompallColors.setFactorColor(false);
		rsScatterTNBandCompallColors.setyRange(new double[] {0, .6});
		rsScatterTNBandCompallColors.setxRange(new double[] {0, .3});
		// Restrictions[] restrictionTNCompallColors = new Restrictions[] { new Restrictions(new
		// String[] { "Exclude" }, new double[] { 0 }, new String[] { "==" }, null) };
		// rsScatterTNBandCompallColors.setRestrictions(restrictionTNCompallColors);
		// rsScatterTNBandComp.setlPosition(LEGEND_POSITION.NONE);
		rsScatterTNBandCompallColors.setxLabel("SD of BAF (0.15-0.85)");
		rsScatterTNBandCompallColors.setyLabel("IQR of BAF (0.15-0.85)");
		rsScatterTNBandCompallColors.setTitle("All Samples");
		rsScatterTNBandCompallColors.execute();
		rsScatterTNBandCompallColors.setFontsize(14);

		rscScatters.add(rsScatterTNBandCompallColors);

		String tumorNormalNoExcludeCompnoEColors = outDir + "tumorNormalCompMethodnEColors.plot";

		RScatter rsScatterTNBandCompallNEColors = new RScatter(	mosFile,
																														tumorNormalNoExcludeCompnoEColors
																																			+ ".rscript",
																														ext.rootOf(tumorNormalNoExcludeCompnoEColors),
																														tumorNormalNoExcludeCompnoEColors + ".jpeg",
																														"SD of BAF (0.15-0.85)",
																														new String[] {"IQR of BAF (0.15-0.85)"},
																														"TUMOR", SCATTER_TYPE.BASIC_POINT,
																														proj.getLog());
		rsScatterTNBandCompallNEColors.setOverWriteExisting(true);
		rsScatterTNBandCompallNEColors.setFactorColor(false);
		rsScatterTNBandCompallNEColors.setyRange(new double[] {0, .6});
		rsScatterTNBandCompallNEColors.setxRange(new double[] {0, .3});
		Restrictions[] restrictionTNCompallNEColors =
																								new Restrictions[] {new Restrictions(	new String[] {"Exclude"},
																																											new double[] {0},
																																											new String[] {"=="},
																																											null)};
		rsScatterTNBandCompallNEColors.setRestrictions(restrictionTNCompallNEColors);
		// rsScatterTNBandComp.setlPosition(LEGEND_POSITION.NONE);
		rsScatterTNBandCompallNEColors.setxLabel("SD of BAF (0.15-0.85)");
		rsScatterTNBandCompallNEColors.setyLabel("IQR of BAF (0.15-0.85)");
		rsScatterTNBandCompallNEColors.setTitle("Non-excluded Samples");
		rsScatterTNBandCompallNEColors.execute();
		rsScatterTNBandCompallNEColors.setFontsize(14);

		rscScatters.add(rsScatterTNBandCompallNEColors);


		String tumorNormalNoExcludeCompall = outDir + "tumorNormalCompMethodall.plot";

		RScatter rsScatterTNBandCompall =
																		new RScatter(	mosFile, tumorNormalNoExcludeCompall + ".rscript",
																									ext.rootOf(tumorNormalNoExcludeCompall),
																									tumorNormalNoExcludeCompall + ".jpeg",
																									"SD of BAF (0.15-0.85)",
																									new String[] {"IQR of BAF (0.15-0.85)"},
																									"BpWeightedAverage", SCATTER_TYPE.BASIC_POINT,
																									proj.getLog());
		rsScatterTNBandCompall.setOverWriteExisting(true);
		rsScatterTNBandCompall.setFactorColor(false);
		rsScatterTNBandCompall.setyRange(new double[] {0, .6});
		rsScatterTNBandCompall.setxRange(new double[] {0, .3});
		Restrictions[] restrictionTNCompall =
																				new Restrictions[] {new Restrictions(	new String[] {"Exclude"},
																																							new double[] {0},
																																							new String[] {"=="},
																																							null)};
		rsScatterTNBandCompall.setRestrictions(restrictionTNCompall);
		// rsScatterTNBandComp.setlPosition(LEGEND_POSITION.NONE);
		rsScatterTNBandCompall.setxLabel("SD of BAF (0.15-0.85)");
		rsScatterTNBandCompall.setyLabel("IQR of BAF (0.15-0.85)");
		rsScatterTNBandCompall.setTitle("Non-excluded Samples");
		rsScatterTNBandCompall.execute();
		rsScatterTNBandCompall.setFontsize(14);

		rscScatters.add(rsScatterTNBandCompall);

		String tumorNormalNoExcludeCompone = outDir + "tumorNormalCompMethod1.plot";

		RScatter rsScatterTNBandCompone =
																		new RScatter(	mosFile, tumorNormalNoExcludeCompone + ".rscript",
																									ext.rootOf(tumorNormalNoExcludeCompone),
																									tumorNormalNoExcludeCompone + ".jpeg",
																									"SD of BAF (0.15-0.85)",
																									new String[] {"IQR of BAF (0.15-0.85)"},
																									"BpWeightedAverage", SCATTER_TYPE.BASIC_POINT,
																									proj.getLog());
		rsScatterTNBandCompone.setOverWriteExisting(true);
		rsScatterTNBandCompone.setFactorColor(false);
		rsScatterTNBandCompone.setyRange(new double[] {0, .6});
		rsScatterTNBandCompone.setxRange(new double[] {0, .3});
		Restrictions[] restrictionTNCompone = new Restrictions[] {mm,
																															new Restrictions(	new String[] {"Exclude"},
																																								new double[] {0},
																																								new String[] {"=="},
																																								null)};
		rsScatterTNBandCompone.setRestrictions(restrictionTNCompone);
		// rsScatterTNBandComp.setlPosition(LEGEND_POSITION.NONE);
		rsScatterTNBandCompone.setxLabel("SD of BAF (0.15-0.85)");
		rsScatterTNBandCompone.setyLabel("IQR of BAF (0.15-0.85)");
		rsScatterTNBandCompone.setTitle("Non-excluded Samples with Mosaic Metric >=" + filter);
		rsScatterTNBandCompone.execute();
		rsScatterTNBandCompone.setFontsize(14);

		rscScatters.add(rsScatterTNBandCompone);

		String tumorNormalNoExcludeComp = outDir + "tumorNormalCompMethod.plot";

		RScatter rsScatterTNBandComp = new RScatter(mosFile, tumorNormalNoExcludeComp + ".rscript",
																								ext.rootOf(tumorNormalNoExcludeComp),
																								tumorNormalNoExcludeComp + ".jpeg",
																								"SD of BAF (0.15-0.85)",
																								new String[] {"IQR of BAF (0.15-0.85)"},
																								"BpWeightedAverage", SCATTER_TYPE.BASIC_POINT,
																								proj.getLog());
		rsScatterTNBandComp.setOverWriteExisting(true);
		rsScatterTNBandComp.setFactorColor(false);
		rsScatterTNBandComp.setyRange(new double[] {0, .6});
		rsScatterTNBandComp.setxRange(new double[] {0, .3});
		Restrictions[] restrictionTNComp = new Restrictions[] {	mm2,
																														new Restrictions(	new String[] {"Exclude"},
																																							new double[] {0},
																																							new String[] {"=="},
																																							null)};
		rsScatterTNBandComp.setRestrictions(restrictionTNComp);
		// rsScatterTNBandComp.setlPosition(LEGEND_POSITION.NONE);
		rsScatterTNBandComp.setxLabel("SD of BAF (0.15-0.85)");
		rsScatterTNBandComp.setyLabel("IQR of BAF (0.15-0.85)");
		rsScatterTNBandComp.setTitle("Non-excluded Samples with Mosaic Metric >=" + filter2);
		rsScatterTNBandComp.execute();
		rsScatterTNBandComp.setFontsize(14);

		rscScatters.add(rsScatterTNBandComp);

		String tumorNormalNoExcludeFocus2arm = outDir + "tumorNormalFocus2arm.plot";
		RScatter rsScatterTNFocus2arm =
																	new RScatter(	mosFile, tumorNormalNoExcludeFocus2arm + ".rscript",
																								ext.rootOf(tumorNormalNoExcludeFocus2arm),
																								tumorNormalNoExcludeFocus2arm + ".jpeg", "mean LRR",
																								new String[] {"BpWeightedAverage"}, "Band",
																								SCATTER_TYPE.BASIC_POINT, proj.getLog());
		rsScatterTNFocus2arm.setOverWriteExisting(true);
		rsScatterTNFocus2arm.setyRange(new double[] {0, 1});
		rsScatterTNFocus2arm.setxRange(new double[] {-.5, .5});
		rsScatterTNFocus2arm.setTitle("Samples with Mosaic Metric >="	+ filter2
																	+ " (Colored by chromosomal Arm)");
		rsScatterTNFocus2arm.setFontsize(14);
		rsScatterTNFocus2arm.setlPosition(LEGEND_POSITION.NONE);

		// rsScatterTNFocus.setgTexts(geomTexts.toArray(new GeomText[geomTexts.size()]));
		// rsScatterTNFocus.setHorizontalLines(new HorizontalLine[] { horizontalLine });
		Restrictions mm2arm = new Restrictions(	new String[] {"BpWeightedAverage"},
																						new double[] {filter2}, new String[] {">="}, null);
		Restrictions[] restrictionTNFocus2arm = new Restrictions[] {mm2arm,
																																new Restrictions(	new String[] {"Exclude"},
																																									new double[] {0},
																																									new String[] {"=="},
																																									null)};
		rsScatterTNFocus2arm.setRestrictions(restrictionTNFocus2arm);
		rsScatterTNFocus2arm.setxLabel("Average LRR of chromosomal arm");
		rsScatterTNFocus2arm.setyLabel("Mosaic Metric of chromosomal arm");

		rsScatterTNFocus2arm.execute();
		rscScatters.add(rsScatterTNFocus2arm);

		String tumorNormalNoExcludeFocus2chr = outDir + "tumorNormalFocus2chr.plot";
		RScatter rsScatterTNFocus2chr =
																	new RScatter(	mosFile, tumorNormalNoExcludeFocus2chr + ".rscript",
																								ext.rootOf(tumorNormalNoExcludeFocus2chr),
																								tumorNormalNoExcludeFocus2chr + ".jpeg", "mean LRR",
																								new String[] {"BpWeightedAverage"}, "Chr",
																								SCATTER_TYPE.BASIC_POINT, proj.getLog());
		rsScatterTNFocus2chr.setOverWriteExisting(true);
		rsScatterTNFocus2chr.setyRange(new double[] {0, 1});
		rsScatterTNFocus2chr.setxRange(new double[] {-.5, .5});
		rsScatterTNFocus2chr.setTitle("Samples with Mosaic Metric >="	+ filter2
																	+ " (Colored by chromosome)");
		rsScatterTNFocus2chr.setFontsize(14);
		// rsScatterTNFocus2chr.setlPosition(LEGEND_POSITION.NONE);
		// rsScatterTNFocus.setgTexts(geomTexts.toArray(new GeomText[geomTexts.size()]));
		// rsScatterTNFocus.setHorizontalLines(new HorizontalLine[] { horizontalLine });
		Restrictions mm2chr = new Restrictions(	new String[] {"BpWeightedAverage"},
																						new double[] {filter2}, new String[] {">="}, null);
		Restrictions[] restrictionTNFocus2chr = new Restrictions[] {mm2chr,
																																new Restrictions(	new String[] {"Exclude"},
																																									new double[] {0},
																																									new String[] {"=="},
																																									null)};
		rsScatterTNFocus2chr.setRestrictions(restrictionTNFocus2chr);
		rsScatterTNFocus2chr.setxLabel("Average LRR of chromosomal arm");
		rsScatterTNFocus2chr.setyLabel("Mosaic Metric of chromosomal arm");

		rsScatterTNFocus2chr.execute();
		rscScatters.add(rsScatterTNFocus2chr);

		String tumorNormalNoExcludeBand = outDir + "tumorNormalBand.plot";

		RScatter rsScatterTNBand = new RScatter(mosFile, tumorNormalNoExcludeBand + ".rscript",
																						ext.rootOf(tumorNormalNoExcludeBand),
																						tumorNormalNoExcludeBand + ".jpeg", "mean LRR",
																						new String[] {"BpWeightedAverage"}, "Band",
																						SCATTER_TYPE.BASIC_POINT, proj.getLog());
		rsScatterTNBand.setOverWriteExisting(true);
		rsScatterTNBand.setyRange(new double[] {0, 1});
		rsScatterTNBand.setxRange(new double[] {-.5, .5});
		Restrictions[] restrictionTNBand = new Restrictions[] {new Restrictions(new String[] {"Exclude"},
																																						new double[] {0},
																																						new String[] {"=="},
																																						null)};
		rsScatterTNBand.setRestrictions(restrictionTNBand);
		rsScatterTNBand.setlPosition(LEGEND_POSITION.NONE);
		rsScatterTNBand.setxLabel("Average LRR of chromosomal arm");
		rsScatterTNBand.setyLabel("Mosaic Metric of chromosomal arm");
		rsScatterTNBand.setTitle("Non-excluded Samples (Colored by chromosomal Arm)");
		rsScatterTNBand.execute();
		rsScatterTNBand.setFontsize(14);

		rscScatters.add(rsScatterTNBand);

		String finalRscript = outDir + "final.rscript";
		String finalout = outDir + "final.pdf";

		RScatters rScatters =
												new RScatters(rscScatters.toArray(new RScatter[rscScatters.size()]),
																			finalRscript, finalout, COLUMNS_MULTIPLOT.COLUMNS_MULTIPLOT_1,
																			PLOT_DEVICE.PDF, proj.getLog());
		rScatters.execute();
	}



	public static void main(String[] args) {
		Project proj = new Project("C:/workspace/Genvisis/projects/Cushings.properties", false);
		// generateGCModel(proj, new ReferenceGenome(proj.getReferenceGenomeFASTAFilename(),
		// proj.getLog()), 1000000);


		String[] mosResults = Files.list(	proj.RESULTS_DIRECTORY.getValue(), "Mosaicism", ".xln", true,
																			false, true);
		for (String mosResult : mosResults) {
			plot(proj, mosResult);
		}

	}

}
