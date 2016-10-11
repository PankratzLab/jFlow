package org.genvisis.one.JL;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.ExtProjectDataParser;
import org.genvisis.cnv.manage.ExtProjectDataParser.ProjectDataParserBuilder;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.stats.Correlation;
import org.genvisis.stats.Rscript.COLUMNS_MULTIPLOT;
import org.genvisis.stats.Rscript.PLOT_DEVICE;
import org.genvisis.stats.Rscript.RScatter;
import org.genvisis.stats.Rscript.RScatters;
import org.genvisis.stats.Rscript.SCATTER_TYPE;

import com.google.common.primitives.Ints;

/**
 * Look at correlation of mito estimates from different pc files, relies on typed/*.gz files
 *
 */
public class CorrelTypePcs {

	private static final String SUB_DIR = "_eval/typed/";

	private static void run(Project proj, String[] fullPathPcFile) throws FileNotFoundException {
		Logger log = new Logger();
		String outDir = proj.PROJECT_DIRECTORY.getValue()	+ "comparePCMethods/"
										+ ext.rootOf(fullPathPcFile[0]) + "/";
		new File(outDir).mkdirs();
		if (fullPathPcFile.length != 2) {
			proj.getLog().reportTimeError("Must be two pc files only");
			return;
		}
		String out = outDir + ext.rootOf(fullPathPcFile[0]) + "_pcs.comp.txt";
		ArrayList<RScatter> rsList = new ArrayList<RScatter>();

		if (!Files.exists(out)) {
			try {
				PrintWriter writer = new PrintWriter(new FileWriter(out));
				ArrayList<String> outHeader = new ArrayList<String>();
				outHeader.add("PC_FILE_1");
				outHeader.add("PC_FILE_2");
				outHeader.add("METHOD_FILE");
				outHeader.add("HAS");
				outHeader.add("PC");
				outHeader.add("SPEARCorrelation");
				outHeader.add("PEARCorrelation");

				writer.println(Array.toStr(Array.toStringArray(outHeader)));

				for (int pcFile1Index = 0; pcFile1Index < fullPathPcFile.length; pcFile1Index++) {
					String getItDir = ext.rootOf(fullPathPcFile[pcFile1Index], false) + SUB_DIR;
					String[] gzippers = Files.list(getItDir, ".gz", false);

					for (String gzipper : gzippers) {
						ProjectDataParserBuilder builderCurrent = new ProjectDataParserBuilder();

						builderCurrent.numericDataTitles(getPCIndices(getItDir + gzipper, log));
						builderCurrent.sampleBased(true);
						builderCurrent.dataKeyColumnName("DNA");
						builderCurrent.treatAllNumeric(false);

						ExtProjectDataParser parser = builderCurrent.build(proj, getItDir + gzipper);
						parser.determineIndicesFromTitles();
						parser.loadData();
						for (int pcFile2Index = pcFile1Index
																		+ 1; pcFile2Index < fullPathPcFile.length; pcFile2Index++) {
							String comp = ext.rootOf(fullPathPcFile[pcFile2Index], false) + SUB_DIR + gzipper;
							if (Files.exists(comp)) {
								ProjectDataParserBuilder builderComp = new ProjectDataParserBuilder();

								builderComp.numericDataTitles(getPCIndices(comp, log));
								builderComp.sampleBased(true);
								builderComp.dataKeyColumnName("DNA");
								builderComp.treatAllNumeric(false);
								System.out.println(Files.getHeaderOfFile(comp, log)[0]);
								ExtProjectDataParser parserComp = builderComp.build(proj, comp);
								parserComp.determineIndicesFromTitles();
								parserComp.loadData();
								int min = Math.min(	parser.getNumericDataTitles().length,
																		parserComp.getNumericDataTitles().length);
								for (int pcIndex = 0; pcIndex < min; pcIndex++) {
									ArrayList<String> outData = new ArrayList<String>();
									outData.add(fullPathPcFile[pcFile1Index]);
									outData.add(fullPathPcFile[pcFile2Index]);
									outData.add(ext.rootOf(gzipper));
									outData.add("TRUE");
									if (parser.getNumericDataTitles()[pcIndex].equals(parserComp.getNumericDataTitles()[pcIndex])) {
										outData.add(parser.getNumericDataTitles()[pcIndex].replaceAll("PC", ""));
										String pc = parser.getNumericDataTitles()[pcIndex];
										double[][] data = new double[][] {parser.getNumericDataForTitle(pc),
																											parserComp.getNumericDataForTitle(pc)};
										ArrayList<Integer> indicesToCorrel = new ArrayList<Integer>();
										for (int i = 0; i < data[0].length; i++) {
											if (!Double.isNaN(data[0][i]) && !Double.isNaN(data[1][i])) {
												indicesToCorrel.add(i);
											}
										}
										// proj.getLog().reportTimeInfo(indicesToCorrel.size() + " valid estimates of "
										// + data[0].length + " total for " + gzippers[gzipIndex]);
										data[0] = Array.subArray(data[0], Ints.toArray(indicesToCorrel));
										data[1] = Array.subArray(data[1], Ints.toArray(indicesToCorrel));

										outData.add(Correlation.Spearman(data)[0] + "");
										outData.add(Correlation.Pearson(data)[0] + "");
										writer.println(Array.toStr(Array.toStringArray(outData)));
									} else {
										writer.close();

										throw new IllegalStateException("Mismatched PC order");
									}
								}
							}
						}
					}

				}

				writer.close();
			} catch (Exception e) {
				log.reportError("Error writing to " + out);
				log.reportException(e);
			}
		}
		String tmp1 = out + "_basic";
		RScatter rsScatter = new RScatter(out, tmp1 + ".rscript", ext.removeDirectoryInfo(tmp1),
																			tmp1 + ".jpeg", "PC", new String[] {"SPEARCorrelation"},
																			"METHOD_FILE", SCATTER_TYPE.POINT, log);
		rsScatter.setTitle(ext.removeDirectoryInfo(fullPathPcFile[0])	+ " vs "
												+ ext.removeDirectoryInfo(fullPathPcFile[1]));
		rsScatter.setOverWriteExisting(true);
		rsScatter.setxLabel("PC");
		rsScatter.setyLabel("SPEARMAN r");
		rsScatter.setyRange(new double[] {-.5, 1});
		rsList.add(rsScatter);

		String tmp2 = out + "_basicPear";
		RScatter rsScatterP = new RScatter(	out, tmp2 + ".rscript", ext.removeDirectoryInfo(tmp2),
																				tmp2 + ".jpeg", "PC", new String[] {"PEARCorrelation"},
																				"METHOD_FILE", SCATTER_TYPE.POINT, log);
		rsScatterP.setTitle(ext.removeDirectoryInfo(fullPathPcFile[0])	+ " vs "
												+ ext.removeDirectoryInfo(fullPathPcFile[1]));
		rsScatterP.setOverWriteExisting(true);
		rsScatterP.setxLabel("PC");
		rsScatterP.setyLabel("PEARSON r");
		rsScatterP.setyRange(new double[] {-.5, 1});
		rsList.add(rsScatterP);

		RScatters rsScatters = new RScatters(	rsList.toArray(new RScatter[rsList.size()]),
																					out + ".rscript", out + ".pdf",
																					COLUMNS_MULTIPLOT.COLUMNS_MULTIPLOT_1, PLOT_DEVICE.PDF,
																					log);
		rsScatters.execute();
		// new RScatter(out, out+".rscript", ext.removeDirectoryInfo(out), out+".jpeg", "PC", new
		// String, sType, log)

		// for (int i = 0; i < fullPathPcFile.length; i++) {
		// String getItDir = ext.rootOf(fullPathPcFile[i], false) + SUB_DIR;
		// String[] gzippers = Files.list(getItDir, ".gz", false);
		//
		// for (int j = 0; j < gzippers.length; j++) {
		//
		// for (int j2 = i + 1; j2 < fullPathPcFile.length; j2++) {
		// String comp = ext.rootOf(fullPathPcFile[j2], false) + SUB_DIR + gzippers[j];
		// if (Files.exists(comp)) {
		// Restrictions restrictions = new Restrictions(new String[]{"METHOD_FILE"}, vals, ops, cond)
		// }
		// }
		// }
		// }

	}

	private static String[] getPCIndices(String file, Logger log) {
		String[] header = Files.getHeaderOfFile(file, log);
		ArrayList<Integer> ind = new ArrayList<Integer>();
		for (int i = 0; i < header.length; i++) {
			if (header[i].startsWith("PC")) {
				ind.add(i);
			}
		}
		return Array.subArray(header, Ints.toArray(ind));

	}

	public static void main(String[] args) {

		ArrayList<String> pcFilesAll = new ArrayList<String>();
		Project proj = new Project("/home/pankrat2/lanej/projects/aric_exome.properties", false);
		pcFilesAll.add("/home/pankrat2/shared/aric_exome_chip/aric_exomeALL_1000PCs_OHW_40_ws15_recomp_gc_corrected.PCs.extrapolated.txt");
		pcFilesAll.add("/home/pankrat2/shared/aric_exome_chip/gc_corrected/aric_exomeALL_1000PCs_OHW_40_ws15_gc_corrected.PCs.extrapolated.txt");
		ArrayList<String> pcFilesW = new ArrayList<String>();
		pcFilesW.add("/home/pankrat2/shared/aric_exome_chip/aric_exomeW_1000PCs_OHW_40_ws15_recomp_gc_corrected.PCs.extrapolated.txt");
		pcFilesW.add("/home/pankrat2/shared/aric_exome_chip/gc_corrected/aric_exomeW_1000PCs_OHW_40_ws15_gc_corrected.PCs.extrapolated.txt");

		try {
			run(proj, Array.toStringArray(pcFilesAll));
			run(proj, Array.toStringArray(pcFilesW));

		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

}
