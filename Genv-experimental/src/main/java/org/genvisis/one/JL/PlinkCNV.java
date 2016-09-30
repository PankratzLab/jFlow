package org.genvisis.one.JL;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.concurrent.Callable;

import org.genvisis.cnv.analysis.FilterCalls;
import org.genvisis.cnv.manage.Resources;
import org.genvisis.cnv.manage.Resources.GENOME_BUILD;
import org.genvisis.common.Array;
import org.genvisis.common.CmdLine;
import org.genvisis.common.ExcelConverter;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.Matrix;
import org.genvisis.common.WorkerHive;
import org.genvisis.common.ext;
import org.genvisis.filesys.CNVariant;
import org.genvisis.filesys.GeneData;
import org.genvisis.filesys.GeneTrack;
import org.genvisis.filesys.LocusSet;
import org.genvisis.filesys.Segment;
import org.genvisis.stats.LeastSquares;

public class PlinkCNV {

	private static void run(String cnvaFile, String sampFile, final String build) {

		final String dir = ext.parseDirectoryOfFile(cnvaFile);
		final Logger log = new Logger(dir + "cnv.log");
		String merge = dir + "penncnv.cnv";
		if (!Files.exists(merge)) {
			log.reportTimeInfo("Merging calls");
			FilterCalls.mergeCNVs(cnvaFile, merge, FilterCalls.DEFAULT_CLEAN_FACTOR, null);
		}
		String[] filters = Files.list(dir, ".crf", false);

		WorkerHive<PlinkResult> hive = new WorkerHive<PlinkResult>(6, 10, log);
		GeneTrack geneTrack = GeneTrack.load(Resources.genome(GENOME_BUILD.valueOf(build), log).getGTrack().get(), false);
		if (Files.exists(dir + "mitoCarta.txt")) {
			String[] mitos = HashVec.loadFileToStringArray(dir + "mitoCarta.txt", false, new int[] { 0 }, true);
			ArrayList<String> mitoList = new ArrayList<String>();
			int count = 0;
			for (int j = 0; j < mitos.length; j++) {
				if (geneTrack.lookupAllGeneData(mitos[j]) != null && geneTrack.lookupAllGeneData(mitos[j]).length > 0) {
					count++;
					GeneData[] gds = geneTrack.lookupAllGeneData(mitos[j]);
					for (int k = 0; k < gds.length; k++) {
						mitoList.add(gds[k].getChr() + "\t" + gds[k].getStart() + "\t" + gds[k].getStop() + "\t" + gds[k].getGeneName().toUpperCase());
					}
				}
			}
			log.reportTimeInfo("Count " + count + " mitos of " + mitos.length);
			Files.writeIterable(mitoList, dir + "mitoList-" + build + ".txt");

		}
		ArrayList<String> glist = new ArrayList<String>();
		for (int j = 0; j < geneTrack.getGenes().length; j++) {
			for (int j2 = 0; j2 < geneTrack.getGenes()[j].length; j2++) {
				GeneData gds = geneTrack.getGenes()[j][j2];
				glist.add(gds.getChr() + "\t" + gds.getStart() + "\t" + gds.getStop() + "\t" + gds.getGeneName().toUpperCase());

			}
		}
		Files.writeIterable(glist, dir + "geneList-" + build + ".txt");

		for (int filt = 0; filt < filters.length; filt++) {// for each filter crf found
			String[][] phenos = HashVec.loadFileToStringMatrix(dir + "pheno.dat", false, null, false);
			final String filtFile = dir + filters[filt];
			String out = dir + ext.rootOf(filters[filt]) + ".cnv";

			if (!Files.exists(out)) {
				FilterCalls.fromParameters(dir + filters[filt], new Logger());
			}
			String sampCNVs = ext.addToRoot(out, "." + ext.rootOf(sampFile));
			String geneSampCnvs = ext.addToRoot(out, "." + ext.rootOf(sampFile) + "_gene");

			HashSet<String> sampSet = HashVec.loadFileToHashSet(sampFile, false);
			if (!Files.exists(sampCNVs) || !Files.exists(geneSampCnvs)) {
				LocusSet<CNVariant> cnvs = CNVariant.loadLocSet(out, log);

				try {
					int numTotal = 0;
					int numSamps = 0;
					PrintWriter writer = new PrintWriter(new FileWriter(sampCNVs)); // subset to samp cnvs
					PrintWriter writerGene = new PrintWriter(new FileWriter(geneSampCnvs));// subset to samp and genes

					writer.println(Array.toStr(CNVariant.PLINK_CNV_HEADER));
					writerGene.println(Array.toStr(CNVariant.PLINK_CNV_HEADER));

					for (int i = 0; i < cnvs.getLoci().length; i++) {
						numTotal++;
						if (sampSet.contains(cnvs.getLoci()[i].getIndividualID())) {
							if (Double.isFinite(cnvs.getLoci()[i].getScore())) {
								writer.println(cnvs.getLoci()[i].toPlinkFormat());
								numSamps++;
								if (geneTrack.getOverlappingGenes(cnvs.getLoci()[i]).length > 0) {
									writerGene.println(cnvs.getLoci()[i].toPlinkFormat());
								}
							}
						}
					}
					// CNVariant cheater = new CNVariant("ALOUD_p_ARIC_batch15_005_affy_GenomeWideSNP_6_E09_236930.CEL.IND.txt", "ALOUD_p_ARIC_batch15_005_affy_GenomeWideSNP_6_E09_236930.CEL.IND.txt", (byte) 6, 155662711, 155675684, 1, 24.5930, 9, -99);
					//
					// log.reportTimeInfo("John remember you are adding this" + cheater.toPlinkFormat());
					// writer.println(cheater.toPlinkFormat());
					// writerGene.println(cheater.toPlinkFormat());

					writer.close();
					writerGene.close();
					log.reportTimeInfo(numTotal + " - > " + numSamps);
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}

			for (int i = 3; i < phenos[0].length; i++) { // test each phenotype

				final String pheno = phenos[0][i];
				final String key = phenos[0][i] + "_" + ext.rootOf(sampFile) + "_" + ext.rootOf(filters[filt]);
				final String opDir = dir + key + "/";
				boolean quant = false;

				new File(opDir).mkdirs();
				String phenoCNV = opDir + pheno + ".cnv";
				String phenoGene = opDir + pheno + "_gene.cnv";
				if (!Files.exists(phenoCNV) || !Files.exists(phenoGene)) {
					Files.copyFile(sampCNVs, phenoCNV);
					Files.copyFile(geneSampCnvs, phenoGene);

				}
				ArrayList<String> fam = new ArrayList<String>();
				ArrayList<String> quants = new ArrayList<String>();
				HashSet<String> types = new HashSet<String>();
				for (int j = 2; j < phenos.length; j++) {
					if (sampSet.contains(phenos[j][1])) {// filter pheno file
						fam.add(phenos[j][0] + "\t" + phenos[j][1] + "\t" + 0 + "\t" + 0 + "\t" + phenos[j][2] + "\t" + phenos[j][i]);
						quants.add(phenos[j][0] + "\t" + phenos[j][1] + "\t" + phenos[j][i]);
					}
					try {
						double val = Double.parseDouble(phenos[j][i]);
						if (val != 1 && val != 2)
							types.add(val + "");
					} catch (NumberFormatException nfe) {

					}
				}
				quant = types.size() > 0;
				if (quant) {
					log.reportTimeInfo("Detected phenotype " + pheno + " is quantatative");
				} else {
					log.reportTimeInfo("Detected phenotype " + pheno + " is discrete");
				}

				// prep .fam and pheno files
				Files.writeIterable(fam, opDir + pheno + ".fam");
				Files.writeIterable(fam, opDir + pheno + "_gene.fam");

				Files.writeIterable(quants, opDir + pheno + ".qPheno.txt");
				Files.writeIterable(quants, opDir + pheno + "_gene.qPheno.txt");
				String map = opDir + pheno + ".cnv.map";
				String mapGene = opDir + pheno + "_gene.cnv.map";

				if (!Files.exists(map) || !Files.exists(mapGene)) { // make the map
					CmdLine.run(dir + "plink --cnv-list " + pheno + ".cnv --cnv-make-map --out " + pheno, opDir);
					CmdLine.run(dir + "plink --cnv-list " + pheno + "_gene.cnv --cnv-make-map --out " + pheno + "_gene", opDir);

				}

				final int perm = 10000;
				final boolean oWrite = false;
				if (quant) {
					Callable<PlinkResult> c1 = new Callable<PlinkResult>() {

						@Override
						public PlinkResult call() throws Exception {
							ArrayList<String> cmd = new ArrayList<String>();
							cmd.add(dir + "plink");
							cmd.add("--cfile");
							cmd.add(opDir + pheno);
							cmd.add("--pheno");
							cmd.add(opDir + pheno + ".qPheno.txt");
							cmd.add("--mperm");
							cmd.add(perm + "");
							cmd.add("--out");
							cmd.add(opDir + pheno);
							// cmd.addAll(covarArgs);
							String out = opDir + pheno + ".cnv.qt.summary.mperm";
							boolean complete = CmdLine.runCommandWithFileChecks(Array.toStringArray(cmd), "", null, new String[] { out }, true, oWrite, true, log);
							return new PlinkResult(filtFile, key, out, Array.toStr(Array.toStringArray(cmd), " "), complete, true);
						}
					};
					Callable<PlinkResult> c2 = new Callable<PlinkResult>() {

						@Override
						public PlinkResult call() throws Exception {

							ArrayList<String> cmd = new ArrayList<String>();
							cmd.add(dir + "plink");
							cmd.add("--cfile");
							cmd.add(opDir + pheno);
							cmd.add("--pheno");
							cmd.add(opDir + pheno + ".qPheno.txt");
							cmd.add("--mperm");
							cmd.add(perm + "");
							cmd.add("--out");
							cmd.add(opDir + pheno + "_mito");
							cmd.add("--cnv-intersect");
							cmd.add(dir + "mitoList-" + build + ".txt");
							cmd.add("--cnv-test-region");
							// cmd.addAll(covarArgs);
							// cmd.add("--cnv-count");
							// cmd.add(dir + "glist-hg18.txt");
							// cmd.add("--cnv-subset");
							// cmd.add(dir + "MitoGenes.txt");
							// cmd.add("--cnv-enrichment-test");
							String out = opDir + pheno + "_mito.cnv.qt.summary.mperm";
							boolean complete = CmdLine.runCommandWithFileChecks(Array.toStringArray(cmd), "", null, new String[] { out }, true, oWrite, true, log);
							return new PlinkResult(filtFile, key + "_enrichment", out, Array.toStr(Array.toStringArray(cmd), " "), complete, true);
						}
					};

					Callable<PlinkResult> c3 = new Callable<PlinkResult>() {

						@Override
						public PlinkResult call() throws Exception {

							ArrayList<String> cmd = new ArrayList<String>();
							cmd.add(dir + "plink");
							cmd.add("--cfile");
							cmd.add(opDir + pheno + "_gene");
							cmd.add("--pheno");
							cmd.add(opDir + pheno + ".qPheno.txt");
							cmd.add("--mperm");
							cmd.add(perm + "");
							cmd.add("--out");
							cmd.add(opDir + pheno + "_mitoBurden");
							cmd.add("--cnv-count");
							cmd.add(dir + "geneList-" + build + ".txt");
							cmd.add("--cnv-subset");
							cmd.add(dir + "mitoCarta.txt");
							cmd.add("--cnv-enrichment-test");
							// cmd.addAll(covarArgs);
							String out = opDir + pheno + "_mitoBurden.cnv.burden.mperm";
							boolean complete = CmdLine.runCommandWithFileChecks(Array.toStringArray(cmd), "", null, new String[] { out }, true, oWrite, true, log);
							return new PlinkResult(filtFile, key + "_burdenenrichment", out, Array.toStr(Array.toStringArray(cmd), " "), complete, true);
						}
					};
					Callable<PlinkResult> c4 = new Callable<PlinkResult>() {

						@Override
						public PlinkResult call() throws Exception {

							ArrayList<String> cmd = new ArrayList<String>();
							cmd.add(dir + "plink");
							cmd.add("--cfile");
							cmd.add(opDir + pheno);
							cmd.add("--pheno");
							cmd.add(opDir + pheno + ".qPheno.txt");
							cmd.add("--mperm");
							cmd.add(perm + "");
							cmd.add("--out");
							cmd.add(opDir + pheno + "_mito_tradBurden");
							cmd.add("--cnv-intersect");
							cmd.add(dir + "mitoList-" + build + ".txt");
							cmd.add("--cnv-indiv-perm");

							String out = opDir + pheno + "_mito_tradBurden.cnv.indiv";
							System.out.println(Files.exists(opDir));
							boolean complete = CmdLine.runCommandWithFileChecks(Array.toStringArray(cmd), "", null, new String[] { out }, true, oWrite, true, log);

							String[] header = Files.getHeaderOfFile(out, "[\\s]+", log);
							ArrayList<String> summary = new ArrayList<String>();

							int[] cols = new int[] { 2, 3, 4, 5 };
							String[][] dataS = HashVec.loadFileToStringMatrix(out, true, cols, "[\\s]+", false, 1000, true);
							double[][] data = Array.toDoubleArrays(dataS, true);
							double[] phe = Matrix.extractColumn(data, 0);
							for (int j = 1; j < cols.length; j++) {
								double[] test = Matrix.extractColumn(data, j);
								LeastSquares ls = new LeastSquares(phe, test);
								String type = header[cols[j]];
								double p = ls.getSigs()[1];
								double b = ls.getBetas()[1];
								summary.add(type + "\t" + b + "\t" + p);
							}
							Files.writeIterable(summary, opDir + pheno + "_mito_tradBurden.cnv.indiv.sigs");

							return new PlinkResult(filtFile, key + "_enrichment", out, Array.toStr(Array.toStringArray(cmd), " "), false, true);
						}
					};

					Callable<PlinkResult> c5 = new Callable<PlinkResult>() {

						@Override
						public PlinkResult call() throws Exception {

							ArrayList<String> cmd = new ArrayList<String>();
							cmd.add(dir + "plink");
							cmd.add("--cfile");
							cmd.add(opDir + pheno);
							cmd.add("--pheno");
							cmd.add(opDir + pheno + ".qPheno.txt");
							cmd.add("--mperm");
							cmd.add(perm + "");
							cmd.add("--out");
							cmd.add(opDir + pheno + "_mito_tradBurden_count");
							cmd.add("--cnv-count");
							cmd.add(dir + "mitoList-" + build + ".txt");

							String out = opDir + pheno + "_mito_tradBurden_count.cnv.qt.summary";
							boolean complete = CmdLine.runCommandWithFileChecks(Array.toStringArray(cmd), "", null, new String[] { out }, true, oWrite, true, log);
							ArrayList<String> summary = new ArrayList<String>();

							String load = opDir + pheno + "_mito_tradBurden_count.cnv.indiv";
							String[] header = Files.getHeaderOfFile(load, "[\\s]+", log);

							int[] cols = new int[] { 2, 3, 4, 5, 6 };
							String[][] dataS = HashVec.loadFileToStringMatrix(load, true, cols, "[\\s]+", false, 1000, true);
							double[][] data = Array.toDoubleArrays(dataS, true);
							double[] phe = Matrix.extractColumn(data, 0);
							for (int j = 1; j < cols.length; j++) {
								double[] test = Matrix.extractColumn(data, j);
								LeastSquares ls = new LeastSquares(phe, test);
								String type = header[cols[j]];
								double p = ls.getSigs()[1];
								double b = ls.getBetas()[1];
								summary.add(type + "\t" + b + "\t" + p);
							}
							Files.writeIterable(summary, opDir + pheno + "_mito_tradBurden_count.cnv.indiv.sigs");

							return new PlinkResult(filtFile, key + "_enrichment", out, Array.toStr(Array.toStringArray(cmd), " "), false, true);
						}
					};

					Callable<PlinkResult> c6 = new Callable<PlinkResult>() {

						@Override
						public PlinkResult call() throws Exception {

							ArrayList<String> cmd = new ArrayList<String>();
							cmd.add(dir + "plink");
							cmd.add("--cfile");
							cmd.add(opDir + pheno);
							cmd.add("--pheno");
							cmd.add(opDir + pheno + ".qPheno.txt");
							cmd.add("--mperm");
							cmd.add(perm + "");
							cmd.add("--out");
							cmd.add(opDir + pheno + "_FulltradBurden");
							cmd.add("--cnv-indiv-perm");

							String out = opDir + pheno + "_FulltradBurden.cnv.indiv";
							System.out.println(Files.exists(opDir));
							boolean complete = CmdLine.runCommandWithFileChecks(Array.toStringArray(cmd), "", null, new String[] { out }, true, oWrite, true, log);

							String[] header = Files.getHeaderOfFile(out, "[\\s]+", log);
							ArrayList<String> summary = new ArrayList<String>();

							int[] cols = new int[] { 2, 3, 4, 5 };
							String[][] dataS = HashVec.loadFileToStringMatrix(out, true, cols, "[\\s]+", false, 1000, true);
							double[][] data = Array.toDoubleArrays(dataS, true);
							double[] phe = Matrix.extractColumn(data, 0);
							for (int j = 1; j < cols.length; j++) {
								double[] test = Matrix.extractColumn(data, j);
								LeastSquares ls = new LeastSquares(phe, test);
								String type = header[cols[j]];
								double p = ls.getSigs()[1];
								double b = ls.getBetas()[1];
								summary.add(type + "\t" + b + "\t" + p);
							}
							Files.writeIterable(summary, opDir + pheno + "_FulltradBurden.cnv.indiv.sigs");

							return new PlinkResult(filtFile, key + "_enrichment", out, Array.toStr(Array.toStringArray(cmd), " "), false, true);
						}
					};

					// plink --cfile mycnv
					//  --cnv-count genes.dat
					//  --cnv-subset pathway.txt
					//  --cnv-enrichment-test

					hive.addCallable(c1);
					hive.addCallable(c2);
					hive.addCallable(c3);
					hive.addCallable(c4);
					hive.addCallable(c5);
					hive.addCallable(c6);

				} else {
					Callable<PlinkResult> c1 = new Callable<PlinkResult>() {

						@Override
						public PlinkResult call() throws Exception {
							ArrayList<String> cmd = new ArrayList<String>();
							cmd.add(dir + "plink");
							cmd.add("--cfile");
							cmd.add(opDir + pheno);
							cmd.add("--cnv-indiv-perm");
							cmd.add("--mperm");
							cmd.add("10000");
							cmd.add("--out");
							cmd.add(opDir + pheno);
							String out = opDir + pheno + ".cnv.summary.mperm";
							boolean complete = CmdLine.runCommandWithFileChecks(Array.toStringArray(cmd), "", null, new String[] { out }, true, false, true, log);

							return new PlinkResult(filtFile, key, out, Array.toStr(Array.toStringArray(cmd), " "), complete, false);
						}
					};
					Callable<PlinkResult> c2 = new Callable<PlinkResult>() {

						@Override
						public PlinkResult call() throws Exception {
							ArrayList<String> cmd = new ArrayList<String>();
							cmd.add(dir + "plink");
							cmd.add("--cfile");
							cmd.add(opDir + pheno);
							cmd.add("--mperm");
							cmd.add("10000");
							cmd.add("--out");
							cmd.add(opDir + pheno + "_position");
							String out = opDir + pheno + "_position.cnv.summary.mperm";

							boolean complete = CmdLine.runCommandWithFileChecks(Array.toStringArray(cmd), "", null, new String[] { out }, true, false, true, log);
							return new PlinkResult(filtFile, key + "_pos", out, Array.toStr(Array.toStringArray(cmd), " "), complete, false);

						}
					};
					Callable<PlinkResult> c3 = new Callable<PlinkResult>() {

						@Override
						public PlinkResult call() throws Exception {
							ArrayList<String> cmd = new ArrayList<String>();
							cmd.add(dir + "plink");
							cmd.add("--cfile");
							cmd.add(opDir + pheno);
							cmd.add("--mperm");
							cmd.add("10000");
							cmd.add("--cnv-test-window");
							cmd.add("200");
							cmd.add("--out");
							cmd.add(opDir + pheno + "_window");
							String out = opDir + pheno + "_window.cnv.summary.mperm";

							boolean complete = CmdLine.runCommandWithFileChecks(Array.toStringArray(cmd), "", null, new String[] { out }, true, false, true, log);
							return new PlinkResult(filtFile, key + "_window", out, Array.toStr(Array.toStringArray(cmd), " "), complete, false);

						}
					};
					hive.addCallable(c1);
					hive.addCallable(c2);
					hive.addCallable(c3);
				}
			}

		}
		hive.execute(true);
		ArrayList<PlinkResult> results = hive.getResults();
		ArrayList<String> filesToCombine = new ArrayList<String>();

		String finalDir = dir + ext.rootOf(sampFile) + "_results/";
		filesToCombine.add(finalDir + "key.txt");
		ArrayList<String> key = new ArrayList<String>();
		ArrayList<String> combo = new ArrayList<String>();
		String comboOut = finalDir + "combo.txt";
		String allSigGenes = finalDir + "sigGenes.txt";
		HashSet<String> allSigGenesList = new HashSet<String>();
		String allSigGenes5k = finalDir + "sigGenes5k.txt";
		HashSet<String> allSigGenesList5k = new HashSet<String>();

		String allSigGenesEMP2 = finalDir + "sigGenesEMP2_0_5.txt";
		HashSet<String> allSigGenesListEMP2 = new HashSet<String>();

		filesToCombine.add(comboOut);
		filesToCombine.add(allSigGenesEMP2);
		filesToCombine.add(allSigGenes);
		filesToCombine.add(allSigGenes5k);

		combo.add("Type\tCHR\tSNP\tEMP1\tEMP2\tUCSC\tGENE\tGENE_5k\tLink_buffer\tAFF\tUNAFF\tNCNV\tM0\tM1\tCMD");
		new File(finalDir).mkdirs();
		int index = 1;
		for (PlinkResult result : results) {
			if (result.complete) {
				String out = finalDir + result.key.substring(0, 26) + "_" + index + ".txt";
				String[][] data = HashVec.loadFileToStringMatrix(result.file, false, null, "[\\s]+", false, 1000, true);
				log.reportTimeInfo("Loading paired " + ext.rootOf(result.file, false));
				String[][] dataCount = HashVec.loadFileToStringMatrix(ext.rootOf(result.file, false), false, null, "[\\s]+", false, 1000, true);
				String allSigGenesSpecific = finalDir + result.key + "_Genes" + index + ".txt";
				HashSet<String> allSigGenesListSpecific = new HashSet<String>();

				ArrayList<String> toReport = new ArrayList<String>();
				toReport.add("Type\t" + Array.toStr(data[0]) + "\tUCSC\tGENE\tGENE_5k\tLink_buffer\tAFF\tUNAFF\tNCNV\tM0\tM1\tCMD");

				for (int i = 1; i < data.length; i++) {
					if (data[i].length > 3 && Double.parseDouble(data[i][result.quant ? 3 : 2]) < 0.05) {
						int loc = Integer.parseInt(data[i][1].replaceAll("p" + data[i][0] + "-", ""));
						Segment seg = new Segment(Byte.parseByte(data[i][0]), loc, loc);
						GeneData[] geneDatas = geneTrack.getOverlappingGenes(seg);
						GeneData[] geneDatasBuff = geneTrack.getOverlappingGenes(seg.getBufferedSegment(5000));
						if (!data[i][1].equals(dataCount[i][1])) {
							throw new IllegalArgumentException("mismatched files");
						}
						String[] dataToReport = data[i];
						if (result.quant) {
							dataToReport = Array.removeFromArray(dataToReport, 2);
						}
						StringBuilder builder = new StringBuilder(result.key + "\t" + Array.toStr(dataToReport) + "\t" + seg.getUCSClocation() + "\t");
						for (int j = 0; j < geneDatas.length; j++) {
							builder.append((j == 0 ? "" : ":") + geneDatas[j].getGeneName());
							allSigGenesList.add(geneDatas[j].getGeneName() + "\t" + geneDatas[j].getUCSClocation() + "\t" + result.key);
							allSigGenesListSpecific.add(geneDatas[j].getGeneName() + "\t" + geneDatas[j].getUCSClocation() + "\t" + result.key);
							if (Double.parseDouble(data[i][result.quant ? 4 : 3]) < 1) {
								allSigGenesListEMP2.add(geneDatas[j].getGeneName() + "\t" + geneDatas[j].getUCSClocation() + "\t" + result.key);
							}
						}
						builder.append("\t");
						for (int j = 0; j < geneDatasBuff.length; j++) {
							builder.append((j == 0 ? "" : ":") + geneDatasBuff[j].getGeneName());
							allSigGenesList5k.add(geneDatasBuff[j].getGeneName() + "\t" + geneDatasBuff[j].getUCSClocation() + "\t" + result.key);
						}
						builder.append("\t" + seg.getBufferedSegment(5000).getUCSCLink(build));
						if (result.quant) {
							builder.append("\tNA\tNA" + "\t" + dataCount[i][3] + "\t" + dataCount[i][4] + "\t" + dataCount[i][5]);

						} else {
							builder.append("\t" + dataCount[i][3] + "\t" + dataCount[i][4] + "\tNA\tNA\tNA");
						}
						builder.append("\t" + result.command);
						toReport.add(builder.toString());
						combo.add(builder.toString());
					}
				}
				if (toReport.size() > 1) {

					key.add(result.key + "\t" + index);
					index++;
					Files.writeIterable(toReport, out);
					Files.writeIterable(allSigGenesListSpecific, allSigGenesSpecific);
					filesToCombine.add(out);
				} else {
					log.reportTimeInfo("Nothing to report for " + result.key);
				}
			} else {
				log.reportTimeError(result.key + " was not complete");
			}
		}

		Files.writeIterable(key, finalDir + "key.txt");
		Files.writeIterable(combo, comboOut);
		Files.writeIterable(allSigGenesList, allSigGenes);
		Files.writeIterable(allSigGenesList5k, allSigGenes5k);
		Files.writeIterable(allSigGenesListEMP2, allSigGenesEMP2);

		ExcelConverter excelConverter = new ExcelConverter(filesToCombine, finalDir + "p_0_05.xlsx", log);
		excelConverter.convert(true);
		System.exit(1);
	}

	private static class PlinkResult {
		private String filtFile;
		private String file;
		private String key;
		private String command;
		private boolean quant;
		private boolean complete;

		public PlinkResult(String filtFile, String key, String file, String command, boolean complete, boolean quant) {
			super();
			this.filtFile = filtFile;
			this.key = key;
			this.file = file;
			this.command = command;
			this.complete = complete;
			this.quant = quant;
		}

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String dir = "C:/data/ARIC/shadowCNVs/";
		String cnvFile = dir + "combinedMF.cnv";
		String build = "hg18";
		String usage = "\n" +
				"one.JL.ARICCNV requires 0-1 arguments\n" +
				"   (1) cnvs (i.e. cnvs=" + cnvFile + " (default))\n" +
				"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				cnvFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("dir=")) {
				dir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("build=")) {
				build = args[i].split("=")[1];
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
			String[] sampFiles = new String[] { dir + "all.txt", dir + "whites.txt", dir + "blacks.txt" };

			for (int i = 0; i < sampFiles.length; i++) {
				run(cnvFile, sampFiles[i], build);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
