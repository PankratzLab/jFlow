package org.genvisis.seq.manage;

import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Hashtable;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.Positions;
import org.genvisis.common.ext;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.manage.VCFOps.VcfPopulation;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import org.genvisis.seq.qc.FilterNGS.RareVariantFilter;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.tribble.Tribble;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.util.LittleEndianOutputStream;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.VariantContextUtils;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class VCF {
	// private static final String SITE_ONLY = ".siteOnly";
	// private static final String BLANK_ANNO = ".";
	private static final Set<String> SAMPLE = new TreeSet<String>();
	private static final String[] testExpress = {"CHROM == 'chr1' && DP > 0 && DP < 100"};
	private static final String[] testExpress2 = {"ExonicFunc.refGene != '.' && ExonicFunc.refGene == 'nonsynonymous_SNV'"};
	private final VCFFileReader vcfFileReader;
	private VariantContextWriter variantContextWriter;
	private final String vcfFile;
	private final boolean fail;
	private boolean extractBams;
	private final Logger log;
	public static final String VCF_INIT = "VCFilt";
	public static final String VCF_DESCRIPTION = "- filter a vcf file using JEXL expressions";
	public static final String VCF_COMMAND = "vcf=";
	public static final String EXPRESSION_COMMAND = "filter=";
	public static final String EXPRESSION_NAME_COMMAND = "name=";
	public static final String ID_FILE_COMMAND = "IDFile=";
	public static final String EXTRACT_ANNOTATION_COMMAND = "-extract";
	public static final String BAM_DIR_COMMAND = "bamDir=";
	public static final String DUMP_COMMAND = "-dump";
	public static final String TO_DUMP_COMMAND = "dump=";

	public VCF(String vcfFile, Logger log) {
		this.log = log;
		// this.fail = verifyIndex(vcfFile, log);
		fail = false;

		this.vcfFile = vcfFile;
		vcfFileReader = new VCFFileReader(new File(vcfFile), true);
		variantContextWriter = null;
		extractBams = false;
	}

	public boolean isFail() {
		return fail;
	}

	public String[] getAvailableAnno() {
		Collection<VCFInfoHeaderLine> vcfInfoHeaderLines = vcfFileReader.getFileHeader()
																																		.getInfoHeaderLines();
		String[] infos = new String[vcfInfoHeaderLines.size()];
		int index = 0;
		for (VCFInfoHeaderLine vcfInfoHeaderLine : vcfInfoHeaderLines) {
			infos[index] = (vcfInfoHeaderLine.getID()	+ "\t" + vcfInfoHeaderLine.getType() + "\t"
											+ vcfInfoHeaderLine.getDescription() + "\t");
			index++;
		}
		return infos;
	}

	public String[] getSamplesInVcf() {
		return Array.toStringArray(vcfFileReader.getFileHeader().getSampleNamesInOrder());
	}

	public void closeReader() {
		vcfFileReader.close();
	}

	public void closeWriter() {
		variantContextWriter.close();
	}

	public boolean hasAllInfos(String[] toDump) {
		if (toDump == null) {
			return true;
		}
		for (int i = 0; i < toDump.length; i++) {
			if (!vcfFileReader.getFileHeader().hasInfoLine(toDump[i])) {
				log.reportTimeError("Could not find common info field " + toDump[i] + " in " + vcfFile);
				return false;
			}
		}
		return true;
	}

	public void filter(	String[] jexpression, VcfPopulation vpop, String[] name, String[] toDump,
											String bamDir, String outputDir,
											Hashtable<String, Vector<String>> IDsToExtract, String segFile,
											int numThreads, int bpBuffer, int mac) {
		BamExtractor.BamSample bamSample = null;
		ArrayList<String> toDumpTmp = new ArrayList<String>();
		ArrayList<String> keys = new ArrayList<String>();
		RareVariantFilter rareVariantFilter = null;
		if (vpop != null) {
			rareVariantFilter = new RareVariantFilter(vpop.getSubPop().get("CASE"),
																								vpop.getSubPop().get("CONTROL"));
			rareVariantFilter.setMafRef(0.01);
			rareVariantFilter.setMacCase(mac);
			rareVariantFilter.initFilters(log);
		}
		Hashtable<String, Integer> histCOUNT = new Hashtable<String, Integer>();
		Hashtable<String, Integer> histALLELE = new Hashtable<String, Integer>();
		Segment[] segs = null;
		if (segFile != null) {
			segs = Segment.loadRegions(segFile, 0, 1, 2, 0, true, true, true, 0);
		}
		if (hasAllInfos(toDump)) {
			if (toDump != null) {
				log.reportTimeInfo("Will be dumping the following annotations " + Array.toStr(toDump, ","));
				toDumpTmp.add("CHR\tPOS\t" + Array.toStr(toDump) + "\tNumAlleles");
			}
			if (outputDir == null) {
				outputDir = ext.parseDirectoryOfFile(vcfFile);
			}
			new File(outputDir).mkdirs();
			String outputVCF = outputDir + ext.removeDirectoryInfo(vcfFile).replaceAll(".vcf", "")
																				.replaceAll(".gz", "");
			outputVCF = ext.addToRoot(outputVCF, "." + Array.toStr(name, "_"));
			if (!outputVCF.endsWith(".vcf.gz")) {
				outputVCF = outputVCF + ".vcf.gz";
			}
			variantContextWriter = initWriter(vcfFileReader, outputVCF, false);
			if (bamDir != null) {
				log.reportTimeInfo("Since a bam directory was provided, we will verify that all samples in the vcf have a corresponding bam file prior to filtering");
				extractBams = true;
				bamSample =
									new BamExtractor.BamSample(Files.listFullPaths(bamDir, ".bam", false), log, true);
				bamSample.generateMap();
				bamSample.getBamSampleMap();
				if ((!bamSample.isFail()) && (bamSample.verify(getSamplesInVcf(), null))) {
				}
			} else {
				extractBams = false;
				log.reportTimeInfo("Since a bam directory was not provided, we will not subset the bam files");
			}
			int count = 0;
			int countPass = 0;
			List<VariantContextUtils.JexlVCMatchExp> jExps = null;
			if (jexpression != null && jexpression[0] != null && name != null) {
				jExps = VariantContextUtils.initializeMatchExps(name, jexpression);
				log.reportTimeInfo("Using " + jExps.size() + " expression(s)");
			}
			// VariantContextUtils.match(vc, g, exps)
			// int report = 0;

			for (VariantContext variantContext : vcfFileReader) {
				// variantContext.getGenotype(1).
				boolean write = true;
				if (jExps != null) {
					for (VariantContextUtils.JexlVCMatchExp jExp : jExps) {
						if (!VariantContextUtils.match(variantContext, jExp)) {
							write = false;
							break;
						}
					}
				}
				if (variantContext.isFiltered()) {
					write = false;
				}
				if (write && rareVariantFilter != null) {
					write = rareVariantFilter.filter(variantContext, log).passed();
				}
				if (write && (IDsToExtract.size() > 0)) {
					write = IDsToExtract.containsKey(variantContext.getID());
					// if (!write && variantContext.getCommonInfo().hasAttribute("snp138")) {
					// write =
					// IDsToExtract.containsKey(variantContext.getCommonInfo().getAttributeAsString("snp138",
					// "."));
					// }
				}
				if (write && segs != null) {
					write = VCOps.isInTheseSegments(variantContext, segs);
				}
				if ((count != 0) && (count % 100000 == 0)) {
					log.report(ext.getTime() + " Info - scanned " + count + " variants");
					log.report(ext.getTime()	+ " Info - " + countPass
											+ " variants have passed the filter thus far...currently on chromosome "
											+ variantContext.getContig());
					if (IDsToExtract.size() > 0) {
						log.report(ext.getTime()	+ " Info - " + countPass + " variants ( " + IDsToExtract.size()
												+ " eligible ) have passed the filter(s) thus far...");
					}
				}
				if (write) {
					countPass++;
					variantContextWriter.add(variantContext);
					byte chr = Positions.chromosomeNumber(variantContext.getContig());
					int start = variantContext.getStart();
					if (extractBams) {
						bamSample.addSegmentToExtract(new Segment(chr, start, start));
					}
					if (toDump != null) {
						if (vpop != null) {
							variantContext = VCOps.getSubset(variantContext, vpop.getSubPop().get("CASE"));
						}

						String tmp = variantContext.getContig() + "\t" + variantContext.getStart();
						int numAlts = variantContext.getHomVarCount() * 2 + variantContext.getHetCount();
						for (String element : toDump) {
							String key = element	+ "\t"
														+ variantContext.getCommonInfo().getAttributeAsString(element, ".");
							try {
								Double.parseDouble(variantContext.getCommonInfo().getAttributeAsString(	element,
																																												"."));
							} catch (NumberFormatException nfe) {
								if (!histCOUNT.containsKey(key)) {
									histALLELE.put(key, Integer.valueOf(0));
									histCOUNT.put(key, Integer.valueOf(0));
									keys.add(key);
								}
								histCOUNT.put(key, Integer.valueOf(histCOUNT.get(key).intValue() + 1));
								histALLELE.put(key, Integer.valueOf(histALLELE.get(key).intValue() + numAlts));
							}
							tmp = tmp + "\t" + variantContext.getCommonInfo().getAttributeAsString(element, ".");
						}
						toDumpTmp.add(tmp + "\t" + numAlts);
					}
				}
				count++;
			}
			log.reportTimeInfo(countPass + " of " + count + " varaints passed the fileter");
			vcfFileReader.close();
			variantContextWriter.close();
			if (toDump != null) {
				Files.writeArray(	Array.toStringArray(toDumpTmp),
													ext.rootOf(outputVCF, false) + ".filteredAnno");
				try {
					PrintWriter writer = new PrintWriter(new FileWriter(ext.rootOf(outputVCF, false)
																															+ ".hist"));
					for (int i = 0; i < keys.size(); i++) {
						writer.println(keys.get(i)	+ "\t" + histCOUNT.get(keys.get(i)) + "\t"
														+ histALLELE.get(keys.get(i)));
					}
					writer.close();
				} catch (Exception e) {
					log.reportError("Error writing to " + ext.rootOf(outputVCF, false) + ".hist");
					log.reportException(e);
				}
			}
			if (extractBams) {
				BamExtractor.extractAll(bamSample, outputDir, bpBuffer, true, true, numThreads, log);
				bamSample = new BamExtractor.BamSample(	Files.listFullPaths(outputDir, ".bam", false), log,
																								true);
				bamSample.generateMap();
				bamSample.dumpToIGVMap(outputVCF, null);
			}
		}
	}

	public void dumpToSiteOnly() {
		variantContextWriter = initWriter(vcfFileReader, ext.addToRoot(vcfFile, ".siteOnly"), true);
		for (VariantContext variantContext : vcfFileReader) {
			variantContextWriter.add(subsetToSamplesWithOriginalAnnotations(variantContext, SAMPLE));
		}
	}

	public void convertToPlinkSet(String vcf) {
		String rootOut = ext.rootOf(vcfFile, false);
		String[] outFiles = PSF.Plink.getPlinkBedBimFam(rootOut);
		String[] plinkCommand = PSF.Plink.getPlinkVCFCommand(vcfFile, rootOut);
		CmdLine.runCommandWithFileChecks(	plinkCommand, "", new String[] {vcfFile}, outFiles, true, true,
																			false, log);
	}

	private static VariantContextWriter initWriter(	VCFFileReader vcfFileReader, String output,
																									boolean siteOnly) {
		VCFHeader inputVcfHeader = siteOnly	? new VCFHeader(vcfFileReader	.getFileHeader()
																																			.getMetaDataInInputOrder())
																				: vcfFileReader.getFileHeader();
		SAMSequenceDictionary sequenceDictionary = inputVcfHeader.getSequenceDictionary();
		VariantContextWriterBuilder builder =
																				new VariantContextWriterBuilder()	.setOutputFile(output)
																																					.setReferenceDictionary(sequenceDictionary);
		builder.setOption(Options.INDEX_ON_THE_FLY);
		VariantContextWriter writer = builder.build();
		writer.writeHeader(siteOnly	? new VCFHeader(inputVcfHeader.getMetaDataInInputOrder(), SAMPLE)
																: inputVcfHeader);
		return writer;
	}

	private static VariantContext subsetToSamplesWithOriginalAnnotations(	VariantContext ctx,
																																				Set<String> samples) {
		VariantContextBuilder builder = new VariantContextBuilder(ctx);
		GenotypesContext newGenotypes = ctx.getGenotypes().subsetToSamples(samples);
		builder.alleles(ctx.getAlleles());
		return builder.genotypes(newGenotypes).make();
	}

	public static boolean verifyIndex(String vcfFile, Logger log) {
		boolean created = false;
		File indexFile = Tribble.indexFile(new File(vcfFile));
		if (indexFile.canRead()) {
			log.report("Info - Loading index file " + indexFile);
			IndexFactory.loadIndex(indexFile.getAbsolutePath());
			created = true;
		} else {
			log.report("Info - creating index file " + indexFile);
			try {
				Index index = IndexFactory.createLinearIndex(new File(vcfFile), new VCFCodec());
				LittleEndianOutputStream stream =
																				new LittleEndianOutputStream(new FileOutputStream(indexFile));
				index.write(stream);
				stream.close();
				created = true;
			} catch (IOException e) {
				log.reportError("Error - could not create index file " + indexFile);
				created = false;
			}
		}
		return created;
	}

	public static String[] getParserParams(String dir) {
		String[] params = new String[23];
		params[0] = "#the full path to a vcf file";
		params[1] = ("vcf=" + dir);
		params[2] = "# extract the annotations available in the vcf";
		params[3] = "-extract";
		params[4] = "# dump the vcf to a site only context";
		params[5] = "#-dump";
		params[6] = "# JEXL formatted filter expression";
		params[7] = "#filter=";
		params[8] = "# Name of the JEXL filter command";
		params[9] = "#name=";
		params[10] = "# Directory of bam files to use ";
		params[11] = "#bamDir=";
		params[12] = "# base pair buffer to use ";
		params[13] = "#bpBuffer=";
		params[14] = "# number of threads ";
		params[15] = "#numthreads=";
		params[16] = "# Below are some example filters... ";
		params[17] = ("#" + testExpress[0]);
		params[18] = ("#" + testExpress2[0]);
		params[19] = "# file name of IDs to extract (filters will be applied to only this subset)";
		params[20] = "#IDFile=";
		params[21] = "# INFO to dump for filtered variants";
		params[22] = "#dump=";
		return params;
	}

	public static String getIgvXmlScript(	String miniSamDir, String chr, String pos,
																				String[] miniSamFilenamesOfOneTrio) {
		return "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n<Session genome=\"hg19\" hasGeneTrack=\"true\" hasSequenceTrack=\"true\" locus=\"chr"
							+ chr + ":" + pos + "\" version=\"8\">" + "\n<Resources>" + "\n<Resource path=\""
						+ miniSamDir + miniSamFilenamesOfOneTrio[1] + "\"/>" + "\n<Resource path=\""
						+ miniSamDir + miniSamFilenamesOfOneTrio[2] + "\"/>" + "\n<Resource path=\""
						+ miniSamDir + miniSamFilenamesOfOneTrio[3] + "\"/>" + "\n</Resources>" + "</Session>";
	}

	public static String getIgvLaunchScript(String fulPathToXml) {
		return "java -Xmx1200m -Dproduction=true -Djava.net.preferIPv4Stack=true -Dsun.java2d.noddraw=true -jar D:/logan/DeNovos/IGV/IGV_2.3.36/igv.jar "
						+ fulPathToXml;
	}

	public static void extractAvaliableAnnotations(String vcfFile, Logger log) {
		VCF vcf = new VCF(vcfFile, new Logger());
		String[] annos = vcf.getAvailableAnno();
		Files.writeArray(annos, ext.rootOf(vcfFile, false) + ".anno");
		vcf.closeReader();
	}

	public static void dumpToSiteOnly(String vcfFile, Logger log) {
		VCF vcf = new VCF(vcfFile, new Logger());
		vcf.dumpToSiteOnly();
		vcf.closeWriter();
		vcf.closeReader();
	}

	public static void filterByExpression(String vcfFile, String popFile, String jexp,
																				String jexpName, String[] toDump, String bamDir,
																				String outputDir, String idFile, String segFile,
																				int numThreads, int bpBuffer, int mac, Logger log) {
		if (jexp != null) {
			log.report(ext.getTime() + " Info - using expression " + jexp + " on " + vcfFile);
		} else {
			log.reportTimeInfo("No filter expressions were supplied");
		}
		VcfPopulation vpop = null;
		if (popFile != null) {
			System.out.println(popFile);
			vpop = VcfPopulation.load(popFile, POPULATION_TYPE.ANY, log);
			vpop.report();
		}
		VCF vcf = new VCF(vcfFile, new Logger(ext.rootOf(vcfFile, false) + ".log"));
		Hashtable<String, Vector<String>> IDsToExtract = new Hashtable<String, Vector<String>>();
		if (idFile != null) {
			IDsToExtract = HashVec.loadFileToHashVec(idFile, 0, new int[1], "\t", false, true);
			log.reportTimeInfo("Subsetting the search to " + IDsToExtract.size() + " ID(s)");
		}
		vcf.filter(	new String[] {jexp == null ? null : jexp}, vpop, new String[] {jexpName}, toDump,
								bamDir, outputDir, IDsToExtract, segFile, numThreads, bpBuffer, mac);
	}

	public static void fromParameters(String filename, Logger log) {
		Vector<String> params = Files.parseControlFile(	filename, "VCFilt",
																										getParserParams(ext.parseDirectoryOfFile(filename)),
																										log);
		if (params != null) {
			main(Array.toStringArray(params));
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String vcfFile = "D:/data/Project_Tsai_Project_021/Variants/joint_genotypes.SNP.recal.INDEL.recal.eff.gatk.vcf";

		String logfile = null;
		String bamDir = null;
		String outputDir = null;
		String idFile = null;

		String[] toDump = null;
		ArrayList<String> filterExpression = new ArrayList<String>();
		String filterName = "filter";
		int numThreads = 1;
		int bpBuffer = 1000;
		boolean extractAnnotation = false;
		String vpop = null;
		boolean dump = false;
		String segFile = null;
		int mac = 2;
		String usage = "\njlDev.VCF requires 0-1 arguments\n";
		usage = usage + "   (1) vcf filename (i.e.vcf=" + vcfFile + " (default))\n";
		usage = usage + "   (2) filter expression to use on the vcf file (i.e. filter= (no default))\n";
		usage = usage + "   (3) filter name (i.e. name=" + filterName + " (no default))\n";
		usage = usage
							+ "   (4) bam directory containing .bam files to match with variants in the vcf, defaults to not matching (i.e. bamDir="
						+ bamDir + " (no default))\n";
		usage =
					usage	+ "   (5) output directory, defualts to directory of the vcf file (i.e. outputdir="
						+ bamDir + " (no default))\n";
		usage = usage	+ "   (6) number of threads to use if extracting bams (i.e. numthreads="
						+ numThreads + " ( default))\n";
		usage = usage	+ "   (7) up and downstream base-pair buffer if extracting bams (i.e. bpBuffer="
						+ bpBuffer + " ( default))\n";
		usage = usage
						+ "   (8) ids to filter from the \"ID\" column (note that any other filters will also be applied) (i.e. IDFile= (no default))\n";
		usage = usage
						+ "   (9) dump these annotations for variants passing the filters (comma separated, must be in INFO column)  (i.e. dump= (no default))\n";
		usage = usage + "   (10) A population File for call rate   (i.e. vpop= (no default))\n";
		usage = usage + "   (11) mac for case    (i.e. mac=" + mac + " ( default))\n";

		usage = usage + "   OR:";
		usage = usage	+ "   (1) extract available annotations from a vcf for filtering (i.e. name="
						+ filterName + " (no default))\n";
		usage = usage + "   OR:";
		usage = usage	+ "   (1) dump to a site only vcf (no genotypes) (i.e. -dump" + filterName
						+ " (no default))\n";
		for (String arg : args) {
			if ((arg.equals("-h"))	|| (arg.equals("-help")) || (arg.equals("/h"))
					|| (arg.equals("/help"))) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("vcf=")) {
				vcfFile = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith("filter=")) {
				filterExpression.add(ext.parseStringArg(arg, ""));
				numArgs--;
			} else if (arg.startsWith("name=")) {
				filterName = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith("-extract")) {
				extractAnnotation = true;
				numArgs--;
			} else if (arg.startsWith("bamDir=")) {
				bamDir = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith("outputdir=")) {
				outputDir = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith("dump=")) {
				toDump = ext.parseStringArg(arg, "").split(",");
				numArgs--;
			} else if (arg.startsWith("IDFile=")) {
				idFile = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith("segs=")) {
				segFile = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith("vpop=")) {
				vpop = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith("numthreads=")) {
				numThreads = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith("bpBuffer=")) {
				bpBuffer = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith("mac=")) {
				mac = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith("-dump")) {
				dump = true;
				numArgs--;
			} else if (arg.startsWith("log=")) {
				logfile = arg.split("=")[1];
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
			Logger log = new Logger(logfile);
			String expression = null;
			if (filterExpression.size() > 0) {
				expression = Array.toStr(	filterExpression.toArray(new String[filterExpression.size()]),
																	"&&");
			}
			if (extractAnnotation) {
				extractAvaliableAnnotations(vcfFile, log);
			}
			if (dump) {
				dumpToSiteOnly(vcfFile, log);
			} else {
				filterByExpression(	vcfFile, vpop, expression, filterName, toDump, bamDir, outputDir, idFile,
														segFile, numThreads, bpBuffer, mac, log);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}

// TODO dynamic

// VariantContext cushings = VCOps.getSubset(variantContext,
// vpop.getSubPop().get("CUSHINGS_JOINT_EUR"));
// VariantContext osteoParents = VCOps.getSubset(variantContext,
// vpop.getSubPop().get("OSTEO_Parents_JOINT_EUR"));
// if (VCOps.getCallRate(cushings, null) >= .90 && VCOps.getCallRate(osteoParents, null) >= .90) {
// double osteoMAF = VCOps.getMAF(osteoParents, null);
// double mafCushingCounts = VCOps.getMAC(cushings, null);
// if (osteoMAF <= 0.01 && mafCushingCounts >= 2) {
// // if (hweFilter.filter(cushings).passed() && hweFilter.filter(variantContext).passed()) {
// write = true;
// if (report == 0) {
// log.reportTimeError("DONT forget about the popFiltering");
// }
// report++;
// if (report % 1000 == 0) {
// double osteoMAFCounts = VCOps.getMAC(osteoParents, null);
// log.reportTimeInfo("NUM CUSHINGS: " + cushings.getNSamples());
// log.reportTimeInfo("NUM OSTEO: " + osteoParents.getNSamples());
//
// log.reportTimeInfo("CALLRATE CUSHINGS: " + VCOps.getCallRate(cushings, null));
// log.reportTimeInfo("CALLRATE OSTEO: " + VCOps.getCallRate(osteoParents, null));
//
// log.reportTimeInfo("MAF Count CUSHINGS: " + mafCushingCounts);
// log.reportTimeInfo("MAF OSTEO: " + osteoMAF);
// log.reportTimeInfo("MAF OSTEO Counts: " + osteoMAFCounts);
//
// }
// } else {
// write = false;
// }
// } else {
// write = false;
// }
