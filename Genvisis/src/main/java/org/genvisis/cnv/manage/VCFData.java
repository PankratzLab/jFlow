package org.genvisis.cnv.manage;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.genvisis.CLI;
import org.genvisis.CLI.Arg;
import org.genvisis.cnv.filesys.ClusterFilterCollection;
import org.genvisis.cnv.filesys.MarkerData;
import org.genvisis.cnv.filesys.MarkerDetailSet;
import org.genvisis.cnv.filesys.MarkerDetailSet.Marker;
import org.genvisis.cnv.filesys.MarkerDetailSet.Marker.RefAllele;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.Positions;
import org.genvisis.common.Sort;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.manage.ReferenceGenome;
import org.genvisis.seq.manage.StrandOps;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.tribble.annotation.Strand;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;

public final class VCFData {

	private VCFData() {
		// private un-instantiable class
	}

	public static void exportGenvisisToVCF(Project proj, String fileOfSamplesToExport,
																				 String fileOfMarkersToExport,
																				 boolean splitChrs, boolean useGRCRefGen,
																				 int[] chrsToExport,
																				 String outputDirAndRoot) {
		String[] samples = null;
		String[] markers = null;

		boolean sampMiss = false;
		boolean markMiss = false;
		if (fileOfSamplesToExport != null && !"".equals(fileOfSamplesToExport)) {
			if (!Files.exists(fileOfSamplesToExport)) {
				sampMiss = true;
				proj.getLog().reportFileNotFound(fileOfSamplesToExport);
			} else {
				samples = HashVec.loadFileToStringArray(fileOfSamplesToExport, false, new int[] {0}, false);
			}
		}
		if (fileOfMarkersToExport != null && !"".equals(fileOfMarkersToExport)) {
			if (!Files.exists(fileOfMarkersToExport)) {
				markMiss = true;
				proj.getLog().reportFileNotFound(fileOfMarkersToExport);
			} else {
				markers = HashVec.loadFileToStringArray(fileOfMarkersToExport, false, new int[] {0}, false);
			}
		}
		if (sampMiss || markMiss) {
			return;
		}

		exportGenvisisToVCF(proj, samples, markers, splitChrs, useGRCRefGen, chrsToExport,
												outputDirAndRoot);
	}

	public static void exportGenvisisToVCF(Project proj, String[] samplesToExport,
																				 String[] markersToExport,
																				 boolean splitChrs, boolean useGRCContigs,
																				 int[] chrsToExport,
																				 String outputDirAndRoot) {
		SampleData sd = proj.getSampleData(false);
		String[] allSamples = proj.getSamples();
		Set<String> idsSet = new HashSet<String>();
		List<String> idsToInclude = new ArrayList<String>();
		for (String s : (samplesToExport == null || samplesToExport.length == 0 ? allSamples
																																						: samplesToExport)) {
			idsSet.add(s);
		}
		Map<String, Integer> idIndexMap = new HashMap<String, Integer>();
		for (int i = 0; i < allSamples.length; i++) {
			String s = allSamples[i];
			if (sd.lookup(s) == null) {
				System.out.println("No SampleData lookup for " + s);
			} else {
				if (idsSet.contains(s) || idsSet.contains(sd.lookup(s)[1])) {
					idsToInclude.add(s);
					idIndexMap.put(s, i);
				}
			}
		}

		// UNUSED - could potentially apply
		// String clusterFilterFileName = proj.CLUSTER_FILTER_COLLECTION_FILENAME.getValue();
		ClusterFilterCollection clusterFilterCollection = null;
		// if (clusterFilterFileName == null || !Files.exists(clusterFilterFileName)) {
		// clusterFilterCollection = null;
		// } else {
		// clusterFilterCollection = ClusterFilterCollection.load(clusterFilterFileName,
		// proj.JAR_STATUS.getValue());
		// }

		ReferenceGenome refGen = useGRCContigs
																					 ? new ReferenceGenome(
																																 Resources.genome(proj.GENOME_BUILD_VERSION.getValue(),
																																									proj.getLog())
																																					.getGRCFASTA()
																																					.getAbsolute(),
																																 proj.getLog())
																					 : new ReferenceGenome(
																																 Resources.genome(proj.GENOME_BUILD_VERSION.getValue(),
																																									proj.getLog())
																																					.getFASTA()
																																					.getAbsolute(),
																																 proj.getLog());

		if (refGen.getIndexedFastaSequenceFile().getSequenceDictionary() == null) {
			proj.getLog()
					.reportError("VCF export requires a valid ReferenceGenome. Please fix your project or configuration and try again.");
			return;
		}

		MarkerDetailSet mds = proj.getMarkerSet();

		HashSet<Integer> allowedChrs = new HashSet<>();
		if (chrsToExport != null && chrsToExport.length != 0) {
			for (int c : chrsToExport) {
				allowedChrs.add(c);
			}
		}

		Map<String, Marker> markerMap = mds.getMarkerNameMap();

		HashMap<Integer, List<String>> markersByChr = new HashMap<>();
		for (String m : (markersToExport == null ? proj.getMarkerNames()
																						 : markersToExport)) {
			Marker mkr = markerMap.get(m);
			if (allowedChrs.isEmpty() || allowedChrs.contains((int) mkr.getChr())) {
				List<String> mkrChr = markersByChr.get((int) mkr.getChr());
				if (mkrChr == null) {
					mkrChr = new ArrayList<>();
					markersByChr.put((int) mkr.getChr(), mkrChr);
				}
				mkrChr.add(m);
			}
		}

		ConcurrentMap<Marker, MATCH_FAIL> failMap = new ConcurrentHashMap<>();

		ArrayList<ActualExporter> runners = new ArrayList<>();
		if (splitChrs) {
			for (int chr : markersByChr.keySet()) {
				String chrName = Positions.chromosomeNumberInverse(chr);
				proj.getLog()
						.report("Exporting " + markersByChr.get(chr).size() + " markers for chr" + chrName);
				List<String> mkrs = getMarkersSorted(markersByChr.get(chr), markerMap);

				String fileOut = outputDirAndRoot + "_chr" + chrName + ".vcf.gz";

				ActualExporter runner = new ActualExporter(proj, refGen, !useGRCContigs, fileOut,
																									 idsToInclude, idIndexMap,
																									 mkrs.toArray(new String[mkrs.size()]),
																									 markerMap,
																									 clusterFilterCollection, failMap);
				runners.add(runner);
			}
		} else {
			proj.getLog()
					.report("Exporting "
									+ (markersToExport == null ? proj.getMarkerNames() : markersToExport).length
									+ " markers");
			List<String> allMarkers = new ArrayList<>();
			Integer[] chrs = markersByChr.keySet().toArray(new Integer[markersByChr.keySet().size()]);
			Arrays.sort(chrs);
			for (int i = 0; i < chrs.length; i++) {
				allMarkers.addAll(getMarkersSorted(markersByChr.get(chrs[i]), markerMap));
			}

			String fileOut = outputDirAndRoot + ".vcf.gz";

			ActualExporter runner = new ActualExporter(proj, refGen, !useGRCContigs, fileOut,
																								 idsToInclude,
																								 idIndexMap,
																								 allMarkers.toArray(new String[allMarkers.size()]),
																								 markerMap, clusterFilterCollection, failMap);
			runners.add(runner);
		}

		if (runners.size() == 1) {
			runners.get(0).run();
		} else {
			boolean b = execute(runners);
			if (!b) {
				System.err.println("Error - Failed to complete properly.  Please check output for errors or missing data/files.");
			}
		}

		int cntS, cntA, cntSA, missSeq, missAll;
		cntS = cntA = cntSA = missSeq = missAll = 0;

		for (MATCH_FAIL mf : failMap.values()) {
			switch (mf) {
				case POSS_ALLELE:
					cntA++;
					break;
				case POSS_STRAND_ALLELE:
					cntSA++;
					break;
				case STRAND_FLIP:
					cntS++;
					break;
				case MISSING_SEQ:
					missSeq++;
					break;
				case MISSING_ALLELES:
					missAll++;
					break;
				default:
					break;
			}
		}

		if (cntS > 0) {
			proj.getLog().reportTimeWarning("Found " + cntS + " strand flipped markers.");
		}

		if (cntA > 0) {
			proj.getLog()
					.reportTimeWarning("Found "
														 + cntA
														 + " possibly allele-flipped markers.  ALT allele is unknown, so these cannot be confirmed.");
		}

		if (cntSA > 0) {
			proj.getLog()
					.reportTimeWarning("Found "
														 + cntSA
														 + " possibly strand- and allele-flipped markers.  ALT allele is unknown, so these cannot be confirmed.");
		}
		if (missSeq > 0) {
			proj.getLog().reportTimeWarning("Found " + missSeq
																			+ " markers missing from the reference genome.");
		}
		if (missAll > 0) {
			proj.getLog()
					.reportTimeWarning("Found "
														 + missAll
														 + " markers missing both alleles in the project data. These markers were dropped.");
		}

	}

	private static boolean execute(List<ActualExporter> runners) {
		int avail = Runtime.getRuntime().availableProcessors();
		ExecutorService executor = Executors.newFixedThreadPool(Math.max(1, avail - 1));

		CountDownLatch cdl = new CountDownLatch(runners.size());
		for (ActualExporter runner : runners) {
			runner.setCDL(cdl);
			runner.setParentThread(Thread.currentThread());
			executor.submit(runner);
		}

		executor.shutdown();
		boolean success = false;
		try {
			cdl.await(Long.MAX_VALUE, TimeUnit.MILLISECONDS);
			success = true;
		} catch (InterruptedException e) {
		}
		return success;
	}

	private static List<String> getMarkersSorted(List<String> chrMarkers,
																							 Map<String, Marker> markerMap) {
		int[] pos = new int[chrMarkers.size()];
		String[] mkr = chrMarkers.toArray(new String[chrMarkers.size()]);
		for (int i = 0; i < mkr.length; i++) {
			pos[i] = markerMap.get(mkr[i]).getPosition();
		}

		int[] indices = Sort.getSortedIndices(pos);

		ArrayList<String> mkrs = new ArrayList<String>();
		for (int i = 0; i < indices.length; i++) {
			// skip if prev (in sorted array) was same position:
			if (i == 0 || pos[indices[i]] != pos[indices[i - 1]]) {
				mkrs.add(mkr[indices[i]]);
			}
		}
		return mkrs;
	}


	private static boolean matches(Allele ref, Allele mkr) {
		if (ref != null) {
			if (ref.basesMatch(mkr)) {
				return true;
			}
		}
		return false;
	}

	enum MATCH_FAIL {
		MISSING_ALLELES,
		MISSING_SEQ,
		STRAND_FLIP,
		POSS_ALLELE,
		POSS_STRAND_ALLELE;
	}

	private static class ActualExporter implements Runnable {

		Project proj;
		ReferenceGenome refGen;
		String fileOut;
		List<String> idsToInclude;
		Map<String, Integer> idIndexMap;
		String[] markers;
		Map<String, Marker> markerMap;
		ClusterFilterCollection clusterFilterCollection;
		boolean useChr;

		CountDownLatch cdl;
		Thread parent;

		public ActualExporter(Project proj, ReferenceGenome refGen, boolean useChr,
													String fileOut, List<String> idsToInclude,
													Map<String, Integer> idIndexMap,
													String[] markers, Map<String, Marker> markerMap,
													ClusterFilterCollection clusterFilterCollection,
													ConcurrentMap<Marker, MATCH_FAIL> failMap) {
			this.proj = proj;
			this.refGen = refGen;
			this.fileOut = fileOut;
			this.idsToInclude = idsToInclude;
			this.idIndexMap = idIndexMap;
			this.markers = markers;
			this.markerMap = markerMap;
			this.clusterFilterCollection = clusterFilterCollection;
			this.fails = failMap;
			this.useChr = useChr;
		}

		protected void setCDL(CountDownLatch cdl) {
			this.cdl = cdl;
		}

		protected void setParentThread(Thread th) {
			this.parent = th;
		}

		@Override
		public void run() {
			try {
				actualRun();
			} catch (Exception e) {
				e.printStackTrace();
				if (parent != null) {
					parent.interrupt();
				}
			}
			if (cdl != null) {
				cdl.countDown();
			}
		}

		private void actualRun() {
			Logger log = proj.getLog();

			VariantContextWriterBuilder builder = new VariantContextWriterBuilder().setOutputFile(fileOut);
			builder.clearOptions();
			builder.setOption(Options.INDEX_ON_THE_FLY);

			HashSet<VCFHeaderLine> lines = new HashSet<VCFHeaderLine>();
			VCFFormatHeaderLine format = new VCFFormatHeaderLine("GT", 1, VCFHeaderLineType.String, "GT");
			lines.add(format);

			VCFHeader vcfHeader = new VCFHeader(lines, idsToInclude);

			SAMSequenceDictionary samSequenceDictionary = refGen.getIndexedFastaSequenceFile()
																													.getSequenceDictionary();

			builder.setReferenceDictionary(samSequenceDictionary);
			vcfHeader.setSequenceDictionary(samSequenceDictionary);

			VariantContextWriter writer = builder.build();
			vcfHeader.hasGenotypingData();
			writer.writeHeader(vcfHeader);

			float gcThreshold = proj.GC_THRESHOLD.getValue().floatValue();

			log.report("Initializing MDL for " + markers.length + " markers...");
			MDL mdl = new MDL(proj, proj.getMarkerSet(), markers);

			while (mdl.hasNext()) {
				MarkerData markerData = mdl.next();

				VariantContextBuilder builderVc = new VariantContextBuilder();
				builderVc.chr((useChr ? "chr" : "")
											+ Positions.chromosomeNumberInverse(markerData.getChr()));
				Marker mkr = markerMap.get(markerData.getMarkerName());
				ArrayList<Allele> a = new ArrayList<Allele>();
				Allele aR = mkr.getRef();
				Allele aA = mkr.getAlt();

				String[] bases = refGen.getSequenceFor(new Segment(mkr.getChr(),
																													 mkr.getPosition(),
																													 mkr.getPosition()));
				if (aR == null && aA == null) {
					missingAlleles(mkr);
					continue;
				}
				if (aR == null) {
					log.reportError("REF allele for marker " + mkr.getName() + " is null!");
					aR = bases == null ? Allele.create("N") : Allele.create(bases[0], true);
				}
				if (aA == null) {
					log.reportError("ALT allele for marker " + mkr.getName() + " is null!");
					aA = Allele.create("N", false);
				}
				a.add(aR);
				a.add(aA);

				if (bases == null || (bases.length == 1 && bases[0].trim().isEmpty())) {
					missingSeq(mkr);
				} else {
					Allele refGenRef = Allele.create(bases[0], true);
					boolean match = matches(refGenRef, aR);
					if (!match) {
						// Cannot fully determine if alleles are mismatched due to lacking ALT allele
						if (matches(refGenRef, aA)) {
							// Possible (but not certain) allele flip
							possibleAlleleFlips(mkr);
						} else {
							Allele mkrRefFlip = Allele.create(StrandOps.flipsIfNeeded(mkr.getRef()
																																					 .getBaseString(),
																																				Strand.NEGATIVE,
																																				true));
							if (matches(refGenRef, mkrRefFlip)) {
								// strand flip
								strandFlip(mkr);
							} else {
								Allele mkrAltFlip = Allele.create(StrandOps.flipsIfNeeded(mkr.getAlt()
																																						 .getBaseString(),
																																					Strand.NEGATIVE,
																																					true));
								if (matches(refGenRef, mkrAltFlip)) {
									// Possible (but not certain) allele and strand flip
									possibleStrandAndAlleleFlip(mkr);
								}
							}
						}
					}
				}

				builderVc.alleles(a);
				builderVc.start(mkr.getPosition());
				builderVc.stop(mkr.getPosition() + mkr.getRef().length() - 1);
				builderVc.id(mkr.getName());
				Collection<Genotype> genos = new ArrayList<Genotype>();
				byte[] genotypes = markerData.getAbGenotypesAfterFilters(clusterFilterCollection,
																																 markerData.getMarkerName(),
																																 gcThreshold, proj.getLog());
				for (int k = 0; k < idsToInclude.size(); k++) {
					int idInd = idIndexMap.get(idsToInclude.get(k));
					String id = idsToInclude.get(k);
					Genotype g;
					List<Allele> all;
					byte indGeno = genotypes[idInd];
					// 0 for A/A, 1 for A/B, 2 for B/B, and -1 for null
					boolean aIsRef = mkr.getRefAllele().equals(RefAllele.A);
					switch (indGeno) {
						case 0:
							all = aIsRef ? Arrays.asList(aR, aR) : Arrays.asList(aA, aA);
							break;
						case 1:
							all = aIsRef ? Arrays.asList(aR, aA) : Arrays.asList(aR, aA);
							break;
						case 2:
							all = aIsRef ? Arrays.asList(aA, aA) : Arrays.asList(aR, aR);
							break;
						case -1:
						default:
							all = Arrays.asList(Allele.NO_CALL, Allele.NO_CALL);
							break;
					}
					g = GenotypeBuilder.create(id, all);

					genos.add(g);
				}

				builderVc.genotypes(genos);
				writer.add(builderVc.make());

			}
			mdl.shutdown();
			mdl = null;

			writer.close();
			writer = null;

			System.gc();

			log.report("Processed " + markers.length + " markers for " + idsToInclude.size()
								 + " samples.");
		}

		Map<Marker, MATCH_FAIL> fails;

		private void missingAlleles(Marker mkr) {
			if (fails != null) {
				fails.put(mkr, MATCH_FAIL.MISSING_ALLELES);
			}
		}

		private void missingSeq(Marker mkr) {
			if (fails != null) {
				fails.put(mkr, MATCH_FAIL.MISSING_SEQ);
			}
		}

		private void possibleStrandAndAlleleFlip(Marker mkr) {
			if (fails != null) {
				fails.put(mkr, MATCH_FAIL.POSS_STRAND_ALLELE);
			}
		}

		private void strandFlip(Marker mkr) {
			if (fails != null) {
				fails.put(mkr, MATCH_FAIL.STRAND_FLIP);
			}
		}

		private void possibleAlleleFlips(Marker mkr) {
			if (fails != null) {
				fails.put(mkr, MATCH_FAIL.POSS_ALLELE);
			}
		}
	}

	private static final String SAMP_ARG = "samples";
	private static final String MARK_ARG = "markers";
	private static final String SPLIT_ARG = "split";
	private static final String GRC_ARG = "grc";
	private static final String CHRS_ARG = "chrs";

	public static void main(String[] args) {
		Project proj = null;
		String samp = null;
		String mark = null;
		boolean split = true;
		int[] chrs = null;
		boolean grc = true;
		String out = null;

		Object[][] argSet = {
												 {CLI.ARG_PROJ, CLI.DESC_PROJ, "example.properties", Arg.STRING},
												 {SAMP_ARG, "List of Sample IDs", "", Arg.STRING},
												 {MARK_ARG, "List of Markers", "", Arg.STRING},
												 {SPLIT_ARG, "Split output by chromosomes", split, Arg.STRING},
												 {GRC_ARG,
													"Use GRC output rather than HG (doesn't export 'chr' in contigs)", grc,
													Arg.STRING},
												 {CHRS_ARG,
													"OPTIONAL: comma-delimited list of specific chromosomes to export", "",
													Arg.STRING},
												 {CLI.ARG_OUTFILE, CLI.DESC_OUTFILE, "", Arg.STRING},
		};

		CLI cli = new CLI(VCFData.class);

		for (Object[] arg : argSet) {
			cli.addArgWithDefault((String) arg[0], (String) arg[1], arg[2] + "", (Arg) arg[3]);
		}

		cli.parseWithExit(args);

		String projFile = cli.get(CLI.ARG_PROJ);
		if ("".equals(projFile) || !Files.exists(projFile)) {
			System.err.println("Error - valid project properties file is required!");
			System.err.println("");
			return;
		}
		proj = new Project(projFile);
		samp = cli.get(SAMP_ARG);
		mark = cli.get(MARK_ARG);
		split = Boolean.parseBoolean(cli.get(SPLIT_ARG));
		grc = Boolean.parseBoolean(cli.get(GRC_ARG));
		String[] chrStrs = cli.get(CHRS_ARG).equals("") ? null : cli.get(CHRS_ARG).split(",");
		chrs = new int[chrStrs.length];
		for (int i = 0; i < chrs.length; i++) {
			chrs[i] = Integer.parseInt(chrStrs[i]);
		}
		out = cli.get(CLI.ARG_OUTFILE);
		if (out == null || "".equals(out)) {
			out = proj.PROJECT_NAME.getValue();
		}

		VCFData.exportGenvisisToVCF(proj, samp, mark, split, grc, chrs, out);
	}
}
