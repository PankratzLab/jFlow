package org.genvisis.cnv.manage;

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
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.genvisis.cnv.filesys.ClusterFilterCollection;
import org.genvisis.cnv.filesys.MarkerData;
import org.genvisis.cnv.filesys.MarkerDetailSet;
import org.genvisis.cnv.filesys.MarkerDetailSet.Marker;
import org.genvisis.cnv.filesys.MarkerDetailSet.Marker.RefAllele;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.Logger;
import org.genvisis.common.Sort;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.manage.ReferenceGenome;
import org.genvisis.seq.manage.StrandOps;

public final class VCFData {

	private VCFData() {
		// private un-instantiable class
	}

	public static void exportGenvisisToVCF(Project proj, String[] samplesToExport,
																				 String[] markersToExport,
																				 boolean splitChrs, boolean useGRCRefGen,
																				 String outputDirAndRoot) {
		SampleData sd = proj.getSampleData(false);
		String[] allSamples = proj.getSamples();
		Set<String> idsSet = new HashSet<String>();
		List<String> idsToInclude = new ArrayList<String>();
		for (String s : (samplesToExport == null ? allSamples : samplesToExport)) {
			idsSet.add(s);
		}
		Map<String, Integer> idIndexMap = new HashMap<String, Integer>();
		for (int i = 0; i < allSamples.length; i++) {
			String s = allSamples[i];
			if (sd.lookup(s) == null) {
				System.out.println("No SampleData lookup for " + s);
			} else {
				if (idsSet.contains(sd.lookup(s)[1])) {
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

		ReferenceGenome refGen = useGRCRefGen
																				 ?
																				 new ReferenceGenome(
																														 Resources.genome(proj.GENOME_BUILD_VERSION.getValue(),
																																							proj.getLog())
																																			.getGRCFASTA().getAbsolute(),
																														 proj.getLog())
																				 : new ReferenceGenome(
																															 Resources.genome(proj.GENOME_BUILD_VERSION.getValue(),
																																								proj.getLog())
																																				.getFASTA().getAbsolute(),
																															 proj.getLog());

		MarkerDetailSet mds = proj.getMarkerSet();

		Map<String, Marker> markerMap = mds.getMarkerNameMap();

		HashMap<Integer, List<String>> markersByChr = new HashMap<>();
		for (String m : markersToExport) {
			Marker mkr = markerMap.get(m);
			List<String> mkrChr = markersByChr.get((int) mkr.getChr());
			if (mkrChr == null) {
				mkrChr = new ArrayList<>();
				markersByChr.put((int) mkr.getChr(), mkrChr);
			}
			mkrChr.add(m);
		}

		ConcurrentMap<Marker, MATCH_FAIL> failMap = new ConcurrentHashMap<>();

		ArrayList<Runnable> runners = new ArrayList<>();
		if (splitChrs) {
			for (int chr : markersByChr.keySet()) {
				proj.getLog()
						.report("Exporting " + markersByChr.get(chr).size() + " markers for chr" + chr);
				List<String> mkrs = getMarkersSorted(markersByChr.get(chr), markerMap);

				String fileOut = outputDirAndRoot + "_chr" + chr + ".vcf.gz";

				Runnable runner = new ActualExporter(proj, refGen, fileOut, idsToInclude, idIndexMap,
																						 mkrs.toArray(new String[mkrs.size()]), markerMap,
																						 clusterFilterCollection, failMap);
				runners.add(runner);
			}
		} else {
			proj.getLog().report("Exporting " + markersToExport.length + " markers");
			List<String> allMarkers = new ArrayList<>();
			Integer[] chrs = markersByChr.keySet().toArray(new Integer[markersByChr.keySet().size()]);
			Arrays.sort(chrs);
			for (int i = 0; i < chrs.length; i++) {
				allMarkers.addAll(getMarkersSorted(markersByChr.get(chrs[i]), markerMap));
			}

			String fileOut = outputDirAndRoot + ".vcf.gz";

			Runnable runner = new ActualExporter(proj, refGen, fileOut, idsToInclude, idIndexMap,
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

		int cntS, cntA, cntSA;
		cntS = cntA = cntSA = 0;

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

	}

	private static boolean execute(List<Runnable> runners) {
		int avail = Runtime.getRuntime().availableProcessors();
		ExecutorService executor = Executors.newFixedThreadPool(avail - 1);

		for (Runnable runner : runners) {
			executor.submit(runner);
		}

		executor.shutdown();
		boolean success = false;
		try {
			boolean s = executor.awaitTermination(Long.MAX_VALUE, TimeUnit.MILLISECONDS);
			success = s;
		} catch (InterruptedException e) {
			e.printStackTrace();
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

		public ActualExporter(Project proj, ReferenceGenome refGen, String fileOut,
													List<String> idsToInclude,
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
		}

		@Override
		public void run() {
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

			float gcThreshold = 0; // include everything

			MDL mdl = new MDL(proj, proj.getMarkerSet(), markers);

			while (mdl.hasNext()) {
				MarkerData markerData = mdl.next();

				VariantContextBuilder builderVc = new VariantContextBuilder();
				builderVc.chr(String.valueOf((int) markerData.getChr()));
				Marker mkr = markerMap.get(markerData.getMarkerName());
				ArrayList<Allele> a = new ArrayList<Allele>();
				Allele aR = mkr.getRef();
				Allele aA = mkr.getAlt();
				a.add(aR);
				a.add(aA);

				String[] bases = refGen.getSequenceFor(new Segment(mkr.getChr(),
																													 mkr.getPosition(),
																													 mkr.getPosition()));
				if (bases == null) {
					log.reportError("Couldn't find sequence for marker in ReferenceGenome");
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
				builderVc.stop(mkr.getPosition());
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
		}

		Map<Marker, MATCH_FAIL> fails;

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

}
