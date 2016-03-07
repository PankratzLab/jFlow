package seq.analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Set;
import java.util.concurrent.Callable;

import cnv.var.LocusSet;
import filesys.Segment;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.tribble.annotation.Strand;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;
import seq.analysis.Blast.BlastResults;
import seq.manage.Adapter;
import seq.manage.BamExtractor;
import seq.manage.CigarOps;
import seq.manage.BamExtractor.WorkerExtractor;
import seq.manage.BamOps;
import seq.manage.SamRecordOps;
import seq.manage.SamRecordOps.SoftClipped;
import seq.manage.StrandOps;
import seq.manage.VCFOps;
import seq.manage.VCOps;
import common.Array;
import common.ExcelConverter;
import common.Files;
import common.HashVec;
import common.Logger;
import common.PSF;
import common.WorkerHive;
import common.WorkerTrain;
import common.WorkerTrain.Producer;
import common.ext;

/**
 * Indel qc
 *
 */
public class Indelathon {

	@SuppressWarnings("unchecked")
	private static void summarizeSoftClippings(String vcf, String bamFilesFile, String outDir, Set<String> variantSets, int buffer, int minSCLenth, int minSCCount, int numThreads) {
		new File(outDir).mkdirs();
		Logger log = new Logger(outDir + "indel.log");
		String[] bams = HashVec.loadFileToStringArray(bamFilesFile, false, new int[] { 0 }, true);
		String[] samps = VCFOps.getSamplesInFile(vcf);
		ArrayList<String> excelFiles = new ArrayList<String>();

		Hashtable<String, String> matchedSamps = BamOps.matchToVcfSamplesToBamFiles(samps, variantSets, bams, numThreads, log);
		String outIndelVCF = outDir + VCFOps.getAppropriateRoot(vcf, true) + ".indels.vcf.gz";
		String outSegSer = outDir + VCFOps.getAppropriateRoot(vcf, true) + ".indels.seg.ser";
		Hashtable<String, ArrayList<Segment>> sampSegs = new Hashtable<String, ArrayList<Segment>>();

		if (!Files.exists(outIndelVCF) || !Files.exists(outSegSer)) {
			extractIndelVariants(vcf, log, outIndelVCF, outSegSer, sampSegs);
		}
		sampSegs = (Hashtable<String, ArrayList<Segment>>) Files.readSerial(outSegSer, false, log, false, true);
		String indelBamDir = outDir + "indel_bams/";
		String[] indelBams = extractIndelSegments(indelBamDir, matchedSamps, sampSegs, buffer, numThreads, log);
		Hashtable<String, String> matchedIndelSamps = BamOps.matchToVcfSamplesToBamFiles(samps, variantSets, indelBams, numThreads, log);

		// extracting soft clips
		SoftClipResultProducer producer = new SoftClipResultProducer(samps, sampSegs, matchedSamps, indelBamDir, log);
		WorkerTrain<SoftClipResult> train = new WorkerTrain<Indelathon.SoftClipResult>(producer, numThreads, 10, log);
		ArrayList<SoftClipResult> results = new ArrayList<Indelathon.SoftClipResult>();
		while (train.hasNext()) {
			results.add(train.next());
		}

		// summarize barcodes
		String outBarCode = outDir + "barCodes.txt";
		String bcPrint = "bamSample\tbam\tBarcode1\tBarcode2";
		for (int i = 0; i < bams.length; i++) {
			ArrayList<String> bc = BamOps.getBarcodesFor(bams[i], log);
			bcPrint += "\n" + BamOps.getSampleName(bams[i]) + "\t" + bams[i] + "\t" + Array.toStr(bc.toArray(new String[bc.size()]));
		}
		Files.write(bcPrint, outBarCode);
		String out = outDir + "countit.txt";
		HashSet<String> allClips = new HashSet<String>();
		ArrayList<Adapter> adapters = Adapter.getCurrentAdapters(Adapter.getBarcodes(BamOps.getAllBarCodes(bams, log)));
		for (SoftClipResult result : results) {
			for (String clip : result.getScAllCounts().keySet()) {
				if (clip.replaceAll("N", "").length() >= minSCLenth && result.getScAllCounts().get(clip) >= minSCCount) {
					allClips.add(clip);
				}
			}
		}
		String blastDir = outDir + "blast/";
		String[] blasts = Adapter.blast(minSCLenth, adapters, allClips.toArray(new String[allClips.size()]), blastDir, "softClip", numThreads, log);
		Hashtable<String, HitSummary> hits = summarizeBlasts(blasts, allClips, log);
		try {

			PrintWriter writer = new PrintWriter(new FileWriter(out));
			writer.print("Clip\t" + HitSummary.getHeader());
			for (SoftClipResult result : results) {
				if (result.getBamFile() != null) {
					writer.print("\t" + result.getVcfSample());
				}
			}
			writer.println();
			for (String clip : allClips) {
				HitSummary hitSummary = null;
				if (hits.containsKey(clip)) {
					hitSummary = hits.get(clip);
				}
				writer.print(clip + "\t" + (hitSummary == null ? HitSummary.getBlank() : hitSummary.getSummary()));
				for (SoftClipResult result : results) {
					if (result.getBamFile() != null) {
						if (result.getScAllCounts().containsKey(clip)) {
							writer.print("\t" + result.getScAllCounts().get(clip));
						} else {
							writer.print("\t0");
						}

					}
				}
				writer.println();
			}

			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + out);
			log.reportException(e);
		}

		excelFiles.add(out);
		excelFiles.add(outBarCode);
		String excelOut = outDir + "summary.xlsx";
		ExcelConverter converter = new ExcelConverter(excelFiles, excelOut, log);
		converter.convert(true);
		for (SoftClipResult result : results) {
			Hashtable<String, Integer> maxCountClipLength = new Hashtable<String, Integer>();
			for (String clip : result.getScAllCounts().keySet()) {
				String len = clip.length() + "";
				if (maxCountClipLength.containsKey(len)) {

				}
				else {
					maxCountClipLength.put(len, result.getScAllCounts().get(clip));
				}
			}
			// result.g
		}

	}

	private static class HitSummary {
		private String clip;
		private String adapterHit;
		private String barcode;
		private Cigar adapterCigar;
		private Cigar barCodeCigar;
		private double minEvaladapter;
		private double minEvalBarcode;
		private String strandClipAdapter;
		private String strandClipBarCode;

		public HitSummary(String clip) {
			super();
			this.clip = clip;
			this.minEvaladapter = Double.NaN;
			this.minEvalBarcode = Double.NaN;
			this.adapterHit = "NONE";
			this.barcode = "NONE";
			this.barCodeCigar = null;
			this.adapterCigar = null;
			this.strandClipAdapter = "NA";
			this.strandClipBarCode = "NA";
		}

		public void setStrandClipAdapter(String strandClipAdapter) {
			this.strandClipAdapter = strandClipAdapter;
		}

		public void setStrandClipBarCode(String strandClipBarCode) {
			this.strandClipBarCode = strandClipBarCode;
		}

		public void setBarCodeCigar(Cigar barCodeCigar) {
			this.barCodeCigar = barCodeCigar;
		}

		public void setAdapterCigar(Cigar adapterCigar) {
			this.adapterCigar = adapterCigar;
		}

		private String getSummary() {
			return strandClipBarCode + "\t" + barcode + "\t" + (barCodeCigar == null ? "NA" : barCodeCigar.toString()) + "\t" + minEvalBarcode + "\t" + strandClipAdapter + "\t" + adapterHit + "\t" + (adapterCigar == null ? "NA" : adapterCigar.toString()) + "\t" + minEvaladapter;
		}

		public static String getBlank() {
			return "NA\tNONE\tNONE\t" + Double.NaN + "\tNA\tNONE\tNONE\t" + Double.NaN;
		}

		public static String getHeader() {
			return "StrandClipBarCode\tBarcode\tBarcodeCigar\tMinEvalBarcode\tStrandClipAdapter\tAdapter\tAdapterCigar\tMinEvalAdapter";

		}

		public double getMinEvaladapter() {
			return minEvaladapter;
		}

		public void setMinEvaladapter(double minEvaladapter) {
			this.minEvaladapter = minEvaladapter;
		}

		public double getMinEvalBarcode() {
			return minEvalBarcode;
		}

		public void setMinEvalBarcode(double minEvalBarcode) {
			this.minEvalBarcode = minEvalBarcode;
		}

		public void setAdapterHit(String adapterHit) {
			this.adapterHit = adapterHit;
		}

		public void setBarcode(String barcode) {
			this.barcode = barcode;
		}

	}

	private static Hashtable<String, HitSummary> summarizeBlasts(String[] blastFiles, Set<String> softClip, Logger log) {
		Hashtable<String, HitSummary> hits = new Hashtable<String, HitSummary>();

		for (int i = 0; i < blastFiles.length; i++) {
			try {
				BufferedReader reader = Files.getAppropriateReader(blastFiles[i]);
				reader.readLine();
				while (reader.ready()) {
					BlastResults result = new BlastResults(reader.readLine().trim().split("\t"), log);
					String clip = result.getQueryID().split("_")[0];
					String strandClip = result.isStrandFlipped() ? StrandOps.flipsIfNeeded(clip, Strand.NEGATIVE, true, true) : clip;
					Cigar cigar = CigarOps.convertBtopToCigar(result, clip.length(), log);
					boolean barcode = result.getSubjectID().split("_")[1].equals(Adapter.BARCODE);

					if (!hits.containsKey(clip)) {
						hits.put(clip, new HitSummary(clip));
					}
					if (barcode) {
						if (Double.isNaN(hits.get(clip).getMinEvalBarcode())) {
							hits.get(clip).setMinEvalBarcode(result.getEvalue());
							hits.get(clip).setBarcode(result.getSubjectID());
							hits.get(clip).setBarCodeCigar(cigar);
							hits.get(clip).setStrandClipBarCode(strandClip);

						} else if (hits.get(clip).getMinEvalBarcode() > result.getEvalue()) {
							hits.get(clip).setMinEvalBarcode(result.getEvalue());
							hits.get(clip).setBarcode(result.getSubjectID());
							hits.get(clip).setBarCodeCigar(cigar);
							hits.get(clip).setStrandClipBarCode(strandClip);

						}
					} else {
						if (Double.isNaN(hits.get(clip).getMinEvaladapter())) {
							hits.get(clip).setAdapterCigar(cigar);
							hits.get(clip).setMinEvaladapter(result.getEvalue());
							hits.get(clip).setAdapterHit(result.getSubjectID());
							hits.get(clip).setStrandClipAdapter(strandClip);

						} else if (hits.get(clip).getMinEvaladapter() > result.getEvalue()) {
							hits.get(clip).setMinEvaladapter(result.getEvalue());
							hits.get(clip).setAdapterHit(result.getSubjectID());
							hits.get(clip).setAdapterCigar(cigar);
							hits.get(clip).setStrandClipAdapter(strandClip);

						}
					}
				}
				reader.close();
			} catch (FileNotFoundException fnfe) {
				log.reportError("Error: file \"" + blastFiles[i] + "\" not found in current directory");
				return null;
			} catch (IOException ioe) {
				log.reportError("Error reading file \"" + blastFiles[i] + "\"");
				return null;
			}
		}
		return hits;

	}

	/**
	 * Stores the soft clip summary from a bam file...ideally one subsetted to indels
	 *
	 */
	private static class SoftClipResult implements Callable<SoftClipResult> {
		private LocusSet<Segment> locs;
		private Hashtable<String, Integer> scAllCounts;
		private Hashtable<String, Integer> scSegCounts;
		private Hashtable<String, ArrayList<String>> segScs;
		private String scAllCountSerFile;
		private String scSegCountSerFile;
		private String segScsSerFile;
		private String vcfSample;
		private String bamFile;
		private Logger log;

		public SoftClipResult(ArrayList<Segment> toQuery, String vcfSample, String bamFile, String outputDir, Logger log) {
			super();
			this.locs = new LocusSet<Segment>(toQuery.toArray(new Segment[toQuery.size()]), true, log) {

				/**
				 * 
				 */
				private static final long serialVersionUID = 1L;
			};
			this.vcfSample = vcfSample;
			this.bamFile = bamFile;
			this.scAllCountSerFile = bamFile == null ? null : outputDir + ext.rootOf(bamFile, true) + ".softclippedAllCounts.ser";
			this.scSegCountSerFile = bamFile == null ? null : outputDir + ext.rootOf(bamFile, true) + ".softclippedSegCounts.ser";
			this.segScsSerFile = bamFile == null ? null : outputDir + ext.rootOf(bamFile, true) + ".softclippedSegmentCounts.ser";
			this.log = log;
		}

		public String getBamFile() {
			return bamFile;
		}

		public String getVcfSample() {
			return vcfSample;
		}

		public Hashtable<String, Integer> getScAllCounts() {
			return scAllCounts;
		}

		@SuppressWarnings("unchecked")
		@Override
		public SoftClipResult call() throws Exception {
			this.scAllCounts = new Hashtable<String, Integer>();
			this.scSegCounts = new Hashtable<String, Integer>();

			this.segScs = new Hashtable<String, ArrayList<String>>();

			if (bamFile != null) {
				if (!Files.exists(scAllCountSerFile) || !Files.exists(segScsSerFile)) {
					SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault();
					samReaderFactory.validationStringency(ValidationStringency.LENIENT);
					SamReader reader = samReaderFactory.open(new File(bamFile));
					int num = 0;
					for (SAMRecord samRecord : reader) {
						num++;
						if (num % 100000 == 0) {
							log.reportTimeInfo("Scanned " + num + " reads, currently on " + SamRecordOps.getDisplayLoc(samRecord) + "from " + bamFile);
							log.reportTimeInfo("Found " + scAllCounts.keySet().size() + " unique soft clipped reads");
						}
						if (!samRecord.getReadUnmappedFlag() && !samRecord.getDuplicateReadFlag()) {
							SoftClipped[] softs = SamRecordOps.getSoftClippedBases(samRecord, log);
							if (softs.length > 0) {
								for (int i = 0; i < softs.length; i++) {
									if (!scAllCounts.containsKey(softs[i].getBases())) {
										scAllCounts.put(softs[i].getBases(), 0);
									}
									scAllCounts.put(softs[i].getBases(), scAllCounts.get(softs[i].getBases()) + 1);
									Segment[] varLocs = locs.getOverLappingLoci(softs[i].getRefSeg());
									if (varLocs != null) {// soft clipped contained in indel call
										if (!scSegCounts.containsKey(softs[i].getBases())) {
											scSegCounts.put(softs[i].getBases(), 0);
										}
										scSegCounts.put(softs[i].getBases(), scSegCounts.get(softs[i].getBases()) + 1);
										for (int j = 0; j < varLocs.length; j++) {
											String key = varLocs[j].getUCSClocation();
											if (!segScs.containsKey(key)) {
												segScs.put(key, new ArrayList<String>());
											}

											segScs.get(key).add(softs[i].getBases());
										}
									}
								}
							}
						}
					}
					Files.writeSerial(scAllCounts, scAllCountSerFile, true);
					Files.writeSerial(scSegCounts, scSegCountSerFile, true);
					Files.writeSerial(segScs, segScsSerFile, true);
				}

				// }
				scAllCounts = (Hashtable<String, Integer>) Files.readSerial(scAllCountSerFile, false, log, false, true);
				segScs = (Hashtable<String, ArrayList<String>>) Files.readSerial(segScsSerFile, false, log, false, true);

			} else {
				log.reportTimeWarning("Skipping sample " + vcfSample + " could not determine bam file");
			}
			return this;
		}

	}

	private static class SoftClipResultProducer implements Producer<SoftClipResult> {
		private String[] samples;
		private Hashtable<String, ArrayList<Segment>> sampSegs;
		private Hashtable<String, String> matchedIndelSamps;
		private String outputDir;
		private Logger log;
		private int index;

		public SoftClipResultProducer(String[] samples, Hashtable<String, ArrayList<Segment>> sampSegs, Hashtable<String, String> matchedSamps, String outputDir, Logger log) {
			super();
			this.samples = samples;
			this.sampSegs = sampSegs;
			this.matchedIndelSamps = matchedSamps;
			this.outputDir = outputDir;
			this.log = log;
			this.index = 0;
		}

		@Override
		public boolean hasNext() {
			return index < samples.length;
		}

		@Override
		public Callable<SoftClipResult> next() {
			String vcfSample = samples[index];
			ArrayList<Segment> toQuery = new ArrayList<Segment>();
			String bamFile = null;
			if (matchedIndelSamps.containsKey(vcfSample)) {
				bamFile = matchedIndelSamps.get(vcfSample);
				toQuery = sampSegs.get(vcfSample);
			}
			SoftClipResult sc = new SoftClipResult(toQuery, vcfSample, bamFile, outputDir, log);
			index++;
			return sc;
		}

		@Override
		public void shutdown() {
			// TODO Auto-generated method stub

		}

	}

	private static void extractIndelVariants(String vcf, Logger log, String outIndelVCF, String outSegSer, Hashtable<String, ArrayList<Segment>> sampSegs) {
		VCFFileReader reader = new VCFFileReader(vcf, true);
		VariantContextWriter writer = VCFOps.initWriterWithHeader(reader, outIndelVCF, VCFOps.DEFUALT_WRITER_OPTIONS, log);
		int numTotal = 0;
		int numIndels = 0;
		for (VariantContext vc : reader) {
			numTotal++;
			if (vc.isIndel()) {
				GenotypesContext gc = vc.getGenotypes();
				for (Genotype g : gc) {
					if (!g.isNoCall() && !g.isHomRef()) {
						if (!sampSegs.containsKey(g.getSampleName())) {
							sampSegs.put(g.getSampleName(), new ArrayList<Segment>());
						}
						sampSegs.get(g.getSampleName()).add(VCOps.getSegment(vc));
					}
				}
				numIndels++;
				writer.add(vc);
			}
		}
		log.reportTimeInfo("Found " + numIndels + " indels for " + numTotal + " variants total");
		writer.close();
		reader.close();
		String[] allSamps = VCFOps.getSamplesInFile(vcf);
		for (int i = 0; i < allSamps.length; i++) {
			if (!sampSegs.containsKey(allSamps[i])) {
				sampSegs.put(allSamps[i], new ArrayList<Segment>());
			}
		}
		Files.writeSerial(sampSegs, outSegSer, true);
	}

	private static String[] extractIndelSegments(String indelBamDir, Hashtable<String, String> matchedSamps, Hashtable<String, ArrayList<Segment>> sampSegs, int buffer, int numThreads, Logger log) {
		WorkerHive<BamExtractor> hive = new WorkerHive<BamExtractor>(numThreads, 10, log);
		ArrayList<String> indelBams = new ArrayList<String>();
		new File(indelBamDir).mkdirs();
		for (String samp : matchedSamps.keySet()) {
			String outbam = indelBamDir + samp + "_indels_" + buffer + "bp.bam";
			indelBams.add(outbam);
			Segment[] segmentsToExtract = new Segment[] {};
			if (sampSegs.containsKey(samp)) {
				segmentsToExtract = sampSegs.get(samp).toArray(new Segment[sampSegs.get(samp).size()]);
			}
			WorkerExtractor extractor = new WorkerExtractor(segmentsToExtract, matchedSamps.get(samp), false, false, buffer, outbam, log);
			hive.addCallable(extractor);
		}
		hive.execute(true);
		return Array.toStringArray(indelBams);
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String vcf = "Indelathon.vcf";
		String outDir = "Indelathon/";
		String bams = null;
		HashSet<String> sets = new HashSet<String>();
		int numThreads = 24;
		int buffer = 100;
		int minSC = 0;
		int minSCCount = 5;

		String usage = "\n" +
				"seq.analysis.Indelathon requires 0-1 arguments\n" +
				"   (1) filename (i.e. vcf=" + vcf + " (default))\n" +
				"   (2) outDir (i.e. out=" + outDir + " (default))\n" +
				"   (3) bams (i.e. bams=" + bams + " (default))\n" +
				"   (4) comma delimited variant sets (i.e. sets= (default))\n" +
				"   (5) min length for soft clipped sequences  (i.e. minSC=" + minSC + " (default))\n" +
				"   (5) min number of occurances for soft clipped sequences per sample (i.e. minSCCount=" + minSCCount + " (default))\n" +

				PSF.Ext.getNumThreadsCommand(4, numThreads) +

				"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("vcf=")) {
				vcf = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("out=")) {
				outDir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("bams=")) {
				bams = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("minSC=")) {
				minSC = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("minSCCount=")) {
				minSCCount = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("sets=")) {
				String[] tmp = ext.parseStringArg(args[i], "").split(",");
				for (int j = 0; j < tmp.length; j++) {
					sets.add(tmp[j]);
				}
				numArgs--;
			} else if (args[i].startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
				numThreads = ext.parseIntArg(args[i]);
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

			summarizeSoftClippings(vcf, bams, outDir, sets, buffer, minSC, minSCCount, numThreads);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
