package org.genvisis.seq.analysis;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Vector;
import java.util.concurrent.Callable;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.WorkerHive;
import org.genvisis.common.ext;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.manage.VCFOps;
import org.genvisis.seq.manage.VCFOps.VcfPopulation;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.AsyncVariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;

/**
 * Class to perform common info queries on local and remote vcf files<br>
 * For spot-checking
 * http://browser.1000genomes.org/Homo_sapiens/Info/Index?db=core;r=1:109818030-109819030;v=rs646776
 * ;vdb=variation;vf=510141
 * http://browser.1000genomes.org/Homo_sapiens/Variation/Population?r=1:109818030-109819030;v=
 * rs646776;vdb=variation;vf=510141
 */
/**
 * @author lane0212
 *
 */
public class VcfQuery {
	/**
	 * This ftp site currently hosts indexed (.tbi) vcf (.vcf.gz) files that can searched<br>
	 * functional Annotations:
	 * ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/supporting/functional_annotation/
	 * unfiltered/ <br>
	 * ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/supporting/functional_annotation/
	 * filtered/
	 */
	private static final String DEFUALT_DIRECTORY =
																								"ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/";
	/**
	 * Continental MAF query
	 */
	private static final String[] CONTINENTAL_QUERY = new String[] {"AF", "EAS_AF", "EUR_AF",
																																	"AFR_AF", "AMR_AF", "SAS_AF"};
	/**
	 * Functional query
	 */
	// private static final String[] FUNTCTIONAL_QUERY = new String[] { "GERP", "SVTYPE", "CSQ" };

	/**
	 * Info that will be reported for every variant
	 */
	private static final String[] DEFAULT_HEADER = new String[] {	"ID", "CHR", "POS", "REF", "ALT",
																																"POP_CALLRATE"};

	public enum Location {
												/**
												 * files are on a local drive
												 */
												LOCAL,
												/**
												 * FTP hosted files
												 */
												REMOTE
	}

	private enum VCF_ORGANIZATION {
																	/**
																	 * Each vcf contains a single chromosome
																	 */
																	ONE_PER_CHROMOSOME,
																	/**
																	 * multiple vcfs per chromosome
																	 */
																	MULTI_CHROMOSOME
	}

	// private static volatile int TRACK_VOLATILE = 0;
	private final VCFFileReader vcfSourceReader;
	private final String vcfFile;
	private final Logger log;

	public VcfQuery(String vcfFile, Logger log) {
		super();
		vcfSourceReader = new VCFSourceReader(vcfFile, true);
		this.vcfFile = vcfFile;
		this.log = log;
	}

	private byte getFirstChr() {
		byte firstChr = 0;
		for (VariantContext variantContext : vcfSourceReader) {
			try {
				firstChr = Byte.parseByte(variantContext.getContig());
			} catch (NumberFormatException nfe) {
				firstChr = 0;
			}
			break;
		}
		return firstChr;
	}

	public String getVcfFile() {
		return vcfFile;
	}

	public VCFFileReader getVcfFileReader() {
		return vcfSourceReader;
	}

	public Logger getLog() {
		return log;
	}

	public void close() {
		vcfSourceReader.close();
	}

	/**
	 * Returns all variants contained within the segment
	 *
	 * @param seg
	 * @return
	 */
	private VariantContext[] queryASegment(Segment seg) {
		CloseableIterator<VariantContext> cdl = vcfSourceReader.query(seg.getChr()	+ "", seg.getStart(),
																																	seg.getStop());
		ArrayList<VariantContext> vcs = new ArrayList<VariantContext>();
		while (cdl.hasNext()) {
			vcs.add(cdl.next());
		}
		return vcs.toArray(new VariantContext[vcs.size()]);
	}

	private static class QueryManager implements Callable<QueryResults> {
		private final String vcfFile;
		private final Segment[] seqsToQuery;
		private final String[] infoToExtract;
		private final Logger log;
		private final String outputDir;
		private final VCF_ORGANIZATION org;
		private final VcfPopulation vpop;
		// private VariantContextWriter writer;

		private QueryManager(	String vcfFile, Segment[] seqsToQuery, QueryParams params,
													String outputDir, VariantContextWriter writer, Logger log) {
			super();
			this.vcfFile = vcfFile;
			this.seqsToQuery = seqsToQuery;
			infoToExtract = params.getInfoToExtract();
			this.outputDir = outputDir;
			org = params.getOrg();
			// this.writer = writer;
			vpop = params.getPopulationFile() == null	? null
																								: VcfPopulation.load(	params.getPopulationFile(),
																																			POPULATION_TYPE.ANY, log);
			this.log = log;
		}

		@Override
		public QueryResults call() throws Exception {
			QueryResults qResults = new QueryResults();
			ArrayList<QueryResult> results = new ArrayList<VcfQuery.QueryResult>();
			String tmpFile = outputDir + ext.rootOf(vcfFile) + ".query";
			VariantContextWriter writer = null;

			// if (!Files.exists(tmpFile) || !Files.headerOfFileContainsAll(tmpFile, infoToExtract, log))
			// {
			try {
				Thread.sleep(10000);
			} catch (InterruptedException ie) {
			}

			log.reportTimeInfo("Beginning query with "	+ vcfFile + "using: "
													+ Thread.currentThread().getName());
			VcfQuery vcfQuery = new VcfQuery(vcfFile, log);
			// Set<String> EUR = vpop.getSuperPop().get("EUR");
			if (vpop != null) {
			}
			byte firstChr = vcfQuery.getFirstChr();
			if (org == VCF_ORGANIZATION.ONE_PER_CHROMOSOME) {
				log.reportTimeWarning("Assuming all variants from "	+ ext.rootOf(vcfFile)
															+ " are on chromosome " + firstChr);
			}

			int thisRound = 0;
			long time = System.currentTimeMillis();

			for (Segment element : seqsToQuery) {
				if (org == VCF_ORGANIZATION.MULTI_CHROMOSOME || element.getChr() == firstChr) {
					VariantContext[] vcs = vcfQuery.queryASegment(element);
					System.out.println(vcs.length);

					for (VariantContext vc : vcs) {
						thisRound++;
						if (writer != null) {
							writer.add(vc);
						}

						QueryResult qResult = QueryResult.getFromVariantContext(vc, infoToExtract, log);

						if (qResult != null) {
							qResult.setVariantContext(vc);
							results.add(qResult);
						}
					}
					if (ext.getTimeSince(time, 's') > 10) {
						time = System.currentTimeMillis();
						log.report("...", false, true);
					}
					if (thisRound > 0 && thisRound % 100 == 0) {
						log.reportTimeInfo("Found "	+ thisRound + " variants in " + vcfFile + " \n using: "
																+ Thread.currentThread().getName());
					}
				}
			}

			vcfQuery.close();
			log.reportTimeInfo("Finished query with "	+ vcfFile + "using: "
													+ Thread.currentThread().getName());
			if (results.size() > 0) {
				log.reportTimeInfo("Reporting temporary results to " + tmpFile);
				qResults.setHasTmpFile(true);
				qResults.setTmpFile(tmpFile);
				qResults.setQueryResults(results.toArray(new QueryResult[results.size()]));
				dumpToTmpFile(qResults, infoToExtract, tmpFile, false, log);
			} else {
				qResults.setHasTmpFile(false);
			}
			// } else {
			// log.reportTimeWarning("File " + tmpFile + " exists so this query will be skipped");
			// qResults.setHasTmpFile(true);
			// qResults.setTmpFile(tmpFile);
			// }

			return qResults;
		}
	}

	private static class QueryResults {
		private QueryResult[] queryResults;
		private boolean hasTmpFile;
		private String tmpFile;

		public QueryResults() {
			super();
		}

		public QueryResult[] getQueryResults() {
			return queryResults;
		}

		public void setQueryResults(QueryResult[] queryResults) {
			this.queryResults = queryResults;
		}

		public boolean isHasTmpFile() {
			return hasTmpFile;
		}

		public void setHasTmpFile(boolean hasTmpFile) {
			this.hasTmpFile = hasTmpFile;
		}

		public String getTmpFile() {
			return tmpFile;
		}

		public void setTmpFile(String tmpFile) {
			this.tmpFile = tmpFile;
		}

	}

	/**
	 * Tries to save any progress made
	 */
	private static void dumpToTmpFile(QueryResults queryResults, String[] infoToExtract, String file,
																		boolean append, Logger log) {
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(file, append));
			if (!append) {
				writer.print(Array.toStr(DEFAULT_HEADER) + "\t" + Array.toStr(infoToExtract));
				writer.println();
			}
			if (queryResults.getQueryResults() != null) {
				for (int i = 0; i < queryResults.getQueryResults().length; i++) {
					writer.println(queryResults.getQueryResults()[i].getDisplayString());
				}
			} else if (queryResults.isHasTmpFile()) {
				Vector<String> tmp = HashVec.loadFileToVec(queryResults.getTmpFile(), true, false, false);
				for (int i = 0; i < tmp.size(); i++) {
					writer.println(tmp.get(i));
				}

			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + file);
			log.reportException(e);
		}
	}

	/**
	 * Stores results for and individual variant
	 *
	 */
	private static class QueryResult {
		private final String id;
		private final String contig;
		private final String ref;
		private final String[] alts;
		private final int start;
		private final String callRate;
		private String[] info;
		private VariantContext variantContext;

		private QueryResult(String id, String contig, String ref, String[] alts, int start,
												String callRate) {
			super();
			this.id = id;
			this.contig = contig;
			this.ref = ref;
			this.alts = alts;
			this.start = start;
			this.callRate = callRate;
			// System.out.println(callRate);
		}

		public VariantContext getVariantContext() {
			return variantContext;
		}

		public void setVariantContext(VariantContext variantContext) {
			this.variantContext = variantContext;
		}

		private String getDisplayString() {
			return id	+ "\t" + contig + "\t" + start + "\t" + ref + "\t" + Array.toStr(alts, ",") + "\t"
							+ callRate + "\t" + Array.toStr(info);
		}

		private void setInfo(String[] info) {
			this.info = info;
		}

		private static QueryResult getFromVariantContext(	VariantContext vc, String[] infoToExtract,
																											Logger log) {
			String name = vc.getID();
			String chr = vc.getContig();
			String ref = vc.getReference().getDisplayString();
			Allele[] alt = vc.getAlternateAlleles().toArray(new Allele[vc.getAlternateAlleles().size()]);
			String callRate = getCallRate(vc) + "";
			String[] alts = new String[alt.length];
			for (int i = 0; i < alts.length; i++) {
				alts[i] = alt[i].getDisplayString();
			}
			int start = vc.getStart();
			if (name.equals(".")) {
				name = vc.getContig() + ":" + vc.getStart();
			}
			QueryResult qResult = new QueryResult(name, chr, ref, alts, start, callRate);
			if (infoToExtract != null) {
				String[] info = new String[infoToExtract.length];
				for (int i = 0; i < info.length; i++) {
					info[i] = vc.getCommonInfo().getAttributeAsString(infoToExtract[i], ".");
				}
				qResult.setInfo(info);
			}
			return qResult;
		}

		private static double getCallRate(VariantContext vc) {
			int noCalls = vc.getNoCallCount();
			int numSamps = vc.getNSamples();
			// System.out.println(numSamps + "\t" + noCalls + "\t" + vc.getHetCount() + "\t" +
			// vc.getHomRefCount() + "\t" + vc.getHomVarCount());
			if (numSamps > 0) {
				return (numSamps - noCalls) / numSamps;
			} else {
				return Double.NaN;
			}
		}
	}

	private static QueryResults[] query(String[] fullPathVCFs, Segment[] segs, QueryParams params,
																			int numThreads, Logger log) {
		QueryManager[] qManagers = new QueryManager[fullPathVCFs.length];
		String tmpDir = ext.parseDirectoryOfFile(params.getOutputFileName());
		new File(tmpDir).mkdirs();
		VariantContextWriter writer =
																VCFOps.initWriter(tmpDir	+ ext.rootOf(fullPathVCFs[0])
																									+ ".query.vcf", null,
																									VCFOps.getSequenceDictionary(new VCFSourceReader(	fullPathVCFs[0],
																																																		true)));

		// VCFOps.copyHeader(new VCFSourceReader(fullPathVCFs[0], true), writer,
		// VcfPopulation.load(params.getPopulationFile(), log).getSuperPop().get("EUR"));

		AsyncVariantContextWriter asWriter = new AsyncVariantContextWriter(writer);
		for (int i = 0; i < fullPathVCFs.length; i++) {
			qManagers[i] = new QueryManager(fullPathVCFs[i], segs, params, tmpDir, writer, log);
		}
		log.reportTimeInfo("Attempting to extract the following:\n"
												+ Array.toStr(params.getInfoToExtract(), "\n"));
		WorkerHive<QueryResults> hive = new WorkerHive<QueryResults>(numThreads, 10, log);
		hive.setReportEvery(1);
		hive.addCallables(qManagers);
		hive.execute(true);
		ArrayList<QueryResults> results = hive.getResults();
		asWriter.close();
		VariantContextWriter writer2 = VCFOps.initWriter(tmpDir	+ "query.vcf", null,
																											VCFOps.getSequenceDictionary(new VCFSourceReader(	fullPathVCFs[0],
																																																				true)));

		// VCFOps.copyHeader(new VCFSourceReader(fullPathVCFs[0], true), writer2,
		// VcfPopulation.load(params.getPopulationFile(), log).getSuperPop().get("EUR"));

		for (int i = 0; i < results.size(); i++) {
			if (results.get(i).getQueryResults() != null) {
				for (int j = 0; j < results.get(i).getQueryResults().length; j++) {
					writer2.add(results.get(i).getQueryResults()[j].getVariantContext());
				}
			}
		}
		writer2.close();
		return results.toArray(new QueryResults[results.size()]);
	}

	/**
	 * Organizes the query
	 *
	 */
	private static class QueryParams {

		private String dir;
		private String[] infoToExtract;
		private Location location;
		private String segFile;
		private String outputFileName;
		private String populationFile;
		private VCF_ORGANIZATION org;

		public QueryParams(	String dir, String[] infoToExtract, Location location,
												VCF_ORGANIZATION org) {
			super();
			this.dir = dir;
			this.infoToExtract = infoToExtract;
			this.location = location;
			this.org = org;
			outputFileName = null;
		}

		public Location getLocation() {
			return location;
		}

		public void setLocation(Location location) {
			this.location = location;
		}

		public String getPopulationFile() {
			return populationFile;
		}

		public void setPopulationFile(String populationFile) {
			this.populationFile = populationFile;
		}

		public VCF_ORGANIZATION getOrg() {
			return org;
		}

		public void setOrg(VCF_ORGANIZATION org) {
			this.org = org;
		}

		public String getDir() {
			return dir;
		}

		public void setDir(String dir) {
			this.dir = dir;
		}

		public String[] getInfoToExtract() {
			return infoToExtract;
		}

		public void setInfoToExtract(String[] infoToExtract) {
			this.infoToExtract = infoToExtract;
		}

		public String getSegFile() {
			return segFile;
		}

		public void setSegFile(String segFile) {
			this.segFile = segFile;
		}

		public String getOutputFileName() {
			return outputFileName;
		}

		public void setOutputFileName(String outputFileName) {
			this.outputFileName = outputFileName;
		}

	}

	public static void queryDir(QueryParams qParams, int numThreads, int chrColumn, int startColumn,
															int stopColumn, Logger log) {
		Segment[] segs = Segment.loadRegions(	qParams.getSegFile(), chrColumn, startColumn, stopColumn,
																					0, true, true, true, 0);
		if (segs == null || segs.length < 1) {
			log.reportError("Did not find any valid segments in file " + qParams.getSegFile());
		} else {
			log.reportTimeInfo("Found " + segs.length + " seqments to query ");
			String[] vcfs = new String[0];
			String[] gzVcfs = new String[0];
			if (qParams.getLocation() == Location.REMOTE) {
				vcfs = Files.parseRemoteFTPFiles(qParams.getDir(), ".vcf.gz", log);
			} else {
				vcfs = Files.listFullPaths(qParams.getDir(), ".vcf", false);
				gzVcfs = Files.listFullPaths(qParams.getDir(), ".vcf.gz", false);
			}
			String[] all = Array.concatAll(vcfs, gzVcfs);
			if (all.length < 1) {
				log.reportError("Did not find any .vcf or .vcf.gz files in directory "
														+ qParams.getDir());
			} else {
				log.reportTimeInfo("Found " + all.length + " files to query from " + qParams.getDir());
				QueryResults[] queryResults = query(all, segs, qParams, numThreads, log);
				for (int i = 0; i < queryResults.length; i++) {
					dumpToTmpFile(queryResults[i], qParams.getInfoToExtract(), qParams.getOutputFileName(),
												i > 0, log);
				}
			}
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		QueryParams params = new QueryParams(	DEFUALT_DIRECTORY, CONTINENTAL_QUERY, Location.REMOTE,
																					VCF_ORGANIZATION.ONE_PER_CHROMOSOME);
		int numThreads = 2;
		int chrColumn = 0;
		int startColumn = 1;
		int stopColumn = 1;

		Logger log;

		String usage = "\n" + "seq.analysis.VcfQuery requires 0-1 arguments\n";
		usage +=
					"   (1) full path to a file of segments (no header,tab-delimited, chr start stop) (i.e. segFile= (no default))\n"
							+ "";
		usage += "	 OPTIONAL:\n";
		usage +=
					"   (2) full path to an output filename (i.e. out= (no default, based off root of segment file))\n"
							+ "";
		usage += "   (3) number of threads for the query (i.e. "	+ PSF.Ext.NUM_THREADS_COMMAND
							+ numThreads + "(default)\n" + "";
		usage += "   (4) full path to a directory containing vcf files to query (i.e. dir= "
							+ DEFUALT_DIRECTORY + " (default)\n" + "";
		usage += "   (5) a comma-delimited list of common info to extract (i.e. info= "
							+ Array.toStr(CONTINENTAL_QUERY, ",") + " (default)\n" + "";
		usage +=
					"   (6) the directory of vcfs to search is a local directory (i.e. -local (not the default)\n"
							+ "";
		usage +=
					"   (7) the vcf files are not one file per chromosome (i.e. -multi (not the default, the first variant's chromosome is used as a filter)\n"
							+ "";
		usage += "   (8) the chromosome column in the segment file (i.e. chr="	+ chrColumn
							+ "(default)\n" + "";
		usage += "   (9) the start column in the segment file (i.e. start="	+ startColumn
							+ "(default)\n" + "";
		usage += "   (10) the stop column in the segment file (i.e. stopColumn="	+ stopColumn
							+ "(default)\n" + "";
		usage += "   (11) full path to a population file (i.e. vpop=null (no default)\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("segFile=")) {
				params.setSegFile(ext.parseStringArg(arg, ""));
				numArgs--;
			} else if (arg.startsWith("dir=")) {
				params.setDir(ext.parseStringArg(arg, ""));
				numArgs--;
			} else if (arg.startsWith("out=")) {
				params.setOutputFileName(ext.parseStringArg(arg, ""));
				numArgs--;
			} else if (arg.startsWith("info=")) {
				params.setInfoToExtract(ext.parseStringArg(arg, "").split(","));
				numArgs--;
			} else if (arg.startsWith("vpop=")) {
				params.setPopulationFile(ext.parseStringArg(arg, ""));
				numArgs--;
			} else if (arg.startsWith("-local")) {
				params.setLocation(Location.LOCAL);
				numArgs--;
			} else if (arg.startsWith("-multi")) {
				params.setOrg(VCF_ORGANIZATION.MULTI_CHROMOSOME);
				numArgs--;
			} else if (arg.startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
				numThreads = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith("chr=")) {
				chrColumn = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith("start=")) {
				startColumn = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith("stop=")) {
				stopColumn = ext.parseIntArg(arg);
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
			if (params.getSegFile() == null) {
				log = new Logger();
				log.reportError("Must provide a valid segment file");
			} else {
				if (params.getOutputFileName() == null) {
					params.setOutputFileName(ext.addToRoot(params.getSegFile(), ".query"));
				}
				log = new Logger(params.getSegFile() + ".log");
				queryDir(params, numThreads, chrColumn, startColumn, stopColumn, log);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}

// bash script for gmaf
// chr=0 start=1 stop=1 dir=/home/pankrat2/public/bin/ref/
//
// echo "usage:"
// echo "-c # sets the chromosome column of the segment file (default $chr)"
// echo "-s # sets the start column of the segment file (default $start)"
// echo "-p # sets the start column of the segment file (default $stop)"
// echo "-d aDir sets the directory of vcf files to search (default $dir)"
// echo "All other arguments are treated as a file of segments"
//
//
// while getopts c:s:p:d: opt; do
// case $opt in
// c)
// chr=$OPTARG
// ;;
// s)
// start=$OPTARG
// ;;
// p)
// stop=$OPTARG
// ;;
// d)
// dir=$OPTARG
// ;;
//
// esac
// done
//
// shift $((OPTIND - 1))
//
//
//
//
// echo "Current Params"
// echo "chromosome column set to $chr"
// echo "start column set to $start"
// echo "stop column set to $stop"
// echo "directory to search set to $dir"
//
//
//
// for file in "$@"
// do
// jcp seq.analysis.VcfQuery -local -multi dir="$dir" segFile="$file" chr="$chr" start="$start"
// stop="$stop"
// done
