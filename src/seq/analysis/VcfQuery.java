package seq.analysis;

import filesys.Segment;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Vector;
import java.util.concurrent.Callable;

import common.Array;
import common.Files;
import common.HashVec;
import common.Logger;
import common.PSF;
import common.WorkerHive;
import common.ext;

/**
 * Class to perform common info queries on local and remote vcf files<br>
 * For spot-checking http://browser.1000genomes.org/Homo_sapiens/Info/Index?db=core;r=1:109818030-109819030;v=rs646776;vdb=variation;vf=510141 http://browser.1000genomes.org/Homo_sapiens/Variation/Population?r=1:109818030-109819030;v=rs646776;vdb=variation;vf=510141
 */
public class VcfQuery {
	/**
	 * This ftp site currently hosts indexed (.tbi) vcf (.vcf.gz) files that can searched
	 */
	private static final String DEFUALT_DIRECTORY = "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/supporting/vcf_with_sample_level_annotation/";
	/**
	 * Continental MAF query
	 */
	private static final String[] CONTINENTAL_QUERY = new String[] { "AF", "EAS_AF", "EUR_AF", "AFR_AF", "AMR_AF", "SAS_AF" };

	/**
	 * Info that will be reported for every variant
	 */
	private static final String[] DEFAULT_HEADER = new String[] { "ID", "CHR", "POS", "REF", "ALT" };

	private enum Location {
		/**
		 * files are on a local drive
		 */
		LOCAL, /**
		 * FTP hosted files
		 */
		REMOTE
	}

	private enum VCF_ORGANIZATION {
		/**
		 * Each vcf contains a single chromosome
		 */
		ONE_PER_CHROMOSOME, /**
		 * multiple vcfs per chromosome
		 */
		MULTI_CHROMOSOME
	}

	// private static volatile int TRACK_VOLATILE = 0;
	private VCFFileReader vcfFileReader;
	private String vcfFile;
	private Logger log;

	public VcfQuery(String vcfFile, Logger log) {
		super();
		this.vcfFileReader = new VCFFileReader(vcfFile, true);
		this.vcfFile = vcfFile;
		this.log = log;
	}

	private byte getFirstChr() {
		byte firstChr = 0;
		for (VariantContext variantContext : vcfFileReader) {
			try {
				firstChr = Byte.parseByte(variantContext.getChr());
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

	public Logger getLog() {
		return log;
	}

	public void close() {
		vcfFileReader.close();
	}

	/**
	 * Returns all variants contained within the segment
	 * 
	 * @param seg
	 * @return
	 */
	private VariantContext[] queryASegment(Segment seg) {
		CloseableIterator<VariantContext> cdl = vcfFileReader.query(seg.getChr() + "", seg.getStart(), seg.getStop());
		ArrayList<VariantContext> vcs = new ArrayList<VariantContext>();
		while (cdl.hasNext()) {
			vcs.add(cdl.next());
		}
		return vcs.toArray(new VariantContext[vcs.size()]);
	}

	private static class QueryManager implements Callable<QueryResults> {
		private String vcfFile;
		private Segment[] seqsToQuery;
		private String[] infoToExtract;
		private Logger log;
		private String outputDir;
		private VCF_ORGANIZATION org;

		private QueryManager(String vcfFile, Segment[] seqsToQuery, String[] infoToExtract, String outputDir, VCF_ORGANIZATION org, Logger log) {
			super();
			this.vcfFile = vcfFile;
			this.seqsToQuery = seqsToQuery;
			this.infoToExtract = infoToExtract;
			this.outputDir = outputDir;
			this.org = org;
			this.log = log;
		}

		@Override
		public QueryResults call() throws Exception {
			QueryResults qResults = new QueryResults();
			ArrayList<QueryResult> results = new ArrayList<VcfQuery.QueryResult>();
			String tmpFile = outputDir + ext.rootOf(vcfFile) + ".query";
			if (!Files.exists(tmpFile) || !Files.headerOfFileContainsAll(tmpFile, infoToExtract, log)) {
				try {
					Thread.sleep(10000);
				} catch (InterruptedException ie) {
				}

				log.reportTimeInfo("Beginning query with " + vcfFile + "using: " + Thread.currentThread().getName());
				VcfQuery vcfQuery = new VcfQuery(vcfFile, log);
				byte firstChr = vcfQuery.getFirstChr();
				if (org == VCF_ORGANIZATION.ONE_PER_CHROMOSOME) {
					log.reportTimeWarning("Assuming all variants from " + ext.rootOf(vcfFile) + " are on chromosome " + firstChr);
				}

				int thisRound = 0;
				long time = System.currentTimeMillis();

				for (int i = 0; i < seqsToQuery.length; i++) {
					if (org == VCF_ORGANIZATION.MULTI_CHROMOSOME || seqsToQuery[i].getChr() == firstChr) {
						VariantContext[] vcs = vcfQuery.queryASegment(seqsToQuery[i]);
						for (int j = 0; j < vcs.length; j++) {
							thisRound++;
							QueryResult qResult = QueryResult.getFromVariantContext(vcs[j], infoToExtract, log);
							if (qResult != null) {
								results.add(qResult);
							}
						}
						if (ext.getTimeSince(time, 's') > 10) {
							time = System.currentTimeMillis();
							log.report("...", false, true);
						}
						if (thisRound > 0 && thisRound % 100 == 0) {
							log.reportTimeInfo("Found " + thisRound + " variants in " + vcfFile + " \n using: " + Thread.currentThread().getName());
						}
					}
				}
				vcfQuery.close();
				log.reportTimeInfo("Finished query with " + vcfFile + "using: " + Thread.currentThread().getName());
				if (results.size() > 0) {
					log.reportTimeInfo("Reporting temporary results to " + tmpFile);
					qResults.setHasTmpFile(true);
					qResults.setTmpFile(tmpFile);
					qResults.setQueryResults(results.toArray(new QueryResult[results.size()]));
					dumpToTmpFile(qResults, infoToExtract, tmpFile, false, log);
				} else {
					qResults.setHasTmpFile(false);
				}
			} else {
				log.reportTimeWarning("File " + tmpFile + " exists so this query will be skipped");
				qResults.setHasTmpFile(true);
				qResults.setTmpFile(tmpFile);
			}

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
	private static void dumpToTmpFile(QueryResults queryResults, String[] infoToExtract, String file, boolean append, Logger log) {
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
		private String id;
		private String contig;
		private String ref;
		private String[] alts;
		private int start;
		private String[] info;

		private QueryResult(String id, String contig, String ref, String[] alts, int start) {
			super();
			this.id = id;
			this.contig = contig;
			this.ref = ref;
			this.alts = alts;
			this.start = start;
		}

		private String getDisplayString() {
			return id + "\t" + contig + "\t" + start + "\t" + ref + "\t" + Array.toStr(alts, ",") + "\t" + Array.toStr(info);
		}

		private void setInfo(String[] info) {
			this.info = info;
		}

		private static QueryResult getFromVariantContext(VariantContext vc, String[] infoToExtract, Logger log) {
			String name = vc.getID();
			String chr = vc.getChr();
			String ref = vc.getReference().getDisplayString();
			Allele[] alt = vc.getAlternateAlleles().toArray(new Allele[vc.getAlternateAlleles().size()]);
			String[] alts = new String[alt.length];
			for (int i = 0; i < alts.length; i++) {
				alts[i] = alt[i].getDisplayString();
			}
			int start = vc.getStart();
			if (name.equals(".")) {
				name = vc.getChr() + ":" + vc.getStart();
			}
			QueryResult qResult = new QueryResult(name, chr, ref, alts, start);
			if (infoToExtract != null) {
				String[] info = new String[infoToExtract.length];
				for (int i = 0; i < info.length; i++) {
					info[i] = vc.getCommonInfo().getAttributeAsString(infoToExtract[i], ".");
				}
				qResult.setInfo(info);
			}
			return qResult;
		}
	}

	/**
	 * Attempts to grab all file paths of a type from an ftp site
	 * 
	 * @param ftpdirAddress
	 *            a remote directory of an ftp site
	 * @param type
	 *            the type of file to collect
	 * @param log
	 * @return
	 */
	private static String[] parseRemoteFTPFiles(String ftpdirAddress, String type, Logger log) {
		ArrayList<String> remoteVcfs = new ArrayList<String>();
		URL url;
		if (!ftpdirAddress.startsWith("ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/")) {
			log.reportTimeWarning("Did not detect that " + ftpdirAddress + " was an ftp address starting with ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/");
			log.reportTimeWarning("\t this parsing method is therefore un-tested");
		}
		try {
			url = new URL(ftpdirAddress);
			BufferedReader in = new BufferedReader(new InputStreamReader(url.openStream()));
			String inputLine;

			while ((inputLine = in.readLine()) != null) {
				String[] line = inputLine.trim().split("[\\s]+");
				int possibleVcfIndex = ext.indexOfEndsWith(type, line, false);
				if (possibleVcfIndex > 0) {
					remoteVcfs.add(ftpdirAddress + line[possibleVcfIndex]);
				}
			}
			in.close();

		} catch (MalformedURLException e) {
			log.reportTimeError("malformed URL " + ftpdirAddress);
			e.printStackTrace();
		} catch (IOException e) {
			log.reportIOException(ftpdirAddress);
			e.printStackTrace();
		}
		return remoteVcfs.toArray(new String[remoteVcfs.size()]);
	}

	private static QueryResults[] query(String[] fullPathVCFs, Segment[] segs, QueryParams params, int numThreads, Logger log) {
		QueryManager[] qManagers = new QueryManager[fullPathVCFs.length];
		String tmpDir = ext.parseDirectoryOfFile(params.getOutputFileName());
		new File(tmpDir).mkdirs();
		for (int i = 0; i < fullPathVCFs.length; i++) {
			qManagers[i] = new QueryManager(fullPathVCFs[i], segs, params.getInfoToExtract(), tmpDir, params.getOrg(), log);
		}
		log.reportTimeInfo("Attempting to extract the following:\n" + Array.toStr(params.getInfoToExtract(), "\n"));
		WorkerHive<QueryResults> hive = new WorkerHive<QueryResults>(numThreads, 10, log);
		hive.setReportEvery(1);
		hive.addCallables(qManagers);
		hive.execute(true);
		ArrayList<QueryResults> results = hive.getResults();
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
		private VCF_ORGANIZATION org;

		public QueryParams(String dir, String[] infoToExtract, Location location, VCF_ORGANIZATION org) {
			super();
			this.dir = dir;
			this.infoToExtract = infoToExtract;
			this.location = location;
			this.org = org;
			this.outputFileName = null;
		}

		public Location getLocation() {
			return location;
		}

		public void setLocation(Location location) {
			this.location = location;
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

	public static void queryDir(QueryParams qParams, int numThreads, Logger log) {
		Segment[] segs = Segment.loadRegions(qParams.getSegFile(), 0, 1, 2, 0, true, true, true);
		if (segs == null || segs.length < 1) {
			log.reportTimeError("Did not find any valid segments in file " + qParams.getSegFile());
		} else {
			log.reportTimeInfo("Found " + segs.length + " seqments to query ");
			String[] vcfs = new String[0];
			String[] gzVcfs = new String[0];
			if (qParams.getLocation() == Location.REMOTE) {
				vcfs = parseRemoteFTPFiles(qParams.getDir(), ".vcf.gz", log);
			} else {
				vcfs = Files.listFullPaths(qParams.getDir(), ".vcf", false);
				gzVcfs = Files.listFullPaths(qParams.getDir(), ".vcf.gz", false);
			}
			String[] all = Array.concatAll(vcfs, gzVcfs);
			if (all.length < 1) {
				log.reportTimeError("Did not find any .vcf or .vcf.gz files in directory " + qParams.getDir());
			} else {
				log.reportTimeInfo("Found " + all.length + " files to query from " + qParams.getDir());
				QueryResults[] queryResults = query(all, segs, qParams, numThreads, log);
				for (int i = 0; i < queryResults.length; i++) {
					dumpToTmpFile(queryResults[i], qParams.getInfoToExtract(), qParams.getOutputFileName(), i > 0, log);
				}
			}
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		QueryParams params = new QueryParams(DEFUALT_DIRECTORY, CONTINENTAL_QUERY, Location.REMOTE, VCF_ORGANIZATION.ONE_PER_CHROMOSOME);
		int numThreads = 2;
		Logger log;

		String usage = "\n" + "seq.analysis.VcfQuery requires 0-1 arguments\n";
		usage += "   (1) full path to a file of segments (no header,tab-delimited, chr start stop) (i.e. segFile= (no default))\n" + "";
		usage += "	 OPTIONAL:";
		usage += "   (2) full path to an output filename (i.e. out= (no default, based off root of segment file))\n" + "";
		usage += "   (3) number of threads for the query (i.e. " + PSF.Ext.NUM_THREADS_COMMAND + numThreads + "(default)\n" + "";
		usage += "   (4) full path to a directory containing vcf files to query (i.e. dir= " + DEFUALT_DIRECTORY + " (default)\n" + "";
		usage += "   (5) a comma-delimited list of common info to extract (i.e. info= " + Array.toStr(CONTINENTAL_QUERY, ",") + " (default)\n" + "";
		usage += "   (6) the directory of vcfs to search is a local directory (i.e. -local (not the default)\n" + "";
		usage += "   (7) the vcf files are not one file per chromosome (i.e. -multi (not the default, the first variant's chromosome is used as a filter)\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("segFile=")) {
				params.setSegFile(ext.parseStringArg(args[i], ""));
				numArgs--;
			} else if (args[i].startsWith("dir=")) {
				params.setDir(ext.parseStringArg(args[i], ""));
				numArgs--;
			} else if (args[i].startsWith("out=")) {
				params.setOutputFileName(ext.parseStringArg(args[i], ""));
				numArgs--;
			} else if (args[i].startsWith("info=")) {
				params.setInfoToExtract(ext.parseStringArg(args[i], "").split(","));
				numArgs--;
			} else if (args[i].startsWith("-local")) {
				params.setLocation(Location.LOCAL);
				numArgs--;
			} else if (args[i].startsWith("-mulit")) {
				params.setOrg(VCF_ORGANIZATION.MULTI_CHROMOSOME);
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
			if (params.getSegFile() == null) {
				log = new Logger();
				log.reportTimeError("Must provide a valid segment file");
			} else {
				if (params.getOutputFileName() == null) {
					params.setOutputFileName(ext.addToRoot(params.getSegFile(), ".query"));
				}
				log = new Logger(params.getSegFile() + ".log");
				queryDir(params, numThreads, log);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
