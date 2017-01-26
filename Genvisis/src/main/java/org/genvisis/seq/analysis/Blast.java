package org.genvisis.seq.analysis;

import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.concurrent.Callable;

import org.genvisis.common.ArrayUtils;
import org.genvisis.common.CmdLine;
import org.genvisis.common.CmdLineProcess;
import org.genvisis.common.CmdLineProcess.ERR_Mode;
import org.genvisis.common.CmdLineProcess.INPUT_Mode;
import org.genvisis.common.CmdLineProcess.OUTPUT_Mode;
import org.genvisis.common.CmdLineProcess.StandardInputProvider;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.Positions;
import org.genvisis.common.ext;
import org.genvisis.filesys.Segment;
import org.genvisis.stats.Histogram.DynamicHistogram;

import htsjdk.tribble.annotation.Strand;

public class Blast {
	public static final String[] DB_EXTs = new String[] {".nsq", ".nin", ".nhr"};
	public static final String[] BLAST_HEADER = new String[] {"query id", "subject id", "% identity",
																														"alignment length", "mismatches",
																														"gap opens", "q. start", "q. end",
																														"s. start", "s. end", "evalue",
																														"bit score", "BTOP"};
	private static final String DB = "-db";
	private static final String IN = "-in";
	private static final String DB_TYPE = "-dbtype";

	private static final String OUT_FMT = "-outfmt";
	private static final String E = "-evalue";

	private static final int DEFAULT_OUT_FMT = 7;
	private static final String WORD_SIZE = "-word_size";
	public static final int DEFAULT_WORD_SIZE = 11;
	public static final double DEFAULT_EVALUE = 10;

	private enum BLAST_COMMANDS {
																MAKE_DB("makeblastdb"), BLASTN("blastn");
		private String command;

		private BLAST_COMMANDS(String command) {
			this.command = command;
		}

		public String getCommand() {
			return command;
		}
	}

	private enum BLAST_DB_TYPE {
															NUCL("nucl");
		private String type;

		private BLAST_DB_TYPE(String type) {
			this.type = type;
		}

		public String getTYPE() {
			return type;
		}
	}

	private final String fastaDb;
	private final Logger log;
	private boolean fail, taxonMode;
	private final boolean overwriteExisting;
	private final boolean verbose;
	private final int blastWordSize, reportWordSize;
	private double evalue;

	public Blast(	String fastaDb, int blastWordSize, int reportWordSize, Logger log,
								boolean overwriteExisting, boolean verbose) {
		super();

		this.fastaDb = fastaDb;
		this.blastWordSize = blastWordSize;
		this.reportWordSize = reportWordSize;
		this.log = log;
		fail = !verify(fastaDb, log);
		this.overwriteExisting = overwriteExisting;
		this.verbose = verbose;
		taxonMode = false;
		evalue = DEFAULT_EVALUE;
		if (!fail) {
			fail = !initDb(BLAST_DB_TYPE.NUCL, fastaDb, log);
		}
	}


	public boolean isOverwriteExisting() {
		return overwriteExisting;
	}


	public void setEvalue(double evalue) {
		this.evalue = evalue;
	}

	public Logger getLog() {
		return log;
	}

	public void setTaxonMode(boolean taxonMode) {
		this.taxonMode = taxonMode;
	}

	public BlastResultsSummary[] blastSequence(FastaEntry[] fastaEntries, PrintWriter tmpFile) {
		BlastResultsSummary[] bSummaries = new BlastResultsSummary[fastaEntries.length];

		if (!Files.isWindows()) {
			if (!fail) {
				String[] command = new String[] {	BLAST_COMMANDS.BLASTN.getCommand(), DB, fastaDb,
																					OUT_FMT, DEFAULT_OUT_FMT	+ " std"
																										+ (taxonMode ? " staxids" : " btop"),
																					WORD_SIZE, blastWordSize + "",
																					evalue != DEFAULT_EVALUE ? E : "",
																					evalue != DEFAULT_EVALUE ? evalue + "" : ""};
				FastaEntryInputStream fStream = new FastaEntryInputStream(fastaEntries, log);
				CmdLineProcess.Builder builder = new CmdLineProcess.Builder();
				builder.STIN(fStream);
				builder.inputMode(INPUT_Mode.STIN);
				builder.outputMode(OUTPUT_Mode.STOUT_CAPTURE_ITERATOR);
				builder.errorMode(ERR_Mode.STERR_CAPTURE_BY_LOG);
				builder.log(log);
				builder.verbose(verbose);
				CmdLineProcess cmdLineProcess = builder.build(command);
				for (int i = 0; i < bSummaries.length; i++) {
					bSummaries[i] = new BlastResultsSummary(fastaEntries[i].getName(), taxonMode,
																									reportWordSize);
				}
				while (cmdLineProcess.hasNext()) {
					String line = cmdLineProcess.next();

					String[] result = line.trim().split("[\\s]+");
					for (int i = 0; i < fastaEntries.length; i++) {
						if (result[0].equals(fastaEntries[i].getName())) {
							BlastResults bResults = new BlastResults(result, log);
							if (bResults.getAlignmentLength() >= reportWordSize) {
								bSummaries[i].addBlastResult(bResults, log);
								if (tmpFile != null) {
									tmpFile.println(ArrayUtils.toStr(bResults.getResults()));
								}
							}
						}
					}
				}
				cmdLineProcess.waitFor();
			}
		} else {
			log.reportError("This command can only be used on *.nix systems, apologies");
		}
		return bSummaries;
	}

	public static boolean initDb(BLAST_DB_TYPE type, String fastaDb, Logger log) {
		boolean dbCreated = true;
		if (!Files.exists("", getDBFiles(fastaDb))) {
			String[] dbCommand = new String[] {	BLAST_COMMANDS.MAKE_DB.getCommand(), IN, fastaDb, DB_TYPE,
																					type.getTYPE()};
			dbCreated = CmdLine.runCommandWithFileChecks(	dbCommand, "", new String[] {fastaDb},
																										getDBFiles(fastaDb), true, false, false, log);
		} else {
			log.reportTimeInfo("Using existing data base files : "	+ ext.rootOf(fastaDb) + " ("
													+ ArrayUtils.toStr(DB_EXTs, ",") + ")");

		}
		return dbCreated;
	}

	private static boolean verify(String fastaDb, Logger log) {
		boolean verified = true;
		if (!Files.exists(fastaDb) && !Files.exists("", getDBFiles(fastaDb))) {

			log.reportError("Could not find fasta file for database " + fastaDb);
			log.reportError(" " + fastaDb);
			verified = false;
		}
		if (!CmdLine.run(BLAST_COMMANDS.BLASTN.getCommand(), "")) {
			log.reportError("It is assumed that the program "	+ BLAST_COMMANDS.BLASTN.getCommand()
													+ " can be found on the system's path, please install before continuing");
			verified = false;
		}
		if (!Files.exists("", getDBFiles(fastaDb))
				&& !CmdLine.run(BLAST_COMMANDS.MAKE_DB.getCommand(), "")) {
			log.reportError("It is assumed that the program "	+ BLAST_COMMANDS.BLASTN.getCommand()
													+ " can be found on the system's path, or the following files are present...");
			log.reportError(ArrayUtils.toStr(getDBFiles(fastaDb), "\n"));
			verified = false;
		}
		return verified;
	}

	private static String[] getDBFiles(String fastaDb) {
		String[] dbFiles = new String[DB_EXTs.length];
		for (int i = 0; i < DB_EXTs.length; i++) {
			dbFiles[i] = fastaDb + DB_EXTs[i];
		}
		return dbFiles;
	}

	public static class FastaEntryInputStream implements StandardInputProvider {
		private final FastaEntry[] fastaEntries;
		private int count;
		private final Logger log;

		public FastaEntryInputStream(FastaEntry[] fastaEntries, Logger log) {
			super();
			this.fastaEntries = fastaEntries;
			count = 0;
			this.log = log;
		}

		public Logger getLog() {
			return log;
		}

		@Override
		public boolean hasNext() {
			return count < fastaEntries.length;
			// TODO Auto-generated method stub
		}

		@Override
		public String next() {
			String entry = fastaEntries[count].getFastaFormat();
			count++;
			return entry;
		}

		@Override
		public void remove() {
			// TODO Auto-generated method stub

		}

	}

	public static class FastaEntry {
		protected String name;
		protected String sequence;

		public FastaEntry(String name, String sequence) {
			super();
			this.name = name;
			this.sequence = sequence;
		}

		public String getFastaFormat() {
			return ">" + name + "\n" + sequence;
		}

		public String getSequence() {
			return sequence;
		}

		public String getSubProcessForLinuxShell() {
			String sub = PSF.Ext.REV_CARROT;
			sub += "(echo -e '" + getFastaFormat() + "')";
			sub = "\"'echo -e \"" + getFastaFormat() + "\"'\"";
			sub = getFastaFormat();

			System.out.println(sub);
			return sub;
		}

		public String getName() {
			return name;
		}

	}

	public static class BlastResults {
		private final String queryID;
		private final String subjectID;
		private final double percentIdentity;
		private final int alignmentLength;
		private final int mismatches;
		private final int gapOpens;
		private final int qstart;
		private final int qstop;
		private final int sstart;
		private final int sstop;
		private final double evalue;
		private final double bitScore;
		private String taxID;
		private String btop;
		private final Logger log;

		public BlastResults(String[] blastLine, Logger log) {
			this(blastLine, false, log);
		}

		public boolean isStrandFlipped() {
			return getSstart() > getSstop();
		}

		public BlastResults(String[] blastLine, boolean taxonMode, Logger log) {
			queryID = blastLine[0];
			subjectID = blastLine[1];
			percentIdentity = tryDouble(blastLine[2], log);
			alignmentLength = tryInt(blastLine[3], log);
			mismatches = tryInt(blastLine[4], log);
			gapOpens = tryInt(blastLine[5], log);
			qstart = tryInt(blastLine[6], log);
			qstop = tryInt(blastLine[7], log);
			sstart = tryInt(blastLine[8], log);
			sstop = tryInt(blastLine[9], log);
			evalue = tryDouble(blastLine[10], log);
			bitScore = tryDouble(blastLine[11], log);
			if (taxonMode && blastLine.length > 12) {
				taxID = blastLine[12];
			} else if (blastLine.length > 12) {
				btop = blastLine[12];
			} else {
				taxID = "NA";
				btop = "NA";
			}
			this.log = log;
		}

		public Segment getSegment() {
			return new Segment(	Positions.chromosomeNumber(subjectID, false, new Logger()),
													Math.min(sstart, sstop), Math.max(sstart, sstop));

		}

		public Strand determineStrand() {
			if (sstart > sstop) {
				return Strand.NEGATIVE;
			} else if (sstop >= sstart) {
				return Strand.POSITIVE;
			} else {
				log.reportTimeWarning("Could not determine strand for " + getResults());
				return Strand.NONE;
			}
		}

		public int getMismatches() {
			return mismatches;
		}

		public int getGapOpens() {
			return gapOpens;
		}

		public double getEvalue() {
			return evalue;
		}

		public String getTaxID() {
			return taxID;
		}

		public int getSstart() {
			return sstart;
		}

		public int getSstop() {
			return sstop;
		}

		public int getQstart() {
			return qstart;
		}

		public String getBtop() {
			return btop;
		}

		public int getQstop() {
			return qstop;
		}

		public String getQueryID() {
			return queryID;
		}

		public String getSubjectID() {
			return subjectID;
		}

		public double getPercentIdentity() {
			return percentIdentity;
		}

		public int getAlignmentLength() {
			return alignmentLength;
		}

		public String[] getResults() {
			return new String[] {	queryID, subjectID, percentIdentity + "", alignmentLength + "",
														mismatches + "", gapOpens + "", qstart + "", qstop + "", sstart + "",
														sstop + "", evalue + "", bitScore + "", btop};
		}

		private double tryDouble(String ad, Logger log) {
			double d = Double.NaN;
			try {
				d = Double.parseDouble(ad);
			} catch (NumberFormatException nfe) {
				log.reportError("Invalid number " + ad + " found ");
			}
			return d;
		}

		private int tryInt(String ai, Logger log) {
			int i = -1;
			try {
				i = Integer.parseInt(ai);
			} catch (NumberFormatException nfe) {
				log.reportError("Invalid number " + ai + " found ");
			}
			return i;
		}

	}

	public static class BlastResultsSummary {
		private final DynamicHistogram percentIdentityHistogram, evalueHistogram;
		private final String name;
		private final Hashtable<String, Integer> hitCounts;
		private Segment perfectMatchSegment;
		private int numPerfectMatches;
		private final boolean taxonMode;

		public BlastResultsSummary(String name, boolean taxonMode, int reportWordSize) {
			super();
			this.name = name;
			this.taxonMode = taxonMode;
			hitCounts = new Hashtable<String, Integer>();
			percentIdentityHistogram = new DynamicHistogram(reportWordSize, 100, 0);
			evalueHistogram = new DynamicHistogram(0, 1, 2);
			numPerfectMatches = 0;
			perfectMatchSegment = null;
		}

		public void addBlastResult(BlastResults blastResults, Logger log) {
			if (!blastResults.getQueryID().equals(name)) {
				log.reportError("Query id and summary name do not match");
			} else {
				percentIdentityHistogram.addDataPointToHistogram(blastResults.getPercentIdentity());
				evalueHistogram.addDataPointToHistogram(blastResults.getEvalue());
				if (!hitCounts.containsKey(blastResults.getSubjectID())) {
					if (taxonMode) {
						hitCounts.put(blastResults.getTaxID(), 0);
					} else {
						hitCounts.put(blastResults.getSubjectID(), 0);
					}
				}

				int cur = hitCounts.get(taxonMode ? blastResults.getTaxID() : blastResults.getSubjectID());
				hitCounts.put(taxonMode ? blastResults.getTaxID() : blastResults.getSubjectID(), cur + 1);

				if (blastResults.getPercentIdentity() == 100) {
					if (numPerfectMatches == 0) {
						if (blastResults.getSubjectID().startsWith("chr")) {
							perfectMatchSegment =
																	new Segment(Positions.chromosomeNumber(blastResults.getSubjectID()),
																							blastResults.getSstart(), blastResults.getSstop());
						}

					} else {
						perfectMatchSegment = null;
					}
					numPerfectMatches++;
				}
			}
		}

		public Hashtable<String, Integer> getHitCounts() {
			return hitCounts;
		}

		public String getName() {
			return name;
		}

		public Segment getPerfectMatchSegment() {
			return perfectMatchSegment;
		}

		public int getNumPerfectMatches() {
			return numPerfectMatches;
		}

		public DynamicHistogram getPercentIdentityHistogram() {
			return percentIdentityHistogram;
		}

		public DynamicHistogram getEvalueHistogram() {
			return evalueHistogram;
		}

	}

	public static class BlastWorker implements Callable<BlastResultsSummary[]> {
		private final Blast blast;
		private final FastaEntry[] fastaEntries;
		private final String tmpFile;

		public BlastWorker(Blast blast, FastaEntry[] fastaEntries, String tmpFile) {
			super();
			this.blast = blast;
			this.fastaEntries = fastaEntries;
			this.tmpFile = tmpFile;
		}

		@Override
		public BlastResultsSummary[] call() throws Exception {
			PrintWriter writer = null;
			if (tmpFile != null) {
				blast.getLog().reportTimeInfo("Output sent to " + tmpFile);
				writer = Files.getAppropriateWriter(tmpFile);
				writer.println(ArrayUtils.toStr(BLAST_HEADER));

			}

			BlastResultsSummary[] blasts = blast.blastSequence(fastaEntries, writer);
			writer.close();
			return blasts;
		}

	}



	public static void test() {
		String fastaDb = "/home/pankrat2/public/bin/ref/hg19_canonical.fa";
		Blast blast = new Blast(fastaDb, 60, 100, new Logger(), true, true);
		FastaEntry fastaEntry =
													new FastaEntry(	"HDSIF",
																					"GAGCCGGAGCACCCTATGTCGCAGTATCTGTCTTTGATTCCTGCCTCATTCTATTATTTATCGCACCTACGTTCAATATTACAGGCGAACATACCTACTAAAGTGTGTTAATTAATTAATGCTTGTAGGACATAATAATAACAATTGAAT");
		blast.blastSequence(new FastaEntry[] {fastaEntry}, null);
	}

	public static void main(String[] args) {
		int numArgs = args.length;

		String filename = "Blast.dat";
		// String logfile = null;
		// Logger log;

		String usage = "\n"	+ "seq.analysis.Blast requires 0-1 arguments\n"
										+ "   (1) filename (i.e. file=" + filename + " (default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("file=")) {
				filename = arg.split("=")[1];
				numArgs--;
			}
			// else if (args[i].startsWith("log=")) {
			// logfile = args[i].split("=")[1];
			// numArgs--;
			// }
			else {
				System.err.println("Error - invalid argument: " + arg);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			// log = new Logger(logfile);
			test();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
//
// else if ((btop.charAt(i) + "").equals("-")) {//gap
// if (currentString != null) {
// btopBroken.add(currentString);
// currentString = null;
// }
// if (currentInt != null) {
// btopBroken.add(currentInt);
// currentInt = null;
// }
// btopBroken.add(btop.charAt(i) + "");
// }
