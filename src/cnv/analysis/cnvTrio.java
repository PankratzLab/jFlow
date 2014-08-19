package cnv.analysis;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import cnv.filesys.Sample;
import cnv.var.CNVariant;
import cnv.var.SampleData;
import common.Array;
import common.Files;
import common.Logger;
import common.Matrix;
import common.ext;

/**
 * Class for filtering denovo calls in offspring by BEAST SCORE and LRR_SD, and a few other metrics Filtering at the default metrics here seems to work OK
 */
//TODO, exclude from sample data
public class cnvTrio {
	private static final String[] SPLITS = { "\t" };
	private static final String[] COMBINED_TRIOS = { "denovo_joint.cnv", "denovo_trio.cnv" };
	private static final String[] TRIOS = { ".jointcnv", ".triocnv" };
	private static final String[] BEAST_OUTPUT = { ".beast.summary", ".beast.cnv", ".filtered" };
	private static final String[] BEAST_HEADER = { "NUM_MARKERS", "IDNA", "ILRR_SD", "IBAF1585_SD", "IHEIGHT", "ILength", "IBEAST", "FADNA", "FALRR_SD", "FABAF1585_SD", "FAHEIGHT", "FALength", "FABEAST", "MODNA", "MOLRR_SD", "MOBAF1585_SD", "MOHEIGHT", "MOLength", "MOBEAST", "MinBeastDiff", "MaxTrioLRR_SD", "MaxTrioBAF1585_SD", "NumDeNovo", "MinBeastDiffQC", "MaxBeastParentQC", "MaxTrioLRR_SDQC", "MaxTrioBAF1585_SDQC", "NumDeNovoQC", "ALLQC", "UCSC", "UCSC_HYPERLINK" };
	public static final float DEFAULT_BEAST_HEIGHT_DIFF_THRESHOLD = 0.25f;
	public static final float DEFAULT_BEAST_HEIGHT_PARENTS_THRESHOLD = 0.20f;
	public static final float DEFAULT_MAXTRIO_LRR_SD_THRESHOLD = 0.30f;
	public static final int DEFAULT_MAXNUM_CALLS_THRESHOLD = 50;

	/**
	 * Steps are:
	 * <p>
	 * 1. load the trios to Trio[]
	 * <p>
	 * 2. Parse the penncnv results and collect in one .cnv file using DeNovoCNV.parsePennCnvResult
	 * <p>
	 * 3. filter with FilterCalls.filter ( only for numProbes and problematic regions)
	 * <p>
	 * 4. assign variants to each trio in Trio[]
	 * <p>
	 * 5. process trios (compute beast and other qc metrics), multithreaded
	 * <p>
	 * 6. summarize the trio calls (calls removed in FilterCalls.filter are not reported)
	 * <p>
	 * Optional: 7. update the projects region and individual list for ease of viewing filtered calls
	 */
	public static void analyze(Project proj, String trioFile, int fileType, Filter filter, String filenameOfProblematicRegions, boolean updateRegionAndList, int numThreads) {
		Logger log = proj.getLog();
		Trio[] trios = Trio.loadTrios(proj, trioFile);
		DeNovoCNV.parsePennCnvResult(proj, proj.getDir(Project.PENNCNV_RESULTS_DIRECTORY), proj.getDir(Project.DATA_DIRECTORY) + trioFile, TRIOS[fileType]);
		String[] no = null;
		FilterCalls.filter(proj.getDir(Project.DATA_DIRECTORY), COMBINED_TRIOS[fileType], COMBINED_TRIOS[fileType] + BEAST_OUTPUT[2], 0, 0, filter.getNumProbes(), 0, filenameOfProblematicRegions, 3, no, false, null, false, 36, log);
		CNVariant[] cnVariants = CNVariant.loadPlinkFile(proj.getDir(Project.DATA_DIRECTORY) + COMBINED_TRIOS[fileType] + BEAST_OUTPUT[2], false);

		parseTrios(cnVariants, trios, proj.getSampleData(0, false), log);
		processTrios(proj, trios, filter, numThreads);
		summarizeTrios(proj, trios, fileType);

		if (updateRegionAndList) {
			writeRegionList(proj, proj.getDir(Project.DATA_DIRECTORY) + COMBINED_TRIOS[fileType] + BEAST_OUTPUT[1]);
		}
	}

	/**
	 * Assigns offspring CNV calls to trios by matching through sampleData
	 * <p>
	 * Essentially parses the array of CNVariants and assigns to the appropriate trio
	 * <p>
	 * If FID/IIDs are not unique (correspond to more than one DNA, the CNV will be dropped due to the ambiguity)
	 * 
	 * @param cnVariants
	 * @param trios
	 * @param sampleData
	 * @param log
	 */

	// TODO, resolve the issue of non-unique FID/IIDs (-> might have to be done manually), this is a problem if there are identical FID/IID combos, and the DNA happens to match. We will not catch it.
	private static void parseTrios(CNVariant[] cnVariants, Trio[] trios, SampleData sampleData, Logger log) {
		Hashtable<String, Integer> trioIndex = hashTrios(trios, sampleData);
		int count = 0;
		for (int i = 0; i < cnVariants.length; i++) {
			String key = cnVariants[i].getFamilyID() + "\t" + cnVariants[i].getIndividualID();
			if (trioIndex.containsKey(key)) {
				if (trios[trioIndex.get(key)].addCNVTrio(sampleData, cnVariants[i], log)) {
					count++;
				}
			} else {
				log.reportError("Warning - the cnv " + cnVariants[i].toPlinkFormat() + " was not defined in the as an offspring in the trio file \n Skipping...");
			}
		}
		log.report(ext.getTime() + " Found " + count + " offspring cnvs");
		if (count != cnVariants.length) {
			log.reportError("Warning - Due to duplicate FID/IID combos or mismatched ids, only " + count + "/" + cnVariants.length + " cnvs are being analyzed");
		}
	}

	/**
	 * To speed things up a bit, we send out the heavy beast lifting to other threads
	 */
	public static void processTrios(Project proj, Trio[] trios, Filter filter, int numThreads) {
		MarkerSet markerSet = proj.getMarkerSet();
		int[][] indi = markerSet.getIndicesByChr();
		ExecutorService executor = Executors.newFixedThreadPool(numThreads);
		for (int i = 0; i < trios.length; i++) {
			Runnable worker = new WorkerThread(proj, trios[i], markerSet.getChrs(), markerSet.getPositions(), Matrix.clone(indi), filter, i);
			executor.execute(worker);
		}
		executor.shutdown();
		try {
			executor.awaitTermination(7, TimeUnit.DAYS);
		} catch (InterruptedException e) {
		}
		proj.getLog().report(ext.getTime() + "Info - completed beast score computation");
	}

	/**
	 * WorkerThreads which processes each trio to compute beast scores and qc metrics
	 */

	public static class WorkerThread implements Runnable {
		private Project proj;
		private Trio trio;
		private byte[] chrs;
		private int[] positions;
		private int[][] indicesByChr;
		private Filter filter;
//		private Logger log;
		private int threadID;

		public WorkerThread(Project proj, Trio trio, byte[] chr, int[] positions, int[][] indicesByChr, Filter filter, int threadID) {
			super();
			this.proj = proj;
			this.trio = trio;
			this.chrs = chr;
			this.positions = positions;
			this.indicesByChr = indicesByChr;
			this.filter = filter;
			this.threadID = threadID;
		}

		@Override
		public void run() {
			Logger log = proj.getLog();
			
			log.report(ext.getTime() + " Info - trio (" + (threadID + 1) + ") " + trio.getTrio());
			trio.computeTrioBeast(proj, chrs, positions, indicesByChr, filter);
		}
	}

	/**
	 * Creates two output files, one with all the qc metrics, and another with just the list of filtered trio calls
	 * 
	 * @param proj
	 * @param trios
	 * @param fileType
	 * @param log
	 */
	private static void summarizeTrios(Project proj, Trio[] trios, int fileType) {
		backup(proj, fileType);
		PrintWriter writerFullSummary = Files.getAppropriateWriter(proj.getDir(Project.DATA_DIRECTORY) + COMBINED_TRIOS[fileType] + BEAST_OUTPUT[0]);
		PrintWriter writerFiltered = Files.getAppropriateWriter(proj.getDir(Project.DATA_DIRECTORY) + COMBINED_TRIOS[fileType] + BEAST_OUTPUT[1]);
		writerFullSummary.print(Array.toStr(CNVariant.PLINK_CNV_HEADER) + "\t" + Array.toStr(BEAST_HEADER) + "\n");
		writerFiltered.println(Array.toStr(CNVariant.PLINK_CNV_HEADER));
		for (int i = 0; i < trios.length; i++) {
			if (trios[i].hasSummary()) {
				writerFullSummary.println(Array.toStr(trios[i].getFullSummary(), "\n"));
				if (trios[i].hasFilteredCNVS()) {
					writerFiltered.println(Array.toStr(trios[i].getFilteredCNVS(), "\n"));
				}
			}
		}
		writerFullSummary.close();
		writerFiltered.close();
	}

	private static void backup(Project proj, int fileType) {
		if (Files.exists(proj.getDir(Project.DATA_DIRECTORY) + COMBINED_TRIOS[fileType] + BEAST_OUTPUT[0])) {
			Files.backup(COMBINED_TRIOS[fileType] + BEAST_OUTPUT[0], proj.getDir(Project.DATA_DIRECTORY), proj.getDir(Project.BACKUP_DIRECTORY));
		}
		if (Files.exists(proj.getDir(Project.DATA_DIRECTORY) + COMBINED_TRIOS[fileType] + BEAST_OUTPUT[1])) {
			Files.backup(COMBINED_TRIOS[fileType] + BEAST_OUTPUT[1], proj.getDir(Project.DATA_DIRECTORY), proj.getDir(Project.BACKUP_DIRECTORY));
		}
	}

	/**
	 * Assign the offspring FID/IIDs to a hash so we know which trio to add a particular CNVariant to
	 * 
	 * @param trios
	 * @param sampleData
	 * @return
	 */
	private static Hashtable<String, Integer> hashTrios(Trio[] trios, SampleData sampleData) {
		Hashtable<String, Integer> trioIndex = new Hashtable<String, Integer>();
		for (int i = 0; i < trios.length; i++) {
			String key = trios[i].getOffspringFIDIID();
			trioIndex.put(key, i);
		}
		return trioIndex;
	}

	/**
	 * Helper class to hold trio related information, and handle the computation/summarization
	 * 
	 */
	public static class Trio {
		private static final String[] PED_TRIO_HEADER = { "fId", "iId", "faId", "moId", "iDna", "faDna", "moDna" };
		private String fID;
		private String iID;
		private String faID;
		private String moID;
		private String iDNA;
		private String faDNA;
		private String moDNA;
		private ArrayList<CNVariant> ICNV;
		private ArrayList<String> fullSummary;
		private ArrayList<String> filteredCNVS;

		public Trio(String fID, String iID, String faID, String moID, String iDNA, String faDNA, String moDNA) {
			super();
			this.fID = fID;
			this.iID = iID;
			this.faID = faID;
			this.moID = moID;
			this.iDNA = iDNA;
			this.faDNA = faDNA;
			this.moDNA = moDNA;
			this.ICNV = new ArrayList<CNVariant>();
			this.fullSummary = new ArrayList<String>();
			this.filteredCNVS = new ArrayList<String>();
		}

		public String getIDNA() {
			return iDNA;
		}

		public String getFADNA() {
			return faDNA;
		}

		public String getMODNA() {
			return moDNA;
		}

		public String getOffspringFIDIID() {
			return fID + "\t" + iID;
		}

		public String getTrio() {
			return getIDNA() + "\t" + getFADNA() + "\t" + getMODNA();
		}

		/**
		 * @return the full summary for all calls passed to the trio
		 */
		public String[] getFullSummary() {
			return fullSummary.toArray(new String[fullSummary.size()]);
		}

		/**
		 * @return the filtered cnv calls
		 */
		public String[] getFilteredCNVS() {
			return filteredCNVS.toArray(new String[filteredCNVS.size()]);
		}

		public boolean hasSummary() {
			return fullSummary.size() > 0;
		}

		public boolean hasFilteredCNVS() {
			return filteredCNVS.size() > 0;
		}

		/***
		 * The main event, grabs LRR_SD and BeastScore and few others. Summarizes QC across trio... max LRR_SD across the trio and minimum beast score (height) distance from offspring to parent Writes to two files: writerFullSummary the full summary of the beast-> writerFullSummary writerFiltered a new filtered cnv file with passing calls, and calls assigned to each of the parents for ease of vis in comp plot->trailer
		 */
		public void computeTrioBeast(Project proj, byte[] chrs, int[] positions, int[][] indicesByChr, Filter filter) {
			CNVariant[] analyzeCNV = ICNV.toArray(new CNVariant[ICNV.size()]);
			Logger log = proj.getLog();
			
			if (analyzeCNV.length == 0) {
				log.reportError("Warning - trio " + getTrio() + " does not have any CNVs defined");
			} else {
				int[][] targetIndices = getCNVIndices(chrs, positions, analyzeCNV, log);
				Sample iSamp = proj.getFullSampleFromRandomAccessFile(iDNA);
				Sample faSamp = proj.getFullSampleFromRandomAccessFile(faDNA);
				Sample moSamp = proj.getFullSampleFromRandomAccessFile(moDNA);

				float iBAF1585SD = getBAF1585SD(iSamp.getBAFs(), log);
				float faBAF1585SD = getBAF1585SD(faSamp.getBAFs(), log);
				float moBAF1585SD = getBAF1585SD(moSamp.getBAFs(), log);
				float maxBAF1585SD = Math.max(Math.max(iBAF1585SD, faBAF1585SD), moBAF1585SD);

				float[] iLrr = iSamp.getLRRs();
				float[] faLrr = faSamp.getLRRs();
				float[] moLrr = moSamp.getLRRs();

				float iStDev = Array.stdev(Array.subArray(iLrr, 0, Array.indexOfByte(chrs, (byte) 23)), true);
				float faStDev = Array.stdev(Array.subArray(faLrr, 0, Array.indexOfByte(chrs, (byte) 23)), true);
				float moStDev = Array.stdev(Array.subArray(moLrr, 0, Array.indexOfByte(chrs, (byte) 23)), true);
				float maxStDev = Math.max(Math.max(iStDev, faStDev), moStDev);

				BeastScore iBeast = getBeast(iLrr, indicesByChr, targetIndices, log);
				BeastScore faBeast = getBeast(faLrr, indicesByChr, targetIndices, log);
				BeastScore moBeast = getBeast(moLrr, indicesByChr, targetIndices, log);

				for (int i = 0; i < analyzeCNV.length; i++) {
					float minBeastDiff = minHeightDist(iBeast.getBeastHeights()[i], faBeast.getBeastHeights()[i], moBeast.getBeastHeights()[i], filter);

					String fSum = "";
					fSum += analyzeCNV[i].toPlinkFormat() + "\t" + targetIndices[i].length;
					fSum += "\t" + iDNA + "\t" + iStDev + "\t" + iBAF1585SD + "\t" + iBeast.getSummaryAt(i);
					fSum += "\t" + faDNA + "\t" + faStDev + "\t" + faBAF1585SD + "\t" + faBeast.getSummaryAt(i);
					fSum += "\t" + moDNA + "\t" + moStDev + "\t" + moBAF1585SD + "\t" + moBeast.getSummaryAt(i);
					fSum += "\t" + minBeastDiff + "\t" + maxStDev + "\t" + maxBAF1585SD + "\t" + analyzeCNV.length;
					fSum += "\t" + getQCString(faBeast.getBeastHeights()[i], moBeast.getBeastHeights()[i], minBeastDiff, maxStDev, maxBAF1585SD, analyzeCNV.length, filter);
					fSum += "\t" + analyzeCNV[i].getUCSClocation() + "\t" + analyzeCNV[i].getUCSCLink();

					fullSummary.add(fSum);

					if (checkAll(faBeast.getBeastHeights()[i], moBeast.getBeastHeights()[i], analyzeCNV.length, maxStDev, maxBAF1585SD, minBeastDiff, filter)) {
						filteredCNVS.add(analyzeCNV[i].toPlinkFormat());
						filteredCNVS.add(getVariant(fID, faID, analyzeCNV[i]).toPlinkFormat());
						filteredCNVS.add(getVariant(fID, moID, analyzeCNV[i]).toPlinkFormat());
					}
				}
			}
		}

		private BeastScore getBeast(float[] lrrs, int[][] indicesByChr, int[][] targetIndices, Logger log) {
			BeastScore iBeast = new BeastScore(lrrs, indicesByChr, targetIndices, log);
			iBeast.computeBeastScores(BeastScore.DEFAULT_ALPHA);
			return iBeast;
		}

		private static float getBAF1585SD(float[] bafs, Logger log) {
			for (int j = 0; j < bafs.length; j++) {
				if (bafs[j] < 0.15 || bafs[j] > 0.85) {
					bafs[j] = Float.NaN;
				}
			}
			return Array.stdev(bafs, true);
		}

		private static boolean checkAll(float faHeight, float moHeight, int numCalls, float maxStDev, float maxBAF1585SD, float minBeastDiff, Filter filter) {
			return faHeight <= filter.getMaxBeastHeightParents() && moHeight <= filter.getMaxBeastHeightParents() && minBeastDiff >= filter.getMinBeastHeightDiffThreshold() && checkLrrSdThreshold(maxStDev, filter) && numCalls <= filter.getMaxNumCallsThreshold() && checkBAF1585SD(maxBAF1585SD, filter);
		}

		/**
		 * 
		 * @return a string with 1s for good and 0s for bad
		 */
		private static String getQCString(float faHeight, float moHeight, float minBeastDiff, float maxStDev, float maxBAF1585SD, int numCalls, Filter filter) {
			String qcString = "";
			if (minBeastDiff >= filter.getMinBeastHeightDiffThreshold()) {
				qcString += "1\t";
			} else {
				qcString += "0\t";
			}
			if (faHeight <= filter.getMaxBeastHeightParents() && moHeight <= filter.getMaxBeastHeightParents()) {
				qcString += "1\t";
			} else {
				qcString += "0\t";
			}
			if (checkLrrSdThreshold(maxStDev, filter)) {
				qcString += "1\t";
			} else {
				qcString += "0\t";
			}
			if (checkBAF1585SD(maxStDev, filter)) {
				qcString += "1\t";
			} else {
				qcString += "0\t";
			}
			if (numCalls <= filter.getMaxNumCallsThreshold()) {
				qcString += "1\t";
			} else {
				qcString += "0\t";
			}
			if (checkAll(faHeight, moHeight, numCalls, maxStDev, maxBAF1585SD, minBeastDiff, filter)) {
				qcString += "1";
			} else {
				qcString += "0";
			}
			return qcString;
		}

		private static boolean checkBAF1585SD(float maxStDev, Filter filter) {
			return maxStDev <= filter.getMaxBAF1585SD();
		}

		private static boolean checkLrrSdThreshold(float maxStDev, Filter filter) {
			return maxStDev <= filter.getMaxLrrSdThreshold();
		}

		/**
		 * If the parents have a height that exceeds the realCNVHeight (MaxBeastHeightParents), we defualt to 0; If the offspring's height does not surpass the realCNVHeight (MaxBeastHeightParents), we defualt to 0; <p>
		 * <p>
		 * Also check if the parents have in opposite signed height as the offspring, if so, we set the parent's height to 0.
		 * Otherwise, we return the minimum distance between parents and offspring
		 * 
		 * @param offSpringHeight
		 *            beast height of the offspring
		 * @param faHeigt
		 *            father beast height
		 * @param moHeight
		 *            mother beast height
		 * @param realCNVHeight
		 *            definition of height for a cnv
		 * @return
		 */
		private static float minHeightDist(float offSpringHeight, float faHeigt, float moHeight, Filter filter) {
			float minHeightDist = 0;
			if (Math.abs(faHeigt) > filter.getMaxBeastHeightParents() || Math.abs(moHeight) > filter.getMaxBeastHeightParents() || Math.abs(offSpringHeight) < filter.getMaxBeastHeightParents()) {
				return minHeightDist;
			}
			faHeigt = checkParentOppositeHeight(offSpringHeight, faHeigt);
			moHeight = checkParentOppositeHeight(offSpringHeight, moHeight);
			return Math.min(Math.abs(offSpringHeight - faHeigt), Math.abs(offSpringHeight - moHeight));
		}

		/** if offSpringHeight and parentHeight are different signs, return 0 for parentHeight
		 * Note, checks for a real cnv should be done prior to using this function
		 * @param offSpringHeight 
		 * @param parentHeight
		 * @return
		 */
		private static float checkParentOppositeHeight(float offSpringHeight , float parentHeight){
			if(offSpringHeight*parentHeight >=0){
				return parentHeight;
			}else{
				return 0;
			}
		}
		/**
		 * Add the variant if lookup matches
		 */
		private boolean addCNVTrio(SampleData sampleData, CNVariant cnVariant, Logger log) {
			if (checkValid(sampleData, cnVariant)) {
				ICNV.add(cnVariant);
				return true;
			} else {
				log.reportError("\nError - The current offspring has DNA: " + iDNA + ", FID: " + fID + ", and IID: " + iID + ", the current cnv call has FID: " + cnVariant.getFamilyID() + " and IID: " + cnVariant.getIndividualID());
				log.reportError("The lookup of the cnv FID/IID combo returned a mismatched DNA: " + sampleData.lookup(cnVariant.getFamilyID() + "\t" + cnVariant.getIndividualID())[0]);
				log.reportError("This can happen if there are identical FID/IID combos in the file coorresponding to different DNA IDs. If so, please make the FID/IID combos unique");
				log.reportError("Skipping CNV " + cnVariant.toPlinkFormat() + " due to unknown source DNA\n");
				return false;
			}
		}

		/**
		 * Make sure the FID/IID lookup returns the correct DNA
		 * 
		 * @param sampleData
		 * @param cnVariant
		 * @return
		 */
		private boolean checkValid(SampleData sampleData, CNVariant cnVariant) {
			String lookup = sampleData.lookup(cnVariant.getFamilyID() + "\t" + cnVariant.getIndividualID())[0];
			String FIDIID = fID + "\t" + iID;
			return iDNA.equals(lookup) && FIDIID.equals(cnVariant.getFamilyID() + "\t" + cnVariant.getIndividualID());
		}

		private static CNVariant getVariant(String FID, String IID, CNVariant cnVariant) {
			return new CNVariant(FID, IID, cnVariant.getChr(), cnVariant.getStart(), cnVariant.getStop(), cnVariant.getCN(), cnVariant.getScore(), cnVariant.getNumMarkers(), cnVariant.getSource());
		}

		/**
		 * @param proj
		 * @param trioFile
		 *            should have header cnvTrio.PED_TRIO_HEADER
		 * @param log
		 * @return
		 */
		public static Trio[] loadTrios(Project proj, String trioFile) {
			ArrayList<Trio> alTrio = new ArrayList<Trio>();
			int trioIndex = 0;
			Logger log = proj.getLog();
			
			try {
				BufferedReader reader = Files.getReader(proj.getDir(Project.DATA_DIRECTORY) + trioFile, false, true, false);
				String[] line = reader.readLine().trim().split(SPLITS[0]);
				int[] indices = ext.indexFactors(PED_TRIO_HEADER, line, true, true);
				while (reader.ready()) {
					line = reader.readLine().trim().split(SPLITS[0]);
					alTrio.add(new Trio(line[indices[0]], line[indices[1]], line[indices[2]], line[indices[3]], line[indices[4]], line[indices[5]], line[indices[6]]));
					trioIndex++;
				}
				reader.close();
			} catch (FileNotFoundException fnfe) {
				log.reportError("Error: file \"" + proj.getProjectDir() + trioFile + "\" not found in current directory");
			} catch (IOException ioe) {
				log.reportError("Error reading file \"" + proj.getProjectDir() + trioFile + "\"");
			}
			if (trioIndex == 0) {
				log.reportError("Error - did not find any trios in " + proj.getDir(Project.DATA_DIRECTORY) + trioFile);
			}
			log.report(ext.getTime() + " Info - found " + (trioIndex > 1 ? trioIndex + " trios" : trioIndex + " trio"));
			return alTrio.toArray(new Trio[alTrio.size()]);
		}
	}

	/**
	 * Get the indices for the markers contained in the cnv
	 * 
	 * @param chr
	 * @param positions
	 * @param cnVariants
	 * @param log
	 * @return
	 */
	private static int[][] getCNVIndices(byte[] chr, int[] positions, CNVariant[] cnVariants, Logger log) {
		int[][] indices = new int[cnVariants.length][];
		for (int i = 0; i < cnVariants.length; i++) {
			indices[i] = BeastScore.getCNVMarkerIndices(chr, positions, cnVariants[i], log);
		}
		return indices;
	}

	/**
	 * @param proj
	 * @param cnvFile
	 *            converts this cnv file to Region list, Individual CNV list, and a bed type file
	 * @param log
	 */
	private static void writeRegionList(Project proj, String cnvFile) {
		CNVariant[] cnvs = CNVariant.loadPlinkFile(cnvFile, false);
		proj.getLog().report("Info - updating regions and list using " + cnvs.length + " cnvs");
		String list = proj.getFilename(Project.INDIVIDUAL_CNV_LIST_FILENAMES).replaceAll(";", "");
		String regionList = proj.getFilename(Project.REGION_LIST_FILENAMES).replaceAll(";", "");
		ArrayList<String> lists = new ArrayList<String>();
		ArrayList<String> bed = new ArrayList<String>();
		ArrayList<String> regions = new ArrayList<String>();
		SampleData sampleData = proj.getSampleData(0, false);
		Hashtable<String, String> track = new Hashtable<String, String>();
		for (int i = 0; i < cnvs.length; i++) {
			String UCSC = cnvs[i].getUCSClocation();
			if (!track.containsKey(UCSC)) {
				bed.add("chr" + cnvs[i].getChr() + "\t" + cnvs[i].getStart() + "\t" + cnvs[i].getStop() + "\t" + cnvs[i].getUCSClocation());
				regions.add(UCSC);
				track.put(UCSC, UCSC);
			}

			lists.add(sampleData.lookup(cnvs[i].getFamilyID() + "\t" + cnvs[i].getIndividualID())[0] + "\t" + UCSC);
		}
		System.out.println(list + "\n" + regionList);
		Files.writeList(lists.toArray(new String[lists.size()]), list);
		Files.writeList(regions.toArray(new String[regions.size()]), regionList);
		Files.writeList(bed.toArray(new String[bed.size()]), regionList + ".bed");
	}

	/**
	 * A small helper class to facilitate passing thresholds around to anyone who wants em
	 * 
	 */
	private static class Filter {
		private float minBeastHeightDiffThreshold;
		private float maxLrrSdThreshold;
		private float maxBeastHeightParents;
		private int maxNumCallsThreshold;
		private int numProbes;
		private float maxBAF1585SD;

		/**
		 * @param minBeastHeightDiffThreshold
		 *            minimum height difference between parents and offspring (also defines what height is a cnv)
		 * @param maxLrrSdThreshold
		 *            max LRR_SD across the trio
		 * @param maxBeastHeightParents
		 *            maximum height of the parents (i.e height of a real cnv)
		 * @param maxNumCallsThreshold
		 *            maximum denovo calls allowed
		 * @param numProbes
		 *            minimum number of probes for a call
		 * @param maxBAF1585SD
		 *            maximum BAF_1585_SD across the trio
		 */
		public Filter(float minBeastHeightDiffThreshold, float maxLrrSdThreshold, float maxBeastHeightParents, int maxNumCallsThreshold, int numProbes, float maxBAF1585SD) {
			super();
			this.minBeastHeightDiffThreshold = minBeastHeightDiffThreshold;
			this.maxLrrSdThreshold = maxLrrSdThreshold;
			this.maxBeastHeightParents = maxBeastHeightParents;
			this.maxNumCallsThreshold = maxNumCallsThreshold;
			this.numProbes = numProbes;
			this.maxBAF1585SD = maxBAF1585SD;
		}

		public float getMinBeastHeightDiffThreshold() {
			return minBeastHeightDiffThreshold;
		}

		public float getMaxLrrSdThreshold() {
			return maxLrrSdThreshold;
		}

		public float getMaxBeastHeightParents() {
			return maxBeastHeightParents;
		}

		public int getMaxNumCallsThreshold() {
			return maxNumCallsThreshold;
		}

		public int getNumProbes() {
			return numProbes;
		}

		public float getMaxBAF1585SD() {
			return maxBAF1585SD;
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "C:/workspace/Genvisis/projects/OSv2.properties";
		String trioFile = "PedigreeOfTrios.txt";
		String logFile = "D:/data/logan/OSv2/trioLog.log";
		String filenameOfProblematicRegions = "D:/data/logan/OSv2/data/problematicRegions_hg18.dat";
		float minBeastHeightDiffThreshold = .5f;
		float maxLrrSdThreshold = DEFAULT_MAXTRIO_LRR_SD_THRESHOLD;
		float maxBeastHeightParents = DEFAULT_BEAST_HEIGHT_PARENTS_THRESHOLD;
		float maxBAF1585SD = 0.8f;
		int maxNumCallsThreshold = DEFAULT_MAXNUM_CALLS_THRESHOLD;
		int numProbes = 5;
		int fileType = 1;
		// TODO, change default numThreads and Region list to false
		int numThreads = 7;
		boolean updateRegionAndList = true;

		String usage = "\n";
		usage += "jlDev.cnvTrio requires 0-13 arguments\n";
		usage += "   (1) project file (i.e. proj=" + filename + " (default))\n";
		usage += "   (2) a file listing pedgiree information for the trios (i.e. trioFile=" + trioFile + " (default))\n";
		usage += "   (3) the type of penncnv trio results to analyze (i.e. fileType=" + TRIOS[fileType] + " (default)), options are " + Array.toStr(TRIOS) + "\n";
		usage += "   (4) the minimum height score difference between offspring and parents (i.e. minBeastThreshold=" + minBeastHeightDiffThreshold + " (default))\n";
		usage += "   (5) the maximum height score for parents  (i.e. maxBeastParents=" + maxBeastHeightParents + " (default))\n";
		usage += "   (6) maximum log R Ratio standard deviation across the trio (i.e. maxLrrSdThreshold=" + maxLrrSdThreshold + " (default))\n";
		usage += "   (7) maximum BAF1585 SD (i.e. maxBAF1585SD=" + maxBAF1585SD + " (default))\n";
		usage += "   (8) maximum number of denovo cnvs in a sample (i.e. maxNumCallsThreshold=" + maxNumCallsThreshold + " (default))\n";
		usage += "   (9) minimum number of probes for a cnv call; this is a hard filter -> cnvs not passing this threshold will not be reported or analyzed (i.e. numProbes=" + numProbes + " (default))\n";
		usage += "   (10) filter out cnvs in known problematicRegions; this is a hard filter -> cnvs in problematic regions will not be reported or analyzed (i.e. filterFile=" + filenameOfProblematicRegions + " (default))\n";
		usage += "   (11) update the cnv lists for the project based on filtered and parsed results ( i.e. -updateRegionAndList (not the default))\n";
		usage += "   (12) log filename (i.e. logFile=" + logFile + " (default))\n";
		usage += "   (13) number of threads to use for the analysis (i.e. numThreads=" + numThreads + " ( no default))\n";
		usage += "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("filename=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("trioFile=")) {
				trioFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("fileType=")) {
				fileType = ext.indexOfStr(args[i].split("=")[1], TRIOS);
				numArgs--;
			} else if (args[i].startsWith("logFile=")) {
				logFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("filterFile=")) {
				filenameOfProblematicRegions = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("minBeastThreshold=")) {
				minBeastHeightDiffThreshold = Float.parseFloat(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("maxLrrSdThreshold=")) {
				maxLrrSdThreshold = Float.parseFloat(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("maxBeastParents=")) {
				maxBeastHeightParents = Float.parseFloat(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("maxBAF1585SD=")) {
				maxBAF1585SD = Float.parseFloat(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("maxNumCallsThreshold=")) {
				maxNumCallsThreshold = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("numProbes=")) {
				numProbes = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("numThreads=")) {
				numThreads = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("-updateRegionAndList=")) {
				updateRegionAndList = true;
				numArgs--;
			} else {
				System.out.println("Invalid argument " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		Project proj = new Project(filename, logFile, false);
		Filter filter = new Filter(minBeastHeightDiffThreshold, maxLrrSdThreshold, maxBeastHeightParents, maxNumCallsThreshold, numProbes, maxBAF1585SD);
		analyze(proj, trioFile, fileType, filter, filenameOfProblematicRegions, updateRegionAndList, numThreads);
	}
}
