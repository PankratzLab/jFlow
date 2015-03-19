package cnv.qc;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Vector;

import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import cnv.var.CNVariant;
import cnv.var.SampleData;
import common.Array;
import common.Files;
import common.HashVec;
import common.Logger;
import common.Positions;
import common.ext;
import filesys.Segment;
import filesys.SnpMarkerSet;

public class CNVFilter {

	/**
	 */
	public static final int GIANT_BP = 10000000;
	public static final int GIANT_NUM_MARKERS = 500;
	public static final int DEFAULT_MIN_NUM_MARKERS = 5;
	public static final double DEFAULT_MIN_SCORE_THRESHOLD = 15;
	public static final boolean DEFAULT_BREAK_UP_CENTROMERES = false;
	public static final boolean DEFAULT_EXCLUDE_SAMPLE_DATA = false;
	public static final boolean DEFAULT_COMMON_IN = false;
	public static final int DEFAULT_BUILD = 36;

	/**
	 * COMMAND_* variables are for building .crf parameters, and can be used for command line usage
	 */
	public static final String COMMAND_MIN_NUM_MARKERS = "minNumMarkers=";//
	public static final String COMMAND_MAX_NUM_MARKERS = "maxNumMarkers=";//
	public static final String COMMAND_MIN_SIZE = "minSize=";//
	public static final String COMMAND_MAX_SIZE = "maxSize=";//
	public static final String COMMAND_MIN_SCORE = "minScore=";//
	public static final String COMMAND_MAX_SCORE = "maxScore=";//
	public static final String COMMAND_COMMON_IN = "commonIn=";//
	public static final String COMMAND_PROBLEM_REGIONS_FILE = "problematicRegionsFile=";//
	public static final String COMMAND_BREAK_UP_CENTROMERES = "breakupCentromeres=";
	public static final String COMMAND_BUILD = "build=";
	public static final String COMMAND_COMMON_REFERENCE = "commonReferenceFile=";
	public static final String COMMAND_INDIVIDUALS_TO_KEEP = "individualsToKeepFile=";
	public static final String COMMAND_EXCLUDE_INDIVIDUALS_FROM_SAMPLE_DATA = "exclude=";
	public static final String COMMAND_PROJECT = "proj=";
	public static final String COMMAND_CNV_FILE = "cnvFile=";
	public static final String COMMAND_CNV_FILE_OUT = "out=";
	public static final String COMMAND_CNV_FILTER = "-filter";

	public static final String COMMAND_CNV_FILTER_CRF = "cnvFilter";
	public static final String COMMAND_CNV_FILTER_DESCRIPTION = "filter a file of cnvs";

	/**
	 * Thes parameters should result in no filtering
	 */
	public static final int NO_FILTER_CN = -1;
	public static final int NO_FILTER_MIN_NUM_MARKERS = -1;
	public static final int NO_FILTER_MAX_NUM_MARKERS = Integer.MAX_VALUE;
	public static final double NO_FILTER_MIN_SCORE = -1;
	public static final double NO_FILTER_MAX_SCORE = Double.MAX_VALUE;
	public static final int NO_FILTER_MIN_SIZE = -1;
	public static final int NO_FILTER_MAX_SIZE = Integer.MAX_VALUE;
	public static final Segment[] NO_FILTER_PROBLEM_REGIONS = new Segment[0];
	public static final Segment[] NO_FILTER_CENTROMERE_MIDPOINTS = new Segment[0];
	public static final Segment[] NO_FILTER_COMMON_REFERENCE = new Segment[0];
	public static final int[][] NO_FILTER_CENTROMERE_BOUNDARIES = null;
	public static final boolean NO_FILTER_BREAK_UP_CENTROMERES = false;
	public static final HashSet<String> NO_FILTER_INDIVIDUAL_HASH = null;

	private int minNumMarkers, maxNumMarkers, minSize, maxSize, build, CN;
	private double minScore, maxScore;
	private Segment[] problemRegions, centromereMidpoints, commonReference;
	private int[][] centromereBoundaries;
	private int[][] positions;
	private boolean breakupCentromeres, commonIn;
	private HashSet<String> indHash;
	private Hashtable<String, String> commandLineFiltersInEffect = new Hashtable<String, String>();
	private Logger log;

	public CNVFilter(int minNumMarkers, int maxNumMarkers, int minSize, int maxSize, double minScore, double maxScore, Segment[] problemRegions, Segment[] centromereMidpoints, Segment[] commonReference, int[][] centromereBoundaries, boolean breakupCentromeres, boolean commonIn, HashSet<String> indHash, int build, int CN, Logger log) {
		super();
		this.minNumMarkers = minNumMarkers;
		this.maxNumMarkers = maxNumMarkers;
		this.minSize = minSize;
		this.maxSize = maxSize;
		this.minScore = minScore;
		this.maxScore = maxScore;
		this.problemRegions = problemRegions;
		this.centromereMidpoints = centromereMidpoints;
		this.commonReference = commonReference;
		this.centromereBoundaries = centromereBoundaries;
		this.breakupCentromeres = breakupCentromeres;
		this.commonIn = commonIn;
		this.indHash = indHash;
		this.build = build;
		this.CN = CN;
		this.log = log;
	}

	public CNVFilter(Logger log) {
		this(NO_FILTER_MIN_NUM_MARKERS, NO_FILTER_MAX_NUM_MARKERS, NO_FILTER_MIN_SIZE, NO_FILTER_MAX_SIZE, NO_FILTER_MIN_SCORE, NO_FILTER_MAX_SCORE, NO_FILTER_PROBLEM_REGIONS, NO_FILTER_CENTROMERE_MIDPOINTS, NO_FILTER_COMMON_REFERENCE, NO_FILTER_CENTROMERE_BOUNDARIES, NO_FILTER_BREAK_UP_CENTROMERES, DEFAULT_COMMON_IN, NO_FILTER_INDIVIDUAL_HASH, DEFAULT_BUILD, NO_FILTER_CN, log);
	}

	public void setCNVDefaults(Project proj) {
		setMinNumMarkers(DEFAULT_MIN_NUM_MARKERS);
		setMinScore(DEFAULT_MIN_SCORE_THRESHOLD);
		setCentromereBoundariesFromProject(proj);
		setBreakupCentromeres(DEFAULT_BREAK_UP_CENTROMERES);
		setIndividualsToKeepFromSampleData(proj);
		setCommonIn(DEFAULT_COMMON_IN);
	}

	/**
	 * @param cnv
	 *            the cnv to be tested against all filters currently set The returned object is a {@link CNVFilterPass}
	 *            <p>
	 *            Can call {@link CNVFilterPass#passedFilter()} to see if the cnv made it through all filters
	 *            <p>
	 *            Can call {@link CNVFilterPass#getReasonNotPassing()} to get a full report of why the cnv was filtered out. This is useful for testing new filters
	 * @return {@link CNVFilterPass}
	 */
	public CNVFilterPass getCNVFilterPass(CNVariant cnv) {
		CNVFilterPass filterPass = new CNVFilterPass();
		if (minSize != NO_FILTER_MIN_SIZE && cnv.getSize() < minSize) {
			filterPass.setFailed("minimum size < " + minSize, ";");
		}
		if (maxSize != NO_FILTER_MAX_SIZE && cnv.getSize() > maxSize) {
			filterPass.setFailed("maximum size > " + maxSize, ";");
		}
		if (minNumMarkers != NO_FILTER_MIN_NUM_MARKERS && cnv.getNumMarkers() < minNumMarkers) {
			filterPass.setFailed("min number of markers < " + minNumMarkers, ";");
		}
		if (maxNumMarkers != NO_FILTER_MAX_NUM_MARKERS && cnv.getNumMarkers() > maxNumMarkers) {
			filterPass.setFailed("max number of markers > " + maxNumMarkers, ";");
		}
		if (minScore != NO_FILTER_MIN_SCORE && cnv.getScore() < minScore) {
			filterPass.setFailed("minimum score < " + minScore, ";");
		}
		if (maxScore != NO_FILTER_MAX_SCORE && cnv.getScore() > maxScore) {
			filterPass.setFailed("maximum score > " + maxScore, ";");
		}
		if (commonReference != NO_FILTER_COMMON_REFERENCE) {
			if (commonIn && !inOneOfTheseRegions(cnv, commonReference)) {
				filterPass.setFailed("commonIn=" + commonIn + " and was not in commonReference", ";");
			} else if (inOneOfTheseRegions(cnv, commonReference)) {
				filterPass.setFailed("commonIn=" + commonIn + " and was in commonReference", ";");
			}
		}
		if (problemRegions != NO_FILTER_PROBLEM_REGIONS) {
			if (inOneOfTheseRegions(cnv, problemRegions)) {
				filterPass.setFailed("problematic regions were defined and was in a problematic region ", ";");
			}
		}
		if (indHash != NO_FILTER_INDIVIDUAL_HASH && !indHash.contains(cnv.getFamilyID() + "\t" + cnv.getIndividualID())) {
			filterPass.setFailed(cnv.getFamilyID() + "\t" + cnv.getIndividualID() + " was not in the individual filter list", ";");
			filterPass.setIndIsExcluded(true);// this is useful if you are computing concordance, and do not want excluded individuals counted against
		}
		// TODO better, and break up everytime if needed without the extra method
		if (centromereBoundaries != NO_FILTER_CENTROMERE_BOUNDARIES && centromereMidpoints != NO_FILTER_CENTROMERE_MIDPOINTS && cnv.overlaps(centromereMidpoints[cnv.getChr()])) {
			// filterPass.prepFail();
			// filterPass.addReasonFailing("breakupCentromeres=" + breakupCentromeres + " and was in a centromere", ";");
			if (!breakupCentromeres) {
				filterPass.setFailed("breakupCentromeres=" + breakupCentromeres + " and was in a centromere", ";");
			}
			filterPass.setCentromeric(true);
		}
		if (cnv.getSize() > GIANT_BP || cnv.getNumMarkers() > GIANT_NUM_MARKERS) {
			filterPass.setGiant(true);
		}
		if (CN != NO_FILTER_CN && cnv.getCN() != CN) {
			filterPass.setFailed("cnv was not copy number" + CN, ";");
		}
		return filterPass;
	}

	/**
	 * @param proj
	 * @param args
	 *            these are usually the args passed to a main method, if an arg is not valid we do not warn or exit, simply skip
	 * @param filter
	 *            a pre-existing filter, can be null
	 * @param defaults
	 *            set the filter to the defualt values as defined by {@link CNVFilter#setCNVDefaults(Project)}
	 * @param log
	 * @return
	 */
	public static CNVFilter setupCNVFilterFromArgs(Project proj, String[] args, CNVFilter filter, boolean defaults, Logger log) {
		if (filter == null) {
			filter = new CNVFilter(log);

		}
		if (defaults) {
			filter.setCNVDefaults(proj);
		}
		filter.setCommandLineFiltersInEffect(new Hashtable<String, String>());
		filter.setCentromereBoundariesFromFile(proj.getFilename(Project.MARKERSET_FILENAME));// resets if neccesary
		for (int i = 0; i < args.length; i++) {
			if (args[i].startsWith(COMMAND_MIN_NUM_MARKERS)) {
				filter.setMinNumMarkers(ext.parseIntArg(args[i]));
				filter.addCommandLineFilter(args[i], COMMAND_MIN_NUM_MARKERS);
			} else if (args[i].startsWith(COMMAND_MAX_NUM_MARKERS)) {
				filter.setMaxNumMarkers(ext.parseIntArg(args[i]));
				filter.addCommandLineFilter(args[i], COMMAND_MAX_NUM_MARKERS);
			} else if (args[i].startsWith(COMMAND_MIN_SIZE)) {
				filter.setMinSize(ext.parseIntArg(args[i]));
				filter.addCommandLineFilter(args[i], COMMAND_MIN_SIZE);
			} else if (args[i].startsWith(COMMAND_MAX_SIZE)) {
				filter.setMaxSize(ext.parseIntArg(args[i]));
				filter.addCommandLineFilter(args[i], COMMAND_MAX_SIZE);
			} else if (args[i].startsWith(COMMAND_MIN_SCORE)) {
				filter.setMinScore(ext.parseDoubleArg(args[i]));
				filter.addCommandLineFilter(args[i], COMMAND_MIN_SCORE);
			} else if (args[i].startsWith(COMMAND_MAX_SCORE)) {
				filter.setMaxScore(ext.parseDoubleArg(args[i]));
				filter.addCommandLineFilter(args[i], COMMAND_MAX_SCORE);
			} else if (args[i].startsWith(COMMAND_PROBLEM_REGIONS_FILE)) {
				filter.setProblemRegionsFromFile(proj.getProjectDir() + ext.parseStringArg(args[i], ""));
				filter.addCommandLineFilter(args[i], COMMAND_PROBLEM_REGIONS_FILE);
			} else if (args[i].startsWith(COMMAND_COMMON_REFERENCE)) {
				filter.setCommonReferenceFromFile(proj.getProjectDir() + ext.parseStringArg(args[i], ""), DEFAULT_COMMON_IN);
				filter.addCommandLineFilter(args[i], COMMAND_COMMON_REFERENCE);
			} else if (args[i].startsWith(COMMAND_INDIVIDUALS_TO_KEEP)) {
				filter.setIndividualsToKeepFromFile(proj.getProjectDir() + ext.parseStringArg(args[i], ""));
				filter.addCommandLineFilter(args[i], COMMAND_INDIVIDUALS_TO_KEEP);
			} else if (args[i].startsWith(COMMAND_BUILD)) {
				filter.setBuild(ext.parseIntArg(args[i]));
				filter.setCentromereBoundariesFromFile(proj.getFilename(Project.MARKERSET_FILENAME));
				filter.addCommandLineFilter(args[i], COMMAND_BUILD);
			} else if (args[i].startsWith(COMMAND_COMMON_IN)) {
				filter.setCommonIn(ext.parseBooleanArg(args[i]));
				filter.addCommandLineFilter(args[i], COMMAND_COMMON_IN);
			} else if (args[i].startsWith(COMMAND_EXCLUDE_INDIVIDUALS_FROM_SAMPLE_DATA)) {
				if (filter.getIndHash() == NO_FILTER_INDIVIDUAL_HASH && ext.parseBooleanArg(args[i])) {
					filter.setIndividualsToKeepFromSampleData(proj);
				}
				filter.addCommandLineFilter(args[i], COMMAND_EXCLUDE_INDIVIDUALS_FROM_SAMPLE_DATA);
			} else if (args[i].startsWith(COMMAND_BREAK_UP_CENTROMERES)) {
				filter.setBreakupCentromeres(ext.parseBooleanArg(args[i]));
				filter.addCommandLineFilter(args[i], COMMAND_BREAK_UP_CENTROMERES);
			}
		}
		return filter;

	}

	public int getBuild() {
		return build;
	}

	public String getHGBuild() {
		return build == 36 ? "hg18" : "hg19";
	}

	public void setBuild(int build) {
		this.build = build;
	}

	public Hashtable<String, String> getCommandLineFiltersInEffect() {
		return commandLineFiltersInEffect;
	}

	public void addCommandLineFilter(String filter, String value) {
		if (commandLineFiltersInEffect == null) {
			this.commandLineFiltersInEffect = new Hashtable<String, String>();
		}
		commandLineFiltersInEffect.put(filter, value);
	}

	public boolean isCommandLineFilterInEffect(String filter) {
		if (commandLineFiltersInEffect == null) {
			return false;
		} else {
			return commandLineFiltersInEffect.containsKey(filter);
		}
	}

	public void setCommandLineFiltersInEffect(Hashtable<String, String> commandLineFiltersInEffect) {
		this.commandLineFiltersInEffect = commandLineFiltersInEffect;
	}

	/**
	 * @return a String array of the default parameters, typically used for creating a .crf
	 */
	public static String[] getDefaultCNVParams() {
		String[] params = new String[27];
		params[0] = "# CNV Specific Filters: ";

		params[1] = "# minimum number of markers";
		params[2] = COMMAND_MIN_NUM_MARKERS + DEFAULT_MIN_NUM_MARKERS;

		params[3] = "# maximum number of markers";
		params[4] = "#" + COMMAND_MAX_NUM_MARKERS;

		params[5] = "# minimum size in bp";
		params[6] = "#" + COMMAND_MIN_SIZE;

		params[7] = "# maximum size in bp";
		params[8] = "#" + COMMAND_MAX_SIZE;

		params[9] = "# minimum score (recomended 15 for PennCNV and 2.5 for cnv.Beast";
		params[10] = "#" + COMMAND_MIN_SCORE + DEFAULT_MIN_SCORE_THRESHOLD;

		params[11] = "# maximum score";
		params[12] = "#" + COMMAND_MAX_SCORE;

		params[13] = "# a path (relative to the project directory) to a file of problematic regions.";
		params[14] = "#" + COMMAND_PROBLEM_REGIONS_FILE;

		params[15] = "# a path (relative to the project directory) to a file common regions";
		params[16] = "#" + COMMAND_COMMON_REFERENCE;

		params[17] = "# if a common reference is provided, keep only variants in the common regions. Defaults to removing (" + DEFAULT_COMMON_IN + ")";
		params[18] = COMMAND_COMMON_IN + DEFAULT_COMMON_IN;

		params[19] = "# a path (relative to the project directory) to a file common individuals to use (note this will override " + COMMAND_EXCLUDE_INDIVIDUALS_FROM_SAMPLE_DATA + ")";
		params[20] = "#" + COMMAND_INDIVIDUALS_TO_KEEP;

		params[21] = "# the genome build to use for centromere locations";
		params[22] = COMMAND_BUILD + DEFAULT_BUILD;

		params[23] = "# break up CNVs spanning centromers, defualts to removing cnvs that span centromeres(" + DEFAULT_BREAK_UP_CENTROMERES + ")";
		params[24] = "#" + COMMAND_BREAK_UP_CENTROMERES;

		params[25] = "# exclude indivudals as defined by sample data";
		params[26] = COMMAND_EXCLUDE_INDIVIDUALS_FROM_SAMPLE_DATA + DEFAULT_EXCLUDE_SAMPLE_DATA;

		return params;
	}

	public boolean passesCNVFilter(CNVariant cnv) {
		return getCNVFilterPass(cnv).passedFilter();
	}

	/**
	 * A bare bones filter
	 */
	public void setBasicCNVFilters(int minNumMarkers, int maxNumMarkers, int minSize, int maxSize, double minScore, double maxScore) {
		this.minNumMarkers = minNumMarkers;
		this.maxNumMarkers = maxNumMarkers;
		this.minSize = minSize;
		this.maxSize = maxSize;
		this.minScore = minScore;
		this.maxScore = maxScore;
	}

	public int getCN() {
		return CN;
	}

	public void setCN(int cN) {
		CN = cN;
	}

	/**
	 * The method {@link CNVFilter#getCNVFilterPass(CNVariant) } does not explicitly break up centromeres, this must be called to do that
	 */
	public CNVariant[] breakUpCentromere(CNVFilterPass filterPass, CNVariant cnv) {
		CNVariant[] cnvCentromere;
		if (breakupCentromeres && filterPass.passedFilter() && centromereMidpoints != NO_FILTER_CENTROMERE_MIDPOINTS && centromereBoundaries != NO_FILTER_CENTROMERE_BOUNDARIES && filterPass.isCentromeric()) {

			int[] bounds = centromereBoundaries[cnv.getChr()];
			boolean startWithin = cnv.getStart() >= bounds[0];
			boolean endWithin = cnv.getStop() <= bounds[1];

			CNVariant newCNV1 = null;
			CNVariant newCNV2 = null;
			if (startWithin && endWithin) {
				cnvCentromere = new CNVariant[] {};
			} else {
				if (startWithin || !endWithin) {
					// if startWithin or (implied: !startWithin and) !endWithin
					int newStart = bounds[1] + 1;
					if (cnv.getStop() - newStart <= 1) {
						cnvCentromere = new CNVariant[] {};
					}

					int secondMarker = Array.binarySearch(positions[cnv.getChr()], cnv.getStop(), true);
					int firstMarker = Array.binarySearch(positions[cnv.getChr()], bounds[1], true);
					int markerCnt = secondMarker - firstMarker + 1;
					newCNV1 = new CNVariant(cnv.getFamilyID(), cnv.getIndividualID(), cnv.getChr(), positions[cnv.getChr()][firstMarker], cnv.getStop(), cnv.getCN(), cnv.getScore(), markerCnt, cnv.getSource());
				}
				if (endWithin || !startWithin) {
					// if endWithin or (implied: !endWithin and) !startWithin
					int newEnd = bounds[0] - 1;
					if (newEnd - cnv.getStart() <= 1) {
						cnvCentromere = new CNVariant[] {};
					}

					int firstMarker = Array.binarySearch(positions[cnv.getChr()], cnv.getStart(), true);
					int secondMarker = Array.binarySearch(positions[cnv.getChr()], bounds[0], true);
					int markerCnt = secondMarker - firstMarker + 1;
					newCNV2 = new CNVariant(cnv.getFamilyID(), cnv.getIndividualID(), cnv.getChr(), cnv.getStart(), positions[cnv.getChr()][secondMarker], cnv.getCN(), cnv.getScore(), markerCnt, cnv.getSource());
				}

				if (!endWithin && !startWithin) {
					cnvCentromere = new CNVariant[] { newCNV1, newCNV2 };
				} else if (startWithin && !endWithin) {
					cnvCentromere = new CNVariant[] { newCNV1 };
				} else if (endWithin && !startWithin) {
					cnvCentromere = new CNVariant[] { newCNV2 };
				} else {
					cnvCentromere = new CNVariant[] {};
				}
			}
		} else {
			cnvCentromere = new CNVariant[] { cnv };
		}
		return cnvCentromere;
	}

	/**
	 * Any or all can be null,
	 * 
	 * @param commonIn
	 *            include (true) or exclude (false) common variants
	 */
	public void setAllAuxillaryRegionsFromFiles(String fullPathToProblematicRegions, String fullPathToMarkerSetFilename, String fullPathToCommonCNPReference, String fullPathToIndividualsToKeepFile, boolean commonIn) {
		setProblemRegionsFromFile(fullPathToProblematicRegions);
		setCentromereBoundariesFromFile(fullPathToMarkerSetFilename);
		setCommonReferenceFromFile(fullPathToCommonCNPReference, commonIn);
		setIndividualsToKeepFromFile(fullPathToIndividualsToKeepFile);
	}

	public void setProblemRegionsFromFile(String fullPathToProblematicRegions) {
		if (fullPathToProblematicRegions == null || fullPathToProblematicRegions.equals("")) {
			this.problemRegions = NO_FILTER_PROBLEM_REGIONS;
		} else {
			this.problemRegions = Segment.loadUCSCregions(fullPathToProblematicRegions, 0, false, log);
		}
	}

	public void setCentromereBoundariesFromProject(Project proj) {
		setCentromereBoundariesFromFile(proj.getFilename(Project.MARKERSET_FILENAME));
	}

	public void setCentromereBoundariesFromFile(String fullPathToMarkerSetFilename) {
		if (fullPathToMarkerSetFilename == null || fullPathToMarkerSetFilename.equals("")) {
			this.centromereBoundaries = NO_FILTER_CENTROMERE_BOUNDARIES;
		} else {
			MarkerSet markerSet = MarkerSet.load(fullPathToMarkerSetFilename, false);
			this.positions = markerSet.getPositionsByChr();
			this.centromereBoundaries = Positions.determineCentromereBoundariesFromMarkerSet(markerSet.getChrs(), markerSet.getPositions(), build, log);
		}
		computeCentromereMidPoints();
	}

	public void computeCentromereMidPoints() {
		this.centromereMidpoints = centromereBoundaries == NO_FILTER_CENTROMERE_BOUNDARIES ? NO_FILTER_CENTROMERE_MIDPOINTS : Positions.computeCentromereMidpoints(centromereBoundaries);
	}

	public void setCommonReferenceFromFile(String fullPathToCommonCNPReference, boolean commonIn) {
		if (fullPathToCommonCNPReference == null || fullPathToCommonCNPReference.equals("")) {
			this.commonReference = NO_FILTER_COMMON_REFERENCE;
		} else {
			this.commonReference = Segment.loadUCSCregions(fullPathToCommonCNPReference, false);
		}
		this.commonIn = commonIn;
	}

	public void setIndividualsToKeepFromFile(String fullPathToIndividualsToKeepFile) {
		if (fullPathToIndividualsToKeepFile == null || fullPathToIndividualsToKeepFile.equals("")) {
			setIndHash(NO_FILTER_INDIVIDUAL_HASH);
		} else {
			setIndHash(HashVec.loadFileToStringArray(fullPathToIndividualsToKeepFile, false, false, new int[] { 0, 1 }, true, false, "\t"));
		}
	}

	public void setIndividualsToKeepFromSampleData(Project proj) {
		String samples[] = Array.subArray(proj.getSamples(), proj.getSamplesToInclude(null));
		SampleData sampleData = proj.getSampleData(0, false);
		for (int i = 0; i < samples.length; i++) {
			samples[i] = sampleData.lookup(samples[i])[1];
		}
		log.report("Info - retaining " + samples.length + (samples.length == 1 ? " sample" : " samples") + " in the cnv filter set");
		setIndHash(samples);
	}

	public boolean isCommonIn() {
		return commonIn;
	}

	public void setCommonIn(boolean commonIn) {
		this.commonIn = commonIn;
	}

	public int getMinNumMarkers() {
		return minNumMarkers;
	}

	public HashSet<String> getIndHash() {
		return indHash;
	}

	public void setIndHash(String[] individualsToKeep) {
		this.indHash = HashVec.loadToHashSet(individualsToKeep);
	}

	public void setIndHash(HashSet<String> indHash) {
		this.indHash = indHash;
	}

	public void setMinNumMarkers(int minNumMarkers) {
		this.minNumMarkers = minNumMarkers;
	}

	public int getMaxNumMarkers() {
		return maxNumMarkers;
	}

	public void setMaxNumMarkers(int maxNumMarkers) {
		this.maxNumMarkers = maxNumMarkers;
	}

	public double getMinScore() {
		return minScore;
	}

	public void setMinScore(double minScore) {
		this.minScore = minScore;
	}

	public double getMaxScore() {
		return maxScore;
	}

	public void setMaxScore(double maxScore) {
		this.maxScore = maxScore;
	}

	public int getMinSize() {
		return minSize;
	}

	public void setMinSize(int minSize) {
		this.minSize = minSize;
	}

	public int getMaxSize() {
		return maxSize;
	}

	public void setMaxSize(int maxSize) {
		this.maxSize = maxSize;
	}

	public Segment[] getProblemRegions() {
		return problemRegions;
	}

	public void setProblemRegions(Segment[] problemRegions) {
		this.problemRegions = problemRegions;
	}

	public Segment[] getCentromereMidpoints() {
		return centromereMidpoints;
	}

	public void setCentromereMidpoints(Segment[] centromereMidpoints) {
		this.centromereMidpoints = centromereMidpoints;
	}

	public Segment[] getCommonReference() {
		return commonReference;
	}

	public void setCommonReference(Segment[] commonReference) {
		this.commonReference = commonReference;
	}

	public int[][] getCentromereBoundaries() {
		return centromereBoundaries;
	}

	public void setCentromereBoundaries(int[][] centromereBoundaries) {
		this.centromereBoundaries = centromereBoundaries;
	}

	public boolean isBreakupCentromeres() {
		return breakupCentromeres;
	}

	public void setBreakupCentromeres(boolean breakupCentromeres) {
		this.breakupCentromeres = breakupCentromeres;
	}

	public Logger getLog() {
		return log;
	}

	public void setLog(Logger log) {
		this.log = log;
	}

	public static boolean inOneOfTheseRegions(CNVariant cnv, Segment[] regions) {
		for (int i = 0; i < regions.length; i++) {
			if (cnv.significantOverlap(regions[i])) {
				return true;
			}
		}
		return false;
	}

	/**
	 * This Class is designed to store filter metrics associated with a single cnv
	 *
	 */
	public static class CNVFilterPass {
		private boolean passes, centromeric, giant;
		private String reasonNotPassing;
		private boolean indIsExcluded;
		private static final String IT_PASSED = "Passed QC";
		private static final String DID_NOT_PASS = "Did not pass QC for the following reasons:";

		public CNVFilterPass() {
			super();
			this.passes = true;
			this.centromeric = false;
			this.giant = false;
			this.reasonNotPassing = IT_PASSED;
			this.indIsExcluded = false;
		}

		public boolean isIndIsExcluded() {
			return indIsExcluded;
		}

		public void setIndIsExcluded(boolean indIsExcluded) {
			this.indIsExcluded = indIsExcluded;
		}

		public boolean isGiant() {
			return giant;
		}

		public void setGiant(boolean giant) {
			this.giant = giant;
		}

		public boolean passedFilter() {
			return passes;
		}

		public boolean isGiantCentromeric() {
			return giant && centromeric;
		}

		public String getReasonNotPassing() {
			return reasonNotPassing;
		}

		public boolean isCentromeric() {
			return centromeric;
		}

		public void setCentromeric(boolean centromeric) {
			this.centromeric = centromeric;
		}

		public void setFailed(String reason, String sep) {
			if (passes)
				passes = false;
			if (DID_NOT_PASS.equals(reasonNotPassing)) {
				reasonNotPassing += reason;
			} else if (IT_PASSED.equals(reasonNotPassing)) {
				reasonNotPassing = DID_NOT_PASS + reason;
			} else {
				reasonNotPassing += sep + reason;
			}
		}

	}

	public static String[] getParserParams() {
		String[] params = new String[7];
		params[0] = "#To intialize the cnv filter, provide the following arguments";
		params[1] = "#the full path to a project properties file";
		params[2] = COMMAND_PROJECT;
		params[3] = "# the path (relative to the project directory) for a cnv file";
		params[4] = COMMAND_CNV_FILE;
		params[5] = "# the path (relative to the project directory) for the filtered output";
		params[6] = COMMAND_CNV_FILE_OUT;
		params = Array.concatAll(params, getDefaultCNVParams());
		return params;
	}

	public static void fromParametersTrio(String filename, Logger log) {
		Vector<String> params;
		params = Files.parseControlFile(filename, CNVTrioFilter.COMMAND_CNV_TRIO_CRF, getParserParams(), log);
		if (params != null) {
			main(Array.toStringArray(params));
		}
	}

	public static void fromParameters(String filename, Logger log) {
		Vector<String> params;
		params = Files.parseControlFile(filename, CNVTrioFilter.COMMAND_CNV_FILTER_CRF, getParserParams(), log);
		params.add(COMMAND_CNV_FILTER);
		if (params != null) {
			main(Array.toStringArray(params));
		}
	}

	public static void filterCNVFile(Project proj, String cnvFile, String out, CNVFilter cnvFilter) {
		CNVariant[] cnvs = CNVariant.loadPlinkFile(proj.getProjectDir() + cnvFile, false);
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(proj.getProjectDir() + out));
			writer.println(Array.toStr(CNVariant.PLINK_CNV_HEADER));
			for (int i = 0; i < cnvs.length; i++) {
				CNVFilterPass pass = cnvFilter.getCNVFilterPass(cnvs[i]);
				if (pass.passedFilter()) {
					writer.println(cnvs[i].toPlinkFormat());
				} else {
					// proj.getLog().report(pass.getReasonNotPassing()+"\t"+cnvs[i].toPlinkFormat());
				}
			}
			writer.close();
		} catch (Exception e) {
			proj.getLog().reportError("Error writing to " + proj.getProjectDir() + out);
			proj.getLog().reportException(e);
		}

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;
		String cnvFile = "Genvisis.cnv";
		String out = "Genvisis.filt.cnv";
		boolean filter = false;
		String logfile = null;
		Project proj;

		String usage = "\n" + "cnv.qc.CNVFilter requires 1 argument\n";
		usage += "   (1) project file name (i.e. " + COMMAND_PROJECT + filename + " (no default))\n" + "";
		usage += "   (2) cnv file name (relative to the project directory) (i.e. " + COMMAND_CNV_FILE + cnvFile + " ( default))\n" + "";
		usage += "   (3) output (relative to the project directory) (i.e. " + COMMAND_CNV_FILE_OUT + out + " ( default))\n" + "";

		if (ext.indexOfStr(COMMAND_PROJECT, args, true, false) >= 0) {
			proj = new Project(ext.parseStringArg(args[ext.indexOfStr(COMMAND_PROJECT, args, true, false)], ""), logfile, false);
		} else {
			proj = new Project(filename, null, false);
		}
		CNVFilter cnvFilter = setupCNVFilterFromArgs(proj, args, null, true, proj.getLog());
		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				logfile = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith(COMMAND_CNV_FILE_OUT)) {
				out = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith(COMMAND_CNV_FILE)) {
				cnvFile = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith(COMMAND_PROJECT)) {
				numArgs--;
			} else if (args[i].startsWith(COMMAND_CNV_FILTER)) {
				filter = true;
				numArgs--;
			} else if (cnvFilter.isCommandLineFilterInEffect(args[i])) {
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
			if (filter) {
				filterCNVFile(proj, cnvFile, out, cnvFilter);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
