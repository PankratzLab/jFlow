package cnv.analysis.pca;

import java.io.FileNotFoundException;
import java.util.ArrayList;

import stats.CrossValidation;
import stats.ICC;
import common.Array;
import common.Files;
import common.Logger;
import cnv.analysis.pca.PrincipalComponentsResiduals.PrincipalComponentsIterator;
import cnv.filesys.Project;
import cnv.manage.ExtProjectDataParser;

public class CorrectionEvaluator {
	private static final String[] STRING_DATA = new String[] {};// CLASS Sex
	private static final String[][] EVAL_MASKS = new String[][] { { "0" }, { "0", "-1" } };
	private static final String[] DOUBLE_DATA = new String[] { "AGE" };
	private static final String[] DOUBLE_DATA_PATTERN = new String[] { "EVAL_DATA" };
	private static final String[] STRING_DATA_PATTERN = new String[] { "EVAL_CLASS","AGE" };

	private Project proj;
	private PrincipalComponentsResiduals pcResiduals;
	private boolean[] samplesToInclude;
	private String[] matchDouble, matchString;
	private ExtProjectDataParser parser;
	private Logger log;

	public CorrectionEvaluator(Project proj, PrincipalComponentsResiduals pcResiduals, boolean[] samplesToExclude) {
		super();
		this.proj = proj;
		this.pcResiduals = pcResiduals;
		this.samplesToInclude = samplesToExclude;
		this.log = proj.getLog();
		this.matchDouble = gatherPatternTitles(proj.SAMPLE_DATA_FILENAME.getValue(), DOUBLE_DATA_PATTERN, log);
		this.matchString = gatherPatternTitles(proj.SAMPLE_DATA_FILENAME.getValue(), STRING_DATA_PATTERN, log);
		loadSampleData();
	}

	private void evaluate(int[] order) {
		PrincipalComponentsIterator iterator = new PrincipalComponentsIterator(pcResiduals, order);
		while (iterator.hasNext()) {
			PrincipalComponentsResiduals tmpResiduals = iterator.next();
			CrossValidation cValidation = tmpResiduals.getCorrectedDataAt(tmpResiduals.getMedians(), null, samplesToInclude, tmpResiduals.getNumComponents(), false, "HFDS", true);

			for (int i = 0; i < DOUBLE_DATA.length; i++) {

			}
			for (int i = 0; i < matchString.length; i++) {
				String[] response = parser.getStringDataForTitle(matchString[i]);
				System.out.println(response.length + "\t" + proj.getSamples().length);
				ICC icc = new ICC(cValidation.getResiduals(), response, EVAL_MASKS[0], null, true, log);
				icc.computeICC();
				System.out.println(matchString[i] + "\t" + icc.getICC());
			}

		}
	}

	private void loadSampleData() {
		log.reportTimeInfo("Found " + matchDouble.length + "(" + Array.toStr(matchDouble) + ") data columns to load matching the patterns defined by " + Array.toStr(DOUBLE_DATA_PATTERN));
		log.reportTimeInfo("Found " + matchString.length + "(" + Array.toStr(matchString) + ") String columns to load matching the patterns defined by " + Array.toStr(STRING_DATA_PATTERN));
		ExtProjectDataParser.Builder builder = new ExtProjectDataParser.Builder();
		builder.sampleBased(true);
		builder.treatAllNumeric(false);
		builder.requireAll(true);
		builder.verbose(true);
		builder.dataKeyColumnName("DNA");
		builder.stringDataTitles(Array.concatAll(STRING_DATA, matchString));
		builder.numericDataTitles(Array.concatAll(DOUBLE_DATA, matchDouble));
		try {
			log.reportTimeInfo("Loading " + proj.SAMPLE_DATA_FILENAME.getValue());
			this.parser = builder.build(proj, proj.SAMPLE_DATA_FILENAME.getValue());
			parser.determineIndicesFromTitles();
			parser.loadData();

			log.reportTimeInfo("Finished loading " + proj.SAMPLE_DATA_FILENAME.getValue());

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}

	private String[] gatherPatternTitles(String dataFile, String[] patterns, Logger log) {

		String[] header = Files.getHeaderOfFile(dataFile, log);
		ArrayList<String> matches = new ArrayList<String>();
		for (int i = 0; i < header.length; i++) {
			for (int j = 0; j < patterns.length; j++) {
				if (header[i].startsWith(patterns[j])) {
					matches.add(header[i]);
				}
			}
		}
		return matches.toArray(new String[matches.size()]);
	}

	public static void run(Project proj, String markesToEvaluate) {
		proj.getLog().reportTimeInfo("Loading " + proj.INTENSITY_PC_FILENAME.getValue());
		PrincipalComponentsResiduals pcResiduals = proj.loadPcResids();
		pcResiduals.setMarkersToAssessFile(markesToEvaluate);
		pcResiduals.computeAssessmentDataMedians();
		proj.getLog().reportTimeInfo("Finished loading " + proj.INTENSITY_PC_FILENAME.getValue());
		boolean[] samplesToInclude = proj.getSamplesToInclude(null);
		CorrectionEvaluator cEvaluator = new CorrectionEvaluator(proj, pcResiduals, samplesToInclude);
		cEvaluator.evaluate(null);
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String proj = null;
		String markers = null;
		String logfile = null;

		String usage = "\n" + "cnv.analysis.pca.CorrectionEvaluator requires 0-1 arguments\n";
		usage += "   (1) project filename (i.e. proj=" + proj + " (default))\n" + "";
		usage += "   (2) markers to Evaluate (i.e. markers=" + markers + " (default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				proj = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("markers=")) {
				markers = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				logfile = args[i].split("=")[1];
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
			run(new Project(proj, false), markers);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
