package cnv.analysis.pca;

import common.Files;
import common.Logger;
import cnv.filesys.Project;
import cnv.var.SampleData;

/**
 * @author lane0212 Class to compute percent race with genotyping pcs, hijacking the {@link PrincipalComponentsResiduals} methods
 */
public class PCPopulation {
	public static final String POP_DEFINITION = "CLASS=POP_DEF";
	private Project proj;
	private SampleData sampleData;
	private PrincipalComponentsResiduals pResiduals;
	private boolean verbose;
	private Logger log;

	public PCPopulation(Project proj, PrincipalComponentsResiduals pResiduals, boolean verbose) {
		super();
		this.proj = proj;
		this.sampleData = proj.getSampleData(0, false);
		this.pResiduals = pResiduals;
		this.log = proj.getLog();
	}

	private void extractDefs() {
		pResiduals.getSamplesInPc();
		String[] sampleDefs = sampleData.getClasses();
		int index = -1;
		for (int i = 0; i < sampleDefs.length; i++) {
			if (sampleDefs[i].startsWith(POP_DEFINITION)) {
				if (verbose) {
					log.reportTimeInfo("Detected sample population definition at index " + i + " (" + sampleDefs[i] + ")");
					index = i;
				}
			}
		}
		if (index < 0) {
			log.reportTimeError("Did not find class " + POP_DEFINITION + " in sample data file " + proj.SAMPLE_DATA_FILENAME.getValue());
		} else {
			String[] samplesInPcs = pResiduals.getSamplesToReport();
			int[] sampleClasses = new int[samplesInPcs.length];
			for (int i = 0; i < samplesInPcs.length; i++) {
				sampleClasses[i] = sampleData.getClassForInd(samplesInPcs[i], index);
			}
		}
	}

	public static void test(Project proj, String genoPCfile) {
		PrincipalComponentsResiduals pResiduals = new PrincipalComponentsResiduals(genoPCfile, Files.getHeaderOfFile(genoPCfile, proj.getLog()).length - 2, proj.getLog());
		PCPopulation pcPopulation = new PCPopulation(proj, pResiduals, true);

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;
		String genoPCfile = null;
		String logfile = null;
		Logger log;

		String usage = "\n" + "cnv.analysis.pca.PCPopulation requires 0-1 arguments\n";
		usage += "   (1) project filename (i.e. proj=" + filename + " (default))\n" + "";
		usage += "   (2) full path to a genotype pc file (i.e. file=" + filename + " (default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("pc=")) {
				genoPCfile = args[i].split("=")[1];
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
			log = new Logger(logfile);
			test(new Project(filename, false), genoPCfile);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
