package cnv.analysis.pca;

import java.io.File;
import java.io.FileNotFoundException;

import common.ext;
import stats.StatsCrossTabs;
import stats.StatsCrossTabs.CORREL_TYPE;
import cnv.filesys.Project;
import cnv.manage.ExtProjectDataParser;
import cnv.manage.ExtProjectDataParser.Builder;
import cnv.qc.SampleQC;

/**
 * Class that may be used for pruning/testing principal components according to correlation with sample QC and other metrics TODO , implement a pc-selector based of thresholds....
 */
public class PrincipalComponentsCrossTabs extends PrincipalComponentsResiduals {

	private SampleQC sampleQC;
	private int numPCs;
	private StatsCrossTabs sTabs;

	public PrincipalComponentsCrossTabs(PrincipalComponentsResiduals pResiduals, int numPCs) {
		super(pResiduals);
		this.sampleQC = SampleQC.loadSampleQC(getProj());
		this.numPCs = numPCs;
	}

	public SampleQC getSampleQC() {
		return sampleQC;
	}

	public void developCrossTabs(String fullPathToAlternateDataFile) {
		double[][] additionalData = null;
		String[] additionalDataTitles = null;
		if (fullPathToAlternateDataFile != null) {
			ExtProjectDataParser extSampleFileParser = loadExternalData(fullPathToAlternateDataFile);
			additionalData = extSampleFileParser.getNumericData();
			additionalDataTitles = extSampleFileParser.getNumericDataTitles();
		}
		log.reportTimeInfo("Developing cross tabs...");
		this.sTabs = getCorrelationTable(this, numPCs, sampleQC, CORREL_TYPE.PEARSON, additionalData, additionalDataTitles);
		log.reportTimeInfo("Finished developing cross tabs...");

	}

	public void dumpTables() {
		String outputDir = getProj().getDir(Project.RESULTS_DIRECTORY);
		proj.getLog().reportTimeInfo("Dumping data to " + outputDir);
		new File(outputDir).mkdirs();
		sTabs.dumpTables(outputDir + "sampleQC_PC.crosstabs");
	}

	private ExtProjectDataParser loadExternalData(String fullPathToAlternateDataFile) {
		ExtProjectDataParser.Builder builder = new Builder();
		ExtProjectDataParser extSampleFileParser = null;
		builder.treatAllNumeric(true);

		try {
			extSampleFileParser = builder.build(getProj(), fullPathToAlternateDataFile);
		} catch (FileNotFoundException e) {
			getProj().getLog().reportTimeError("Could not read file " + fullPathToAlternateDataFile);
			e.printStackTrace();
		}
		extSampleFileParser.loadData();
		return extSampleFileParser;
	}

	private static StatsCrossTabs getCorrelationTable(PrincipalComponentsCrossTabs pCorrelation, int numPcs, SampleQC sampleQC, StatsCrossTabs.CORREL_TYPE cType, double[][] additionalData, String[] additionalDataTitles) {
		int numCorrels = sampleQC.getQctitles().length + numPcs + (additionalData == null ? 0 : additionalData.length);
		double[][] data = new double[numCorrels][];
		String[] titles = new String[numCorrels];
		int dataIndex = 0;
		for (int i = 0; i < sampleQC.getQctitles().length; i++) {
			titles[dataIndex] = sampleQC.getQctitles()[i];
			data[dataIndex] = sampleQC.getDataFor(titles[dataIndex]);
			dataIndex++;
		}
		for (int i = 0; i < numPcs; i++) {
			titles[dataIndex] = "PC" + (i + 1);
			data[dataIndex] = pCorrelation.getBasisAt((i + 1));
			dataIndex++;
		}
		if (additionalData != null) {
			if (additionalDataTitles == null || additionalDataTitles.length != additionalData.length) {
				pCorrelation.getProj().getLog().reportTimeError("Additional data must have identical lengths");
				return null;
			} else {
				for (int i = 0; i < additionalData.length; i++) {
					titles[dataIndex] = additionalDataTitles[i];
					data[dataIndex] = additionalData[i];
					dataIndex++;
				}
			}
		}
		StatsCrossTabs statsCrossTabs = new StatsCrossTabs(data, titles, cType, pCorrelation.getProj().getLog());
		statsCrossTabs.computeTable();
		return statsCrossTabs;
	}

	public static void crossTabulate(Project proj, int numPCs, String alternateDataFile) {
		PrincipalComponentsResiduals pcResiduals = proj.loadPcResids();
		numPCs = numPCs > 0 ? numPCs : pcResiduals.getNumComponents();
		PrincipalComponentsCrossTabs pcCorrelation = new PrincipalComponentsCrossTabs(pcResiduals, numPCs);
		pcCorrelation.developCrossTabs(alternateDataFile == null ? null : proj.getProjectDir() + alternateDataFile);
		pcCorrelation.dumpTables();

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;
		int numPCs = 0;
		String alternateDataFile = null;

		String usage = "\n" + "cnv.analysis.pca.PrincipalComponentsCrossTabs requires 0-1 arguments\n";
		usage += "   (1) filename (i.e. proj=" + filename + " (no default))\n" + "";
		usage += "   (2) number of principal components to cross tabulate (i.e. numPCs=" + numPCs + " (defaults to all PCs in the file))\n" + "";
		usage += "   (3) an optional data file to add, must have a column with \"DNA\" and numeric data in all other columns.  (i.e. alternateDataFile=" + alternateDataFile + " (no default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("alternateDataFile=")) {
				alternateDataFile = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("numPCs=")) {
				numPCs = ext.parseIntArg(args[i]);
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
			Project proj = new Project(filename, false);
			crossTabulate(proj, numPCs, alternateDataFile);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
