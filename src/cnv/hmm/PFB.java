package cnv.hmm;

import java.io.FileNotFoundException;

import cnv.filesys.Project;
import cnv.manage.ExtProjectDataParser;
import cnv.manage.ExtProjectDataParser.Builder;

/**
 * @author lane0212 Handles the pfb data for the hmm in reduced format
 */
public class PFB {
	private Project proj;
	private double[] pfbs;

	private PFB(Project proj, double[] pfbs) {
		super();
		this.proj = proj;
		this.pfbs = pfbs;
		if (pfbs.length != proj.getMarkerNames().length) {
			String error = "Found " + pfbs.length + " pfb entries, but the project has" + pfbs.length + " markers";
			proj.getLog().reportTimeError(error);
			throw new IllegalArgumentException(error);
		} else {
			this.proj.getLog().reportTimeInfo("Loaded " + pfbs.length + " pfb entries");
		}
	}

	public double[] getPfbs() {
		return pfbs;
	}

	public static PFB loadPFB(Project proj) {
		return loadPFB(proj, proj.CUSTOM_PFB_FILENAME.getValue());
	}

	public static PFB loadPFB(Project proj, String fullPathToPfb) {

		Builder builder = new Builder();
		builder.dataKeyColumnName("Name");
		builder.numericDataTitles(new String[] { "PFB" });
		builder.sampleBased(false);
		builder.requireAll(true);
		builder.treatAllNumeric(false);

		try {
			ExtProjectDataParser extProjectDataParser = builder.build(proj, fullPathToPfb);
			extProjectDataParser.determineIndicesFromTitles();
			extProjectDataParser.loadData();
			double[] pfbs = extProjectDataParser.getNumericDataForTitle("PFB");
			return new PFB(proj, pfbs);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			proj.getLog().reportFileNotFound(fullPathToPfb);
			return null;
		}

	}

}
