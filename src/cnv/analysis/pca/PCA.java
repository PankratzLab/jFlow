package cnv.analysis.pca;

import java.util.ArrayList;

import cnv.filesys.Project;
import cnv.filesys.SampleList;
import cnv.qc.GcAdjustorParameter;
import cnv.qc.GcAdjustorParameter.GcAdjustorParameters;
import common.Files;
import common.Logger;
import common.ext;

/**
 * <p>
 * One class to bring other PC classes together
 * 
 *
 */
public class PCA {
	public static final String[] FILE_EXTs = { ".PCs.extrapolated.txt", ".PCs.summary.txt" };
	private static final String EVALUATION_FILENAME = "Evaluated_PCs.txt";
	public static final String PCA_SAMPLES = ".samples.USED_PC.txt";

	/**
	 * @param proj
	 *            project to use
	 * @param excludeSamples
	 *            exlude samples as defined in sample data
	 * @param numComponents
	 *            number of principal components to compute
	 * @param printFullData
	 *            print all the data used to compute principal components (used only for testing)
	 * @param center
	 *            center each marker about it's mean value across samples
	 * @param reportMarkerLoadings
	 *            report the loadings for each marker/component
	 * @param reportSingularValues
	 *            report the singular values
	 * @param imputeMeanForNaN
	 *            for samples with NaN values, impute with the mean of the marker
	 * @param recomputeLRR
	 *            recompute the log R ratios on the fly
	 * @param useFile
	 *            a subset of individuals to use
	 * @param output
	 *            a base name
	 * @param log
	 *            Warning - if the principal component, singular values, and marker loadings files already exist, the {@link PrincipalComponentsCompute} object returned will only have the existing filenames populated. The U,V,and W matrices will not be computed again
	 */
	public static PrincipalComponentsCompute computePrincipalComponents(Project proj, boolean excludeSamples, int numComponents, boolean printFullData, boolean center, boolean reportMarkerLoadings, boolean reportSingularValues, boolean imputeMeanForNaN, boolean recomputeLRR, String useFile, String output, GcAdjustorParameters parameters) {
		return PrincipalComponentsCompute.getPrincipalComponents(proj, excludeSamples, numComponents, printFullData, center, reportMarkerLoadings, reportSingularValues, imputeMeanForNaN, recomputeLRR, useFile, output, parameters);
	}

	/**
	 * @param proj
	 *            as above
	 * @param numComponents
	 *            number of components to apply, must be less than or equal to the number computed
	 * @param singularFile
	 * @param markerLoadingFile
	 * @param useFile
	 * @param excludeSamples
	 * @param imputeMeanForNaN
	 * @param recomputeLRR
	 *            recompute the log R ratios on the fly
	 * @param output
	 * @param log
	 *            Warning - if the extrapolated principal component file already exists, the {@link PrincipalComponentsApply} object returned will only have the extrapolated pc file populated
	 */
	public static PrincipalComponentsApply applyLoadings(Project proj, int numComponents, String singularFile, String markerLoadingFile, String useFile, boolean excludeSamples, boolean imputeMeanForNaN, boolean recomputeLRR, String output, GcAdjustorParameters params) {
		// first retrieve the samples to apply marker loadings to
		boolean[] samplesToUse = PrincipalComponentsCompute.getSamples(proj, excludeSamples, useFile);
		PrincipalComponentsApply pcApply = new PrincipalComponentsApply(proj, numComponents, singularFile, markerLoadingFile, samplesToUse, imputeMeanForNaN, recomputeLRR);
		if (pcApply.outputExists(proj.PROJECT_DIRECTORY.getValue() + ext.rootOf(output) + FILE_EXTs[0], true)) {// warn about existence when checking
			pcApply.setExtrapolatedPCsFile(ext.rootOf(output) + FILE_EXTs[0]);
		} else {
			pcApply.setParams(params);
			pcApply.applyLoadings();
			pcApply.reportExtropolatedPCs(ext.rootOf(output) + FILE_EXTs[0]);
		}
		return pcApply;
	}

	/**
	 * @param proj
	 *            as above
	 * @param pcFile
	 *            either a raw or extrapolated principal component file
	 * @param markersToAssessFile
	 *            the markers to use for median values/resdiual reporting
	 * @param numComponents
	 *            number of components to regress out
	 * @param printFull
	 *            print all data associated with the median markers
	 * @param gcThreshold
	 *            threshold to use for median value computation
	 * @param homozygousOnly
	 *            only compute using homozygous calls
	 * @param recomputeLRR
	 *            recompute the log R ratios on the fly
	 * @param output
	 * @param log
	 * @return
	 */
	public static PrincipalComponentsResiduals computeResiduals(Project proj, String pcFile, String markersToAssessFile, int numComponents, boolean printFull, float gcThreshold, boolean homozygousOnly, boolean recomputeLRR, String output, GcAdjustorParameters params) {
		PrincipalComponentsResiduals pcResids = new PrincipalComponentsResiduals(proj, pcFile, markersToAssessFile, numComponents, printFull, gcThreshold, homozygousOnly, recomputeLRR, ext.rootOf(output));
		String file = proj.PROJECT_DIRECTORY.getValue() + ext.rootOf(ext.rootOf(output) + FILE_EXTs[1]) + PrincipalComponentsResiduals.MT_REPORT_EXT[0];
		if (!Files.exists(file)) {
			pcResids.setParams(params);// params can be null
			pcResids.computeAssessmentDataMedians();
			pcResids.computeResiduals();
			pcResids.computeInverseTransformResiduals();
			pcResids.summarize(ext.rootOf(output) + FILE_EXTs[1]);
		} else {
			pcResids.setResidOutput(ext.removeDirectoryInfo(file));
			proj.getLog().report("Detected that the following report file already exists:\n" + output + "\n");
			proj.getLog().report("Skipping the report generation and using this file instead");
			proj.getLog().report("If this is incorrect (using a different number of components, new samples, etc...),  please remove or change the name of the file listed above.\n Alternatively, specify a new analysis name");
		}
		return pcResids;
	}
	
	public static PrincipalComponentsApply generateFullPCA(Project proj, int numComponents, String outputBase, boolean recomputeLRR_PCs, boolean imputeMeanForNaN,GcAdjustorParameters params, Logger log) {
		log.report("\nReady to perform the principal components analysis (PCA)\n");
		PrincipalComponentsCompute pcs = PCA.computePrincipalComponents(proj, false, numComponents, false, false, true, true, imputeMeanForNaN, recomputeLRR_PCs, proj.PROJECT_DIRECTORY.getValue() + outputBase + PCA_SAMPLES, outputBase, params);

		if (pcs == null) {
			return null;
		}
		// apply PCs to everyone, we set useFile to null and excludeSamples to false to get all samples in the current project.
		// TODO, if we ever want to apply to only a subset of the project, we can do that here.....
		log.report("\nApplying the loadings from the principal components analysis to all samples\n");
		PrincipalComponentsApply pcApply = PCA.applyLoadings(proj, numComponents, pcs.getSingularValuesFile(), pcs.getMarkerLoadingFile(), null, false, imputeMeanForNaN, recomputeLRR_PCs, outputBase, null);
		// Compute Medians for (MT) markers and compute residuals from PCs for everyone
		log.report("\nComputing residuals after regressing out " + numComponents + " principal component" + (numComponents == 1 ? "" : "s") + "\n");
		return pcApply;
	}

	public static void main(String[] args) {
		String filename = null;
		String logfile = "PCA.log";
		String useFile = null;
		String markersToassessFile = "MT_Markers.txt";
		String markerLoadingFile = "loadings.txt";
		String singularValueFile = "singularValues.txt";
		String useApplyfile = null;
		String pcFile = "PCs.txt";
		String evalOut = EVALUATION_FILENAME;
		int numArgs = args.length;
		boolean excludeSamples = false;
		boolean printFullData = false;
		boolean applyLoadings = false;
		boolean computeResiduals = false;
		boolean imputeMeanForNaN = false;
		boolean homozygousOnly = false;
		boolean center = false;
		boolean recomputeLRR = false;
		float gcThreshold = 0f;
		int numComponents = 100;

		String usage = // String usage = "\n"+
		"cnv.analysis.PCA requires 1 argument\n" +
				"   To generate principal components, use the following options \n" +
				"   (1) project (i.e. proj=" + filename + " (default))\n" +
				"  OPTIONAL \n" +
				"   (2) output base name (i.e. out=" + evalOut + " (default))\n" +
				"   (3) exclude samples as defined in Sample Data for PCA (i.e. -exclude (not the default))\n" +
				"   (4) supply a file with a list of samples (one per line) to use for PCA computation (i.e. useFile=" + useFile + " (default , use all samples))\n" +
				"   (5) number of principal components to compute  (i.e. components=" + numComponents + " (default))\n" +
				"   (6) mean center the data (each marker) prior to PC computation (i.e. -center (not the default))\n" +
				"   (7) print full marker data used for PCA (i.e. -printFull (not the default))\n" +
				"   (8) impute the mean marker value for samples with NaN Log R Ratios when computing PCs(i.e. -impute (not the default))\n" +
				"   (9) logfile (i.e. log=" + logfile + " (default))\n" +
				"   (10) recompute Log R Ratios for each marker from genotypes/intensities (i.e. -recomputeLRR (not the default))\n" +

				"  NOTE\n" +
				"      Markers in the target marker file will be used for the PCA \n" +
				"      Imputing mean marker values for NaN intensities is not the default, instead, the marker is dropped from the PCA computation if any sample has a NaN intensity\n" +
				"  OR \n" +
				"   To apply marker loadings to a new set of project samples (and extrapolate new PCs), use the following options \n" +
				"   (1) apply loadings (i.e. -apply (not the default))\n" +
				"   (2) project (i.e. proj=" + filename + " (default))\n" +
				"   (3) output file baseName for extrapolation (i.e. out=" + evalOut + " (default))\n" +
				"   (4) a marker loading file to use (i.e. markerLoadingFile=" + markerLoadingFile + " (default))\n" +
				"   (5) a file containing singular values (i.e. singularValueFile=" + singularValueFile + " (default))\n" +
				"   (6) supply a file with a list of samples (one per line) to use for PCA extrapolation (i.e. useSampleFile=" + useFile + " (default , apply to all samples))\n" +
				"   (7) number of principal components (using marker loadings to apply)   (i.e. components=" + numComponents + " (default))\n" +
				"   (8) exclude samples from pc extrapolation as defined in Sample Data for PCA extrapolation (i.e. -exclude (not the default))\n" +
				"   (9) impute the mean marker value for samples with NaN Log R Ratios when applying marker loadings (i.e. -impute (not the default))\n" +
				"   (10) logfile (i.e. log=" + logfile + " (default))\n" +
				"  NOTE\n" +
				" Imputing mean marker values for NaN intensities is not the default, instead, the marker's loadings are not applied if any sample has a NaN intensity\n" +
				"  OR \n" +
				"   To compute median LRR values for a set of markers using PCs, use the following options \n" +
				"   (1) compute residuals  (i.e. -residuals (not the default))\n" +
				"   (2) project (i.e. proj=" + filename + " (default))\n" +
				"   (4) output file baseName for residuals (i.e. out=" + evalOut + " (default))\n" +
				"   (5) supply a file with a list of markers (one per line) to use to compute median LRR values and residuals from PCs (i.e. markers=" + markersToassessFile + " (default))\n" +
				"   (6) pc file to use for residuals  (i.e. pcFile=" + pcFile + " (default))\n" +
				"   (7) number of principal components to use for residuals (i.e. components=" + numComponents + " (default))\n" +
				"   (8) gcThreshold filter for computing median LRR values (i.e. gcThreshold=" + gcThreshold + " (default, no filtering))\n" +
				"   (9) compute median LRR values off homozygous calls only (i.e. -homozygousOnly (not the default)\n" +
				"   (10) logfile (i.e. log=" + logfile + " (default))\n" +
				"   (11) print full marker data used for median LRR computations (i.e. -printFull (not the default))\n" +
				"   (12) recompute Log R Ratios for each marker from genotypes/intensities (i.e. -recomputeLRR (not the default))\n" +
				"  NOTE\n" +
				"  Markers with NaN values for any sample will not be included in either the PCA computation or PCA extrapolation unless -impute is flagged\n" +
				"  However, markers with NaN values used to compute the median Log R Ratio will only be excluded from that sample\n" +
				"  NOTE\n" +
				"  All files must be located in the project directory\n" +
				"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("pcFile=")) {
				pcFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				logfile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("out=")) {
				evalOut = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("-exclude")) {
				excludeSamples = true;
				numArgs--;
			} else if (args[i].startsWith("-printFull")) {
				printFullData = true;
				numArgs--;
			} else if (args[i].startsWith("-residuals")) {
				computeResiduals = true;
				numArgs--;
			} else if (args[i].startsWith("-apply")) {
				applyLoadings = true;
				numArgs--;
			} else if (args[i].startsWith("components=")) {
				numComponents = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("gcThreshold=")) {
				gcThreshold = Float.parseFloat(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("useSampleFile=")) {
				useFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("markers=")) {
				markersToassessFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("markerLoadingFile=")) {
				markerLoadingFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("singularValueFile=")) {
				singularValueFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("useApplyfile=")) {
				useApplyfile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("-impute")) {
				imputeMeanForNaN = true;
				numArgs--;
			} else if (args[i].startsWith("-center")) {
				center = true;
				numArgs--;
			} else if (args[i].startsWith("-homozygousOnly")) {
				homozygousOnly = true;
				numArgs--;
			} else if (args[i].startsWith("-recomputeLRR")) {
				recomputeLRR = true;
				numArgs--;
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		Project proj = new Project(filename, logfile, false);

		if (applyLoadings) {
			applyLoadings(proj, numComponents, singularValueFile, markerLoadingFile, useFile, excludeSamples, imputeMeanForNaN, recomputeLRR, evalOut, null);
		} else if (computeResiduals) {
			computeResiduals(proj, pcFile, markersToassessFile, numComponents, printFullData, gcThreshold, homozygousOnly, recomputeLRR, evalOut, null);
		} else {
			PrincipalComponentsCompute pc = computePrincipalComponents(proj, excludeSamples, numComponents, printFullData, center, true, true, imputeMeanForNaN, recomputeLRR, useFile, evalOut, null);
			applyLoadings(proj, numComponents, pc.getSingularValuesFile(), pc.getMarkerLoadingFile(), useApplyfile, excludeSamples, imputeMeanForNaN, recomputeLRR, evalOut, null);
		}
	}
}
