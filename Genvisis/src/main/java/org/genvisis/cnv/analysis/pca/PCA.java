package org.genvisis.cnv.analysis.pca;

import org.genvisis.cnv.analysis.pca.PrincipalComponentsCompute.PRE_PROCESSING_METHOD;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.qc.GcAdjustorParameter.GcAdjustorParameters;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

/**
 * <p>
 * One class to bring other PC classes together
 */
public class PCA {

  public static final String[] FILE_EXTs = {".PCs.extrapolated.txt", ".PCs.summary.txt"};
  private static final String EVALUATION_FILENAME = "Evaluated_PCs.txt";
  public static final String PCA_SAMPLES = ".samples.USED_PC.txt";

  /**
   * @param proj project to use
   * @param excludeSamples exlude samples as defined in sample data
   * @param numComponents number of principal components to compute
   * @param printFullData print all the data used to compute principal components (used only for
   *          testing)
   * @param center center each marker about it's mean value across samples
   * @param reportMarkerLoadings report the loadings for each marker/component
   * @param reportSingularValues report the singular values
   * @param imputeMeanForNaN for samples with NaN values, impute with the mean of the marker
   * @param recomputeLRR recompute the log R ratios on the fly
   * @param useFile a subset of individuals to use
   * @param output a base name
   * @param log Warning - if the principal component, singular values, and marker loadings files
   *          already exist, the {@link PrincipalComponentsCompute} object returned will only have
   *          the existing filenames populated. The U,V,and W matrices will not be computed again
   */
  public static PrincipalComponentsCompute computePrincipalComponents(Project proj,
                                                                      boolean excludeSamples,
                                                                      int numComponents,
                                                                      boolean printFullData,
                                                                      PRE_PROCESSING_METHOD method,
                                                                      boolean reportMarkerLoadings,
                                                                      boolean reportSingularValues,
                                                                      boolean imputeMeanForNaN,
                                                                      boolean recomputeLRR,
                                                                      String useFile, String output,
                                                                      GcAdjustorParameters parameters) {
    return PrincipalComponentsCompute.getPrincipalComponents(proj, excludeSamples, numComponents,
                                                             printFullData, method,
                                                             reportMarkerLoadings,
                                                             reportSingularValues, imputeMeanForNaN,
                                                             recomputeLRR, useFile, output,
                                                             parameters);
  }

  /**
   * @param proj as above
   * @param numComponents number of components to apply, must be less than or equal to the number
   *          computed
   * @param singularFile
   * @param markerLoadingFile
   * @param useFile
   * @param excludeSamples
   * @param imputeMeanForNaN
   * @param recomputeLRR recompute the log R ratios on the fly
   * @param output
   * @param log Warning - if the extrapolated principal component file already exists, the
   *          {@link PrincipalComponentsApply} object returned will only have the extrapolated pc
   *          file populated
   */
  public static PrincipalComponentsApply applyLoadings(Project proj, int numComponents,
                                                       String singularFile,
                                                       String markerLoadingFile, String useFile,
                                                       boolean excludeSamples,
                                                       boolean imputeMeanForNaN,
                                                       boolean recomputeLRR, String output,
                                                       GcAdjustorParameters params) {
    // first retrieve the samples to apply marker loadings to
    boolean[] samplesToUse = PrincipalComponentsCompute.getSamples(proj, excludeSamples, useFile);
    PrincipalComponentsApply pcApply = new PrincipalComponentsApply(proj, numComponents,
                                                                    singularFile, markerLoadingFile,
                                                                    samplesToUse, imputeMeanForNaN,
                                                                    recomputeLRR);
    if (pcApply.outputExists(proj.PROJECT_DIRECTORY.getValue() + ext.rootOf(output) + FILE_EXTs[0],
                             true)) {// warn about existence when checking
      pcApply.setExtrapolatedPCsFile(ext.rootOf(output) + FILE_EXTs[0]);
    } else {
      pcApply.setParams(params);
      pcApply.applyLoadings();
      pcApply.reportExtropolatedPCs(ext.rootOf(output) + FILE_EXTs[0]);
    }
    return pcApply;
  }

  /**
   * @param proj as above
   * @param pcFile either a raw or extrapolated principal component file
   * @param markersToAssessFile the markers to use for median values/resdiual reporting
   * @param numComponents number of components to regress out
   * @param printFull print all data associated with the median markers
   * @param gcThreshold threshold to use for median value computation
   * @param homozygousOnly only compute using homozygous calls
   * @param recomputeLRR recompute the log R ratios on the fly
   * @param output
   * @param log
   * @return
   */
  public static PrincipalComponentsResiduals computeResiduals(Project proj, String pcFile,
                                                              String markersToAssessFile,
                                                              int numComponents, boolean printFull,
                                                              float gcThreshold,
                                                              boolean homozygousOnly,
                                                              boolean recomputeLRR, String output,
                                                              GcAdjustorParameters params) {
    PrincipalComponentsResiduals pcResids = new PrincipalComponentsResiduals(proj, pcFile,
                                                                             markersToAssessFile,
                                                                             numComponents,
                                                                             printFull, gcThreshold,
                                                                             homozygousOnly,
                                                                             recomputeLRR,
                                                                             ext.rootOf(output));
    String file = proj.PROJECT_DIRECTORY.getValue() + ext.rootOf(ext.rootOf(output) + FILE_EXTs[1])
                  + PrincipalComponentsResiduals.MT_REPORT_EXT[0];
    if (!Files.exists(file)) {
      pcResids.setParams(params);// params can be null
      pcResids.computeAssessmentDataMedians();
      pcResids.computeResiduals();
      pcResids.computeInverseTransformResiduals();
      pcResids.summarize(ext.rootOf(output) + FILE_EXTs[1]);
    } else {
      pcResids.setResidOutput(ext.removeDirectoryInfo(file));
      proj.getLog()
          .report("Detected that the following report file already exists:\n" + output + "\n");
      proj.getLog().report("Skipping the report generation and using this file instead");
      proj.getLog()
          .report("If this is incorrect (using a different number of components, new samples, etc...),  please remove or change the name of the file listed above.\n Alternatively, specify a new analysis name");
    }
    return pcResids;
  }

  public static PrincipalComponentsApply generateFullPCA(Project proj, int numComponents,
                                                         String outputBase,
                                                         boolean recomputeLRR_PCs,
                                                         boolean imputeMeanForNaN,
                                                         GcAdjustorParameters params,
                                                         PRE_PROCESSING_METHOD method, Logger log) {
    log.report("\nReady to perform the principal components analysis (PCA)\n");
    PrincipalComponentsCompute pcs = PCA.computePrincipalComponents(proj, false, numComponents,
                                                                    false, method, true, true,
                                                                    imputeMeanForNaN,
                                                                    recomputeLRR_PCs,
                                                                    proj.PROJECT_DIRECTORY.getValue()
                                                                                      + outputBase
                                                                                      + PCA_SAMPLES,
                                                                    outputBase, params);

    if (pcs == null) {
      return null;
    }
    // apply PCs to everyone, we set useFile to null and excludeSamples to false to get all samples
    // in the current project.
    // TODO, if we ever want to apply to only a subset of the project, we can do that here.....
    log.report("\nApplying the loadings from the principal components analysis to all samples\n");
    PrincipalComponentsApply pcApply = PCA.applyLoadings(proj, numComponents,
                                                         pcs.getSingularValuesFile(),
                                                         pcs.getMarkerLoadingFile(), null, false,
                                                         imputeMeanForNaN, recomputeLRR_PCs,
                                                         outputBase, params);
    // Compute Medians for (MT) markers and compute residuals from PCs for everyone
    log.report("\nComputing residuals after regressing out " + numComponents
               + " principal component" + (numComponents == 1 ? "" : "s") + "\n");
    return pcApply;
  }

}
