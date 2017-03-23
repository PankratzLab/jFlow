package org.genvisis.gwas;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.genvisis.CLI;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.gwas.MarkerQC.QC_METRIC;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Maps;

public class FurtherAnalysisQc extends Qc {

	private static final Map<QC_METRIC, String> DEFAULT_ILLUMINA_MARKER_QC_THRESHOLDS;
	private static final Map<QC_METRIC, String> DEFAULT_AFFY_MARKER_QC_THRESHOLDS;

	public static final String FURTHER_ANALYSIS_DIR = "further_analysis_QC/";
	public static final String FURTHER_ANALYSIS_QC_PLINK_SUFFIX = "_QCd";

	public static final String SAMPLE_QC_DROPS = "mind_drops.dat";

	public static final String ARG_UNRELATEDS = "unrelateds";
	public static final String ARG_EUROPEANS = "europeans";

	private static final String BONFERRONI_CORRECTED_P_THRESHOLD = "<1E-7";

	static {
		Map<QC_METRIC, String> defaultIlluminaMetricThresholds = Maps.newEnumMap(QC_METRIC.class);
		defaultIlluminaMetricThresholds.put(QC_METRIC.CHR, "<1");
		defaultIlluminaMetricThresholds.put(QC_METRIC.MAF, "<0");
		defaultIlluminaMetricThresholds.put(QC_METRIC.CALLRATE,
																				MarkerQC.DEFAULT_ILLUMINA_CALLRATE_THRESHOLD);
		defaultIlluminaMetricThresholds.put(QC_METRIC.HWE, BONFERRONI_CORRECTED_P_THRESHOLD);
		defaultIlluminaMetricThresholds.put(QC_METRIC.MISHAP_HETERO, BONFERRONI_CORRECTED_P_THRESHOLD);
		defaultIlluminaMetricThresholds.put(QC_METRIC.MISHAP_MIN, BONFERRONI_CORRECTED_P_THRESHOLD);
		defaultIlluminaMetricThresholds.put(QC_METRIC.P_MISS, BONFERRONI_CORRECTED_P_THRESHOLD);
		defaultIlluminaMetricThresholds.put(QC_METRIC.P_GENDER, BONFERRONI_CORRECTED_P_THRESHOLD);
		defaultIlluminaMetricThresholds.put(QC_METRIC.P_GENDER_MISS, BONFERRONI_CORRECTED_P_THRESHOLD);

		DEFAULT_ILLUMINA_MARKER_QC_THRESHOLDS = Collections.unmodifiableMap(defaultIlluminaMetricThresholds);
	}
	static {
		Map<QC_METRIC, String> defaultAffyMetricThresholds = Maps.newEnumMap(DEFAULT_ILLUMINA_MARKER_QC_THRESHOLDS);
		defaultAffyMetricThresholds.put(QC_METRIC.CALLRATE, MarkerQC.DEFAULT_AFFY_CALLRATE_THRESHOLD);

		DEFAULT_AFFY_MARKER_QC_THRESHOLDS = Collections.unmodifiableMap(defaultAffyMetricThresholds);
	}

	private final String unrelatedsFile;
	private final String europeansFile;

	/**
	 * 
	 * @param sourceDir @see Qc#Qc(String, String, Map, Logger)
	 * @param plinkPrefix @see Qc#Qc(String, String, Map, Logger)
	 * @param markerQCThresholds @see Qc#Qc(String, String, Map, Logger)
	 * @param unrelatedsFile filename of list of unrelated FID/IID pairs to use for marker QC,
	 *        {@code null} to use all samples
	 * @param europeansFile filename of list of European samples to use for Hardy-Weinberg equilibrium
	 *        tests, {@code null} to use all unrelateds
	 * @param log @see Qc#Qc(String, String, Map, Logger)
	 */
	public FurtherAnalysisQc(String sourceDir, String plinkPrefix,
													 Map<QC_METRIC, String> markerQCThresholds,
													 String unrelatedsFile, String europeansFile,
													 Logger log) {
		super(sourceDir, plinkPrefix, markerQCThresholds, log);
		this.unrelatedsFile = unrelatedsFile;
		this.europeansFile = europeansFile;
	}

	public void runFurtherAnalysisQC() {
		final String subDir = FurtherAnalysisQc.FURTHER_ANALYSIS_DIR;
		if (!markerQc(subDir, unrelatedsFile, europeansFile))
			return;
		final String sourcePlink = sourceDir + plinkroot;
		final String markerDrops = qcDir + subDir + MARKER_QC_DROPS;
		final String plinkQCd = qcDir + subDir + plinkroot
														+ FurtherAnalysisQc.FURTHER_ANALYSIS_QC_PLINK_SUFFIX;
		List<String> applyQCCommand = ImmutableList.of("plink2", "--noweb", "--bfile", sourcePlink,
																									 "--exclude", markerDrops, "--mind", "0.05",
																									 "--make-bed", "--out", plinkQCd);
		Set<String> requiredInputs = PSF.Plink.getPlinkBedBimFamSet(sourcePlink);
		requiredInputs.add(markerDrops);
		Set<String> requiredOutputs = PSF.Plink.getPlinkBedBimFamSet(plinkQCd);
		if (CmdLine.runCommandWithFileChecks(applyQCCommand, "", requiredInputs,
																				 requiredOutputs,
																				 true, false, true, log)) {
			File iremFile = new File(plinkQCd + ".irem");
			File sampleDrops = new File(qcDir + subDir + SAMPLE_QC_DROPS);
			if (!iremFile.exists()) {
				// PLINK only generates .irem if there were samples to drop
				try {
					iremFile.createNewFile();
				} catch (IOException e) {
					log.reportError("Could not generate " + iremFile.getAbsolutePath());
					return;
				}
			}
			if (!iremFile.renameTo(sampleDrops)) {
				log.reportError("Could not move " + iremFile.getAbsolutePath() + " to "
												+ sampleDrops.getAbsolutePath());
			}
		}
	}

	public static Map<QC_METRIC, String> getDefaultMarkerQCThresholds(Project.ARRAY arrayType) {
		switch (arrayType) {
			case AFFY_GW6:
			case AFFY_GW6_CN:
				return DEFAULT_AFFY_MARKER_QC_THRESHOLDS;
			case ILLUMINA:
				return DEFAULT_ILLUMINA_MARKER_QC_THRESHOLDS;
			default:
				throw new IllegalArgumentException("Undefined for " + arrayType.getClass().getName() + ": "
																					 + arrayType.toString());
		}
	}

	public static void main(String[] args) {
		CLI c = new CLI(FurtherAnalysisQc.class);
		c.addArgWithDefault(CLI.ARG_INDIR, "directory with binary plink dataset", "./");
		c.addArgWithDefault(CLI.ARG_PLINKROOT, CLI.DESC_PLINKROOT, Qc.DEFAULT_PLINKROOT);
		c.addArg(ARG_UNRELATEDS,
						 "file of unrelated samples to use for marker QC, one FID/IID pair per line",
						 "unrelateds.txt");
		c.addArg(ARG_EUROPEANS,
						 "file of european samples to use for Hardy-Weinberg Equilibirum filtering, one FID/IID pair per line",
						 "europeans.txt");
		for (QC_METRIC metric : QC_METRIC.values()) {
			String defaultThreshold = DEFAULT_ILLUMINA_MARKER_QC_THRESHOLDS.get(metric);
			c.addArgWithDefault(metric.getKey(), metric.getCLIDescription(), defaultThreshold);
		}

		c.addArgWithDefault(CLI.ARG_LOG, CLI.DESC_LOG, "furtherAnalysisQC.log");

		c.parseWithExit(args);

		String dir = c.get(CLI.ARG_INDIR);
		String inputPlinkroot = c.get(CLI.ARG_PLINKROOT);
		String unrelatedsFile = c.get(ARG_UNRELATEDS);
		if (unrelatedsFile != null) {
			unrelatedsFile = new File(unrelatedsFile).getAbsolutePath();
		}
		String europeansFile = c.get(ARG_EUROPEANS);
		if (europeansFile != null) {
			europeansFile = new File(europeansFile).getAbsolutePath();
		}
		Map<QC_METRIC, String> markerQCThresholds = Maps.newEnumMap(QC_METRIC.class);
		for (QC_METRIC metric : QC_METRIC.values()) {
			markerQCThresholds.put(metric, c.get(metric.getKey()));
		}
		Logger log = new Logger(dir + c.get(CLI.ARG_LOG));


		try {
			new FurtherAnalysisQc(dir, inputPlinkroot, markerQCThresholds, unrelatedsFile, europeansFile,
														log).runFurtherAnalysisQC();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
