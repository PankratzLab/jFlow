package org.genvisis.gwas;

import java.io.File;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.genvisis.CLI;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.gwas.MarkerQC.QC_METRIC;
import org.genvisis.stats.Maths;

import com.google.common.base.Joiner;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Maps;

public class FurtherAnalysisQc extends Qc {

	public static final Map<QC_METRIC, String> DEFAULT_MARKER_QC_THRESHOLDS;

	public static final String FURTHER_ANALYSIS_DIR = "further_analysis_QC/";
	public static final String FURTHER_ANALYSIS_QC_PLINK_SUFFIX = "_QCd";

	public static final String ARG_UNRELATEDS = "unrelateds";
	public static final String ARG_EUROPEANS = "europeans";

	private static final String BONFERRONI_CORRECTED_P_THRESHOLD = "<1E-7";

	static {
		Map<QC_METRIC, String> defaultMetricThresholds = Maps.newEnumMap(QC_METRIC.class);
		defaultMetricThresholds.put(QC_METRIC.CHR, "<1");
		defaultMetricThresholds.put(QC_METRIC.MAF, "<0");
		defaultMetricThresholds.put(QC_METRIC.CALLRATE, "<0.98");
		defaultMetricThresholds.put(QC_METRIC.HWE, BONFERRONI_CORRECTED_P_THRESHOLD);
		defaultMetricThresholds.put(QC_METRIC.MISHAP_HETERO, BONFERRONI_CORRECTED_P_THRESHOLD);
		defaultMetricThresholds.put(QC_METRIC.MISHAP_MIN, BONFERRONI_CORRECTED_P_THRESHOLD);
		defaultMetricThresholds.put(QC_METRIC.P_MISS, BONFERRONI_CORRECTED_P_THRESHOLD);
		defaultMetricThresholds.put(QC_METRIC.P_GENDER, BONFERRONI_CORRECTED_P_THRESHOLD);
		defaultMetricThresholds.put(QC_METRIC.P_GENDER_MISS, BONFERRONI_CORRECTED_P_THRESHOLD);

		DEFAULT_MARKER_QC_THRESHOLDS = Collections.unmodifiableMap(defaultMetricThresholds);
	}

	private final Map<QC_METRIC, String> markerQCThresholds;
	private final String unrelatedsFile;
	private final String europeansFile;

	/**
	 * 
	 * @param dir @see Qc.#Qc(String, String, Logger)
	 * @param plinkPrefix @see Qc.#Qc(String, String, Logger)
	 * @param unrelatedsFile filename of list of unrelated FID/IID pairs to use for marker QC,
	 *        {@code null} to use all samples
	 * @param europeansFile filename of list of European samples to use for Hardy-Weinberg equilibrium
	 *        tests, {@code null} to use all unrelateds
	 * @param markerQCThresholds thresholds to apply for each desired marker QC metric, {@code null}
	 *        for defaults
	 * @param log @see Qc.#Qc(String, String, Logger)
	 */
	public FurtherAnalysisQc(String dir, String plinkPrefix, String unrelatedsFile,
													 String europeansFile, Map<QC_METRIC, String> markerQCThresholds,
													 Logger log) {
		super(dir, plinkPrefix, log);
		this.markerQCThresholds = markerQCThresholds;
		this.unrelatedsFile = unrelatedsFile;
		this.europeansFile = europeansFile;
	}

	private void runFurtherAnalysisQC() {
		final String subDir = FurtherAnalysisQc.FURTHER_ANALYSIS_DIR;
		final String plinkQCd = plink + FurtherAnalysisQc.FURTHER_ANALYSIS_QC_PLINK_SUFFIX;
		if (!markerQc(subDir, markerQCThresholds, unrelatedsFile, europeansFile))
			return;
		List<String> applyQCCommand = ImmutableList.of("plink2", "--noweb", "--bfile", plink,
																									 "--exclude", "miss_drops.dat", "--mind", "0.05",
																									 "--make-bed", "--out", plinkQCd);
		Set<String> requiredInputs = PSF.Plink.getPlinkBedBimFamSet(plink);
		requiredInputs.add("miss_drops.dat");
		Set<String> requiredOutputs = PSF.Plink.getPlinkBedBimFamSet(plinkQCd);
		CmdLine.runCommandWithFileChecks(applyQCCommand, dir + subDir, requiredInputs, requiredOutputs,
																		 true, false, true, log);
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
			StringBuilder descriptionBuilder = new StringBuilder();
			descriptionBuilder.append("threshold to reject ").append(metric.getKey()).append(", using: ");
			descriptionBuilder.append(Joiner.on(", ").join(Maths.OPERATORS));
			descriptionBuilder.append(" (<0 to not filter)");
			String description = descriptionBuilder.toString();
			String defaultThreshold = DEFAULT_MARKER_QC_THRESHOLDS.get(metric);
			c.addArgWithDefault(metric.getKey(), description, defaultThreshold);
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
			new FurtherAnalysisQc(dir, inputPlinkroot, unrelatedsFile, europeansFile, markerQCThresholds,
														log).runFurtherAnalysisQC();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
