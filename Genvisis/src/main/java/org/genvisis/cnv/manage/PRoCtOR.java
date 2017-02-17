package org.genvisis.cnv.manage;

import java.io.File;
import java.io.FilenameFilter;

import org.genvisis.cnv.analysis.PennCNVPrep;
import org.genvisis.cnv.analysis.pca.PCA;
import org.genvisis.cnv.analysis.pca.PCAPrep;
import org.genvisis.cnv.analysis.pca.PrincipalComponentsApply;
import org.genvisis.cnv.analysis.pca.PrincipalComponentsIntensity.CHROMOSOME_X_STRATEGY;
import org.genvisis.cnv.analysis.pca.PrincipalComponentsIntensity.CORRECTION_TYPE;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.ext;
import org.genvisis.stats.LeastSquares.LS_TYPE;

// PRincipal COmponents Residuals - PR [o] C[t]O R
public class PRoCtOR {

	private static final String SHADOW_PREP_DIR = "shadowPrep/";

	private static long getMaxSampleSize(Project proj) {
		File[] sampleFiles = (new File(proj.SAMPLE_DIRECTORY.getValue())).listFiles(new FilenameFilter() {
			@Override
			public boolean accept(File dir, String name) {
				return name.endsWith(Sample.SAMPLE_FILE_EXTENSION);
			}
		});
		long max = 0;
		for (File f : sampleFiles) {
			if (f.length() > max) {
				max = f.length();
			}
		}
		return max;
	}

	private static int getSampleChunks(Project proj, int numThreads) {
		long mem = Runtime.getRuntime().maxMemory();
		long samp = getMaxSampleSize(proj);
		// With the assumption that all samples in a given chunk can be open when a chunk is processed
		// by a thread.
		// Thus we want a number of chunks such that we will not run out of memory if each thread is
		// simultaneously
		// processing chunks of (chunk size * max_chunk_size) files.
		// A fraction of max memory is used to account for additional files which may need to be in
		// memory during
		// this process (e.g. marker data)
		double sampleChunks = (0.65 * mem) / (numThreads * samp);
		return (int) sampleChunks;
	}

	public static String shadow(Project proj, String tmpDir, String outputBase,
															double markerCallRateFilter, boolean recomputeLRR_PCs,
															CORRECTION_TYPE correctionType, CHROMOSOME_X_STRATEGY strategy,
															int numComponents, int totalThreads) {
		int numMarkerThreads = 1;
		int numThreads = (int) Math.ceil((double) totalThreads / (double) numMarkerThreads);
		boolean markerQC = true;
		String useFile = null;
		int sampleChunks = getSampleChunks(proj, numThreads);
		proj.getLog().report("Using " + sampleChunks + " sample chunks");

		int retCode = PCAPrep.prepPCA(proj, numThreads, outputBase, markerQC, markerCallRateFilter,
																	useFile, proj.getSampleList(), proj.getLog());
		if (retCode != 42) {
			return PCAPrep.errorMessage(retCode);
		}
		PrincipalComponentsApply pcApply = PCA.generateFullPCA(proj, numComponents, outputBase,
																													 recomputeLRR_PCs, true, null,
																													 proj.getLog());
		proj.getLog().reportTime("Setting PCs file: " + pcApply.getExtrapolatedPCsFile());
		proj.INTENSITY_PC_FILENAME.setValue(pcApply.getExtrapolatedPCsFile());

		if (correctionType == CORRECTION_TYPE.GENERATE_PCS_ONLY) {
			return "";
		}
		PennCNVPrep.prepExport(proj, SHADOW_PREP_DIR, tmpDir, numComponents, null, numThreads,
													 numMarkerThreads, LS_TYPE.REGULAR, false, correctionType, strategy);
		PennCNVPrep.exportSpecialPennCNV(proj, SHADOW_PREP_DIR, tmpDir, numComponents, null, numThreads,
																		 numMarkerThreads, true, LS_TYPE.REGULAR, sampleChunks, false,
																		 false, correctionType, strategy);
		return "";
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "D:/projects/Poynter.properties";
		String tempDir = null;
		String outputBase = MitoPipeline.FILE_BASE;
		double callrate = MitoPipeline.DEFAULT_MKR_CALLRATE_FILTER;
		boolean recomputeLRR = false;
		int numComponents = MitoPipeline.DEFAULT_NUM_COMPONENTS;
		CORRECTION_TYPE correctionType = CORRECTION_TYPE.XY;
		CHROMOSOME_X_STRATEGY strategy = CHROMOSOME_X_STRATEGY.BIOLOGICAL;
		int numThreads = Runtime.getRuntime().availableProcessors();

		String usage = "\n" + "cnv.manage.PRoCtOR requires 0-1 arguments\n"
									 + "   (1) project properties filename (i.e. proj=" + filename + " (default))\n"
									 + "   (2) Number of principal components for correction (i.e. numComponents="
									 + numComponents + " (default))\n"
									 + "   (3) Output file full path and baseName for principal components correction files (i.e. outputBase="
									 + outputBase + " (default))\n"
									 + "   (4) Call-rate filter for determining high-quality markers (i.e. callrate="
									 + callrate + " (default))\n"
									 + "   (5) Flag specifying whether or not to re-compute Log-R Ratio values (usually false if LRRs already exist) (i.e. recomputeLRR="
									 + recomputeLRR + " (default))\n"
									 + "   (6) Type of correction.  Options include: "
									 + ArrayUtils.toStr(CORRECTION_TYPE.values(), ", ") + " (i.e. type="
									 + correctionType + " (default))\n"

									 + "   (7) Chromosome X correction strategy.  Options include: "
									 + ArrayUtils.toStr(CHROMOSOME_X_STRATEGY.values(), ", ") + " (i.e. sexStrategy="
									 + strategy + " (default))\n"

									 + "   (8) Total number of threads to use (i.e. numThreads=" + numThreads
									 + " (default))\n"

									 + "   (8) OPTIONAL: temp directory for intermediate files (which tend to be very large) (i.e. tmp="
									 + tempDir + " (default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("proj=")) {
				filename = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("numComponents=")) {
				numComponents = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith("outputBase=")) {
				outputBase = ext.parseStringArg(arg, outputBase);
				numArgs--;
			} else if (arg.startsWith("callrate=")) {
				callrate = ext.parseDoubleArg(arg);
				numArgs--;
			} else if (arg.startsWith("recomputeLRR=")) {
				recomputeLRR = ext.parseBooleanArg(arg);
				numArgs--;
			} else if (arg.startsWith("numThreads=")) {
				numThreads = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith("type=")) {
				correctionType = CORRECTION_TYPE.valueOf(ext.parseStringArg(arg,
																																		correctionType.toString()));
				numArgs--;
			} else if (arg.startsWith("sexStrategy=")) {
				strategy = CHROMOSOME_X_STRATEGY.valueOf(ext.parseStringArg(arg, strategy.toString()));
				numArgs--;
			} else if (arg.startsWith("tmp=")) {
				tempDir = ext.parseStringArg(arg, null);
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + arg);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			Project proj = new Project(filename, false);
			String err = shadow(proj, tempDir, outputBase, callrate, recomputeLRR, correctionType,
													strategy, numComponents, numThreads);
			if (!"".equals(err)) {
				System.err.println("Error - " + err);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
