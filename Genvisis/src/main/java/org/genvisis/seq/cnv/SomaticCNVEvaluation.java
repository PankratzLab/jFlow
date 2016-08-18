package org.genvisis.seq.cnv;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Set;
import java.util.concurrent.Callable;

import org.genvisis.cnv.analysis.BeastScore;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.filesys.MarkerSet.PreparedMarkerSet;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.WorkerTrain;
import org.genvisis.common.ext;
import org.genvisis.common.WorkerTrain.AbstractProducer;
import org.genvisis.filesys.CNVariant;
import org.genvisis.filesys.GeneData;
import org.genvisis.filesys.GeneTrack;
import org.genvisis.filesys.LocusSet;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.manage.BamImport.NGS_MARKER_TYPE;
import org.genvisis.seq.manage.VCFOps.VcfPopulation;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;

import com.google.common.primitives.Ints;

/**
 * Class to refine somatic cnv calls
 */
public class SomaticCNVEvaluation {

	private static class TNTrack {
		private String tumorFidIid;
		private String tumorSample;
		private String normalSample;

		private TNTrack(String tumorFidIid, String tumorSample, String normalSample) {
			super();
			this.tumorFidIid = tumorFidIid;
			this.tumorSample = tumorSample;
			this.normalSample = normalSample;
		}

		private String getTumorFidIid() {
			return tumorFidIid;
		}

		private String getTumorSample() {
			return tumorSample;
		}

		private String getNormalSample() {
			return normalSample;
		}

	}

	private static void filter(Project proj, String vpopFile, String cnvFile, double normalCutoff, double diffCutoff, int numThreads) {
		Logger log = proj.getLog();
		VcfPopulation vpop = VcfPopulation.load(vpopFile, POPULATION_TYPE.TUMOR_NORMAL, log);
		vpop.report();
		// String[] markFilter = proj.MARKER_COLOR_KEY_FILENAMES.getValue();
		// log.reportTimeInfo(markFilter.length + " filter files detected");

		LocusSet<CNVariant> cnvs = CNVariant.loadLocSet(cnvFile, log);
		Hashtable<String, LocusSet<CNVariant>> inds = CNVariant.breakIntoInds(cnvs, log);
		Set<String> tumors = vpop.getTumorSamples();
		ArrayList<TNTrack> trackers = new ArrayList<TNTrack>();
		SampleData sampleData = proj.getSampleData(0, false);

		for (String tnPair : vpop.getSubPop().keySet()) {
			Set<String> pair = vpop.getSubPop().get(tnPair);
			String tumor = null;
			String normal = null;
			for (String samp : pair) {
				if (tumors.contains(samp)) {
					tumor = samp;
				} else {
					normal = samp;
				}
			}

			String fidIid = sampleData.lookup(tumor)[1];
			if (inds.containsKey(fidIid)) {
				trackers.add(new TNTrack(fidIid, tumor, normal));
			}
		}
		// ArrayList<ColorManager<String>> managers = new ArrayList<ColorExt.ColorManager<String>>();
		// if (markFilter != null) {
		// for (int i = 0; i < markFilter.length; i++) {
		// proj.getLog().reportTimeInfo("Loading " + markFilter[i]);
		// managers.add(ColorExt.getColorManager(proj, markFilter[i]));
		//
		// }
		// }
		BeastFilt beastFilt = new BeastFilt(normalCutoff, diffCutoff);
		TNCNVProducer producer = new TNCNVProducer(proj, inds, trackers, beastFilt);
		WorkerTrain<TNCNV> train = new WorkerTrain<SomaticCNVEvaluation.TNCNV>(producer, numThreads, numThreads, proj.getLog());
		String outDir = proj.PROJECT_DIRECTORY.getValue() + "SomaticCNV/";
		new File(outDir).mkdirs();
		String outFile = outDir + ext.rootOf(cnvFile) + ".somaticEvals.txt";
		if (!Files.exists(outFile)) {
			System.exit(1);
			try {
				PrintWriter writer = new PrintWriter(new FileWriter(outFile));
				writer.println(Array.toStr(CNVariant.PLINK_CNV_HEADER) + "\t" + Array.toStr(SomaticEvaluation.HEADER));
				while (train.hasNext()) {
					TNCNV current = train.next();
					for (int i = 0; i < current.getTumorCnvs().getLoci().length; i++) {
						CNVariant currentCNV = current.getTumorCnvs().getLoci()[i];
						writer.println(currentCNV.toPlinkFormat() + "\t" + Array.toStr(current.getSomaticEvaluations()[i].getSummary()) + "\t" + Array.toStr(currentCNV.toTrailerFormat()));
					}
				}
				writer.close();
			} catch (Exception e) {
				log.reportError("Error writing to " + outFile);
				log.reportException(e);
			}
		}
		String[][] results = HashVec.loadFileToStringMatrix(outFile, false, null, false);
		String geneOutFile = ext.addToRoot(outFile, ".genes");
		try {
			GeneTrack geneTrack = GeneTrack.load(proj.getGeneTrackFilename(true), false);
			PrintWriter writer = new PrintWriter(new FileWriter(geneOutFile));
			writer.println(Array.toStr(results[0]) + "\tGENE");
			for (int i = 1; i < results.length; i++) {
				Segment current = new Segment(Byte.parseByte(results[i][2]), Integer.parseInt(results[i][3]), Integer.parseInt(results[i][4]));
				GeneData[] geneDatas = geneTrack.getOverlappingGenes(current);
				HashSet<String> written = new HashSet<String>();
				for (int j = 0; j < geneDatas.length; j++) {
					if (!written.contains(geneDatas[j].getGeneName())) {
						writer.println(Array.toStr(results[i]) + "\t" + geneDatas[j].getGeneName());
						written.add(geneDatas[j].getGeneName());
					}
				}
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + outFile);
			log.reportException(e);
		}
	}

	private static class BeastFilt {
		private double normalCutoff;
		private double diffCutoff;

		private double getNormalCutoff() {
			return normalCutoff;
		}

		private double getDiffCutoff() {
			return diffCutoff;
		}

		public BeastFilt(double normalCutoff, double diffCutoff) {
			super();
			this.normalCutoff = normalCutoff;
			this.diffCutoff = diffCutoff;
		}

	}

	private static class TNCNVProducer extends AbstractProducer<TNCNV> {
		private Project proj;

		private Hashtable<String, LocusSet<CNVariant>> inds;
		private ArrayList<TNTrack> tncnvs;
		private int index;
		private BeastFilt beastFilt;

		public TNCNVProducer(Project proj, Hashtable<String, LocusSet<CNVariant>> inds, ArrayList<TNTrack> tncnvs, BeastFilt beastFilt) {
			super();
			this.proj = proj;
			this.inds = inds;
			this.tncnvs = tncnvs;
			this.index = 0;
			this.beastFilt = beastFilt;
		}

		@Override
		public boolean hasNext() {
			// TODO Auto-generated method stub
			return index < tncnvs.size();
		}

		@Override
		public Callable<TNCNV> next() {
			TNTrack trackCurrent = tncnvs.get(index);
			TNCNV current = new TNCNV(proj, inds.get(trackCurrent.getTumorFidIid()), trackCurrent.getTumorSample(), trackCurrent.getNormalSample(), beastFilt);
			index++;
			return current;
		}
	}

	private static class SomaticEvaluation {
		private static final String[] HEADER = new String[] { "TYPE", "BEAST_HEIGHT_TUMOR", "BEAST_HEIGHT_NORMAL", "BEAST_HEIGHT_DIFF", "HQ" };
		private String type;
		private double bhTumor;
		private double bhNormal;
		private double heightDiff;
		private boolean hq;

		public SomaticEvaluation(String type, double bhTumor, double bhNormal) {
			super();
			this.type = type;
			this.bhTumor = bhTumor;
			this.bhNormal = bhNormal;
		}

		private void determineHQ(double normalCutoff, double diffCutoff) {
			this.hq = true;
			if (Math.abs(bhNormal) > normalCutoff) {// might remove since somatic can get weird
				hq = false;
			}
			boolean posN = bhNormal > 0;
			boolean posT = bhTumor > 0;

			boolean same = posN == posT;

			if (same) {
				this.heightDiff = Math.abs(bhTumor) - Math.abs(bhNormal);
				if (heightDiff < diffCutoff) {
					hq = false;
				}

			} else {
				this.heightDiff = Math.abs(bhTumor);
				if (heightDiff < diffCutoff) {
					hq = false;

				}
			}

		}

		private String[] getSummary() {
			ArrayList<String> summary = new ArrayList<String>();
			summary.add(type);
			summary.add(bhTumor + "");
			summary.add(bhNormal + "");
			summary.add(heightDiff + "");
			summary.add(hq + "");
			return Array.toStringArray(summary);
		}

	}

	private static class TNCNV implements Callable<TNCNV> {
		private Project proj;
		private LocusSet<CNVariant> tumorCnvs;
		private String tumorSample;
		private String normalSample;
		private PreparedMarkerSet markerSet;
		private SomaticEvaluation[] somaticEvaluations;
		private BeastFilt beastFilt;

		public TNCNV(Project proj, LocusSet<CNVariant> tumorCnvs, String tumorSample, String normalSample, BeastFilt beastFilt) {
			super();
			this.proj = proj;
			this.tumorCnvs = tumorCnvs;
			this.tumorSample = tumorSample;
			this.normalSample = normalSample;
			this.beastFilt = beastFilt;
		}

		private LocusSet<CNVariant> getTumorCnvs() {
			return tumorCnvs;
		}

		private SomaticEvaluation[] getSomaticEvaluations() {
			return somaticEvaluations;
		}

		@Override
		public TNCNV call() throws Exception {
			this.somaticEvaluations = new SomaticEvaluation[tumorCnvs.getLoci().length];

			markerSet = PreparedMarkerSet.getPreparedMarkerSet(proj.getMarkerSet());

			proj.getLog().reportTimeInfo("Loading" + tumorSample);

			Sample tumorSamp = proj.getFullSampleFromRandomAccessFile(tumorSample);
			proj.getLog().reportTimeInfo("Loading" + normalSample);

			Sample normalSamp = proj.getFullSampleFromRandomAccessFile(normalSample);

			int[][] cnvIndices = new int[tumorCnvs.getLoci().length][];
			Hashtable<String, Integer> track = proj.getMarkerIndices();
			for (int i = 0; i < tumorCnvs.getLoci().length; i++) {
				String[] namesIn = markerSet.getMarkersIn(tumorCnvs.getLoci()[i], markerSet.getIndicesByChr());
				ArrayList<Integer> nonVariant = new ArrayList<Integer>();

				for (int j = 0; j < namesIn.length; j++) {
					NGS_MARKER_TYPE type = NGS_MARKER_TYPE.getType(namesIn[j]);
					if (type != NGS_MARKER_TYPE.VARIANT_SITE) {
						nonVariant.add(track.get(namesIn[j]));
					}
				}
				cnvIndices[i] = Ints.toArray(nonVariant);
			}

			proj.getLog().reportTimeInfo("Computing scores for " + tumorSample);
			BeastScore beastScoreTumor = new BeastScore(tumorSamp.getLRRs(), markerSet.getIndicesByChr(), cnvIndices, proj.getLog());
			// beastScoreTumor.setUse(useForMedian);
			beastScoreTumor.computeBeastScores();
			proj.getLog().reportTimeInfo("Computing scores for " + normalSample);
			BeastScore beastScoreNormal = new BeastScore(normalSamp.getLRRs(), markerSet.getIndicesByChr(), cnvIndices, proj.getLog());
			beastScoreNormal.computeBeastScores();
			this.markerSet = null;
			for (int i = 0; i < cnvIndices.length; i++) {
				double bhTumor = beastScoreTumor.getBeastHeights()[i];
				double bhNormal = beastScoreNormal.getBeastHeights()[i];
				SomaticEvaluation somaticEvaluation = new SomaticEvaluation("ALL_MARKERS", bhTumor, bhNormal);
				somaticEvaluation.determineHQ(beastFilt.getNormalCutoff(), beastFilt.getDiffCutoff());
				somaticEvaluations[i] = somaticEvaluation;
			}
			return this;
		}
	}

	public static void main(String[] args) {
		Project proj = new Project("C:/workspace/Genvisis/projects/CushingsPCCorrected.properties", false);
		String[] cnvFiles = proj.CNV_FILENAMES.getValue();
		String vpopFile = proj.PROJECT_DIRECTORY.getValue() + "TN.vpop";
		int numthreads = 4;
		double normalCutoff = .25;
		double diffCutoff = .5;
		filter(proj, vpopFile, cnvFiles[0], normalCutoff, diffCutoff, numthreads);

		// for (int i = 0; i < cnvFiles.length; i++) {
		// }

	}

}
