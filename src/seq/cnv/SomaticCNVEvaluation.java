package seq.cnv;

import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Set;
import java.util.concurrent.Callable;

import common.Logger;
import common.WorkerTrain;
import common.WorkerTrain.Producer;
import seq.manage.VCFOps.VcfPopulation;
import seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import cnv.analysis.BeastScore;
import cnv.filesys.Project;
import cnv.filesys.Sample;
import cnv.filesys.MarkerSet.PreparedMarkerSet;
import cnv.var.CNVariant;
import cnv.var.LocusSet;
import cnv.var.SampleData;

/**
 * Class to refine somatic cnv calls
 */
public class SomaticCNVEvaluation {

	private static void filter(Project proj, String vpopFile, String cnvFile, int numThreads) {
		Logger log = proj.getLog();
		VcfPopulation vpop = VcfPopulation.load(vpopFile, POPULATION_TYPE.TUMOR_NORMAL, log);
		vpop.report();
		// String[] markFilter = proj.MARKER_COLOR_KEY_FILENAMES.getValue();
		// log.reportTimeInfo(markFilter.length + " filter files detected");

		LocusSet<CNVariant> cnvs = CNVariant.loadLocSet(cnvFile, log);
		Hashtable<String, LocusSet<CNVariant>> inds = CNVariant.breakIntoInds(cnvs, log);
		Set<String> tumors = vpop.getTumorSamples();
		ArrayList<TNCNV> tncnvs = new ArrayList<SomaticCNVEvaluation.TNCNV>();
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

				tncnvs.add(new TNCNV(proj, inds.get(fidIid), tumor, normal));
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
		TNCNVProducer producer = new TNCNVProducer(PreparedMarkerSet.getPreparedMarkerSet(proj.getMarkerSet()), tncnvs);
		WorkerTrain<TNCNV> train = new WorkerTrain<SomaticCNVEvaluation.TNCNV>(producer, numThreads, numThreads, proj.getLog());
		while (train.hasNext()) {
			System.out.println("DHFSD");
			train.next();
		}
	}

	private static class TNCNVProducer implements Producer<TNCNV> {
		private PreparedMarkerSet markerSet;
		private ArrayList<TNCNV> tncnvs;
		private int index;

		public TNCNVProducer(PreparedMarkerSet markerSet, ArrayList<TNCNV> tncnvs) {
			super();
			this.tncnvs = tncnvs;
		}

		@Override
		public boolean hasNext() {
			// TODO Auto-generated method stub
			return index < tncnvs.size();
		}

		@Override
		public Callable<TNCNV> next() {
			TNCNV current = tncnvs.get(index);
			index++;
			return current;
		}

		@Override
		public void shutdown() {
			// TODO Auto-generated method stub

		}

	}

	private static class TNCNV implements Callable<TNCNV> {
		private Project proj;
		private LocusSet<CNVariant> tumorCnvs;
		private String tumorSample;
		private String normalSample;
		private PreparedMarkerSet markerSet;

		public TNCNV(Project proj, LocusSet<CNVariant> tumorCnvs, String tumorSample, String normalSample) {
			super();
			this.proj = proj;
			this.tumorCnvs = tumorCnvs;
			this.tumorSample = tumorSample;
			this.normalSample = normalSample;
		}

		@Override
		public TNCNV call() throws Exception {
			markerSet = PreparedMarkerSet.getPreparedMarkerSet(proj.getMarkerSet());

			proj.getLog().reportTimeInfo("Loading" + tumorSample);

			Sample tumorSamp = proj.getFullSampleFromRandomAccessFile(tumorSample);
			proj.getLog().reportTimeInfo("Loading" + normalSample);

			Sample normalSamp = proj.getFullSampleFromRandomAccessFile(normalSample);

			int[][] cnvIndices = new int[tumorCnvs.getLoci().length][];
			Hashtable<String, Integer> track = proj.getMarkerIndices();
			for (int i = 0; i < tumorCnvs.getLoci().length; i++) {
				String[] namesIn = markerSet.getMarkersIn(tumorCnvs.getLoci()[i], markerSet.getIndicesByChr());
				cnvIndices[i] = new int[namesIn.length];
				for (int j = 0; j < namesIn.length; j++) {
					cnvIndices[i][j] = track.get(namesIn[j]);
				}
			}

			proj.getLog().reportTimeInfo("Computing scores for " + tumorSample);
			BeastScore beastScoreTumor = new BeastScore(tumorSamp.getLRRs(), markerSet.getIndicesByChr(), cnvIndices, proj.getLog());
			//beastScoreTumor.setUse(useForMedian);
			beastScoreTumor.computeBeastScores();
			proj.getLog().reportTimeInfo("Computing scores for " + normalSample);
			BeastScore beastScoreNormal = new BeastScore(normalSamp.getLRRs(), markerSet.getIndicesByChr(), cnvIndices, proj.getLog());
			beastScoreNormal.computeBeastScores();
			this.markerSet = null;
			return this;
		}

	}

	public static void main(String[] args) {
		Project proj = new Project("C:/workspace/Genvisis/projects/Cushings.properties", false);
		String[] cnvFiles = proj.CNV_FILENAMES.getValue();
		String vpopFile = proj.PROJECT_DIRECTORY.getValue() + "TN.vpop";
		int numthreads = 5;

		for (int i = 0; i < cnvFiles.length; i++) {
			filter(proj, vpopFile, cnvFiles[i], numthreads);

		}

	}

}
