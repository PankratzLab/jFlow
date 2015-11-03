package one.JL;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Set;
import java.util.concurrent.Callable;

import common.Array;
import common.WorkerTrain;
import common.WorkerTrain.Producer;
import common.ext;
import cnv.analysis.MosaicismDetect;
import cnv.analysis.MosaicismDetect.MosaicBuilder;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import cnv.filesys.Sample;
import cnv.var.CNVariant;
import cnv.var.LocusSet;
import cnv.var.MosaicRegion;
import cnv.var.SampleData;

public class CNVMosaic {

	private static void computeCNVMosaic(Project proj) {
		String[] cnvFiles = proj.CNV_FILENAMES.getValue();
		for (int i = 0; i < cnvFiles.length; i++) {
			LocusSet<CNVariant> cnvs = CNVariant.loadLocSet(cnvFiles[i], proj.getLog());
			Hashtable<String, LocusSet<CNVariant>> indSets = CNVariant.breakIntoInds(cnvs, proj.getLog());
			MosaicForceProducer producer = new MosaicForceProducer(proj, indSets);
			WorkerTrain<MosaicRegion[]> train = new WorkerTrain<MosaicRegion[]>(producer, proj.NUM_THREADS.getValue(), 2, proj.getLog());

			String output = ext.addToRoot(cnvFiles[i], ".mosaicMetrics");
			try {
				PrintWriter writer = new PrintWriter(new FileWriter(output));
				writer.println(Array.toStr(CNVariant.PLINK_CNV_HEADER) + "\t" + Array.toStr(MosaicRegion.ADD_HEADER));
				int index = 0;
				while (train.hasNext()) {
					MosaicRegion[] regions = train.next();
					for (int j = 0; j < regions.length; j++) {
						writer.println(regions[i].toAnalysisString());
					}
					index++;
					proj.getLog().reportTimeInfo(index + " of " + indSets.size());
				}
				writer.close();
			} catch (Exception e) {
				proj.getLog().reportError("Error writing to " + output);
				proj.getLog().reportException(e);
			}

		}
	}

	private static class MosaicForceProducer implements Producer<MosaicRegion[]> {
		private Project proj;
		private Hashtable<String, LocusSet<CNVariant>> indSets;
		private int index;
		private String[] fidsIids;
		private SampleData sampleData;
		private MarkerSet markerSet;

		public MosaicForceProducer(Project proj, Hashtable<String, LocusSet<CNVariant>> indSets) {
			super();
			this.proj = proj;
			this.indSets = indSets;
			this.index = 0;
			Set<String> tmp = indSets.keySet();
			this.fidsIids = new String[tmp.size()];
			int i = 0;
			for (String fidIid : tmp) {
				fidsIids[i] = fidIid;
				i++;
			}
			this.sampleData = proj.getSampleData(0, false);
			this.markerSet = proj.getMarkerSet();

		}

		@Override
		public boolean hasNext() {

			return index < fidsIids.length;
		}

		@Override
		public Callable<MosaicRegion[]> next() {
			final String sample = sampleData.lookup(fidsIids[index])[0];
			final LocusSet<CNVariant> sampCNVs = indSets.get(fidsIids[index]);

			Callable<MosaicRegion[]> callable = new Callable<MosaicRegion[]>() {

				@Override
				public MosaicRegion[] call() throws Exception {

					Sample samp = proj.getFullSampleFromRandomAccessFile(sample);
					ArrayList<MosaicRegion> all = new ArrayList<MosaicRegion>();
					MosaicBuilder builderMosaic = new MosaicBuilder();
					builderMosaic.verbose(true);
					MosaicismDetect md = builderMosaic.build(proj, sample, markerSet, Array.toDoubleArray(samp.getBAFs()));
					for (int i = 0; i < sampCNVs.getLoci().length; i++) {
						LocusSet<MosaicRegion> mosSet = md.callMosaic(sampCNVs.getLoci()[i], true);
						if (mosSet.getLoci().length > 0) {

							for (int j = 0; j < mosSet.getLoci().length; j++) {

								MosaicRegion recaptured = new MosaicRegion(sampCNVs.getLoci()[i], mosSet.getLoci()[j]);
								all.add(recaptured);
							}
						} else {
							all.add(new MosaicRegion(sampCNVs.getLoci()[i], Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN));
						}
					}
					if (all.size() != sampCNVs.getLoci().length) {
						throw new IllegalStateException("Did not call all cnvs for sample " + sample);
					}
					return all.toArray(new MosaicRegion[all.size()]);
				}
			};
			index++;
			return callable;
		}

		@Override
		public void shutdown() {
			// TODO Auto-generated method stub

		}

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;

		String usage = "\n" + "one.JL.CNVMosaic requires 0-1 arguments\n" + "   (1) project  (i.e. proj=" + filename + " (default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
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
			computeCNVMosaic(proj);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
