package org.genvisis.one.JL.cushing;

import org.genvisis.cnv.filesys.MarkerSet.PreparedMarkerSet;

import java.io.File;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.ArrayList;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.genvisis.cnv.filesys.MarkerData;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.manage.MDL;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.SerializedFiles;
import org.genvisis.common.ext;
import org.genvisis.filesys.Segment;

import com.google.common.primitives.Ints;

/**
 * Specific to Cables1 cnv detection
 *
 */

public class Cables1 {
	private Cables1() {

	}

	private static class CNDataPoint {
		private String sample;
		private String marker;
		private double lrr;
		private double sampPdense;
		private double popPdense;
		private double popMean;
		private double popSd;
		private double sampMean;
		private double sampSd;
		private double sampPopLogProbDel;
		private double sampPopLogProbDup;

		private static final String[] Header = new String[] { "SAMPLE", "MARKER", "LRR", "SAMP_PROB", "POP_PROB",
				"POP_MEAN", "POP_SD", "SAMP_MEAN", "SAMP_SD", "SAMP_POP_LOG_DEL", "SAMP_POP_LOG_DUP" };

		/**
		 * @param sample
		 *            sample represented
		 * @param marker
		 *            current marker
		 * @param lrr
		 *            lrr value
		 * @param sampPdense
		 *            sample-relative probability
		 * @param popPdense
		 *            population-relative probabiliyt
		 * @param popMean
		 *            population mean of marker
		 * @param popSd
		 *            population sd of marker
		 * @param sampMean
		 *            sample mean lrr
		 * @param sampSd
		 *            sample sd
		 */
		public CNDataPoint(String sample, String marker, double lrr, double sampPdense, double popPdense,
				double popMean, double popSd, double sampMean, double sampSd) {
			super();
			this.sample = sample;
			this.marker = marker;
			this.lrr = lrr;
			this.sampPdense = sampPdense;
			this.popPdense = popPdense;
			this.popMean = popMean;
			this.popSd = popSd;
			this.sampMean = sampMean;
			this.sampSd = sampSd;
			this.sampPopLogProbDel = Math.log(1 - sampPdense) + Math.log(1 - popPdense);
			this.sampPopLogProbDup = Math.log(sampPdense) + Math.log(popPdense);

		}

		@Override
		public String toString() {
			return sample + "\t" + marker + "\t" + lrr + "\t" + sampPdense + "\t" + popPdense + "\t" + popMean + "\t"
					+ popSd + "\t" + sampMean + "\t" + sampSd + "\t" + sampPopLogProbDel + "\t" + sampPopLogProbDup;
		}

	}

	private static class SampleDistParams implements Serializable {
		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;
		// private String sampleName;
		private double mean;
		private double sd;

		private SampleDistParams(String sampleName, double mean, double sd) {
			super();
			// this.sampleName = sampleName;
			this.mean = mean;
			this.sd = sd;
		}

	}

	private static SampleDistParams[] generateSampleParams(Project proj, PreparedMarkerSet preparedMarkerSet,
			String outputDir, String[] pattersToExclude) {
		String serFile = outputDir + "sampleLrrParams.ser";
		if (Files.exists(serFile)) {
			return (SampleDistParams[]) SerializedFiles.readSerial(serFile, false, proj.getLog(), false, true);
		} else {
			String[] names = preparedMarkerSet.getMarkerNames();
			ArrayList<Integer> toUse = new ArrayList<>();
			for (int i = 0; i < names.length; i++) {
				boolean use = true;
				if (preparedMarkerSet.getChrs()[i] > 0 && preparedMarkerSet.getChrs()[i] < 23) {
					for (int j = 0; j < pattersToExclude.length; j++) {
						if (names[i].contains(pattersToExclude[j])) {
							use = false;
							break;
						}
					}
					if (use) {
						toUse.add(i);
					}
				}
			}
			proj.getLog().reportTimeInfo("Developing dists from " + toUse.size() + " markers");
			int[] indices = Ints.toArray(toUse);
			String[] samples = proj.getSamples();
			SampleDistParams[] sampleDistParams = new SampleDistParams[samples.length];
			StringBuilder raw = new StringBuilder("Sample\tLRR_mean\tLrrSD");
			for (int i = 0; i < samples.length; i++) {
				proj.getLog().reportTimeInfo("Sample " + i);
				Sample samp = proj.getFullSampleFromRandomAccessFile(samples[i]);
				float[] lrrs = Array.subArray(samp.getLRRs(), indices);
				float mean = Array.mean(lrrs, true);
				float sd = Array.stdev(lrrs, true);
				sampleDistParams[i] = new SampleDistParams(samples[i], mean, sd);
				raw.append("\n" + samples[i] + "\t" + mean + "\t" + sd);
				proj.getLog().reportTimeInfo("Sample " + i + "\t" + samples[i] + "\t" + mean + "\t" + sd);

			}
			SerializedFiles.writeSerial(sampleDistParams, serFile, true);
			Files.write(raw.toString(), ext.rootOf(serFile, false) + ".txt");
			return sampleDistParams;

		}

	}

	public static void main(String[] args) {
		Project proj = new Project(
				"/Users/Kitty/workspace.other/Genvisis/Genvisis/projects/cushings_corrected.properties", false);
		proj.verifyAndGenerateOutliers(true);
		Segment cables1Loc = new Segment("chr18:20,714,528-20,840,434");
		Segment alk = new Segment("chr2:28,961,923-31,735,067");
		Segment test = new Segment("chr6:306,447-338,866");
		Segment test2 = new Segment("chr3:141,874,465-142,094,208");
		Segment bai1 =new Segment("chr8:143,545,377-143,626,368");
		ArrayList<Segment> segs = new ArrayList<>();
		segs.add(cables1Loc);
		segs.add(alk);
		segs.add(test);
		segs.add(test2);
		segs.add(bai1);

		String outDir = proj.PROJECT_DIRECTORY.getValue() + "CablesCNVs/";

		new File(outDir).mkdirs();
		PreparedMarkerSet preparedMarkerSet = PreparedMarkerSet.getPreparedMarkerSet(proj.getMarkerSet());

		for (Segment seg : segs) {
			SampleDistParams[] sampleDistParams = generateSampleParams(proj, preparedMarkerSet, outDir,
					new String[] { "VARIANT_SITE", "OFF_TARGET" });
			String[] markersInSeg = preparedMarkerSet.getMarkersIn(seg, preparedMarkerSet.getIndicesByChr());

			MDL mdl = new MDL(proj, preparedMarkerSet, markersInSeg);

			ArrayList<CNDataPoint> points = new ArrayList<>();

			String[] samples = proj.getSamples();
			NormalDistribution[] ndsample = new NormalDistribution[samples.length];
			for (int i = 0; i < ndsample.length; i++) {
				ndsample[i] = new NormalDistribution(sampleDistParams[i].mean, sampleDistParams[i].sd);
			}

			while (mdl.hasNext()) {
				MarkerData md = mdl.next();
				if (!md.getMarkerName().contains("VARIANT_SITE") && !md.getMarkerName().contains("OFF_TARGET")) {
					proj.getLog().reportTimeInfo("Marker " + md.getMarkerName());

					NormalDistribution nd = new NormalDistribution(Array.mean(md.getLRRs(), true),
							Array.stdev(md.getLRRs(), true));

					for (int i = 0; i < samples.length; i++) {

						double sampProb = ndsample[i].density(md.getLRRs()[i]);
						double popProb = nd.density(md.getLRRs()[i]);
//						sampProb *= ndsample[i].getStandardDeviation();
//						popProb *= nd.getStandardDeviation();
						// if (md.getLRRs()[i] > 0) {
						// sampProb = 1 - sampProb;
						// popProb = 1 - popProb;
						// }

						CNDataPoint cnDataPointNorm = new CNDataPoint(samples[i], md.getMarkerName(), md.getLRRs()[i],
								sampProb, popProb, nd.getMean(), nd.getStandardDeviation(), ndsample[i].getMean(),
								ndsample[i].getStandardDeviation());

						points.add(cnDataPointNorm);
					}
				}
			}
			PrintWriter writer = Files.getAppropriateWriter(outDir + seg.getUCSClocation() + "cnvs.txt");
			writer.println(Array.toStr(CNDataPoint.Header));
			for (CNDataPoint cDataPoint : points) {
				writer.println(cDataPoint.toString());
			}
			writer.close();
		}
	}

}
