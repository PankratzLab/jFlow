package org.genvisis.one.JL.shadow;

import java.util.ArrayList;

import org.genvisis.cnv.analysis.pca.PrincipalComponentsIntensity;
import org.genvisis.cnv.analysis.pca.PrincipalComponentsResiduals;
import org.genvisis.cnv.analysis.pca.PrincipalComponentsIntensity.PcCorrectionProducer;
import org.genvisis.cnv.filesys.MarkerData;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.TransposeData;
import org.genvisis.common.Files;
import org.genvisis.common.WorkerTrain;
import org.genvisis.stats.LeastSquares.LS_TYPE;

/**
 * Draft of a possibly a more friendly,less temporary version of shadowing
 *
 * Needs a few magic methods
 */
public class ShadowRework {
	private ShadowRework() {

	}

	/**
	 * @param proj
	 *            Project to correct
	 * @param principalComponentsResiduals
	 *            PCs to do the correcting
	 * @param preserveBafs
	 *            preserve BAF values (NGS specific), you likely want false here
	 * @param sampleSex
	 *            for Sex specific clustering
	 * @param samplesToUseCluster
	 *            samples to seed correction
	 * @param lType
	 * @param numComponents
	 *            number of PCs to correct for
	 * @param numCorrectionThreads
	 *            number of threads within a marker (max of 6 can be utilized)
	 * @param numMarkerThreads
	 *            number of markers corrected at once
	 */
	public static void correctProject(Project proj, PrincipalComponentsResiduals principalComponentsResiduals,
			boolean preserveBafs, int[] sampleSex, boolean[] samplesToUseCluster, int numComponents,
			int numCorrectionThreads, int numMarkerThreads) {

		Project shadowProject = new Project();
		// TODO update shadow project for new location of files,
		// transposed/samples dirs, etc

		String[] markers = proj.getMarkerNames(); // Correct the entire thing
		PcCorrectionProducer producer = new PcCorrectionProducer(principalComponentsResiduals, numComponents, sampleSex,
				samplesToUseCluster, LS_TYPE.REGULAR, numCorrectionThreads, 1, proj.getMarkerNames());
		WorkerTrain<PrincipalComponentsIntensity> train = new WorkerTrain<PrincipalComponentsIntensity>(producer,
				numMarkerThreads, 10, proj.getLog());
		ArrayList<String> notCorrected = new ArrayList<String>();
		int index = 0;

		// Magic marker data writer, should handle indexing, and verifying
		// order, properly splitting files, writes .mdRAF files

		// MarkerDataWriter markerDataWriter = new MarkerDataWriter(Project
		// proj);

		while (train.hasNext()) {
			PrincipalComponentsIntensity principalComponentsIntensity = train.next();
			MarkerData markerData = principalComponentsIntensity.getCentroidCompute().getMarkerData();
			if (principalComponentsIntensity.isFail()) {
				notCorrected.add(markers[index]);
			} else {
				byte[] abGenotypes = principalComponentsIntensity.getGenotypesUsed();// for
				// now
				float[][] correctedXY = principalComponentsIntensity
						.getCorrectedIntensity(PrincipalComponentsIntensity.XY_RETURN, true);
				float[][] correctedLRRBAF = principalComponentsIntensity
						.getCorrectedIntensity(PrincipalComponentsIntensity.BAF_LRR_RETURN, true);
				markerData = new MarkerData(markerData.getMarkerName(), markerData.getChr(), markerData.getPosition(),
						markerData.getFingerprint(), markerData.getGCs(), null, null, correctedXY[0], correctedXY[1],
						null, null, preserveBafs ? markerData.getBAFs() : correctedLRRBAF[0], correctedLRRBAF[1],
						abGenotypes, abGenotypes);
			}

			// TODO manage outliers

			// Magic method to write to mdRAF
			// markerDataWriter.write(markerData)

			// NOTE, we could actually write to all sampRAF files here as well,
			// but might be extra complexity and not worth it. Would just allow
			// to skip the reverse transpose

			index++;
		}
		if (!notCorrected.isEmpty()) {
			Files.writeArray(notCorrected.toArray(new String[notCorrected.size()]),
					shadowProject.PROJECT_DIRECTORY.getValue() + notCorrected.size()
							+ "_markersThatFailedCorrection.txt");
		}

		// Magic method to reverse transpose - if needed, not sure if this works
		TransposeData.reverseTranspose(proj);
	}

}
