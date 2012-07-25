package cnv.plots;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Vector;

import javax.swing.JOptionPane;
import stats.LogisticRegression;
import common.Array;
import common.Files;
import cnv.filesys.MarkerData;
import cnv.filesys.MarkerDataCollection;
import cnv.filesys.MarkerLookup;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import cnv.var.SampleData;

public class ScanForProblemSnp {
	public static final int DEFAULT_GC_THRESHOLD = 15;
	private Project proj;
	private boolean jar;
//	private float gcThreshold;
//	private String[] samples;
//	private SampleData sampleData;
	private MarkerData[] markerData;
	private String[] markerList;
//	private String[] commentList;
	private MarkerLookup markerLookup;

	public ScanForProblemSnp (Project proj, String filename) {
		PrintWriter writer;
//		String[] samples;
//		Vector<double[]> xys, baflrrs; // raw_xys, thetars, 
//		float[] xs, ys, lrrs, bafs;
//		LogisticRegression lr;
		String output;
//        String[] files;
        int count, numModes, numGT5;
//        SampleData sampleData;
//        int[] sexes;
        
//        sampleData = proj.getSampleData(false);
//        samples = proj.getSamples();
//        sexes = new int[samples.length];
//        for (int i = 0; i < samples.length; i++) {
//        	sexes[i] = sampleData.getSexForIndividual(samples[i]);
//		}

        loadMarkListFromTxtFile();
        loadMarkerData(proj, markerList);
		try {
			writer = new PrintWriter(new FileWriter(proj.getDir(Project.RESULTS_DIRECTORY)+filename));
			writer.println("Marker\tChromosome\tCNP\tBlurryMargin\tAA\tAB\tBB");
			count = 0;
			numGT5 = 0;
			for (int j = 0; j < markerData.length; j++) {
				if (count++ % 1000 == 0) {
					System.out.println(count+"\t"+numGT5);
				}
				numModes = markerData[j].detectCNP();
				output = markerData[j].getMarkerName()
						+ "\t" + markerData[j].getChr()
						+ "\t" + numModes
						+ "\t" + Array.toStr(markerData[j].getAB_GenotypeCount())
						+ "\t.";
				writer.println(output);
				writer.flush();
				if (numModes > 5) {
					numGT5++;
				}
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing results");
			e.printStackTrace();
		}
        
        /*
        files = Files.list(proj.getDir(Project.PLOT_DIRECTORY), ".scat", false);		//zx
		try {
			writer = new PrintWriter(new FileWriter(proj.getDir(Project.RESULTS_DIRECTORY)+filename));
			writer.println("Marker\tChromosome\tCNP\tBlurryMargin\tAA\tAB\tBB");
			count = 0;
			numGT5 = 0;
			for (int i=0; i<files.length; i++) {
				markerData = MarkerDataCollection.load(proj.getDir(Project.PLOT_DIRECTORY)+files[i], proj.getJarStatus()).getCollection();
				for (int j = 0; j < markerData.length; j++) {
					if (count++ % 1000 == 0) {
						System.out.println(count+"\t"+numGT5);
					}
					numModes = markerData[j].detectCNP();
					output = markerData[j].getMarkerName()
							+ "\t" + markerData[j].getChr()
							+ "\t" + numModes
							+ "\t" + Array.toStr(markerData[j].getAB_GenotypeCount())
							+ "\t.";
					writer.println(output);
					writer.flush();
					if (numModes > 5) {
						numGT5++;
					}
				}
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing results");
			e.printStackTrace();
		}
		*/
		System.out.println("done");
	}

	public void loadMarkDataFromScatFiles() {
        String[] files;

        files = Files.list(proj.getDir(Project.PLOT_DIRECTORY), ".scat", false);
		try {
			for (int i=0; i<files.length; i++) {
				markerData = MarkerDataCollection.load(proj.getDir(Project.PLOT_DIRECTORY)+files[i], proj.getJarStatus()).getCollection();
			}
		} catch (Exception e) {
			System.err.println("Error writing results");
			e.printStackTrace();
		}
		System.out.println("done");
	}

	public void loadMarkListFromTxtFile() {
		BufferedReader reader;
		Vector<String> markerNames = new Vector<String>();
		Vector<String> markerComments = new Vector<String>();
		String[] line;
		String filename;
		Vector<String> missingMarkers;
		
		missingMarkers = new Vector<String>();
		filename = proj.getFilename(Project.DISPLAY_MARKERS_FILENAME);
		//System.out.println("filename: "+filename);//zx
		try {
			try {
				reader = Files.getReader(filename, jar, true, false);
				if (reader == null) {
					JOptionPane.showMessageDialog(null, "Failed to load '"+filename+"'", "Error", JOptionPane.ERROR_MESSAGE);
					return;
				}
				while (reader.ready()) {
					line = reader.readLine().trim().split("\t", -1);
					if (markerLookup.contains(line[0])) {
						markerNames.add(line[0]);
						markerComments.add(line.length>1?line[1]:"");
					} else {
//						System.err.println("Error - could not find "+line[0]+" in the lookup table");
						missingMarkers.add(line[0]);
					}
				}
				reader.close();
				if (missingMarkers.size() > 0) {
					JOptionPane.showMessageDialog(null, "Error - the following markers were not found in the MarkerSet:\n"+Array.toStr(Array.toStringArray(missingMarkers), "\n"), "Error", JOptionPane.ERROR_MESSAGE);
				}
			} catch (FileNotFoundException fnfe) {
				JOptionPane.showMessageDialog(null, "Error - could not find \""+filename+"\"", "Error", JOptionPane.ERROR_MESSAGE);
				// handle error by closing window
			} catch (Exception e) {
				System.err.println("Error reading file \""+filename+"\"");
				e.printStackTrace();
				System.exit(2);
			}

			markerList = Array.toStringArray(markerNames);
		} catch (Exception e) {
			System.err.println("Error loading: "+filename);
			e.printStackTrace();
		}
	}

	public void loadMarkerData(Project proj, String[] markerList) {
		markerData = MarkerSet.loadFromList(proj, markerList);
	}
	
//	public void main () {
//		loadData();
//		for (int i=0; i<; i++) {
//			clusterKMeans(marker[i]);
//			isCnvBimodal(marker[i]);
//			isCnvOther(marker[i]);
//		}
//	}
}
