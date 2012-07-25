package cnv.plots;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.Vector;

import javax.swing.JOptionPane;

import cnv.filesys.MarkerData;
import cnv.filesys.MarkerLookup;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;

import common.Array;
import common.Files;

public class Test_CnpDetection {
	private MarkerData[] markerData;
	private String[] markerList;
	String[] actualIsCnp;
	String[] actualClusters;
	private MarkerLookup markerLookup;

	public Test_CnpDetection(Project proj) {
		PrintWriter writer1, writer2;
		int numModes;
		int fpCount, fnCount, tpCount, tnCount;
		boolean isPositive;
		this.markerLookup = proj.getMarkerLookup();
		loadMarkerList(proj);
		loadMarkerData(proj, markerList);
		try {
			writer1 = new PrintWriter(new FileWriter(proj.getDir(Project.RESULTS_DIRECTORY)+"testResult_detail.txt"));
			writer1.println("a\tb\tMarker\tChromosome\tactualIsCnp\tmodelOuputIsCnp\tactualClusters\tmodelOuputnumModes\tBlurryMargin\tAA\tAB\tBB");
			writer2 = new PrintWriter(new FileWriter(proj.getDir(Project.RESULTS_DIRECTORY)+"testResult_summary.txt"));
//			writer2.println("Threshold1\tThreshold2\tFP_Rate\tnumFalseCases\tTN_Rate\tnumTrueCases");
			writer2.println("Threshold1\tThreshold2\tFP_Rate\tTN_Rate");
			for (double threshold1=0; threshold1<=0.5; threshold1=threshold1+.05) {
				System.out.print(".");
				for (double threshold2=0; threshold2<=0.5; threshold2=threshold2+.05) {
					fpCount = 0;
					fnCount = 0;
					tpCount = 0;
					tnCount = 0;
					isPositive = false;
					for (int i=0; i<markerList.length; i++) {
						numModes = markerData[i].detectCNP(threshold1, threshold2);
						//Customize Test Criteria Here
						if (actualIsCnp[i].equals("0")) {
							if (numModes>3
								|| (numModes>2 && actualClusters!=null && actualClusters.length>0 && actualClusters.equals("2"))
								|| (numModes>1 && actualClusters!=null && actualClusters.length>0 && actualClusters.equals("1"))) {
								isPositive = true;
								fpCount ++;
							} else {
								isPositive = false;
								fnCount ++;
							}
						} else {
							if (numModes>4
								|| (numModes>3 && actualClusters!=null && actualClusters.length>0 && actualClusters.equals("2"))
								|| (numModes>2 && actualClusters!=null && actualClusters.length>0 && actualClusters.equals("1"))) {
								isPositive = true;
								tpCount ++;
							} else {
								isPositive = false;
								tnCount ++;
							}
						}
						writer1.println((new DecimalFormat("#.##")).format(threshold1) + "\t" + (new DecimalFormat("#.##")).format(threshold2)
								+ "\t" + markerList[i]
								+ "\t" + markerData[i].getChr()
								+ "\t" + actualIsCnp[i]
								+ "\t" + (isPositive?"1":"0")
								+ "\t" + actualClusters[i]
								+ "\t" + numModes
								+ "\t" + "."
								+ "\t" + Array.toStr(markerData[i].getAB_GenotypeCount())
								);
						writer1.flush();
					}
					writer2.println((new DecimalFormat("#.##")).format(threshold1)
									+ "\t" + (new DecimalFormat("#.##")).format(threshold2)
									+ "\t" + (new DecimalFormat("#.####")).format((double)fpCount/(double)(fpCount+fnCount))
//									+ "\t" + (fpCount+fnCount)
									+ "\t" + (new DecimalFormat("#.####")).format((double)tnCount/(double)(tpCount+tnCount)));
//									+ "\t" + (tpCount+tnCount)
					writer2.flush();
				}
			}
			writer1.close();
			writer2.close();
			System.out.println("done");
		} catch(Exception e) {
					System.err.println("Error writing results");
					e.printStackTrace();
		}
	}

	public void loadMarkerList(Project proj) {
		BufferedReader reader;
		Vector<String> markerListTemp = new Vector<String>();
		Vector<String> actualIsCnpTemp = new Vector<String>();
		Vector<String> actualClustersTemp = new Vector<String>();
//		Vector<String> markerComments = new Vector<String>();
		String[] line;
		String filename;
		Vector<String> missingMarkers;
		
		missingMarkers = new Vector<String>();
//		filename = proj.getFilename(Project.DISPLAY_MARKERS_FILENAME);
		filename = proj.getProjectDir()+"data/testScanForCnp.txt";
		//System.out.println("filename: "+filename);//zx
		try {
			try {
				reader = Files.getReader(filename, proj.getJarStatus(), true, false);
				if (reader == null) {
					JOptionPane.showMessageDialog(null, "Failed to load '"+filename+"'", "Error", JOptionPane.ERROR_MESSAGE);
					return;
				}
				reader.readLine();
				while (reader.ready()) {
					line = reader.readLine().trim().split("\t", -1);
					if (markerLookup.contains(line[0])) {
						markerListTemp.add(line[0]);
						if (line.length>1) {
							actualIsCnpTemp.add(line[1]);
							if (line.length>2) {
								actualClustersTemp.add(line[2]);
							}
						}
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

			markerList = Array.toStringArray(markerListTemp);
			if (actualIsCnpTemp.size()>0) {
				actualIsCnp = Array.toStringArray(actualIsCnpTemp);
				if (actualClustersTemp.size()>0) {
					actualClusters = Array.toStringArray(actualClustersTemp);
				}
			}
//			markerData = MarkerSet.loadFromList(proj, markerList);
			//System.out.println("markerName: "+markerData.+);//zx
		} catch (Exception e) {
			System.err.println("Error loading: "+filename);
			e.printStackTrace();
		}
	}
	
	public void loadMarkerData(Project proj, String[] markerList) {
		markerData = MarkerSet.loadFromList(proj, markerList);
	}

}
