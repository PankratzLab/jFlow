import java.awt.FlowLayout;
import java.awt.Font;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.Vector;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.WindowConstants;

import common.Array;
import cnv.filesys.Centroids;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;

public class lab {
	
	private static void mockupGUI() {
		JFrame frame = new JFrame();
		frame.setLayout(new FlowLayout());
		JCheckBox chkBx = new JCheckBox("Derive from Centroids:");
		chkBx.setFont(new Font("Arial", 0, 14));
		JComboBox<String> comboBx = new JComboBox<String>(new String[]{"sexSpecific_Male", "sexSpecific_Female"});
		comboBx.setFont(new Font("Arial", 0, 14));
		frame.add(chkBx);
		frame.add(comboBx);
		frame.add(new JButton());
		frame.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
		frame.setVisible(true);
	}
	
	
	private static void fakeAutoCentFile(Project proj, String centFile) {
		byte[] markerChrs;
		String[] allMarkers;
		MarkerSet ms;
		PrintWriter writer;
		Vector<String> markerList;

		ms = proj.getMarkerSet();
		allMarkers = ms.getMarkerNames();
		markerChrs = ms.getChrs();
		markerList = new Vector<String>();
		
		Centroids centroid = Centroids.load(centFile, false);
		float[][][] newCentroids = new float[allMarkers.length][][];
		
		int breakInd = 0;
		for (int i = 0; i < markerChrs.length; i++) {
			switch(markerChrs[i]) {
				case 23:
				case 24:
				case 25:
				case 26:
					if (breakInd == 0) breakInd = i;
					break;
				default:
					markerList.add(allMarkers[i]);
					break;
			}
		}
		
		for (int i = 0; i < breakInd; i++) {
			newCentroids[i] = new float[][]{{Float.NaN, Float.NaN}, {Float.NaN, Float.NaN}, {Float.NaN, Float.NaN}};
		}
		for (int i = breakInd; i < newCentroids.length; i++) {
			newCentroids[i] = centroid.getCentroids()[i - breakInd];
		}
		
		Centroids newCentObj = new Centroids(newCentroids, ms.getFingerprint());
		newCentObj.serialize(centFile + ".faked");
		String dir = proj.getProjectDir();
		String relativeFile = centFile.substring(dir.length());
		Centroids.exportToText(proj, relativeFile + ".faked", relativeFile + ".faked.txt");
		
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		Project proj;
		String filename = "lab.dat";
		String logfile = null;
		String centFile = null;
		
		
		String usage = "";
		
		if (numArgs == 0) {
			mockupGUI();
			return;
		}
		
		
		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				logfile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("cent=")) {
				centFile = args[i].split("=")[1];
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
			proj = new Project(filename, logfile, false);
			if (centFile != null) {
				fakeAutoCentFile(proj, centFile);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
}
