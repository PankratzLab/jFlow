package cnv.filesys;

import java.awt.BorderLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.*;
import java.util.*;

import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;

import common.Array;
import common.Files;
import common.HashVec;
import common.IntVector;
import common.ext;

import cnv.plots.GenericRectangle;

/**
 * This is data structure to hold a group of filters to screen data points with.
 * @author npankratz and zxu
 *
 */
public class ClusterFilterCollection implements Serializable {
	private static final long serialVersionUID = 1L;

	/**
	 * The structure of "hash" is Hashtable<String markerName, ArrayList<ClusterFilter> clusterFilter>.
	 */
	private Hashtable<String,ArrayList<ClusterFilter>> hash;
	
	// constructor, get(String markerName), load, serialize
	public ClusterFilterCollection() {
		hash = new Hashtable<String, ArrayList<ClusterFilter>>();
	}
	
	public void addClusterFilter(String markerName, ClusterFilter filter) {
//		ArrayList<ClusterFilter> filters;
//		
//		if (hash.containsKey(markerName)) {
//			filters = hash.get(markerName);
//		} else {
//			hash.put(markerName, filters = new ArrayList<ClusterFilter>());
//		}
//		filters.add(filter);

		if (!hash.containsKey(markerName)) {
			hash.put(markerName, new ArrayList<ClusterFilter>());
		}
		hash.get(markerName).add(filter);
	}

	public String[] getMarkerNames () {
		return HashVec.getKeys(hash);
	}

	public byte getSize (String markerName) {
		return (byte) (hash.get(markerName)==null?0:hash.get(markerName).size());//???
	}

	public int getSize () {
		return (hash==null?0:hash.size());
	}

	//??? How to select the last filter???
	public byte getGenotype(String markerName, byte index) {
		if (hash.containsKey(markerName) && hash.get(markerName).size() > index) {
			return hash.get(markerName).get(index).getNewGenotype();
		} else {
			System.err.println("Error - Trying to get a ClusterFilter that does not exist.");
			return (byte)-1;
		}
	}
		
	// filterIndividual, worry about this later!!

	public ArrayList<ClusterFilter> getClusterFilters(String markerName) {
		return hash.get(markerName);
	}

	public void deleteClusterFilter(String markerName, byte index) {
		if (hash.containsKey(markerName) && hash.get(markerName).size()>index) {
			hash.get(markerName).remove(index);
			if (hash.get(markerName).size()==0) {
				hash.remove(markerName);
			}
		} else {
			System.err.println("Error deleting the cluster filter: either no cluster filter associate with this marker name, or the index for the cluster filter to be deleted does not exist.");
		}
	}

	public void updateGenotype(String markerName, byte index, byte newGenotype) {
		hash.get(markerName).get(index).setNewGenotype(newGenotype);
	}
		
	public byte[] filterMarker(MarkerData markerData, float gcThreshold) {
		byte[] result, original;
		float[] realX;
		float[] realY;
		ArrayList<ClusterFilter> clusterFilters;
		
		original = markerData.getAB_Genotypes(); 
		clusterFilters = hash.get(markerData.getMarkerName());
		result = new byte[original.length];
		for (int j=0; j<original.length; j++) {
			if (markerData.getGCs()[j]<gcThreshold) {
				result[j]=(byte)-1;
			} else {
				result[j] = original[j];
			}
		}
		
		for (int i=0; clusterFilters != null && i < clusterFilters.size(); i++) {
			// get appropriate data (X/Y Theta/R LRR/BAF)
			switch(clusterFilters.get(i).getPlotType()) {
			case 0:
				realX = markerData.getX_Raws();
				realY = markerData.getY_Raws();
				break;
			case 1:
				realX = markerData.getXs();
				realY = markerData.getYs();
				break;
			case 2:
				realX = markerData.getThetas();
				realY = markerData.getRs();
				break;
			case 3:
				realX = markerData.getBAFs();
				realY = markerData.getLRRs();
				break;
			default:
				realX = markerData.getXs();
				realY = markerData.getYs();
			}
			// iterate through all samples
			for (int j=0; j<markerData.getAB_Genotypes().length; j++) {
				if (realX[j]>=clusterFilters.get(i).getXMin()
						&& realY[j]>=clusterFilters.get(i).getYMin()
						&& realX[j]<=clusterFilters.get(i).getXMax()
						&& realY[j]<=clusterFilters.get(i).getYMax()) {
					result[j]=clusterFilters.get(i).getNewGenotype();
					if (result[j] < -1) {
						System.err.println("Error - result["+j+"]="+result[j]);
					}
				}
			}
		}
		
		return result;
	}

	public GenericRectangle[] getRectangles(String markerName, byte plotType, byte thickness, boolean fill, boolean roundedCorners, byte color, byte layer) {
		ArrayList<ClusterFilter> clusterFilters;
		GenericRectangle[] result = new GenericRectangle[getSize(markerName)];
//		ArrayList<GenericRectangle> rectangles = new ArrayList<GenericRectangle>();
		
		clusterFilters = hash.get(markerName);
		for (int i = 0; clusterFilters != null && i < clusterFilters.size(); i++) {
			if (clusterFilters.get(i).getPlotType()==plotType) {
				result[i] = new GenericRectangle(clusterFilters.get(i).getXMin(),
												 clusterFilters.get(i).getYMin(),
												 clusterFilters.get(i).getXMax(),
												 clusterFilters.get(i).getYMax(),
												 thickness,
												 fill,
												 roundedCorners,
												 color,
												 layer);
			} else {
				result[i] = new GenericRectangle(clusterFilters.get(i).getXMin(),
												 clusterFilters.get(i).getYMin(),
												 clusterFilters.get(i).getXMax(),
												 clusterFilters.get(i).getYMax(),
												 thickness,
												 fill,
												 roundedCorners,
												 color,
												 (byte) -1);
			}
		}
		return result;
	}


	/*
	public GenericRectangle[] getRectangles(String markerName, byte plotType, byte thickness, boolean fill, boolean roundedCorners, byte color, byte layer) {
		ArrayList<ClusterFilter> clusterFilters;
//		GenericRectangle[] result = new GenericRectangle[getSize(markerName)];
		GenericRectangle[] result;
		ArrayList<GenericRectangle> rectangles = new ArrayList<GenericRectangle>();
		
		clusterFilters = hash.get(markerName);
		for (int i = 0; clusterFilters != null && i < clusterFilters.size(); i++) {
			if (clusterFilters.get(i).getPlotType()==plotType) {
//				result[i] = new GenericRectangle(clusterFilters.get(i).getXYBoundary()[0],
				rectangles.add(new GenericRectangle(clusterFilters.get(i).getXYBoundary()[0],
												 clusterFilters.get(i).getXYBoundary()[1],
												 clusterFilters.get(i).getXYBoundary()[2],
												 clusterFilters.get(i).getXYBoundary()[3],
												 thickness,
												 fill,
												 roundedCorners,
												 color,
												 layer)
						   );
			}
		}
		
		if (rectangles.size()>0) {
			result = new GenericRectangle[rectangles.size()];
			for (int i=0; i<rectangles.size(); i++) {
				result[i]=rectangles.get(i);
			}
		} else {
			result = null;
		}
		return result;
	}
	*/

	public void serialize(String filename) {
		Files.writeSerial(this, filename);
	}

	public static ClusterFilterCollection load(String filename, boolean jar) {
		return (ClusterFilterCollection)Files.readSerial(filename, jar, true);
	}
	
	public static String getClusterFilterFilenameSelection(Project proj) {
		/*
		return (String) JOptionPane.showInputDialog(null,
													"Please select a cluster filter file:",
													"Apply Cluster Filters",
													JOptionPane.QUESTION_MESSAGE,
													null,
													Files.list(proj.getDir(Project.DATA_DIRECTORY), null, proj.getProperty(Project.CLUSTER_FILTER_COLLECTION_FILENAME), false, proj.getJarStatus()),
													proj.getProperty(Project.CLUSTER_FILTER_COLLECTION_FILENAME));
		*/

		String result;
		result = (String)JOptionPane.showInputDialog(null,
													"Please select a cluster filter file:",
													"Apply Cluster Filters",
													JOptionPane.QUESTION_MESSAGE,
													null,
													Array.addStrToArray("(--Do not apply any cluster filter--)", Files.list(proj.getDir(Project.DATA_DIRECTORY), null, proj.getProperty(Project.CLUSTER_FILTER_COLLECTION_FILENAME), false, proj.getJarStatus())),
													proj.getProperty(Project.CLUSTER_FILTER_COLLECTION_FILENAME));
		if (result==null) {
			result = "cancel";
		} else if (result.equals("(--Do not apply any cluster filter--)")) {
			result = null;
		}
		return result;
	}

	public static String getGenotypeLookupTableSelection(Project proj) {
		/*
		return (String) JOptionPane.showInputDialog(null,
													"Please select a cluster filter file:",
													"Apply Cluster Filters",
													JOptionPane.QUESTION_MESSAGE,
													null,
													Files.list(proj.getDir(Project.DATA_DIRECTORY), null, proj.getProperty(Project.CLUSTER_FILTER_COLLECTION_FILENAME), false, proj.getJarStatus()),
													proj.getProperty(Project.CLUSTER_FILTER_COLLECTION_FILENAME));
		*/

		String result;
		result = (String)JOptionPane.showInputDialog(null,
													"Please select a AB genotyp lookup table:",
													"Select AB Genotype Lookup Table",
													JOptionPane.QUESTION_MESSAGE,
													null,
													new String[] {"Lookup Table 1", "Lookup Table 2", "Lookup Table 3"},
													proj.getProperty(Project.CLUSTER_FILTER_COLLECTION_FILENAME));
		return result;
	}

	public static String getClusterFilterFilenameSelection_Old(Project proj) {
		JFrame frame;
		JLabel label;
		JComboBox filenameComboBox;
		Vector<String> options;
		String[] filenames;
//		String filenameSelected;
		byte index;
		FilenameSelectionListener filenameSelectionListener;
		
		frame = new JFrame("Alert");
		frame.setSize(500, 150);
		frame.setVisible(true);
		frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		
//		create frame with question as label: Which cluster filter file would you like to apply (if any) to the genotypes before exporting?
		label = new JLabel("Which cluster filter file would you like to\n apply (if any) to the genotypes before exporting?");
		frame.add(label, BorderLayout.NORTH);

		// get all files with "clusterFilter.ser" suffix in data directory and add to options
		filenames = Files.list(proj.getDir(Project.DATA_DIRECTORY), null, "clusterFilters.ser", false, false);// 1) Property; 2) jar;
		options = new Vector<String> (Arrays.asList(filenames));
		// include None as an option
        options.add(0, "-");

		// include options in comboBox
		filenameComboBox = new JComboBox(options);

		// if default cluster file filter file is not included in options and it exists, then include as first option
//		System.out.println(proj.getFilename(Project.CLUSTER_FILTER_COLLECTION_FILENAME, Project.DATA_DIRECTORY, false, false));
//        if (options.contains(proj.getFilename(Project.CLUSTER_FILTER_COLLECTION_FILENAME, Project.DATA_DIRECTORY, false, false))) {
//        	comboBox.setSelectedIndex(options.indexOf(proj.getFilename(Project.CLUSTER_FILTER_COLLECTION_FILENAME, Project.DATA_DIRECTORY, false, false)));
        if (options.contains("clusterFilters.ser")) {
        	filenameComboBox.setSelectedIndex(options.indexOf("clusterFilters.ser"));
        } else {
        	filenameComboBox.setSelectedIndex(0);
        }

		// determine answer, should be null if None
        /*
        ActionListener filenameComboBoxListener = new ActionListener() {
    		public void actionPerformed(ActionEvent e) {
    			String filenameSelected = (String)filenameComboBox.getSelectedItem();
    		}
    		public void comboBoxChanged(ActionEvent e) {
    		}
    	};
		filenameComboBox.addActionListener(this);
		frame.add(filenameComboBox, BorderLayout.SOUTH);
		index=(byte) options.indexOf(listenerResult);
		*/

		filenameSelectionListener = new FilenameSelectionListener();
        filenameComboBox.addActionListener(filenameSelectionListener);
		frame.add(filenameComboBox, BorderLayout.SOUTH);
		index=(byte) options.indexOf(filenameSelectionListener.getResult());
		while (filenameSelectionListener.getResult()==null) {
			index=(byte) options.indexOf(filenameSelectionListener.getResult());
		}

		return index<=0?null:options.get(index);
	}
	
	public static void merge(String[] filenames, String outfile) {
		String[] markerNames;
		ClusterFilterCollection master, trav;
		Vector<String> v = new Vector<String>();
		ArrayList<ClusterFilter> masterArray, travArray;
		
		master = new ClusterFilterCollection();
		for (int i = 0; i < filenames.length; i++) {
			trav = load(filenames[i], false);
			markerNames = trav.getMarkerNames();
			for (int j = 0; j < markerNames.length; j++) {
				masterArray = master.getClusterFilters(markerNames[j]);
				travArray = trav.getClusterFilters(markerNames[j]);
				if (masterArray != null) {
					HashVec.addIfAbsent(markerNames[j], v);
				}
				for (int k = 0; k < travArray.size(); k++) {
					master.addClusterFilter(markerNames[j], travArray.get(k));
				}
			}
		}
		
		master.serialize(outfile);
		
		if (v.size() > 0) {
			System.out.println("The following markers had cluster filters in multiple files:");
			System.out.println(Array.toStr(Array.toStringArray(v), "\n"));
		}
			
	}
	
	public static void describe(String filename) {
		PrintWriter writer;
		String[] markerNames;
		ClusterFilterCollection trav;
		int count;
		
		trav = load(filename, false);

		count = 0;
		try {
			writer = new PrintWriter(new FileWriter(ext.rootOf(filename)+"_described.xln"));
			writer.println("MarkerName\t#filters");
			markerNames = trav.getMarkerNames();
			for (int i = 0; i < markerNames.length; i++) {
				writer.println(markerNames[i]+"\t"+trav.getClusterFilters(markerNames[i]).size());
				count += trav.getClusterFilters(markerNames[i]).size();
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + ext.rootOf(filename)+"_described.xln");
			e.printStackTrace();
		}
		
		System.out.println("There were a total of "+count+" clusterFilters across "+trav.getSize()+" markers");
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String[] filenames = null;
		String out = "merged.ser";
		String filename = "clusterFilters.ser";

		String usage = "\n" + 
				"cnv.filesys.ClusterFilterCollection requires 0-1 arguments\n" +
				"   (1) names of clusterFilter files to merge (i.e. files=file1.ser,file2.ser,file3.ser (not the default))\n" + 
				"   (2) output filename (i.e. out="+out+" (default))\n" +
				" OR:\n" +
				"   (1) names of clusterFilter file to describe (i.e. files="+filename+" (default))\n" + 
				"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("files=")) {
				filenames = args[i].split("=")[1].split(",");
				numArgs--;
			} else if (args[i].startsWith("out=")) {
				out = args[i].split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
//		filenames = new String[] {"clusterFilters_part1.ser", "clusterFilters_part2.ser"};
		try {
			if (filenames != null) {
				merge(filenames, out);
			} else {
				describe(filename);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
}

class FilenameSelectionListener implements ActionListener {
	private String result;
	
	public FilenameSelectionListener() {
		result = null;
	}

	public void actionPerformed(ActionEvent e) {
		JComboBox comboBox = (JComboBox)e.getSource();
		result = (String)comboBox.getSelectedItem();
//		System.out.println("Result in ClusterFilterCollection.actionPerformed(): "+result);
	}
	
	public void comboBoxChanged(ActionEvent e) {
//		JComboBox comboBox = (JComboBox)e.getSource();
//		result = (String)comboBox.getSelectedItem();
	}

	public String getResult() {
//		System.out.println("Result in ClusterFilterCollection.getValue(): "+result);
		return result;
	}
}
