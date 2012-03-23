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

	public byte getSize () {
		return (byte) (hash==null?0:hash.size());//???
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
