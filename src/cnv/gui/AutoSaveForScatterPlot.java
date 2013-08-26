package cnv.gui;

import java.io.File;

import javax.swing.JOptionPane;

import common.Files;

import cnv.filesys.AnnotationCollection;
import cnv.filesys.ClusterFilterCollection;

public class AutoSaveForScatterPlot implements Runnable {

	private ClusterFilterCollection clusterFilters;
	private String clusterFilterFilename;
	private AnnotationCollection annotations;
	private String annotationFilename;
	private int period;
	private boolean isKilled;
	private boolean isClusterFiltersUpdated;
	private boolean isAnnotationsUpdated;
	
	public AutoSaveForScatterPlot(ClusterFilterCollection clusterFilterCollection, String clusterFilterFilename, AnnotationCollection AnnotationCollection, String annotaionFilename, int periodInSeconds) {
		this.clusterFilters = clusterFilterCollection;
		this.clusterFilterFilename = clusterFilterFilename;
		this.annotations = AnnotationCollection;
		this.annotationFilename = annotaionFilename;
		this.period = periodInSeconds;
		this.isKilled = false;
		this.isClusterFiltersUpdated = false;
		this.isAnnotationsUpdated = false;
	}

//	public AutoSaveForScatterPlot(AnnotationCollection collection, String annotaionFilename, int periodInSeconds) {
//		this.clusterFilters = null;
//		this.clusterFilterFilename = null;
//		this.annotations = collection;
//		this.annotationFilename = annotaionFilename;
//		this.period = periodInSeconds;
//		this.isKilled = false;
//	}

	public void addToAutoSave(AnnotationCollection collection, String annotaionFilename) {
		this.annotations = collection;
		this.annotationFilename = annotaionFilename;
	}

	public void addToAutoSave(ClusterFilterCollection collection, String clusterFilterFilename) {
		this.clusterFilters = collection;
		this.clusterFilterFilename = clusterFilterFilename;
	}

	public boolean isClusterFilterNull() {
		return clusterFilters == null;
	}

	public boolean isAnnotationNull() {
		return annotations == null;
	}

	public void setClusterFilterUpdated(boolean isUpdated) {
		isClusterFiltersUpdated = isUpdated;
	}

	public void setAnnotationUpdated(boolean isUpdated) {
		isAnnotationsUpdated = isUpdated;
	}

	public void saveNow() {
//		if (clusterFilters != null) {
		if (isClusterFiltersUpdated) {
			clusterFilters.serialize(clusterFilterFilename);
			isClusterFiltersUpdated = false;
		}
		
//		if (annotations != null) {
		if (isAnnotationsUpdated) {
			Files.writeSerial(annotations, annotationFilename);
			isAnnotationsUpdated = false;
		}
	}
	
	public void kill() {
		if (clusterFilters != null) {
			new File(clusterFilterFilename).delete();
		}
		if (annotations != null) {
			new File(annotationFilename).delete();
		}
		
		isKilled = true;
	}
	
	public void run() {
		do {
			try {
				Thread.sleep(period * 1000);
			} catch (InterruptedException ie) {}
			
			saveNow();

//		} while (! isKilled && (Files.exists(clusterFilterFilename, false) || Files.exists(annotationFilename, false)));
//		} while (! isKilled && ((isClusterFiltersUpdated && !Files.exists(clusterFilterFilename, false)) || (isAnnotationsUpdated)));
		} while (! isKilled);
		
//		if (! isKilled ) {
//			JOptionPane.showMessageDialog(null, "Temporary file(s) containing new ClusterFilters and/or Annotations have been deleted; Auto-save has been turned off", "Error", JOptionPane.ERROR_MESSAGE);
//		}
	}
}
