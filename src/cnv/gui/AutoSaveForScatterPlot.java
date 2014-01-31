package cnv.gui;

import java.io.File;
import cnv.filesys.AnnotationCollection;
import cnv.filesys.ClusterFilterCollection;

public class AutoSaveForScatterPlot implements Runnable {

	private ClusterFilterCollection clusterFilters;
	private String clusterFilterFilename;
	private AnnotationCollection annotations;
	private String annotationFilename;
	private int period;
	private boolean killed;
	private boolean isClusterFiltersUpdated;
	private boolean isAnnotationsUpdated;
	
	public AutoSaveForScatterPlot(ClusterFilterCollection clusterFilterCollection, String clusterFilterFilename, AnnotationCollection AnnotationCollection, String annotaionFilename, int periodInSeconds) {
		this.clusterFilters = clusterFilterCollection;
		this.clusterFilterFilename = clusterFilterFilename;
		this.annotations = AnnotationCollection;
		this.annotationFilename = annotaionFilename;
		this.period = periodInSeconds;
		this.killed = false;
		this.isClusterFiltersUpdated = false;
		this.isAnnotationsUpdated = false;
	}

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
		if (isClusterFiltersUpdated) {
			clusterFilters.serialize(clusterFilterFilename);
			isClusterFiltersUpdated = false;
		}
		
		if (isAnnotationsUpdated) {
			annotations.serialize(annotationFilename);
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
		
		killed = true;
	}
	
	public void run() {
		while (!killed) {
			saveNow();
			try {
				Thread.sleep(period * 1000);
			} catch (InterruptedException ie) {}
		}
	}
}
