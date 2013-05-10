package cnv.filesys;

import java.io.Serializable;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Vector;

import javax.swing.JOptionPane;

import common.Array;
import common.Files;
import common.HashVec;

public class AnnotationCollection implements Serializable {
	private static final long serialVersionUID = 1L;
	
	private Hashtable<Character, String> commentsHash;
	private Hashtable<String, Vector<String>> markerAnnotations;	//markerName -> AnnotationA, AnnotationB, ...
	private Hashtable<String, Vector<String>> annotationMarkerLists;	//Annotation -> marker1Name, marker2Name, ...  
	
	public AnnotationCollection() {
		commentsHash = new Hashtable<Character, String>();
		markerAnnotations = new Hashtable<String, Vector<String>>();
		annotationMarkerLists = new Hashtable<String, Vector<String>>();
	}

	public void addAnnotation(char c, String comment) {
		commentsHash.put(c, comment);
		annotationMarkerLists.put(c + "", new Vector<String>());
	}

	public void removeAnnotation(char c) {
		int response;
		Vector<String> markers;
		
		response = JOptionPane.showConfirmDialog(null, "This will remove the annotaion '" + commentsHash.get(c) + "' from all markers (n="+annotationMarkerLists.get(c+"").size() + ") from the annoation database", "Warning", JOptionPane.ERROR_MESSAGE);
		if (response != 1) {
			commentsHash.remove(c);
			markers = annotationMarkerLists.get(c);
			for (int i=0; markers != null && i < markers.size(); i++) {
				markerAnnotations.remove(markers.elementAt(i));
			}
			annotationMarkerLists.remove(c);
		}
		
	}

	public void renameAnnotation(char c, String newComment) {
		
	}
	
	public void addAnnotationForMarker(String markerName, char c) {
		HashVec.addToHashVec(markerAnnotations, markerName, c+"", true);
		HashVec.addToHashVec(annotationMarkerLists, c+"", markerName, true);
	}
	
	public void removeAnnotationForMarker(String markerName, char c) {
		if (markerAnnotations.containsKey(markerName)) {
			if (markerAnnotations.get(markerName).contains(c+"")) {
				markerAnnotations.get(markerName).remove(c+"");
				annotationMarkerLists.get(c+"").remove(markerName);
			} else {
				System.err.println("Error - cannot remove "+c+" from "+markerName+" if it's not already checked");
			}
		} else {
			System.err.println("Error - cannot remove "+c+" from "+markerName+" if the marker does not have any annotation");
		}
	}
	
	public boolean markerHasAnnotation(String markerName, char annotation) {
		if (markerAnnotations.containsKey(markerName)) {
			return markerAnnotations.get(markerName).contains(annotation + "");
		} else {
			return false;
		}
	}

	public char[] annotationsForMarker(String markerName) {
		if (markerAnnotations.containsKey(markerName)) {
			return Array.toCharArray(Array.toStringArray(markerAnnotations.get(markerName)));
		} else {
			return new char[0];
		}
	}

	public void dumpLists(Project proj) {
		String[] list, keys;

		list = new String[markerAnnotations.size()];
		keys = HashVec.getKeys(markerAnnotations);
		for (int i = 0; i < list.length; i++) {
			list[i] = keys[i]+"\t"+Array.toStr(markerAnnotations.get(keys[i]));
		}
		Files.writeList(list, proj.getProjectDir()+"markerAnnoations.xln");

		keys = HashVec.getKeys(annotationMarkerLists);
		for (int i = 0; i < keys.length; i++) {
			Files.writeList(Array.toStringArray(annotationMarkerLists.get(keys[i])), proj.getProjectDir()+"annoation_"+keys[i]+"_"+commentsHash.get(keys[i])+".xln");
		}
	}
	
	public char[] getKeys() {
//		return Array.toCharArray(HashVec.getKeys(commentsHash));

		char[] result;
		Enumeration<Character> keys;
		result = new char[commentsHash.size()];
		keys = commentsHash.keys();
		for (int i=0; i<result.length; i++) {
			result[i] = keys.nextElement();
		}
		return result;
	}
	
	public String[] getMarkerLists() {
		String[] result;
		Enumeration<String> keys;

		result = new String[markerAnnotations.size()];
		keys = markerAnnotations.keys();
		for (int i=0; i<result.length; i++) {
			result[i] = keys.nextElement();
		}
		return result;
	}

	public String getDescriptionForComment(char c, boolean shortcuts) {
		return (shortcuts?"'"+c+"' ":"") + commentsHash.get(c) + " (n=" + annotationMarkerLists.get(c+"").size() + ")";
	}
}
