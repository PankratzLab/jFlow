package cnv.filesys;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.Serializable;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Vector;

import javax.swing.JOptionPane;

import common.Array;
import common.Files;
import common.HashVec;
import common.Logger;
import common.Sort;
import common.ext;

public class AnnotationCollection implements Serializable {
	private static final long serialVersionUID = 1L;
	
	private Hashtable<Character, String> commentsHash;		//annotation keys
	private Hashtable<String, Vector<String>> markerAnnotations;	//annotations organized by markers.   markerName -> AnnotationA, AnnotationB, ...
	private Hashtable<String, Vector<String>> annotationMarkerLists;	//annotations organized by annotation keys.   Annotation -> marker1Name, marker2Name, ...
	
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
//		if (response != 1) {
		if (response == 0) {
			commentsHash.remove(c);
			markers = annotationMarkerLists.get(c+"");
			for (int i=0; markers != null && i < markers.size(); i++) {
				markerAnnotations.remove(markers.elementAt(i));
			}
			System.out.println(ext.listWithCommas(HashVec.getKeys(annotationMarkerLists, true, false)));
			annotationMarkerLists.remove(c+"");
			System.out.println(ext.listWithCommas(HashVec.getKeys(annotationMarkerLists, true, false)));
			System.out.println();
		}
		
	}

	public void renameAnnotation(char c, String newComment) {
		commentsHash.put(c, newComment);
	}
	
	public void renameAnnotationKey(char oldShortcut, char newShortcut) {
		//TODO
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
			if (markerAnnotations.get(markerName).size() == 0) {
				markerAnnotations.remove(markerName);
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

	public boolean markerHasAnyAnnotation(String markerName) {
		return markerAnnotations.containsKey(markerName);
	}
	
	public void dumpLists(Project proj) {
		String[] list, keys;
		Vector<String> annotationsVector;
		String annotationsOfTheMarker;

		list = new String[markerAnnotations.size()];
		keys = HashVec.getKeys(markerAnnotations);
		for (int i = 0; i < list.length; i++) {
//			list[i] = keys[i]+"\t"+Array.toStr(markerAnnotations.get(keys[i]));
			annotationsVector = markerAnnotations.get(keys[i]);
			annotationsOfTheMarker = "";
			for (int j=0; j<annotationsVector.size(); j++) {
				annotationsOfTheMarker += (j==0?"":";") + commentsHash.get(annotationsVector.elementAt(j).toCharArray()[0]);
			}
			list[i] = keys[i] + "\t" + annotationsOfTheMarker;
		}
		Files.writeList(list, proj.getProjectDir() + "annoations.xln");
//		Files.writeList(list, proj.getDir(Project.ANNOTATION_DIRECTORY) + "annoations.xln");

		keys = HashVec.getKeys(annotationMarkerLists);
		for (int i = 0; i < keys.length; i++) {
			Files.writeList(Array.toStringArray(annotationMarkerLists.get(keys[i])), proj.getProjectDir() + "annoation_" + keys[i] + "_" + commentsHash.get(keys[i]) + ".xln");
		}
	}
	
	public static AnnotationCollection loadFromLists(String filename, Logger log) {
		AnnotationCollection annotationCollection;
		annotationCollection = new AnnotationCollection();
		appendFromLists(filename, annotationCollection, null, log);
		return annotationCollection;
	}
	
	public static void appendFromLists(String filename, AnnotationCollection currentAnnotationCollection, String[] markerNamesProj, Logger log) {
		BufferedReader reader;
		String[] line;
		char key;
		char[] keys;
		Hashtable<String, Character> annotationKeys;
		boolean found;

		annotationKeys = new Hashtable<String, Character>();
		keys = currentAnnotationCollection.getKeys();
		for (int i = 0; i < keys.length; i++) {
			annotationKeys.put(currentAnnotationCollection.getDescriptionForComment(keys[i], false, false), keys[i]);
		}
		try {
			reader = new BufferedReader(new FileReader(filename));
			while (reader.ready()) {
				key = 0;
				line = reader.readLine().split("[\t;]");
				found = false;
				if (markerNamesProj != null) {
					for (int i = 0; i < markerNamesProj.length; i++) {
						if (markerNamesProj[i].equals(line[0])) {
							found = true;
						}
					}
				}
				if (markerNamesProj != null && ! found) {
					log.reportError("Skipped importing annotations for marker " + line[0] + ", which is found in the annotation file " + filename + ", but not in the marker name list of the project.");
				} else {
					for (int i = 1; i < line.length; i++) {
						if (annotationKeys.containsKey(line[i])) {
							key = annotationKeys.get(line[i]);
						} else {
//							letters = line[i].toLowerCase().toCharArray();
//							for (int j = 0; j < letters.length; j++) {
//								if (letters[j] >= 97 && letters[j] <= 122 && ! currentAnnotationCollection.containsKey(letters[j])) {
//									key = letters[j];
//									break;
//								}
//							}
//							if (key == 0) {
//								for (int j = 97; j <= 122; j++) {
//									if (! currentAnnotationCollection.containsKey((char) j)) {
//										key = (char) j;
//										break;
//									}
//								}
//							}
							key = assignKey(line[i], currentAnnotationCollection);
							if (key == 0) {
								if (log == null) {
									System.out.println("cannot automatically assign a shortcut to the annotation key '" + line[i] + "'. Skipped importing this annotation for marker " + line[0] + ".");
								} else {
									log.reportError("cannot automatically assign a shortcut to the annotation key '" + line[i] + "'. Skipped importing this annotation for marker " + line[0] + ".");
								}
							} else {
								annotationKeys.put(line[i], key);
								currentAnnotationCollection.addAnnotation(key, line[i]);
							}
						}
						currentAnnotationCollection.addAnnotationForMarker(line[0], key);
					}
				}
			}
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}


	public static char assignKey(String comment, AnnotationCollection currentAnnotationCollection) {
		char key;
		char[] letters;

		key = 0;
		letters = comment.toLowerCase().toCharArray();
		for (int j = 0; j < letters.length; j++) {
			if (letters[j] >= 97 && letters[j] <= 122 && ! currentAnnotationCollection.containsKey(letters[j])) {
				key = letters[j];
				break;
			}
		}
		if (key == 0) {
			for (int j = 97; j <= 122; j++) {
				if (! currentAnnotationCollection.containsKey((char) j)) {
					key = (char) j;
					break;
				}
			}
		}

		return key;
	}


	public char[] getKeys() {
//		String[] result;
//		Enumeration<Character> keys;
//		result = new String[commentsHash.size()];
//		keys = commentsHash.keys();
//		for (int i=0; i<result.length; i++) {
//			result[i] = keys.nextElement()+"";
//		}
//		return Array.toCharArray(Sort.putInOrder(result));

		String[] commentDescriptions;
		char[] result;
		Enumeration<Character> keys;
		char key;
		Hashtable<String, Character> temp;

		temp = new Hashtable<String, Character>(commentsHash.size() + 1, (float) 1.00);
		keys = commentsHash.keys();
		commentDescriptions = new String[commentsHash.size()];
		for (int i=0; i<commentDescriptions.length; i++) {
			key = keys.nextElement();
			commentDescriptions[i] = commentsHash.get(key);
			temp.put(commentDescriptions[i], key);
		}
		commentDescriptions = Sort.putInOrder(commentDescriptions);

		result = new char[commentsHash.size()];
		for (int i=0; i<commentDescriptions.length; i++) {
			result[i] = temp.get(commentDescriptions[i]);
		}
		return result;
	}
	
	public boolean containsKey(char c) {
		return commentsHash.containsKey(c);
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

	public String[] getMarkerLists(char key) {
		return annotationMarkerLists.get(key + "").toArray(new String[0]);
	}

	public String getDescriptionForComment(char c, boolean includeShortcuts, boolean includeNumbers) {
		return (includeShortcuts? "'" + c + "' " : "") + commentsHash.get(c) + (includeNumbers? " (n=" + annotationMarkerLists.get(c+"").size() + ")" : "");
	}
	
	public void serialize(String filename) {
		Files.writeSerial(this, filename);
	}

	public static AnnotationCollection load(String filename, boolean jar) {
		return (AnnotationCollection)Files.readSerial(filename, jar, true);
	}
}
