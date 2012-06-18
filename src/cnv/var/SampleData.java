package cnv.var;

import java.io.*;
import java.util.*;

import javax.swing.JOptionPane;

import common.*;
import cnv.filesys.*;

public class SampleData {
	public static final String[] BASIC_CLASSES = {"All", "Genotype"};
	public static final String[][][] KEYS_FOR_BASIC_CLASSES = {{{"0", "All"}}, {{"1", "A/A"}, {"2", "A/B"}, {"3", "B/B"}}};
	
//	public static final String[] BASIC_FILTERS = {"GC"};

	private String[] filters;
	private String[] covars;
	private String[] classes;
	private String[][][] classColorKeys;
	private String[] cnvClasses;
//	private Hashtable<String,String> sampleLookup;
//	private Hashtable<String,String> famIndLookup;
//	private Hashtable<String,String> indLookup;
	private Hashtable<String,String[]> lookup;
	private Hashtable<String,IndiPheno> sampleHash;
	private boolean failedToLoad;
	private int sexClassIndex;
	private String clusterUserSpecified;
	
	public SampleData(Project proj, boolean loadCNVs) {
		this(proj, loadCNVs?proj.getFilenames(Project.CNV_FILENAMES):null);
	}
	
	public SampleData(Project proj, String[] cnvFilesnames) {
		BufferedReader reader;
		String[] line, header;
		IntVector filterIs = new IntVector();
		IntVector covarIs = new IntVector();
		IntVector classIs = new IntVector();
		IntVector iv;
		DoubleVector dv;
		IndiPheno indi;
		int dnaIndex, famIndex, indIndex;
//		Hashtable<String,IndiPheno> sampleHash; 
		String filename;
		CountVector sexCountHash;
		String[] sexValues;
		int[] sexCounts;
		String[] ids;

		failedToLoad = true;
		if (cnvFilesnames == null) {
			cnvFilesnames = new String[0];
		}
		
		try {
			filename = proj.getFilename(Project.SAMPLE_DATA_FILENAME);
			if (!Files.exists(filename, proj.getJarStatus())) {
    			JOptionPane.showMessageDialog(null, "SampleData file does not exist: "+filename, "Error", JOptionPane.ERROR_MESSAGE);
				return;
			}			
			reader = Files.getReader(filename, proj.getJarStatus(), true, true); // to do, don't kill?
			header = reader.readLine().split("\t");
			dnaIndex = ext.indexOfStr("DNA", header);
			famIndex = ext.indexOfStr("FID", header);
			indIndex = ext.indexOfStr("IID", header);
			if (dnaIndex == -1) {
				System.err.println("Error - 'DNA' was not a header in the SampleData file; assuming lookup with first column");
				dnaIndex = 0;
			}
			if (cnvFilesnames.length > 0 && famIndex == -1) {
				System.err.println("Error - 'FID' was not a header in the SampleData file; cannot lookup cnv data");
				cnvFilesnames = new String[0];
			}
			if (cnvFilesnames.length > 0 && indIndex == -1) {
				System.err.println("Error - 'IID' was not a header in the SampleData file; cannot lookup cnv data");
				cnvFilesnames = new String[0];
			}
			for (int i = 1; i<header.length; i++) {
				if (header[i].toUpperCase().startsWith("FILTER=")) {
					filterIs.add(i);
				} else if (header[i].toUpperCase().startsWith("CLASS=")) {
					classIs.add(i);
				} else if (header[i].toUpperCase().startsWith("COVAR=")) {
					covarIs.add(i);
				}
			}
			filters = new String[filterIs.size()];
			for (int i = 0; i<filters.length; i++) {
				filters[i] = header[filterIs.elementAt(i)].split("=")[1];
			}
			covars = new String[covarIs.size()];
			for (int i = 0; i<covars.length; i++) {
				covars[i] = header[covarIs.elementAt(i)].split("=")[1];
			}
			classes = new String[classIs.size()];
			classColorKeys = new String[classIs.size()][][];
			for (int i = 0; i<classes.length; i++) {
				line = header[classIs.elementAt(i)].split(";");
				classes[i] = line[0].split("=")[1];
				classColorKeys[i] = new String[line.length-1][];
				for (int j = 1; j<line.length; j++) {
					classColorKeys[i][j-1] = line[j].split("=");
					if (classColorKeys[i][j-1].length != 2) {
						System.err.println("Error - invalid key for class '"+classes[i]+"'; must use format #=Key (not '"+line[j]+"'), separated by semicolons");
						classColorKeys[i][j-1] = new String[0];
					}
                }
			}
			sexClassIndex = ext.indexFactors(new String[][] {{"Sex", "CLASS=Sex", "Gender", "CLASS=Gender"}}, classes, false, true, true, false)[0];
			System.out.println(Array.toStr(classes));

			sexCountHash = new CountVector();
			sampleHash = new Hashtable<String,IndiPheno>();
			lookup = new Hashtable<String, String[]>();
			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				indi = new IndiPheno();
				
				ids = new String[] {line[dnaIndex], line[famIndex]+"\t"+line[indIndex], line[indIndex]};
				lookup.put(line[dnaIndex], ids);
				lookup.put(line[famIndex]+"\t"+line[indIndex], ids);
				lookup.put(line[indIndex], ids);
				
				dv = new DoubleVector();
				for (int i = 0; i<filterIs.size(); i++) {
					dv.add(Double.parseDouble(line[filterIs.elementAt(i)]));
				}
				indi.setFilters(dv.toArray());

				dv = new DoubleVector();
				for (int i = 0; i<covarIs.size(); i++) {
					dv.add(line[covarIs.elementAt(i)].equals(".")?Double.NaN:Double.parseDouble(line[covarIs.elementAt(i)]));
				}
				indi.setCovars(dv.toArray());

				iv = new IntVector();
				for (int i = 0; i<classIs.size(); i++) {
					iv.add(line[classIs.elementAt(i)].equals(".")||Integer.parseInt(line[classIs.elementAt(i)])<0?Integer.MIN_VALUE:Integer.parseInt(line[classIs.elementAt(i)]));
				}
				indi.setClasses(iv.toArray());
				if (sexClassIndex != -1) {
					sexCountHash.add(indi.getClasses()[sexClassIndex]+"");					
				}

				sampleHash.put(line[dnaIndex].toLowerCase(), indi);
			}
			reader.close();
			
			if (sexClassIndex != -1) {
				sexValues = sexCountHash.getValues();
				sexValues = Sort.putInOrder(sexValues, Sort.quicksort(sexValues, Sort.DESCENDING));
				if (!sexValues[0].equals("2")) {
					System.err.println("Error - warning no females listed in SampleData file; make sure 1=male and 2=female in the coding");
					JOptionPane.showMessageDialog(null, "descending "+ Array.toStr(sexValues, " ")+"\tError - warning no females listed in SampleData file; make sure 1=male and 2=female in the coding", "Error", JOptionPane.ERROR_MESSAGE);
				}
//
//				sexCountHash.sort(true);
//				sexValues = sexCountHash.getValues();
//				sexCounts = sexCountHash.getCounts();
//				if (!sexValues[0].equals("2")) {
//					System.err.println("Error - warning no females listed in SampleData file; make sure 1=male and 2=female in the coding");
//					JOptionPane.showMessageDialog(null, "ascending "+ Array.toStr(sexValues, " ")+"\tError - warning no females listed in SampleData file; make sure 1=male and 2=female in the coding", "Error", JOptionPane.ERROR_MESSAGE);
//				}
//			
			} else {
				JOptionPane.showMessageDialog(null, "Error - variable names 'Sex' was found in the SampleData file; also make sure 1=male and 2=female in the coding", "Error", JOptionPane.ERROR_MESSAGE);
			}
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+proj.getFilename(Project.SAMPLE_DATA_FILENAME)+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+proj.getFilename(Project.SAMPLE_DATA_FILENAME)+"\"");
			System.exit(2);
		}
		
		
		
		if (cnvFilesnames.length > 0) {
			loadCNVs(cnvFilesnames, proj.getJarStatus());
		} else {
			cnvClasses =  new String[0];
		}
		
		failedToLoad = false;
	}
	
	public boolean failedToLoad() {
		return failedToLoad;
	}
	
	public int getSexClassIndex() {
		return sexClassIndex;
	}
	
	public int getSexForIndividual(String id) {
		IndiPheno indi;
		String[] ids;
		
		indi = sampleHash.get(id);
//		indi = sampleHash.get("S_"+id);
//		sampleLookup.put("FI_"+line[famIndex]+"\t"+line[indIndex], "S_"+line[dnaIndex]);
//		sampleLookup.put("I_"+line[indIndex], "S_"+line[dnaIndex]);
//		famIndLookup.put("I_"+line[indIndex], "FI_"+line[famIndex]+"\t"+line[indIndex]);
//		famIndLookup.put("S_"+line[dnaIndex], "FI_"+line[famIndex]+"\t"+line[indIndex]);
//		indLookup.put("FI_"+line[famIndex]+"\t"+line[dnaIndex], "I_"+line[indIndex]);
//		indLookup.put("S_"+line[dnaIndex], "I_"+line[indIndex]);
		if (indi == null) {
			ids = lookup.get(id);
			if (ids != null) {
				indi = sampleHash.get(ids[0]);
			}
		}
		
		if (indi == null) {
			System.err.println("Error - id '"+id+"' was not present in the SampleData");
			return -1;
		}
		
		if (indi.getClasses()[sexClassIndex] == Integer.MIN_VALUE) {
			return -1;
		} else {
			return indi.getClasses()[sexClassIndex];
		}
	}
	
	public void loadCNVs(String[] files, boolean jar) {
		Vector<Hashtable<String,CNVariant[]>> finalHashes;
		CNVariantHash[] cnvhs;
		IndiPheno indi;
		String[] inds;
		String trav;
		long time;
		
		time = new Date().getTime();
		cnvClasses = new String[files.length];
		cnvhs = new CNVariantHash[files.length];
		for (int i = 0; i<files.length; i++) {
			cnvClasses[i] = ext.rootOf(files[i]);
			cnvhs[i] = CNVariantHash.load(files[i], CNVariantHash.CONSTRUCT_BY_IND, jar);
		}
		System.out.println("Read in CNV data in "+ext.getTimeElapsed(time));

		time = new Date().getTime();
		inds = HashVec.getKeys(sampleHash);
		for (int i = 0; i<inds.length; i++) {
			indi = sampleHash.get(inds[i]);
//			trav = lookup.get(inds[i]);
//			trav = famIndLookup.get("S_"+inds[i]);
			trav = lookup.get(inds[i])[1];

			finalHashes = new Vector<Hashtable<String,CNVariant[]>>();
			for (int j = 0; j<files.length; j++) {
				finalHashes.add(cnvhs[j].getDataFor(trav));
			}
			indi.setCNVclasses(finalHashes);
        }
		System.out.println("Added CNV data in "+ext.getTimeElapsed(time));
	}
	
	public String[] getFilters() {
		if (filters == null) {
			return new String[0];
		} else {
			return filters;
		}
	}

	public String[] getCovars() {
		if (covars == null) {
			return new String[0];
		} else {
			return covars;
		}
	}

	public String[] getClasses() {
		if (classes == null) {
			return new String[0];
		} else {
			return classes;
		}
	}

	public String[] getCnvClasses() {
		return cnvClasses;
	}

	
	public String[] lookup(String str) {
		return lookup.get(str);
	}
	
	public IndiPheno getIndiFromSampleHash(String sampleID) {
		return sampleHash.get(sampleID.toLowerCase());
	}
	
	public String[] getListOfSamples() {
		return HashVec.getKeys(sampleHash);
	}

	public IndiPheno getIndiPheno(String sample) {
		return sampleHash.get(sample);
	}
	
	public byte getClassForInd(String sample, int currentClass) {
		int[] classes;
		IndiPheno indi;
		
		indi = sampleHash.get(sample);
		if (indi == null) {
			return 0;
		}
		classes = indi.getClasses();
		if (currentClass == -1 || classes[currentClass] == Integer.MIN_VALUE) {
			return -1;
		} else {
			return (byte)classes[currentClass];
		}
	}
	
	public int getNumClasses() {
		return BASIC_CLASSES.length+classes.length+cnvClasses.length;
	}

	public int getNumActualClasses() {
		if (classes == null) {
			return -1;
		}
		return classes.length;
	}

	public int[] getClassCategoryAndIndex(int index) {
		int[] indices = new int[2];

		if (index<BASIC_CLASSES.length) {
			indices[0] = 0;
			indices[1] = index;
		} else if (index<BASIC_CLASSES.length+classes.length) {
			indices[0] = 1;
			indices[1] = index-BASIC_CLASSES.length;
		} else if (index<BASIC_CLASSES.length+classes.length+cnvClasses.length) {
			indices[0] = 2;
			indices[1] = index-BASIC_CLASSES.length-classes.length;
		} else {
			System.err.println("Error - invalid class index");
		}
		
		return indices;		
	}
	
	public String getClassName(int index) {
		int[] indices = getClassCategoryAndIndex(index);
		
		switch (indices[0]) {
		case 0:
			return BASIC_CLASSES[indices[1]];
		case 1:
			return classes[indices[1]];
		case 2:
			return cnvClasses[indices[1]];
		default:
			return null;
		}
	}	

	public String getActualClassName(int index) {
		return classes[index];
	}	

	public String[][] getActualClassColorKey(int index) {
		return classColorKeys[index];
	}	
}
