package cnv.var;

import java.io.*;
import java.util.*;
import common.*;

public class CNVariantHash implements Serializable {
	public static final long serialVersionUID = 1L;
	
	public static final int CONSTRUCT_BY_IND = 1;
	public static final int CONSTRUCT_ALL = 2;

	private Hashtable<String,Hashtable<String,CNVariant[]>> hashes;

	public CNVariantHash(String filename, int structureType, boolean jar) {
		Hashtable<String,Hashtable<String,Vector<CNVariant>>> vHashes;
		Hashtable<String,Vector<CNVariant>> vHash;
		Vector<CNVariant> v;
		Hashtable<String,CNVariant[]> finalHash;
		CNVariant cnv;
		String trav, temp;
		String[] inds, chrs;
//		long time;
		BufferedReader reader;
		String[] line;
		ProgressBarDialog prog;
		int count;

		prog = new ProgressBarDialog("Converting CNVs", 0, Files.getSize(filename, jar), 800, 200);
		
//		time = new Date().getTime();
		vHashes = new Hashtable<String, Hashtable<String,Vector<CNVariant>>>();
		count = 0;
		try {
			reader = Files.getReader(filename, jar, true, true);

			reader.mark(1000);
			temp = reader.readLine();
			count += temp.length();
			line = temp.trim().split("[\\s]+");
			if (!line[2].toLowerCase().equals("chr")&&Positions.chromosomeNumber(line[2])!=-1) {
				reader.reset();
			}
			while (reader.ready()) {
				temp = reader.readLine();
				count += temp.length();
				prog.setProgress(count);
				cnv = new CNVariant(temp.trim().split("[\\s]+"));
				trav = cnv.getFamilyID()+"\t"+cnv.getIndividualID();
				if (structureType == CONSTRUCT_ALL) {
					cnv.setFamilyID(null);
					cnv.setIndividualID(null);
					trav = "all";					
				}

				if (vHashes.containsKey(trav)) {
					vHash = vHashes.get(trav);
				} else {
					vHashes.put(trav, vHash = new Hashtable<String,Vector<CNVariant>>());
				}
				if (vHash.containsKey(cnv.getChr()+"")) {
					v = vHash.get(cnv.getChr()+"");
				} else {
					vHash.put(cnv.getChr()+"", v = new Vector<CNVariant>());
				}
				v.add(cnv);
				
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+filename+"\" not found in current directory");
			fnfe.printStackTrace();
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+filename+"\"");
			ioe.printStackTrace();
		}
		prog.close();

		hashes = new Hashtable<String,Hashtable<String,CNVariant[]>>();
//		time = new Date().getTime();
		
		inds = HashVec.getKeys(vHashes);
		prog = new ProgressBarDialog("Serializing CNVs", 0, inds.length, 800, 200);
		for (int i = 0; i<inds.length; i++) {
			prog.setProgress(i);
			vHash = vHashes.get(inds[i]);
			finalHash = new Hashtable<String,CNVariant[]>();
			chrs = HashVec.getKeys(vHash);
			for (int j = 0; j<chrs.length; j++) {
				finalHash.put(chrs[j], CNVariant.toArray(vHash.get(chrs[j])));
            }
			hashes.put(inds[i], finalHash);
        }
		prog.close();
	}
	
	public Hashtable<String,CNVariant[]> getDataFor(String famID_indID) {
		if (hashes.containsKey(famID_indID)) {
			return hashes.get(famID_indID);
		} else {
			return new Hashtable<String,CNVariant[]>(); 
		}
	}
	
	public void serialize(String filename) {
		Files.writeSerial(this, filename);
	}

	public static CNVariantHash load(String filename, int structureType, boolean jar) {
		CNVariantHash hashes;
		String suffix;
		
		if (structureType == CONSTRUCT_BY_IND) {
			suffix = ".ind.ser";
		} else if (structureType == CONSTRUCT_ALL) {
			suffix = ".all.ser";
		} else {
			System.err.println("Error - invalid CONSTRUCT type");
			suffix = null;
		}
		
		if (Files.exists(filename+suffix, jar)) {
			hashes = (CNVariantHash)Files.readSerial(filename+suffix, jar, false);
		} else {
			hashes = new CNVariantHash(filename, structureType, jar);
			hashes.serialize(filename+suffix);
		}
		
		return hashes;
	}

	public static void main(String[] args) {
	    int numArgs = args.length;
	    String filename = "C:\\Documents and Settings\\npankrat\\My Documents\\CNV_PD\\data\\penncnv_1SNP.cnv";

	    String usage = "\n"+"cnv.var.CNVariantHash requires 0-1 arguments\n"+"   (1) filename (i.e. file="+filename+" (default))\n"+"";

	    for (int i = 0; i<args.length; i++) {
		    if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
			    System.err.println(usage);
			    System.exit(1);
		    } else if (args[i].startsWith("file=")) {
			    filename = args[i].split("=")[1];
			    numArgs--;
		    }
	    }
	    if (numArgs!=0) {
		    System.err.println(usage);
		    System.exit(1);
	    }
	    try {
		    new CNVariantHash(filename, 1, false);
	    } catch (Exception e) {
		    e.printStackTrace();
	    }
    }
}
