// -Xms1024M -Xmx1024M
package cnv.analysis;

import java.io.*;
import java.util.*;

import javax.swing.JOptionPane;

import cnv.filesys.*;
import cnv.var.CNVariant;
import cnv.var.IndiPheno;
import cnv.var.SampleData;
import common.*;
import filesys.Segment;

public class Mosaicism {
	public static final String[] HEADER = {"Sample", "Arm", "#CNVs", "Summed_Size", "%covered", "Custom_metric", "LRR_SD", "LRR_SD_flag", "Flagged_via_Mosaicism", "Mosaicism_level", "Mosaicism description"};
	public static final double LOWER_BOUND = 0.15;
	public static final double UPPER_BOUND = 0.85;

	public static void findOutliers(Project proj) {
		PrintWriter writer;
		String[] samples;
		Sample samp;

		int chr;
		long time;
		Hashtable<String,String> hash = new Hashtable<String,String>();
		String[] markerNames;
		byte[] chrs;
		int[] positions;
		boolean[] snpDropped;
		int[][] chrBoundaries;
		float baf;
		FloatVector lrrVec, bafVec;
		float[] lrrs, bafs;
		MarkerSet markerSet;

		hash = proj.getFilteredHash();

		chrBoundaries = new int[27][3];
		snpDropped = null;
		for (int i = 0; i<chrBoundaries.length; i++) {
			chrBoundaries[i][0] = chrBoundaries[i][1] = chrBoundaries[i][2] = -1;
		}
		time = new Date().getTime();
		markerSet = proj.getMarkerSet();
		markerNames = markerSet.getMarkerNames();
		chrs = markerSet.getChrs();
		positions = markerSet.getPositions();
		snpDropped = new boolean[markerNames.length];
		System.out.println("Mosacism will be estimated using "+markerNames.length+" markers");
		chr = 0;
		for (int i = 0; i<markerNames.length; i++) {
			snpDropped[i] = hash.containsKey(markerNames[i]);
			if (positions[i]>Positions.CENTROMERE_MIDPOINTS[chr]&&chrBoundaries[chr][1]==-1) {
				chrBoundaries[chr][1] = i;
			}
			if (chrs[i]>chr||i==markerNames.length-1) {
				if (chr!=0) {
					chrBoundaries[chr][2] = i-1;
				}
				chr = chrs[i];
				chrBoundaries[chr][0] = i;
			}
		}
		chrBoundaries[0][0] = 0;
		chrBoundaries[0][2] = markerNames.length-1;
		
		for (int i = 0; i<chrBoundaries.length; i++) {
			if (chrBoundaries[i][0] == -1 || chrBoundaries[i][2] == -1) {
				System.err.println("Error - no data for chromosome '"+i+"'");
			}
		}
		

		samples = proj.getSamples();
		try {
			writer = new PrintWriter(new FileWriter(proj.getDir(Project.RESULTS_DIRECTORY, true)+"Mosaicism.xln"));
			writer.println("Sample\tBand\tLRR N\tmean LRR\tBAF N\tSD of BAF (0.15-0.85)\tIQR of BAF (0.15-0.85)\t%Homo");
			for (int i = 0; i<samples.length; i++) {
				System.out.println((i+1)+" of "+samples.length+" in "+ext.getTimeElapsed(time));
				time = new Date().getTime();
				samp = proj.getSample(samples[i]);
				if (samp.getFingerprint()!=markerSet.getFingerprint()) {
					System.err.println("Error - cannot estimate mosaics if MarkerSet and Sample ("+samples[i]+") don't use the same markers");
					return;
				}
				lrrs = samp.getLRRs();
				bafs = samp.getBAFs();
				for (int j = 1; j<=23; j++) {
					for (int arm = 0; arm<2; arm++) {
						lrrVec = new FloatVector();
						bafVec = new FloatVector();
						for (int k = (arm==0?chrBoundaries[j][0]:chrBoundaries[j][1]); k<(arm==0?chrBoundaries[j][1]:chrBoundaries[j][2]+1); k++) {
							if (!snpDropped[k]) {
								if (!Float.isNaN(lrrs[k])) {
									lrrVec.add(lrrs[k]);
								}
								baf = bafs[k];
								if (baf>LOWER_BOUND&&baf<UPPER_BOUND) {
									bafVec.add(baf);
								}
							}
						}
						if (lrrVec.size()>100) {
							writer.println(samples[i]+"\t"+"chr"+j+(arm==0?"p":"q")+"\t"+lrrVec.size()+"\t"+ext.formDeci(Array.mean(lrrVec.toArray()), 5)+"\t"+bafVec.size()+(bafVec.size()>10?"\t"+ext.formDeci(Array.stdev(bafVec.toArray(), true), 5)+"\t"+ext.formDeci(Array.iqr(bafVec.toArray()), 5):"\t.\t.")+"\t"+ext.formDeci((double)(lrrVec.size()-bafVec.size())/(double)lrrVec.size(), 5));
						}
					}
				}
				writer.flush();
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to summary file");
			e.printStackTrace();
		}
	}

	public static void checkForOverlap(Project proj) {
		BufferedReader reader;
        PrintWriter writer;
        String[] line;
        Vector<String> v;
        Vector<String[]> list;
        int count;
        String listOfMosaicArms;
        SampleData sampleData;
        int chr, sum;
        Segment arm;
        IndiPheno indiPheno;
        CNVariant[] cnvs;
        String[] sampleList;
        String[][] listOfArms;
        Hashtable<String,String> lrrsdHash, mosaicHash;
        double proportion;
        long time;
        String[] cnvFiles;
        
        time = new Date().getTime();
        listOfMosaicArms = proj.getFilename(Project.MOSAIC_ARMS_FILENAME, false, false);
        cnvFiles = proj.getFilenames(Project.CNV_FILENAMES);
        if (cnvFiles.length == 0) {
        	System.err.println("Error - need to specify the name of a CNV file in the project properties file before running Mosaicism.checkForOverlap()");
        	return;
        }
        sampleData = new SampleData(proj, new String[] {cnvFiles[0]});
        if (Files.exists(proj.getProjectDir()+"lrr_sd.xln", proj.getJarStatus())) {
        	lrrsdHash = HashVec.loadFileToHashString(proj.getProjectDir()+"lrr_sd.xln", false);
        } else {
        	System.err.println("Warning - could not find 'lrr_sd.xln' in project directory; no flags will be generated");
        	lrrsdHash = new Hashtable<String,String>();
        }
        if (Files.exists(proj.getFilename(Project.MOSAIC_COLOR_CODES_FILENAME, false, false), proj.getJarStatus())) {
        	mosaicHash = HashVec.loadFileToHashString(proj.getFilename(Project.MOSAIC_COLOR_CODES_FILENAME, false, false), new int[] {0,1}, new int[] {2,3}, false, "\t", true, proj.getJarStatus(), true);
        } else {
        	System.err.println("Warning - could not find "+proj.getFilename(Project.MOSAIC_COLOR_CODES_FILENAME, false, false)+"; no annotation possible");
        	mosaicHash = new Hashtable<String,String>();
        }
        
    	if (listOfMosaicArms.toLowerCase().endsWith("/all")) {
            sampleList = sampleData.getListOfSamples();
            listOfArms = new String[sampleList.length*22*2][2];
            for (int i = 0; i<sampleList.length; i++) {
            	for (int j = 0; j<22; j++) {
            		listOfArms[i*22*2+j*2+0][0] = sampleList[i];
            		listOfArms[i*22*2+j*2+0][1] = "chr"+(j+1)+"p";
            		listOfArms[i*22*2+j*2+1][0] = sampleList[i];
            		listOfArms[i*22*2+j*2+1][1] = "chr"+(j+1)+"q";
                }
            }    		
    	} else {
    		list = new Vector<String[]>();
            try {
	    		reader = new BufferedReader(new FileReader(listOfMosaicArms));
		        line = reader.readLine().trim().split("[\\s]+");
		        if (!ext.checkHeader(line, new String[] {"Sample", "Arm"}, false)) {
		        	return;
		        }

		        while (reader.ready()) {
		        	line = reader.readLine().trim().split("[\\s]+");
		        	list.add(new String[] {line[0], line[1]});
		        }
		        reader.close();
            
            } catch (FileNotFoundException fnfe) {
    	        System.err.println("Error: file \""+listOfMosaicArms+"\" not found in current directory");
    	        return;
            } catch (IOException ioe) {
    	        System.err.println("Error reading file \""+listOfMosaicArms+"\"");
    	        return;
            }
            listOfArms = Matrix.toStringArrays(list);
    	}
        
        v = new Vector<String>();
        try {
	        writer = new PrintWriter(new FileWriter(ext.rootOf(listOfMosaicArms, false)+"_counts.xln"));
	        writer.println(Array.toStr(HEADER));
	        for (int i = 0; i<listOfArms.length; i++) {
	            indiPheno = sampleData.getIndiPheno(listOfArms[i][0]);

	    		if (indiPheno == null) {
	        		HashVec.addIfAbsent(listOfArms[i][0], v);
	        	} else {
		        	chr = Integer.parseInt(listOfArms[i][1].substring(3, listOfArms[i][1].length()-1));
		        	if (listOfArms[i][1].charAt(listOfArms[i][1].length()-1) == 'p') {
		        		arm = new Segment((byte)chr, 0, Positions.CENTROMERE_MIDPOINTS[chr]);
		        	} else {
		        		arm = new Segment((byte)chr, Positions.CENTROMERE_MIDPOINTS[chr], Positions.CHROMOSOME_LENGTHS[chr]);
		        	}
		        	
		        	cnvs = indiPheno.getCNVs(0, chr);
					if (cnvs == null) {
						cnvs = new CNVariant[0];
					}

					count = 0;
		        	sum = 0;
		        	for (int j = 0; j<cnvs.length; j++) {
	        			if (cnvs[j].overlaps(arm)) {
	        				count++;
	        				sum += cnvs[j].amountOfOverlapInBasepairs(arm);
	        			}
                    }
		        	proportion = (double)sum/(double)arm.getSize();
		        	writer.print(listOfArms[i][0]+"\t"+listOfArms[i][1]);
		        	writer.print("\t"+count+"\t"+sum+"\t"+proportion+"\t"+(400*proportion+count));
		        	if (lrrsdHash.containsKey(listOfArms[i][0])) {
			        	writer.print("\t"+lrrsdHash.get(listOfArms[i][0])+"\t"+(Double.parseDouble(lrrsdHash.get(listOfArms[i][0]))>0.28?1:0));
		        	} else {
			        	writer.print("\t.\t.");
		        	}
		        	if (mosaicHash.containsKey(listOfArms[i][0]+"\t"+listOfArms[i][1])) {
			        	writer.print("\t1\t"+mosaicHash.get(listOfArms[i][0]+"\t"+listOfArms[i][1]));
		        	} else {
			        	writer.print("\t0\t.\t.");
		        	}
		        	writer.println();
	        	}
	        }
            writer.close();
            
            if (v.size() > 0) {
                writer = new PrintWriter(new FileWriter(proj.getProjectDir()+"SAMPLES_IN_CNVFILE_NOT_IN_SAMPLE_DATA.txt"));
                for (int i = 0; i<v.size(); i++) {
                	writer.println(v.elementAt(i));
                }
                writer.close();
    			JOptionPane.showMessageDialog(null, "There were "+v.size()+" samples not present in the SampleData file; check file in the project directory for a list", "Error", JOptionPane.ERROR_MESSAGE);

            }
        } catch (IOException ioe) {
	        System.err.println("Error writing to file \""+ext.rootOf(listOfMosaicArms, false)+"_counts.xln"+"\"");
	        return;
        }
        System.out.println("Finished in "+ext.getTimeElapsed(time));
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = Project.DEFAULT_PROJECT;
		Project proj;
		boolean check = false;

		String usage = "\n"+
		"filesys.ParseIllumina requires 0-1 arguments\n"+
		"   (1) project file (i.e. proj="+filename+" (default))\n"+
		"   (2) check for overlap between mosaic arms and CNV calls in the first CNV file listed in the project file (i.e. -check (not the default))\n"+
		"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
		        return;
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("-check")) {
				check = true;
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
	        return;
		}
		
		try {
			proj = new Project(filename, false);
			if (check) {
				checkForOverlap(proj);
			} else {
				findOutliers(proj);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
