// make sure to use the scratch on the actual node for the large temp files Beagle creates, otherwise it hangs the whole cluster 
package gwas;

import java.io.*;
import java.util.*;

import common.*;
import filesys.*;

public class Beagle {
	public static final double[] PI_HAT_THRESHOLDS = {0.3, 0.5};
	public static final String[] SEGMENT_HEADER = {"FID1", "IID1", "FID2", "IID2", "PHE", "CHR", "BP1", "BP2", "SNP1", "SNP2", "NSNP", "KB"};
	public static final String[] ALLELES = {"0", "A", "C", "G", "T"};

	public static final int[] FILTER_CM_THRESHOLDS = {1,2,3,4,5};
	public static final double FILTER_MAX_THRESHOLD = 0.90;
	
	public static void pairUp(String filename) {
        PrintWriter writer;
        String[] files;
        String[][] ids;
        String root;
        
        files = filename.split(",");
        if (files.length > 2) {
        	System.err.println("Error - can only handle 1 or 2 files at a time");
        	return;
        }
        ids = new String[files.length][];
        
        root = "";
        for (int i = 0; i<files.length; i++) {
            if (!new File(files[i]).exists()) {
            	System.err.println("Error - '"+files[i]+"' does not exist");
            	return;
            }
            ids[i] = HashVec.loadFileToStringArray(files[i], false, new int[] {1}, false);
            
            root += (i==0?"":"_")+ext.rootOf(files[i]);
        }        

		try {
	        writer = new PrintWriter(new FileWriter(root+".list"));
	        if (files.length == 1) {
		        for (int i = 0; i<ids[0].length; i++) {
		        	for (int j = i+1; j<ids[0].length; j++) {
		        		writer.println(ids[0][i]+"\t"+ids[0][j]);
	                }
	            }
	        } else {
		        for (int i = 0; i<ids[0].length; i++) {
		        	for (int j = 0; j<ids[1].length; j++) {
		        		writer.println(ids[0][i]+"\t"+ids[1][j]);
	                }
	            }
	        }
	        writer.close();
        } catch (Exception e) {
	        System.err.println("Error writing to "+root+".list");
	        e.printStackTrace();
        }
	}
	
	public static void prepFiles(String list, int step) {
		new File("lists/").mkdirs();
    	Files.splitFile(list, (int)Math.ceil((double)Files.countLines(list, true)/(double)step), 0, step, "lists/", ".list", true);
		writeScript(list);
	}

	public static void writeScript(String filename) {
		PrintWriter writer;
		String root;
		
		root = ext.rootOf(filename);
		try {
	        writer = new PrintWriter(new FileWriter("script."+root));
	        writer.println("#!/bin/bash");
	        writer.println("#$ -cwd");
	        writer.println("#$ -S /bin/bash");
	        writer.println("");
	        writer.println("");
	        writer.println("echo \"start at: \" `date`");
	        writer.println("/bin/hostname");
	        writer.println(">plug");
	        writer.println("while [ -e \"plug\" ]; do ");
	        writer.println("	if [ -e \"wait\" ]; then");
	        writer.println("		sleep $(($RANDOM % 5 + 1))");
	        writer.println("	else");
	        writer.println("		>wait");
	        writer.println("		chr=1");
	        writer.println("		rep=1");
	        writer.println("		while [ -e \"chr$chr/$rep.log\" ] || [ -e \"chr$chr/$rep.ibds.pre_phase.bgl.ibd.gz\" ]; do");
	        writer.println("			chr=$(expr $chr + 1)");
	        writer.println("		        if [ ! -d \"chr$chr\" ]; then");
	        writer.println("		                chr=1");
	        writer.println("		                rep=$(expr $rep + 1)");
	        writer.println("		        fi");
	        writer.println("		done");
	        writer.println("		echo \"\" > chr$chr/$rep.log");
	        writer.println("		echo \"Next is chr=$chr and rep=$rep\"");
	        writer.println("		rm wait");
	        writer.println("		if [ ! -e \"lists/$rep.list\" ]; then");
	        writer.println("			echo \"$rep.list does not exist; exiting\"");
	        writer.println("			rm plug");
	        writer.println("			rep=$(expr $rep - 1)");
	        writer.println("		else");
	        writer.println("			cd chr$chr/");
	        writer.println("			echo \"Starting rep $rep at \" `date` > $rep_times.log");
	        writer.println("	        /bin/hostname > $rep.host");
	        writer.println("			tmpdir=temp.$rep");
	        writer.println("			mkdir $tmpdir");
	        writer.println("#			"+Files.JCP+"park.temp chr$chr/$tmpdir 4 $rep");
	        writer.println("			"+Files.JAVA+" -Djava.io.tmpdir=$tmpdir -d64 -Xmx1600M -jar /home/npankrat/bin/beagle.jar unphased=ibds.pre_phase.bgl ibdpairs=../lists/$rep.list out=$rep missing=0 markers=markers.dat seed=1234 > $rep.out");
	        writer.println("			rm -r $tmpdir");
	        writer.println("			echo \"Finished rep $rep at \" `date` >> $rep_times.log");
	        writer.println("			echo \"Starting to zip...\" >> $rep_times.log");
	        writer.println("			gzip $rep.ibds.pre_phase.bgl.ibd");
	        writer.println("			echo \"Finished zipping at \" `date` >> $rep_times.log");
	        writer.println("	        	if [ $rep -gt 1 ]; then");
	        writer.println("		        	rm $rep.ibds.pre_phase.bgl.phased.gz");
	        writer.println("		        	rm $rep.ibds.pre_phase.bgl.gprobs.gz");
	        writer.println("		        	rm $rep.ibds.pre_phase.bgl.r2");
	        writer.println("	        	fi");
	        writer.println("			cd ..");
	        writer.println("		fi");
	        writer.println("	fi		");
	        writer.println("done");
	        writer.println("echo \"Ran until \" `date`");
	        writer.println("echo \"Produced $rep replicates\"");
	        writer.close();
	        Files.chmod("script."+root);
        } catch (Exception e) {
	        System.err.println("Error writing to "+"script."+root);
	        e.printStackTrace();
        }
	}
	
	public static void batchIBD(String filename, int step) {
        PrintWriter writer;
        String root;
        String[] files;
        
        files = filename.split(",");
        if (files.length > 2) {
        	System.err.println("Error - can only handle 1 or 2 files at a time");
        	return;
        }
        
        root = "";
        for (int i = 0; i<files.length; i++) {
            if (!new File(files[i]).exists()) {
            	System.err.println("Error - '"+files[i]+"' does not exist");
            	return;
            }
            
            root += (i==0?"":"_")+ext.rootOf(files[i]);
        }        
        
        try {
	        writer = new PrintWriter(new FileWriter("batch."+root));
	        writer.println("#!/bin/bash");
	        writer.println();
	        writer.println("mkdir "+root);
	        for (int i = 0; i<files.length; i++) {
		        writer.println("cp "+files[i]+" "+root+"/");
            }
	        writer.println("cd "+root+"/");
	        writer.println("java -cp /home/npankrat/park.jar gwas.Beagle pair="+filename);
	        if (files.length == 1) {
		        writer.println("plink --bfile ../plink --keep "+filename+" --make-bed");
	        } else {
	        	Files.cat(files, root+".txt", Array.intArray(files.length, 0), new Logger());
		        writer.println("cp "+root+".txt "+root+"/");
		        writer.println("plink --bfile ../plink --keep "+root+".txt"+" --make-bed");
	        }
	        writer.println("splitByChrAlt plink &");
	        writer.println("splitByChrAlt2 plink");
	        writer.println("cp "+root+".list split/");
	        writer.println("cd split");
	        writer.println("java -cp /home/npankrat/park.jar gwas.Beagle -prepFiles list="+root+".list step="+step);
	        writer.println("for chr in {1..23}");
	        writer.println("do");
	        writer.println("echo \"Transposing chromosome $chr...\"");
	        writer.println("cd chr$chr/");
	        writer.println("/home/npankrat/bin/mine_ped_to_bgl plink.ped plink.map > ibds.pre_phase.bgl");
	        writer.println(Files.JCP+"parse.GenParser plink.bim out=markers.dat 1 2 4 5 noHeader");
	        writer.println("cd ..");
	        writer.println("done");
//	        writer.println("./script."+root);
	        writer.close();
	        Files.chmod("batch."+root);
        } catch (Exception e) {
	        System.err.println("Error writing to "+"batch."+root);
	        e.printStackTrace();
        }
	}
	
	public static void deleteFailedRuns() {
		String[] files;
		String trav;
		
		for (int chr = 1; chr<=22; chr++) {
			files = Files.list("chr"+chr+"/", ".log", false);
			for (int i = 0; i<files.length; i++) {
				trav = files[i].substring(0, files[i].lastIndexOf("."));
				if (!new File("chr"+chr+"/"+trav+".ibds.pre_phase.bgl.ibd.gz").exists()) {
					new File("chr"+chr+"/"+trav+".log").delete();
					new File("chr"+chr+"/"+trav+"_times.log").delete();
					new File("chr"+chr+"/"+trav+".host").delete();
				}
            }
			files = Files.list("chr"+chr+"/", ".ibds.pre_phase.bgl.phased.gz", false);
			for (int i = 0; i<files.length; i++) {
				if (!files[i].startsWith("1.")) {
					new File("chr"+chr+"/"+files[i]).delete();
				}
            }
        }
	}
	
	public static void checkProgress() {
		long[][] sizes;
		int[][] percentages;
		int numLists;
		double max;
		int compl, count, prob, miss;
		double average;
		
		numLists = Files.list("lists/", ".list", false).length;
		for (int rep = 1; rep<=numLists; rep++) {
			if (!new File("lists/"+rep+".list").exists()) {
				System.err.println("Error - abnormal list directory: missing '"+rep+".list'");
				return;
			}
        }
		
		sizes = new long[22][numLists];
		for (int chr = 1; chr<=22; chr++) {
			for (int rep = 1; rep<=numLists; rep++) {
				if (new File("chr"+chr+"/"+rep+".ibds.pre_phase.bgl.ibd.gz").exists()) {
					sizes[chr-1][rep-1] = new File("chr"+chr+"/"+rep+".ibds.pre_phase.bgl.ibd.gz").length();
				} else if (new File("chr"+chr+"/"+rep+".ibds.pre_phase.bgl.ibd").exists()) {
					sizes[chr-1][rep-1] = -3;
				} else if (new File("chr"+chr+"/"+rep+".log").exists()) {
					if (!new File("chr"+chr+"/temp."+rep).exists()) {
						sizes[chr-1][rep-1] = -5;
					} else if (new Date().getTime() - Files.getMostRecentUpdate("chr"+chr+"/temp."+rep) > 1800000) {
//						System.out.println("Current time: "+new Date().getTime()+"\tLast updated: "+Files.getMostRecentUpdate("chr"+chr+"/temp."+rep)+"\tDiff: "+(new Date().getTime() - Files.getMostRecentUpdate("chr"+chr+"/temp."+rep)));
						sizes[chr-1][rep-1] = -4;
					} else if (new File("chr"+chr+"/"+rep+".ibds.pre_phase.bgl.phased.gz").exists()) {
						sizes[chr-1][rep-1] = -1;
					} else {
						sizes[chr-1][rep-1] = -2;
					}
				} else {
					sizes[chr-1][rep-1] = -9;
				}
            }
        }
		
		percentages = new int[22][numLists];
		for (int chr = 1; chr<=22; chr++) {
			max = Math.max(Array.max(Array.toDoubleArray(sizes[chr-1])), 1);
			for (int rep = 1; rep<=numLists; rep++) {
				percentages[chr-1][rep-1] = (int)Math.round((double)sizes[chr-1][rep-1]/max*100);
			}
		}

		System.out.println("Reps:\t"+Array.toStr(Array.stringArraySequence(numLists, ""))+"\tNum remaining");
		for (int chr = 1; chr<=22; chr++) {
			compl = count = prob = miss = 0;
			max = Array.max(Array.toDoubleArray(sizes[chr-1]));
			System.out.print("chr"+chr);
			for (int rep = 1; rep<=numLists; rep++) {
//				average = 0;
				if (sizes[chr-1][rep-1] == -9) {
					System.out.print("\t.");
					miss++;
				} else if (sizes[chr-1][rep-1] == -5) {
					System.out.print("\tterm");
				} else if (sizes[chr-1][rep-1] == -4) {
					System.out.print("\tlost");
				} else if (sizes[chr-1][rep-1] == -3) {
					System.out.print("\ttrans");
				} else if (sizes[chr-1][rep-1] == -2) {
					System.out.print("\tprep");
				} else if (sizes[chr-1][rep-1] == -1) {
					System.out.print("\trun");
				} else {
					average = count = 0;
					for (int i = 1; i<=22; i++) {
						if (i != chr && percentages[i-1][rep-1] > 0) {
							average += percentages[i-1][rep-1];
							count++;
						}
                    }
					average /= (double)count;
					if (average - (double)percentages[chr-1][rep-1] > 2) {
						if (new File("chr"+chr+"/"+rep+".ibds.pre_phase.bgl.ibd").exists()) {
							System.out.print("\t("+percentages[chr-1][rep-1]+")");
						} else {
							System.out.print("\t"+percentages[chr-1][rep-1]+"**");
							prob++;
						}
					} else {
						System.out.print("\t"+percentages[chr-1][rep-1]);
					}
					compl++;
				}
//				System.out.print("\tp: "+ext.formDeci(percentages[chr-1][rep-1], 1)+"\ta: "+ext.formDeci(average, 1));
			}
			System.out.print("\t"+miss+"/"+(numLists-compl));
			if (prob > 0) {
				System.out.print("**");
			}
			System.out.println();
		}
	}

	public static void summarizeFile(String filename, String list) {
		BufferedReader reader;
        PrintWriter writer;
        PrintWriter[][][] writers;
        String[] line, markerNames;
        String[][] ids, pedIDs;
        String dir;
        int count;
        long time;
        double[][] avgs;
        double[][][] segAvgs, segMaxes;
        Logger log;
        int[][][] starts;
        String[][][][] segsInfo;
        double d;
        SnpMarkerSet markerSet;
        byte[] chrs, affs;
        int[] positions;
        double[] centiMorgans;
        Pedfile pedfile;
        FamilyStructure famStruct;
        Hashtable<String,String> hash, famids;
        boolean problem;
        int end;
        int[][] alleleCounts;
        
        log = new Logger(filename+"_sum.log");
        time = new Date().getTime();
        ids = HashVec.loadFileToStringMatrix(list, false, new int[] {0, 1}, false);
        dir = ext.parseDirectoryOfFile(filename);
        if (new File(dir+"plink.map").exists() || new File(dir+"plink.bim").exists()) {
        	markerSet = new SnpMarkerSet(new File(dir+"plink.map").exists()?dir+"plink.map":dir+"plink.bim", false, new Logger());
        	markerNames = markerSet.getMarkerNames();
        	chrs = markerSet.getChrs();
        	positions = markerSet.getPositions();
        	centiMorgans = markerSet.getCentiMorgans();
        } else {
        	log.reportError("Error - could not find plink.map or plink.bim in the same directory; required to define segments and perform a crucial datacheck");
        	return;
        }
        if (new File(dir+"plink.fam").exists() || new File(dir+"plink.ped").exists()) {
        	pedfile = new Pedfile(new File(dir+"plink.fam").exists()?dir+"plink.fam":dir+"plink.ped");
        	famStruct = pedfile.getFamilyStructure();
        	pedIDs = famStruct.getIds();
        	affs = famStruct.getAffections();
        	famids = new Hashtable<String,String>();
        	for (int i = 0; i<pedIDs.length; i++) {
        		famids.put(pedIDs[i][1], pedIDs[i][0]);        		
            }
        	hash = new Hashtable<String,String>();
        	for (int i = 0; i<pedIDs.length; i++) {
        		hash.put(pedIDs[i][1], affs[i]+"");        		
            }
        	problem = false;
        	for (int i = 0; i<ids.length; i++) {
        		affs = new byte[] {(byte)-999, (byte)-999};
        		for (int j = 0; j<2; j++) {
            		if (!hash.containsKey(ids[i][j])) {
            			log.reportError("Error - indiviudal "+ids[i][j]+" was not present in family stucture");
            			problem = true;
            		} else {
            			affs[j] = Byte.parseByte(hash.get(ids[i][j]));
            			if (affs[j] != 1 && affs[j] != 2) {
            				log.reportError("Error - indiviudal "+ids[i][j]+" has an invalid affection status: '"+hash.get(ids[i][j])+"'");
            				problem = true;
            			}
            		}
                }
        		hash.put(ids[i][0]+"\t"+ids[i][1], (affs[0]+affs[1]-3)+"");
            }
        	if (problem) {
        		return;
        	}
        } else {
        	log.reportError("Error - could not find plink.fam or plink.ped in the same directory; required to define segments and perform a crucial datacheck");
        	return;
        }
        
        
        log.report("Found "+markerNames.length+" markers to parse");
        
        avgs = new double[markerNames.length][2];
        alleleCounts = new int[markerNames.length][ALLELES.length];
        segAvgs = new double[ids.length][PI_HAT_THRESHOLDS.length][2];
        segMaxes = new double[ids.length][PI_HAT_THRESHOLDS.length][2];
        try {
	        reader = new BufferedReader(new FileReader(filename));
	        count = 0;
	        for (int i = 0; i<2; i++) {
		        line = reader.readLine().trim().split("[\\s]+");
		        for (int j = 0; j<ids.length; j++) {
		        	for (int k = 0; k<4; k++) {
			        	if (!line[j*4+k+1].equals(ids[j][i])) {
			        		log.reportError("Error - mismatched ids: expecting "+ids[j][i]+" in column "+(j*4+k+2)+", but found "+line[j*4+k+1]);
			        	}
                    }
                }
            }
	        
	        starts = new int[ids.length][][];
	        segsInfo = new String[ids.length][PI_HAT_THRESHOLDS.length][2][4];
	        for (int i = 0; i<ids.length; i++) {
		        starts[i] = Matrix.intMatrix(PI_HAT_THRESHOLDS.length, 2, -1);
            }
	        writers = new PrintWriter[PI_HAT_THRESHOLDS.length][2][2];
	        for (int j = 0; j<PI_HAT_THRESHOLDS.length; j++) {
	        	writers[j][0][0] = new PrintWriter(new FileWriter(filename+"_PI_"+PI_HAT_THRESHOLDS[j]+"_normal.segment"));
	        	writers[j][0][1] = new PrintWriter(new FileWriter(filename+"_PI_"+PI_HAT_THRESHOLDS[j]+"_normal.seginfo"));
//	        	writers[j][0].println(Array.toStr(SEGMENT_HEADER)+"\tSCORE");
	        	writers[j][1][0] = new PrintWriter(new FileWriter(filename+"_PI_"+PI_HAT_THRESHOLDS[j]+"_nuanced.segment"));
	        	writers[j][1][1] = new PrintWriter(new FileWriter(filename+"_PI_"+PI_HAT_THRESHOLDS[j]+"_nuanced.seginfo"));
//	        	writers[j][1].println(Array.toStr(SEGMENT_HEADER)+"\tSCORE");
            }
	        while (reader.ready()) {
	        	line = reader.readLine().trim().split("[\\s]+");
	        	if (line.length == ids.length*4+1) {
	        		for (int i = 0; i<ids.length; i++) {
	        			alleleCounts[count][ext.indexOfStr(line[i*4+2+1], ALLELES)]++; 
	        			for (int strict = 0; strict<2; strict++) {
		        			d = Double.parseDouble(line[i*4+strict+1]);
		        			avgs[count][strict] += d;
		        			for (int j = 0; j<PI_HAT_THRESHOLDS.length; j++) {
		        				if (d>PI_HAT_THRESHOLDS[j] && count < markerNames.length -1) {
		        					if (starts[i][j][strict] == -1) {
		        						starts[i][j][strict] = count;
		        						segAvgs[i][j][strict] = 0;
		        						segMaxes[i][j][strict] = 0;
		        						for (int k = 0; k<4; k++) {
			        						segsInfo[i][j][strict][k] = "";
										}
		        					}
	        						segAvgs[i][j][strict] += d;
	        						if (d > segMaxes[i][j][strict]) {
	        							segMaxes[i][j][strict] = d;
	        						}
	        						for (int k = 0; k<4; k++) {
		        						segsInfo[i][j][strict][k] += (segsInfo[i][j][strict][k].equals("")?"":"\t")+line[i*4+k+1];
									}
		        				} else if (starts[i][j][strict] != -1) {
		        					writers[j][strict][0].print(famids.get(ids[i][0])+"\t"+ids[i][0]+"\t"+famids.get(ids[i][1])+"\t"+ids[i][1]); // FID1 IID1 FID2 IID2
		        					writers[j][strict][0].print("\t"+hash.get(ids[i][0]+"\t"+ids[i][1])); // PHE  Phenotype concordance: -1,0,1
		        					writers[j][strict][0].print("\t"+chrs[starts[i][j][strict]]); // CHR
		        					end = count < markerNames.length -1?count - 1:count;
		        					writers[j][strict][0].print("\t"+positions[starts[i][j][strict]]+"\t"+positions[end]); // BP1 BP2
		        					writers[j][strict][0].print("\t"+markerNames[starts[i][j][strict]]+"\t"+markerNames[end]); // SNP1 SNP2
		        					writers[j][strict][0].print("\t"+(end-starts[i][j][strict]+1)); // NSNP  Number of SNPs in this segment
		        					writers[j][strict][0].print("\t"+(positions[end]-positions[starts[i][j][strict]]+1)); // KB  Physical length of segment (kb)
		        					writers[j][strict][0].print("\t"+ext.formDeci(segAvgs[i][j][strict]/(double)(end-starts[i][j][strict]+1), 5, true)); // extra quality score not used/allowed by PLINK
		        					writers[j][strict][0].print("\t"+ext.formDeci(segMaxes[i][j][strict], 3, true)); // extra max score not used/allowed by PLINK
		        					writers[j][strict][0].print("\t"+ext.formDeci(centiMorgans[end]-centiMorgans[starts[i][j][strict]], 4, true)); // cM  Genetic length of segment (cM)
		        					writers[j][strict][0].println();

		        					writers[j][strict][1].println(ids[i][0]+"\t"+ids[i][1]+"\t"+chrs[starts[i][j][strict]]+"\t"+starts[i][j][strict]+"\t"+end+"\tchr"+chrs[starts[i][j][strict]]+":"+positions[starts[i][j][strict]]+"-"+positions[end]);
		        					for (int k = 0; k<4; k++) {
		        						writers[j][strict][1].println(segsInfo[i][j][strict][k]);
		        					}
		        					
		        					starts[i][j][strict] = -1;
		        				}
	                        }
                        }
	        		}
	        		for (int strict = 0; strict<2; strict++) {
		        		avgs[count][strict] /= ids.length;
                    }
		        	count++;
	        	} else {
	        		log.reportError("Error - misamatched number of columns in row "+(count+3));
	        	}
	        }
	        if (count == markerNames.length) {
	        	log.report("Successfully processed all "+count+" lines");
	        } else {
	        	log.report("Error - only processed "+count+" of "+markerNames.length+" markers");
	        	reader.close();
	        	return;
	        }
	        reader.close();
	        for (int j = 0; j<PI_HAT_THRESHOLDS.length; j++) {
	        	writers[j][0][0].close();
	        	writers[j][0][1].close();
	        	writers[j][1][0].close();
	        	writers[j][1][1].close();
            }

        } catch (FileNotFoundException fnfe) {
        	log.reportError("Error: file \""+filename+"\" not found in current directory");
	        return;
        } catch (IOException ioe) {
        	log.reportError("Error reading file \""+filename+"\"");
        	log.reportException(ioe);
	        return;
        }
        
        try {
	        writer = new PrintWriter(new FileWriter(filename+"_sum.xln"));
	        for (int i = 0; i<markerNames.length; i++) {
	        	writer.println(markerNames[i]+"\t"+avgs[i][0]+"\t"+avgs[i][1]+"\t"+Array.toStr(alleleCounts[i]));
            }
	        writer.close();
        } catch (Exception e) {
        	log.reportError("Error writing to "+filename+"_sum.xln");
        	log.reportException(e);
        }
        log.report("Checked and summed in "+ext.getTimeElapsed(time));
        
	}
	
	public static void summarizeAllFiles(boolean del) {
		BufferedReader reader;
        PrintWriter writer;
        String[] line;
        String[][] markerNames;
        int count;
        double[][] avgs;
        Logger log;
        int rep;
        IntVector iv;
        IntVector[] missings;
        int n;
        SnpMarkerSet markerSet;
        byte[][] chrs;
        int[][] positions;
        String trav;
        IntVector missingChrs;
        int[][] alleleCounts;
        int[] actualAlleleCounts;
        String[][] markerFreq;
        String conAllele;
        double conFreq;
        
        log = new Logger("allSummary.log");
        iv = new IntVector();
        rep = 1;
        while (new File("lists/"+rep+".list").exists()) {
        	iv.add(Files.countLines("lists/"+rep+".list", false));
        	rep++;
        }
        n = Array.sum(iv.toArray());
        log.report("Detected "+iv.size()+" reps to parse for each chromosome");
        missings = IntVector.newIntVectors(22);
        markerNames = new String[22][];
        chrs = new byte[22][];
        positions = new int[22][];
        missingChrs = new IntVector();
        for (int chr = 1; chr<=22; chr++) {
            if (new File("chr"+chr+"/").exists()) {
                if (new File("chr"+chr+"/plink.map").exists() || new File("chr"+chr+"/plink.bim").exists()) {
                	markerSet = new SnpMarkerSet(new File("chr"+chr+"/plink.map").exists()?"chr"+chr+"/plink.map":"chr"+chr+"/plink.bim", false, new Logger());
                	markerNames[chr-1] = markerSet.getMarkerNames();
                	chrs[chr-1] = markerSet.getChrs();
                	positions[chr-1] = markerSet.getPositions();
                } else {
                	log.reportError("Error - could not find plink.map or plink.bim in chr"+chr+"/ -- required to define segments and perform a crucial datacheck");
                	return;
                }
            	for (int i = 0; i<iv.size(); i++) {
            		if (!new File("chr"+chr+"/"+(i+1)+".ibds.pre_phase.bgl.ibd_sum.xln").exists()) {
            			missings[chr-1].add(i+1);
            			if (del) {
                    		new File("chr"+chr+"/"+(i+1)+".ibds.pre_phase.bgl.ibd_sum.log").delete();
            			}
            		}
                }
            } else {
            	missingChrs.add(chr);
            }
        }
        if (missingChrs.size() > 0) {
        	log.reportError("Warning - missing data for the following chromosomes: "+ext.listRanges(missingChrs.toArray()));
        }
        for (int i = 0; i<missings.length; i++) {
        	if (missings[i].size() > 0) {
        		log.reportError("For chr"+(i+1)+", missing rep"+(missings[i].size()>1?"s":"")+": "+ext.listRanges(missings[i].toArray()));
        	}
        }

        log.report(ext.getTime()+"\tAveraging IBDs");
        try {
	        writer = new PrintWriter(new FileWriter("avgIBD.xln"));
	        writer.println("MarkerName\tChr\tPosition\t"+ext.replaceDirectoryCharsWithUnderscore(new File(".").getCanonicalPath(), 1)+"\t"+Array.toStr(ALLELES)+"\tConsensus\tFreq\tActualFreq\tpval");
	        writer.println("\t\t\t"+n);
	        for (int chr = 1; chr<=22; chr++) {
	            if (new File("chr"+chr+"/").exists()) {
		        	log.report(".", false, true);
		        	avgs = new double[markerNames[chr-1].length][2];
		        	alleleCounts = new int[markerNames[chr-1].length][ALLELES.length];
		        	for (rep = 1; rep<=iv.size(); rep++) {
		        		try {
		                    reader = new BufferedReader(new FileReader("chr"+chr+"/"+rep+".ibds.pre_phase.bgl.ibd_sum.xln"));
		                    count = 0;
		                    while (reader.ready()) {
		                    	line = reader.readLine().trim().split("[\\s]+");
		                    	if (!line[0].equals(markerNames[chr-1][count])) {
		                    		log.reportError("Error - mismatch in '"+"chr"+chr+"/"+rep+".ibds.pre_phase.bgl.ibd_sum.xln"+"' at line "+(count+1)+"; expecting "+markerNames[chr-1][count]+" and found "+line[0]);
		                    	}
		                    	avgs[count][0] += Double.parseDouble(line[1])*iv.elementAt(rep-1);
//		                    	avgs[count][1] += Double.parseDouble(line[2])*iv.elementAt(rep-1);
		                    	if (line.length > 3) {
			                    	for (int i = 0; i < ALLELES.length; i++) {
			                    		String trav2 = line[3+i];
			                    		alleleCounts[count][i] += Integer.parseInt(trav2);
									}
		                    	}
		                    	count++;
		                    }
		                    if (count != markerNames[chr-1].length) {
	                    		log.reportError("Error - mismatch in '"+"chr"+chr+"/"+rep+".ibds.pre_phase.bgl.ibd_sum.xln"+"' at line "+(count+1)+"; only processed "+count+" of "+markerNames.length+" markers");
			    	        }
		                    reader.close();
	                    } catch (FileNotFoundException fnfe) {
	                    	log.reportError("Error: file \""+"chr"+chr+"/"+rep+".ibds.pre_phase.bgl.ibd_sum.xln"+"\" not found in current directory");
	                    } catch (IOException ioe) {
	                    	log.reportError("Error reading file \""+"chr"+chr+"/"+rep+".ibds.pre_phase.bgl.ibd_sum.xln"+"\"");
	                    	log.reportException(ioe);
	                    }
	                }
		        	CmdLine.run("plink --bfile plink --freq", "chr"+chr+"/");
		        	markerFreq = HashVec.loadFileToStringMatrix("chr"+chr+"/plink.frq", true, new int[] {1,2,3,4,5}, false);
		        	for (int i = 0; i<markerNames[chr-1].length; i++) {
		        		avgs[i][0] /= n; 
		        		avgs[i][1] /= n;
		        		if (!markerNames[chr-1][i].equals(markerFreq[i][0])) {
		        			System.err.println("Error - mismatch in freq file order");
		        		}
		        		actualAlleleCounts = Array.subArray(alleleCounts[i], 1, alleleCounts[i].length);
		        		conAllele = ALLELES[ext.indexOfStr(Array.max(actualAlleleCounts)+"", Array.toStr(actualAlleleCounts).split("\t"))+1];
		        		writer.print(markerNames[chr-1][i]+"\t"+chrs[chr-1][i]+"\t"+positions[chr-1][i]+"\t"+avgs[i][0]+"\t"+Array.toStr(alleleCounts[i])+"\t"+conAllele+"\t"+ext.formDeci((double)Array.max(actualAlleleCounts)/(double)Array.sum(actualAlleleCounts), 5));
		        		if (conAllele.equals("0")) {
		        			writer.println("\t.\t.");
		        		} else {
		        			if (conAllele.equals(markerFreq[i][2])) {
		        				conFreq = 1-Double.parseDouble(markerFreq[i][3]);
		        			} else {
		        				conFreq = Double.parseDouble(markerFreq[i][3]);
		        			}
//	        				writer.println("\t"+ext.formDeci(conFreq, 5)+"\t"+ProbDist.ChiDist(ContingencyTable.ChiSquare(new double[][] {{Array.max(actualAlleleCounts), Array.sum(actualAlleleCounts) - Array.max(actualAlleleCounts)}, {conFreq*Double.parseDouble(markerFreq[i][4]), Double.parseDouble(markerFreq[i][4]) - conFreq*Double.parseDouble(markerFreq[i][4])}}, false, true), 1));
	        				writer.println("\t"+ext.formDeci(conFreq, 5)+"\t"+ext.formDeci(((double)Array.max(actualAlleleCounts)/(double)Array.sum(actualAlleleCounts))/conFreq, 5));
		        		}
		        			
	                }
		        	writer.flush();	            	
	            }
	        }
	        writer.close();
	        log.report("");
        } catch (Exception e) {
        	log.reportError("Error writing to "+"avgIBD.xln");
        	log.reportException(e);
        	return;
        }        

        log.report(ext.getTime()+"\tConcatenating segments");
        for (int i = 0; i<PI_HAT_THRESHOLDS.length; i++) {
        	for (int strict = 0; strict<2; strict++) {
        		for (int suffix = 0; suffix < 2; suffix++) {
            		trav = "PI_"+PI_HAT_THRESHOLDS[i]+"_"+(strict==0?"normal":"nuanced")+(suffix==0?".segment":".seginfo");
            		try {
    	                writer = new PrintWriter(new FileWriter(trav+(suffix==0?"Plus":"")));
    	                if (suffix == 0) {
    	                	writer.println(Array.toStr(SEGMENT_HEADER)+"\tMeanIBD\tMaxIBD\tGeneticLength");
    	                }
    		        	for (int chr = 1; chr<=22; chr++) {
    			            if (new File("chr"+chr+"/").exists()) {
    			        		for (rep = 1; rep<=iv.size(); rep++) {
    			        			try {
    		                            reader = new BufferedReader(new FileReader("chr"+chr+"/"+rep+".ibds.pre_phase.bgl.ibd_"+trav));
    		                            while (reader.ready()) {
    		                            	writer.println(reader.readLine());
    		                            }
    		                            reader.close();
    	                            } catch (FileNotFoundException fnfe) {
    	                            	log.reportError("Error: file \""+"chr"+chr+"/"+rep+".ibds.pre_phase.bgl.ibd_"+trav+"\" not found");
    	                            } catch (IOException ioe) {
    	            	            	log.reportException(ioe);
    	                            }
    	                        }
    			            }
                        }
    	                writer.close();
                    } catch (Exception e) {
    	            	log.reportError("Error writing to "+trav+"Plus");
    	            	log.reportException(e);
    	            	return;
                    }
                    
                    if (suffix == 0) {
	                    try {
	    	                writer = new PrintWriter(new FileWriter(trav));
	    	                try {
	    	                    reader = new BufferedReader(new FileReader(trav+"Plus"));
	    	                    while (reader.ready()) {
	    	                    	line = reader.readLine().trim().split("[\\s]+");
	    	                    	for (int j = 0; j<SEGMENT_HEADER.length; j++) {
	    	                    		writer.print((j==0?"":"\t")+line[j]);
	                                }
	    	                    	writer.println();
	    	                    }
	    	                    reader.close();
	                        } catch (FileNotFoundException fnfe) {
	    	                    System.err.println("Error: file \""+trav+"Plus"+"\" not found in current directory");
	    	                    System.exit(1);
	                        } catch (IOException ioe) {
	    	                    System.err.println("Error reading file \""+trav+"Plus"+"\"");
	    	                    System.exit(2);
	                        }
	    	                writer.close();
	                    } catch (Exception e) {
	    	                System.err.println("Error writing to "+trav);
	    	                e.printStackTrace();
	                    }
                    }
				}
            }
        }
	}

	public static void summarizeGermline() {
		BufferedReader reader;
        PrintWriter writer;
        String[] line;
        Logger log;
        Pedfile pedfile;
        FamilyStructure famStruct;
        String[][] pedIDs;
        byte[] affs;
        Hashtable<String,String> hash;
        boolean problem;

        log = new Logger("allGermline.log");

        if (new File("plink.fam").exists() || new File("plink.ped").exists()) {
        	pedfile = new Pedfile(new File("plink.fam").exists()?"plink.fam":"plink.ped");
        	famStruct = pedfile.getFamilyStructure();
        	pedIDs = famStruct.getIds();
        	affs = famStruct.getAffections();
        	hash = new Hashtable<String,String>();
        	for (int i = 0; i<pedIDs.length; i++) {
    			if (affs[i] != 1 && affs[i] != 2) {
    				log.reportError("Error - indiviudal "+pedIDs[i][0]+","+pedIDs[i][1]+" has an invalid affection status: '"+affs[i]+"'");
    				problem = true;
    			}
        		hash.put(pedIDs[i][0]+"\t"+pedIDs[i][1], affs[i]+"");        		
            }
        	problem = false;
        	if (problem) {
        		return;
        	}
        } else {
        	log.reportError("Error - could not find plink.fam or plink.ped in the same directory; required to define segments and perform a crucial datacheck");
        	return;
        }
        
        log.report(ext.getTime()+"\tConcatenating .match files");
		try {
            writer = new PrintWriter(new FileWriter("germline.segmentPlus"));
        	writer.println(Array.toStr(SEGMENT_HEADER)+"\tGeneticLength");
        	for (int chr = 1; chr<=22; chr++) {
    			try {
                    reader = new BufferedReader(new FileReader("chr"+chr+"/chr"+chr+".match"));
                    while (reader.ready()) {
                    	line = reader.readLine().trim().split("[\\s]+");
                    	writer.print(line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]);
                		affs = new byte[] {(byte)-999, (byte)-999};
                		for (int j = 0; j<2; j++) {
                    		if (!hash.containsKey(line[j*2+0]+"\t"+line[j*2+1])) {
                    			log.reportError("Error - indiviudal "+line[j*2+0]+","+line[j*2+1]+" was not present in family stucture");
                    			problem = true;
                    		} else {
                    			affs[j] = Byte.parseByte(hash.get(line[j*2+0]+"\t"+line[j*2+1]));
                    		}
                        }
                		writer.print("\t"+(affs[0]+affs[1]-3));
                    	writer.print("\t"+line[4]+"\t"+line[5]+"\t"+line[6]+"\t"+line[7]+"\t"+line[8]+"\t"+line[9]);
                		writer.print("\t"+(Integer.parseInt(line[6])+Integer.parseInt(line[5])+1));
                		writer.println("\t"+line[10]);
                    }
                    reader.close();
                } catch (FileNotFoundException fnfe) {
                	log.reportError("Error: file \""+"chr"+chr+"/chr"+chr+".match"+"\" not found");
                } catch (IOException ioe) {
	            	log.reportException(ioe);
                }
            }
            writer.close();
        } catch (Exception e) {
        	log.reportError("Error writing to "+"germline.segmentPlus");
        	log.reportException(e);
        	return;
        }
        
        try {
            writer = new PrintWriter(new FileWriter("germline.segment"));
            try {
                reader = new BufferedReader(new FileReader("germline.segmentPlus"));
                while (reader.ready()) {
                	line = reader.readLine().trim().split("[\\s]+");
                	for (int j = 0; j<SEGMENT_HEADER.length; j++) {
                		writer.print((j==0?"":"\t")+line[j]);
                    }
                	writer.println();
                }
                reader.close();
            } catch (FileNotFoundException fnfe) {
                System.err.println("Error: file \""+"germline.segmentPlus"+"\" not found in current directory");
                System.exit(1);
            } catch (IOException ioe) {
                System.err.println("Error reading file \""+"germline.segmentPlus"+"\"");
                System.exit(2);
            }
            writer.close();
        } catch (Exception e) {
            System.err.println("Error writing to "+"germline.segment");
            e.printStackTrace();
        }
	}
	
	// at some point, if we run into confusion, we might need to read in the phased data and flag haplotype strings instead 
	public static void align(String segInfoFile, String intervalFile, String segFile) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line, keys, alleles, probs, data, merge;
		Vector<String> v;
		Vector<String[]> dataV;
		Hashtable<String, Hashtable<String, Vector<String[]>>> hashes;
		Hashtable<String, Vector<String[]>> hash;
		Hashtable<String, Vector<String>> includes;
		Hashtable<String, Segment[]> includeSegs;
		String[][] intervalStartAndStopMarkers;
		SnpMarkerSet markerSet;
		String[] markerNames;
		int[][] intervalIndices;
		Segment[] intervals, segs;
		Segment seg;
		int[] positions, order;
		String dir;
		byte[] chrs;
		byte chr;
		int start, stop;
		boolean merged;
		
		chr = -1;
		dir = ext.parseDirectoryOfFile(segInfoFile);
		intervalStartAndStopMarkers = HashVec.loadFileToStringMatrix(intervalFile, false, new int[] {0,1}, false);
		if (segFile != null) {
			includes = HashVec.loadFileToHashVec(segFile, new int[] {1,3}, new int[] {5,6,7}, "\t", true, false);
			includeSegs = new Hashtable<String, Segment[]>();
			keys = HashVec.getKeys(includes);
			for (int i = 0; i < keys.length; i++) {
				v = includes.get(keys[i]);
				segs = new Segment[v.size()];
				for (int j = 0; j < segs.length; j++) {
					line = v.elementAt(j).split("[\\s]+");
					segs[j] = new Segment(Byte.parseByte(line[0]), Integer.parseInt(line[1]), Integer.parseInt(line[2]));
				}
				includeSegs.put(keys[i], segs);
			}
		} else {
			includes = null;
			includeSegs = null;
		}
		
		if (!new File(dir+"plink.map").exists()) {
			System.err.println("Error - need a map file to align segments and plink.map does not exist in directory '"+dir+"'");
			return;
		}
		markerSet = new SnpMarkerSet(dir+"plink.map");
		markerNames = markerSet.getMarkerNames();
		chrs = markerSet.getChrs();
		positions = markerSet.getPositions();
		
		intervalIndices = new int[intervalStartAndStopMarkers.length][2];
		intervals = new Segment[intervalStartAndStopMarkers.length];
		System.out.println("Aligning to the following segments:");
		for (int i = 0; i<intervalStartAndStopMarkers.length; i++) {
			for (int j = 0; j<2; j++) {
				intervalIndices[i][j] = ext.indexOfStr(intervalStartAndStopMarkers[i][j], markerNames);
				if (intervalIndices[i][j] == -1) {
					System.err.println("Error - marker '"+intervalStartAndStopMarkers[i][j]+"' not found in plink.map; aborting");
					System.exit(1);
				}
				if (chr == -1) {
					chr = chrs[intervalIndices[i][j]];
				} else if (chr != chrs[intervalIndices[i][j]]) {
					System.err.println("Error - markers are on different chromosomes; can only process one chromosome at a time");
				}
			}
			intervals[i] = new Segment(chr, positions[intervalIndices[i][0]], positions[intervalIndices[i][1]]);
			System.out.println(intervals[i].getUCSClocation());
		}
		order = Segment.quicksort(intervals);
		intervals = Segment.putInOrder(intervals, order);
		intervalIndices = Matrix.putInOrder(intervalIndices, order);
		
		hashes = new Hashtable<String, Hashtable<String,Vector<String[]>>>();
		for (int i = 0; i < intervals.length; i++) {
			hashes.put(i+"", new Hashtable<String, Vector<String[]>>());
		}
		try {
			reader = new BufferedReader(new FileReader(segInfoFile));
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				seg = new Segment(line[5]);
				if (Segment.overlapsAny(seg, intervals) && (includes == null || (includeSegs.containsKey(line[0]+"\t"+line[1]) && Segment.contains(seg, includeSegs.get(line[0]+"\t"+line[1]))))) {
					reader.readLine(); // normal confidence
					probs = reader.readLine().trim().split("[\\s]+"); // nuanced confidence
					alleles = reader.readLine().trim().split("[\\s]+"); // alleles
					reader.readLine(); // genotypes
					
					for (int i = 0; i < intervals.length; i++) {
						start = Integer.parseInt(line[3]);
						stop = Integer.parseInt(line[4]);
						
						if (seg.overlaps(intervals[i])) {
							hash = hashes.get(i+"");
							
							if (positions[start] != seg.getStart()) {
								System.err.println("Error - different start; must have used a different map to generate segInfo file");
							} else if (positions[stop] != seg.getStop()) {
								System.err.println("Error - different stop; must have used a different map to generate segInfo file");
							} else {
								System.out.print(".");
							}
							
							data = new String[intervalIndices[i][1]-intervalIndices[i][0]+2];
							for (int j = start; j < stop; j++) {
								if (j >= intervalIndices[i][0] && j <= intervalIndices[i][1]) {
									if (Double.parseDouble(probs[j-start]) >= 0.500) {
										data[j-intervalIndices[i][0]+1] = alleles[j-start];
									}
								}
							}
							
							for (int j = 0; j < 2; j++) {
								data[0] = line[1-j];
								HashVec.addToHashArrayVec(hash, line[j], Array.clone(data));
							}
						}
					}
				} else {
					reader.readLine();
					reader.readLine();
					reader.readLine();
					reader.readLine();
					System.out.print("X");
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + segInfoFile + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + segInfoFile + "\"");
			System.exit(2);
		}
		
		for (int i = 0; i < intervals.length; i++) {
			try {
				writer = new PrintWriter(new FileWriter(dir+"SegAlignment_"+ext.replaceAllWith(intervals[i].getUCSClocation(), ":", "@")+".xln"));
				writer.print("IID\tPairing(s)");
				for (int j = intervalIndices[i][0]; j <= intervalIndices[i][1]; j++) {
					writer.print("\t"+markerNames[j]);
				}
				writer.println();
				hash = hashes.get(i+"");
				keys = HashVec.getKeys(hash);
				for (int j = 0; j < keys.length; j++) {
					dataV = hash.get(keys[j]);
					do {
						merged = false;
						for (int i2 = 0; i2 < dataV.size(); i2++) {
							for (int j2 = i2+1; j2 < dataV.size(); j2++) {
								merge = Array.merge(dataV.elementAt(i2), dataV.elementAt(j2), dataV.elementAt(i2).length/10);
								if (merge != null) {
									dataV.removeElementAt(j2);
									dataV.removeElementAt(i2);
									dataV.insertElementAt(merge, i2);
									merged = true;
								}
							}
						}						
					} while (merged);
					for (int k = 0; k < dataV.size(); k++) {
						writer.println(keys[j]+"\t"+Array.toStr(dataV.elementAt(k), null, "\t", "."));
					}
				}
				writer.close();
			} catch (Exception e) {
				System.err.println("Error writing to " + dir+"SegAlignment_"+ext.replaceAllWith(intervals[i].getUCSClocation(), ":", "@")+".xln");
				e.printStackTrace();
			}
			Files.transpose(dir+"SegAlignment_"+ext.replaceAllWith(intervals[i].getUCSClocation(), ":", "@")+".xln", "\t", null, "\t", new Logger());
		}
		
		
		// delineate positions and indices of markers from mapfile
		// create a hashtable of all regions overlapping these intervals, one record per person
		// records need to fill the interval (sheer and fill in blanks)
		// records can be added to, if there is significant deviation, may need to delineate both alleles in the file
		// add/display an individual only if they overlap the region, display minimum number of alleles 
	}
	
	public static void filter(String filename) {
		BufferedReader reader;
        PrintWriter writer;
        String[] line;
        Logger log;
        
        log = new Logger("filter.log");
		
        log.report(ext.getTime()+"\tFiltering segments");
        for (int i = 0; i<FILTER_CM_THRESHOLDS.length; i++) {
        	for (int strict = 0; strict<2; strict++) {
        		try {
					reader = new BufferedReader(new FileReader(filename));
					writer = new PrintWriter(new FileWriter(FILTER_CM_THRESHOLDS[i]+"cM"+(strict==1?"_max"+(int)(FILTER_MAX_THRESHOLD*100):"")+".segment"));
					line = reader.readLine().trim().split("[\\s]+");
					writer.println(Array.toStr(Array.subArray(line, 0, line.length-3)));
					while (reader.ready()) {
						line = reader.readLine().trim().split("[\\s]+");
						if (Double.parseDouble(line[14]) >= (double)FILTER_CM_THRESHOLDS[i] && (strict==0 || Double.parseDouble(line[13]) >= (double)FILTER_MAX_THRESHOLD)) {
							writer.println(Array.toStr(Array.subArray(line, 0, line.length-3)));
						}
					}
					writer.close();
					reader.close();
				} catch (FileNotFoundException fnfe) {
					System.err.println("Error: file \"" + filename + "\" not found in current directory");
					System.exit(1);
				} catch (IOException ioe) {
					System.err.println("Error reading file \"" + filename + "\"");
					System.exit(2);
				}
            }
        }
	}

	public static void checkRelatedness(String groupsFile, String plinkFile, String segFile) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		Hashtable<String, String> groupsHash, checks, extras;
		String[] ids;
		double[][] means;
		int[][] meanCounts;
		CountVector cv;
		String[] values;
		int[] counts, order;
		int index;
		
		groupsHash = HashVec.loadFileToHashString(groupsFile, new int[] {0,1}, new int[] {2}, false, "\t", false, false, false);
		groupsHash.remove("FID\tIID");
		ids = HashVec.getKeys(groupsHash);
		cv = new CountVector(groupsHash);
		values = cv.getValues();
		counts = cv.getCounts();
		
		order = Sort.quicksort(values);
		values = Sort.putInOrder(values, order);
		counts = Sort.putInOrder(counts, order);
		
		means = new double[2][values.length+1];
		meanCounts = new int[2][values.length+1];
		
		if (plinkFile != null) {
			checks = new Hashtable<String, String>();
			extras = new Hashtable<String, String>();
			for (int i = 0; i < ids.length; i++) {
				for (int j = i+1; j < ids.length; j++) {
					checks.put(ids[i]+"\t"+ids[j], "");
				}
			}
			try {
				reader = new BufferedReader(new FileReader(plinkFile));
				ext.checkHeader(reader.readLine().trim().split("[\\s]+"), Plink.CLUSTER_HEADER, true);
				while (reader.ready()) {
					line = reader.readLine().trim().split("[\\s]+");
					if (groupsHash.containsKey(line[0]+"\t"+line[1]) && groupsHash.containsKey(line[2]+"\t"+line[3])) {
						if (checks.containsKey(line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]) || checks.containsKey(line[2]+"\t"+line[3]+"\t"+line[0]+"\t"+line[1])) {
							checks.remove(line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]);
							checks.remove(line[2]+"\t"+line[3]+"\t"+line[0]+"\t"+line[1]);
							
							if (groupsHash.get(line[0]+"\t"+line[1]).equals(groupsHash.get(line[2]+"\t"+line[3]))) {
								index = ext.indexOfStr(groupsHash.get(line[0]+"\t"+line[1]), values);
								means[0][index] += Double.parseDouble(line[9]);
								meanCounts[0][index]++;
							} else {
								means[0][values.length] += Double.parseDouble(line[9]);
								meanCounts[0][values.length]++;
							}
							
						} else {
							System.err.println("Error that should not happen - pair present in cluster file twice ("+line[2]+"-"+line[3]+" and "+line[0]+"-"+line[1]+")");
						}
					} else if (!groupsHash.containsKey(line[0]+"\t"+line[1])) {
						extras.put(line[0]+"\t"+line[1], "");
					} else if (!groupsHash.containsKey(line[2]+"\t"+line[3])) {
						extras.put(line[2]+"\t"+line[3], "");
					} else {
						System.err.println("Error that cannot happen");
					}
				}
				reader.close();

				if (extras.size() > 0) {
					Files.writeList(HashVec.getKeys(extras), ext.parseDirectoryOfFile(groupsFile)+"SAMPLES_IN_"+ext.removeDirectoryInfo(plinkFile)+"_BUT_NOT_IN_"+ext.removeDirectoryInfo(groupsFile)+".TXT");
				}
				if (checks.size() > 0) {
					Files.writeList(HashVec.getKeys(extras), ext.parseDirectoryOfFile(groupsFile)+"SAMPLE_PAIRS_FROM_"+ext.removeDirectoryInfo(groupsFile)+"_BUT_NOT_IN_"+ext.removeDirectoryInfo(plinkFile)+".TXT");
				}
			} catch (FileNotFoundException fnfe) {
				System.err.println("Error: file \"" + plinkFile + "\" not found in current directory");
				System.exit(1);
			} catch (IOException ioe) {
				System.err.println("Error reading file \"" + plinkFile + "\"");
				System.exit(2);
			}
		}

		try {
			writer = new PrintWriter(new FileWriter(ext.rootOf(groupsFile, false)+"_ibd.xln"));
			writer.println("Group\tN\tnCr(2)\tmeanNumSegments\tmeanTotalSegmentsSize");
			for (int i = 0; i < values.length+1; i++) {
				writer.print(i==values.length?"across\t.":values[i]+"\t"+counts[i]);
				if (plinkFile == null) {
					writer.print("\t.\t.");
				} else {
					writer.print("\t"+meanCounts[0][i]+"\t"+ext.formDeci(means[0][i]/(double)meanCounts[0][i], 5));
				}
				if (segFile == null) {
					writer.print("\t.\t.");
				} else {
					
				}
				writer.println();
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + ext.rootOf(groupsFile, false)+"_ibd.xln");
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args) {
	    int numArgs = args.length;
	    String pairup = "";
	    boolean prepFiles = false;
	    String list = "all.txt";
	    int step = 10000;
	    boolean batchIBD = false;
	    boolean writeScript = false;
	    boolean del = false;
	    boolean check = false;
	    String sumFile = "";
	    boolean sumAll = false;
	    boolean sumGermline = false;
	    String segInfoFile = "";
	    String intervalFile = "";
	    String segFile = null;
	    String filter = "";
	    String groupsFile = null;
	    String plinkFile = null;
	    
	    String usage = "\n"+
	    "gwas.Beagle requires 0-1 arguments\n"+
	    "   (1) pair up IIDs from a PLINK keep file (i.e. pair=keeps.txt (not the default))\n"+
	    "  OR\n"+
	    "   (1) generate list files and script for IBD run (i.e. -prepFiles (not the default))\n"+
	    "   (2) list of pairs (i.e. list="+list+" (default))\n"+
	    "   (3) number of pairs to do at a time (i.e. step="+step+" (default))\n"+
	    "  OR\n"+
	    "   (1) writeScript (i.e. -writeScript (not the default))\n"+
	    "  OR\n"+
	    "   (1) create batch file for entire IBD process (i.e. -batchIBD (not the default))\n"+
	    "  OR\n"+
	    "   (1) delete .log files for failed runs (i.e. -del (not the default))\n"+
	    "  OR\n"+
	    "   (1) check progress on IBD runs (i.e. -check (not the default))\n"+
	    "  OR\n"+
	    "   (1) summarize the IBD estimates of a single file (i.e. sumFile=10.ibds.pre_phase.bgl.ibd (not the default))\n"+
	    "  OR\n"+
	    "   (1) sumamrize all of the summarized IBD files (i.e. -sumAll (not the default))\n"+
	    "  OR\n"+
	    "   (1) sumamrize all GERMLINE .match files (i.e. -sumGermline (not the default))\n"+
	    "  OR\n"+
	    "   (1) filter a .segmentPlus (i.e. filter=PI_0.5_normal.segmentPlus (not the default))\n"+
	    "  OR\n"+
	    "   (1) align .seginfo files at specific intervals (i.e. align=1.ibds.pre_phase.bgl.ibd_PI_0.5_normal.seginfo (not the default))\n"+
	    "   (2) list of those intervals containing two columns of start and stop markerNames  (i.e. ints=intervals.txt (not the default))\n"+
	    "   (3) (optional) only use segments in the following file (i.e. seg=PI_0.5_normal.segment (not the default))\n"+
	    "  OR\n"+
	    "   (1) check relatedness between groups of individuals (3 cols: FID IID group) (i.e. related=ids_to_check.txt (not the default))\n"+
	    "   (2) plink cluste file to check (i.e. plink=cluster.genome (not the default))\n"+
	    "   (3) segments to use in calculation (i.e. seg=PI_0.5_normal.segment (not the default))\n"+
	    "";

	    for (int i = 0; i<args.length; i++) {
		    if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
			    System.err.println(usage);
			    System.exit(1);
		    } else if (args[i].startsWith("pair=")) {
			    pairup = args[i].split("=")[1];
			    numArgs--;
		    } else if (args[i].startsWith("-prepFiles")) {
		    	prepFiles = true;
			    numArgs--;
		    } else if (args[i].startsWith("list=")) {
			    list = args[i].split("=")[1];
			    numArgs--;
		    } else if (args[i].startsWith("step=")) {
			    step = ext.parseIntArg(args[i]);
			    numArgs--;
		    } else if (args[i].startsWith("-batchIBD")) {
			    batchIBD = true;
			    numArgs--;
		    } else if (args[i].startsWith("-writeScript")) {
		    	writeScript = true;
			    numArgs--;
		    } else if (args[i].startsWith("-del")) {
			    del = true;
			    numArgs--;
		    } else if (args[i].startsWith("-check")) {
			    check = true;
			    numArgs--;
		    } else if (args[i].startsWith("sumFile=")) {
			    sumFile = args[i].split("=")[1];
			    numArgs--;
		    } else if (args[i].startsWith("-sumAll")) {
		    	sumAll = true;
			    numArgs--;
		    } else if (args[i].startsWith("-sumGermline")) {
		    	sumGermline = true;
			    numArgs--;
		    } else if (args[i].startsWith("filter=")) {
			    filter = args[i].split("=")[1];
			    numArgs--;
		    } else if (args[i].startsWith("align=")) {
			    segInfoFile = args[i].split("=")[1];
			    numArgs--;
		    } else if (args[i].startsWith("ints=")) {
			    intervalFile = args[i].split("=")[1];
			    numArgs--;
		    } else if (args[i].startsWith("seg=")) {
			    segFile = args[i].split("=")[1];
			    numArgs--;
		    } else if (args[i].startsWith("related=")) {
			    groupsFile = ext.parseStringArg(args[i], null);
			    numArgs--;
		    } else if (args[i].startsWith("plink=")) {
			    groupsFile = ext.parseStringArg(args[i], null);
			    numArgs--;
		    } else {
		    	System.err.println("Error - don't know what to do with argument: "+args[i]);
		    }
	    }
	    if (numArgs!=0) {
		    System.err.println(usage);
		    System.exit(1);
	    }

//	    sumAll = true;
	    
//	    segInfoFile = "D:\\Segments\\targetedRegions\\chr22seg\\chr1\\reparse\\1.ibds.pre_phase.bgl.ibd_PI_0.3_normal.seginfo";
//	    intervalFile = "D:\\Segments\\targetedRegions\\chr22seg\\chr1\\reparse\\intervalStartAndStopMarkers.txt";

//	    segInfoFile = "D:\\Segments\\targetedRegions\\chr22\\PI_0.5_normal.seginfo";
//	    intervalFile = "D:\\Segments\\targetedRegions\\chr22\\intervals.txt";
//	    segFile = "D:\\Segments\\targetedRegions\\chr22\\5cM.segment";
	    
//	    groupsFile = "ancestryAffection.xln";
//	    plinkFile = "filtered_cluster.genome";
	    
	    try {
	    	if (batchIBD) {
	    		batchIBD(list, step);
	    	} else if (!pairup.equals("")) {
	    		pairUp(pairup);
	    	} else if (prepFiles) {
	    		prepFiles(list, step);
	    	} else if (writeScript) {
	    		writeScript(list);
	    	} else if (check) {
		    	checkProgress();
	    	} else if (!sumFile.equals("")) {
	    		summarizeFile(sumFile, list);
	    	} else if (sumAll) {
	    		summarizeAllFiles(del);
	    	} else if (sumGermline) {
	    		summarizeGermline();
	    	} else if (del) {
	    		deleteFailedRuns();
	    	} else if (!filter.equals("")) {
	    		filter(filter);
	    	} else if (!segInfoFile.equals("")) {
	    		align(segInfoFile, intervalFile, segFile);
	    	} else if (groupsFile != null) {
	    		checkRelatedness(groupsFile, plinkFile, segFile);
	    	}
	    } catch (Exception e) {
		    e.printStackTrace();
	    }
    }
}
