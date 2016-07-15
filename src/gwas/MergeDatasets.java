// you have to force A1 onto the same strand if you want to use the HWE shortcut 
// -Xms1024M -Xmx1024M
package gwas;

import java.io.*;
import java.util.*;

import stats.ContingencyTable;
import stats.ProbDist;
import stats.Rscript;
import common.*;
import filesys.SerialHash;

public class MergeDatasets {
	public static final double HOMOGENEITY_THRESHOLD = 0.001;
//	public static final double LOOSE_HOMOGENEITY_THRESHOLD = 0.015;
	public static final double LOOSE_HOMOGENEITY_THRESHOLD = 0.05;
	
	public static final String CHI_SQUARE_DROPS_FILENAME = "lackOfHomogeneity.dat";
	public static final String FISHER_OR_CHI_SQUARE_DROPS_FILENAME = "FisherOrChiSquareDrops.dat";
	
	public static void checkForHomogeneity(String dir) {
		checkForHomogeneity(dir, null, null, "UNAFF", new Logger());
	}

	public static void checkForHomogeneity(String dir, String[] dirs, String outputDir,String hweCountFlag, Logger log) {
        Hashtable<String,Hashtable<String,String>> hashes;
        Vector<String> v = new Vector<String>();
        Hashtable<String,String> hash;
        BufferedReader reader;
		PrintWriter writer, writer2, writer3, writerProblematic;
        String[] line, markerNames, keys, refAlleles, gCounts;
        int[][] alleleCounts, genotypeCounts;
        double[] freqs;
        int index, temp;
        double p;
		String filename;
//		System.out.println(Array.toStr(dirs));
		if (dirs == null) {
			dirs = Files.listDirectories(dir, false);
		} else if (dir == null) {
			dir = "";
		} else {
			log.reportTimeError("Can only provide a single top directory, or multiple data dirs");
			return;
		}
		outputDir = outputDir == null ? dir : outputDir;
        for (int i = 0; i<dirs.length; i++) {
        	dirs[i] = ext.verifyDirFormat(dirs[i]);
        	System.out.println(dir+dirs[i]+"hardy.hwe");
        	if (new File(dir+dirs[i]+"hardy.hwe").exists() || new File(dir+dirs[i]+"plink.hwe").exists()) {
        		v.add(dirs[i]);
        	}
        }
        dirs = Array.toStringArray(v);
        
        hashes = new Hashtable<String,Hashtable<String,String>>();
        for (int i = 0; i<dirs.length; i++) {
			filename = (Files.exists(dir+dirs[i]+"hardy.hwe")?"hardy.hwe":"plink.hwe");
    		try {
    			reader = new BufferedReader(new FileReader(dir+dirs[i]+filename));
    			System.out.println(ext.getTime()+"\tLoading "+filename+" in "+dirs[i]);
    			ext.checkHeader(reader.readLine().trim().split("[\\s]+"), MarkerQC.HWE_HEADER, true);
    			while (reader.ready()) {
    				line = reader.readLine().trim().split("[\\s]+");
    				if (line[2].equals(hweCountFlag) || line[2].equals("ALL(NP)") || line[2].equals("ALL(QT)")) {
    					HashVec.addToHashHash(hashes, line[1].toLowerCase(), i+"", line[3]+" "+line[4]+" "+line[5]);
    				}
    			}
    			reader.close();
    		} catch (FileNotFoundException fnfe) {
    			System.err.println("Error: file \""+dirs[i]+filename+"\" not found in "+dir);
    			System.exit(1);
    		} catch (IOException ioe) {
    			System.err.println("Error reading file \""+dirs[i]+filename+"\"");
    			System.exit(2);
    		}
        }
        
        markerNames = HashVec.getKeys(hashes, false, false);
        try {
	        writer = new PrintWriter(new FileWriter(outputDir+CHI_SQUARE_DROPS_FILENAME));
	        writer2 = new PrintWriter(new FileWriter(outputDir+"homogeneityTests.xln"));
	        
	        writer2.println("SNP\t"+Array.toStr(dirs)+"\t"+Array.toStr(dirs)+"\t2allele_p-value\t3genotype_p-value");
	        writer3 = new PrintWriter(new FileWriter(outputDir+"homo.R"));
	        
			writerProblematic = new PrintWriter(new FileWriter(outputDir + "problematic.dat"));
	        System.out.println(ext.getTime()+"\tComputing tests of homogeneity...");
	        Hashtable<String, String> problematicTrack = new Hashtable<String, String>();
	        for (int i = 0; i<markerNames.length; i++) {
	        	hash = hashes.get(markerNames[i]);
	        	keys = HashVec.getKeys(hash, false, false);
	        	alleleCounts = new int[keys.length][2];
	        	genotypeCounts = new int[keys.length][2];
	        	gCounts = Array.stringArray(dirs.length);
	        	freqs = Array.doubleArray(dirs.length, Double.MIN_VALUE);
	        	refAlleles = new String[2];
	        	for (int j = 0; j<keys.length; j++) {
	        		index = Integer.parseInt(keys[j]);
	        		line = hash.get(keys[j]).split("[\\s]+");
	        		genotypeCounts[j] = Array.toIntArray(line[2].split("/"));

	        		switch (Metal.determineStrandConfig(new String[] {line[0], line[1]}, refAlleles)) {
	        		case Metal.STRAND_CONFIG_SAME_ORDER_FLIPPED_STRAND:
//	        			countFlipped[i]++;
	        		case Metal.STRAND_CONFIG_SAME_ORDER_SAME_STRAND:
	                    break;
	        		case Metal.STRAND_CONFIG_OPPOSITE_ORDER_FLIPPED_STRAND:
//	        			countFlipped[i]++;
	        		case Metal.STRAND_CONFIG_OPPOSITE_ORDER_SAME_STRAND:
	        			temp = genotypeCounts[j][0];
	        			genotypeCounts[j][0] = genotypeCounts[j][2];
	        			genotypeCounts[j][2] = temp;
	                    break;
	        		case Metal.STRAND_CONFIG_BOTH_NULL:
	        			break;
	        		case Metal.STRAND_CONFIG_DIFFERENT_ALLELES:
						if (!problematicTrack.containsKey(markerNames[i] + "\t" + Metal.STRAND_CONFIG_DIFFERENT_ALLELES)) {
							writerProblematic.println(markerNames[i] + "\t" + Metal.STRAND_CONFIG_DIFFERENT_ALLELES + "\t" + dirs[index]);
							problematicTrack.put(markerNames[i] + "\t" + Metal.STRAND_CONFIG_DIFFERENT_ALLELES, markerNames[i] + "\t" + Metal.STRAND_CONFIG_DIFFERENT_ALLELES);
						}
						System.err.println("Error - for marker " + markerNames[i] + " " + dirs[index] + " has different alleles (" + line[0] + "/" + line[1] + ") than the rest (" + refAlleles[0] + "/" + refAlleles[1] + ")");
						break;
					case Metal.STRAND_CONFIG_SPECIAL_CASE:
						if (!problematicTrack.containsKey(markerNames[i] + "\t" + Metal.STRAND_CONFIG_SPECIAL_CASE)) {
							writerProblematic.println(markerNames[i] + "\t" + Metal.STRAND_CONFIG_SPECIAL_CASE + "\t" + dirs[index]);
							problematicTrack.put(markerNames[i] + "\t" + Metal.STRAND_CONFIG_SPECIAL_CASE, markerNames[i] + "\t" + Metal.STRAND_CONFIG_SPECIAL_CASE);
						}
	        			//System.err.println("Warning - marker "+markerNames[i]+" has a special case starting with "+dirs[index]+": alleles ("+line[0]+"/"+line[1]+") where previous had only ("+refAlleles[0]+"/"+refAlleles[1]+")");
	                    break;
                    default:
                    	System.err.println("Error - unknown determineStrandConfig return code");
                    	break;
                    }
	        		
		        	gCounts[index] = Array.toStr(genotypeCounts[j], "/");
	        		alleleCounts[j][0] = genotypeCounts[j][0]*2+genotypeCounts[j][1];
	        		alleleCounts[j][1] = genotypeCounts[j][1]+genotypeCounts[j][2]*2;
	        		freqs[index] = (double)(alleleCounts[j][0])/(double)(alleleCounts[j][0]+alleleCounts[j][1]);
                }
	        	if (keys.length > 1) {
		        	writer2.print(markerNames[i]+"\t"+Array.toStr(gCounts));
		        	for (int j = 0; j<freqs.length; j++) {
		        		writer2.print("\t");
		        		if (freqs[j] > Double.MIN_VALUE) {
			        		writer2.print(ext.formDeci(freqs[j], 3, true));
		        		}
	                }
		        	alleleCounts = Matrix.prune(alleleCounts);
		        	if (alleleCounts.length < 2 || alleleCounts[0].length < 2) {
			        	writer2.print("\tX");
		        	} else {
		        		writer2.print("\t"+ProbDist.ChiDist(ContingencyTable.ChiSquare(Matrix.toDoubleArrays(alleleCounts), false, false), (alleleCounts.length-1)*(alleleCounts[0].length-1)));
		        	}
	
		        	genotypeCounts = Matrix.prune(genotypeCounts);
		        	if (genotypeCounts.length < 2 || genotypeCounts[0].length < 2) {
			        	writer2.print("\tX");
			        	p = Double.NaN;
		        	} else {
		        		p = ProbDist.ChiDist(ContingencyTable.ChiSquare(Matrix.toDoubleArrays(genotypeCounts), false, false), (genotypeCounts.length-1)*(genotypeCounts[0].length-1)); 
		        		writer2.print("\t"+p);
		        		if (p<LOOSE_HOMOGENEITY_THRESHOLD) {
			        		writer3.println("#"+markerNames[i]);
			        		writer3.print("alleles <- matrix(c(");
			        		for (int j = 0; j<genotypeCounts.length; j++) {
				        		writer3.print((j==0?"":", ")+Array.toStr(genotypeCounts[j], ", "));
	                        }
			        		writer3.println("), nr="+genotypeCounts[0].length+")");
			        		writer3.println("chisq.test(alleles)");
			        		writer3.println("fisher.test(alleles, simulate.p.value = TRUE, B=1e6)");
		        		}
		        	}
	        		writer2.println();
		        	
		        	if (p < HOMOGENEITY_THRESHOLD) {
		        		writer.println(markerNames[i]);
		        	}
	        	}
            }
	        writer.close();
	        writer2.close();
	        writer3.close();
	        writer3 = new PrintWriter(new FileWriter(outputDir+"split.crf"));
	        writer3.println("split\nhomo.R numFiles=12 sizeOfHeader=0 blockSize=4 root=homo ext=.R");
            writer3.close();
//	        writer3 = new PrintWriter(new FileWriter(dir+"master"));
//	        for (int i = 1; i<=12; i++) {
//		        writer3.println("R CMD BATCH homo"+i+".R &");
//            }
//            writer3.close();
        } catch (Exception e) {
	        System.err.println("Error writing to " + outputDir + CHI_SQUARE_DROPS_FILENAME);
	        e.printStackTrace();
        }
        System.out.println("Splitting up file");
    	Files.splitFile(outputDir+"homo.R", 12, 0, 4, outputDir+"homo", ".R", false);
		Files.qsub(outputDir + "runHomo", null, -1, Rscript.getRExecutable(new Logger()) + " CMD BATCH " + outputDir + "homo[%0].R " + outputDir + "homo[%0].Rout", Matrix.toMatrix(Array.stringArraySequence(12, "")), 1000, 12);
        System.out.println("To run Fisher's exact with permutations in R use either ./master or ./master.runHomo");
        System.out.println(ext.getTime()+"\tDone!");
	}
	
	public static void batchMerge(String dir, String[] subdirs, String commonPlinkRoot) {
        PrintWriter writer;
        String trav, next;

        for (int i = 0; i<subdirs.length; i++) {
        	subdirs[i] = ext.verifyDirFormat(subdirs[i]);
        	if (!new File(dir+subdirs[i]).exists() || !new File(dir+subdirs[i]).isDirectory()) {
        		System.err.println("Error - '"+subdirs[i]+"' is not a valid directory");
        		return;
        	}
        	if (!new File(dir+subdirs[i]+commonPlinkRoot+".bed").exists() || !new File(dir+subdirs[i]+commonPlinkRoot+".bim").exists() || !new File(dir+subdirs[i]+commonPlinkRoot+".fam").exists()) {
        		System.err.println("Error - '"+subdirs[i]+"' does not contain a full set of "+commonPlinkRoot+".bed/bim/fam files");
        		return;
        	}
        	System.out.println(subdirs[i]+" checks out...");
        }
        
        try {
	        writer = new PrintWriter(new FileWriter(dir+"batchMerge"+(Files.isWindows()?".bat":"")));
	        for (int i = 0; i<subdirs.length; i++) {
	        	writer.println("java -cp /home/npankrat/" + common.PSF.Java.GENVISIS + " filesys.SnpMarkerSet file="+subdirs[i]+commonPlinkRoot+".bim");
	        }
	        writer.println();
	        if (Files.isWindows()) {
		        writer.println("start /B batchMerge.unambiguous.bat");
		        writer.println("start /B batchMerge.ambiguous.bat");
	        } else {
		        writer.println("./batchMerge.unambiguous &");
		        writer.println("./batchMerge.ambiguous &");
	        }
	        writer.close();
	        
	        if (!Files.isWindows()) {
		        Files.chmod("batchMerge");
			}
        } catch (Exception e) {
	        System.err.println("Error writing to "+"batchMerge");
	        e.printStackTrace();
        }
        
        try {
	        writer = new PrintWriter(new FileWriter(dir+"batchMerge.unambiguous"+(Files.isWindows()?".bat":"")));
	        
        	writer.println("plink --bfile "+subdirs[0]+"plink --extract "+subdirs[0]+commonPlinkRoot+".bim_unambiguous.txt --make-bed --out "+subdirs[0]+"unambiguous");
        	writer.println();
	        trav = subdirs[0]+"unambiguous";
	        for (int i = 1; i<subdirs.length; i++) {
	        	next = "merged_1-"+(i+1);

	        	writer.println("plink --bfile "+subdirs[i]+commonPlinkRoot+" --extract "+subdirs[i]+commonPlinkRoot+".bim_unambiguous.txt --make-bed --out "+subdirs[i]+"unambiguous");
	        	writer.println("plink --bfile "+trav+" --bmerge "+subdirs[i]+"unambiguous.bed "+subdirs[i]+"unambiguous.bim "+subdirs[i]+"unambiguous.fam --make-bed --out "+next);
	        	writer.println("plink --bfile "+subdirs[i]+"unambiguous --flip "+next+".missnp --make-bed --out "+subdirs[i]+"unambiguousFlipped");
	        	writer.println("plink --bfile "+trav+" --bmerge "+subdirs[i]+"unambiguousFlipped.bed "+subdirs[i]+"unambiguousFlipped.bim "+subdirs[i]+"unambiguousFlipped.fam --make-bed --out "+next);
	        	writer.println();

	        	trav = next;
            }
	        
	        writer.println("mkdir unambiguous/");
	        writer.println("cp "+trav+".bed unambiguous/plink.bed");
	        writer.println("cp "+trav+".bim unambiguous/plink.bim");
	        writer.println("cp "+trav+".fam unambiguous/plink.fam");
        	writer.println();

	        writer.close();
	        
	        if (Files.isWindows()) {
				new File(dir+"unambiguous/");
			} else {
				Files.chmod("batchMerge.unambiguous");
			}
        } catch (Exception e) {
	        System.err.println("Error writing to "+"batchMerge.unambiguous");
	        e.printStackTrace();
        }

        try {
	        writer = new PrintWriter(new FileWriter(dir+"batchMerge.ambiguous"+(Files.isWindows()?".bat":"")));
	        
        	writer.println("plink --bfile "+subdirs[0]+commonPlinkRoot+" --exclude "+subdirs[0]+commonPlinkRoot+".bim_unambiguous.txt --make-bed --out "+subdirs[0]+"ambiguous");
        	writer.println();
	        trav = subdirs[0]+"ambiguous";
	        for (int i = 1; i<subdirs.length; i++) {
	        	next = "mergedAmbiguous_1-"+(i+1);

	        	writer.println("plink --bfile "+subdirs[i]+commonPlinkRoot+" --exclude "+subdirs[i]+commonPlinkRoot+".bim_unambiguous.txt --make-bed --out "+subdirs[i]+"ambiguous");
	        	writer.println("plink --bfile "+trav+" --bmerge "+subdirs[i]+"ambiguous.bed "+subdirs[i]+"ambiguous.bim "+subdirs[i]+"ambiguous.fam --make-bed --out "+next);
	        	writer.println("plink --bfile "+subdirs[i]+"ambiguous --flip "+next+".missnp --make-bed --out "+subdirs[i]+"ambiguousFlipped");
	        	writer.println("plink --bfile "+trav+" --bmerge "+subdirs[i]+"ambiguousFlipped.bed "+subdirs[i]+"ambiguousFlipped.bim "+subdirs[i]+"ambiguousFlipped.fam --make-bed --out "+next);
	        	writer.println();

	        	trav = next;
            }

	        writer.println("mkdir ambiguous/");
	        writer.println("cp "+trav+".bed ambiguous/plink.bed");
	        writer.println("cp "+trav+".bim ambiguous/plink.bim");
	        writer.println("cp "+trav+".fam ambiguous/plink.fam");
        	writer.println();

	        writer.close();
	        
	        if (Files.isWindows()) {
				new File(dir+"ambiguous/");
			} else {
				Files.chmod("batchMerge.ambiguous");
			}
        } catch (Exception e) {
	        System.err.println("Error writing to "+"batchMerge.unambiguous");
	        e.printStackTrace();
        }
	
	}

	public static void mergeMaps(String dir) {
		BufferedReader reader;
        PrintWriter writer, writer2;
        String[] line, files, keys, datasetNames;
        String temp;
        Hashtable<String,Vector<String>> hash = new Hashtable<String,Vector<String>>();
        Vector<String> v, values;
        Vector<Vector<String>> datasets;
        int[] counts, order;
        int index;
        
        ext.timestamp();
        
        files = Files.list(dir, ".bim", false);
        datasetNames = new String[files.length];
        for (int i = 0; i<files.length; i++) {
        	System.out.println("Loading "+files[i]);
        	try {
        		datasetNames[i] = files[i].substring(0, files[i].length()-4);
	            reader = new BufferedReader(new FileReader(dir+files[i]));
	            while (reader.ready()) {
	            	line = reader.readLine().trim().split("[\\s]+");
	            	HashVec.addToHashVec(hash, line[1], i+"\t"+line[0]+"\t"+line[3], false);
	            }
	            reader.close();
            } catch (FileNotFoundException fnfe) {
	            System.err.println("Error: file \""+dir+files[i]+"\" not found in current directory");
	            System.exit(1);
            } catch (IOException ioe) {
	            System.err.println("Error reading file \""+dir+files[i]+"\"");
	            System.exit(2);
            }
        }
        
        keys = HashVec.getKeys(hash, false, false);
        Files.writeList(keys, dir+"list.0.snps");
        try {
	        writer = new PrintWriter(new FileWriter(dir+"unanimousPositions.xln"));
	        writer2 = new PrintWriter(new FileWriter(dir+"discrepancies.xln"));
	        writer.println("SNP\tConsensus?\tMajor Count\tTotal Count\tPercentage\tChr\tPosition\tMajor datasets\tAlt_chr\tAlt_pos\tAlt datasets...");
	        writer2.println("SNP\tConsensus?\tMajor Count\tTotal Count\tPercentage\tChr\tPosition\tMajor datasets\tAlt_chr\tAlt_pos\tAlt datasets...");
//        	ext.timestamp();
//        	System.out.println("Sorting marker names");
        	ext.timestamp();
        	System.out.println("Writing consensus map");
	        for (int i = 0; i<keys.length; i++) {
	        	v = hash.get(keys[i]);
	        	values = new Vector<String>();
	        	datasets = new Vector<Vector<String>>();
	        	for (int j = 0; j<v.size(); j++) {
	        		line = v.elementAt(j).split("[\\s]+");
	        		index = values.indexOf(line[1]+"\t"+line[2]);
	        		if (index == -1) {
	        			values.add(line[1]+"\t"+line[2]);
	        			index = datasets.size(); 
	        			datasets.add(new Vector<String>());
	        		}
	        		datasets.elementAt(index).add(datasetNames[Integer.parseInt(line[0])]);
                }
	        	counts = new int[datasets.size()];
	        	for (int j = 0; j<counts.length; j++) {
	        		counts[j] = datasets.elementAt(j).size();
                }
	        	order = Sort.quicksort(counts, Sort.DESCENDING);
	        	temp = keys[i]+"\t"+(order.length==1?1:0)+"\t"+counts[order[0]]+"\t"+Array.sum(counts)+"\t"+((double)counts[order[0]]/(double)Array.sum(counts));
	        	for (int j = 0; j<order.length; j++) {
	        		temp += "\t"+values.elementAt(order[j])+"\t"+ext.listWithCommas(Array.toStringArray(datasets.elementAt(order[j])));
                }
	        	if (order.length>1) {
	        		writer2.println(temp);
	        	} else {
		        	writer.println(temp);
	        	}
	        }
	        writer.close();
        } catch (Exception e) {
	        System.err.println("Error writing to "+dir+"consensusMap.xln");
	        e.printStackTrace();
        }
        
        ext.timestamp();
	}
	
	public static void parseHomo(String dir) {
		BufferedReader reader;
        PrintWriter writer, writer2;
        String[] files;
        String temp;
        int count, trav;
        long time;
        String[] record;
        
        if (!Files.exists(dir, false)) {
        	System.err.println("Error - directory not found: "+dir);
        	System.err.println("      - using current directory");
        	dir = "./";
        }
        files = Files.list(dir, ".Rout", false);
        System.out.println("Found "+files.length+" files with the extension .Rout");
        
        try {
	        writer = new PrintWriter(new FileWriter(dir+"FisherResults.xln"));
	        writer2 = new PrintWriter(new FileWriter(dir + FISHER_OR_CHI_SQUARE_DROPS_FILENAME));
	        writer.println("Marker\tPearson\tFisher10K\tminPvalue");
            time = new Date().getTime();
	        count = 0;
	        trav = -2;
	        record = null;
	        for (int i = 0; i<files.length; i++) {
	        	try {
	                reader = new BufferedReader(new FileReader(dir+files[i]));
	                while (reader.ready()) {
	                	temp = reader.readLine();
	                	if (temp.startsWith("> #")) {
	                		record = new String[] {temp.substring(3), "error", "error", "error"};
	                		trav = 0;
	                	} else if (temp.indexOf("Pearson's Chi-squared test") > 0) {
	                		trav = 1;
 	                	} else if (temp.indexOf("Fisher's Exact Test") > 0) {
	                		trav = 2;
	                	} else if (temp.indexOf("p-value = ") >= 0 || temp.indexOf("p-value < ") >= 0) {
	                		try {
		                		record[trav] = temp.substring(temp.indexOf("p-value ")+10);
		                		if (trav == 2) {
			                		if (Double.parseDouble(record[1]) < Double.parseDouble(record[2])) {
			                			record[3] = record[1];
			                		} else {
			                			record[3] = record[2];
			                		}
		                			writer.println(Array.toStr(record));
		                			if (Double.parseDouble(record[3]) < HOMOGENEITY_THRESHOLD) {
		                				writer2.println(record[0]);
		                			}
		                			trav = -1;
		                			count++;
		                		}
	                		} catch (Exception e) {
	                			System.err.println("Error - p-values are not being parsed correctly, aborting at marker "+record[0]+" in file "+files[i]);
	                			return;
	                		}
	                	}
	                }
	                reader.close();
                } catch (FileNotFoundException fnfe) {
	                System.err.println("Error: file \""+dir+files[i]+"\" not found in current directory");
	                System.exit(1);
                } catch (IOException ioe) {
	                System.err.println("Error reading file \""+dir+files[i]+"\"");
	                System.exit(2);
                }
	        }
	        writer.close();
	        System.out.println("Found results for "+count+" markers in "+ext.getTimeElapsed(time));
	        System.out.println("There were "+Files.countLines(dir+FISHER_OR_CHI_SQUARE_DROPS_FILENAME, 0)+" markers that had a minimum p-value less than "+HOMOGENEITY_THRESHOLD);
        } catch (Exception e) {
	        System.err.println("Error writing to "+dir+"FisherResults.xln");
	        e.printStackTrace();
        }
	}
	
	public static void updateMap(String fileToUpdate, String mergedMap) {
		BufferedReader reader;
        PrintWriter writer;
        String[] line, loc;
        Hashtable<String,String[]> hash;
        
        if (Files.exists(mergedMap+".ser", false)) {
        	System.out.print("Loading serialized "+mergedMap);
        	hash = SerialHash.loadSerializedStringArrayHash(mergedMap+".ser");
            System.out.println("...done");
        } else {
        	System.out.print("Serializing "+mergedMap);
        	hash = new Hashtable<String,String[]>();
        	try {
	            reader = new BufferedReader(new FileReader(mergedMap));
	            while (reader.ready()) {
	            	line = reader.readLine().trim().split("[\\s]+");
	            	hash.put(line[0], new String[] {line[1], line[2]});
	            }
	            reader.close();
            } catch (FileNotFoundException fnfe) {
	            System.err.println("Error: file \""+mergedMap+"\" not found in current directory");
	            System.exit(1);
            } catch (IOException ioe) {
	            System.err.println("Error reading file \""+mergedMap+"\"");
	            System.exit(2);
            }
            SerialHash.createSerializedStringArrayHash(mergedMap+".ser", hash);
            System.out.println("...done");
        }
        
        new File(fileToUpdate).renameTo(new File(fileToUpdate+"_"));
        try {
	        reader = new BufferedReader(new FileReader(fileToUpdate+"_"));
	        writer = new PrintWriter(new FileWriter(fileToUpdate));
	        while (reader.ready()) {
	        	line = reader.readLine().trim().split("[\\s]+");
	        	loc = hash.get(line[1]);
	        	if (loc == null) {
	        		System.err.println("Error - merged map '"+mergedMap+"' does not contain a position for '"+line[1]+"'; aborting");
	        		writer.close();
	        		new File(fileToUpdate).delete();
	                new File(fileToUpdate+"_").renameTo(new File(fileToUpdate));
	        	}
	        	line[0] = loc[0];
	        	line[3] = loc[1];
	        	writer.println(Array.toStr(line));
	        }
	        reader.close();
            writer.close();
        } catch (FileNotFoundException fnfe) {
	        System.err.println("Error: file \""+fileToUpdate+"_"+"\" not found in current directory");
	        System.exit(1);
        } catch (IOException ioe) {
	        System.err.println("Error reading file \""+fileToUpdate+"_"+"\"");
	        System.exit(2);
        }
	}

	public static void main(String[] args) {
	    int numArgs = args.length;
	    String dir = "";
//	    String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\tWork\\Consortium\\00src\\Miami\\";
//	    String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\tWork\\Consortium\\00src\\Miami\\cnvis\\";
//	    String dir = "D:/tWork/Consortium/00src/";
	    String batch = "";
//	    String batch = "cidr/,miami/,ngrc/";
	    boolean checkHomo = false;
	    boolean consens = false;
	    boolean parseHomo = false;
	    String update = "";
	    String map = "";
	    String plinkRoot = "final";
	    
	    String usage = "\n"+
	    "gwas.MergeDatasets requires 0-1 arguments\n"+
	    "   (1) directory (i.e. dir="+dir+" (default))\n"+
	    "   (2) check for homogeneity among control frequencies (i.e. -checkHomo (not the default))\n"+
	    " OR:\n"+
	    "   (2) parse test of homogeneity results (i.e. -parseHomo (not the default))\n"+
	    " OR:\n"+
	    "   (2) set up batch merge (i.e. batch=dir1/,dir2/,lastDir/ (not the default))\n"+
	    "   (3) root of plink files (i.e. root="+plinkRoot+" (default))\n"+
	    " OR:\n"+
	    "   (2) update indiviudal map with mergedMap (i.e. update=plink.bim (not the default))\n"+
	    " OR:\n"+
	    "   (2) mergedMap filename (i.e. map=allSNPs.xln (not the default))\n"+
	    "";

	    for (int i = 0; i<args.length; i++) {
		    if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
			    System.err.println(usage);
			    System.exit(1);
		    } else if (args[i].startsWith("dir=")) {
			    dir = args[i].split("=")[1];
			    numArgs--;
		    } else if (args[i].startsWith("batch=")) {
			    batch = args[i].split("=")[1];
			    numArgs--;
		    } else if (args[i].startsWith("-checkHomo")) {
			    checkHomo = true;
			    numArgs--;
		    } else if (args[i].startsWith("-parseHomo")) {
			    parseHomo = true;
			    numArgs--;
		    } else if (args[i].startsWith("-findConsensusPositions")) {
			    consens = true;
			    numArgs--;
		    } else if (args[i].startsWith("update=")) {
			    update = args[i].split("=")[1];
			    numArgs--;
		    } else if (args[i].startsWith("map=")) {
			    map = args[i].split("=")[1];
			    numArgs--;
		    }
	    }
	    if (numArgs!=0) {
		    System.err.println(usage);
		    System.exit(1);
	    }
	    
//	    checkHomo = true;
//	    parseHomo = true;
	    if (!Files.exists(dir, false)) {
	    	System.err.println("Error - using current directory instead of the one that does not exist: "+dir);
	    	dir = "./";
	    }
	    
//	    parseHomo = true;
//	    dir = "D:/data/aric_affy/ancestry26K/homogeneityWithPDs/";
	    
		try {
	    	if (checkHomo) {
	    		checkForHomogeneity(dir);
	    	} else if (parseHomo) {
	    		parseHomo(dir);
	    	} else if (!update.equals("")) {
	    		updateMap(update, map);
	    	} else if (!batch.equals("")) {
	    		batchMerge(dir, batch.split(","), plinkRoot);
	    	} else if (consens) {
	    		mergeMaps(dir);
	    	}
	    } catch (Exception e) {
		    e.printStackTrace();
	    }
    }
}
