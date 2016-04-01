// -Xms1024M -Xmx1024M
package cnv.analysis;

import java.io.*;
import java.util.*;
import java.util.concurrent.Callable;

import javax.swing.JOptionPane;

import cnv.analysis.MosaicismDetect.MosaicBuilder;
import cnv.filesys.*;
import cnv.filesys.MarkerSet.PreparedMarkerSet;
import cnv.filesys.Project.ARRAY;
import cnv.plots.MosaicPlot;
import cnv.var.CNVariant;
import cnv.var.IndiPheno;
import cnv.var.LocusSet;
import cnv.var.MosaicRegion;
import cnv.var.SampleData;
import common.*;
import common.WorkerTrain.Producer;
import filesys.Segment;

public class Mosaicism {
	public static final String[] HEADER = {"Sample", "Arm", "#CNVs", "Summed_Size", "%covered", "Custom_metric", "LRR_SD", "LRR_SD_flag", "Flagged_via_Mosaicism", "Mosaicism_level", "Mosaicism description"};
	public static final double LOWER_BOUND = 0.15;
	public static final double UPPER_BOUND = 0.85;

	public static void findOutliers(Project proj) {
		findOutliers(proj, -1);
	}

	public static void findOutliers(Project proj, int numthreads) {
		PrintWriter writer;
		String[] samples;
		int chr;
		Hashtable<String,String> hash = new Hashtable<String,String>();
		String[] markerNames;
		byte[] chrs;
		int[] positions;
		boolean[] snpDropped;
		int[][] chrBoundaries;
		MarkerSet markerSet;

		hash = proj.getFilteredHash();

		chrBoundaries = new int[27][3];
		snpDropped = null;
		for (int i = 0; i<chrBoundaries.length; i++) {
			chrBoundaries[i][0] = chrBoundaries[i][1] = chrBoundaries[i][2] = -1;
		}
		markerSet = proj.getMarkerSet();
		markerNames = markerSet.getMarkerNames();
		chrs = markerSet.getChrs();
		positions = markerSet.getPositions();
		snpDropped = new boolean[markerNames.length];
		int[][] indicesByChr= markerSet.getIndicesByChr();
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
		
        if (Thread.currentThread().isInterrupted()) { throw new RuntimeException(new InterruptedException()); }
		samples = proj.getSamples();
		try {
			writer = new PrintWriter(new FileWriter(proj.MOSAIC_RESULTS_FILENAME.getValue()));
			writer.println(Array.toStr(MosaicPlot.MOSAICISM_HEADER));
//			samples = new String[] { "7355066051_R03C01", "7330686030_R02C01", "7159911135_R01C02" };
//			samples = new String[] { "7355066051_R03C01" };

			MosaicResultProducer producer = new MosaicResultProducer(proj, samples, snpDropped, chrBoundaries, markerSet, indicesByChr);
			WorkerTrain<String[]> train = new WorkerTrain<String[]>(producer, numthreads > 0 ? numthreads : proj.NUM_THREADS.getValue(), 2, proj.getLog());
			int index =0;
			long timePer = System.currentTimeMillis();
			long time = System.currentTimeMillis();

			while (train.hasNext()) {
		        if (Thread.currentThread().isInterrupted()) { writer.close(); throw new RuntimeException(new InterruptedException()); }
				try {
					String[] results = train.next();
					index++;
					if (index % numthreads == 0) {
						proj.getLog().reportTimeInfo((index) + " of " + samples.length + " in " + ext.getTimeElapsed(timePer) + ", total time at " + ext.getTimeElapsed(time));
						timePer = System.currentTimeMillis();
					}

					for (int i = 0; i < results.length; i++) {
						writer.println(results[i]);
					}
				} catch (Exception e) {
					proj.getLog().reportException(e);
					producer.shutdown();
				}

			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to summary file");
			e.printStackTrace();
		}
	}
	
	private static class MosaicResultProducer implements Producer<String[]>{
		private Project proj;
		private String[] samples;
		private boolean[]snpDropped;
		private int[][] chrBoundaries;
		private PreparedMarkerSet markerSet;
		private int[][] indicesByChr;
		private int index;

		public MosaicResultProducer(Project proj, String[] samples, boolean[] snpDropped, int[][] chrBoundaries, MarkerSet markerSet, int[][] indicesByChr) {
			super();
			this.proj = proj;
			this.samples = samples;
			this.snpDropped = snpDropped;
			this.chrBoundaries = chrBoundaries;
			this.indicesByChr = indicesByChr;
			this.markerSet = PreparedMarkerSet.getPreparedMarkerSet(markerSet);
			this.index = 0;
		}

		@Override
		public boolean hasNext() {
			return index < samples.length;
		}

		@Override
		public Callable<String[]> next() {
			final String currentSample =samples[index];
			Callable<String[]> callable = new Callable<String[]>() {
				@Override
				public String[] call() throws Exception {

					return getMosaicResults(proj, currentSample,  snpDropped, chrBoundaries, markerSet, indicesByChr);
				}
			};
			index++;
			return callable;
		}

		@Override
		public void shutdown() {
			index = samples.length;
		}
		
		
	}

	private static String[] getMosaicResults(Project proj, String sample, boolean[] snpDropped, int[][] chrBoundaries, MarkerSet markerSet, int[][] indicesByChr) {
		Sample samp;
		float baf;
		float[] lrrs;
		float[] bafs;
		samp = proj.getPartialSampleFromRandomAccessFile(sample);
		ArrayList<String> results = new ArrayList<String>();
		if (samp.getFingerprint() != markerSet.getFingerprint()) {
			String error = "Error - cannot estimate mosaics if MarkerSet and Sample (" + sample + ") don't use the same markers";
			throw new IllegalArgumentException(error);
		}
		lrrs = samp.getLRRs();
		bafs = samp.getBAFs();
		MosaicBuilder builder = new MosaicBuilder();
		builder.indicesByChr(indicesByChr);
		builder.verbose(false);
		builder.markerIndices(proj.getMarkerIndices());
		if (proj.getArrayType() == ARRAY.NGS) {
			proj.getLog().reportTimeWarning("Masking non-variant sites for project type " + ARRAY.NGS);
			boolean[] use = new boolean[markerSet.getMarkerNames().length];
			int numMasked = 0;
			for (int i = 0; i < use.length; i++) {
				boolean useit = !proj.getArrayType().isCNOnly(markerSet.getMarkerNames()[i]);
				if (!useit) {
					numMasked++;
				}
				use[i] = useit;
			}
			proj.getLog().reportTimeInfo(numMasked + " markers were masked");
			builder.use(use);
		}
		MosaicismDetect md = builder.build(proj, sample, markerSet, Array.toDoubleArray(bafs));
		int[] positions = markerSet.getPositions();
		byte[]chrs = markerSet.getChrs();
		for (int j = 1; j<=23; j++) {
			for (int arm = 0; arm<2; arm++) {
				int startIndex =  (arm==0?chrBoundaries[j][0]:chrBoundaries[j][1]);
				int stopIndex = (arm == 0 ? chrBoundaries[j][1] : chrBoundaries[j][2] + 1);
				

				ArrayList<Float> lrrAl = new ArrayList<Float>(stopIndex+10-startIndex);
				ArrayList<Float> bafAl = new ArrayList<Float>(stopIndex+10-startIndex);
				for (int k = startIndex; k < stopIndex; k++) {
					if (!snpDropped[k] && !(md.getUse() != null && !md.getUse()[k])) {
						if (!Float.isNaN(lrrs[k])) {
							lrrAl.add(lrrs[k]);
						}
						baf = bafs[k];
						if (baf>LOWER_BOUND&&baf<UPPER_BOUND) {
							bafAl.add(baf);
						}
					}
				}
				if (lrrAl.size()>100) {
					if (chrs[startIndex] != (byte) j || chrs[stopIndex - 1] != (byte) j) {

						throw new IllegalStateException("Internal Error, mismatched chromosome indices, start =" + chrs[startIndex] + "\tstop = " + chrs[stopIndex] + "\tarm = " + arm);
					}
					Segment armSeg = new Segment((byte) j, positions[startIndex], positions[stopIndex - 1]);
					MosaicMetric mosaicMetrics = getMosiacMetric(md, armSeg, proj.getLog());

					int bafSize = bafAl.size();
					int lrrSize = lrrAl.size();
					float[] bafTmp = Array.toFloatArray(bafAl);
					String result = sample + "\t" + "chr" + j + (arm == 0 ? "p" : "q") + "\t" + lrrSize + "\t" + ext.formDeci(Array.mean(Array.toFloatArray(lrrAl)), 5) + "\t" + bafAl.size() + (bafSize > 10 ? "\t" + ext.formDeci(Array.stdev(bafTmp, true), 5) + "\t" + ext.formDeci(Array.iqr(bafTmp), 5) : "\t.\t.") + "\t" + ext.formDeci((double) (lrrSize - bafSize) / (double) lrrSize, 5) + "\t" + mosaicMetrics.getPercentBandMosaic() + "\t" + mosaicMetrics.getBpWeightedAverage() + "\t" + mosaicMetrics.getMosaicRegionsDetected();
					results.add(result);
				}
			}
		}
		return Array.toStringArray(results);
	}

	private static class MosaicMetric {
		private double percentBandMosaic;
		private double bpWeightedAverage;
		private int mosaicRegionsDetected;

		public MosaicMetric(double percentBandMosaic, double mosaiceMetric,int mosaicRegionsDetected) {
			super();
			this.percentBandMosaic = percentBandMosaic;
			this.bpWeightedAverage = mosaiceMetric;
		}

		private int getMosaicRegionsDetected() {
			return mosaicRegionsDetected;
		}

		private void setMosaicRegionsDetected(int mosaicRegionsDetected) {
			this.mosaicRegionsDetected = mosaicRegionsDetected;
		}

		private double getPercentBandMosaic() {
			return percentBandMosaic;
		}

		private void setPercentBandMosaic(double percentBandMosaic) {
			this.percentBandMosaic = percentBandMosaic;
		}

		private double getBpWeightedAverage() {
			return bpWeightedAverage;
		}

		private void setBpWeightedAverage(double bpWeightedAverage) {
			this.bpWeightedAverage = bpWeightedAverage;
		}

	}

	private static MosaicMetric getMosiacMetric(MosaicismDetect md, Segment seg, Logger log) {
		LocusSet<MosaicRegion> mosSet = md.callMosaic(seg, true);
		MosaicMetric mosaicMetric = new MosaicMetric(-1, -1, -1);
		if (mosSet.getLoci().length != 1 || !seg.equals(mosSet.getLoci()[0])) {
			log.reportTimeError("Mosaic caller not in force call mode");
			log.reportTimeError(seg.getUCSClocation() + " went in, and " + mosSet.getLoci()[0].getUCSClocation() + " came out");
		} else if (seg.getChr() < 23) {// can't call chr23 yet
			mosaicMetric.setPercentBandMosaic(100 * mosSet.getLoci()[0].getScore());
		}

		LocusSet<MosaicRegion> tmp = md.callMosaic(seg, false);
		mosaicMetric.setMosaicRegionsDetected(tmp.getLoci().length);
		if (tmp.getLoci().length < 1 || seg.getChr() >= 23) {

		} else {
			double bpWeightedAverage = 0;
			int numRegions =0;
			for (int i = 0; i < tmp.getLoci().length; i++) {
				if (tmp.getLoci()[i].getNumMarkers() > 2 * md.getMovingFactor()) {
					bpWeightedAverage += (double) tmp.getLoci()[i].getSize() * tmp.getLoci()[i].getScore();
					numRegions++;
				}
			}
			bpWeightedAverage = (double) bpWeightedAverage / seg.getSize();
			mosaicMetric.setMosaicRegionsDetected(numRegions);
			mosaicMetric.setBpWeightedAverage(bpWeightedAverage);
		}
		return mosaicMetric;
	}

	public static void checkForOverlap(Project proj, String listOfMosaicArms) {
		BufferedReader reader;
        PrintWriter writer;
        String[] line;
        Vector<String> v;
        Vector<String[]> list;
        int count;
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
        cnvFiles = proj.CNV_FILENAMES.getValue();
        if (cnvFiles.length == 0) {
        	System.err.println("Error - need to specify the name of a CNV file in the project properties file before running Mosaicism.checkForOverlap()");
        	return;
        }
        sampleData = proj.getSampleData(2, new String[] {cnvFiles[0]});
        if (Files.exists(proj.PROJECT_DIRECTORY.getValue()+"lrr_sd.xln", proj.JAR_STATUS.getValue())) {
        	lrrsdHash = HashVec.loadFileToHashString(proj.PROJECT_DIRECTORY.getValue()+"lrr_sd.xln", false);
        } else {
        	System.err.println("Warning - could not find 'lrr_sd.xln' in project directory; no flags will be generated");
        	lrrsdHash = new Hashtable<String,String>();
        }
        if (Files.exists(proj.MOSAIC_COLOR_CODES_FILENAME.getValue(false, false), proj.JAR_STATUS.getValue())) {
        	mosaicHash = HashVec.loadFileToHashString(proj.MOSAIC_COLOR_CODES_FILENAME.getValue(false, false), new int[] {0,1}, new int[] {2,3}, false, "\t", true, proj.JAR_STATUS.getValue(), true);
        } else {
        	System.err.println("Warning - could not find "+proj.MOSAIC_COLOR_CODES_FILENAME.getValue(false, false)+"; no annotation possible");
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
		        	reader.close();
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
		        		arm = new Segment((byte)chr, Positions.CENTROMERE_MIDPOINTS[chr], Positions.CHROMOSOME_LENGTHS_B36_HG18[chr]);
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
                writer = new PrintWriter(new FileWriter(proj.PROJECT_DIRECTORY.getValue()+"SAMPLES_IN_CNVFILE_NOT_IN_SAMPLE_DATA.txt"));
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
		String filename = null;
		String arms = "MosaicArms.txt";
		Project proj;
		boolean check = false;
		int numthreads = 24;

		String usage = "\n"+
				"filesys.ParseIllumina requires 0-1 arguments\n" +
				"   (1) project properties filename (i.e. proj=" + cnv.Launch.getDefaultDebugProjectFile(false) + " (default))\n" +
				"   (2) check for overlap between mosaic arms and CNV calls in the first CNV file listed in the project file (i.e. -check (not the default))\n" +
				"   (3) mosaic arms file (i.e. arms=MosaicArms.txt (default))\n" +
				PSF.Ext.getNumThreadsCommand(4, numthreads) +

				"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
		        return;
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("arms=")) {
			    arms = args[i].split("=")[1];
			    numArgs--;
			} else if (args[i].startsWith("-check")) {
				check = true;
				numArgs--;
			}else if (args[i].startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
				numthreads = ext.parseIntArg(args[i]);
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
				checkForOverlap(proj, arms);
			} else {
				findOutliers(proj);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
