package affy;


import java.io.FileWriter;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.Arrays;
import java.util.Date;


import cnv.filesys.Centroids;
import cnv.filesys.MarkerData;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import cnv.filesys.Sample;
import cnv.filesys.SampleList;
import cnv.manage.MarkerDataLoader;
import stats.Maths;
import common.Array;
import common.Files;
import common.Sort;
import common.ext;
public class AffyCentroids implements Serializable {
	public static final long serialVersionUID = 1L;
	public static final String[] AFFY_CENTROID_SUFFIXES = {"Name", "AA T Mean", "AA R Mean", "AB T Mean", "AB R Mean", "BB T Mean", "BB R Mean"};

	private float[][][] AffyCentroids;  // marker, genotype (0=AA, 1=AB, 2=BB), coordinates (0=Mean Theta, 1=Mean R) (a.k.a. follows the suffix order above)
	private long fingerprint;
	private static int starter =0;
	private static int stopper = 1855448;

	public AffyCentroids(float[][][] AffyCentroids, long fingerprint) {
		this.AffyCentroids = AffyCentroids;
		this.fingerprint = fingerprint;
	}

	public float[][][] getCentroids() {
		return AffyCentroids;
	}

	public long getFingerprint() {
		return fingerprint;
	}

	public void serialize(String filename) {
		Files.writeSerial(this, filename);
	}

	public static AffyCentroids load(String filename, boolean jar) {
		return (AffyCentroids)Files.readSerial(filename, jar, true);
	}

	public static float calcR(float x, float y) {
		return (float)(Math.max(x, 0.0001)+Math.max(y, 0.0001));
	}

	public static float calcTheta(float x, float y) {
		return (float)(Math.atan2(Math.max(y, 0.0001),Math.max(x, 0.0001))/(Math.PI/2));
	}

	public static float calcBAF(float theta, float[][] AffyCentroids) {
		if (AffyCentroids[0] != null && theta < AffyCentroids[0][0]) {
			return 0;
		} else if (AffyCentroids[1] != null && theta < AffyCentroids[1][0]) {
			if (AffyCentroids[0] == null) {
				return 0.50f;
			} else {
				return 0.5f*(theta-AffyCentroids[0][0])/(AffyCentroids[1][0]-AffyCentroids[0][0]);
			}
		} else if (AffyCentroids[2] != null && theta < AffyCentroids[2][0]) {
			if (AffyCentroids[1] == null) {
				return 1.0f;
			} else {
				return 0.5f+0.5f*(theta-AffyCentroids[1][0])/(AffyCentroids[2][0]-AffyCentroids[1][0]);
			}
		} else {
			if (AffyCentroids[2] == null) {
				return 0.50f;
			} else {
				return 1;
			}
		}
	}

	public static float calcLRR(float theta, float r, float[][] AffyCentroids) {
		float estimatedR;
		if(AffyCentroids[0] == null){
			estimatedR = Float.NaN;			
		}
		else if(theta < AffyCentroids[1][0]){
			if(AffyCentroids[1][0] -AffyCentroids[0][0] !=0){
				estimatedR = AffyCentroids[0][1]+(theta-AffyCentroids[0][0])*(AffyCentroids[1][1]-AffyCentroids[0][1])/(AffyCentroids[1][0]-AffyCentroids[0][0]);
			}
			else{
				estimatedR = AffyCentroids[0][1];
			}
			if(estimatedR< 0){
				estimatedR = r;
				System.out.println("Warning - estimatedR < 0 ");
			}
		}
		else{
			if(AffyCentroids[2][0] -AffyCentroids[1][0] !=0){
				estimatedR = AffyCentroids[1][1]+(theta-AffyCentroids[1][0])*(AffyCentroids[2][1]-AffyCentroids[1][1])/(AffyCentroids[2][0]-AffyCentroids[1][0]);
			}
			else{
				estimatedR = AffyCentroids[1][1];
			}
			if(estimatedR< 0){
				estimatedR = r;
				System.out.println("Warning - estimatedR < 0 ");
			}
		}
		return (float)Maths.log2(r/estimatedR);

	}

	public static float[] getAFFYBAF( String[] markerNames, float[][][] affyCents, float[] Xs , float[] Ys , int i){
		float[] AFFYBAFs;
		AFFYBAFs =new float[markerNames.length];
		for (int k = starter; k<stopper; k++) {
			if(markerNames[k].startsWith("CN_")){
				AFFYBAFs[k] = 0;
			}
			else if (markerNames[k].startsWith("SNP_")){
				AFFYBAFs[k] = calcBAF(calcTheta(Xs[k],Ys[k]), affyCents[k]);
			}
		}
		return AFFYBAFs;	
	}


	public static float[] getAFFYLRR( String[] markerNames, float[][][] affyCents, float[] Xs , float[] Ys, int i){
		float[] AFFYLRRs;
		AFFYLRRs =new float[Xs.length];
		for (int k = starter; k<stopper; k++) {
			if(markerNames[k].startsWith("CN_")){
				//if a copy number probeset, take the log2 intensity value - median log2 intensity
				AFFYLRRs[k] = (float) ( Maths.log2(Xs[k]) - affyCents[k][0][1]);
			}
			else if (markerNames[k].startsWith("SNP_")){
				AFFYLRRs[k] = calcLRR(calcTheta(Xs[k],Ys[k]), calcR(Xs[k],Ys[k]),affyCents[k] );
			}
		}
		return AFFYLRRs;
	}

	public static void recompute(Project proj, String centroidsFile) {
		MarkerSet markerSet;
		Centroids affyCentroids;
		Sample original, sample;
		String[] samples, markerNames;
		float[][][] affyCents;
		float[] Xs;
		float[] Ys; 
		float[] AFFYBAFs;
		float[] AFFYLRRs;
		markerSet = proj.getMarkerSet();
		markerNames = markerSet.getMarkerNames();
		affyCentroids = Centroids.load(centroidsFile, proj.getJarStatus());
		if (affyCentroids.getFingerprint() != markerSet.getFingerprint()) {
			System.err.println("Error - fingerprint for Centroids file '"+centroidsFile+"' does not match the fingerprint for the current MarkerSet");
		}
		affyCents = affyCentroids.getCentroids(); 
		samples = proj.getSamples();
		for (int i = 0; i<samples.length; i++) {
			System.out.println(samples[i]);
			original = proj.getFullSampleFromRandomAccessFile(samples[i]);
			Xs = original.getXs();
			Ys = original.getYs();
			AFFYBAFs = getAFFYBAF(markerNames , affyCents , Xs ,Ys ,i );
			AFFYLRRs = getAFFYLRR(markerNames , affyCents , Xs ,Ys ,i );
			sample = new Sample(original.getSampleName(), original.getFingerprint(), original.getGCs(), original.getXs(), original.getYs(), AFFYBAFs, AFFYLRRs, original.getForwardGenotypes(), original.getAB_Genotypes(), original.getCanXYBeNegative());
			sample.saveToRandomAccessFile(proj.getDir(Project.SAMPLE_DIRECTORY) + original.getSampleName() + Sample.SAMPLE_DATA_FILE_EXTENSION);
		}
	}

	public static void parseCentroids(Project proj, boolean[] samplesToBeUsed, double missingnessThreshold) {
		String[] samples, markerNames;
		float[][][] centroids;
		int count;
		SampleList sampleList;
		MarkerDataLoader markerDataLoader;
		MarkerData markerData;
		long time;
		MarkerSet markerSet;
		time = new Date().getTime();
		System.out.println("Computing centroids from intensity means");
		sampleList = proj.getSampleList();
		samples = sampleList.getSamples();
		if (samples.length != samplesToBeUsed.length) {
			System.err.println("Error - mismatched number of samples in project versus sample mask");
			return;
		}
		markerSet = proj.getMarkerSet();
		markerNames = markerSet.getMarkerNames();
		centroids = new float[markerNames.length][][];
		markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, Arrays.copyOfRange(markerNames,starter ,stopper));

		time = new Date().getTime();
		for (int i = 0; i < stopper-starter; i++) {
			markerData = markerDataLoader.requestMarkerData(i);
			centroids[i] = computeCluster(proj ,samplesToBeUsed, missingnessThreshold, samples,  markerData, i, 0.99);
			markerDataLoader.releaseIndex(i);
		}
		count = 0;
		for (int i = 0; i<centroids.length; i++) {
			if (centroids[i] == null) {
				if (count == 0) {
					System.out.println("The following marker(s) could not be computed:");
				}
				count++;
			}
		}
		if (count > 0) {
			System.out.println("Computed mean genotyped centroids for "+(centroids.length-count)+" of "+centroids.length+" markers, "+count+" missing");
		} else {
			System.out.println("Computed mean genotyped centroids for all "+centroids.length+" markers");
		}
		new Centroids(centroids, markerSet.getFingerprint()).serialize(proj.getFilename(Project.GENOTYPE_CENTROIDS_FILENAME));
		System.out.println("Computation took "+ext.getTimeElapsed(time));
	}

	private static float[][] computeCluster(Project proj, boolean[] samplesToBeUsed, double missingnessThreshold, String[] samples,  MarkerData markerData, int i ,double confThreshold) {
		float[][] centroid;
		float[] thetas;
		float[] rs;
		float[] confs;
		double[] meanThetas;
		double[] meanRs;
		byte[] genotypes;
		int[] counts;
		PrintWriter writer =null;
		String markerName;
		markerName = markerData.getMarkerName();
		meanThetas = new double[5];
		meanRs = new double[5];
		counts = new int[5];
		boolean use = true;

		if(markerName.startsWith("SNP_")){
			for (int k = 0; k<samples.length; k++) {
				//use if X chromosome and sex is female, dont use if X chromsome and sex is male
				//use if Ychromosome and sex is male, dont use if Ychromsome and sex is female
				use = checkSexMarker(proj, samples[k], markerData);
				if(!use || !samplesToBeUsed[k] ){
					continue;
				}
				confs =markerData.getGCs();
				genotypes = markerData.getAB_Genotypes();
				thetas = markerData.getThetas();
				rs = markerData.getRs();
				//maybe alter the confidence checks
				if ( !Float.isNaN(thetas[k]) && !Float.isNaN(rs[k])&& !Float.isNaN(confs[k])&& confs[k]>=confThreshold) {
					meanThetas[0] += thetas[k];
					meanRs[0] += rs[k];
					counts[0]++;
					meanThetas[genotypes[k]+2] += thetas[k];
					meanRs[genotypes[k]+2] += rs[k];
					counts[genotypes[k]+2]++;
				}
			}
			for (int k = 0; k<5; k++) {
				meanThetas[k] /= counts[k];
				meanRs[k] /= counts[k];
			}
		}
		else if(markerName.startsWith("CN_")){
			float[] log2Xs = log2Array(proj, markerData, samplesToBeUsed,samples);
			meanRs[2] = median(log2Xs);
			meanRs[3] = Array.mean(log2Xs );
			meanThetas[2] =0.5;
			meanThetas[3] =0.5;
			counts[2]++;
			counts[3]++;
			counts[0]++;
		}
		centroid = new float[3][];
		if (counts[1] >= counts[0]*missingnessThreshold) {
			for (int k = 0; k<3; k++) {
				centroid[k] = new float[] { (float)meanThetas[0], (float)meanRs[0] };  	
			}
		} else {
			processCentroid(centroid, meanThetas, meanRs, counts, writer, markerName);
		}
		return centroid;
	}

	private static boolean checkSexMarker(Project proj, String sample, MarkerData markerData) {
		int sex =0;
		byte chr;
		boolean use =true;
		sex = proj.getSampleData(0, false).getSexForIndividual(sample);
		chr = markerData.getChr();
		if((int)chr == 23 && sex != 2){
			use =false;
		}
		if((int)chr == 24 && sex != 1){
			use =false;
		}
		return use;
	}

	private static void processCentroid(float[][] centroid, double[] meanThetas, double[] meanRs, int[] counts, PrintWriter writer, String markerName) {
		if(markerName.startsWith("CN_")){
			centroid[0] = new float[] { (float)meanThetas[2], (float)meanRs[2] };
			centroid[1] = new float[] { (float)meanThetas[3], (float)meanRs[3] };
			centroid[2] = new float[] { (float)meanThetas[2], (float)meanRs[2] };

		}
		else if(counts[2] ==0 && counts[3] == 0 && counts[4] == 0){
			nullify(centroid);
		}
		else if (counts[2] ==0 && counts[3] == 0){
			nullify(centroid);
		}
		else if (counts[2] ==0 && counts[4] == 0){
			nullify(centroid);
		}

		else if (counts[3] ==0 && counts[4] == 0){
			nullify(centroid);
		}
		else if(counts[2] ==0){	
			centroid[0] = new float[] { (float)(meanThetas[3] -0.3), getAltRs(centroid ,meanRs , counts , 2 , 3 ,4 ) };
			centroid[1] = new float[] { (float)meanThetas[3], (float)meanRs[3] };
			centroid[2] = new float[] { (float)meanThetas[4], (float)meanRs[4] };
			//System.out.println("Estimating the AA cluster for  " + markerName);
		}
		else if(counts[3] ==0){
			centroid[0] = new float[] { (float)meanThetas[2], (float)meanRs[2] };
			if(counts[2]  >0 && counts[4] >0){
				centroid[1] = new float[] { (float)((meanThetas[2]+meanThetas[4])/2), getAltRs(centroid ,meanRs , counts , 3, 2 ,4 ) };
				//System.out.println("Estimating the AB cluster for  " + markerName);
			}
			else{
				centroid[1] = new float[] { (float)(0.5), getAltRs(centroid ,meanRs , counts , 3, 2 ,4 ) };
				//System.out.println("Estimating the AB cluster for  " + markerName);
			}
			centroid[2] = new float[] { (float)meanThetas[4], (float)meanRs[4] };
		}
		else if(counts[4] ==0){
			centroid[0] = new float[] { (float)meanThetas[2], (float)meanRs[2] };
			centroid[1] = new float[] { (float)meanThetas[3], (float)meanRs[3] };
			centroid[2] = new float[] { (float)(meanThetas[3] +0.3), getAltRs(centroid ,meanRs , counts , 4 , 3 ,2 ) };
			//System.out.println("Estimating the BB cluster for  " + markerName);
		}

		else{	
			for (int k = 0; k<3; k++) {
				if (counts[k+2] > 0) {
					//writer.println(markerName+ "\t"+(float)meanThetas[k+2] +"\t" +(float)meanRs[k+2] );
					centroid[k] = new float[] { (float)meanThetas[k+2], (float)meanRs[k+2] };
				} else {
					System.out.println("this should not happen , but does it?");
				}
			}
		}
		if(centroid[2]!=null && ( centroid[0][0] > centroid[1][0] ||centroid[2][0] <centroid[1][0]) ){
			nullify(centroid);
		}
	}



	private static float getAltRs(float[][] centroid, double[] meanRs, int[] counts ,int checkIndex , int primaryAlt, int secondaryAlt){
		//AssignAltRs(centroid ,meanRs , counts , 2 , 3 ,4 );
		float altR;
		if ( counts[primaryAlt] > 0){
			altR =  (float)meanRs[primaryAlt];
		}
		else if(counts[secondaryAlt] > 0){
			altR = (float)meanRs[secondaryAlt];
		}
		else {
			altR =0;
			System.err.println("This should not happen when assigning altRs");
		}
		return altR;
	}

	private static void nullify(float[][] centroid) {
		centroid[0] = null;
		centroid[1] = null;
		centroid[2] = null;

	}

	public static float median(float[] array) {
		return (quant(array, (float)0.50));
	}

	public static float quant(float[] array, float q) {
		int keys[] = Sort.quicksort(array);
		try {
			if (q>1||q<0) {
				return (0);
			} else {
				double index = (array.length+1)*q;
				if (index-(int)index==0) {
					return array[keys[(int)index-1]];
				} else {
					return q*array[keys[(int)Math.floor(index)-1]]+(1-q)*array[keys[(int)Math.ceil(index)-1]];
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
			return -1234567890;
		}
	}

	//only get samps to be used
	public static float[] log2Array(Project proj ,MarkerData markerData ,boolean[] samplesToBeUsed, String[] samples){
		float[] xs = markerData.getXs();
		float[] log2Xs;
		int count = 0;
		boolean use;
		for (int k = 0; k<samplesToBeUsed.length; k++) { 
			if (samplesToBeUsed[k]){
				count++;
			}
		}
		//System.out.println(count);
		log2Xs = new float[count] ;
		for (int k = 0; k<xs.length; k++) {
			if(samplesToBeUsed[k]){
				use = checkSexMarker(proj, samples[k], markerData);
				if(!use){
					continue;
				}
				log2Xs[k] =(float)Maths.log2(xs[k]);
			}
		}
		return log2Xs;
	}

	public static void exportToText(Project proj, String centFilename, String exportFilename) {
		PrintWriter writer;
		Centroids centObject;
		float[][][] centroids;
		MarkerSet markerSet;
		String[] markerNames;
		String dir;

		dir = proj.getProjectDir();
		markerSet = proj.getMarkerSet();
		markerNames = markerSet.getMarkerNames();
		centObject = Centroids.load(dir+centFilename, false);
		centroids = centObject.getCentroids();

		if (markerNames.length != centroids.length) {
			System.err.println("Error - mismatched number of markers in the project's marker set and the imported centroids file ("+centFilename+"); aborting");
			return;
		}

		if (markerSet.getFingerprint() != centObject.getFingerprint()) {
			System.err.println("Error - mismatched marker fingerprints in the project's marker set and the imported centroids file ("+centFilename+"); aborting");
			return;
		}
		try {
			writer = new PrintWriter(new FileWriter(dir+exportFilename));
			writer.println("marker_fingerprint="+centObject.getFingerprint());
			writer.println("MarkerName\tAA_Theta_Mean\tAA_R_Mean\tAB_Theta_Mean\tAB_R_Mean\tBB_Theta_Mean\tBB_R_Mean");
			for (int i = 0; i < markerNames.length; i++) {
				writer.print(markerNames[i]);
				for (int j = 0; j < 3; j++) {
					if (centroids[i][j] == null) {
						writer.print("\t.\t.");
					} else {
						writer.print("\t"+centroids[i][j][0]+"\t"+centroids[i][j][1]);
					}
				}
				writer.println();
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + exportFilename);
			e.printStackTrace();
		}
	}


	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "C:/workspace/Genvisis/projects/tableTest.properties";
		//		String centFile = "ForNathan_table.csv";
		//		String centFile = "SNP Table Myers raw dataset final 022908.csv";
		//		String centFile = "Myers_final_042208_ReclusteredCNV_SNP_Table2.csv";
		//		String centFile = "Theta_R_mean_dev_550.csv";
		String centFile = "CentroidExample.csv";
		String intensityFlags = "";
		//		String intensityFlags = "intesityFlag.txt";
		//		String clusteredCentroids = Project.DEFAULT_ORIGINAL_CENTROIDS_FILENAME;
		//		String unclusteredCentroids = Project.DEFAULT_GENOTYPE_CENTROIDS_FILENAME;
		boolean fromGenotypes = false;
		Project proj;
		String compute = "";
		String importFile = null;
		String exportFile = null;

		String usage = "\n"+
				"cnv.filesys.Centroids requires 0-1 arguments\n"+
				"   (1) project (i.e. proj=" + filename + " (default))\n"+
				"   (2) filename (i.e. file=" + centFile + " (default))\n"+
				" OR\n"+
				"   (2) generate centroids from genotypes (i.e. -fromGenotypes (not the default))\n"+
				" OR\n"+
				"   (2) file with intensity only flags (i.e. flags=intensityFlags.dat (not the default))\n"+
				"   (3) centroid file for clustered markers (see " + Project.GENOTYPE_CENTROIDS_FILENAME + " in the Project properties file)\n"+
				"   (4) centroid file for intensity only markers (see " + Project.GENOTYPE_CENTROIDS_FILENAME + " in the Project properties file)\n"+
				" OR\n"+
				"   (2) recompute BAF/LRR and generate new Sample files using these centroids (i.e. compute=genotype.cent (not the default))\n"+
				"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help")
					|| args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("file=")) {
				centFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("import=")) {
				importFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("export=")) {
				exportFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("-fromGenotypes")) {
				fromGenotypes = true;
				numArgs--;
			} else if (args[i].startsWith("flags=")) {
				intensityFlags = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("compute=")) {
				compute = args[i].split("=")[1];
				numArgs--;
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}

		proj = new Project(filename, false);
		//		fromGenotypes = true;
		////		compute = "genotype.cent";
		//		
		centFile = "C:/data/AFFYtable/data/genotype.cent";
		exportFile = "data/genotype.cent.xln";
		try {
			if (exportFile != null) {
				//exportToText(proj, centFile, exportFile);
			} //else if (importFile != null) {
			//				importFromText(proj, importFile, centFile);
			//			} else if (fromGenotypes) {
			//				parseCentroidsFromGenotypes(proj, Array.booleanArray(proj.getSamples().length, true), 1);
			//			} else if (!compute.equals("")) {
			//				recompute(proj, compute);
			//			} else if (!intensityFlags.equals("")) {
			//				generateChimeraCentroids(proj, intensityFlags);
			//			} else { 
			//				parseIlluminaCentroidsFromCSV(proj, centFile);
			//			}
			//exportToText(proj, centFile, exportFile);
			//parseCentroids(proj, Array.booleanArray(proj.getSamples().length, true), 1);
			recompute(proj ,centFile );
		} catch (Exception e) {
			e.printStackTrace();
		}
	}


}
//private static float median(float[] array){
//float[] findMedian=array.clone();
//
//Arrays.sort(findMedian);
//float median;
//int middle = ((findMedian.length) / 2);
//if(array.length % 2 == 0){
//	float medianA = findMedian[middle];
//	float medianB = findMedian[middle-1];
//	median = (medianA + medianB) / 2;
//} else{
//	median = findMedian[middle + 1];
//}
//return median;
//}



//print $psid, "\t", sprintf ("%.4f", median (\@sig)), "\t", sprintf ("%.4f", mean (\@sig)), "\t", sprintf ("%.4f", sd (\@sig)), "\n";
//private static float[][] computeCN(MarkerData markerData , String Chr , String[] samples){
//float medianCN;
//float meanCN;
//float[] xs;
//float[][] cnCentroid;
//xs = markerData.getXs();		
//medianCN = median(xs);
//meanCN = Array.mean(xs);
//Array.stdev(xs, false);
//cnCentroid = new float[3][];
//cnCentroid[0][1] =medianCN;
//cnCentroid[1][1] = meanCN;
//
//for (int k = 0; k<samples.length; k++) {
//	//TODO check sex and get sex specific cluster
//	//if chr =X, if chr =Y etc
//	//else
//
//}
//return cnCentroid;
//}
//
//try {
//	writer = new PrintWriter(new FileWriter("C:/data/AFFYtable/JLtest.txt",true));	
//}catch (FileNotFoundException fnfe) {
//	//System.err.println("Error: file \""+Project.GENOTYPE_CENTROIDS_FILENAME +"JLtest"+"\" could not be written to (it's probably open)");
//	System.exit(1);
//} catch (IOException ioe) {
//	//System.err.println("Error reading file \""+Project.GENOTYPE_CENTROIDS_FILENAME +"JLtest"+"\"");
//	System.exit(2);
//}
//for (int k = 0; k<log2Xs.length; k++) {
//	writer.println(samples[k]+"\t"+markerName +"\t" +(log2Xs[k] - meanRs[2]));
//}

//if (AffyCentroids[0] != null && theta < AffyCentroids[0][0] || (AffyCentroids[1] == null && AffyCentroids[2] == null)) {
//			estimatedR = AffyCentroids[0][1];
//		} else if (AffyCentroids[1] != null && theta < AffyCentroids[1][0]) {
//			if (AffyCentroids[0] == null) {
//				estimatedR = AffyCentroids[1][1];
//			} else {
//				estimatedR = AffyCentroids[0][1]+(theta-AffyCentroids[0][0])*(AffyCentroids[1][1]-AffyCentroids[0][1])/(AffyCentroids[1][0]-AffyCentroids[0][0]);
//			}
//		} else if (AffyCentroids[2] != null && theta < AffyCentroids[2][0]) {
//			if (AffyCentroids[1] == null) {
//				estimatedR = AffyCentroids[2][1];
//			} else {
//				estimatedR = AffyCentroids[1][1]+(theta-AffyCentroids[1][0])*(AffyCentroids[2][1]-AffyCentroids[1][1])/(AffyCentroids[2][0]-AffyCentroids[1][0]);
//			}
//		} else {
//			if (AffyCentroids[2] == null) {
//				estimatedR = AffyCentroids[1][1];
//			} else {
//				estimatedR = AffyCentroids[2][1];
//			}
//		}
//		if (AffyCentroids[1] != null && theta < AffyCentroids[1][0]|| (AffyCentroids[1] == null && AffyCentroids[2] == null)) {
//			if(AffyCentroids[0] != null && AffyCentroids[1][0] - AffyCentroids[0][0] !=0){
//				estimatedR = AffyCentroids[0][1]+(theta-AffyCentroids[0][0])*(AffyCentroids[1][1]-AffyCentroids[0][1])/(AffyCentroids[1][0]-AffyCentroids[0][0]);
//			}
//			else{
//				estimatedR = AffyCentroids[0][1];
//						
//			}
//			if(estimatedR <0 ){
//				estimatedR =r;
//				System.out.println("WARNING: estimatedR < 0");
//			}
//			
//		}
//		else{
//			if (AffyCentroids[2] != null && AffyCentroids[1] != null && AffyCentroids[2][0]-AffyCentroids[1][0] != 0){
//				estimatedR = AffyCentroids[1][1]+(theta-AffyCentroids[1][0])*(AffyCentroids[2][1]-AffyCentroids[1][1])/(AffyCentroids[2][0]-AffyCentroids[1][0]);			
//			}
//			else{
//				estimatedR = AffyCentroids[1][1];
//						
//			}
//			if(estimatedR <0 ){
//				estimatedR =r;
//				System.out.println("WARNING: estimatedR < 0");
//			}		
//		}




//		//		if (AffyCentroids[2][1] < 0.0000001) {
//		//			AffyCentroids[2][0] = 1;
//		//		}
//		if (AffyCentroids[0] != null && theta < AffyCentroids[0][0] || (AffyCentroids[1] == null && AffyCentroids[2] == null)) {
//			estimatedR = AffyCentroids[0][1];
//		} else if (AffyCentroids[1] != null && theta < AffyCentroids[1][0]) {
//			if (AffyCentroids[0] == null) {
//				estimatedR = AffyCentroids[1][1];
//			} else {
//				estimatedR = AffyCentroids[0][1]+(theta-AffyCentroids[0][0])*(AffyCentroids[1][1]-AffyCentroids[0][1])/(AffyCentroids[1][0]-AffyCentroids[0][0]);
//			}
//		} else if (AffyCentroids[2] != null && theta < AffyCentroids[2][0]) {
//			if (AffyCentroids[1] == null) {
//				estimatedR = AffyCentroids[2][1];
//			} else {
//				estimatedR = AffyCentroids[1][1]+(theta-AffyCentroids[1][0])*(AffyCentroids[2][1]-AffyCentroids[1][1])/(AffyCentroids[2][0]-AffyCentroids[1][0]);
//			}
//		} else {
//			if (AffyCentroids[2] == null) {
//				estimatedR = AffyCentroids[1][1];
//			} else {
//				estimatedR = AffyCentroids[2][1];
//			}
//		}
