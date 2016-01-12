package cnv.manage;

import java.io.File;
import java.util.ArrayList;
import java.util.Hashtable;

import common.Files;
import common.ext;
import cnv.analysis.PennCNVPrep;
import cnv.analysis.PennCNVPrep.ShadowSample;
import cnv.analysis.pca.PrincipalComponentsIntensity;
import cnv.analysis.pca.PrincipalComponentsResiduals;
import cnv.filesys.MarkerData;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import cnv.filesys.Sample;

// PRincipal COmponents Residuals - PR [o] C[t]O R
public class PRoCtOR {
    
    // relevant code from PennCNVPrep
    /*
     
        PrincipalComponentsResiduals principalComponentsResiduals = loadPcResids(proj, numComponents);
        if (principalComponentsResiduals == null) {
            return;
        }
        int[] sex = getSampleSex(proj);
        if (sex == null && proj.getAutosomalMarkers().length != proj.getMarkerNames().length) {
    
            proj.getLog().reportError("Error - missing sex codes");
            return;
        }
        if (markerFile == null) {
            markers = proj.getMarkerNames();
            proj.getLog().report("Info - a file of markers was not provided, exporting all in this batch");
        } else {
            markers = HashVec.loadFileToStringArray(markerFile, false, new int[] { 0 }, false);
            proj.getLog().report("Info - loaded " + markers.length + " markers from " + markerFile + " to export");
    
        }
        PennCNVPrep specialPennCNVFormat = new PennCNVPrep(proj, principalComponentsResiduals, null, proj.getSamplesToInclude(null), sex, markers, numComponents, dir, svdRegression, numThreads, numMarkerThreads);
        specialPennCNVFormat.exportSpecialMarkerDataMoreThreads(tmpDir, preserveBafs);
     
     */
    
    private String[] getMarkerNames(Project proj) {
        return proj.getMarkerNames();
    }
    
    private int getNumThreads() {
        return 4;
    }
    
    private int getMarkerBuffer() {
        return 20000;
    }
    
    private int getNumComponents() {
        return 40;
    }
    
    private boolean getSVD() {
        return false;
    }
    
    private static final String SHADOW_DIR = "shadowSamples/";
    
    void run(Project proj) {
        System.gc();
        String[] markerNames = getMarkerNames(proj);
        int numThreads = getNumThreads();
        int mkrBuffer = getMarkerBuffer();
        int numComponents = getNumComponents();
        boolean svdRegression = getSVD();
        
        String[] samples = proj.getSamples();
        String dir = proj.PROJECT_DIRECTORY.getValue() + SHADOW_DIR;
        new File(dir).mkdirs();
        
        proj.getLog().report("Info - checking for existing files in " + dir + "...");
        
        boolean allExist = true;
        for (int i = 0; i < samples.length; i++) {
            if (!Files.exists(dir + samples[i] + Sample.SAMPLE_DATA_FILE_EXTENSION)) {
                allExist = false;
            }
        }
        
        long msFingerprint = proj.getMarkerSet().getFingerprint();
        Hashtable<String, Float> allOutliers = new Hashtable<String, Float>();
		MarkerSet markerSet = proj.getMarkerSet();
        int numMarkers = markerSet.getPositions().length;
        ArrayList<String> notCorrected = new ArrayList<String>();
        PrincipalComponentsResiduals principalComponentsResiduals = PennCNVPrep.loadPcResids(proj, numComponents);
        if (principalComponentsResiduals == null) {
            return;
        }
        int[] sex = PennCNVPrep.getSampleSex(proj);
        if (sex == null && proj.getAutosomalMarkers().length != proj.getMarkerNames().length) {
            proj.getLog().reportError("Error - missing sex codes");
            return;
        }
        boolean[] samplesToUseCluster = proj.getSamplesToInclude(null);
        if (!allExist) {// if all files exist we skip this export
            proj.getLog().report("Info - detected that not all files exist");

            proj.getLog().report(ext.getTime() + "]\tLoading all markers into memory...");
            MarkerData[] allMarkers = new MarkerData[numMarkers];
			MDL dataLoader = new MDL(proj, markerSet, markerNames, numThreads, mkrBuffer);
            dataLoader.setReportEvery(50000);
            int index = 0;
            while (dataLoader.hasNext()) {
                MarkerData markerData = dataLoader.next();
                MarkerData markerDataToStore;
                PrincipalComponentsIntensity principalComponentsIntensity = new PrincipalComponentsIntensity(principalComponentsResiduals, markerData, true, sex, samplesToUseCluster, 1, 0, null, true, svdRegression, 2, 5, PrincipalComponentsIntensity.DEFAULT_RESID_STDV_FILTER, PrincipalComponentsIntensity.DEFAULT_CORRECTION_RATIO, numThreads, false, null);
                principalComponentsIntensity.correctXYAt(numComponents);
                if (principalComponentsIntensity.isFail()) {
                    notCorrected.add(markerData.getMarkerName());
                    markerDataToStore = markerData;
                } else {
                    byte[] abGenotypes = principalComponentsIntensity.getGenotypesUsed();
                    float[][] correctedXY = principalComponentsIntensity.getCorrectedIntensity(PrincipalComponentsIntensity.XY_RETURN, true);
                    float[][] correctedLRRBAF = principalComponentsIntensity.getCorrectedIntensity(PrincipalComponentsIntensity.BAF_LRR_RETURN, true);// for now
                    markerDataToStore = new MarkerData(markerData.getMarkerName(), markerData.getChr(), markerData.getPosition(), markerData.getFingerprint(), markerData.getGCs(), null, null, correctedXY[0], correctedXY[1], null, null, correctedLRRBAF[0], correctedLRRBAF[1], abGenotypes, abGenotypes);
                }
                allMarkers[index++] = markerDataToStore;
            }
            dataLoader.shutdown();
            dataLoader = null;
            System.gc();
            proj.getLog().report(ext.getTime() + "]\tExporting " + samples.length + " sample files...");
            for (int i = 0; i < samples.length; i++) {
                ShadowSample shadowSample = new ShadowSample(samples[i], proj.getMarkerNames());
                for (int m = 0; m < allMarkers.length; m++) {
                    shadowSample.addData(samples[i], m, allMarkers[m].getMarkerName(), allMarkers[m].getXs()[i], allMarkers[m].getYs()[i], allMarkers[m].getGCs()[i], allMarkers[m].getBAFs()[i], allMarkers[m].getLRRs()[i], allMarkers[m].getAbGenotypes()[i], proj.getLog());
                }
                shadowSample.writeShadow(proj, dir, msFingerprint, allOutliers);
                shadowSample = null;
                System.gc();
                if ((i + 1) % 50000 == 0) {
                    float usedMemory = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
                    float freeMemory = Runtime.getRuntime().maxMemory() - usedMemory;
                    float maxMemory = Runtime.getRuntime().maxMemory();
                    proj.getLog().report(ext.getTime() + "]\tSample Written = " + Math.round(((double) i / (double) samples.length * 100.0)) + "%\tFree memory: " + Math.round(((double) freeMemory / (double) maxMemory * 100.0)) + "%");
                }
            }
            String outlierFile = proj.PROJECT_DIRECTORY.getValue() + SHADOW_DIR + "outliers.ser";
            Files.writeSerial(allOutliers, outlierFile);
        } else {
            proj.getLog().report("Info - detected that all " + samples.length + " shadow samples exist in " + dir + ", skipping export...");
        }
    }

    
    
    public static void main(String[] args) {
        int numArgs = args.length;
        String filename = "D:/projects/Poynter.properties";

        String usage = "\n" + 
                       "cnv.manage.PRoCtOR requires 0-1 arguments\n" + 
                       "   (1) project properties filename (i.e. proj=" + filename + " (default))\n" + "";

        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
                System.err.println(usage);
                System.exit(1);
            } else if (args[i].startsWith("proj=")) {
                filename = args[i].split("=")[1];
                numArgs--;
            } else {
                System.err.println("Error - invalid argument: " + args[i]);
            }
        }
        if (numArgs != 0) {
            System.err.println(usage);
            System.exit(1);
        }
        try {
            Project proj = new Project(filename, false);
            new PRoCtOR().run(proj);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
