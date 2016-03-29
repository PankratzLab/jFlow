package cnv.manage;

import gwas.Qc;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;

import cnv.Launch;
import cnv.analysis.AnalysisFormats;
import cnv.analysis.Mosaicism;
import cnv.analysis.pca.PCA;
import cnv.analysis.pca.PrincipalComponentsApply;
import cnv.analysis.pca.PrincipalComponentsCompute;
import cnv.analysis.pca.PrincipalComponentsResiduals;
import cnv.filesys.ABLookup;
import cnv.filesys.Centroids;
import cnv.filesys.MarkerData;
import cnv.filesys.Pedigree;
import cnv.filesys.Project;
import cnv.filesys.Sample;
import cnv.gui.GenvisisWorkflowGUI;
import cnv.hmm.CNVCaller;
import cnv.qc.GcAdjustor;
import cnv.qc.GcAdjustorParameter;
import cnv.qc.GcAdjustor.GCAdjustorBuilder;
import cnv.qc.GcAdjustorParameter.GcAdjustorParameters;
import cnv.qc.SampleQC;
import cnv.var.SampleData;
import common.Aliases;
import common.Array;
import common.Files;
import common.Logger;
import common.ext;

public class GenvisisWorkflow {
    
    Project proj;
    Logger log;
    private GenvisisWorkflowGUI gui;
    private Launch launch;
    
    static final STEP S1I_CREATE_MKR_POS = new STEP("Create Marker Positions (if not already exists)", 
                      "", 
                      new String[][]{{"An Illumina SNP_map file."}}, 
                      new RequirementInputType[][]{{RequirementInputType.FILE}}) {
        
        @Override
        public void setNecessaryPreRunProperties(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            // not needed for step
        }
        
        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            proj.getLog().report("Generating marker positions file");
            String filename = variables.get(this).get(0);
            cnv.manage.Markers.generateMarkerPositions(proj, filename);
        }
        
        @Override
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, Boolean> stepSelections, HashMap<STEP, ArrayList<String>> variables) {
            return new boolean[][]{{Files.exists(variables.get(this).get(0))}};
        }
        
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            String filename = proj.getLocationOfSNP_Map(false);
            if (filename == null) {
                filename = "";
            }
            return new Object[]{filename};
        }

        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            boolean mkrPosFile = Files.exists(proj.MARKER_POSITION_FILENAME.getValue(false, false));
//            boolean mkrSetFile = Files.exists(proj.MARKERSET_FILENAME.getValue(false, false));
            return mkrPosFile/* && mkrSetFile*/;
        };
        
        @Override
        public String getCommandLine(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String projFile = proj.getPropertyFilename();
            String filename = variables.get(this).get(0);
            return "jcp cnv.manage.Markers proj=" + projFile + " snps=" + filename;
        }
        
    };
    
    static final STEP S2I_PARSE_SAMPLES = new STEP("Parse Illumina Sample Files", 
                     "", 
                     new String[][]{{"[Create Marker Positions] step must be selected and valid.", "Parsed markerPositions file must already exist."}, {"Number of Threads to Use"}}, 
                     new RequirementInputType[][]{{RequirementInputType.NONE, RequirementInputType.FILE}, {RequirementInputType.INT}}) {

        @Override
        public void setNecessaryPreRunProperties(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String projFile = proj.MARKER_POSITION_FILENAME.getValue(false, false);
            String mkrFile = variables.get(this).get(0);
            mkrFile = ext.verifyDirFormat(mkrFile);
            mkrFile = mkrFile.substring(0, mkrFile.length() - 1);
            if (!mkrFile.equals(projFile)) {
                proj.MARKER_POSITION_FILENAME.setValue(mkrFile);
            }
            int numThreads = proj.NUM_THREADS.getValue();
            try {
                numThreads = Integer.parseInt(variables.get(this).get(1));
            } catch (NumberFormatException e) {}
            if (numThreads != proj.NUM_THREADS.getValue()) {
                proj.NUM_THREADS.setValue(numThreads);
            }
        }
        
        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            int numThreads = proj.NUM_THREADS.getValue();
//            left as unneeded/pre-run checks ensure int value
//            try {
//                numThreads = Integer.parseInt(variables.get(this).get(1));
//            } catch (NumberFormatException e) {}
            proj.getLog().report("Parsing sample files");
            int retCode = cnv.manage.SourceFileParser.createFiles(proj, numThreads);
            switch (retCode) {
            case 0:
                this.setFailed();
                this.failReasons.add("Operation failure, please check log for more information.");
                break;
            case 1:
                break;
            case 6:
                this.failReasons.add("ABLookup required but wasn't found.");
                break;
            }
        }
        
        @Override
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, Boolean> stepSelections, HashMap<STEP, ArrayList<String>> variables) {
        	int numThreads = -1;
        	try {
        		numThreads = Integer.parseInt(variables.get(this).get(1));
        	} catch (NumberFormatException e) {}
            return new boolean[][]{
                    { stepSelections.get(S1I_CREATE_MKR_POS) && S1I_CREATE_MKR_POS.hasRequirements(proj, stepSelections, variables),
                    	Files.exists(ext.verifyDirFormat(variables.get(this).get(0))),},
                	{numThreads > 0}
            };
        }
        
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            return new Object[]{proj.MARKER_POSITION_FILENAME.getValue(false, false), proj.NUM_THREADS.getValue()};
        }

        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String sampleDirectory = proj.SAMPLE_DIRECTORY.getValue(false, false);
            boolean mkrSetFile = Files.exists(proj.MARKERSET_FILENAME.getValue(false, false));
            return mkrSetFile && Files.exists(sampleDirectory) && Files.list(sampleDirectory, Sample.SAMPLE_DATA_FILE_EXTENSION, false).length > 0 && proj.getSampleList() != null && proj.getSampleList().getSamples().length > 0;
        }
        
        @Override
        public String getCommandLine(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String projPropFile = proj.getPropertyFilename();
            StringBuilder kvCmd = new StringBuilder("jcp cnv.filesys.Project proj=").append(projPropFile);
            StringBuilder kvPairs = new StringBuilder();
            String projFile = proj.MARKER_POSITION_FILENAME.getValue(false, false);
            String mkrFile = variables.get(this).get(0);
            mkrFile = ext.verifyDirFormat(mkrFile);
            mkrFile = mkrFile.substring(0, mkrFile.length() - 1);
            if (!mkrFile.equals(projFile)) {
                kvPairs.append(" MARKER_POSITION_FILENAME=").append(mkrFile);
            }
            int numThreads = proj.NUM_THREADS.getValue();
            try {
                numThreads = Integer.parseInt(variables.get(this).get(1));
            } catch (NumberFormatException e) {}
            if (numThreads != proj.NUM_THREADS.getValue()) {
                kvPairs.append(" NUM_THREADS=").append(numThreads);
            }
            StringBuilder command = new StringBuilder();
            if (kvPairs.length() != 0) {
                command.append(kvCmd).append(kvPairs).append("\n");
            }
            command.append("jcp cnv.manage.SourceFileParser proj=").append(projPropFile).append(" threads=").append(numThreads);
            return command.toString();
        }
        
    };
    
    static final STEP S2A_PARSE_SAMPLES = new STEP("Parse Sample Files", 
            "", 
            new String[][]{{"markerPositions file must already exist."}, {"Number of Threads to Use"}}, 
            new RequirementInputType[][]{{RequirementInputType.FILE}, {RequirementInputType.INT}}) {
        
        @Override
        public void setNecessaryPreRunProperties(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String projFile = proj.MARKER_POSITION_FILENAME.getValue(false, false);
            String mkrFile = variables.get(this).get(0);
            mkrFile = ext.verifyDirFormat(mkrFile);
            mkrFile = mkrFile.substring(0, mkrFile.length() - 1);
            if (!mkrFile.equals(projFile)) {
                proj.MARKER_POSITION_FILENAME.setValue(mkrFile);
            }
            int numThreads = proj.NUM_THREADS.getValue();
            try {
                numThreads = Integer.parseInt(variables.get(this).get(1));
            } catch (NumberFormatException e) {}
            if (numThreads != proj.NUM_THREADS.getValue()) {
                proj.NUM_THREADS.setValue(numThreads);
            }
        }
        
        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            int numThreads = proj.NUM_THREADS.getValue();
            numThreads = Integer.parseInt(variables.get(this).get(1));
            proj.getLog().report("Parsing sample files");
            int retCode = cnv.manage.SourceFileParser.createFiles(proj, numThreads);
            switch (retCode) {
            case 0:
                this.setFailed();
                this.failReasons.add("Operation failure, please check log for more information.");
                break;
            case 1:
                break;
            case 6:
                this.failReasons.add("ABLookup required but wasn't found.");
                break;
            }
        }
        
        @Override
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, Boolean> stepSelections, HashMap<STEP, ArrayList<String>> variables) {
            int numThreads = -1;
            try {
                numThreads = Integer.parseInt(variables.get(this).get(1));
            } catch (NumberFormatException e) {}
            return new boolean[][]{
                    {Files.exists(ext.verifyDirFormat(variables.get(this).get(0))),},
                    {numThreads > 0}
            };
        }
        
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            return new Object[]{proj.MARKER_POSITION_FILENAME.getValue(false, false), proj.NUM_THREADS.getValue()};
        }
        
        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String sampleDirectory = proj.SAMPLE_DIRECTORY.getValue(false, false);
            boolean mkrSetFile = Files.exists(proj.MARKERSET_FILENAME.getValue(false, false));
            return mkrSetFile && Files.exists(sampleDirectory) && Files.list(sampleDirectory, Sample.SAMPLE_DATA_FILE_EXTENSION, false).length > 0 && proj.getSampleList() != null && proj.getSampleList().getSamples().length > 0;
        }
        
        @Override
        public String getCommandLine(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String projPropFile = proj.getPropertyFilename();
            StringBuilder kvCmd = new StringBuilder("jcp cnv.filesys.Project proj=").append(projPropFile);
            StringBuilder kvPairs = new StringBuilder();
            String projFile = proj.MARKER_POSITION_FILENAME.getValue(false, false);
            String mkrFile = variables.get(this).get(0);
            mkrFile = ext.verifyDirFormat(mkrFile);
            mkrFile = mkrFile.substring(0, mkrFile.length() - 1);
            if (!mkrFile.equals(projFile)) {
                kvPairs.append(" MARKER_POSITION_FILENAME=").append(mkrFile);
            }
            int numThreads = proj.NUM_THREADS.getValue();
            try {
                numThreads = Integer.parseInt(variables.get(this).get(1));
            } catch (NumberFormatException e) {}
            if (numThreads != proj.NUM_THREADS.getValue()) {
                kvPairs.append(" NUM_THREADS=").append(numThreads);
            }
            StringBuilder command = new StringBuilder();
            if (kvPairs.length() != 0) {
                command.append(kvCmd).append(kvPairs).append("\n");
            }
            command.append("jcp cnv.manage.SourceFileParser proj=").append(projPropFile).append(" threads=").append(numThreads);
            return command.toString();
        }
        
    };
    
    static final STEP S3_CREATE_SAMPLEDATA = new STEP("Create SampleData.txt File", 
                         "", 
                         new String[][]{{"[Parse Sample Files] step must be selected and valid (will create a minimal SampleData.txt file)", "Parsed sample files must already exist (will create a minimal SampleData.txt file)", "A tab-delimited .PED format file with header \"" + Array.toStr(MitoPipeline.PED_INPUT, ", ") + "\"", "A Sample_Map.csv file, with at least two columns having headers \"" + MitoPipeline.SAMPLEMAP_INPUT[1] + "\" and \"" + MitoPipeline.SAMPLEMAP_INPUT[2] + "\""}}, 
                         new RequirementInputType[][]{{RequirementInputType.NONE, RequirementInputType.DIR, RequirementInputType.FILE, RequirementInputType.FILE}}) {

        @Override
        public void setNecessaryPreRunProperties(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String projDir = proj.SAMPLE_DIRECTORY.getValue(false, false);
            String setDir = variables.get(this).get(0);
            if (!ext.verifyDirFormat(setDir).equals(projDir)) {
                proj.SAMPLE_DIRECTORY.setValue(setDir);
            }
        }
        
        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String pedFile = variables.get(this).get(1);
            String sampleMapCsv = variables.get(this).get(2);
            proj.getLog().report("Creating SampleData.txt");
            /*int retStat = */SampleData.createSampleData(pedFile, sampleMapCsv, proj);
        }
        
        @Override
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, Boolean> stepSelections, HashMap<STEP, ArrayList<String>> variables) {
            String sampDir = variables.get(this).get(0);
            STEP parseStep = stepSelections.containsKey(S2I_PARSE_SAMPLES) ? S2I_PARSE_SAMPLES : S2A_PARSE_SAMPLES;
            boolean checkStepParseSamples = stepSelections.get(parseStep) && parseStep.hasRequirements(proj, stepSelections, variables);
            return new boolean[][]{
                    {checkStepParseSamples,
                        (Files.exists(sampDir) && Files.list(sampDir, ".sampRAF", proj.JAR_STATUS.getValue()).length > 0), 
                            Files.exists(variables.get(this).get(1)), 
                            Files.exists(variables.get(this).get(2))}};
        }
        
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            return new Object[]{proj.SAMPLE_DIRECTORY.getValue(false, false), "", ""};
        }

        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            return Files.exists(proj.SAMPLE_DATA_FILENAME.getValue(false, false));
        }
        @Override
        public String getCommandLine(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String projPropFile = proj.getPropertyFilename();
            String pedFile = variables.get(this).get(1);
            String sampleMapCsv = variables.get(this).get(2);
            StringBuilder cmd = new StringBuilder();
            cmd.append("jcp cnv.var.SampleData proj=").append(projPropFile);
            if (!"".equals(pedFile)) {
                cmd.append(" ped=").append(pedFile);
            }
            if (!"".equals(sampleMapCsv)) {
                cmd.append(" sampleMap=").append(sampleMapCsv);
            }
            return cmd.toString();
        }
        
    };
    
    static final STEP S4_TRANSPOSE_TO_MDF = new STEP("Transpose Data into Marker-Dominant Files", 
                        "", 
                        new String[][]{{"[Parse Sample Files] step must be selected and valid.", "Parsed sample files must already exist."}}, 
                        new RequirementInputType[][]{{RequirementInputType.NONE, RequirementInputType.DIR}}) {

        @Override
        public void setNecessaryPreRunProperties(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String projDir = proj.SAMPLE_DIRECTORY.getValue(false, false);
            String setDir = variables.get(this).get(0);
            if (!ext.verifyDirFormat(setDir).equals(projDir)) {
                proj.SAMPLE_DIRECTORY.setValue(setDir);
            }
        }
        
        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            proj.getLog().report("Transposing data");
            TransposeData.transposeData(proj, 2000000000, false); // compact if no LRR was provided
        }
        
        @Override
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, Boolean> stepSelections, HashMap<STEP, ArrayList<String>> variables) {
            String sampDir = variables.get(this).get(0);
            STEP parseStep = stepSelections.containsKey(S2I_PARSE_SAMPLES) ? S2I_PARSE_SAMPLES : S2A_PARSE_SAMPLES;
            boolean checkStepParseSamples = stepSelections.get(parseStep) && parseStep.hasRequirements(proj, stepSelections, variables);
            return new boolean[][]{
                    {checkStepParseSamples,
                    (Files.exists(sampDir) && Files.list(sampDir, ".sampRAF", proj.JAR_STATUS.getValue()).length > 0),}
            };
        }
        
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            return new Object[]{proj.SAMPLE_DIRECTORY.getValue(false, false)};
        }

        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            return Files.list(proj.MARKER_DATA_DIRECTORY.getValue(false, false), MarkerData.MARKER_DATA_FILE_EXTENSION, false).length > 0;
        }
        
        @Override
        public String getCommandLine(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String kvCmd = "";
            String projDir = proj.SAMPLE_DIRECTORY.getValue(false, false);
            String setDir = variables.get(this).get(0);
            if (!ext.verifyDirFormat(setDir).equals(projDir)) {
                kvCmd += " SAMPLE_DIRECTORY=" + setDir;
            }
            String projPropFile = proj.getPropertyFilename();
            StringBuilder cmd = new StringBuilder();
            if (kvCmd.length() > 0) {
                cmd.append("jcp cnv.filesys.Project proj=" + projPropFile).append(kvCmd).append("\n");
            }
            return cmd.append("jcp cnv.manage.TransposeData -transpose proj=" + projPropFile + " max=" + 2000000000).toString();
        }
        
    };
    static final STEP S5_MARKER_QC = new STEP("Run Marker QC Metrics", "",
            new String[][]{{"Marker Call-Rate Filter Threshold"},
    						{"[Parse Sample Files] step must be selected and valid.", "A MarkerSet file must already exist."}, 
    						{"[Parse Sample Files] step must be selected and valid.", "A SampleList file must already exist."},
    						{"Export all markers in project.", "A targetMarkers files listing the markers to QC."},
    						{"Number of threads to use."}},
            new RequirementInputType[][]{{RequirementInputType.INT},
    									 {RequirementInputType.NONE, RequirementInputType.FILE}, 
    									 {RequirementInputType.NONE, RequirementInputType.FILE},
    									 {RequirementInputType.NONE, RequirementInputType.FILE},
    									 {RequirementInputType.INT}}
            ) {
    
        @Override
        public void setNecessaryPreRunProperties(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String mkrPosProj = proj.MARKERSET_FILENAME.getValue(false, false);
            String mkrPosFile = variables.get(this).get(1);
            String setSampList = proj.SAMPLELIST_FILENAME.getValue(false, false);
            String sampList = variables.get(this).get(2);
            String setTgtFile = proj.TARGET_MARKERS_FILENAMES.getValue().length > 0 ? proj.TARGET_MARKERS_FILENAMES.getValue()[0] : "";
            String tgtFile = variables.get(this).get(3);
            if (!mkrPosProj.equals(mkrPosFile)) {
                proj.MARKERSET_FILENAME.setValue(mkrPosFile);
            }
            if (!ext.verifyDirFormat(setSampList).equals(sampList)) {
                proj.SAMPLELIST_FILENAME.setValue(sampList);
            }
            if (!"".equals(tgtFile) && !setTgtFile.equals(tgtFile)) {
                // String[] arr = proj.TARGET_MARKERS_FILENAMES.getValue();
                // arr[0] = tgtFile;
                proj.TARGET_MARKERS_FILENAMES.addValue(tgtFile);
            }
            int numThreads = proj.NUM_THREADS.getValue();
            try {
                numThreads = Integer.parseInt(variables.get(this).get(4));
            } catch (NumberFormatException e) {}
            if (numThreads != proj.NUM_THREADS.getValue()) {
                proj.NUM_THREADS.setValue(numThreads);
            }
        }
        
        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            double markerCallRateFilter = Double.parseDouble(variables.get(this).get(0));
            String tgtFile = variables.get(this).get(3);
            int numThreads = proj.NUM_THREADS.getValue();
            try {
                numThreads = Integer.parseInt(variables.get(this).get(4));
            } catch (NumberFormatException e) {}
            String markersToQC = proj.PROJECT_DIRECTORY.getValue() + MitoPipeline.FILE_BASE + "_" + MitoPipeline.MARKERS_TO_QC_FILE;
            String markersABCallrate = proj.PROJECT_DIRECTORY.getValue() + MitoPipeline.FILE_BASE + "_" + MitoPipeline.MARKERS_FOR_ABCALLRATE;
            MitoPipeline.qcMarkers(proj, "".equals(tgtFile) ? null : tgtFile, markersToQC, markersABCallrate, markerCallRateFilter, numThreads);
            // TODO Replace with MarkerMetrics.fullQC
        }
        
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            return new Object[]{0.98, proj.MARKERSET_FILENAME.getValue(), proj.SAMPLELIST_FILENAME.getValue(), ""/*proj.TARGET_MARKERS_FILENAMES.getValue()[0]*/, proj.NUM_THREADS.getValue()};
        }
        
        @Override
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, Boolean> stepSelections, HashMap<STEP, ArrayList<String>> variables) {
            double mkr = -1;
            try {
                mkr = Double.parseDouble(variables.get(this).get(0));
            } catch (NumberFormatException e) {}
            String mkrPosFile = variables.get(this).get(1);
            String sampList = variables.get(this).get(2);
            String tgtFile = variables.get(this).get(3);
        	int numThreads = -1;
        	try {
        		numThreads = Integer.parseInt(variables.get(this).get(4));
        	} catch (NumberFormatException e) {}
            boolean step12 = Files.exists(mkrPosFile);
            boolean step2 = stepSelections.get(S2I_PARSE_SAMPLES) && S2I_PARSE_SAMPLES.hasRequirements(proj, stepSelections, variables);
            boolean step22 = Files.exists(sampList);
            boolean step3 = Files.exists(tgtFile);
            return new boolean[][]{{mkr != -1}, {step2, step12}, {step2, step22}, {true, step3}, {numThreads > 0}};
        }
        
        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String markersForABCallRate = proj.PROJECT_DIRECTORY.getValue() + MitoPipeline.FILE_BASE + "_" + MitoPipeline.MARKERS_FOR_ABCALLRATE;
            if (!Files.exists(markersForABCallRate)) {
                return false;
            }
            return true;
        }
        @Override
        public String getCommandLine(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            return "## << Marker QC >> Not Implemented For Command Line Yet ##"; // TODO
        }
    };
//    static final STEP S6_EXTRACT_LRRSD_AND_FILTER = new STEP("Extract Lrrsd.xln File and Filter Samples by CallRate", 
//            "", 
//            new String[][]{
//                  {"[Parse Sample Files] step must be selected and valid.", "Parsed sample files must already exist."}, 
//                  {"Number of threads to use."},
//                  {"Sample call rate filter threshold."},
//              }, 
//            new RequirementInputType[][]{{RequirementInputType.NONE, RequirementInputType.DIR}, {RequirementInputType.INT}, {RequirementInputType.STRING}}) {
//
//        @Override
//        public void setNecessaryPreRunProperties(Project proj, HashMap<STEP, ArrayList<String>> variables) {
//            String projDir = proj.SAMPLE_DIRECTORY.getValue(false, false);
//            String setDir = variables.get(this).get(0);
//            if (!ext.verifyDirFormat(setDir).equals(projDir)) {
//                proj.SAMPLE_DIRECTORY.setValue(setDir);
//            }
//            int numThreads = proj.NUM_THREADS.getValue();
//            try {
//                numThreads = Integer.parseInt(variables.get(this).get(1));
//            } catch (NumberFormatException e) {
//            }
//            if (numThreads != proj.NUM_THREADS.getValue()) {
//                proj.NUM_THREADS.setValue(numThreads);
//            }
//        }
//
//        @Override
//        public void run(Project proj, HashMap<STEP, ArrayList<String>> variables) {
//            proj.getLog().report("Running LrrSd");
//            int numThreads = proj.NUM_THREADS.getValue();
//            try {
//                numThreads = Integer.parseInt(variables.get(this).get(1));
//            } catch (NumberFormatException e) {
//            }
//            String callRate = variables.get(this).get(2);
//            String markersForAB = Files.exists(proj.PROJECT_DIRECTORY.getValue() + MitoPipeline.FILE_BASE + "_" + MitoPipeline.MARKERS_FOR_ABCALLRATE) ? proj.PROJECT_DIRECTORY.getValue() + MitoPipeline.FILE_BASE + "_" + MitoPipeline.MARKERS_FOR_ABCALLRATE : proj.PROJECT_DIRECTORY.getValue() + MitoPipeline.FILE_BASE + "_" + MitoPipeline.MARKERS_TO_QC_FILE;
//            String markersForEverything = proj.PROJECT_DIRECTORY.getValue() + MitoPipeline.FILE_BASE + "_" + MitoPipeline.MARKERS_TO_QC_FILE;
//            
//            markersForAB = Files.exists(markersForAB) ? markersForAB : null;
//            markersForEverything = Files.exists(markersForEverything) ? markersForEverything : null;
//            
//            LrrSd.filterSamples(proj, MitoPipeline.FILE_BASE, markersForAB, markersForEverything, numThreads, callRate, null);
//        }
//
//        @Override
//        public boolean[][] checkRequirements(Project proj, HashMap<STEP, Boolean> stepSelections, HashMap<STEP, ArrayList<String>> variables) {
//            int numThreads = -1;
//            try {
//                numThreads = Integer.parseInt(variables.get(this).get(1));
//            } catch (NumberFormatException e) {
//            }
//            String sampDir = variables.get(this).get(0);
//            STEP parseStep = stepSelections.containsKey(S2I_PARSE_SAMPLES) ? S2I_PARSE_SAMPLES : S2A_PARSE_SAMPLES;
//            boolean checkStepParseSamples = stepSelections.get(parseStep) && parseStep.hasRequirements(proj, stepSelections, variables);
//            boolean validCallRate = false;
//            try {
//                Double.parseDouble(variables.get(this).get(2));
//                validCallRate = true;
//            } catch (NumberFormatException e) {
//            }
//            return new boolean[][] { { checkStepParseSamples, (Files.exists(sampDir) && Files.list(sampDir, ".sampRAF", proj.JAR_STATUS.getValue()).length > 0), }, { numThreads > 0 }, { validCallRate } };
//        }
//
//        @Override
//        public Object[] getRequirementDefaults(Project proj) {
//            return new Object[] { proj.SAMPLE_DIRECTORY.getValue(false, false), proj.NUM_THREADS.getValue(), "0.95" };
//        }
//
//        @Override
//        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variables) {
//            return Files.exists(proj.SAMPLE_QC_FILENAME.getValue(false, false)) || (Files.exists(MitoPipeline.FILE_BASE + MitoPipeline.PCA_SAMPLES) && Files.exists(MitoPipeline.FILE_BASE + MitoPipeline.PCA_SAMPLES_SUMMARY));
//        }
//
//        @Override
//        public String getCommandLine(Project proj, HashMap<STEP, ArrayList<String>> variables) {
//            proj.getLog().report("Running LrrSd");
//            int numThreads = proj.NUM_THREADS.getValue();
//            try {
//                numThreads = Integer.parseInt(variables.get(this).get(1));
//            } catch (NumberFormatException e) {
//            }
//            String callRate = variables.get(this).get(2);
//            String markersForAB = Files.exists(proj.PROJECT_DIRECTORY.getValue() + MitoPipeline.FILE_BASE + "_" + MitoPipeline.MARKERS_FOR_ABCALLRATE) ? proj.PROJECT_DIRECTORY.getValue() + MitoPipeline.FILE_BASE + "_" + MitoPipeline.MARKERS_FOR_ABCALLRATE : proj.PROJECT_DIRECTORY.getValue() + MitoPipeline.FILE_BASE + "_" + MitoPipeline.MARKERS_TO_QC_FILE;
//            String markersForEverything = proj.PROJECT_DIRECTORY.getValue() + MitoPipeline.FILE_BASE + "_" + MitoPipeline.MARKERS_TO_QC_FILE;
//            
//            markersForAB = Files.exists(markersForAB) ? markersForAB : null;
//            markersForEverything = Files.exists(markersForEverything) ? markersForEverything : null;
//            
//            String projPropFile = proj.getPropertyFilename();
//            StringBuilder cmd = new StringBuilder();
//            cmd.append("jcp cnv.qc.LrrSd -filter")
//                    .append(" proj=").append(projPropFile)
//                    .append(" outBase=").append(MitoPipeline.FILE_BASE);
//            if (markersForAB != null) {
//                cmd.append(" callRateMarkers=").append(markersForAB);
//            }
//            if (markersForEverything != null) {
//                cmd.append(" otherMarkers=").append(markersForEverything);
//            }
//            cmd.append(" threads=").append(numThreads)
//                .append(" callRate=").append(callRate);
//            return cmd.toString();
//        }
//    };

    static final STEP S6_SEX_CHECKS = new STEP("Run Sex Checks", 
                  "", 
                  new String[][]{{"[Parse Sample Files] step must be selected and valid.", "Parsed sample files must already exist."},
                                 {"[Create SampleData.txt File] step must be selected and valid.", "SampleData.txt file must already exist."}}, 
                  new RequirementInputType[][]{{RequirementInputType.NONE, RequirementInputType.DIR}, 
                                               {RequirementInputType.NONE, RequirementInputType.FILE}}) {

        @Override
        public void setNecessaryPreRunProperties(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String projDir = proj.SAMPLE_DIRECTORY.getValue(false, false);
            String setDir = variables.get(this).get(0);
            if (!ext.verifyDirFormat(setDir).equals(projDir)) {
                proj.SAMPLE_DIRECTORY.setValue(setDir);
            }
        }
        
        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            proj.getLog().report("Running SexCheck");
            cnv.qc.SexChecks.sexCheck(proj);
        }
        
        @Override
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, Boolean> stepSelections, HashMap<STEP, ArrayList<String>> variables) {
            String sampDir = variables.get(this).get(0);
            STEP parseStep = stepSelections.containsKey(S2I_PARSE_SAMPLES) ? S2I_PARSE_SAMPLES : S2A_PARSE_SAMPLES;
            boolean checkStepParseSamples = stepSelections.get(parseStep) && parseStep.hasRequirements(proj, stepSelections, variables);
            return new boolean[][]{
                    {checkStepParseSamples,
                     (Files.exists(sampDir) && Files.list(sampDir, ".sampRAF", proj.JAR_STATUS.getValue()).length > 0)},
                    {stepSelections.get(S3_CREATE_SAMPLEDATA) && S3_CREATE_SAMPLEDATA.hasRequirements(proj, stepSelections, variables),
                     Files.exists(proj.SAMPLE_DATA_FILENAME.getValue(false, false))}, 
            };
        }
        
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            return new Object[]{proj.SAMPLE_DIRECTORY.getValue(false, false), proj.SAMPLE_DATA_FILENAME.getValue(false, false)};
        }

        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            return Files.exists(proj.SEXCHECK_RESULTS_FILENAME.getValue());
        }
        
        @Override
        public String getCommandLine(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String projPropFile = proj.getPropertyFilename();
            String kvCmd = "";
            String projDir = proj.SAMPLE_DIRECTORY.getValue(false, false);
            String setDir = variables.get(this).get(0);
            if (!ext.verifyDirFormat(setDir).equals(projDir)) {
                kvCmd += " SAMPLE_DIRECTORY=" + setDir;
            }
            StringBuilder cmd = new StringBuilder();
            if (kvCmd.length() > 0) {
                cmd.append("jcp cnv.filesys.Project proj=" + projPropFile).append(kvCmd).append("\n");
            }
            return cmd.append("jcp cnv.qc.SexChecks -check proj=" + projPropFile).toString();
        }
        
    };
    
    static final STEP S7_RUN_PLINK = new STEP("Create PLINK Files", 
                 "", 
                 new String[][]{{"[Parse Sample Files] step must be selected and valid.", "Parsed sample files must already exist."}, {"A pedigree.dat file is must exist.", "Create a minimal pedigree.dat file."}}, 
                 new RequirementInputType[][]{{RequirementInputType.NONE, RequirementInputType.DIR}, {RequirementInputType.FILE, RequirementInputType.BOOL}}) {

        @Override
        public void setNecessaryPreRunProperties(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String projDir = proj.SAMPLE_DIRECTORY.getValue(false, false);
            String setDir = variables.get(this).get(0);
            if (!ext.verifyDirFormat(setDir).equals(projDir)) {
                proj.SAMPLE_DIRECTORY.setValue(setDir);
            }
            if (!Boolean.valueOf(variables.get(this).get(2))) {
                String projPedFile = proj.PEDIGREE_FILENAME.getValue(false, false);
                String pedFile = variables.get(this).get(1);
                if (!pedFile.equals(projPedFile)) {
                    proj.PEDIGREE_FILENAME.setValue(pedFile);
                }
            }
        }
        
        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            if (Boolean.valueOf(variables.get(this).get(2))) {
                proj.getLog().report("Creating Pedigree File");
                Pedigree.build(proj, null, false);
            }
            if (!Files.exists(proj.PEDIGREE_FILENAME.getValue())) {
                setFailed();
                this.failReasons.add("Creation of Pedigree file in [Create/Run PLINK Files] step failed.");
                return;
            }
            
            proj.getLog().report("Running PLINK");
            
            boolean create = PlinkData.saveGenvisisToPlinkBedSet(proj, "plink/plink", null, null, -1, true);
            if (!create) {
                setFailed();
                this.failReasons.add("Creation of initial PLINK files failed.");
            }
        }
        
        @Override
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, Boolean> stepSelections, HashMap<STEP, ArrayList<String>> variables) {
            String sampDir = variables.get(this).get(0);
            String pedFile = variables.get(this).get(1);
            STEP parseStep = stepSelections.containsKey(S2I_PARSE_SAMPLES) ? S2I_PARSE_SAMPLES : S2A_PARSE_SAMPLES;
            boolean checkStepParseSamples = stepSelections.get(parseStep) && parseStep.hasRequirements(proj, stepSelections, variables);
            return new boolean[][]{
                    {checkStepParseSamples,
                        (Files.exists(sampDir) && Files.list(sampDir, ".sampRAF", proj.JAR_STATUS.getValue()).length > 0)},
                    {Files.exists(pedFile),
                            Boolean.parseBoolean(variables.get(this).get(2))}
            };
        }
        
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            return new Object[]{proj.SAMPLE_DIRECTORY.getValue(false, false), proj.PEDIGREE_FILENAME.getValue(false, false), false};
        }

        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variables) {
//            String fileCheck1 = proj.PROJECT_DIRECTORY.getValue()+"gwas.map";
            String fileCheck2 = proj.PROJECT_DIRECTORY.getValue()+"plink/plink.bed";
//            String fileCheck3 = proj.PROJECT_DIRECTORY.getValue()+"genome/";
            boolean pedCheck = Boolean.valueOf(variables.get(this).get(2)) ? Files.exists(proj.PEDIGREE_FILENAME.getValue()) : true;
            return/* Files.exists(fileCheck1) &&*/ Files.exists(fileCheck2) /*&& Files.exists(fileCheck3)*/ /*&& Files.list(fileCheck3, ".bed", false).length > 0*/ && pedCheck;
        }
        
        @Override
        public String getCommandLine(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String kvCmd = "";
            
            String projDir = proj.SAMPLE_DIRECTORY.getValue(false, false);
            String setDir = variables.get(this).get(0);
            if (!ext.verifyDirFormat(setDir).equals(projDir)) {
                kvCmd += " SAMPLE_DIRECTORY=" + setDir;
            }
            if (!Boolean.valueOf(variables.get(this).get(2))) {
                String projPedFile = proj.PEDIGREE_FILENAME.getValue(false, false);
                String pedFile = variables.get(this).get(1);
                if (!pedFile.equals(projPedFile)) {
                    kvCmd += " PEDIGREE_FILENAME=" + pedFile;
                }
            }
            
            String projPropFile = proj.getPropertyFilename();
            StringBuilder cmd = new StringBuilder();
            if (kvCmd.length() > 0) {
                cmd.append("jcp cnv.filesys.Project proj=").append(projPropFile).append(kvCmd).append("\n");
            }
            if (Boolean.valueOf(variables.get(this).get(2))) {
                cmd.append("jcp cnv.filesys.Pedigree proj=").append(projPropFile).append("\n");
            }
            cmd.append("jcp cnv.manage.PlinkData -genvisisToBed plinkdata=plink/plink gcthreshold=-1 proj=").append(proj.getPropertyFilename());
            return cmd.toString();
        }
    };
    
    static final STEP S8_GWAS_QC = new STEP("Run GWAS QC", 
               "", 
               new String[][]{{"[Create/Run PLINK Files] step must be selected and valid.", "PLINK files must already exist in the following directory."}, {"Keep genome info for unrelateds only?"}},
               new RequirementInputType[][]{{RequirementInputType.NONE, RequirementInputType.DIR}, {RequirementInputType.BOOL}}) {

        @Override
        public void setNecessaryPreRunProperties(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            // not needed for step
        }
        
        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String dir = getPlinkDir(proj, variables);
            boolean keepUnrelatedsOnly = Boolean.valueOf(variables.get(this).get(1));
            Qc.fullGamut(dir, null, keepUnrelatedsOnly, proj.getLog());
        }
        
        @Override
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, Boolean> stepSelections, HashMap<STEP, ArrayList<String>> variables) {
            String dir = getPlinkDir(proj, variables);
            boolean files = Files.exists(dir);
            return new boolean[][]{
                    {(stepSelections.get(S7_RUN_PLINK) && S7_RUN_PLINK.hasRequirements(proj, stepSelections, variables)), 
                        files},
                    {true}
            };
        }
        
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            return new Object[]{"plink/", ""};
        }

        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String dir = getPlinkDir(proj, variables);
            boolean allExist = true;
            folders: for (int i = 0; i < gwas.Qc.FOLDERS_CREATED.length; i++) {
                for (int j = 0; j < gwas.Qc.FILES_CREATED[i].length; j++) {
                    if (!Files.exists(dir + gwas.Qc.FOLDERS_CREATED[i] + gwas.Qc.FILES_CREATED[i][j])) {
                        allExist = false;
                        break folders;
                    }
                }
            }
            return allExist;
        }
        
        @Override
        public String getCommandLine(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String dir = getPlinkDir(proj, variables);
            boolean keepUnrelatedsOnly = Boolean.valueOf(variables.get(this).get(1));
            return "jcp gwas.Qc dir=" + dir + " keepGenomeInfoForRelatedsOnly=" + keepUnrelatedsOnly;
        }
        
        private String getPlinkDir(Project proj, HashMap<STEP, ArrayList<String>> variables){
        	// This directory is also used in S9_ANNOTATE_SAMPLE_DATA, any changes should be reflected there
        	String dir = variables.get(this).get(0);
            if (!dir.startsWith("/") && !dir.contains(":")) {
                dir = ext.verifyDirFormat(proj.PROJECT_DIRECTORY.getValue() + dir);
            }
            return dir;
        }
        
    };
    
    static final STEP S9_ANNOTATE_SAMPLE_DATA = new STEP("Annotate Sample Data File", 
            "", 
            new String[][]{{"Sample QC File must exist."},
    					   {"[" + S3_CREATE_SAMPLEDATA.stepName + "] step must be selected and valid.", "SampleData.txt file must already exist."},
    					   {"Skip identifying duplicates?","[" + S8_GWAS_QC.stepName + "] step must be selected and valid.","Duplicates Set file must already exist."},
    					   {"Do not use GC corrected LRR SD?","GC Corrected LRR SD must exist in Sample QC File"},
    					   {"LRR SD Threshold"},
    					   {"Callrate Threshold"},
    					   {"Number of Quantiles to Generate"},
    					   {"Replace FID and IID with data from Pedigree"}},
            new RequirementInputType[][]{{RequirementInputType.FILE},
    						   			 {RequirementInputType.NONE, RequirementInputType.FILE},
    						   			 {RequirementInputType.BOOL, RequirementInputType.NONE, RequirementInputType.FILE},
    						   			 {RequirementInputType.BOOL, RequirementInputType.NONE},
    						   			 {RequirementInputType.INT},
    						   			 {RequirementInputType.INT},
    						   			 {RequirementInputType.INT},
    						   			 {RequirementInputType.BOOL}}) {

     @Override
     public void setNecessaryPreRunProperties(Project proj, HashMap<STEP, ArrayList<String>> variables) {
    	 String projSampleQCFile = proj.SAMPLE_QC_FILENAME.getValue();
    	 String sampleQCFile = variables.get(this).get(0);
    	 String projSampleDataFile = proj.SAMPLE_DATA_FILENAME.getValue();
    	 String sampleDataFile = variables.get(this).get(1);
    	 double projLrrSdThreshold = proj.LRRSD_CUTOFF.getValue();
    	 double lrrSdThreshold = Double.parseDouble(variables.get(this).get(5));
    	 double projCallrateThreshold = proj.SAMPLE_CALLRATE_THRESHOLD.getValue();
    	 double callrateThreshold = Double.parseDouble(variables.get(this).get(6));
    	 
         if (!projSampleQCFile.equals(sampleQCFile) && Files.exists(sampleQCFile)) {
        	 proj.SAMPLE_QC_FILENAME.setValue(sampleQCFile);
         }
         if (!projSampleDataFile.equals(sampleDataFile) && Files.exists(sampleDataFile)) {
        	 proj.SAMPLE_DATA_FILENAME.setValue(sampleDataFile);
         }
         if (projLrrSdThreshold != lrrSdThreshold) {
        	 proj.LRRSD_CUTOFF.setValue(lrrSdThreshold);
         }
         if (projCallrateThreshold != callrateThreshold) {
        	 proj.SAMPLE_CALLRATE_THRESHOLD.setValue(callrateThreshold);
         }
     }
     
     @Override
     public void run(Project proj, HashMap<STEP, ArrayList<String>> variables) {
    	 boolean checkDuplicates = !Boolean.valueOf(variables.get(this).get(2));
         String duplicatesSetFile = null;
         if (checkDuplicates) {
        	 duplicatesSetFile = variables.get(this).get(3);
        	 if (!Files.exists(duplicatesSetFile)){
        		 String dir = variables.get(S8_GWAS_QC).get(0);
                 if (!dir.startsWith("/") && !dir.contains(":")) {
                     dir = ext.verifyDirFormat(proj.PROJECT_DIRECTORY.getValue() + dir);
                 }
        		 duplicatesSetFile = proj.PROJECT_DIRECTORY.getValue() + dir + "/quality_control/genome/plink.genome_duplicatesSet.dat";
        	 }
         }
         boolean gcCorrectedLrrSd = !Boolean.valueOf(variables.get(this).get(4));
         int numQ = Integer.parseInt(variables.get(this).get(7));
         boolean correctFidIids = Boolean.valueOf(variables.get(this).get(8));
         SampleQC.parseAndAddToSampleData(proj, numQ, 0, false, gcCorrectedLrrSd, duplicatesSetFile, correctFidIids);
     }
     
     @Override
     public boolean[][] checkRequirements(Project proj, HashMap<STEP, Boolean> stepSelections, HashMap<STEP, ArrayList<String>> variables) {
    	 String sampleQCFile = variables.get(this).get(0);
    	 boolean sampleQCFileExists = Files.exists(sampleQCFile);
    	 String sampleDataFile = variables.get(this).get(1);
    	 boolean checkDuplicates = !Boolean.valueOf(variables.get(this).get(2));
    	 String duplicatesSetFile = variables.get(this).get(3);
    	 boolean gcCorrectedLrrSd = !Boolean.valueOf(variables.get(this).get(4));
    	 boolean gcCorrectedLrrSdExists = false;
    	 if (sampleQCFileExists && ext.indexOfStr("LRR_SD_Post_Correction", Files.getHeaderOfFile(sampleQCFile, proj.getLog())) != -1) gcCorrectedLrrSdExists = true;
    	 double lrrSdThreshold = -1.0;
    	 try { lrrSdThreshold = Double.parseDouble(variables.get(this).get(5)); } catch (NumberFormatException nfe) {}
    	 double callrateThreshold = -1.0;
    	 try { callrateThreshold = Double.parseDouble(variables.get(this).get(6)); } catch (NumberFormatException nfe) {}
    	 int numQ = -1;
    	 try { numQ = Integer.parseInt(variables.get(this).get(7)); } catch (NumberFormatException nfe) {}
         return new boolean[][]{
                 {sampleQCFileExists},
                 {stepSelections.get(S3_CREATE_SAMPLEDATA) && S3_CREATE_SAMPLEDATA.hasRequirements(proj, stepSelections, variables),
                  Files.exists(sampleDataFile)},
                 {!checkDuplicates,
                  stepSelections.get(S8_GWAS_QC) && S8_GWAS_QC.hasRequirements(proj, stepSelections, variables),
                  Files.exists(duplicatesSetFile)},
                 {!gcCorrectedLrrSd, gcCorrectedLrrSdExists},
                 {lrrSdThreshold > proj.LRRSD_CUTOFF.getMinValue() && lrrSdThreshold < proj.LRRSD_CUTOFF.getMaxValue()},
                 {callrateThreshold > proj.SAMPLE_CALLRATE_THRESHOLD.getMinValue() && callrateThreshold < proj.SAMPLE_CALLRATE_THRESHOLD.getMaxValue()},
                 {numQ > 0},
                 {true}
         };
     }
     
     @Override
     public Object[] getRequirementDefaults(Project proj) {
    	 String sampleQCFile = proj.SAMPLE_QC_FILENAME.getValue();
         return new Object[]{sampleQCFile,
        		 			 proj.SAMPLE_DATA_FILENAME.getValue(),
        		 			 "false",
        		 			 proj.PROJECT_DIRECTORY.getValue() + "plink/quality_control/genome/plink.genome_duplicatesSet.dat",
        		 			 "false",
        		 			 proj.LRRSD_CUTOFF.getValue(),
        		 			 proj.SAMPLE_CALLRATE_THRESHOLD.getValue(),
        		 			 10,
        		 			 "false"};
     }

     @Override
     public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variables) {
    	 String sampleDataFile = variables.get(this).get(1);
    	 if (!Files.exists(sampleDataFile)) return false;
    	 boolean checkDuplicates = !Boolean.valueOf(variables.get(this).get(2));
    	 String[] header = Files.getHeaderOfFile(sampleDataFile, proj.getLog());
    	 if (checkDuplicates && ext.indexOfStr("DuplicateId", header, false, true) == -1) return false;
    	 if (ext.indexOfStr("Class=Exclude", header, false, true) == -1) return false;
    	 if (ext.indexOfStr("ExcludeNote", header, false, true) == -1) return false;
    	 if (ext.indexOfStr("Use", header, false, true) == -1) return false;
    	 if (ext.indexOfStr("UseNote", header, false, true) == -1) return false;
    	 if (ext.indexOfStr("Use_cnv", header, false, true) == -1) return false;
    	 if (ext.indexOfStr("Use_cnvNote", header, false, true) == -1) return false;
    	 return true;
     }
     
     @Override
     public String getCommandLine(Project proj, HashMap<STEP, ArrayList<String>> variables) {

    	 String projSampleQCFile = proj.SAMPLE_QC_FILENAME.getValue();
    	 String sampleQCFile = variables.get(this).get(0);
    	 String projSampleDataFile = proj.SAMPLE_DATA_FILENAME.getValue();
    	 String sampleDataFile = variables.get(this).get(1);
    	 double projLrrSdThreshold = proj.LRRSD_CUTOFF.getValue();
    	 double lrrSdThreshold = Double.parseDouble(variables.get(this).get(5));
    	 double projCallrateThreshold = proj.SAMPLE_CALLRATE_THRESHOLD.getValue();
    	 double callrateThreshold = Double.parseDouble(variables.get(this).get(6));

    	 String projPropFile = proj.getPropertyFilename();

    	 boolean checkDuplicates = !Boolean.valueOf(variables.get(this).get(2));
    	 String duplicatesSetFile = null;
    	 if (checkDuplicates) {
    		 duplicatesSetFile = variables.get(this).get(3);
    		 if (!Files.exists(duplicatesSetFile)){
    			 String dir = variables.get(S8_GWAS_QC).get(0);
    			 if (!dir.startsWith("/") && !dir.contains(":")) {
    				 dir = ext.verifyDirFormat(proj.PROJECT_DIRECTORY.getValue() + dir);
    			 }
    			 duplicatesSetFile = proj.PROJECT_DIRECTORY.getValue() + dir + "/quality_control/genome/plink.genome_duplicatesSet.dat";
    		 }
    	 }
    	 boolean gcCorrectedLrrSd = !Boolean.valueOf(variables.get(this).get(4));
    	 int numQ = Integer.parseInt(variables.get(this).get(7));
    	 boolean correctFidIids = Boolean.valueOf(variables.get(this).get(8));
         
    	 String kvCmd = "";
         
         if (!projSampleQCFile.equals(sampleQCFile)) {
             kvCmd += " SAMPLE_QC_FILENAME=" + sampleQCFile;
         }
         if (!projSampleDataFile.equals(sampleDataFile)) {
             kvCmd += " SAMPLE_DATA_FILENAME=" + sampleDataFile;
         }
         if (projLrrSdThreshold != lrrSdThreshold) {
             kvCmd += " LRRSD_CUTOFF=" + lrrSdThreshold;
         }
         if (projCallrateThreshold != callrateThreshold) {
             kvCmd += " SAMPLE_CALLRATE_THRESHOLD=" + callrateThreshold;
         }
         
         StringBuilder cmd = new StringBuilder();
         if (kvCmd.length() > 0) {
             cmd.append("jcp cnv.filesys.Project proj=" + projPropFile).append(kvCmd).append("\n");
         }
         cmd.append("jcp cnv.qc.SampleQC proj=" + projPropFile +
        		 					   " numQ=" + numQ +
        		 					   " justQuantiles=false" +
        							   " gcCorrectedLrrSd=" + gcCorrectedLrrSd +
        		 					   " duplicatesSetFile=" + duplicatesSetFile +
        		 					   " correctFidIids=" + correctFidIids);
         return cmd.toString();
     }
     
 };
    
    static final STEP S10_GENERATE_ABLOOKUP = new STEP("Generate AB Lookup File", "", 
            new String[][]{{"[Parse Sample Files] step must be selected and valid.", "A MarkerSet file must already exist."}, 
                           {"[Parse Sample Files] step must be selected and valid.", "Parsed sample files must already exist."}},
            new RequirementInputType[][]{{RequirementInputType.NONE, RequirementInputType.FILE}, {RequirementInputType.NONE, RequirementInputType.DIR}}) {

        @Override
        public void setNecessaryPreRunProperties(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String mkrPosProj = proj.MARKERSET_FILENAME.getValue(false, false);
            String mkrPosFile = variables.get(this).get(0);
            String setDir = proj.SAMPLE_DIRECTORY.getValue(false, false);
            String sampDir = variables.get(this).get(1);
            if (!mkrPosProj.equals(mkrPosFile)) {
                proj.MARKERSET_FILENAME.setValue(mkrPosFile);
            }
            if (!ext.verifyDirFormat(setDir).equals(sampDir)) {
                proj.SAMPLE_DIRECTORY.setValue(sampDir);
            }
        }
        
        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            ABLookup abLookup;
            String filename;
            
            filename = proj.PROJECT_DIRECTORY.getValue()+ext.addToRoot(ABLookup.DEFAULT_AB_FILE, "_parsed");
            if (!Files.exists(filename)) {
                abLookup = new ABLookup();
                abLookup.parseFromAnnotationVCF(proj);
                abLookup.writeToFile(filename, proj.getLog());
            }
            if (Files.exists(filename)) {
                ABLookup.fillInMissingAlleles(proj, filename, proj.getLocationOfSNP_Map(true), false);
                ABLookup.applyABLookupToFullSampleFiles(proj);
                proj.AB_LOOKUP_FILENAME.setValue(filename);
                proj.saveProperties(new Project.Property[]{proj.AB_LOOKUP_FILENAME});
            } else {
                setFailed();
            }
        }
        
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            return new Object[]{proj.MARKERSET_FILENAME.getValue(false, false), proj.SAMPLE_DIRECTORY.getValue(false, false)};
        }
        
        @Override
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, Boolean> stepSelections, HashMap<STEP, ArrayList<String>> variables) {
            String mkrPosFile = variables.get(this).get(0);
            String sampDir = variables.get(this).get(1);
            boolean checkStepParseSamples = stepSelections.get(S2I_PARSE_SAMPLES) && S2I_PARSE_SAMPLES.hasRequirements(proj, stepSelections, variables);
            return new boolean[][]{
                    {checkStepParseSamples, 
                        Files.exists(mkrPosFile)},
                    {checkStepParseSamples,
                        (Files.exists(sampDir) && Files.list(sampDir, ".sampRAF", proj.JAR_STATUS.getValue()).length > 0)},
            };
        }
        
        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            return Files.exists(proj.AB_LOOKUP_FILENAME.getValue(false, false));
        }
        @Override
        public String getCommandLine(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            return "## << Generate ABLookup >> Not Implemented For Command Line Yet ##"; // TODO
        }
    };

    static final STEP S11_COMPUTE_PFB = new STEP("Compute Population BAF files", "", new String[][]{
            {"[Parse Sample Files] step must be selected and valid.", "A SampleList file must already exist.", "A Sample subset file must exist."}, 
            {"PFB (population BAF) output file must be specified."}},
            new RequirementInputType[][]{
                    {RequirementInputType.NONE, RequirementInputType.FILE, RequirementInputType.FILE},
                    {RequirementInputType.FILE}}
            ) {
    
        @Override
        public void setNecessaryPreRunProperties(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String setSampList = proj.SAMPLELIST_FILENAME.getValue();
            String sampListFile = variables.get(this).get(0);
            String setSubSampFile = proj.SAMPLE_SUBSET_FILENAME.getValue();
            String subSampFile = variables.get(this).get(1);
            String setPFBFile = proj.CUSTOM_PFB_FILENAME.getValue();
            String pfbOutputFile = variables.get(this).get(2);
            
            if (!ext.verifyDirFormat(setSampList).equals(sampListFile)) {
                proj.SAMPLELIST_FILENAME.setValue(sampListFile);
            }
            if (!ext.verifyDirFormat(setSubSampFile).equals(subSampFile)) {
                proj.SAMPLE_SUBSET_FILENAME.setValue(subSampFile);
            }
            if (!ext.verifyDirFormat(setPFBFile).equals(pfbOutputFile)) {
                proj.CUSTOM_PFB_FILENAME.setValue(pfbOutputFile);
            }
        }
        
        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            cnv.analysis.PennCNV.populationBAF(proj);
        }
        
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            return new Object[]{proj.SAMPLELIST_FILENAME.getValue(), proj.SAMPLE_SUBSET_FILENAME.getValue(), Files.exists(proj.SAMPLE_SUBSET_FILENAME.getValue()) ? ext.rootOf(proj.SAMPLE_SUBSET_FILENAME.getValue()) + ".pfb" : proj.CUSTOM_PFB_FILENAME.getValue()};
        }
        
        @Override
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, Boolean> stepSelections, HashMap<STEP, ArrayList<String>> variables) {
            String sampListFile = variables.get(this).get(0);
            String subSampFile = variables.get(this).get(1);
            String pfbOutputFile = variables.get(this).get(2);
            
            STEP parseStep = stepSelections.containsKey(S2I_PARSE_SAMPLES) ? S2I_PARSE_SAMPLES : S2A_PARSE_SAMPLES;
            boolean checkStepParseSamples = stepSelections.get(parseStep) && parseStep.hasRequirements(proj, stepSelections, variables);
            boolean step12 = Files.exists(sampListFile);
            boolean step13 = Files.exists(subSampFile);
            boolean step21 = !Files.exists(pfbOutputFile);
            
            return new boolean[][]{{checkStepParseSamples, step12, step13}, {step21}};
        }
        
        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String subSampFile = variables.get(this).get(1);
            String pfbOutputFile = variables.get(this).get(2);
            boolean pfbExists = Files.exists(pfbOutputFile) || Files.exists(ext.rootOf(subSampFile) + ".pfb");
            return pfbExists;
        }
        
        @Override
        public String getCommandLine(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String kvCmd = "";
            
            String setSampList = proj.SAMPLELIST_FILENAME.getValue();
            String sampListFile = variables.get(this).get(0);
            String setSubSampFile = proj.SAMPLE_SUBSET_FILENAME.getValue();
            String subSampFile = variables.get(this).get(1);
            String setPFBFile = proj.CUSTOM_PFB_FILENAME.getValue();
            String pfbOutputFile = variables.get(this).get(2);
            
            if (!ext.verifyDirFormat(setSampList).equals(sampListFile)) {
                kvCmd += " SAMPLELIST_FILENAME=" + sampListFile;
            }
            if (!ext.verifyDirFormat(setSubSampFile).equals(subSampFile)) {
                kvCmd += " SAMPLE_SUBSET_FILENAME=" + subSampFile;
            }
            if (!ext.verifyDirFormat(setPFBFile).equals(pfbOutputFile)) {
                kvCmd += " CUSTOM_PFB_FILENAME=" + pfbOutputFile;
            }
            
            String projPropFile = proj.getPropertyFilename();
            StringBuilder cmd = new StringBuilder();
            if (kvCmd.length() > 0) {
                cmd.append("jcp cnv.filesys.Project proj=" + projPropFile).append(kvCmd).append("\n");
            }
            return cmd.append("jcp cnv.analysis.PennCNV -pfb proj=" + proj.getPropertyFilename() + " log=" + proj.getLog().getFilename()).toString();
        }
        
    };
    static final STEP S12_COMPUTE_GCMODEL = new STEP("Compute GCMODEL File", 
                    "", 
                    new String[][]{
                                {"A GC Base file must exist."}, 
                                {"GCModel output file must be specified."}},
                    new RequirementInputType[][]{
                            {RequirementInputType.FILE}, 
                            {RequirementInputType.FILE}
                }) {
    
        @Override
        public void setNecessaryPreRunProperties(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String setGCOutputFile = proj.GC_MODEL_FILENAME.getValue();
            String gcOutputFile = variables.get(this).get(1);
            if (!ext.verifyDirFormat(setGCOutputFile).equals(gcOutputFile)) {
                proj.GC_MODEL_FILENAME.setValue(gcOutputFile);
            }
        }
        
        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String gcBaseFile = variables.get(this).get(0);
            String gcOutputFile = variables.get(this).get(1);
            cnv.analysis.PennCNV.gcModel(proj, gcBaseFile, gcOutputFile, 100);
        }
        @Override
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, Boolean> stepSelections, HashMap<STEP, ArrayList<String>> variables) {
            String gcBaseFile = variables.get(this).get(0);
            String gcOutputFile = variables.get(this).get(1);
            return new boolean[][]{{Files.exists(gcBaseFile)},{!Files.exists(gcOutputFile)}};
        }
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            return new Object[]{Files.firstPathToFileThatExists(Aliases.REFERENCE_FOLDERS, "gc5Base.txt", true, false, proj.getLog()), proj.GC_MODEL_FILENAME.getValue()};
        }
        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String gcOutputFile = variables.get(this).get(1);
            boolean gcExists = Files.exists(gcOutputFile);
            return gcExists;
        }
        
        @Override
        public String getCommandLine(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String kvCmd = "";
            
            String setGCOutputFile = proj.GC_MODEL_FILENAME.getValue();
            String gcOutputFile = variables.get(this).get(1);
            if (!ext.verifyDirFormat(setGCOutputFile).equals(gcOutputFile)) {
                kvCmd += " GC_MODEL_FILENAME=" + gcOutputFile;
            }
            
            String projPropFile = proj.getPropertyFilename();
            StringBuilder cmd = new StringBuilder();
            if (kvCmd.length() > 0) {
                cmd.append("jcp cnv.filesys.Project proj=" + projPropFile).append(kvCmd).append("\n");
            }
            String gcBaseFile = variables.get(this).get(0);
            return cmd.append("jcp cnv.analysis.PennCNV proj=" + proj.getPropertyFilename() + " log=" + proj.getLog().getFilename() + " gc5base=" + gcBaseFile).toString();
        }
        
    };
    static final STEP S13_CREATE_PCS = new STEP("Create Principal Components File and Mitochondrial Copy-Number Estimates File", 
                  "", 
                  new String[][]{
                            {"[Transpose Data into Marker-Dominant Files] step must be selected and valid.", "Parsed marker data files must already exist."}, 
                            {"MedianMarkers file must exist."}, 
                            {"FASTA Reference Genome file"},
                            {"Compute PCs with samples passing QC only?"},
                            {"Should impute mean value for NaN?"}, 
                            {"Should recompute Log-R ratio for PC markers?"}, 
                            {"Should recompute Log-R ratio for median markers?"}, 
                            {"Homozygous only?"}, 
                            {"Base-pair bins for the GC model generated from the reference"},
                            {"Regression distance for the GC adjustment"},
                            {"Number of principal components."}, 
                            {"Number of threads to use"},
                            },
                  new RequirementInputType[][]{
                            {RequirementInputType.NONE, RequirementInputType.DIR}, 
                            {RequirementInputType.FILE},
                            {RequirementInputType.FILE},
                            {RequirementInputType.BOOL}, 
                            {RequirementInputType.BOOL}, 
                            {RequirementInputType.BOOL}, 
                            {RequirementInputType.BOOL}, 
                            {RequirementInputType.BOOL}, 
                            {RequirementInputType.INT}, 
                            {RequirementInputType.INT}, 
                            {RequirementInputType.INT}, 
                            {RequirementInputType.INT}, 
                    }) {

        @Override
        public void setNecessaryPreRunProperties(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            // not needed for step
        }
        
        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String medianMarkers = variables.get(this).get(1);
            String refGenomeFasta = variables.get(this).get(2);
            boolean gcCorrect = Boolean.valueOf(variables.get(this).get(3));
            boolean imputeMeanForNaN = Boolean.valueOf(variables.get(this).get(4));
            boolean recomputeLRR_PCs = Boolean.valueOf(variables.get(this).get(5));
            boolean recomputeLRR_Median = Boolean.valueOf(variables.get(this).get(6));
            boolean homozygousOnly = Boolean.valueOf(variables.get(this).get(7));
            int bpGcModel = Integer.parseInt(variables.get(this).get(8));
            int regressionDistance = Integer.parseInt(variables.get(this).get(9));
            int numComponents = Integer.parseInt(variables.get(this).get(10));
            int numThreads = Integer.parseInt(variables.get(this).get(11));
            String outputBase = proj.PROJECT_DIRECTORY.getValue() + MitoPipeline.FILE_BASE;

            GcAdjustorParameters params = null;
            if (gcCorrect) {
                if ((refGenomeFasta != null && !Files.exists(refGenomeFasta)) && Files.exists(proj.REFERENCE_GENOME_FASTA_FILENAME.getValue())) {
                    proj.getLog().reportTimeWarning("Command line reference genome did not exist or was not provided, using default " + proj.REFERENCE_GENOME_FASTA_FILENAME.getValue());
                    refGenomeFasta = proj.REFERENCE_GENOME_FASTA_FILENAME.getValue();
                }
                if (Files.exists(refGenomeFasta) || Files.exists(proj.GC_MODEL_FILENAME.getValue())) {// TODO, after evaluating reference genome based gc model files, will demand a refGenome
                    if (refGenomeFasta != null && Files.exists(refGenomeFasta)) {
                        proj.REFERENCE_GENOME_FASTA_FILENAME.setValue(refGenomeFasta);
                    }
                    GCAdjustorBuilder gAdjustorBuilder = new GCAdjustorBuilder();
                    gAdjustorBuilder.regressionDistance(regressionDistance);
                    params = GcAdjustorParameter.generate(proj, "GC_ADJUSTMENT/", refGenomeFasta, gAdjustorBuilder, recomputeLRR_Median || recomputeLRR_PCs, bpGcModel, numThreads);
                    if ((recomputeLRR_Median || recomputeLRR_PCs) && params.getCentroids() == null) {
                        throw new IllegalStateException("Internal error, did not recieve centroids");
                    } else if ((!recomputeLRR_Median && !recomputeLRR_PCs) && params.getCentroids() != null) {
                        throw new IllegalStateException("Internal error, should not have recieved centroids");
                    }
                    recomputeLRR_Median = false;// recomputed if params has centroid
                    recomputeLRR_PCs = false;
                } else {
                    proj.getLog().reportTimeError("Can not gc correct values without a valid reference genome");
                    proj.getLog().reportTimeError("please supply a valid reference genome (full path) with the \"ref=\" argument");
                }
            }
            proj.getLog().report("\nReady to perform the principal components analysis (PCA)\n");
            PrincipalComponentsCompute pcs = PCA.computePrincipalComponents(proj, false, numComponents, false, false, true, true, imputeMeanForNaN, recomputeLRR_PCs, outputBase + MitoPipeline.PCA_SAMPLES, outputBase, params);
            if (pcs == null) {
                setFailed();
                this.failReasons.add("# of Principal Components is greater than either the # of samples or the # of markers.  Please lower the # of PCs and try again.");
                return;
            }
            // apply PCs to everyone, we set useFile to null and excludeSamples to false to get all samples in the current project.
            // TODO, if we ever want to apply to only a subset of the project, we can do that here.....
            proj.getLog().report("\nApplying the loadings from the principal components analysis to all samples\n");
            PrincipalComponentsApply pcApply = PCA.applyLoadings(proj, numComponents, pcs.getSingularValuesFile(), pcs.getMarkerLoadingFile(), null, false, imputeMeanForNaN, recomputeLRR_PCs, outputBase, params);
            // Compute Medians for (MT) markers and compute residuals from PCs for everyone
            proj.getLog().report("\nComputing residuals after regressing out " + numComponents + " principal component" + (numComponents == 1 ? "" : "s") + "\n");
            PrincipalComponentsResiduals pcResids = PCA.computeResiduals(proj, pcApply.getExtrapolatedPCsFile(), ext.removeDirectoryInfo(medianMarkers), numComponents, true, 0f, homozygousOnly, recomputeLRR_Median, outputBase, params);
            MitoPipeline.generateFinalReport(proj, outputBase, pcResids.getResidOutput());
            proj.setProperty(proj.INTENSITY_PC_FILENAME, pcApply.getExtrapolatedPCsFile());
            proj.setProperty(proj.INTENSITY_PC_NUM_COMPONENTS, numComponents);
            proj.saveProperties();
        }

        @Override
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, Boolean> stepSelections, HashMap<STEP, ArrayList<String>> variables) {
            String markerDir = variables.get(this).get(0);
            String medianMkrs = variables.get(this).get(1);
            String fastaFile = variables.get(this).get(2);
            int numComponents = -1;
            int bpGcModel = -1;
            int regressionDistance = -1;
            int numThreads = -1;
            try { bpGcModel = Integer.parseInt(variables.get(this).get(8)); } catch (NumberFormatException e) {}
            try { regressionDistance = Integer.parseInt(variables.get(this).get(9)); } catch (NumberFormatException e) {}
            try { numComponents = Integer.parseInt(variables.get(this).get(10)); } catch (NumberFormatException e) {}
            try { numThreads = Integer.parseInt(variables.get(this).get(11)); } catch (NumberFormatException e) {}
            return new boolean[][]{
                    {(stepSelections.get(S4_TRANSPOSE_TO_MDF) && S4_TRANSPOSE_TO_MDF.hasRequirements(proj, stepSelections, variables)), Files.exists(markerDir)},
                    {Files.exists(medianMkrs)},
                    {Files.exists(fastaFile)},
                    {true}, // TRUE or FALSE are both valid selections
                    {true}, 
                    {true}, 
                    {true}, 
                    {true}, 
                    {bpGcModel > 0},
                    {regressionDistance > 0},
                    {numComponents > 0},
                    {numThreads > 0},
            };
        }
        
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            return new Object[]{
                    proj.MARKER_DATA_DIRECTORY.getValue(false, false),
                    "",
                    proj.REFERENCE_GENOME_FASTA_FILENAME.getValue(false, false),
                    "true", 
                    "true", 
                    "true", 
                    "true", 
                    "true",
                    GcAdjustor.GcModel.DEFAULT_GC_MODEL_BIN_FASTA,
                    GcAdjustor.DEFAULT_REGRESSION_DISTANCE[0],
                    proj.INTENSITY_PC_NUM_COMPONENTS.getValue().toString(), 
                    proj.NUM_THREADS.getValue()/*,""*/
            };
        }

        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String outputBase = proj.PROJECT_DIRECTORY.getValue() + MitoPipeline.FILE_BASE;
            String finalReport = outputBase + PCA.FILE_EXTs[0];//PrincipalComponentsResiduals.MT_REPORT_EXT[0];
//            boolean mkrFiles = true;
//            for (String file : PrincipalComponentsResiduals.MT_REPORT_MARKERS_USED) {
//                if (!Files.exists(outputBase + file)) {
//                    mkrFiles = false;
//                    break;
//                }
//            }
            return Files.exists(finalReport) /*&& mkrFiles*/;
        }
        @Override
        public String getCommandLine(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            return "## << Create PrincipalComponents File >> Not Implemented For Command Line Yet ##"; // TODO
        }
        
    };
//    
//    static final STEP S14_CREATE_MT_CN_EST = new STEP("Create Mitochondrial Copy-Number Estimates File", 
//                        "", 
//                        new String[][]{
//                                    {"[Transpose Data into Marker-Dominant Files] step must be selected and valid.", "Parsed marker data files must already exist."}, 
//                                    {"[Create Principal Components File] step must be selected and valid.", "Extrapolated PCs file must already exist."},
//                                    {"MedianMarkers file must exist."}, 
//                                    {"FASTA Reference Genome file"},
//                                    {"Number of principal components."}, 
//                                    {"Should recompute Log-R ratio median?"}, 
//                                    {"Should recompute Log-R ratio median?"}, 
//                                    {"Homozygous only?"}, 
//                                    {"GC correct?"},
//                                    {"Regression distance for the GC adjustment"},
//                                    {"Number of threads to use"},
//                                    
//                                },
//                        new RequirementInputType[][]{
//                                    {RequirementInputType.NONE, RequirementInputType.DIR},
//                                    {RequirementInputType.NONE, RequirementInputType.FILE},
//                                    {RequirementInputType.FILE}, 
//                                    {RequirementInputType.FILE}, 
//                                    {RequirementInputType.INT}, 
//                                    {RequirementInputType.BOOL}, 
//                                    {RequirementInputType.BOOL}, 
//                                    {RequirementInputType.BOOL}, 
//                                    {RequirementInputType.BOOL}, 
//                                    {RequirementInputType.INT}, 
//                                    {RequirementInputType.INT}, 
//                        }) {
//
//        @Override
//        public void setNecessaryPreRunProperties(Project proj, HashMap<STEP, ArrayList<String>> variables) {
//            // not needed for step
//        }
//        
//        @Override
//        public void run(Project proj, HashMap<STEP, ArrayList<String>> variables) {
//            String medianMarkers = variables.get(this).get(2);
//            String refGenomeFasta = variables.get(this).get(3);
//            int numComponents = Integer.parseInt(variables.get(this).get(4));
//            boolean recomputeLRR_Median = Boolean.valueOf(variables.get(this).get(5));
//            boolean recomputeLRR_PCs = Boolean.valueOf(variables.get(this).get(6));
//            boolean homozygousOnly = Boolean.valueOf(variables.get(this).get(7));
//            boolean gcCorrect = Boolean.valueOf(variables.get(this).get(8));
//            int regressionDistance = Integer.parseInt(variables.get(this).get(9));
//            int numThreads = Integer.parseInt(variables.get(this).get(10));
//            
//            String outputBase = proj.PROJECT_DIRECTORY.getValue() + MitoPipeline.FILE_BASE;//variables.get(this).get(6);
//            String extrapolatedPCsFile = MitoPipeline.FILE_BASE + PCA.FILE_EXTs[0];
//            
//            
//            GcAdjustorParameters params = null;
//            if (gcCorrect) {
//                if ((refGenomeFasta != null && !Files.exists(refGenomeFasta)) && Files.exists(proj.REFERENCE_GENOME_FASTA_FILENAME.getValue())) {
//                    proj.getLog().reportTimeWarning("Command line reference genome did not exist or was not provided, using default " + proj.REFERENCE_GENOME_FASTA_FILENAME.getValue());
//                    refGenomeFasta = proj.REFERENCE_GENOME_FASTA_FILENAME.getValue();
//                }
//                if (Files.exists(refGenomeFasta) || Files.exists(proj.GC_MODEL_FILENAME.getValue())) {// TODO, after evaluating reference genome based gc model files, will demand a refGenome
//                    if (refGenomeFasta != null && Files.exists(refGenomeFasta)) {
//                        proj.REFERENCE_GENOME_FASTA_FILENAME.setValue(refGenomeFasta);
//                    }
//                    // try {
//                    GCAdjustorBuilder gAdjustorBuilder = new GCAdjustorBuilder();
//                    gAdjustorBuilder.regressionDistance(regressionDistance);
//                    params = GcAdjustorParameter.generate(proj, "GC_ADJUSTMENT/", refGenomeFasta, gAdjustorBuilder, recomputeLRR_Median || recomputeLRR_PCs, bpGcModel, numThreads);
//                    if ((recomputeLRR_Median || recomputeLRR_PCs) && params.getCentroids() == null) {
//                        throw new IllegalStateException("Internal error, did not recieve centroids");
//                    } else if ((!recomputeLRR_Median && !recomputeLRR_PCs) && params.getCentroids() != null) {
//                        throw new IllegalStateException("Internal error, should not have recieved centroids");
//                    }
//                    recomputeLRR_Median = false;// recomputed if params has centroid
//                    recomputeLRR_PCs = false;
//                    // } catch (IllegalStateException e) {
//                    //
//                    // proj.getLog().reportTimeError("GC adjustment was flagged, but could not generate neccesary files");
//                    // }
//                } else {
//                    proj.getLog().reportTimeError("Can not gc correct values without a valid reference genome");
//                    proj.getLog().reportTimeError("please supply a valid reference genome (full path) with the \"ref=\" argument");
//                }
//            }
//
//            proj.getLog().report("\nComputing residuals after regressing out " + numComponents + " principal component" + (numComponents == 1 ? "" : "s") + "\n");
//            PrincipalComponentsResiduals pcResids = PCA.computeResiduals(proj, extrapolatedPCsFile, ext.removeDirectoryInfo(medianMarkers), numComponents, true, 0f, homozygousOnly, recomputeLRR_Median, outputBase, params);
//            MitoPipeline.generateFinalReport(proj, outputBase, pcResids.getResidOutput());
//            proj.setProperty(proj.INTENSITY_PC_FILENAME, extrapolatedPCsFile);
//            proj.setProperty(proj.INTENSITY_PC_NUM_COMPONENTS, numComponents);
//            proj.saveProperties(new Project.Property[]{proj.INTENSITY_PC_FILENAME});
//        }
//        @Override
//        public boolean[][] checkRequirements(Project proj, HashMap<STEP, Boolean> stepSelections, HashMap<STEP, ArrayList<String>> variables) {
//            String markerDir = variables.get(this).get(0);
//            String extrapPCFile = variables.get(this).get(1);
//            String medianMarkers = variables.get(this).get(2);
//            int numComponents = -1;
//            try {
//                numComponents = Integer.parseInt(variables.get(this).get(3));
//            } catch (NumberFormatException e) {}
//            return new boolean[][]{
//                    {(stepSelections.get(S4_TRANSPOSE_TO_MDF) && S4_TRANSPOSE_TO_MDF.hasRequirements(proj, stepSelections, variables)), Files.exists(markerDir)},
//                    {(stepSelections.get(S13_CREATE_PCS) && S13_CREATE_PCS.hasRequirements(proj, stepSelections, variables)), Files.exists(extrapPCFile)},
//                    {Files.exists(medianMarkers)},
//                    {numComponents > 0},
//                    {true}, // TRUE or FALSE are both valid selections
//                    {true}, 
//            };
//        }
//        @Override
//        public Object[] getRequirementDefaults(Project proj) {
//            return new String[]{proj.MARKER_DATA_DIRECTORY.getValue(false, false), proj.INTENSITY_PC_FILENAME.getValue(), "", proj.INTENSITY_PC_NUM_COMPONENTS.getValue().toString(), "", ""/*, ""*/};
//        }
//        @Override
//        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variables) {
//            // existing files will be backed up if re-run
//            return false;
//        }
//        @Override
//        public String getCommandLine(Project proj, HashMap<STEP, ArrayList<String>> variables) {
//            return "## << Create Mitochondrial Copy-Number Estimates File >> Not Implemented For Command Line Yet ##"; // TODO
//        }
//    };
    
    static final STEP S15_MOSAIC_ARMS = new STEP("Create Mosaic Arms File", 
                                                 "", 
                                                 new String[][]{
                                                        {"[Parse Sample Files] step must be selected and valid.","A MarkerSet file must already exist."}, 
                                                        {"Number of threads to use."}},
                                                 new RequirementInputType[][]{
                                                        {RequirementInputType.NONE, RequirementInputType.FILE}, 
                                                        {RequirementInputType.INT}}) {

        @Override
        public void setNecessaryPreRunProperties(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String mkrPosProj = proj.MARKERSET_FILENAME.getValue(false, false);
            String mkrPosFile = variables.get(this).get(0);
            if (!mkrPosProj.equals(mkrPosFile)) {
                proj.MARKERSET_FILENAME.setValue(mkrPosFile);
            }
            int numThreads = proj.NUM_THREADS.getValue();
            try {
                numThreads = Integer.parseInt(variables.get(this).get(1));
            } catch (NumberFormatException e) {}
            if (numThreads != proj.NUM_THREADS.getValue()) {
                proj.NUM_THREADS.setValue(numThreads);
            }
        }
        
        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            Mosaicism.findOutliers(proj);
        }
        
        @Override
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, Boolean> stepSelections, HashMap<STEP, ArrayList<String>> variables) {
            boolean checkStepParseSamples = stepSelections.get(S2I_PARSE_SAMPLES) && S2I_PARSE_SAMPLES.hasRequirements(proj, stepSelections, variables);
            String mkrPosFile = variables.get(this).get(0);
            boolean step11 = Files.exists(mkrPosFile);
            int numThreads = -1;
            try {
                numThreads = Integer.parseInt(variables.get(this).get(1));
            } catch (NumberFormatException e) {}
            return new boolean[][]{{checkStepParseSamples, step11}, {numThreads > 0}};
        }
        
        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            boolean outputCheck = Files.exists(proj.RESULTS_DIRECTORY.getValue(false, false) + "Mosaicism.xln");
            return outputCheck;
        }
        
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            return new Object[]{proj.MARKERSET_FILENAME.getValue(), proj.NUM_THREADS.getValue()};
        }
        
        @Override
        public String getCommandLine(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String kvCmd = "";
            
            String mkrPosProj = proj.MARKERSET_FILENAME.getValue(false, false);
            String mkrPosFile = variables.get(this).get(0);
            if (!mkrPosProj.equals(mkrPosFile)) {
                kvCmd += " MARKERSET_FILENAME=" + mkrPosFile;
            }
            int numThreads = proj.NUM_THREADS.getValue();
            try {
                numThreads = Integer.parseInt(variables.get(this).get(1));
            } catch (NumberFormatException e) {}
            if (numThreads != proj.NUM_THREADS.getValue()) {
                kvCmd += " NUM_THREADS=" + numThreads;
            }
            
            String projPropFile = proj.getPropertyFilename();
            StringBuilder cmd = new StringBuilder();
            if (kvCmd.length() > 0) {
                cmd.append("jcp cnv.filesys.Project proj=" + projPropFile).append(kvCmd).append("\n");
            }
            return cmd.append("jcp cnv.analysis.Mosaicism proj=" + proj.getPropertyFilename()).toString();
        }
        
    };
    
    static final STEP S16_SEX_CENTROIDS_PFB_GCMODEL = new STEP("Create Sex-Specific Centroids; Filter PFB and GCMODEL Files", 
                            "", 
                            new String[][]{
                                {"Full GC Model File."},
                                {"Number of threads to use."},
                            },
                            new RequirementInputType[][]{
                                {RequirementInputType.FILE},
                                {RequirementInputType.INT}
                            }) {
        @Override
        public void setNecessaryPreRunProperties(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            int numThreads = proj.NUM_THREADS.getValue();
            try {
                numThreads = Integer.parseInt(variables.get(this).get(1));
            } catch (NumberFormatException e) {}
            if (numThreads != proj.NUM_THREADS.getValue()) {
                proj.NUM_THREADS.setValue(numThreads);
            }
        }
        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String pennData, sexDir, malePFB, femalePFB, centFilePathM, centFilePathF, newGCFile;
            pennData = proj.getProperty(proj.PENNCNV_DATA_DIRECTORY);
            sexDir = pennData + "sexSpecific/";
            newGCFile = sexDir + "sexSpecific.gcModel";
            malePFB = sexDir + "males.pfb";
            femalePFB = sexDir + "females.pfb";
            centFilePathM = sexDir + "sexSpecific_Male.cent";
            centFilePathF = sexDir + "sexSpecific_Female.cent";
            int numThreads = proj.NUM_THREADS.getValue();
            try {
                numThreads = Integer.parseInt(variables.get(this).get(1));
            } catch (NumberFormatException e) {}
            String gcModelFile = variables.get(this).get(0);
            Centroids.computeSexSpecificCentroids(proj, AnalysisFormats.getChromosomalMarkersOnly(proj), new String[]{malePFB, femalePFB}, new String[]{centFilePathM, centFilePathF}, true, numThreads);

            AnalysisFormats.filterSexSpecificGCModel(proj, gcModelFile, newGCFile);
        }
        @Override
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, Boolean> stepSelections, HashMap<STEP, ArrayList<String>> variables) {
            int numThreads = -1;
            try {
                numThreads = Integer.parseInt(variables.get(this).get(1));
            } catch (NumberFormatException e) {}
            String gcModelFile = variables.get(this).get(0);
            return new boolean[][]{
                {Files.exists(gcModelFile)},
                {numThreads > 0}
            };
        }
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            return new Object[]{proj.GC_MODEL_FILENAME.getValue(), proj.NUM_THREADS.getValue()};
        }
        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String pennData, sexDir, malePFB, femalePFB, centFilePathM, centFilePathF, newGCFile;
            pennData = proj.getProperty(proj.PENNCNV_DATA_DIRECTORY);
            sexDir = pennData + "sexSpecific/";
            malePFB = sexDir + "males.pfb";
            femalePFB = sexDir + "females.pfb";
            centFilePathM = sexDir + "sexSpecific_Male.cent";
            centFilePathF = sexDir + "sexSpecific_Female.cent";
            newGCFile = sexDir + "sexSpecific.gcModel";
            return Files.exists(malePFB) && Files.exists(femalePFB) && Files.exists(centFilePathM) && Files.exists(centFilePathF) && Files.exists(newGCFile);
        }
        @Override
        public String getCommandLine(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            int numThreads = proj.NUM_THREADS.getValue();
            try {
                numThreads = Integer.parseInt(variables.get(this).get(1));
            } catch (NumberFormatException e) {}
            String mainCmd = "jcp cnv.filesys.Centroids proj=" + proj.getPropertyFilename() + " threads=" + numThreads;
            String gcModelFile = variables.get(this).get(0);
            String gcCmd = "jcp cnv.analysis.AnalysisFormats proj=" + proj.getPropertyFilename() + " gcmodel=" + gcModelFile;
            return mainCmd + "\n" + gcCmd; 
        }
        
    };
    
    static final STEP S17_CNV_CALLING = new STEP("Call CNVs", 
                            "", 
                            new String[][]{
                                {"Hidden Markov Model File Must Exist"},
                                {"[Compute Population BAF File] step must be selected and valid", "PFB File Must Exist"},
                                {"[Compute GCMODEL File] step must be selected and valid", "GCMODEL File Must Exist"},
                                {"Number of threads To use."},
                                {"Output filename."}
                            },
                            new RequirementInputType[][]{
                                {RequirementInputType.FILE},
                                {RequirementInputType.NONE, RequirementInputType.FILE},
                                {RequirementInputType.NONE, RequirementInputType.FILE},
                                {RequirementInputType.INT},
                                {RequirementInputType.FILE},
                            }) {
        @Override
        public void setNecessaryPreRunProperties(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String hmm_P = proj.HMM_FILENAME.getValue();
            String hmm_G = variables.get(this).get(0);
            if (!hmm_P.equals(hmm_G)) {
                proj.HMM_FILENAME.setValue(hmm_G);
            }
            String pfb_P = proj.CUSTOM_PFB_FILENAME.getValue();
            String pfb_G = variables.get(this).get(1);
            if (!pfb_P.equals(pfb_G)) {
                proj.CUSTOM_PFB_FILENAME.setValue(pfb_G);
            }
            String gcm_P = proj.GC_MODEL_FILENAME.getValue();
            String gcm_G = variables.get(this).get(2);
            if (!gcm_P.equals(gcm_G)) {
                proj.GC_MODEL_FILENAME.setValue(gcm_G);
            }
            int numThreads = proj.NUM_THREADS.getValue();
            try {
                numThreads = Integer.parseInt(variables.get(this).get(3));
            } catch (NumberFormatException e) {}
            if (numThreads != proj.NUM_THREADS.getValue()) {
                proj.NUM_THREADS.setValue(numThreads);
            }
        }
        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            int numThreads = proj.NUM_THREADS.getValue();
            try {
                numThreads = Integer.parseInt(variables.get(this).get(3));
            } catch (NumberFormatException e) {}
            if (numThreads != proj.NUM_THREADS.getValue()) {
                proj.NUM_THREADS.setValue(numThreads);
            }
            String output = variables.get(this).get(4); // gets PROJ_DIR prepended, so NOT ABSOLUTE
            (new File(ext.parseDirectoryOfFile(proj.PROJECT_DIRECTORY.getValue() + output))).mkdirs();
			CNVCaller.callAutosomalCNVs(proj, output, proj.getSamples(), null, null, numThreads, 1);// TODO, sex specific centroids,etc
            proj.CNV_FILENAMES.addValue(proj.PROJECT_DIRECTORY.getValue() + output);
            proj.saveProperties(new Project.Property[]{proj.CNV_FILENAMES});
        }
        @Override
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, Boolean> stepSelections, HashMap<STEP, ArrayList<String>> variables) {
            boolean checkHMM = Files.exists(variables.get(this).get(0)); 
            boolean checkPFB1 = stepSelections.get(S11_COMPUTE_PFB).booleanValue();
            boolean checkPFB2 = Files.exists(variables.get(this).get(1));
            boolean checkGC1 = stepSelections.get(S12_COMPUTE_GCMODEL).booleanValue();
            boolean checkGC2 = Files.exists(variables.get(this).get(2));
            int numThreads = -1;
            try {
                numThreads = Integer.parseInt(variables.get(this).get(3));
            } catch (NumberFormatException e) {}
            return new boolean[][]{
                    {checkHMM},
                    {checkPFB1, checkPFB2},
                    {checkGC1, checkGC2},
                    {numThreads > 0},
                    {!Files.exists(variables.get(this).get(4))}, 
            };
        }
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            return new Object[]{proj.HMM_FILENAME.getValue(), proj.CUSTOM_PFB_FILENAME.getValue(), proj.GC_MODEL_FILENAME.getValue(), proj.NUM_THREADS.getValue(), "cnvs/genvisis.cnv"};
        }
        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String output = variables.get(this).get(4);
            return Files.exists(proj.PROJECT_DIRECTORY.getValue() + output);
        }
        @Override
        public String getCommandLine(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String kvCmd = "";
            
            String hmm_P = proj.HMM_FILENAME.getValue();
            String hmm_G = variables.get(this).get(0);
            if (!hmm_P.equals(hmm_G)) {
                kvCmd += " HMM_FILENAME=" + hmm_G;
            }
            String pfb_P = proj.CUSTOM_PFB_FILENAME.getValue();
            String pfb_G = variables.get(this).get(1);
            if (!pfb_P.equals(pfb_G)) {
                kvCmd += " CUSTOM_PFB_FILENAME=" + pfb_G;
            }
            String gcm_P = proj.GC_MODEL_FILENAME.getValue();
            String gcm_G = variables.get(this).get(2);
            if (!gcm_P.equals(gcm_G)) {
                kvCmd += " GC_MODEL_FILENAME=" + gcm_G;
            }
            int numThreads = proj.NUM_THREADS.getValue();
            try {
                numThreads = Integer.parseInt(variables.get(this).get(3));
            } catch (NumberFormatException e) {}
            if (numThreads != proj.NUM_THREADS.getValue()) {
                kvCmd += " NUM_THREADS=" + numThreads;
            }
            String projPropFile = proj.getPropertyFilename();
            StringBuilder cmd = new StringBuilder();
            if (kvCmd.length() > 0) {
                cmd.append("jcp cnv.filesys.Project proj=" + projPropFile).append(kvCmd).append("\n");
            }
            return "jcp cnv.hmm.CNVCaller proj=" + projPropFile + " out=" + variables.get(this).get(4) + " numthreads=" + numThreads;
        }
        
    };
    
    static final STEP S99_SHADOW_SAMPLES = new STEP("Create 'Shadow' Project", 
                "", 
                new String[][]{},
                new RequirementInputType[][]{}) {
        @Override
        public void setNecessaryPreRunProperties(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            // not needed for step
        }

        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            // TODO Auto-generated method stub
        }

        @Override
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, Boolean> stepSelections, HashMap<STEP, ArrayList<String>> variables) {
            // TODO Auto-generated method stub
            return new boolean[][] {};
        }

        @Override
        public Object[] getRequirementDefaults(Project proj) {
            // TODO Auto-generated method stub
            return null;
        }

        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            // TODO Auto-generated method stub
            return false;
        }

        @Override
        public String getCommandLine(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            // TODO Auto-generated method stub
            return null;
        }
    };
    static final STEP S00_TEMPLATE = new STEP("", 
                            "", 
                            new String[][]{
                                
                            },
                            new RequirementInputType[][]{
                                
                            }) {
        @Override
        public void setNecessaryPreRunProperties(Project proj, HashMap<STEP, ArrayList<String>> variables) {
        }
        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<String>> variables) {
        }
        @Override
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, Boolean> stepSelections, HashMap<STEP, ArrayList<String>> variables) {
            return new boolean[][]{};
        }
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            return null;
        }
        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            return false;
        }
        @Override
        public String getCommandLine(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            return null;
        }
        
    };

    public static abstract class STEP {
        public String stepName;
        public String stepDesc;
        public String[][] reqs;
        private boolean failed = false;
        protected ArrayList<String> failReasons = new ArrayList<String>();
        public RequirementInputType[][] reqTypes;
        public boolean getFailed() { return failed; }
        protected void setFailed() { failed = true; }
        public ArrayList<String> getFailureMessages() { return failReasons; }
        public abstract void setNecessaryPreRunProperties(Project proj, HashMap<STEP, ArrayList<String>> variables);
        public abstract void run(Project proj, HashMap<STEP, ArrayList<String>> variables);
//        public abstract void gracefulDeath();
        public void gracefulDeath(Project proj) {
            return;
        }
        public abstract boolean[][] checkRequirements(Project proj, HashMap<STEP, Boolean> stepSelections, HashMap<STEP, ArrayList<String>> variables);
        public boolean hasRequirements(Project proj, HashMap<STEP, Boolean> stepSelections, HashMap<STEP, ArrayList<String>> variables) {
            if (stepSelections.get(this) == null || variables.get(this) == null) {
                return false;
            }
            int sum = 0;
            boolean[][] reqs = checkRequirements(proj, stepSelections, variables);
            for (boolean[] req : reqs) {
                sum += Array.booleanArraySum(req) > 0 ? 1 : 0;
            }
            return sum == reqs.length;
        }
        public String[][] getRequirements() { return reqs; }
        public RequirementInputType[][] getRequirementInputTypes() { return this.reqTypes; }
        public abstract boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variables);
        public void resetRun() {
            failed = false;
            failReasons.clear();
        }
        /**
         * Get the default values for requirements in the order they're set [skipping NONE input types]
         * @param proj
         * @return
         */
        public abstract Object[] getRequirementDefaults(Project proj);
        public abstract String getCommandLine(Project proj, HashMap<STEP, ArrayList<String>> variables);
        
        STEP(String name, String desc, String[][] requirements, RequirementInputType[][] reqTypes) {
            this.stepName = name;
            this.stepDesc = desc;
            this.reqs = requirements;
            this.reqTypes = reqTypes;
        }

        
    }
    
    public enum RequirementInputType {
        NONE(),
        FILE(),
        DIR(),
        STRING(),
        INT(),
        BOOL()
    }
    
    public GenvisisWorkflow(Project project, Launch launch) {
        this.proj = project;
        this.log = project == null ? new Logger() : project.getLog();
        this.launch = launch;
    }
    
    public void showDialogAndRun() {
        gui = new GenvisisWorkflowGUI(this.proj, this.launch);
        if (!gui.getCancelled()) {
            gui.setModal(true);
            gui.setVisible(true);
            
            if (gui.getCancelled() == true) {
                return;
            }
        }
    }

    private static STEP[] ILLUMINA_STEPS = {
        S1I_CREATE_MKR_POS,
        S2I_PARSE_SAMPLES,
        S3_CREATE_SAMPLEDATA,
        S4_TRANSPOSE_TO_MDF,
        S5_MARKER_QC,
        S6_SEX_CHECKS,
        S7_RUN_PLINK,
        S8_GWAS_QC,
        S9_ANNOTATE_SAMPLE_DATA,
        S10_GENERATE_ABLOOKUP,
        S11_COMPUTE_PFB,
        S12_COMPUTE_GCMODEL,
        S13_CREATE_PCS,
//        S14_CREATE_MT_CN_EST,
        S15_MOSAIC_ARMS,
        S16_SEX_CENTROIDS_PFB_GCMODEL,
        S17_CNV_CALLING,
    };
    private static STEP[] AFFY_STEPS = {
        S2A_PARSE_SAMPLES,
        S3_CREATE_SAMPLEDATA,
        S4_TRANSPOSE_TO_MDF,
        S5_MARKER_QC,
        S6_SEX_CHECKS,
        S7_RUN_PLINK,
        S8_GWAS_QC,
        S9_ANNOTATE_SAMPLE_DATA,
        S10_GENERATE_ABLOOKUP,
        S11_COMPUTE_PFB,
        S12_COMPUTE_GCMODEL,
        S13_CREATE_PCS,
//        S14_CREATE_MT_CN_EST,
        S15_MOSAIC_ARMS,
        S16_SEX_CENTROIDS_PFB_GCMODEL,
        S17_CNV_CALLING,
    };
    
    public static STEP[] getStepsForProject(Project proj) {
        switch (proj.ARRAY_TYPE.getValue()) {
            case AFFY_GW6: 
            case AFFY_GW6_CN:
                return AFFY_STEPS;
            case ILLUMINA:
            default:
                return ILLUMINA_STEPS;
        }
    }
    
}
