package cnv.manage;

import gwas.Qc;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;

import cnv.Launch;
import cnv.analysis.Mosaicism;
import cnv.analysis.pca.PCA;
import cnv.analysis.pca.PrincipalComponentsApply;
import cnv.analysis.pca.PrincipalComponentsCompute;
import cnv.analysis.pca.PrincipalComponentsResiduals;
import cnv.filesys.ABLookup;
import cnv.filesys.MarkerData;
import cnv.filesys.Project;
import cnv.filesys.Sample;
import cnv.gui.GenvisisPipelineGUI;
import common.Aliases;
import common.Array;
import common.CmdLine;
import common.Files;
import common.Logger;
import common.ext;

public class GenvisisPipeline {
    
    Project proj;
    Logger log;
    private GenvisisPipelineGUI gui;
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
            boolean mkrSetFile = Files.exists(proj.MARKERSET_FILENAME.getValue(false, false));
            return mkrPosFile && mkrSetFile;
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
                // TODO ab genotypes not in samples files, show warning when parsing sample file headers
                this.failReasons.add("ABLookup ");
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
                	{numThreads != -1 && numThreads > 0}
            };
        }
        
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            return new Object[]{proj.MARKER_POSITION_FILENAME.getValue(false, false), proj.NUM_THREADS.getValue()};
        }

        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String sampleDirectory = proj.SAMPLE_DIRECTORY.getValue(false, false);
            // TODO strict check for #files == #samples?
            return Files.exists(sampleDirectory) && Files.list(sampleDirectory, Sample.SAMPLE_DATA_FILE_EXTENSION, false).length > 0 && proj.getSampleList() != null && proj.getSampleList().getSamples().length > 0;
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
                // TODO ab genotypes not in samples files, show warning when parsing sample file headers
                this.failReasons.add("ABLookup ");
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
                    {numThreads != -1 && numThreads > 0}
            };
        }
        
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            return new Object[]{proj.MARKER_POSITION_FILENAME.getValue(false, false), proj.NUM_THREADS.getValue()};
        }
        
        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String sampleDirectory = proj.SAMPLE_DIRECTORY.getValue(false, false);
            // TODO strict check for #files == #samples?
            return Files.exists(sampleDirectory) && Files.list(sampleDirectory, Sample.SAMPLE_DATA_FILE_EXTENSION, false).length > 0 && proj.getSampleList() != null && proj.getSampleList().getSamples().length > 0;
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
            /*int retStat = */MitoPipeline.createSampleData(pedFile, sampleMapCsv, proj);
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
            // TODO don't forget setting project properties 
            return "## << Create SampleData >> Not Implemented For Command Line Yet ##";
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
    
    static final STEP S5_SEX_CHECKS = new STEP("Run Sex Checks", 
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
            return Files.exists(proj.PROJECT_DIRECTORY.getValue()+"sexCheck.xln");
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
    
    static final STEP S6_RUN_PLINK = new STEP("Create/Run PLINK Files", 
                 "", 
                 new String[][]{{"[Parse Sample Files] step must be selected and valid.", "Parsed sample files must already exist."}, {"A pedigree.dat file is must exist."}}, 
                 new RequirementInputType[][]{{RequirementInputType.NONE, RequirementInputType.DIR}, {RequirementInputType.FILE}}) {

        @Override
        public void setNecessaryPreRunProperties(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String projDir = proj.SAMPLE_DIRECTORY.getValue(false, false);
            String setDir = variables.get(this).get(0);
            if (!ext.verifyDirFormat(setDir).equals(projDir)) {
                proj.SAMPLE_DIRECTORY.setValue(setDir);
            }
            String projPedFile = proj.PEDIGREE_FILENAME.getValue(false, false);
            String pedFile = variables.get(this).get(1);
            if (!pedFile.equals(projPedFile)) {
                proj.PEDIGREE_FILENAME.setValue(pedFile);
            }
        }
        
        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            
            proj.getLog().report("Running PLINK");
            
            boolean create = PlinkData.saveGenvisisToPlinkBedSet(proj, "gwas", null, null, -1, true);
            if (create) {
                File genDir = new File(proj.PROJECT_DIRECTORY.getValue()+"genome/");
                create = genDir.exists() || genDir.mkdirs();
                if (create) {
                    proj.getLog().report(ext.getTime() + "]\tCalculating frequencies...");
                    boolean run1 = CmdLine.run("plink --bfile ../plink --freq", null, proj.PROJECT_DIRECTORY.getValue()+"genome/", System.out, System.err, proj.getLog(), false);
                    if (run1) {
                        proj.getLog().report(ext.getTime() + "]\tCalculating missingness...");
                        boolean run2 = CmdLine.run("plink --bfile ../plink --missing", null, proj.PROJECT_DIRECTORY.getValue()+"genome/", System.out, System.err, proj.getLog(), false);
                        if (!run2) {
                            setFailed();
                            this.failReasons.add("Error occured while running missingness analysis.");
                        }
                    } else {
                        setFailed();
                        this.failReasons.add("Error occured while running frequency analysis.");
                    }
                } else {
                    setFailed();
                    this.failReasons.add("Could not create genome/ folder in " + proj.PROJECT_DIRECTORY.getValue());
                }
            } else {
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
                    {Files.exists(pedFile)}
            };
        }
        
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            return new Object[]{proj.SAMPLE_DIRECTORY.getValue(false, false), proj.PEDIGREE_FILENAME.getValue(false, false)};
        }

        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String fileCheck1 = proj.PROJECT_DIRECTORY.getValue()+"gwas.map";
            String fileCheck2 = proj.PROJECT_DIRECTORY.getValue()+"plink.bed";
            String fileCheck3 = proj.PROJECT_DIRECTORY.getValue()+"genome/";
            return Files.exists(fileCheck1) && Files.exists(fileCheck2) && Files.exists(fileCheck3) /*&& Files.list(fileCheck3, ".bed", false).length > 0*/;
        }
        @Override
        public String getCommandLine(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            return "## << Run Plink >> Not Implemented For Command Line Yet ##";
        }
    };
    
    static final STEP S7_GWAS_QC = new STEP("Run GWAS QC", 
               "", 
               new String[][]{{"[Create/Run PLINK Files] step must be selected and valid.", "PLINK files must already exist."}, {"Keep genome info for unrelateds only?"}},
               new RequirementInputType[][]{{RequirementInputType.NONE, RequirementInputType.DIR}, {RequirementInputType.BOOL}}) {

        @Override
        public void setNecessaryPreRunProperties(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            // not needed for step
        }
        
        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String dir = variables.get(this).get(0);
            boolean keepUnrelatedsOnly = Boolean.valueOf(variables.get(this).get(1));
            Qc.fullGamut(dir, keepUnrelatedsOnly, proj.getLog());
        }
        
        @Override
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, Boolean> stepSelections, HashMap<STEP, ArrayList<String>> variables) {
            return new boolean[][]{
                    {(stepSelections.get(S6_RUN_PLINK) && S6_RUN_PLINK.hasRequirements(proj, stepSelections, variables)), 
                        Files.exists(variables.get(this).get(0))},
                    {true}
            };
        }
        
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            return new Object[]{"", ""};
        }

        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String projDir = variables.get(this).get(0);
            boolean allExist = true;
            folders: for (int i = 0; i < gwas.Qc.FOLDERS_CREATED.length; i++) {
                for (int j = 0; j < gwas.Qc.FILES_CREATED[i].length; j++) {
                    if (!Files.exists(projDir + gwas.Qc.FOLDERS_CREATED[i] + gwas.Qc.FILES_CREATED[i][j])) {
                        allExist = false;
                        break folders;
                    }
                }
            }
            return allExist;
        }
        
        @Override
        public String getCommandLine(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String dir = variables.get(this).get(0);
            boolean keepUnrelatedsOnly = Boolean.valueOf(variables.get(this).get(1));
            String logFile = proj.getLog().getFilename();
            return "jcp gwas.Qc dir=" + dir + " log=" + logFile + " keepGenomeInfoForRelatedsOnly=" + keepUnrelatedsOnly;
        }
        
    };

    /**
     * Near-duplicate of S9A, including the requirement of Step 1 (parse Snp_map)
     * CAUTION: any changes to S9I must be reflected in S9A!
     */
    static final STEP S8I_GENERATE_ABLOOKUP = new STEP("Generate AB Lookup File", "", 
            new String[][]{{"[Create Marker Positions] step must be selected and valid.", "A MarkerSet file must already exist."}, 
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
                    {stepSelections.get(S1I_CREATE_MKR_POS) && S1I_CREATE_MKR_POS.hasRequirements(proj, stepSelections, variables), 
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
            return "## << Generate ABLookup >> Not Implemented For Command Line Yet ##";
        }
    };

    /**
     * Near-duplicate of S9I, excluding the requirement of Step 1 (parse Snp_map)
     * CAUTION: any changes to S9A must be reflected in S9I!
     */
    static final STEP S8A_GENERATE_ABLOOKUP = new STEP("Generate AB Lookup File", "", 
            new String[][]{{"A MarkerSet file must already exist."}, 
                            {"[Parse Sample Files] step must be selected and valid.", "Parsed sample files must already exist."}},
            new RequirementInputType[][]{{RequirementInputType.FILE}, {RequirementInputType.NONE, RequirementInputType.DIR}}) {

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
            boolean checkStepParseSamples = stepSelections.get(S2A_PARSE_SAMPLES) && S2A_PARSE_SAMPLES.hasRequirements(proj, stepSelections, variables);
            return new boolean[][]{
                    {Files.exists(mkrPosFile)},
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
            return "## << Generate ABLookup >> Not Implemented For Command Line Yet ##";
        }
    };
    
    /**
     * Near-duplicate of S10A, including the requirement of Step 1 (parse Snp_map)
     * CAUTION: any changes to S10I must be reflected in S10A!
     */
    static final STEP S9I_MARKER_QC = new STEP("Run Marker QC Metrics", "",
            new String[][]{{"Marker Call-Rate Filter Threshold"},
    						{"[Create Marker Positions] step must be selected and valid.", "A MarkerSet file must already exist."}, 
    						{"[Parse Sample Files] step must be selected and valid (will create a SampleList file)", "A SampleList file must already exist."},
    						{"Export all markers in project", "A targetMarkers files listing the markers to QC"},
    						{"Number of Threads to Use"}},
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
            String setTgtFile = proj.TARGET_MARKERS_FILENAMES.getValue()[0];
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
            
            MitoPipeline.qcMarkers(proj, "".equals(tgtFile) ? null : tgtFile, markerCallRateFilter, numThreads);
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
            boolean step11 = stepSelections.get(S1I_CREATE_MKR_POS) && S1I_CREATE_MKR_POS.hasRequirements(proj, stepSelections, variables);
            boolean step12 = Files.exists(mkrPosFile);
            boolean step21 = stepSelections.get(S2I_PARSE_SAMPLES) && S2I_PARSE_SAMPLES.hasRequirements(proj, stepSelections, variables);
            boolean step22 = Files.exists(sampList);
            boolean step3 = Files.exists(tgtFile);
            return new boolean[][]{{mkr != -1}, {step11, step12}, {step21, step22}, {true, step3}, {numThreads != -1 && numThreads > 0}};
        }
        
        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String markersForABCallRate = proj.PROJECT_DIRECTORY.getValue() + MitoPipeline.MARKERS_FOR_ABCALLRATE;
            if (!Files.exists(markersForABCallRate)) {
                return false;
            }
            return true;
        }
        @Override
        public String getCommandLine(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            return "## << Marker QC >> Not Implemented For Command Line Yet ##";
        }
    };
    
    /**
     * Near-duplicate of S10I, excluding the requirement of Step 1 (parse Snp_map)
     * CAUTION: any changes to S10A must be reflected in S10I!
     */
    static final STEP S9A_MARKER_QC = new STEP("Run Marker QC Metrics", "",
                    new String[][]{{"Marker Call-Rate Filter Threshold"},
                                    {"A MarkerSet file must already exist."}, 
                                    {"[Parse Sample Files] step must be selected and valid (will create a SampleList file)", "A SampleList file must already exist."},
                                    {"Export all markers in project", "A targetMarkers files listing the markers to QC"},
                                    {"Number of Threads to Use"}},
                    new RequirementInputType[][]{{RequirementInputType.INT},
                                                    {RequirementInputType.FILE}, 
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
            String setTgtFile = proj.TARGET_MARKERS_FILENAMES.getValue()[0];
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
            
            MitoPipeline.qcMarkers(proj, "".equals(tgtFile) ? null : tgtFile, markerCallRateFilter, numThreads);
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
//            boolean step11 = checkBoxes.get(S1I_CREATE_MKR_POS).isSelected() && S1I_CREATE_MKR_POS.hasRequirements(proj, checkBoxes, variables);
            boolean step12 = Files.exists(mkrPosFile);
            boolean step21 = stepSelections.get(S2A_PARSE_SAMPLES) && S2A_PARSE_SAMPLES.hasRequirements(proj, stepSelections, variables);
            boolean step22 = Files.exists(sampList);
            boolean step3 = Files.exists(tgtFile);
            return new boolean[][]{{mkr != -1}, {step12}, {step21, step22}, {true, step3}, {numThreads != -1 && numThreads > 0}};
        }
        
        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String markersForABCallRate = proj.PROJECT_DIRECTORY.getValue() + MitoPipeline.MARKERS_FOR_ABCALLRATE;
            if (!Files.exists(markersForABCallRate)) {
                return false;
            }
            return true;
        }
        @Override
        public String getCommandLine(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            return "## << Marker QC >> Not Implemented For Command Line Yet ##";
        }
    };
    
    static final STEP S10_EXTRACT_LRRSD_AND_FILTER = new STEP("Extract Lrrsd.xln File and Filter Samples by CallRate", 
                          "", 
                          new String[][]{
                                {"[Parse Sample Files] step must be selected and valid.", "Parsed sample files must already exist."}, 
                                {"Number of Threads to Use"},
                                {"Sample CallRate Threshold"},
                            }, 
                          new RequirementInputType[][]{{RequirementInputType.NONE, RequirementInputType.DIR}, {RequirementInputType.INT}, {RequirementInputType.STRING}}) {

        @Override
        public void setNecessaryPreRunProperties(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String projDir = proj.SAMPLE_DIRECTORY.getValue(false, false);
            String setDir = variables.get(this).get(0);
            if (!ext.verifyDirFormat(setDir).equals(projDir)) {
                proj.SAMPLE_DIRECTORY.setValue(setDir);
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
            proj.getLog().report("Running LrrSd");
            int numThreads = proj.NUM_THREADS.getValue();
            try {
                numThreads = Integer.parseInt(variables.get(this).get(1));
            } catch (NumberFormatException e) {}
            String callRate = variables.get(this).get(2);
            String markersForAB = Files.exists(proj.PROJECT_DIRECTORY.getValue() + MitoPipeline.MARKERS_FOR_ABCALLRATE) ? proj.PROJECT_DIRECTORY.getValue() + MitoPipeline.MARKERS_FOR_ABCALLRATE : proj.PROJECT_DIRECTORY.getValue() + MitoPipeline.MARKERS_TO_QC_FILE;
            String markersForEverything =  proj.PROJECT_DIRECTORY.getValue() + MitoPipeline.MARKERS_TO_QC_FILE;
            MitoPipeline.filterSamples(proj, "PCA_GENVISIS", markersForAB, markersForEverything, numThreads, callRate, null);
        }
        
        @Override
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, Boolean> stepSelections, HashMap<STEP, ArrayList<String>> variables) {
        	int numThreads = -1;
        	try {
        		numThreads = Integer.parseInt(variables.get(this).get(1));
        	} catch (NumberFormatException e) {}
            String sampDir = variables.get(this).get(0);
            STEP parseStep = stepSelections.containsKey(S2I_PARSE_SAMPLES) ? S2I_PARSE_SAMPLES : S2A_PARSE_SAMPLES;
            boolean checkStepParseSamples = stepSelections.get(parseStep) && parseStep.hasRequirements(proj, stepSelections, variables);
            boolean validCallRate = false;
            try {
                Double.parseDouble(variables.get(this).get(2));
                validCallRate = true;
            } catch (NumberFormatException e) {}
            return new boolean[][]{
                    {checkStepParseSamples,
                    	(Files.exists(sampDir) && Files.list(sampDir, ".sampRAF", proj.JAR_STATUS.getValue()).length > 0),},
                    {numThreads != -1 && numThreads > 0},
                    {validCallRate}
            };
        }
        
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            return new Object[]{proj.SAMPLE_DIRECTORY.getValue(false, false), proj.NUM_THREADS.getValue(), "0.95"};
        }
    
        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            return Files.exists(proj.SAMPLE_QC_FILENAME.getValue(false, false)) || (Files.exists("PCA_GENVISIS" + MitoPipeline.PCA_SAMPLES) && Files.exists("PCA_GENVISIS" + MitoPipeline.PCA_SAMPLES_SUMMARY));
        }
        @Override
        public String getCommandLine(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            return "## << Extract LRRSD and Filter Samples by CallRate >> Not Implemented For Command Line Yet ##";
        }
    };
    
    static final STEP S11_COMPUTE_PFB = new STEP("Compute Population BAF files", "", new String[][]{
            {"[Parse Sample Files] step must be selected and valid (will create a SampleList file)", "A SampleList file must already exist.", "A Sample subset file must exist."}, 
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
    static final STEP S13_CREATE_PCS = new STEP("Create Principal Components File", 
                  "", 
                  new String[][]{
                            {"[Transpose Data into Marker-Dominant Files] step must be selected and valid.", "Parsed marker data files must already exist."}, 
                            {"Number of Principal Components"}, 
                            {"Should impute mean value for NaN?"}, 
                            {"Should recompute Log-R ratio?"}, 
                            },
                  new RequirementInputType[][]{
                            {RequirementInputType.NONE, RequirementInputType.DIR}, 
                            {RequirementInputType.INT}, 
                            {RequirementInputType.BOOL}, 
                            {RequirementInputType.BOOL}, 
                    }) {

        @Override
        public void setNecessaryPreRunProperties(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            // not needed for step
        }
        
        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            int numComponents = Integer.parseInt(variables.get(this).get(1));
            boolean imputeMeanForNaN = Boolean.valueOf(variables.get(this).get(2));
            boolean recomputeLRR_PCs = Boolean.valueOf(variables.get(this).get(3));
            String outputBase = proj.PROJECT_DIRECTORY.getValue() + "PCA_GENVISIS";//variables.get(this).get(4);
            
            proj.getLog().report("\nReady to perform the principal components analysis (PCA)\n");
			//TODO, load gc params as needed instead of passing null...
            PrincipalComponentsCompute pcs = PCA.computePrincipalComponents(proj, false, numComponents, false, false, true, true, imputeMeanForNaN, recomputeLRR_PCs, outputBase + MitoPipeline.PCA_SAMPLES, outputBase, null);
            if (pcs == null) {
                setFailed();
                this.failReasons.add("# of Principal Components is greater than either the # of samples or the # of markers.  Please lower the # of PCs and try again.");
                return;
            }
            // apply PCs to everyone, we set useFile to null and excludeSamples to false to get all samples in the current project.
            // TODO, if we ever want to apply to only a subset of the project, we can do that here.....
            proj.getLog().report("\nApplying the loadings from the principal components analysis to all samples\n");
			PrincipalComponentsApply pcApply = PCA.applyLoadings(proj, numComponents, pcs.getSingularValuesFile(), pcs.getMarkerLoadingFile(), null, false, imputeMeanForNaN, recomputeLRR_PCs, outputBase, null);
            // Compute Medians for (MT) markers and compute residuals from PCs for everyone
            proj.setProperty(proj.INTENSITY_PC_FILENAME, pcApply.getExtrapolatedPCsFile());
            proj.setProperty(proj.INTENSITY_PC_NUM_COMPONENTS, numComponents);
            proj.saveProperties();
        }
        
        @Override
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, Boolean> stepSelections, HashMap<STEP, ArrayList<String>> variables) {
            String markerDir = variables.get(this).get(0);
            int numComponents = -1;
            try {
                numComponents = Integer.parseInt(variables.get(this).get(1));
            } catch (NumberFormatException e) {}
            return new boolean[][]{
                    {(stepSelections.get(S4_TRANSPOSE_TO_MDF) && S4_TRANSPOSE_TO_MDF.hasRequirements(proj, stepSelections, variables)), Files.exists(markerDir)},
                    {numComponents != -1},
                    {true}, // TRUE or FALSE are both valid selections
                    {true}, 
            };
        }
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            return new String[]{proj.MARKER_DATA_DIRECTORY.getValue(false, false),proj.INTENSITY_PC_NUM_COMPONENTS.getValue().toString(), "true", "true"/*,""*/};
        }

        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String outputBase = proj.PROJECT_DIRECTORY.getValue() + "PCA_GENVISIS";//ext.rootOf(variables.get(this).get(4));
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
            return "## << Create PrincipalComponents File >> Not Implemented For Command Line Yet ##";
        }
        
    };
    
    static final STEP S14_CREATE_MT_CN_EST = new STEP("Create Mitochondrial Copy-Number Estimates File", 
                        "", 
                        new String[][]{
                                    {"[Transpose Data into Marker-Dominant Files] step must be selected and valid.", "Parsed marker data files must already exist."}, 
                                    {"[Create Principal Components File] step must be selected and valid.", "Extrapolated PCs file must already exist."},
                                    {"Median Markers file"}, 
                                    {"Number of Principal Components"}, 
                                    {"Should recompute Log-R ratio median?"}, 
                                    {"Homozygous only?"}, 
                                },
                        new RequirementInputType[][]{
                                    {RequirementInputType.NONE, RequirementInputType.DIR},
                                    {RequirementInputType.NONE, RequirementInputType.FILE},
                                    {RequirementInputType.FILE}, 
                                    {RequirementInputType.INT}, 
                                    {RequirementInputType.BOOL}, 
                                    {RequirementInputType.BOOL}, 
                        }) {

        @Override
        public void setNecessaryPreRunProperties(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            // not needed for step
        }
        
        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            String extrapolatedPCsFile = "";
            int numComponents = Integer.parseInt(variables.get(this).get(2));
            String medianMarkers = variables.get(this).get(3);
            boolean recomputeLRR_Median = Boolean.valueOf(variables.get(this).get(4));
            boolean homozygousOnly = Boolean.valueOf(variables.get(this).get(5));
            String outputBase = proj.PROJECT_DIRECTORY.getValue() + "PCA_GENVISIS";//variables.get(this).get(6);
            
            proj.getLog().report("\nComputing residuals after regressing out " + numComponents + " principal component" + (numComponents == 1 ? "" : "s") + "\n");
            PrincipalComponentsResiduals pcResids = PCA.computeResiduals(proj, extrapolatedPCsFile, ext.removeDirectoryInfo(medianMarkers), numComponents, true, 0f, homozygousOnly, recomputeLRR_Median, outputBase, null);
            MitoPipeline.generateFinalReport(proj, outputBase, pcResids.getResidOutput());
        }
        @Override
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, Boolean> stepSelections, HashMap<STEP, ArrayList<String>> variables) {
            String markerDir = variables.get(this).get(0);
            String extrapPCFile = variables.get(this).get(1);
            String medianMarkers = variables.get(this).get(2);
            int numComponents = -1;
            try {
                numComponents = Integer.parseInt(variables.get(this).get(3));
            } catch (NumberFormatException e) {}
            return new boolean[][]{
                    {(stepSelections.get(S4_TRANSPOSE_TO_MDF) && S4_TRANSPOSE_TO_MDF.hasRequirements(proj, stepSelections, variables)), Files.exists(markerDir)},
                    {(stepSelections.get(S13_CREATE_PCS) && S13_CREATE_PCS.hasRequirements(proj, stepSelections, variables)), Files.exists(extrapPCFile)},
                    {Files.exists(medianMarkers)},
                    {numComponents != -1},
                    {true}, // TRUE or FALSE are both valid selections
                    {true}, 
            };
        }
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            return new String[]{proj.MARKER_DATA_DIRECTORY.getValue(false, false), proj.INTENSITY_PC_FILENAME.getValue(), "", proj.INTENSITY_PC_NUM_COMPONENTS.getValue().toString(), "", ""/*, ""*/};
        }
        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            // existing files will be backed up if re-run
            return false;
        }
        @Override
        public String getCommandLine(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            return "## << Create Mitochondrial Copy-Number Estimates File >> Not Implemented For Command Line Yet ##";
        }
    };
    
    static final STEP S15_MOSAIC_ARMS = new STEP("Create Mosaic Arms File", 
                                                 "", 
                                                 new String[][]{
                                                        {"A MarkerSet file must already exist."}, 
                                                        {"Number of Threads to Use"}},
                                                 new RequirementInputType[][]{
                                                        {RequirementInputType.FILE}, 
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
            String mkrPosFile = variables.get(this).get(0);
            boolean step11 = Files.exists(mkrPosFile);
            int numThreads = -1;
            try {
                numThreads = Integer.parseInt(variables.get(this).get(1));
            } catch (NumberFormatException e) {}
            return new boolean[][]{{step11}, {numThreads != -1 && numThreads > 0}};
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
    
    static final STEP S16_SHADOW_SAMPLES = new STEP("Create 'Shadow' Sample Files", 
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
            return new boolean[][]{};
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
    
    static final STEP S17_SHADOW_MARKERS = new STEP("Create 'Shadow' Marker Files", 
                            "", 
                            new String[][]{
                                
                            },
                            new RequirementInputType[][]{
                                
                            }) {
        @Override
        public void setNecessaryPreRunProperties(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            // TODO Auto-generated method stub
        }
        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<String>> variables) {
            // TODO Auto-generated method stub
        }
        @Override
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, Boolean> stepSelections, HashMap<STEP, ArrayList<String>> variables) {
            // TODO Auto-generated method stub
            return new boolean[][]{};
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
        public abstract void run(Project proj, HashMap<STEP, ArrayList<String>> variables);
        public abstract void setNecessaryPreRunProperties(Project proj, HashMap<STEP, ArrayList<String>> variables);
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
    
    public GenvisisPipeline(Project project, Launch launch) {
        this.proj = project;
        this.log = project == null ? new Logger() : project.getLog();
        this.launch = launch;
    }
    
    public void showDialogAndRun() {
        gui = new GenvisisPipelineGUI(this.proj, this.launch);
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
        S5_SEX_CHECKS,
        S6_RUN_PLINK,
        S7_GWAS_QC,
        S8I_GENERATE_ABLOOKUP,
        S9I_MARKER_QC,
        S10_EXTRACT_LRRSD_AND_FILTER,
        S11_COMPUTE_PFB,
        S12_COMPUTE_GCMODEL,
        S13_CREATE_PCS,
        S14_CREATE_MT_CN_EST,
        S15_MOSAIC_ARMS,
        S16_SHADOW_SAMPLES
    };
    private static STEP[] AFFY_STEPS = {
        S2A_PARSE_SAMPLES,
        S3_CREATE_SAMPLEDATA,
        S4_TRANSPOSE_TO_MDF,
        S5_SEX_CHECKS,
        S6_RUN_PLINK,
        S7_GWAS_QC,
        S8A_GENERATE_ABLOOKUP,
        S9A_MARKER_QC,
        S10_EXTRACT_LRRSD_AND_FILTER,
        S11_COMPUTE_PFB,
        S12_COMPUTE_GCMODEL,
        S13_CREATE_PCS,
        S14_CREATE_MT_CN_EST,
        S15_MOSAIC_ARMS,
        S16_SHADOW_SAMPLES
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
