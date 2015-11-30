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
        public void run(Project proj, HashMap<STEP, ArrayList<String>> variableFields) {
            proj.getLog().report("Generating marker positions file");
            String filename = variableFields.get(this).get(0);
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
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variableFields) {
            boolean mkrPosFile = Files.exists(proj.MARKER_POSITION_FILENAME.getValue(false, false));
            boolean mkrSetFile = Files.exists(proj.MARKERSET_FILENAME.getValue(false, false));
            return mkrPosFile && mkrSetFile;
        };
        
    };
    
    static final STEP S2I_PARSE_SAMPLES = new STEP("Parse Illumina Sample Files", 
                     "", 
                     new String[][]{{"[Create Marker Positions] step must be selected and valid.", "Parsed markerPositions file must already exist."}, {"Number of Threads to Use"}}, 
                     new RequirementInputType[][]{{RequirementInputType.NONE, RequirementInputType.FILE}, {RequirementInputType.INT}}) {
        
        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<String>> variableFields) {
            proj.getLog().report("Parsing sample files");
            String projFile = proj.MARKER_POSITION_FILENAME.getValue(false, false);
            String mkrFile = variableFields.get(this).get(0);
            mkrFile = ext.verifyDirFormat(mkrFile);
            mkrFile = mkrFile.substring(0, mkrFile.length() - 1);
            if (!mkrFile.equals(projFile)) {
                proj.MARKER_POSITION_FILENAME.setValue(mkrFile);
            }
            int numThreads = proj.NUM_THREADS.getValue();
        	try {
        		numThreads = Integer.parseInt(variableFields.get(this).get(1));
        	} catch (NumberFormatException e) {}
            if (numThreads != proj.NUM_THREADS.getValue()) {
            	proj.NUM_THREADS.setValue(numThreads);
            }
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
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, Boolean> stepSelections, HashMap<STEP, ArrayList<String>> variableFields) {
        	int numThreads = -1;
        	try {
        		numThreads = Integer.parseInt(variableFields.get(this).get(1));
        	} catch (NumberFormatException e) {}
            return new boolean[][]{
                    { stepSelections.get(S1I_CREATE_MKR_POS) && S1I_CREATE_MKR_POS.hasRequirements(proj, stepSelections, variableFields),
                    	Files.exists(ext.verifyDirFormat(variableFields.get(this).get(0))),},
                	{numThreads != -1 && numThreads > 0}
            };
        }
        
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            return new Object[]{proj.MARKER_POSITION_FILENAME.getValue(false, false), proj.NUM_THREADS.getValue()};
        }

        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variableFields) {
            String sampleDirectory = proj.SAMPLE_DIRECTORY.getValue(false, false);
            // TODO strict check for #files == #samples?
            return Files.exists(sampleDirectory) && Files.list(sampleDirectory, Sample.SAMPLE_DATA_FILE_EXTENSION, false).length > 0 && proj.getSampleList() != null && proj.getSampleList().getSamples().length > 0;
        }
        
    };
    
    static final STEP S2A_PARSE_SAMPLES = new STEP("Parse Sample Files", 
            "", 
            new String[][]{{"markerPositions file must already exist."}, {"Number of Threads to Use"}}, 
            new RequirementInputType[][]{{RequirementInputType.FILE}, {RequirementInputType.INT}}) {
        
        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<String>> variableFields) {
            proj.getLog().report("Parsing sample files");
            String projFile = proj.MARKER_POSITION_FILENAME.getValue(false, false);
            String mkrFile = variableFields.get(this).get(0);
            mkrFile = ext.verifyDirFormat(mkrFile);
            mkrFile = mkrFile.substring(0, mkrFile.length() - 1);
            if (!mkrFile.equals(projFile)) {
                proj.MARKER_POSITION_FILENAME.setValue(mkrFile);
            }
            int numThreads = proj.NUM_THREADS.getValue();
            try {
                numThreads = Integer.parseInt(variableFields.get(this).get(1));
            } catch (NumberFormatException e) {}
            if (numThreads != proj.NUM_THREADS.getValue()) {
                proj.NUM_THREADS.setValue(numThreads);
            }
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
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, Boolean> stepSelections, HashMap<STEP, ArrayList<String>> variableFields) {
            int numThreads = -1;
            try {
                numThreads = Integer.parseInt(variableFields.get(this).get(1));
            } catch (NumberFormatException e) {}
            return new boolean[][]{
                    {Files.exists(ext.verifyDirFormat(variableFields.get(this).get(0))),},
                    {numThreads != -1 && numThreads > 0}
            };
        }
        
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            return new Object[]{proj.MARKER_POSITION_FILENAME.getValue(false, false), proj.NUM_THREADS.getValue()};
        }
        
        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variableFields) {
            String sampleDirectory = proj.SAMPLE_DIRECTORY.getValue(false, false);
            // TODO strict check for #files == #samples?
            return Files.exists(sampleDirectory) && Files.list(sampleDirectory, Sample.SAMPLE_DATA_FILE_EXTENSION, false).length > 0 && proj.getSampleList() != null && proj.getSampleList().getSamples().length > 0;
        }
        
    };
    
    static final STEP S3_CREATE_SAMPLEDATA = new STEP("Create SampleData.txt File", 
                         "", 
                         new String[][]{{"[Parse Sample Files] step must be selected and valid (will create a minimal SampleData.txt file)", "Parsed sample files must already exist (will create a minimal SampleData.txt file)", "A tab-delimited .PED format file with header \"" + Array.toStr(MitoPipeline.PED_INPUT, ", ") + "\"", "A Sample_Map.csv file, with at least two columns having headers \"" + MitoPipeline.SAMPLEMAP_INPUT[1] + "\" and \"" + MitoPipeline.SAMPLEMAP_INPUT[2] + "\""}}, 
                         new RequirementInputType[][]{{RequirementInputType.NONE, RequirementInputType.DIR, RequirementInputType.FILE, RequirementInputType.FILE}}) {
        
        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<String>> variableFields) {
            String projDir = proj.SAMPLE_DIRECTORY.getValue(false, false);
            String setDir = variableFields.get(this).get(0);
            if (!ext.verifyDirFormat(setDir).equals(projDir)) {
                proj.SAMPLE_DIRECTORY.setValue(setDir);
            }
            String pedFile = variableFields.get(this).get(1);
            String sampleMapCsv = variableFields.get(this).get(2);
            proj.getLog().report("Creating SampleData.txt");
            /*int retStat = */MitoPipeline.createSampleData(pedFile, sampleMapCsv, proj);
        }
        
        @Override
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, Boolean> stepSelections, HashMap<STEP, ArrayList<String>> variableFields) {
            String sampDir = variableFields.get(this).get(0);
            STEP parseStep = stepSelections.containsKey(S2I_PARSE_SAMPLES) ? S2I_PARSE_SAMPLES : S2A_PARSE_SAMPLES;
            boolean checkStepParseSamples = stepSelections.get(parseStep) && parseStep.hasRequirements(proj, stepSelections, variableFields);
            return new boolean[][]{
                    {checkStepParseSamples,
                        (Files.exists(sampDir) && Files.list(sampDir, ".sampRAF", proj.JAR_STATUS.getValue()).length > 0), 
                            Files.exists(variableFields.get(this).get(1)), 
                            Files.exists(variableFields.get(this).get(2))}};
        }
        
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            return new Object[]{proj.SAMPLE_DIRECTORY.getValue(false, false), "", ""};
        }

        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variableFields) {
            return Files.exists(proj.SAMPLE_DATA_FILENAME.getValue(false, false));
        }
        
    };
    
    static final STEP S4_TRANSPOSE_TO_MDF = new STEP("Transpose Data into Marker-Dominant Files", 
                        "", 
                        new String[][]{{"[Parse Sample Files] step must be selected and valid.", "Parsed sample files must already exist."}}, 
                        new RequirementInputType[][]{{RequirementInputType.NONE, RequirementInputType.DIR}}) {
        
        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<String>> variableFields) {
            String projDir = proj.SAMPLE_DIRECTORY.getValue(false, false);
            String setDir = variableFields.get(this).get(0);
            if (!ext.verifyDirFormat(setDir).equals(projDir)) {
                proj.SAMPLE_DIRECTORY.setValue(setDir);
            }
            proj.getLog().report("Transposing data");
            TransposeData.transposeData(proj, 2000000000, false); // compact if no LRR was provided
        }
        
        @Override
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, Boolean> stepSelections, HashMap<STEP, ArrayList<String>> variableFields) {
            String sampDir = variableFields.get(this).get(0);
            STEP parseStep = stepSelections.containsKey(S2I_PARSE_SAMPLES) ? S2I_PARSE_SAMPLES : S2A_PARSE_SAMPLES;
            boolean checkStepParseSamples = stepSelections.get(parseStep) && parseStep.hasRequirements(proj, stepSelections, variableFields);
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
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variableFields) {
            return Files.list(proj.MARKER_DATA_DIRECTORY.getValue(false, false), MarkerData.MARKER_DATA_FILE_EXTENSION, false).length > 0;
        }
        
    };
    
    static final STEP S5_EXTRACT_LRRSD = new STEP("Extract Sample Data to Lrrsd.xln File", 
                          "", 
                          new String[][]{{"[Parse Sample Files] step must be selected and valid.", "Parsed sample files must already exist."}, {"Number of Threads to Use"}}, 
                          new RequirementInputType[][]{{RequirementInputType.NONE, RequirementInputType.DIR}, {RequirementInputType.INT}}) {
        
        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<String>> variableFields) {
            proj.getLog().report("Running LrrSd");
            String projDir = proj.SAMPLE_DIRECTORY.getValue(false, false);
            String setDir = variableFields.get(this).get(0);
            if (!ext.verifyDirFormat(setDir).equals(projDir)) {
                proj.SAMPLE_DIRECTORY.setValue(setDir);
            }
            int numThreads = proj.NUM_THREADS.getValue();
        	try {
        		numThreads = Integer.parseInt(variableFields.get(this).get(1));
        	} catch (NumberFormatException e) {}
            if (numThreads != proj.NUM_THREADS.getValue()) {
            	proj.NUM_THREADS.setValue(numThreads);
            }
            cnv.qc.LrrSd.init(proj, null, null, numThreads);
        }
        
        @Override
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, Boolean> stepSelections, HashMap<STEP, ArrayList<String>> variableFields) {
        	int numThreads = -1;
        	try {
        		numThreads = Integer.parseInt(variableFields.get(this).get(1));
        	} catch (NumberFormatException e) {}
            String sampDir = variableFields.get(this).get(0);
            STEP parseStep = stepSelections.containsKey(S2I_PARSE_SAMPLES) ? S2I_PARSE_SAMPLES : S2A_PARSE_SAMPLES;
            boolean checkStepParseSamples = stepSelections.get(parseStep) && parseStep.hasRequirements(proj, stepSelections, variableFields);
            return new boolean[][]{
                    {checkStepParseSamples,
                    	(Files.exists(sampDir) && Files.list(sampDir, ".sampRAF", proj.JAR_STATUS.getValue()).length > 0),},
                    {numThreads != -1 && numThreads > 0}
            };
        }
        
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            return new Object[]{proj.SAMPLE_DIRECTORY.getValue(false, false), proj.NUM_THREADS.getValue()};
        }

        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variableFields) {
            return Files.exists(proj.SAMPLE_QC_FILENAME.getValue(false, false));
        }
    };
    
    static final STEP S6_SEX_CHECKS = new STEP("Run Sex Checks", 
                  "", 
                  new String[][]{{"[Parse Sample Files] step must be selected and valid.", "Parsed sample files must already exist."},
                                 {"[Create SampleData.txt File] step must be selected and valid.", "SampleData.txt file must already exist."}}, 
                  new RequirementInputType[][]{{RequirementInputType.NONE, RequirementInputType.DIR}, 
                                               {RequirementInputType.NONE, RequirementInputType.FILE}}) {
        
        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<String>> variableFields) {
            proj.getLog().report("Running SexCheck");
            String projDir = proj.SAMPLE_DIRECTORY.getValue(false, false);
            String setDir = variableFields.get(this).get(0);
            if (!ext.verifyDirFormat(setDir).equals(projDir)) {
                proj.SAMPLE_DIRECTORY.setValue(setDir);
            }
            cnv.qc.SexChecks.sexCheck(proj);
        }
        
        @Override
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, Boolean> stepSelections, HashMap<STEP, ArrayList<String>> variableFields) {
            String sampDir = variableFields.get(this).get(0);
            STEP parseStep = stepSelections.containsKey(S2I_PARSE_SAMPLES) ? S2I_PARSE_SAMPLES : S2A_PARSE_SAMPLES;
            boolean checkStepParseSamples = stepSelections.get(parseStep) && parseStep.hasRequirements(proj, stepSelections, variableFields);
            return new boolean[][]{
                    {checkStepParseSamples,
                     (Files.exists(sampDir) && Files.list(sampDir, ".sampRAF", proj.JAR_STATUS.getValue()).length > 0)},
                    {stepSelections.get(S3_CREATE_SAMPLEDATA) && S3_CREATE_SAMPLEDATA.hasRequirements(proj, stepSelections, variableFields),
                     Files.exists(proj.SAMPLE_DATA_FILENAME.getValue(false, false))}, 
            };
        }
        
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            return new Object[]{proj.SAMPLE_DIRECTORY.getValue(false, false), proj.SAMPLE_DATA_FILENAME.getValue(false, false)};
        }

        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variableFields) {
            return Files.exists(proj.PROJECT_DIRECTORY.getValue()+"sexCheck.xln");
        }
        
    };
    
    static final STEP S7_RUN_PLINK = new STEP("Create/Run PLINK Files", 
                 "", 
                 new String[][]{{"[Parse Sample Files] step must be selected and valid.", "Parsed sample files must already exist."}, {"A pedigree.dat file is must exist."}}, 
                 new RequirementInputType[][]{{RequirementInputType.NONE, RequirementInputType.DIR}, {RequirementInputType.FILE}}) {
        
        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<String>> variableFields) {
            String projDir = proj.SAMPLE_DIRECTORY.getValue(false, false);
            String setDir = variableFields.get(this).get(0);
            if (!ext.verifyDirFormat(setDir).equals(projDir)) {
                proj.SAMPLE_DIRECTORY.setValue(setDir);
            }
            String projPedFile = proj.PEDIGREE_FILENAME.getValue(false, false);
            String pedFile = variableFields.get(this).get(1);
            if (!pedFile.equals(projPedFile)) {
                proj.PEDIGREE_FILENAME.setValue(pedFile);
            }
            
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
                            this.failReasons.add("Error running missingness analysis.");
                        }
                    } else {
                        setFailed();
                        this.failReasons.add("Error running frequency analysis.");
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
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, Boolean> stepSelections, HashMap<STEP, ArrayList<String>> variableFields) {
            String sampDir = variableFields.get(this).get(0);
            String pedFile = variableFields.get(this).get(1);
            STEP parseStep = stepSelections.containsKey(S2I_PARSE_SAMPLES) ? S2I_PARSE_SAMPLES : S2A_PARSE_SAMPLES;
            boolean checkStepParseSamples = stepSelections.get(parseStep) && parseStep.hasRequirements(proj, stepSelections, variableFields);
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
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variableFields) {
            String fileCheck1 = proj.PROJECT_DIRECTORY.getValue()+"gwas.map";
            String fileCheck2 = proj.PROJECT_DIRECTORY.getValue()+"plink.bed";
            String fileCheck3 = proj.PROJECT_DIRECTORY.getValue()+"genome/";
            return Files.exists(fileCheck1) && Files.exists(fileCheck2) && Files.exists(fileCheck3) /*&& Files.list(fileCheck3, ".bed", false).length > 0*/;
        }
    };
    
    static final STEP S8_GWAS_QC = new STEP("Run GWAS QC", 
               "", 
               new String[][]{{"[Create/Run PLINK Files] step must be selected and valid.", "PLINK files must already exist."}, {"Keep genome info for unrelateds only?"}},
               new RequirementInputType[][]{{RequirementInputType.NONE, RequirementInputType.DIR}, {RequirementInputType.BOOL}}) {
        
        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<String>> variableFields) {
            String dir = variableFields.get(this).get(0);
            boolean keepUnrelatedsOnly = Boolean.getBoolean(variableFields.get(this).get(1));
            Qc.fullGamut(dir, keepUnrelatedsOnly, proj.getLog());
        }
        
        @Override
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, Boolean> stepSelections, HashMap<STEP, ArrayList<String>> variableFields) {
            return new boolean[][]{
                    {(stepSelections.get(S7_RUN_PLINK) && S7_RUN_PLINK.hasRequirements(proj, stepSelections, variableFields)), 
                        Files.exists(variableFields.get(this).get(0))},
                    {true}
            };
        }
        
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            return new Object[]{"", ""};
        }

        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variableFields) {
            String projDir = variableFields.get(this).get(0);
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
        
    };

    /**
     * Near-duplicate of S9A, including the requirement of Step 1 (parse Snp_map)
     * CAUTION: any changes to S9I must be reflected in S9A!
     */
    static final STEP S9I_GENERATE_ABLOOKUP = new STEP("Generate AB Lookup File", "", 
            new String[][]{{"[Create Marker Positions] step must be selected and valid.", "A MarkerSet file must already exist."}, 
                           {"[Parse Sample Files] step must be selected and valid.", "Parsed sample files must already exist."}},
            new RequirementInputType[][]{{RequirementInputType.NONE, RequirementInputType.FILE}, {RequirementInputType.NONE, RequirementInputType.DIR}}) {
        
        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<String>> variableFields) {
            String mkrPosProj = proj.MARKERSET_FILENAME.getValue(false, false);
            String mkrPosFile = variableFields.get(this).get(0);
            String setDir = proj.SAMPLE_DIRECTORY.getValue(false, false);
            String sampDir = variableFields.get(this).get(1);
            if (!mkrPosProj.equals(mkrPosFile)) {
                proj.MARKERSET_FILENAME.setValue(mkrPosFile);
            }
            if (!ext.verifyDirFormat(setDir).equals(sampDir)) {
                proj.SAMPLE_DIRECTORY.setValue(sampDir);
            }
            
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
//            if (MitoPipeline.generateABLookup(proj, proj.getLog()) == 0) {
                setFailed();
            }
        }
        
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            return new Object[]{proj.MARKERSET_FILENAME.getValue(false, false), proj.SAMPLE_DIRECTORY.getValue(false, false)};
        }
        
        @Override
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, Boolean> stepSelections, HashMap<STEP, ArrayList<String>> variableFields) {
            String mkrPosFile = variableFields.get(this).get(0);
            String sampDir = variableFields.get(this).get(1);
            boolean checkStepParseSamples = stepSelections.get(S2I_PARSE_SAMPLES) && S2I_PARSE_SAMPLES.hasRequirements(proj, stepSelections, variableFields);
            return new boolean[][]{
                    {stepSelections.get(S1I_CREATE_MKR_POS) && S1I_CREATE_MKR_POS.hasRequirements(proj, stepSelections, variableFields), 
                        Files.exists(mkrPosFile)},
                    {checkStepParseSamples,
                        (Files.exists(sampDir) && Files.list(sampDir, ".sampRAF", proj.JAR_STATUS.getValue()).length > 0)},
            };
        }
        
        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variableFields) {
            return Files.exists(proj.AB_LOOKUP_FILENAME.getValue(false, false));
        }
    };

    /**
     * Near-duplicate of S9I, excluding the requirement of Step 1 (parse Snp_map)
     * CAUTION: any changes to S9A must be reflected in S9I!
     */
    static final STEP S9A_GENERATE_ABLOOKUP = new STEP("Generate AB Lookup File", "", 
            new String[][]{{"A MarkerSet file must already exist."}, 
                            {"[Parse Sample Files] step must be selected and valid.", "Parsed sample files must already exist."}},
            new RequirementInputType[][]{{RequirementInputType.FILE}, {RequirementInputType.NONE, RequirementInputType.DIR}}) {
        
        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<String>> variableFields) {
            String mkrPosProj = proj.MARKERSET_FILENAME.getValue(false, false);
            String mkrPosFile = variableFields.get(this).get(0);
            String setDir = proj.SAMPLE_DIRECTORY.getValue(false, false);
            String sampDir = variableFields.get(this).get(1);
            if (!mkrPosProj.equals(mkrPosFile)) {
                proj.MARKERSET_FILENAME.setValue(mkrPosFile);
            }
            if (!ext.verifyDirFormat(setDir).equals(sampDir)) {
                proj.SAMPLE_DIRECTORY.setValue(sampDir);
            }
            if (MitoPipeline.generateABLookup(proj, proj.getLog()) == 0) {
                setFailed();
            }
        }
        
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            return new Object[]{proj.MARKERSET_FILENAME.getValue(false, false), proj.SAMPLE_DIRECTORY.getValue(false, false)};
        }
        
        @Override
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, Boolean> stepSelections, HashMap<STEP, ArrayList<String>> variableFields) {
            String mkrPosFile = variableFields.get(this).get(0);
            String sampDir = variableFields.get(this).get(1);
            boolean checkStepParseSamples = stepSelections.get(S2A_PARSE_SAMPLES) && S2A_PARSE_SAMPLES.hasRequirements(proj, stepSelections, variableFields);
            return new boolean[][]{
                    {Files.exists(mkrPosFile)},
                    {checkStepParseSamples,
                        (Files.exists(sampDir) && Files.list(sampDir, ".sampRAF", proj.JAR_STATUS.getValue()).length > 0)},
            };
        }
        
        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variableFields) {
            return Files.exists(proj.AB_LOOKUP_FILENAME.getValue(false, false));
        }
    };
    
    /**
     * Near-duplicate of S10A, including the requirement of Step 1 (parse Snp_map)
     * CAUTION: any changes to S10I must be reflected in S10A!
     */
    static final STEP S10I_MARKER_QC = new STEP("Run Marker QC Metrics", "",
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
        public void run(Project proj, HashMap<STEP, ArrayList<String>> variableFields) {
            double markerCallRateFilter = Double.parseDouble(variableFields.get(this).get(0));
            String mkrPosProj = proj.MARKERSET_FILENAME.getValue(false, false);
            String mkrPosFile = variableFields.get(this).get(1);
            String setSampList = proj.SAMPLELIST_FILENAME.getValue(false, false);
            String sampList = variableFields.get(this).get(2);
            String setTgtFile = proj.TARGET_MARKERS_FILENAMES.getValue()[0];
            String tgtFile = variableFields.get(this).get(3);
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
        		numThreads = Integer.parseInt(variableFields.get(this).get(4));
        	} catch (NumberFormatException e) {}
            if (numThreads != proj.NUM_THREADS.getValue()) {
            	proj.NUM_THREADS.setValue(numThreads);
            }
            
            MitoPipeline.qcMarkers(proj, "".equals(tgtFile) ? null : tgtFile, markerCallRateFilter, numThreads);
        }
        
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            return new Object[]{0.98, proj.MARKERSET_FILENAME.getValue(), proj.SAMPLELIST_FILENAME.getValue(), ""/*proj.TARGET_MARKERS_FILENAMES.getValue()[0]*/, proj.NUM_THREADS.getValue()};
        }
        
        @Override
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, Boolean> stepSelections, HashMap<STEP, ArrayList<String>> variableFields) {
            double mkr = -1;
            try {
                mkr = Double.parseDouble(variableFields.get(this).get(0));
            } catch (NumberFormatException e) {}
            String mkrPosFile = variableFields.get(this).get(1);
            String sampList = variableFields.get(this).get(2);
            String tgtFile = variableFields.get(this).get(3);
        	int numThreads = -1;
        	try {
        		numThreads = Integer.parseInt(variableFields.get(this).get(4));
        	} catch (NumberFormatException e) {}
            boolean step11 = stepSelections.get(S1I_CREATE_MKR_POS) && S1I_CREATE_MKR_POS.hasRequirements(proj, stepSelections, variableFields);
            boolean step12 = Files.exists(mkrPosFile);
            boolean step21 = stepSelections.get(S2I_PARSE_SAMPLES) && S2I_PARSE_SAMPLES.hasRequirements(proj, stepSelections, variableFields);
            boolean step22 = Files.exists(sampList);
            boolean step3 = Files.exists(tgtFile);
            return new boolean[][]{{mkr != -1}, {step11, step12}, {step21, step22}, {true, step3}, {numThreads != -1 && numThreads > 0}};
        }
        
        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variableFields) {
            String markersForABCallRate = proj.PROJECT_DIRECTORY.getValue() + MitoPipeline.MARKERS_FOR_ABCALLRATE;
            if (!Files.exists(markersForABCallRate)) {
                return false;
            }
            return true;
        }
    };
    
    /**
     * Near-duplicate of S10I, excluding the requirement of Step 1 (parse Snp_map)
     * CAUTION: any changes to S10A must be reflected in S10I!
     */
    static final STEP S10A_MARKER_QC = new STEP("Run Marker QC Metrics", "",
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
        public void run(Project proj, HashMap<STEP, ArrayList<String>> variableFields) {
            double markerCallRateFilter = Double.parseDouble(variableFields.get(this).get(0));
            String mkrPosProj = proj.MARKERSET_FILENAME.getValue(false, false);
            String mkrPosFile = variableFields.get(this).get(1);
            String setSampList = proj.SAMPLELIST_FILENAME.getValue(false, false);
            String sampList = variableFields.get(this).get(2);
            String setTgtFile = proj.TARGET_MARKERS_FILENAMES.getValue()[0];
            String tgtFile = variableFields.get(this).get(3);
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
                numThreads = Integer.parseInt(variableFields.get(this).get(4));
            } catch (NumberFormatException e) {}
            if (numThreads != proj.NUM_THREADS.getValue()) {
                proj.NUM_THREADS.setValue(numThreads);
            }
            
            MitoPipeline.qcMarkers(proj, "".equals(tgtFile) ? null : tgtFile, markerCallRateFilter, numThreads);
            
            // TODO new step for this: requires lrrsd and marker qc
            String markersForAB = Files.exists(proj.PROJECT_DIRECTORY.getValue() + MitoPipeline.MARKERS_FOR_ABCALLRATE) ? proj.PROJECT_DIRECTORY.getValue() + MitoPipeline.MARKERS_FOR_ABCALLRATE : proj.PROJECT_DIRECTORY.getValue() + MitoPipeline.MARKERS_TO_QC_FILE;
            MitoPipeline.filterSamples(proj, "PCA_GENVISIS", markersForAB, proj.PROJECT_DIRECTORY.getValue() + MitoPipeline.MARKERS_TO_QC_FILE, numThreads, "0.95", null);
        }
        
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            return new Object[]{0.98, proj.MARKERSET_FILENAME.getValue(), proj.SAMPLELIST_FILENAME.getValue(), ""/*proj.TARGET_MARKERS_FILENAMES.getValue()[0]*/, proj.NUM_THREADS.getValue()};
        }
        
        @Override
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, Boolean> stepSelections, HashMap<STEP, ArrayList<String>> variableFields) {
            double mkr = -1;
            try {
                mkr = Double.parseDouble(variableFields.get(this).get(0));
            } catch (NumberFormatException e) {}
            String mkrPosFile = variableFields.get(this).get(1);
            String sampList = variableFields.get(this).get(2);
            String tgtFile = variableFields.get(this).get(3);
            int numThreads = -1;
            try {
                numThreads = Integer.parseInt(variableFields.get(this).get(4));
            } catch (NumberFormatException e) {}
//            boolean step11 = checkBoxes.get(S1I_CREATE_MKR_POS).isSelected() && S1I_CREATE_MKR_POS.hasRequirements(proj, checkBoxes, variableFields);
            boolean step12 = Files.exists(mkrPosFile);
            boolean step21 = stepSelections.get(S2A_PARSE_SAMPLES) && S2A_PARSE_SAMPLES.hasRequirements(proj, stepSelections, variableFields);
            boolean step22 = Files.exists(sampList);
            boolean step3 = Files.exists(tgtFile);
            return new boolean[][]{{mkr != -1}, {step12}, {step21, step22}, {true, step3}, {numThreads != -1 && numThreads > 0}};
        }
        
        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variableFields) {
            String markersForABCallRate = proj.PROJECT_DIRECTORY.getValue() + MitoPipeline.MARKERS_FOR_ABCALLRATE;
            if (!Files.exists(markersForABCallRate)) {
                return false;
            }
            return true;
        }
    };
    
    static final STEP S11_FILTER_SAMPLES = new STEP("Filter Samples by Call Rate", 
            "", 
            new String[][]{
                    
                },
            new RequirementInputType[][]{
                      
                }) {

        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<String>> variableFields) {
        }
        @Override
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, Boolean> stepSelections, HashMap<STEP, ArrayList<String>> variableFields) {
            return null;
        }
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            return null;
        }
        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variableFields) {
            return false;
        }
        
    };
    
    static final STEP S11_CREATE_PCS = new STEP("Create Principal Components File", 
                  "", 
                  new String[][]{
                            {"[Transpose Data into Marker-Dominant Files] step must be selected and valid.", "Parsed marker data files must already exist."}, 
                            {"Number of Principal Components"}, 
                            {"Should impute mean value for NaN?"}, 
                            {"Should recompute Log-R ratio?"}, 
                            {"Output directory"}},
                  new RequirementInputType[][]{
                            {RequirementInputType.NONE, RequirementInputType.DIR}, 
                            {RequirementInputType.INT}, 
                            {RequirementInputType.BOOL}, 
                            {RequirementInputType.BOOL}, 
                            {RequirementInputType.DIR}
                    }) {
        
        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<String>> variableFields) {
            int numComponents = Integer.parseInt(variableFields.get(this).get(1));
            boolean imputeMeanForNaN = Boolean.getBoolean(variableFields.get(this).get(2));
            boolean recomputeLRR_PCs = Boolean.getBoolean(variableFields.get(this).get(3));
            String outputBase = variableFields.get(this).get(4);
            
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
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, Boolean> stepSelections, HashMap<STEP, ArrayList<String>> variableFields) {
            String markerDir = variableFields.get(this).get(0);
//            String medianMarkers = variableFields.get(this).get(1);
            int numComponents = -1;
            try {
                numComponents = Integer.parseInt(variableFields.get(this).get(1));
            } catch (NumberFormatException e) {}
//                boolean imputeMeanForNaN = ((JCheckBox)variableFields.get(this).get(2)).isSelected();
//                boolean recomputeLRR_PCs = ((JCheckBox)variableFields.get(this).get(3)).isSelected();
            String outputBase = variableFields.get(this).get(4);
            boolean outputCheck = (Files.exists(outputBase) && (new File(outputBase)).isDirectory());
            File fil = new File(outputBase);
            boolean exists = fil.exists();
            boolean write = fil.canWrite();
            while(!exists) {
                fil = fil.getParentFile();
                exists = fil.exists();
                write = fil.canWrite();
            }
            return new boolean[][]{
                    {(stepSelections.get(S4_TRANSPOSE_TO_MDF) && S4_TRANSPOSE_TO_MDF.hasRequirements(proj, stepSelections, variableFields)), Files.exists(markerDir)},
//                    {Files.exists(medianMarkers)},
                    {numComponents != -1},
//                    {true}, // TRUE or FALSE are both valid selections
//                    {true}, 
                    {true}, 
                    {true}, 
                    {outputCheck || (exists && write)}
            };
        }
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            return new String[]{proj.MARKER_DATA_DIRECTORY.getValue(false, false),proj.INTENSITY_PC_NUM_COMPONENTS.getValue().toString(),"", "",""};
        }

        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variableFields) {
            String outputBase = ext.rootOf(variableFields.get(this).get(4));
            String finalReport = outputBase + PrincipalComponentsResiduals.MT_REPORT_EXT[0];
            boolean mkrFiles = true;
            for (String file : PrincipalComponentsResiduals.MT_REPORT_MARKERS_USED) {
                if (!Files.exists(outputBase + file)) {
                    mkrFiles = false;
                    break;
                }
            }
            return Files.exists(finalReport) && mkrFiles;
        }
        
    };
    
    static final STEP S12_CREATE_MT_CN_EST = new STEP("Create Mitochondrial Copy-Number Estimates File", 
                            "", 
                            new String[][]{
                                        {"[Transpose Data into Marker-Dominant Files] step must be selected and valid.", "Parsed marker data files must already exist."}, 
                                        {"[Create Principal Components File] step must be selected and valid.", "Extrapolated PCs file must already exist."},
                                        {"Median Markers file"}, 
                                        {"Number of Principal Components"}, 
                                        {"Should recompute Log-R ratio median?"}, 
                                        {"Homozygous only?"}, 
                                        {"Output directory"},
                                    },
                            new RequirementInputType[][]{
                                        {RequirementInputType.NONE, RequirementInputType.DIR},
                                        {RequirementInputType.NONE, RequirementInputType.FILE},
                                        {RequirementInputType.FILE}, 
                                        {RequirementInputType.INT}, 
                                        {RequirementInputType.BOOL}, 
                                        {RequirementInputType.BOOL}, 
                                        {RequirementInputType.DIR},
                            }) {
            @Override
            public void run(Project proj, HashMap<STEP, ArrayList<String>> variableFields) {
                String extrapolatedPCsFile = "";
                int numComponents = Integer.parseInt(variableFields.get(this).get(2));
                String medianMarkers = variableFields.get(this).get(3);
                boolean recomputeLRR_Median = Boolean.getBoolean(variableFields.get(this).get(4));
                boolean homozygousOnly = Boolean.getBoolean(variableFields.get(this).get(5));
                String outputBase = variableFields.get(this).get(6);
                
                proj.getLog().report("\nComputing residuals after regressing out " + numComponents + " principal component" + (numComponents == 1 ? "" : "s") + "\n");
                PrincipalComponentsResiduals pcResids = PCA.computeResiduals(proj, extrapolatedPCsFile, ext.removeDirectoryInfo(medianMarkers), numComponents, true, 0f, homozygousOnly, recomputeLRR_Median, outputBase,null);
                MitoPipeline.generateFinalReport(proj, outputBase, pcResids.getResidOutput());
            }
            @Override
            public boolean[][] checkRequirements(Project proj, HashMap<STEP, Boolean> stepSelections, HashMap<STEP, ArrayList<String>> variableFields) {
                String markerDir = variableFields.get(this).get(0);
                String extrapPCFile = variableFields.get(this).get(1);
                String medianMarkers = variableFields.get(this).get(2);
                int numComponents = -1;
                try {
                    numComponents = Integer.parseInt(variableFields.get(this).get(3));
                } catch (NumberFormatException e) {}
    //                boolean recomputeLRR_Median = ((JCheckBox)variableFields.get(this).get(4)).isSelected();
    //                boolean homozygousOnly = ((JCheckBox)variableFields.get(this).get(5)).isSelected();
                String outputBase = variableFields.get(this).get(6);
                boolean outputCheck = (Files.exists(outputBase) && (new File(outputBase)).isDirectory());
                File fil = new File(outputBase);
                boolean exists = fil.exists();
                boolean write = fil.canWrite();
                while(!exists) {
                    fil = fil.getParentFile();
                    exists = fil.exists();
                    write = fil.canWrite();
                }
                return new boolean[][]{
                        {(stepSelections.get(S4_TRANSPOSE_TO_MDF) && S4_TRANSPOSE_TO_MDF.hasRequirements(proj, stepSelections, variableFields)), Files.exists(markerDir)},
                        {(stepSelections.get(S11_CREATE_PCS) && S11_CREATE_PCS.hasRequirements(proj, stepSelections, variableFields)), Files.exists(extrapPCFile)},
                        {Files.exists(medianMarkers)},
                        {numComponents != -1},
                        {true}, // TRUE or FALSE are both valid selections
                        {true}, 
                        {outputCheck || (exists && write)}
                };
            }
            @Override
            public Object[] getRequirementDefaults(Project proj) {
                return new String[]{proj.MARKER_DATA_DIRECTORY.getValue(false, false), proj.INTENSITY_PC_FILENAME.getValue(), "", proj.INTENSITY_PC_NUM_COMPONENTS.getValue().toString(), "", "", ""};
            }
            @Override
            public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variableFields) {
                // existing files will be backed up if re-run
                return false;
            }
        };
    static final STEP S13_COMPUTE_PFB = new STEP("Compute Population BAF files", "", new String[][]{
            {"[Parse Sample Files] step must be selected and valid (will create a SampleList file)", "A SampleList file must already exist.", "A Sample subset file must exist."}, 
            {"PFB (population BAF) output file must be specified."}},
            new RequirementInputType[][]{
                    {RequirementInputType.NONE, RequirementInputType.FILE, RequirementInputType.FILE},
                    {RequirementInputType.FILE}}
            ) {
        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<String>> variableFields) {
            String setSampList = proj.SAMPLELIST_FILENAME.getValue();
            String sampListFile = variableFields.get(this).get(0);
            String setSubSampFile = proj.SAMPLE_SUBSET_FILENAME.getValue();
            String subSampFile = variableFields.get(this).get(1);
            String setPFBFile = proj.CUSTOM_PFB_FILENAME.getValue();
            String pfbOutputFile = variableFields.get(this).get(2);
            
            if (!ext.verifyDirFormat(setSampList).equals(sampListFile)) {
                proj.SAMPLELIST_FILENAME.setValue(sampListFile);
            }
            if (!ext.verifyDirFormat(setSubSampFile).equals(subSampFile)) {
                proj.SAMPLE_SUBSET_FILENAME.setValue(subSampFile);
            }
            if (!ext.verifyDirFormat(setPFBFile).equals(pfbOutputFile)) {
                proj.CUSTOM_PFB_FILENAME.setValue(pfbOutputFile);
            }
            cnv.analysis.PennCNV.populationBAF(proj);
        }
        
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            return new Object[]{proj.SAMPLELIST_FILENAME.getValue(), proj.SAMPLE_SUBSET_FILENAME.getValue(), Files.exists(proj.SAMPLE_SUBSET_FILENAME.getValue()) ? ext.rootOf(proj.SAMPLE_SUBSET_FILENAME.getValue()) + ".pfb" : proj.CUSTOM_PFB_FILENAME.getValue()};
        }
        
        @Override
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, Boolean> stepSelections, HashMap<STEP, ArrayList<String>> variableFields) {
            String sampListFile = variableFields.get(this).get(0);
            String subSampFile = variableFields.get(this).get(1);
            String pfbOutputFile = variableFields.get(this).get(2);
            
            STEP parseStep = stepSelections.containsKey(S2I_PARSE_SAMPLES) ? S2I_PARSE_SAMPLES : S2A_PARSE_SAMPLES;
            boolean checkStepParseSamples = stepSelections.get(parseStep) && parseStep.hasRequirements(proj, stepSelections, variableFields);
            boolean step12 = Files.exists(sampListFile);
            boolean step13 = Files.exists(subSampFile);
            boolean step21 = !Files.exists(pfbOutputFile);
            
            return new boolean[][]{{checkStepParseSamples, step12, step13}, {step21}};
        }
        
        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variableFields) {
            String subSampFile = variableFields.get(this).get(1);
            String pfbOutputFile = variableFields.get(this).get(2);
            boolean pfbExists = Files.exists(pfbOutputFile) || Files.exists(ext.rootOf(subSampFile) + ".pfb");
            return pfbExists;
        }
    };
    
    static final STEP S14_COMPUTE_BAF_GCMODEL = new STEP("Compute GCMODEL File", 
                    "", 
                    new String[][]{
                                {"A GC Base file must exist."}, 
                                {"GCModel output file must be specified."}},
                    new RequirementInputType[][]{
                            {RequirementInputType.FILE}, 
                            {RequirementInputType.FILE}
                }) {
        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<String>> variableFields) {
            String gcBaseFile = variableFields.get(this).get(0);
            String setGCOutputFile = proj.GC_MODEL_FILENAME.getValue();
            String gcOutputFile = variableFields.get(this).get(1);
            if (!ext.verifyDirFormat(setGCOutputFile).equals(gcOutputFile)) {
                proj.GC_MODEL_FILENAME.setValue(gcOutputFile);
            }
            cnv.analysis.PennCNV.gcModel(proj, gcBaseFile, gcOutputFile, 100);
        }
        @Override
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, Boolean> stepSelections, HashMap<STEP, ArrayList<String>> variableFields) {
            String gcBaseFile = variableFields.get(this).get(0);
            String gcOutputFile = variableFields.get(this).get(1);
            return new boolean[][]{{Files.exists(gcBaseFile)},{!Files.exists(gcOutputFile)}};
        }
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            return new Object[]{Files.firstPathToFileThatExists(Aliases.REFERENCE_FOLDERS, "gc5Base.txt", true, false, proj.getLog()), proj.GC_MODEL_FILENAME.getValue()};
        }
        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variableFields) {
            String gcOutputFile = variableFields.get(this).get(1);
            boolean gcExists = Files.exists(gcOutputFile);
            return gcExists;
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
        public void run(Project proj, HashMap<STEP, ArrayList<String>> variableFields) {
            String mkrPosProj = proj.MARKERSET_FILENAME.getValue(false, false);
            String mkrPosFile = variableFields.get(this).get(0);
            if (!mkrPosProj.equals(mkrPosFile)) {
                proj.MARKERSET_FILENAME.setValue(mkrPosFile);
            }
            int numThreads = proj.NUM_THREADS.getValue();
            try {
                numThreads = Integer.parseInt(variableFields.get(this).get(1));
            } catch (NumberFormatException e) {}
            if (numThreads != proj.NUM_THREADS.getValue()) {
                proj.NUM_THREADS.setValue(numThreads);
            }
            Mosaicism.findOutliers(proj);
        }
        
        @Override
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, Boolean> stepSelections, HashMap<STEP, ArrayList<String>> variableFields) {
            String mkrPosFile = variableFields.get(this).get(0);
            boolean step11 = Files.exists(mkrPosFile);
            int numThreads = -1;
            try {
                numThreads = Integer.parseInt(variableFields.get(this).get(1));
            } catch (NumberFormatException e) {}
            return new boolean[][]{{step11}, {numThreads != -1 && numThreads > 0}};
        }
        
        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variableFields) {
            boolean outputCheck = Files.exists(proj.RESULTS_DIRECTORY.getValue(false, false) + "Mosaicism.xln");
            return outputCheck;
        }
        
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            return new Object[]{proj.MARKERSET_FILENAME.getValue(), proj.NUM_THREADS.getValue()};
        }
        
    };
    
    static final STEP S16_SHADOW_SAMPLES = new STEP("Create 'Shadow' Sample Files", 
                       "", 
                       new String[][]{},
                       new RequirementInputType[][]{}) {
        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<String>> variableFields) {
            // TODO Auto-generated method stub
        }
        @Override
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, Boolean> stepSelections, HashMap<STEP, ArrayList<String>> variableFields) {
            // TODO Auto-generated method stub
            return new boolean[][]{};
        }
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            // TODO Auto-generated method stub
            return null;
        }
        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variableFields) {
            // TODO Auto-generated method stub
            return false;
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
        public abstract boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<String>> variableFields);
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
        S5_EXTRACT_LRRSD,
        S6_SEX_CHECKS,
        S7_RUN_PLINK,
        S8_GWAS_QC,
        S9I_GENERATE_ABLOOKUP,
        S10I_MARKER_QC,
        S11_CREATE_PCS,
        S12_CREATE_MT_CN_EST,
        S13_COMPUTE_PFB,
        S14_COMPUTE_BAF_GCMODEL,
        S15_MOSAIC_ARMS,
        S16_SHADOW_SAMPLES
    };
    private static STEP[] AFFY_STEPS = {
        S2A_PARSE_SAMPLES,
        S3_CREATE_SAMPLEDATA,
        S4_TRANSPOSE_TO_MDF,
        S5_EXTRACT_LRRSD,
        S6_SEX_CHECKS,
        S7_RUN_PLINK,
        S8_GWAS_QC,
        S9A_GENERATE_ABLOOKUP,
        S10A_MARKER_QC,
        S11_CREATE_PCS,
        S12_CREATE_MT_CN_EST,
        S13_COMPUTE_PFB,
        S14_COMPUTE_BAF_GCMODEL,
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
