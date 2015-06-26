package cnv.manage;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;

import javax.swing.JCheckBox;
import javax.swing.JComponent;
import javax.swing.JTextField;

import cnv.Launch;
import cnv.analysis.pca.PCA;
import cnv.analysis.pca.PrincipalComponentsApply;
import cnv.analysis.pca.PrincipalComponentsCompute;
import cnv.analysis.pca.PrincipalComponentsResiduals;
import cnv.filesys.MarkerData;
import cnv.filesys.Project;
import cnv.filesys.Sample;
import cnv.gui.KitAndKaboodleGUI;
import common.Array;
import common.CmdLine;
import common.Files;
import common.Logger;
import common.ext;

public class KitAndKaboodle {
    
    Project proj;
    Logger log;
    private KitAndKaboodleGUI gui;
    private Launch launch;

    static final STEP S1_CREATE_MKR_POS = new STEP("Create Illumina Marker Positions (if not already exists)", 
                      "", 
                      new String[][]{{"An Illumina SNP_map file."}}, 
                      new RequirementInputType[][]{{RequirementInputType.FILE}}) {
        
        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<? extends JComponent>> variableFields) {
            proj.getLog().report("Generating marker positions file");
            String filename = ((JTextField)variableFields.get(this).get(0)).getText().trim();
            cnv.manage.Markers.generateMarkerPositions(proj, filename);
        }
        
        @Override
        public boolean[][] checkRequirements(Project proj, HashMap<STEP,JCheckBox> checkBoxes, HashMap<STEP, ArrayList<? extends JComponent>> variableFields) {
            return new boolean[][]{{Files.exists(((JTextField)variableFields.get(this).get(0)).getText().trim())}};
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
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<? extends JComponent>> variableFields) {
            return Files.exists(proj.MARKER_POSITION_FILENAME.getValue(false, false));
        };
        
    };
    
    static final STEP S2_PARSE_SAMPLES = new STEP("Parse Illumina Sample Files", 
                     "", 
                     new String[][]{{"Option 1 must be selected and valid.", "Parsed markerPositions file must already exist."}}, 
                     new RequirementInputType[][]{{RequirementInputType.NONE, RequirementInputType.FILE}}) {
        
        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<? extends JComponent>> variableFields) {
            proj.getLog().report("Parsing sample files");
            String projFile = proj.MARKER_POSITION_FILENAME.getValue(false, false);
            String mkrFile = ((JTextField)variableFields.get(this).get(0)).getText().trim();
            mkrFile = ext.verifyDirFormat(mkrFile);
            mkrFile = mkrFile.substring(0, mkrFile.length() - 1);
            if (!mkrFile.equals(projFile)) {
                proj.MARKER_POSITION_FILENAME.setValue(mkrFile);
            }
            int retCode = cnv.manage.ParseIllumina.createFiles(proj, proj.NUM_THREADS.getValue());
            switch (retCode) {
            case 0:
                this.setFailed();
                this.failReasons.add("Operation failure, please check log for more information.");
                break;
            case 1:
                break;
            case 6:
                this.failReasons.add("ABLookup ");
                break;
            }
        }
        
        @Override
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, JCheckBox> checkBoxes, HashMap<STEP, ArrayList<? extends JComponent>> variableFields) {
            return new boolean[][]{
                    { checkBoxes.get(S1_CREATE_MKR_POS).isSelected() && S1_CREATE_MKR_POS.hasRequirements(proj, checkBoxes, variableFields),
                    Files.exists(ext.verifyDirFormat(((JTextField)variableFields.get(this).get(0)).getText().trim())),}
            };
        }
        
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            return new Object[]{proj.MARKER_POSITION_FILENAME.getValue(false, false)};
        }

        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<? extends JComponent>> variableFields) {
            String sampleDirectory = proj.SAMPLE_DIRECTORY.getValue(false, false);
            // TODO strict check for #files == #samples?
            return Files.exists(sampleDirectory) && Files.list(sampleDirectory, Sample.SAMPLE_DATA_FILE_EXTENSION, false).length > 0 && proj.getSampleList() != null && proj.getSampleList().getSamples().length > 0;
        }
        
    };
    
    static final STEP S3_CREATE_SAMPLEDATA = new STEP("Create SampleData.txt File", 
                         "", 
                         new String[][]{{"Option 2 must be selected and valid (will create a minimal SampleData.txt file)", "Parsed sample files must already exist (will create a minimal SampleData.txt file)", "A tab-delimited .PED format file with header \"" + Array.toStr(MitoPipeline.PED_INPUT, ", ") + "\"", "A Sample_Map.csv file, with at least two columns having headers \"" + MitoPipeline.SAMPLEMAP_INPUT[1] + "\" and \"" + MitoPipeline.SAMPLEMAP_INPUT[2] + "\""}}, 
                         new RequirementInputType[][]{{RequirementInputType.NONE, RequirementInputType.DIR, RequirementInputType.FILE, RequirementInputType.FILE}}) {
        
        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<? extends JComponent>> variableFields) {
            String projDir = proj.SAMPLE_DIRECTORY.getValue(false, false);
            String setDir = ((JTextField)variableFields.get(this).get(0)).getText().trim();
            if (!ext.verifyDirFormat(setDir).equals(projDir)) {
                proj.SAMPLE_DIRECTORY.setValue(setDir);
            }
            String pedFile = ((JTextField)variableFields.get(this).get(1)).getText().trim();
            String sampleMapCsv = ((JTextField)variableFields.get(this).get(2)).getText().trim();
            proj.getLog().report("Creating SampleData.txt");
            /*int retStat = */MitoPipeline.createSampleData(pedFile, sampleMapCsv, proj);
        }
        
        @Override
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, JCheckBox> checkBoxes, HashMap<STEP, ArrayList<? extends JComponent>> variableFields) {
            String sampDir = ((JTextField)variableFields.get(this).get(0)).getText().trim();
            return new boolean[][]{
                    {checkBoxes.get(S2_PARSE_SAMPLES).isSelected() && S2_PARSE_SAMPLES.hasRequirements(proj, checkBoxes, variableFields),
                        (Files.exists(sampDir) && Files.list(sampDir, ".sampRAF", proj.JAR_STATUS.getValue()).length > 0), 
                            Files.exists(((JTextField)variableFields.get(this).get(1)).getText().trim()), 
                            Files.exists(((JTextField)variableFields.get(this).get(2)).getText().trim())}};
        }
        
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            return new Object[]{proj.SAMPLE_DIRECTORY.getValue(false, false), "", ""};
        }

        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<? extends JComponent>> variableFields) {
            return Files.exists(proj.SAMPLE_DATA_FILENAME.getValue(false, false));
        }
        
    };
    
    static final STEP S4_TRANSPOSE_TO_MDF = new STEP("Transpose Data into Marker-Dominant Files", 
                        "", 
                        new String[][]{{"Option 2 must be selected and valid.", "Parsed sample files must already exist."}}, 
                        new RequirementInputType[][]{{RequirementInputType.NONE, RequirementInputType.DIR}}) {
        
        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<? extends JComponent>> variableFields) {
            String projDir = proj.SAMPLE_DIRECTORY.getValue(false, false);
            String setDir = ((JTextField)variableFields.get(this).get(0)).getText().trim();
            if (!ext.verifyDirFormat(setDir).equals(projDir)) {
                proj.SAMPLE_DIRECTORY.setValue(setDir);
            }
            proj.getLog().report("Transposing data");
            TransposeData.transposeData(proj, 2000000000, false); // compact if no LRR was provided
        }
        
        @Override
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, JCheckBox> checkBoxes, HashMap<STEP, ArrayList<? extends JComponent>> variableFields) {
            String sampDir = ((JTextField)variableFields.get(this).get(0)).getText().trim();
            return new boolean[][]{
                    {checkBoxes.get(S2_PARSE_SAMPLES).isSelected() && S2_PARSE_SAMPLES.hasRequirements(proj, checkBoxes, variableFields),
                    (Files.exists(sampDir) && Files.list(sampDir, ".sampRAF", proj.JAR_STATUS.getValue()).length > 0),}
            };
        }
        
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            return new Object[]{proj.SAMPLE_DIRECTORY.getValue(false, false)};
        }

        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<? extends JComponent>> variableFields) {
            return Files.list(proj.MARKER_DATA_DIRECTORY.getValue(false, false), MarkerData.MARKER_DATA_FILE_EXTENSION, false).length > 0;
        }
        
    };
    
    static final STEP S5_EXTRACT_LRRSD = new STEP("Extract Sample Data to Lrrsd.xln File", 
                          "", 
                          new String[][]{{"Option 2 must be selected and valid.", "Parsed sample files must already exist."}}, 
                          new RequirementInputType[][]{{RequirementInputType.NONE, RequirementInputType.DIR}}) {
        
        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<? extends JComponent>> variableFields) {
            proj.getLog().report("Running LrrSd");
            String projDir = proj.SAMPLE_DIRECTORY.getValue(false, false);
            String setDir = ((JTextField)variableFields.get(this).get(0)).getText().trim();
            if (!ext.verifyDirFormat(setDir).equals(projDir)) {
                proj.SAMPLE_DIRECTORY.setValue(setDir);
            }
            cnv.qc.LrrSd.init(proj, null, null, proj.getProperty(proj.NUM_THREADS));
        }
        
        @Override
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, JCheckBox> checkBoxes, HashMap<STEP, ArrayList<? extends JComponent>> variableFields) {
            String sampDir = ((JTextField)variableFields.get(this).get(0)).getText().trim();
            return new boolean[][]{
                    {checkBoxes.get(S2_PARSE_SAMPLES).isSelected() && S2_PARSE_SAMPLES.hasRequirements(proj, checkBoxes, variableFields),
                    (Files.exists(sampDir) && Files.list(sampDir, ".sampRAF", proj.JAR_STATUS.getValue()).length > 0),}
            };
        }
        
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            return new Object[]{proj.SAMPLE_DIRECTORY.getValue(false, false)};
        }

        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<? extends JComponent>> variableFields) {
            return Files.exists(proj.SAMPLE_QC_FILENAME.getValue(false, false));
        }
    };
    
    static final STEP S6_SEX_CHECKS = new STEP("Run Sex Checks", 
                  "", 
                  new String[][]{{"Option 2 must be selected and valid.", "Parsed sample files must already exist."},
                                 {"Option 3 must be selected and valid.", "SampleData.txt file must already exist."}}, 
                  new RequirementInputType[][]{{RequirementInputType.NONE, RequirementInputType.DIR}, 
                                               {RequirementInputType.NONE, RequirementInputType.FILE}}) {
        
        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<? extends JComponent>> variableFields) {
            proj.getLog().report("Running SexCheck");
            String projDir = proj.SAMPLE_DIRECTORY.getValue(false, false);
            String setDir = ((JTextField)variableFields.get(this).get(0)).getText().trim();
            if (!ext.verifyDirFormat(setDir).equals(projDir)) {
                proj.SAMPLE_DIRECTORY.setValue(setDir);
            }
            cnv.qc.SexChecks.sexCheck(proj);
        }
        
        @Override
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, JCheckBox> checkBoxes, HashMap<STEP, ArrayList<? extends JComponent>> variableFields) {
            String sampDir = ((JTextField)variableFields.get(this).get(0)).getText().trim();
            return new boolean[][]{
                    {checkBoxes.get(S2_PARSE_SAMPLES).isSelected() && S2_PARSE_SAMPLES.hasRequirements(proj, checkBoxes, variableFields),
                     (Files.exists(sampDir) && Files.list(sampDir, ".sampRAF", proj.JAR_STATUS.getValue()).length > 0)},
                    {checkBoxes.get(S3_CREATE_SAMPLEDATA).isSelected() && S3_CREATE_SAMPLEDATA.hasRequirements(proj, checkBoxes, variableFields),
                     Files.exists(proj.SAMPLE_DATA_FILENAME.getValue(false, false))}, 
            };
        }
        
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            return new Object[]{proj.SAMPLE_DIRECTORY.getValue(false, false), proj.SAMPLE_DATA_FILENAME.getValue(false, false)};
        }

        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<? extends JComponent>> variableFields) {
            return Files.exists(proj.PROJECT_DIRECTORY.getValue()+"sexCheck.xln");
        }
        
    };
    
    static final STEP S7_RUN_PLINK = new STEP("Create/Run PLINK Files", 
                 "", 
                 new String[][]{{"Option 2 must be selected and valid.", "Parsed sample files must already exist."}}, 
                 new RequirementInputType[][]{{RequirementInputType.NONE, RequirementInputType.DIR}}) {
        
        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<? extends JComponent>> variableFields) {
            String projDir = proj.SAMPLE_DIRECTORY.getValue(false, false);
            String setDir = ((JTextField)variableFields.get(this).get(0)).getText().trim();
            if (!ext.verifyDirFormat(setDir).equals(projDir)) {
                proj.SAMPLE_DIRECTORY.setValue(setDir);
            }
            proj.getLog().report("Running PLINK");
            // TODO Requires pedigree file exist
            
            boolean create = cnv.manage.PlinkFormat.createPlink(proj, "gwas", null);
            if (create) {
                CmdLine.run("plink --file gwas --make-bed --out plink", proj.PROJECT_DIRECTORY.getValue());
                create = new File(proj.PROJECT_DIRECTORY.getValue()+"genome/").mkdirs();
                if (create) {
                    CmdLine.run("plink --bfile ../plink --freq", proj.PROJECT_DIRECTORY.getValue()+"genome/");
                    CmdLine.run("plink --bfile ../plink --missing", proj.PROJECT_DIRECTORY.getValue()+"genome/");
                } else {
                    setFailed();
                    this.failReasons.add("Could not create genome/ folder in " + proj.PROJECT_DIRECTORY.getValue());
                }
            } else {
                setFailed();
                this.failReasons.add("Creation of initial PLINK files failed. Please check log for more information and try again.");
            }
        }
        
        @Override
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, JCheckBox> checkBoxes, HashMap<STEP, ArrayList<? extends JComponent>> variableFields) {
            String sampDir = ((JTextField)variableFields.get(this).get(0)).getText().trim();
            return new boolean[][]{{
                    checkBoxes.get(S2_PARSE_SAMPLES).isSelected() && S2_PARSE_SAMPLES.hasRequirements(proj, checkBoxes, variableFields),
                    (Files.exists(sampDir) && Files.list(sampDir, ".sampRAF", proj.JAR_STATUS.getValue()).length > 0)}
            };
        }
        
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            return new Object[]{proj.SAMPLE_DIRECTORY.getValue(false, false)};
        }

        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<? extends JComponent>> variableFields) {
            String fileCheck1 = proj.PROJECT_DIRECTORY.getValue()+"gwas.map";
            String fileCheck2 = proj.PROJECT_DIRECTORY.getValue()+"plink.bed";
            String fileCheck3 = proj.PROJECT_DIRECTORY.getValue()+"genome/";
            return Files.exists(fileCheck1) && Files.exists(fileCheck2) && Files.exists(fileCheck3) && Files.list(fileCheck3, ".bed", false).length > 0;
        }
    };
    
    static final STEP S8_GWAS_QC = new STEP("Run GWAS QC", 
               "", 
               new String[][]{{"Option 7 must be selected", "PLINK files must already exist"}, {"Keep genome info for unrelateds only?"}},
               new RequirementInputType[][]{{RequirementInputType.NONE, RequirementInputType.DIR}, {RequirementInputType.BOOL}}) {
        
        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<? extends JComponent>> variableFields) {
            String dir = ((JTextField)variableFields.get(this).get(0)).getText().trim();
            boolean keepUnrelatedsOnly = ((JCheckBox)variableFields.get(this).get(1)).isSelected();
            gwas.Qc.fullGamut(dir, keepUnrelatedsOnly, proj.getLog());
        }
        
        @Override
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, JCheckBox> checkBoxes, HashMap<STEP, ArrayList<? extends JComponent>> variableFields) {
            return new boolean[][]{
                    {(checkBoxes.get(S7_RUN_PLINK).isSelected() && S7_RUN_PLINK.hasRequirements(proj, checkBoxes, variableFields)), Files.exists(((JTextField)variableFields.get(this).get(0)).getText())},
                    {true}
            };
        }
        
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            return new Object[]{"", ""};
        }

        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<? extends JComponent>> variableFields) {
            String projDir = ((JTextField)variableFields.get(this).get(0)).getText().trim();
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
    
    static final STEP S9_CREATE_PCS = new STEP("Create Principal Components File", 
                  "", 
                  new String[][]{{"Option 4 must be selected.", "Parsed marker data files must already exist."}, {"Median Markers file"}, {"Number of Principal Components"}, {"Should impute mean value for NaN?"}, {"Should recompute Log-R ratio?"}, {"Should recompute Log-R ratio median?"}, {"Homozygous only?"}, {"Output directory"}},
                  new RequirementInputType[][]{{RequirementInputType.NONE, RequirementInputType.DIR}, {RequirementInputType.FILE}, {RequirementInputType.INT}, {RequirementInputType.BOOL}, {RequirementInputType.BOOL}, {RequirementInputType.BOOL}, {RequirementInputType.BOOL}, {RequirementInputType.DIR}}) {
        
        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<? extends JComponent>> variableFields) {
            String medianMarkers = ((JTextField)variableFields.get(this).get(1)).getText().trim();
            int numComponents = Integer.parseInt(((JTextField) variableFields.get(this).get(2)).getText().trim());
            boolean imputeMeanForNaN = ((JCheckBox) variableFields.get(this).get(3)).isSelected();
            boolean recomputeLRR_PCs = ((JCheckBox) variableFields.get(this).get(4)).isSelected();
            boolean recomputeLRR_Median = ((JCheckBox) variableFields.get(this).get(5)).isSelected();
            boolean homozygousOnly = ((JCheckBox)variableFields.get(this).get(6)).isSelected();
            String outputBase = ((JTextField)variableFields.get(this).get(7)).getText().trim();
            
            proj.getLog().report("\nReady to perform the principal components analysis (PCA)\n");
            PrincipalComponentsCompute pcs = PCA.computePrincipalComponents(proj, false, numComponents, false, false, true, true, imputeMeanForNaN, recomputeLRR_PCs, outputBase + MitoPipeline.PCA_SAMPLES, outputBase);
            if (pcs == null) {
                setFailed();
                this.failReasons.add("# of Principal Components is greater than either the # of samples or the # of markers.  Please lower the # of PCs and try again.");
                return;
            }
            // apply PCs to everyone, we set useFile to null and excludeSamples to false to get all samples in the current project.
            // TODO, if we ever want to apply to only a subset of the project, we can do that here.....
            proj.getLog().report("\nApplying the loadings from the principal components analysis to all samples\n");
            PrincipalComponentsApply pcApply = PCA.applyLoadings(proj, numComponents, pcs.getSingularValuesFile(), pcs.getMarkerLoadingFile(), null, false, imputeMeanForNaN, recomputeLRR_PCs, outputBase);
            // Compute Medians for (MT) markers and compute residuals from PCs for everyone
            proj.getLog().report("\nComputing residuals after regressing out " + numComponents + " principal component" + (numComponents == 1 ? "" : "s") + "\n");
            PrincipalComponentsResiduals pcResids = PCA.computeResiduals(proj, pcApply.getExtrapolatedPCsFile(), ext.removeDirectoryInfo(medianMarkers), numComponents, true, 0f, homozygousOnly, recomputeLRR_Median, outputBase);
            MitoPipeline.generateFinalReport(proj, outputBase, pcResids.getResidOutput());
            proj.setProperty(proj.INTENSITY_PC_FILENAME, pcApply.getExtrapolatedPCsFile());
            proj.setProperty(proj.INTENSITY_PC_NUM_COMPONENTS, numComponents);
            proj.saveProperties();
        }
        
        @Override
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, JCheckBox> checkBoxes, HashMap<STEP, ArrayList<? extends JComponent>> variableFields) {
            String markerDir = ((JTextField)variableFields.get(this).get(0)).getText().trim();
            String medianMarkers = ((JTextField)variableFields.get(this).get(1)).getText().trim();
            int numComponents = -1;
            try {
                numComponents = Integer.parseInt(((JTextField)variableFields.get(this).get(2)).getText().trim());
            } catch (NumberFormatException e) {}
//                boolean imputeMeanForNaN = ((JCheckBox)variableFields.get(this).get(3)).isSelected();
//                boolean recomputeLRR_PCs = ((JCheckBox)variableFields.get(this).get(4)).isSelected();
//                boolean recomputeLRR_Median = ((JCheckBox)variableFields.get(this).get(5)).isSelected();
//                boolean homozygousOnly = ((JCheckBox)variableFields.get(this).get(6)).isSelected();
            String outputBase = ((JTextField)variableFields.get(this).get(7)).getText().trim();
            return new boolean[][]{
                    {(checkBoxes.get(S4_TRANSPOSE_TO_MDF).isSelected() && S4_TRANSPOSE_TO_MDF.hasRequirements(proj, checkBoxes, variableFields)), Files.exists(markerDir)},
                    {Files.exists(medianMarkers)},
                    {numComponents != -1},
                    {true}, // TRUE or FALSE are both valid selections
                    {true}, 
                    {true}, 
                    {true}, 
                    {Files.exists(outputBase)}
            };
        }
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            return new String[]{proj.MARKER_DATA_DIRECTORY.getValue(false, false),"",proj.INTENSITY_PC_NUM_COMPONENTS.getValue().toString(),"","", "", "", ""};
        }

        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<? extends JComponent>> variableFields) {
            String outputBase = ext.rootOf(((JTextField)variableFields.get(this).get(7)).getText().trim());
            
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
    
    static final STEP S10_CREATE_MT_CN_EST = new STEP("Create Mitochondrial Copy-Number Estimates File", 
                        "", 
                        new String[][]{},
                        new RequirementInputType[][]{}) {
        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<? extends JComponent>> variableFields) {
            // TODO Auto-generated method stub
        }
        @Override
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, JCheckBox> checkBoxes, HashMap<STEP, ArrayList<? extends JComponent>> variableFields) {
            // TODO Auto-generated method stub
            return new boolean[][]{};
        }
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            // TODO Auto-generated method stub
            return null;
        }
        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<? extends JComponent>> variableFields) {
            // TODO Auto-generated method stub
            return false;
        }
    };
    
    static final STEP S11_PENNCNV = new STEP("Export PennCNV Files", 
                "", 
                new String[][]{},
                new RequirementInputType[][]{}) {
        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<? extends JComponent>> variableFields) {
            // TODO Auto-generated method stub
        }
        @Override
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, JCheckBox> checkBoxes, HashMap<STEP, ArrayList<? extends JComponent>> variableFields) {
            // TODO Auto-generated method stub
            return new boolean[][]{};
        }
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            // TODO Auto-generated method stub
            return null;
        }
        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<? extends JComponent>> variableFields) {
            // TODO Auto-generated method stub
            return false;
        }
    };
    
    static final STEP S12_SHADOW_SAMPLES = new STEP("Create 'Shadow' Sample Files", 
                       "", 
                       new String[][]{},
                       new RequirementInputType[][]{}) {
        @Override
        public void run(Project proj, HashMap<STEP, ArrayList<? extends JComponent>> variableFields) {
            // TODO Auto-generated method stub
        }
        @Override
        public boolean[][] checkRequirements(Project proj, HashMap<STEP, JCheckBox> checkBoxes, HashMap<STEP, ArrayList<? extends JComponent>> variableFields) {
            // TODO Auto-generated method stub
            return new boolean[][]{};
        }
        @Override
        public Object[] getRequirementDefaults(Project proj) {
            // TODO Auto-generated method stub
            return null;
        }
        @Override
        public boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<? extends JComponent>> variableFields) {
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
        public abstract void run(Project proj, HashMap<STEP, ArrayList<? extends JComponent>> variableFields);
        public abstract boolean[][] checkRequirements(Project proj, HashMap<STEP, JCheckBox> checkBoxes, HashMap<STEP, ArrayList<? extends JComponent>> variableFields);
        public boolean hasRequirements(Project proj, HashMap<STEP, JCheckBox> checkBoxes, HashMap<STEP, ArrayList<? extends JComponent>> variableFields) {
            int sum = 0;
            boolean[][] reqs = checkRequirements(proj, checkBoxes, variableFields);
            for (boolean[] req : reqs) {
                sum += Array.booleanArraySum(req) > 0 ? 1 : 0;
            }
            return sum == reqs.length;
        }
        public String[][] getRequirements() { return reqs; }
        public RequirementInputType[][] getRequirementInputTypes() { return this.reqTypes; }
        public abstract boolean checkIfOutputExists(Project proj, HashMap<STEP, ArrayList<? extends JComponent>> variableFields);
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
    
    public KitAndKaboodle(Project project, Launch launch) {
        this.proj = project;
        this.log = project == null ? new Logger() : project.getLog();
        this.launch = launch;
    }
    
    public void showDialogAndRun() {
        gui = new KitAndKaboodleGUI(this.proj, this.launch);
        if (!gui.getCancelled()) {
            gui.setModal(true);
            gui.setVisible(true);
            
            if (gui.getCancelled() == true) {
                return;
            }
        }
    }

    private static STEP[] ILLUMINA_STEPS = {
        S1_CREATE_MKR_POS,
        S2_PARSE_SAMPLES,
        S3_CREATE_SAMPLEDATA,
        S4_TRANSPOSE_TO_MDF,
        S5_EXTRACT_LRRSD,
        S6_SEX_CHECKS,
        S7_RUN_PLINK,
        S8_GWAS_QC,
        S9_CREATE_PCS,
        S10_CREATE_MT_CN_EST,
        S11_PENNCNV,
        S12_SHADOW_SAMPLES
    };
    private static STEP[] AFFY_STEPS = {};
    private static STEP[] AFFY_NOCN_STEPS = {};
    
    public static STEP[] getStepsForProject(Project proj) {
        switch (proj.ARRAY_TYPE.getValue()) {
            case AFFY_GW6: 
                return AFFY_NOCN_STEPS;
            case AFFY_GW6_CN:
                return AFFY_STEPS;
            case ILLUMINA:
            default:
                return ILLUMINA_STEPS;
        }
    }
    
}
