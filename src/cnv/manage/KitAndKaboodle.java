package cnv.manage;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;

import javax.swing.JCheckBox;
import javax.swing.JComponent;
import javax.swing.JTextField;

import cnv.analysis.pca.PCA;
import cnv.analysis.pca.PrincipalComponentsApply;
import cnv.analysis.pca.PrincipalComponentsCompute;
import cnv.analysis.pca.PrincipalComponentsResiduals;
import cnv.filesys.Project;
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

    public static enum STEPS {

        S1_CREATE_MKR_POS("Create Illumina Marker Positions (if not already exists)", 
                          "", 
                          new String[][]{{"An Illumina SNP_map file."}}, 
                          new RequirementInputType[][]{{RequirementInputType.FILE}}) {
            
            @Override
            public void run(Project proj, HashMap<STEPS, ArrayList<? extends JComponent>> variableFields) {
                proj.getLog().report("Generating marker positions file");
                String filename = ((JTextField)variableFields.get(this).get(0)).getText().trim();
                cnv.manage.Markers.generateMarkerPositions(proj, filename);
            }
            
            @Override
            public boolean[][] checkRequirements(Project proj, HashMap<STEPS,JCheckBox> checkBoxes, HashMap<STEPS, ArrayList<? extends JComponent>> variableFields) {
                return new boolean[][]{{Files.exists(((JTextField)variableFields.get(this).get(0)).getText().trim())}};
            }
            
            @Override
            public Object[] getRequirementDefaults(Project proj) {
                String filename = proj.getLocationOfSNP_Map(false);
                if (filename == null) {
                    filename = "";
                }
                return new Object[]{filename};
            };
            
        },
        
        S2_PARSE_SAMPLES("Parse Illumina Sample Files", 
                         "", 
                         new String[][]{{"Option 1 must be selected and valid.", "Parsed markerPositions file must already exist."}}, 
                         new RequirementInputType[][]{{RequirementInputType.NONE, RequirementInputType.FILE}}) {
            
            @Override
            public void run(Project proj, HashMap<STEPS, ArrayList<? extends JComponent>> variableFields) {
                proj.getLog().report("Parsing sample files");
                String projFile = proj.MARKER_POSITION_FILENAME.getValue(false, false);
                String mkrFile = ((JTextField)variableFields.get(this).get(0)).getText().trim();
                mkrFile = ext.verifyDirFormat(mkrFile);
                mkrFile = mkrFile.substring(0, mkrFile.length() - 1);
                if (!mkrFile.equals(projFile)) {
                    proj.MARKER_POSITION_FILENAME.setValue(mkrFile);
                }
                cnv.manage.ParseIllumina.createFiles(proj, proj.NUM_THREADS.getValue());
            }
            
            @Override
            public boolean[][] checkRequirements(Project proj, HashMap<STEPS, JCheckBox> checkBoxes, HashMap<STEPS, ArrayList<? extends JComponent>> variableFields) {
                return new boolean[][]{
                        { checkBoxes.get(STEPS.S1_CREATE_MKR_POS).isSelected() && STEPS.S1_CREATE_MKR_POS.hasRequirements(proj, checkBoxes, variableFields),
                        Files.exists(ext.verifyDirFormat(((JTextField)variableFields.get(this).get(0)).getText().trim())),}
                };
            }
            
            @Override
            public Object[] getRequirementDefaults(Project proj) {
                return new Object[]{proj.MARKER_POSITION_FILENAME.getValue(false, false)};
            }
            
        },
        
        S3_CREATE_SAMPLEDATA("Create SampleData.txt File", 
                             "", 
                             new String[][]{{"Option 2 must be selected and valid (will create a minimal SampleData.txt file)", "Parsed sample files must already exist (will create a minimal SampleData.txt file)", "A tab-delimited .PED format file with header \"" + Array.toStr(MitoPipeline.PED_INPUT, ", ") + "\"", "A Sample_Map.csv file, with at least two columns having headers \"" + MitoPipeline.SAMPLEMAP_INPUT[1] + "\" and \"" + MitoPipeline.SAMPLEMAP_INPUT[2] + "\""}}, 
                             new RequirementInputType[][]{{RequirementInputType.NONE, RequirementInputType.DIR, RequirementInputType.FILE, RequirementInputType.FILE}}) {
            
            @Override
            public void run(Project proj, HashMap<STEPS, ArrayList<? extends JComponent>> variableFields) {
                String pedFile = ((JTextField)variableFields.get(this).get(0)).getText().trim();
                String sampleMapCsv = ((JTextField)variableFields.get(this).get(1)).getText().trim();
                proj.getLog().report("Creating SampleData.txt");
                /*int retStat = */MitoPipeline.createSampleData(pedFile, sampleMapCsv, proj);
            }
            
            @Override
            public boolean[][] checkRequirements(Project proj, HashMap<STEPS, JCheckBox> checkBoxes, HashMap<STEPS, ArrayList<? extends JComponent>> variableFields) {
                String sampDir = ((JTextField)variableFields.get(this).get(0)).getText().trim();
                return new boolean[][]{
                        {checkBoxes.get(STEPS.S2_PARSE_SAMPLES).isSelected() && STEPS.S2_PARSE_SAMPLES.hasRequirements(proj, checkBoxes, variableFields),
                            (Files.exists(sampDir) && Files.list(sampDir, ".sampRAF", proj.JAR_STATUS.getValue()).length > 0), 
                                Files.exists(((JTextField)variableFields.get(this).get(1)).getText().trim()), 
                                Files.exists(((JTextField)variableFields.get(this).get(2)).getText().trim())}};
            }
            
            @Override
            public Object[] getRequirementDefaults(Project proj) {
                return new Object[]{proj.SAMPLE_DIRECTORY.getValue(false, false), "", ""};
            }
            
        },
        
        S4_TRANSPOSE_TO_MDF("Transpose Data into Marker-Dominant Files", 
                            "", 
                            new String[][]{{"Option 2 must be selected and valid.", "Parsed sample files must already exist."}}, 
                            new RequirementInputType[][]{{RequirementInputType.NONE, RequirementInputType.DIR}}) {
            
            @Override
            public void run(Project proj, HashMap<STEPS, ArrayList<? extends JComponent>> variableFields) {
                String projDir = proj.SAMPLE_DIRECTORY.getValue(false, false);
                String setDir = ((JTextField)variableFields.get(this).get(0)).getText().trim();
                if (!ext.verifyDirFormat(setDir).equals(projDir)) {
                    proj.SAMPLE_DIRECTORY.setValue(setDir);
                }
                proj.getLog().report("Transposing data");
                TransposeData.transposeData(proj, 2000000000, false); // compact if no LRR was provided
            }
            
            @Override
            public boolean[][] checkRequirements(Project proj, HashMap<STEPS, JCheckBox> checkBoxes, HashMap<STEPS, ArrayList<? extends JComponent>> variableFields) {
                String sampDir = ((JTextField)variableFields.get(this).get(0)).getText().trim();
                return new boolean[][]{
                        {checkBoxes.get(STEPS.S2_PARSE_SAMPLES).isSelected() && STEPS.S2_PARSE_SAMPLES.hasRequirements(proj, checkBoxes, variableFields),
                        (Files.exists(sampDir) && Files.list(sampDir, ".sampRAF", proj.JAR_STATUS.getValue()).length > 0),}
                };
            }
            
            @Override
            public Object[] getRequirementDefaults(Project proj) {
                return new Object[]{proj.SAMPLE_DIRECTORY.getValue(false, false)};
            }
            
        },
        
        S5_EXTRACT_LRRSD("Extract Sample Data to Lrrsd.xln File", 
                              "", 
                              new String[][]{{"Option 2 must be selected and valid.", "Parsed sample files must already exist."}}, 
                              new RequirementInputType[][]{{RequirementInputType.NONE, RequirementInputType.DIR}}) {
            
            @Override
            public void run(Project proj, HashMap<STEPS, ArrayList<? extends JComponent>> variableFields) {
                proj.getLog().report("Running LrrSd");
                String projDir = proj.SAMPLE_DIRECTORY.getValue(false, false);
                String setDir = ((JTextField)variableFields.get(this).get(0)).getText().trim();
                if (!ext.verifyDirFormat(setDir).equals(projDir)) {
                    proj.SAMPLE_DIRECTORY.setValue(setDir);
                }
                cnv.qc.LrrSd.init(proj, null, null, proj.getProperty(proj.NUM_THREADS));
            }
            
            @Override
            public boolean[][] checkRequirements(Project proj, HashMap<STEPS, JCheckBox> checkBoxes, HashMap<STEPS, ArrayList<? extends JComponent>> variableFields) {
                String sampDir = ((JTextField)variableFields.get(this).get(0)).getText().trim();
                return new boolean[][]{
                        {checkBoxes.get(STEPS.S2_PARSE_SAMPLES).isSelected() && STEPS.S2_PARSE_SAMPLES.hasRequirements(proj, checkBoxes, variableFields),
                        (Files.exists(sampDir) && Files.list(sampDir, ".sampRAF", proj.JAR_STATUS.getValue()).length > 0),}
                };
            }
            
            @Override
            public Object[] getRequirementDefaults(Project proj) {
                return new Object[]{proj.SAMPLE_DIRECTORY.getValue(false, false)};
            }
        },
        
        S6_SEX_CHECKS("Run Sex Checks", 
                      "", 
                      new String[][]{{"Option 2 must be selected and valid.", "Parsed sample files must already exist."},
                                     {"Option 3 must be selected and valid.", "SampleData.txt file must already exist."}}, 
                      new RequirementInputType[][]{{RequirementInputType.NONE, RequirementInputType.DIR}, 
                                                   {RequirementInputType.NONE, RequirementInputType.FILE}}) {
            
            @Override
            public void run(Project proj, HashMap<STEPS, ArrayList<? extends JComponent>> variableFields) {
                proj.getLog().report("Running SexCheck");
                String projDir = proj.SAMPLE_DIRECTORY.getValue(false, false);
                String setDir = ((JTextField)variableFields.get(this).get(0)).getText().trim();
                if (!ext.verifyDirFormat(setDir).equals(projDir)) {
                    proj.SAMPLE_DIRECTORY.setValue(setDir);
                }
                cnv.qc.SexChecks.sexCheck(proj);
            }
            
            @Override
            public boolean[][] checkRequirements(Project proj, HashMap<STEPS, JCheckBox> checkBoxes, HashMap<STEPS, ArrayList<? extends JComponent>> variableFields) {
                String sampDir = ((JTextField)variableFields.get(this).get(0)).getText().trim();
                return new boolean[][]{
                        {checkBoxes.get(STEPS.S2_PARSE_SAMPLES).isSelected() && STEPS.S2_PARSE_SAMPLES.hasRequirements(proj, checkBoxes, variableFields),
                         (Files.exists(sampDir) && Files.list(sampDir, ".sampRAF", proj.JAR_STATUS.getValue()).length > 0)},
                        {checkBoxes.get(STEPS.S3_CREATE_SAMPLEDATA).isSelected() && STEPS.S3_CREATE_SAMPLEDATA.hasRequirements(proj, checkBoxes, variableFields),
                         Files.exists(proj.SAMPLE_DATA_FILENAME.getValue(false, false))}, 
                };
            }
            
            @Override
            public Object[] getRequirementDefaults(Project proj) {
                return new Object[]{proj.SAMPLE_DIRECTORY.getValue(false, false), proj.SAMPLE_DATA_FILENAME.getValue(false, false)};
            }
            
        },
        
        S7_RUN_PLINK("Create/Run PLINK Files", 
                     "", 
                     new String[][]{{"Option 2 must be selected and valid.", "Parsed sample files must already exist."}}, 
                     new RequirementInputType[][]{{RequirementInputType.NONE, RequirementInputType.DIR}}) {
            
            @Override
            public void run(Project proj, HashMap<STEPS, ArrayList<? extends JComponent>> variableFields) {
                String projDir = proj.SAMPLE_DIRECTORY.getValue(false, false);
                String setDir = ((JTextField)variableFields.get(this).get(0)).getText().trim();
                if (!ext.verifyDirFormat(setDir).equals(projDir)) {
                    proj.SAMPLE_DIRECTORY.setValue(setDir);
                }
                proj.getLog().report("Running PLINK");
                cnv.manage.PlinkFormat.createPlink(proj, "gwas", null);
                CmdLine.run("plink --file gwas --make-bed --out plink", proj.PROJECT_DIRECTORY.getValue());
                new File(proj.PROJECT_DIRECTORY.getValue()+"genome/").mkdirs();
                CmdLine.run("plink --bfile ../plink --freq", proj.PROJECT_DIRECTORY.getValue()+"genome/");
                CmdLine.run("plink --bfile ../plink --missing", proj.PROJECT_DIRECTORY.getValue()+"genome/");
            }
            
            @Override
            public boolean[][] checkRequirements(Project proj, HashMap<STEPS, JCheckBox> checkBoxes, HashMap<STEPS, ArrayList<? extends JComponent>> variableFields) {
                String sampDir = ((JTextField)variableFields.get(this).get(0)).getText().trim();
                return new boolean[][]{{
                        checkBoxes.get(STEPS.S2_PARSE_SAMPLES).isSelected() && STEPS.S2_PARSE_SAMPLES.hasRequirements(proj, checkBoxes, variableFields),
                        (Files.exists(sampDir) && Files.list(sampDir, ".sampRAF", proj.JAR_STATUS.getValue()).length > 0)}
                };
            }
            
            @Override
            public Object[] getRequirementDefaults(Project proj) {
                return new Object[]{proj.SAMPLE_DIRECTORY.getValue(false, false)};
            }
        },
        
        S8_GWAS_QC("Run GWAS QC", 
                   "", 
                   new String[][]{{"Option 7 must be selected", "PLINK files must already exist"}, {"Keep genome info for unrelateds only?"}},
                   new RequirementInputType[][]{{RequirementInputType.NONE, RequirementInputType.DIR}, {RequirementInputType.BOOL}}) {
            
            @Override
            public void run(Project proj, HashMap<STEPS, ArrayList<? extends JComponent>> variableFields) {
                String dir = ((JTextField)variableFields.get(this).get(0)).getText().trim();
                boolean keepUnrelatedsOnly = ((JCheckBox)variableFields.get(this).get(1)).isSelected();
                gwas.Qc.fullGamut(dir, keepUnrelatedsOnly, proj.getLog());
            }
            
            @Override
            public boolean[][] checkRequirements(Project proj, HashMap<STEPS, JCheckBox> checkBoxes, HashMap<STEPS, ArrayList<? extends JComponent>> variableFields) {
                return new boolean[][]{
                        {(checkBoxes.get(S7_RUN_PLINK).isSelected() && STEPS.S7_RUN_PLINK.hasRequirements(proj, checkBoxes, variableFields)), Files.exists(((JTextField)variableFields.get(this).get(0)).getText())},
                        {true}
                };
            }
            
            @Override
            public Object[] getRequirementDefaults(Project proj) {
                return new Object[]{"", ""};
            }
        },
        
        S9_CREATE_PCS("Create Principal Components File", 
                      "", 
                      new String[][]{{"Option 4 must be selected.", "Parsed marker data files must already exist."}, {"Median Markers file"}, {"Number of Principal Components"}, {"Should impute mean value for NaN?"}, {"Should recompute Log-R ratio?"}, {"Should recompute Log-R ratio median?"}, {"Homozygous only?"}, {"Output directory"}},
                      new RequirementInputType[][]{{RequirementInputType.NONE, RequirementInputType.DIR}, {RequirementInputType.FILE}, {RequirementInputType.INT}, {RequirementInputType.BOOL}, {RequirementInputType.BOOL}, {RequirementInputType.BOOL}, {RequirementInputType.BOOL}, {RequirementInputType.DIR}}) {
            
            @Override
            public void run(Project proj, HashMap<STEPS, ArrayList<? extends JComponent>> variableFields) {
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
                    // TODO computePrincipalComponents returns null when #PCs is > #samples or #markers
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
            public boolean[][] checkRequirements(Project proj, HashMap<STEPS, JCheckBox> checkBoxes, HashMap<STEPS, ArrayList<? extends JComponent>> variableFields) {
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
                        {(checkBoxes.get(S4_TRANSPOSE_TO_MDF).isSelected() && STEPS.S4_TRANSPOSE_TO_MDF.hasRequirements(proj, checkBoxes, variableFields)), Files.exists(markerDir)},
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
        },
        S10_CREATE_MT_CN_EST("Create Mitochondrial Copy-Number Estimates File", 
                            "", 
                            new String[][]{},
                            new RequirementInputType[][]{}) {
            @Override
            public void run(Project proj, HashMap<STEPS, ArrayList<? extends JComponent>> variableFields) {
                // TODO Auto-generated method stub
            }
            @Override
            public boolean[][] checkRequirements(Project proj, HashMap<STEPS, JCheckBox> checkBoxes, HashMap<STEPS, ArrayList<? extends JComponent>> variableFields) {
                // TODO Auto-generated method stub
                return new boolean[][]{};
            }
            @Override
            public Object[] getRequirementDefaults(Project proj) {
                // TODO Auto-generated method stub
                return null;
            }
        },
        S11_PENNCNV("Export PennCNV Files", 
                    "", 
                    new String[][]{},
                    new RequirementInputType[][]{}) {
            @Override
            public void run(Project proj, HashMap<STEPS, ArrayList<? extends JComponent>> variableFields) {
                // TODO Auto-generated method stub
            }
            @Override
            public boolean[][] checkRequirements(Project proj, HashMap<STEPS, JCheckBox> checkBoxes, HashMap<STEPS, ArrayList<? extends JComponent>> variableFields) {
                // TODO Auto-generated method stub
                return new boolean[][]{};
            }
            @Override
            public Object[] getRequirementDefaults(Project proj) {
                // TODO Auto-generated method stub
                return null;
            }
        },
        S12_SHADOW_SAMPLES("Create 'Shadow' Sample Files", 
                           "", 
                           new String[][]{},
                           new RequirementInputType[][]{}) {
            @Override
            public void run(Project proj, HashMap<STEPS, ArrayList<? extends JComponent>> variableFields) {
                // TODO Auto-generated method stub
            }
            @Override
            public boolean[][] checkRequirements(Project proj, HashMap<STEPS, JCheckBox> checkBoxes, HashMap<STEPS, ArrayList<? extends JComponent>> variableFields) {
                // TODO Auto-generated method stub
                return new boolean[][]{};
            }
            @Override
            public Object[] getRequirementDefaults(Project proj) {
                // TODO Auto-generated method stub
                return null;
            }
            
        };
        
        public String stepName;
        public String stepDesc;
        public String[][] reqs;
        public RequirementInputType[][] reqTypes;
//        public abstract boolean checkIfProbablyAlreadyRan();
//        public  void preRunChecks() {}
        public abstract void run(Project proj, HashMap<STEPS, ArrayList<? extends JComponent>> variableFields);
//        public  void postRunChecks() {}
        public abstract boolean[][] checkRequirements(Project proj, HashMap<STEPS, JCheckBox> checkBoxes, HashMap<STEPS, ArrayList<? extends JComponent>> variableFields);
        public boolean hasRequirements(Project proj, HashMap<STEPS, JCheckBox> checkBoxes, HashMap<STEPS, ArrayList<? extends JComponent>> variableFields) {
            int sum = 0;
            boolean[][] reqs = checkRequirements(proj, checkBoxes, variableFields);
            for (boolean[] req : reqs) {
                sum += Array.booleanArraySum(req) > 0 ? 1 : 0;
            }
            return sum == reqs.length;
        }
        public String[][] getRequirements() { return reqs; }
        public RequirementInputType[][] getRequirementInputTypes() { return this.reqTypes; }
        /**
         * Get the default values for requirements in the order they're set [skipping NONE input types]
         * @param proj
         * @return
         */
        public abstract Object[] getRequirementDefaults(Project proj);
        /**
         * Placeholder/Kludge for retrieving a new project object from S0_CREATE_PROJECT
         * @return
         */
        public Project getProject() { return null; }
        
        STEPS(String name, String desc, String[][] requirements, RequirementInputType[][] reqTypes) {
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
    
    public KitAndKaboodle(Project project) {
        this.proj = project;
        this.log = project == null ? new Logger() : project.getLog();
    }
    
    public void showDialogAndRun() {
        gui = new KitAndKaboodleGUI(this.proj);
        if (!gui.getCancelled()) {
            gui.setModal(true);
            gui.setVisible(true);
            
            if (gui.getCancelled() == true) {
                return;
            }
        }
    }
    
}
