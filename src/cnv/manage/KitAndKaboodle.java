package cnv.manage;

import java.io.File;
import java.util.HashMap;

import javax.swing.JCheckBox;

import cnv.filesys.Project;
import cnv.gui.KitAndKaboodleGUI;
import common.Array;
import common.CmdLine;
import common.Files;
import common.Logger;

public class KitAndKaboodle {
    
    Project proj;
    Logger log;

public static enum STEPS {
        
        S1_CREATE_MKR_POS("Create Illumina Marker Positions (if not already exists)", "", new String[]{"An Illumina SNP_map file."}){
            @Override
            public void run(Project proj) {
                proj.getLog().report("Generating marker positions file");
                String filename = proj.getLocationOfSNP_Map(false);
                if (filename == null) {
                    return;
                }
                cnv.manage.Markers.generateMarkerPositions(proj, filename);
            }
            public boolean[] hasRequirements(Project proj, HashMap<STEPS,JCheckBox> checkBoxes) {
                return new boolean[]{proj.getLocationOfSNP_Map(false) != null};
            };
        },
        S2_PARSE_SAMPLES("Parse Illumina Sample Files", "", new String[]{"Option 1 must be selected and valid.", "Parsed markerPositions file must already exist."}) {
            @Override
            public void run(Project proj) {
                proj.getLog().report("Parsing sample files");
                cnv.manage.ParseIllumina.createFiles(proj, proj.NUM_THREADS.getValue());
            }
            @Override
            public boolean[] hasRequirements(Project proj, HashMap<STEPS, JCheckBox> checkBoxes) {
                return new boolean[]{
                        checkBoxes.get(STEPS.S1_CREATE_MKR_POS).isSelected() && Array.booleanArraySum(STEPS.S1_CREATE_MKR_POS.hasRequirements(proj, checkBoxes)) > 0,
                        Files.exists(proj.MARKER_POSITION_FILENAME.getValue(false, false)),
                };
            }
        },
        S3_TRANSPOSE_TO_MDF("Transpose Data into Marker-dominant Files", "", new String[]{"Option 2 must be selected and valid.", "Parsed sample files must already exist."}) {
            @Override
            public void run(Project proj) {
                proj.getLog().report("Transposing data");
                TransposeData.transposeData(proj, 2000000000, false); // compact if no LRR was provided
            }
            @Override
            public boolean[] hasRequirements(Project proj, HashMap<STEPS, JCheckBox> checkBoxes) {
                String sampDir = proj.SAMPLE_DIRECTORY.getValue(false, false);
                return new boolean[]{
                        checkBoxes.get(STEPS.S2_PARSE_SAMPLES).isSelected() && Array.booleanArraySum(STEPS.S2_PARSE_SAMPLES.hasRequirements(proj, checkBoxes)) > 0,
                        (Files.exists(sampDir) && Files.list(sampDir, ".sampRAF", proj.JAR_STATUS.getValue()).length > 0),
                };
            }
        },
        S4_EXTRACT_SAMPLEDATA("Extract Sample Data to lrrsd.xln File", "", new String[]{"Option 2 must be selected and valid.", "Parsed sample files must already exist."}) {
            @Override
            public void run(Project proj) {
                proj.getLog().report("Running LrrSd");
                cnv.qc.LrrSd.init(proj, null, null, proj.getProperty(proj.NUM_THREADS));
            }
            @Override
            public boolean[] hasRequirements(Project proj, HashMap<STEPS, JCheckBox> checkBoxes) {
                String sampDir = proj.SAMPLE_DIRECTORY.getValue(false, false);
                return new boolean[]{
                        checkBoxes.get(STEPS.S2_PARSE_SAMPLES).isSelected() && Array.booleanArraySum(STEPS.S2_PARSE_SAMPLES.hasRequirements(proj, checkBoxes)) > 0,
                        (Files.exists(sampDir) && Files.list(sampDir, ".sampRAF", proj.JAR_STATUS.getValue()).length > 0),
                };
            }
        },
        S5_SEX_CHECKS("Run Sex Checks", "", new String[]{"Option 2 must be selected and valid.", "Parsed sample files must already exist."}) {
            @Override
            public void run(Project proj) {
                proj.getLog().report("Running SexCheck");
                cnv.qc.SexChecks.sexCheck(proj);
            }
            @Override
            public boolean[] hasRequirements(Project proj, HashMap<STEPS, JCheckBox> checkBoxes) {
                String sampDir = proj.SAMPLE_DIRECTORY.getValue(false, false);
                return new boolean[]{
                        checkBoxes.get(STEPS.S2_PARSE_SAMPLES).isSelected() && Array.booleanArraySum(STEPS.S2_PARSE_SAMPLES.hasRequirements(proj, checkBoxes)) > 0,
                        (Files.exists(sampDir) && Files.list(sampDir, ".sampRAF", proj.JAR_STATUS.getValue()).length > 0),
                };
            }
        },
        S6_RUN_PLINK("Create/Run PLINK Files", "", new String[]{"Option 2 must be selected and valid.", "Parsed sample files must already exist."}) {
            @Override
            public void run(Project proj) {
                proj.getLog().report("Running PLINK");
                cnv.manage.PlinkFormat.createPlink(proj, "gwas", null);
                CmdLine.run("plink --file gwas --make-bed --out plink", proj.PROJECT_DIRECTORY.getValue());
                new File(proj.PROJECT_DIRECTORY.getValue()+"genome/").mkdirs();
                CmdLine.run("plink --bfile ../plink --freq", proj.PROJECT_DIRECTORY.getValue()+"genome/");
                CmdLine.run("plink --bfile ../plink --missing", proj.PROJECT_DIRECTORY.getValue()+"genome/");
            }
            @Override
            public boolean[] hasRequirements(Project proj, HashMap<STEPS, JCheckBox> checkBoxes) {
                String sampDir = proj.SAMPLE_DIRECTORY.getValue(false, false);
                return new boolean[]{
                        checkBoxes.get(STEPS.S2_PARSE_SAMPLES).isSelected() && Array.booleanArraySum(STEPS.S2_PARSE_SAMPLES.hasRequirements(proj, checkBoxes)) > 0,
                        (Files.exists(sampDir) && Files.list(sampDir, ".sampRAF", proj.JAR_STATUS.getValue()).length > 0),
                };
            }
        },
        S7_GWAS_QC("Run GWAS QC", "", new String[]{}) {
            @Override
            public void run(Project proj) {
                // TODO Auto-generated method stub
            }
            @Override
            public boolean[] hasRequirements(Project proj, HashMap<STEPS, JCheckBox> checkBoxes) {
                // TODO Auto-generated method stub
                return new boolean[]{false};
            }
        },
        S8_CREATE_PCS("Create Principal Components File", "", new String[]{}) {
            @Override
            public void run(Project proj) {
                // TODO Auto-generated method stub
            }
            @Override
            public boolean[] hasRequirements(Project proj, HashMap<STEPS, JCheckBox> checkBoxes) {
                // TODO Auto-generated method stub
                return new boolean[]{false};
            }
        },
        S9_CREATE_MT_CN_EST("Create Mitochondrial Copy-Number Estimates File", "", new String[]{}) {
            @Override
            public void run(Project proj) {
                // TODO Auto-generated method stub
            }
            @Override
            public boolean[] hasRequirements(Project proj, HashMap<STEPS, JCheckBox> checkBoxes) {
                // TODO Auto-generated method stub
                return new boolean[]{false};
            }
        },
        S10_PENNCNV("Export PennCNV Files", "", new String[]{}) {
            @Override
            public void run(Project proj) {
                // TODO Auto-generated method stub
            }
            @Override
            public boolean[] hasRequirements(Project proj, HashMap<STEPS, JCheckBox> checkBoxes) {
                // TODO Auto-generated method stub
                return new boolean[]{false};
            }
        },
        S11_SHADOW_SAMPLES("Create 'Shadow' Sample Files", "", new String[]{}) {
            @Override
            public void run(Project proj) {
                // TODO Auto-generated method stub
            }
            @Override
            public boolean[] hasRequirements(Project proj, HashMap<STEPS, JCheckBox> checkBoxes) {
                // TODO Auto-generated method stub
                return new boolean[]{false};
            }
        };
        
        public String stepName;
        public String stepDesc;
        public String[] reqs;
        public abstract void run(Project proj);
        public abstract boolean[] hasRequirements(Project proj, HashMap<STEPS, JCheckBox> checkBoxes);
        public String[] getRequirements() {
            return reqs; 
        }
        STEPS(String name, String desc, String[] requirements) {
            this.stepName = name;
            this.stepDesc = desc;
            this.reqs = requirements;
        }

        
    }
    
    public KitAndKaboodle(Project project) {
        this.proj = project;
        this.log = project.getLog();
    }
    
    public void showDialogAndRun() {
        KitAndKaboodleGUI kkGUI = new KitAndKaboodleGUI(this.proj);
        kkGUI.setModal(true);
        kkGUI.setVisible(true);
        
        if (kkGUI.getCancelled() == true) {
            return;
        }
        
        boolean[] options = kkGUI.getSelectedOptions();
        
        for (int i = 0; i < options.length; i++) {
            if (options[i]) {
                STEPS.values()[i].run(proj);
            }
        }
        
    }
    
    
    
}
