
package cnv.filesys;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.EventQueue;
import java.awt.Font;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;
import java.awt.event.KeyEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.EventObject;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;

import javax.swing.AbstractAction;
import javax.swing.ActionMap;
import javax.swing.DefaultCellEditor;
import javax.swing.InputMap;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.KeyStroke;
import javax.swing.ListSelectionModel;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;
import javax.swing.border.EmptyBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableCellEditor;
import javax.swing.table.TableCellRenderer;

import net.miginfocom.swing.MigLayout;
import cnv.LaunchProperties;
import cnv.filesys.Project.DoubleProperty;
import cnv.filesys.Project.FileProperty;
import cnv.filesys.Project.IntegerProperty;
import cnv.filesys.Project.Property;
import cnv.filesys.Project.StringListProperty;
//import cnv.filesys.Project.MultiFileProperty;
import common.Array;
import common.Files;
import common.Grafik;
import common.ext;

public class Configurator extends JFrame {
    
    public static final String[][] ALL_PROPERTY_SETS = new String[][]{
            {
                "Project Name and Locations",
                "PROJECT_NAME",
                "PROJECT_DIRECTORY",
                "DATA_DIRECTORY",
                "SAMPLE_DATA_FILENAME",
                "SAMPLE_DIRECTORY",
                "MARKER_DATA_DIRECTORY",
                "RESULTS_DIRECTORY",
                "DEMO_DIRECTORY",
                "BACKUP_DIRECTORY",
                "ARRAY_TYPE",
                "XY_SCALE_FACTOR"
            },
            {
                "Import",
                "SOURCE_DIRECTORY",
                "SOURCE_FILENAME_EXTENSION",
                "LONG_FORMAT",
                "SOURCE_FILE_DELIMITER",
                "ID_HEADER",
                "PARSE_AT_AT_SYMBOL",
                "MARKER_POSITION_FILENAME",
                "SAMPLE_ALIAS",
                "FID_ALIAS",
                "IID_ALIAS"
            },
            {
                "Global",
                "NUM_THREADS",
                "LOG_LEVEL",
                "CUSTOM_COLOR_SCHEME_FILENAME",
                "CLUSTER_FILTER_COLLECTION_FILENAME",
                "ANNOTATION_FILENAME",
                "AB_LOOKUP_FILENAME",
                "GENETRACK_FILENAME",
                "GC_MODEL_FILENAME",
                "GC_CORRECTION_PARAMETERS_FILENAMES",
                "REFERENCE_GENOME_FASTA_FILENAME"
            },
            {
                "Centroids",
                "SEX_CENTROIDS_MALE_FILENAME",
                "SEX_CENTROIDS_FEMALE_FILENAME",
                "ORIGINAL_CENTROIDS_FILENAME",
                "GENOTYPE_CENTROIDS_FILENAME",
                "CHIMERA_CENTROIDS_FILENAME",
                "CUSTOM_CENTROIDS_FILENAME"
            },
            {
                "DataExport",
                "PEDIGREE_FILENAME",
                "FILTERED_MARKERS_FILENAME",
                "SAMPLE_SUBSET_FILENAME",
                "TARGET_MARKERS_FILENAMES",
                "GC_THRESHOLD",
            },
            {
                "MosaicPlot",
                "MOSAIC_RESULTS_FILENAME",
                "MOSAIC_COLOR_CODES_FILENAME",
            },
            {
                "Data Cleaning",
                "SAMPLE_QC_FILENAME",
                "SEXCHECK_RESULTS_FILENAME",
                "LRRSD_CUTOFF",
                "MARKER_METRICS_FILENAME",
                "MARKER_EXCLUSION_CRITERIA_FILENAME",
                "MARKER_REVIEW_CRITERIA_FILENAME",
                "MARKER_COMBINED_CRITERIA_FILENAME",
                "GENOME_CLUSTER_FILENAME"
            },
            {
                "CNV Files",
                "CNV_FILENAMES",
                "HMM_FILENAME",
                "CUSTOM_PFB_FILENAME",
            },
            {
                "CompPlot",
                "REGION_LIST_FILENAMES"
            },
            {
                "Trailer",
                "INDIVIDUAL_CNV_LIST_FILENAMES",
                "WINDOW_AROUND_SNP_TO_OPEN_IN_TRAILER"
            },
            {
                "ScatterPlot",
                "DISPLAY_MARKERS_FILENAMES",
                "SHIFT_SEX_CHR_COLORS_YESNO",
                "BLAST_ANNOTATION_FILENAME",
                "BLAST_PROPORTION_MATCH_FILTER",
            },
            {
                "TwoDPlot",
                "TWOD_LOADED_FILENAMES",
                "TWOD_LOADED_VARIABLES"
            },
            {
                "ForestPlot",
                "FOREST_PLOT_FILENAMES"
            },
            {
                "QQ-plot",
                "QQ_FILENAMES",
                "DISPLAY_QUANTILES",
                "DISPLAY_STANDARD_QQ",
                "DISPLAY_ROTATED_QQ",
                "QQ_MAX_NEG_LOG10_PVALUE"
            },
            {
                "PennCNV",
                "PENNCNV_EXECUTABLE_DIRECTORY",
                "PENNCNV_DATA_DIRECTORY",
                "PENNCNV_RESULTS_DIRECTORY",
                "PENNCNV_GZIP_YESNO",
            }, 
            {
                "CytoSpecific",
                "UNREPORTED_CNP_FILENAME",
                "COMMON_CNP_FILENAME",
                "REPORTED_CNP_FILENAME"
            },
            {
                "PC Intensity Correction",
                "INTENSITY_PC_FILENAME",
                "INTENSITY_PC_NUM_COMPONENTS"
            },
            {
                "Optimization parameters",
                "MAX_MEMORY_USED_TO_LOAD_MARKER_DATA",
                "MAX_MARKERS_LOADED_PER_CYCLE"
            },
            {
                "Plink Directory/Filename Roots (edit to remove extension)",
                "PLINK_DIR_FILEROOTS"
            }
    };
    String[] hiddenProperties = new String[]{
            "PROJECT_PROPERTIES_FILENAME",
            "MARKERSET_FILENAME",
            "MARKERLOOKUP_FILENAME",
            "SAMPLELIST_FILENAME",
            "JAR_STATUS",
            "STRATIFICATION_RESULTS_FILENAMES"
    };
    String[] uneditableProperties = new String[] {
            "ARRAY_TYPE",
            "XY_SCALE_FACTOR"
    };
    
	private static final long serialVersionUID = 1L;

	private JPanel contentPane;
	
	private JTable table;

	private Project proj;
	
	private ArrayList<Integer> labelRows = new ArrayList<Integer>();
	
	private abstract class InputValidator {
	    abstract Object processNewValue(Object newValue, Object oldValue);
	    abstract boolean acceptNewValue(Object newValue);
	}
	
	private HashMap<String, InputValidator> validators = new HashMap<String, Configurator.InputValidator>() {
        private static final long serialVersionUID = 1L;
    {
	    put("PLINK_DIR_FILEROOTS", new InputValidator() {
            @Override
            Object processNewValue(Object newValue, Object oldValue) {
                if (newValue instanceof File[]) {
                    if (((File[])newValue).length == 0 || (((File[])newValue).length == 1 && ((File[])newValue)[0].getPath().equals(""))) {
                        return newValue;
                    }
                    File[] newFiles = new File[((File[])newValue).length];
                    for (int i = 0; i < ((File[])newValue).length; i++) {
                        String newPath = ext.rootOf(((File[])newValue)[i].getAbsolutePath(), false);
                        newFiles[i] = new File(newPath);
                    }
                    return newFiles;
                } else {
                    System.out.println("ERROR - new value should be a File[] (array)!");
                    return oldValue;
                }
            }
            @Override
            boolean acceptNewValue(Object newValue) {
                if (newValue instanceof File[]) {
                    File[] files = (File[]) newValue;
                    if (files.length == 0 || (files.length == 1 && files[0].getPath().equals(""))) {
                        return true;
                    }
                    for (File f : files) {
                        if (!hasAllPLINKFiles(f)) {
                            proj.message("Error - couldn't find all files {.bim/.fam/.bed} for PLINK root [" + ext.rootOf(f.getPath(), true) + "]");
                            return false;
                        }
                    }
                    return true;
                } else {
                    System.out.println("ERROR - new value should be a File[] (array)!");
                    return false;
                }
            }
            private boolean hasAllPLINKFiles(File f) {
                String[] files = f.getParentFile().list(new FilenameFilter() {
                    @Override
                    public boolean accept(File dir, String name) {
                        return name.endsWith(".fam") || name.endsWith(".bim") || name.endsWith(".bed");
                    }
                });
                if (files.length < 3) return false;
                boolean hasBed = false, hasFam = false, hasBim = false;
                for (String file : files) {
                    if (file.endsWith(".fam")) {
                        hasFam = true;
                    } else if (file.endsWith(".bim")) {
                        hasBim = true;
                    } else if (file.endsWith(".bed")) {
                        hasBed = true;
                    }
                }
                return hasBed && hasFam && hasBim;
            }
        });
	}};
	
	/**
	 * Launch the application.
	 */
	public static void main(String[] args) {
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					Configurator frame = new Configurator(new Project(), ALL_PROPERTY_SETS);
					frame.setVisible(true);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		});
	}

	/**
	 * Create the frame.
	 */
	public Configurator(Project project, String[][] propertySets) {
		setTitle("Genvisis - " + project.PROJECT_NAME.getValue() + " - Project Configuration");
		setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
		setBounds(100, 100, 700, 800);
		contentPane = new JPanel();
		contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
		contentPane.setLayout(new BorderLayout(0, 0));
		setContentPane(contentPane);
		this.proj = project;
		
		JPanel panel = new JPanel();
		getContentPane().add(panel, BorderLayout.SOUTH);
		
        JButton notepad = new JButton("Edit with Notepad");
        notepad.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                proj.getLog().report("Launching notepad...");
                try {
                    /*Process p = */Runtime.getRuntime().exec("C:\\Windows\\System32\\Notepad.exe \""+proj.getPropertyFilename()+"\"");
                    Configurator.this.setVisible(false);
                    // TODO update properties in Project and Configurator - they may have changed
                } catch (IOException ioe) {
                    proj.getLog().reportError("Error - failed to open Notepad");
                }
            }
        });
        panel.add(notepad);
        
        
        boolean includeNotepad = false;
        if (System.getProperty("os.name").startsWith("Windows")) {
            includeNotepad = Files.programExists("notepad.exe");
        }
	    notepad.setVisible(includeNotepad);
        
		JButton btnSave = new JButton("Save");
		btnSave.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				Configurator.this.save();
			}
		});
		panel.add(btnSave);
		
		JButton btnSaveAndClose = new JButton("Save and Close");
		btnSaveAndClose.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
			    doClose(true, false);
			}
		});
		panel.add(btnSaveAndClose);
		
		JButton btnClose = new JButton("Close");
		btnClose.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
			    doClose(false, true);
			}
		});
		panel.add(btnClose);
		
		JScrollPane scrollPane = new JScrollPane();
		getContentPane().add(scrollPane, BorderLayout.CENTER);
		
		final JCheckBox rendererChkBx = new JCheckBox();
		rendererChkBx.setBackground(Color.WHITE);
		rendererChkBx.setHorizontalAlignment(SwingConstants.TRAILING);
		rendererChkBx.setFont(rendererChkBx.getFont().deriveFont(16));
		Grafik.scaleCheckBoxIcon(rendererChkBx);
		rendererChkBx.setBorder(null);
		final JCheckBox editorChkBx = new JCheckBox();  
		editorChkBx.setBackground(Color.WHITE);
		editorChkBx.setHorizontalAlignment(SwingConstants.TRAILING);
		editorChkBx.setFont(editorChkBx.getFont().deriveFont(16));
		Grafik.scaleCheckBoxIcon(editorChkBx);
		final DefaultCellEditor boolEditor = new DefaultCellEditor(editorChkBx);
		final JSpinner rendererSpinner = new JSpinner();
		rendererSpinner.setBorder(null);
		final SpinnerEditor numberEditor = new SpinnerEditor();
		final JPanel fileRenderer = new JPanel(new MigLayout("ins -1 0 0 0, hidemode 3", "[grow][]", ""));//new BorderLayout());
		final JLabel fileLabel = new JLabel();
		final JButton fileBtn2 = new JButton("...");
		final JButton fileAddBtn2 = new JButton(" + ");
		fileAddBtn2.setFont(fileAddBtn2.getFont().deriveFont(Font.BOLD));
		fileBtn2.setMargin(new Insets(1, 1, 0, 1));
		fileAddBtn2.setMargin(new Insets(1, 1, 0, 1));
		fileLabel.setBackground(Color.WHITE);
		fileRenderer.setBackground(Color.WHITE);
		fileRenderer.add(fileLabel, "cell 0 0");
		fileRenderer.add(fileAddBtn2, "cell 1 0, width 21px, split 1");
		fileRenderer.add(fileBtn2, "cell 1 0, width 21px");
		final FileChooserCellEditor fileEditor = new FileChooserCellEditor(true, proj.getProperty(proj.PROJECT_DIRECTORY));
		final DefaultCellEditor stringEditor = new DefaultCellEditor(new JTextField()) {
			private static final long serialVersionUID = 1L;
			@Override
			public Object getCellEditorValue() {
				return ((String) super.getCellEditorValue()).trim();
			}
			{
				this.editorComponent.addFocusListener(new FocusListener() {
					@Override
					public void focusLost(FocusEvent e) {
						stopCellEditing();
					}
					@Override
					public void focusGained(FocusEvent e) {
					}
				});
			}
		};
		final DefaultCellEditor stringListEditor = new DefaultCellEditor(new JTextField()) {
		    private static final long serialVersionUID = 1L;
		    @Override
		    public Object getCellEditorValue() {
		        return ((String) super.getCellEditorValue()).trim().split(";");
		    }
		    @Override
		    public Component getTableCellEditorComponent(JTable table, Object value, boolean isSelected, int row, int column) {
		        String[] val = (String[]) value;
		        String valueString = Array.toStr(val, ";");
		        return super.getTableCellEditorComponent(table, valueString, isSelected, row, column);
		    }
		    {
		        this.editorComponent.addFocusListener(new FocusListener() {
		            @Override
		            public void focusLost(FocusEvent e) {
		                stopCellEditing();
		            }
		            @Override
		            public void focusGained(FocusEvent e) {
		            }
		        });
		    }
		};

		final DefaultTableCellRenderer renderer = new DefaultTableCellRenderer() {
			private static final long serialVersionUID = 1L;

			@Override
			public Component getTableCellRendererComponent(final JTable table, Object value, boolean isSelected, boolean hasFocus, final int row, final int column) {
				String projDir = proj.getProperty(proj.PROJECT_DIRECTORY);
				String tempKey = (String) table.getValueAt(row, 0);
				Component returnComp;
				if (labelRows.contains(row)) {
				    returnComp = super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, column);
                    ((JComponent) returnComp).setToolTipText(null);
                    ((JComponent) returnComp).setFont(((JComponent) returnComp).getFont().deriveFont(Font.BOLD, 14f));
                    return returnComp;
				}
				String desc = proj.getProperty(tempKey).getDescription();
				if (value instanceof Boolean) {
					rendererChkBx.setSelected(((Boolean) value).booleanValue());
					returnComp = rendererChkBx;
				} else if (value instanceof File) {
					String txt = ((File) value).getPath();
					if (((File) value).isDirectory() || ((FileProperty)proj.getProperty(tempKey)).isDir) {
						txt = ext.verifyDirFormat(txt);
					} else {
						txt = ext.replaceAllWith(txt, "\\", "/");
					}
					if (txt.startsWith(projDir)) {
						txt = txt.substring(projDir.length());
						if (txt.trim().length() == 0) {
						    txt = projDir;
						}
					}
					boolean isDefault = txt.equals(proj.getProperty(tempKey).getDefaultValueString());
					if (!isDefault && !Files.exists(((File) value).getPath())) {
					    txt = "<font color='red'>" + txt + "</font>";
					} else {
					    txt = "<font>" + txt + "</font>";
					}
					fileLabel.setText("<html><nobr>"+txt+"</nobr></html>");
					fileAddBtn2.setVisible(false);
                    fileBtn2.setVisible(true);
					returnComp = fileRenderer;
				} else if (value instanceof File[]) {
					StringBuilder sb = new StringBuilder("<html><nobr>");
					boolean exists = false;
					String defVal = proj.getProperty(tempKey).getDefaultValueString();
					if (((File[])value).length > 0) {
						String txt = ((File[])value)[0].getPath();
						exists = Files.exists(txt);
						if (((File[]) value)[0].isDirectory()) {
							txt = ext.verifyDirFormat(txt);
						} else {
							txt = ext.replaceAllWith(txt, "\\", "/");
						}
						if (txt.startsWith(projDir)) {
							txt = txt.substring(projDir.length());
						}
                        if (!exists && !defVal.contains(txt)) {
						    sb.append("<font color='red'>");
						} else {
						    exists = true;
						}
						sb.append(txt);
                        if (!exists) {
                            sb.append("</font>");
                        }
						for (int i = 1; i < ((File[])value).length; i++) {
							txt = ((File[])value)[i].getPath();
							exists = Files.exists(txt);
							if (((File[]) value)[i].isDirectory()) {
								txt = ext.verifyDirFormat(txt);
							} else {
								txt = ext.replaceAllWith(txt, "\\", "/");
							}
							if (txt.startsWith(projDir)) {
								txt = txt.substring(projDir.length());
							}
							sb.append(";");
                            if (!exists && !defVal.contains(txt)) {
                                sb.append("<font color='red'>");
                            } else {
                                exists = true;
                            }
                            sb.append(txt);
                            if (!exists) {
                                sb.append("</font>");
                            }
						}
						sb.append("</nobr></html>");
					}
                    fileAddBtn2.setVisible(true);
                    fileBtn2.setVisible(false);
					fileLabel.setText(sb.toString());
					returnComp = fileRenderer;
				} else if (value instanceof String[]) {
				    StringListProperty prop = proj.getProperty((String) table.getValueAt(row, 0));
                    StringBuilder sb = new StringBuilder("<html><nobr>");
                    String[] values = (String[]) value;
                    for (int i = 0; i < values.length; i++) {
                        if (values[i].equals("")) {
                            continue;
                        }
                        boolean exists = (!prop.isFile && !prop.isDir) || Files.exists(values[i]);
                        if (!exists) {
                            sb.append("<font color='red'>");
                        }
                        if (values[i].startsWith(projDir)) {
                            sb.append(values[i].substring(projDir.length()));
                        } else {
                            sb.append(values[i]);
                        }
                        if (!exists) {
                            sb.append("</font>");
                        }
                        if (i < values.length - 1) {
                            sb.append(";");
                        }
                    }
                    sb.append("</nobr></html>");
                    if (prop.isFile || prop.isDir) {
                        fileLabel.setText(sb.toString());
                        fileAddBtn2.setVisible(true);
                        fileBtn2.setVisible(false);
                        returnComp = fileRenderer;
                    } else {
                        returnComp = super.getTableCellRendererComponent(table, sb.toString(), isSelected, hasFocus, row, column);
                    }
			    } else if (value instanceof Number) {
					String propKey = (String) table.getModel().getValueAt(row, 0);
					Object propVal = table.getModel().getValueAt(row, 1);
					setByPropertyKey(rendererSpinner, proj, propKey, propVal);
					returnComp = rendererSpinner;
			    } else if (value instanceof Enum<?>) {
//			        System.out.println("Found ENUM; values: [" + Array.toStr(((Enum) value).getDeclaringClass(), ",") + "]");
			        // use string renderer for enums
                    returnComp = super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, column);
				} else {
				    if (column == 0) {
				        returnComp = super.getTableCellRendererComponent(table, "            "+value, isSelected, hasFocus, row, column);
					    ((JComponent) returnComp).setFont(((JComponent) returnComp).getFont().deriveFont(Font.PLAIN, 12f));
					} else {
					    returnComp = super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, column);
					}
				}
				if (!"".equals(desc)) {
					((JComponent) returnComp).setToolTipText(desc);
				} else {
					((JComponent) returnComp).setToolTipText(null);
				}
				return returnComp;
			}
		};

        final HashSet<String> disabledKeys = new HashSet<String>();
        for (String key : uneditableProperties) {
            disabledKeys.add(key);
        }
        
		final Color bgColor2 = new Color(184,207,229);
		table = new JTable() {
			private static final long serialVersionUID = 1L;
			volatile private int editRow = -1;
			private void setEditRow(int row) {
				editRow = row;
			}
			private int getEditRow() {
				return editRow;
			}
			public javax.swing.table.TableCellEditor getCellEditor(int row, int column) {
				if (column == 0) {
					return super.getCellEditor(row, column);
				}
				
				String propKey = (String) this.getModel().getValueAt(row, 0);
				Object propVal = this.getModel().getValueAt(row, 1);
				
				javax.swing.table.TableCellEditor editor;
				if (propVal instanceof File || propVal instanceof File[]) {
					if (getEditRow() != row) {
						setEditRow(row);
						fileEditor.reset();
					}
					fileEditor.setValue(this.getValueAt(row, column));
					editor = fileEditor;
				} else if (propVal instanceof Boolean) {
					 editor = boolEditor;
				} else if (propVal instanceof Number) {
					setByPropertyKey(numberEditor.spinner, proj, propKey, propVal);
					editor = numberEditor;
				} else if (propVal instanceof String) {
					editor = stringEditor;
				} else if (propVal instanceof String[]) {
				    StringListProperty prop = ((StringListProperty)proj.getProperty((String) table.getValueAt(row, 0)));
				    if (prop.isFile || prop.isDir) {
				        String[] valStrs = (String[]) propVal;
				        File[] vals = new File[valStrs.length];
				        for (int i = 0; i < valStrs.length; i++) {
				            vals[i] = new File(valStrs[i]);
				        }
				        if (getEditRow() != row) {
	                        setEditRow(row);
	                        fileEditor.reset();
	                    }
	                    fileEditor.setValue(vals);
	                    editor = fileEditor;
				    } else {
				        editor = stringListEditor;
				    }
				} else if (propVal instanceof Enum<?>) {
				    @SuppressWarnings("rawtypes")
                    Object[] values = ((Enum)propVal).getDeclaringClass().getEnumConstants();
                    DefaultCellEditor enumEditor = new DefaultCellEditor(new JComboBox<Object>(values));
				    editor = enumEditor;
				} else {
//				    System.out.println("Not found: Class<" + propVal.getClass().getName() + ">");
				    editor = super.getCellEditor(row, column);
				}
				
				return editor;
			}
			@Override
			public TableCellRenderer getCellRenderer(int row, int column) {
				return renderer;
			}
			@Override
			public Component prepareRenderer(TableCellRenderer renderer, int row, int column) {
				Component superComp = super.prepareRenderer(renderer, row, column);
				if (super.isCellSelected(row, column)) {
					superComp.setBackground(bgColor2);
				} else {
					superComp.setBackground(Color.WHITE);
				}
				if (disabledKeys.contains(getValueAt(row, 0))) {
				    superComp.setEnabled(false);
				} else {
				    superComp.setEnabled(true);
				}
				return superComp;
			}
			@Override
		    public void editingStopped(ChangeEvent e) {
		        // Take in the new value
		        TableCellEditor editor = getCellEditor();
		        if (editor != null) {
		            Object value = editor.getCellEditorValue();
	                String propKey = (String) this.getModel().getValueAt(this.editingRow, 0);
		            Object oldValue = getValueAt(editingRow, editingColumn);
	                if (acceptNewValue(propKey, value)) {
    		            value = processNewValue(propKey, value, oldValue);
    		            setValueAt(value, editingRow, editingColumn);
	                } else {
	                    setValueAt(oldValue, editingRow, editingColumn);
	                }
		            removeEditor();
		        }
		    }
		};

		DefaultTableModel model = new DefaultTableModel(new String[]{"Property Name", "Property Value"}, 0) {
			private static final long serialVersionUID = 1L;
			@Override
			public boolean isCellEditable(int row, int column) { 
			    return column != 0 && !labelRows.contains(row) && super.isCellEditable(row, column) && !disabledKeys.contains(getValueAt(row, 0)); 
		    }
			
		};
		
		int count = 0;
		HashSet<String> allKeys = new HashSet<String>();
		for (String key : proj.getPropertyKeys()) {
		    allKeys.add(key);
		}
		for (String[] propKeySet : propertySets) {
		    String setName = propKeySet[0];
		    model.addRow(new Object[]{setName, ""});
		    labelRows.add(count);
		    count++;
		    for (int i = 1; i < propKeySet.length; i++) {
		        boolean removed = allKeys.remove(propKeySet[i]);
		        if (removed) {
    		        Object[] values = parseProperty(proj, propKeySet[i]);
    	            model.addRow(values);
    	            count++;
		        } else {
		            proj.getLog().reportError("Unknown key found: " + propKeySet[i]);
		        }
		    }
		}
		HashSet<String> excludedKeys = new HashSet<String>();
		for (String key : hiddenProperties) {
		    excludedKeys.add(key);
		}
		ArrayList<String> leftovers = new ArrayList<String>();
		for (String leftoverKey : allKeys) {
		    if (!excludedKeys.contains(leftoverKey)) {
		        leftovers.add(leftoverKey);
		    }
		}
		if (!leftovers.isEmpty()) {
		    proj.getLog().report("Found " + leftovers.size() + " unknown keys: " + leftovers.toString());
		}
		
		table.setModel(model);
		table.getColumnModel().getColumn(0).setPreferredWidth(200);
		table.getColumnModel().getColumn(1).setPreferredWidth(200);
		table.setFillsViewportHeight(true);
		table.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
		table.setRowHeight(24);
		table.setRowMargin(2);
		table.setShowVerticalLines(false);
		table.setShowHorizontalLines(false);
		
		
		InputMap inMap = table.getInputMap(JTable.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT);
		ActionMap actMap = table.getActionMap();
		KeyStroke spaceKey = KeyStroke.getKeyStroke(KeyEvent.VK_SPACE, 0);
		AbstractAction fileSelectAction = new AbstractAction() {
			private static final long serialVersionUID = 1L;
			public void actionPerformed(ActionEvent evt) {
				fileEditor.keyed = true;
		        table.changeSelection(table.getSelectedRow(), 1, false, false);
		        if (!table.editCellAt(table.getSelectedRow(), 1)) {
//		            JOptionPane.showMessageDialog(table, "Failed to start cell editing");
		        }
				fileEditor.keyed = false;
		    }
		};
		inMap.put(spaceKey, "Action.spacebar");
		actMap.put("Action.spacebar", fileSelectAction);

		table.setSurrendersFocusOnKeystroke(true);

		scrollPane.setViewportView(table);
		
		addWindowListener(new WindowAdapter() {
		    @Override
		    public void windowClosing(WindowEvent e) {
		        doClose(false, true);
		    }
        });
		
	}
	
	protected void doClose(boolean save, boolean promptChanges) {
	    HashMap<String, String> changes = extract();
	    if (changes.size() > 0) {
	        if (promptChanges) {
    	        StringBuilder message = new StringBuilder("The following properties have been changed.  Would you like to save your changes?");
    	        int cnt = 0;
    	        for (String key : changes.keySet()) {
    	            if (cnt == 10) {
    	                message.append("\n... and ").append(changes.size() - cnt).append(" additional changes...");
    	                break;
    	            }
    	            message.append("\n").append(key);
    	            cnt++;
    	        }
    	        String title = "Save changes?";
    	        int opt = JOptionPane.showConfirmDialog(Configurator.this, message, title, JOptionPane.YES_NO_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE);
    	        if (opt == JOptionPane.CLOSED_OPTION || opt == JOptionPane.CANCEL_OPTION) {
    	            return;
    	        } else if (opt == JOptionPane.YES_OPTION) {
    	            save(changes);
    	        }
	        } else if (save) {
	            save(changes);
	        }
	    }
        Configurator.this.setVisible(false);
        Configurator.this.dispose();
    }

    private HashMap<String, String> extract() {
	    table.editingStopped(new ChangeEvent(table));
        String projectsDir = new LaunchProperties(LaunchProperties.DEFAULT_PROPERTIES_FILE).getProperty(LaunchProperties.PROJECTS_DIR);
        String currProjDir = proj.getProperty(proj.PROJECT_DIRECTORY);
        int rowCount = table.getRowCount();
        
        HashMap<String, String> newValues = new HashMap<String, String>();
        
        for (int i = 0; i < rowCount; i++) {
            if (labelRows.contains(i)) {
                continue;
            }
            String key = ((String) table.getValueAt(i, 0)).trim();
            Object rawValue = table.getValueAt(i, 1);
            
            String value = "";
            if (rawValue instanceof File[]) {
                File[] set = (File[]) rawValue;
                if (set.length > 0) {
                    value = set[0].getPath();
                    if (!set[0].exists()) {
                        value = ((StringListProperty)proj.getProperty(key)).isDir ? ext.verifyDirFormat(value) : ext.replaceAllWith(value, "\\", "/"); 
                    } else {
                        value = set[0].isDirectory() ? ext.verifyDirFormat(value) : ext.replaceAllWith(value, "\\", "/");
                    }
                    value = set[0].isDirectory() ? ext.verifyDirFormat(value) : ext.replaceAllWith(value, "\\", "/");
                    if (value.startsWith(projectsDir)) {
                        value = value.substring(projectsDir.length());
                    } else if (value.startsWith(currProjDir)) {
                        value = value.substring(currProjDir.length());
                    }
                    for (int k = 1; k < set.length; k++) {
                        String fNm = set[k].getPath();
                        if (!set[k].exists()) {
                            fNm = ((StringListProperty)proj.getProperty(key)).isDir ? ext.verifyDirFormat(fNm) : ext.replaceAllWith(fNm, "\\", "/");
                        } else {
                            fNm = set[k].isDirectory() ? ext.verifyDirFormat(fNm) : ext.replaceAllWith(fNm, "\\", "/");
                        }
                        if (fNm.startsWith(projectsDir)) {
                            fNm = fNm.substring(projectsDir.length());
                        } else if (fNm.startsWith(currProjDir)) {
                            fNm = fNm.substring(currProjDir.length());
                        }
                        value += ";" + fNm;
                    }
                }
            } else if (rawValue instanceof File) { 
                File set = (File) rawValue;
                value = set.getPath();
                if (!set.exists()) {
                    value = ((FileProperty)proj.getProperty(key)).isDir ? ext.verifyDirFormat(value) : ext.replaceAllWith(value, "\\", "/");
                } else {
                    value = set.isDirectory() ? ext.verifyDirFormat(value) : ext.replaceAllWith(value, "\\", "/");
                }
                if (!key.equals(proj.SOURCE_DIRECTORY.getName()) && !key.equals(proj.PROJECT_DIRECTORY.getName())) {
                    if (value.startsWith(projectsDir)) {
                        value = value.substring(projectsDir.length());
                    } else if (value.startsWith(currProjDir)) {
                        value = value.substring(currProjDir.length());
                    }
                }
            } else if (rawValue instanceof String[]) {
                value = "";//((String[])rawValue)[0];
                for (int k = 0; k < ((String[])rawValue).length; k++) {
                    String tempValue = ((String[])rawValue)[k]; 
                    if (tempValue.startsWith(projectsDir)) {
                        value += tempValue.substring(projectsDir.length());
                    } else if (tempValue.startsWith(currProjDir)) {
                        value += tempValue.substring(currProjDir.length());
                    } else {
                        value += tempValue;
                    }
                    if (k < ((String[])rawValue).length - 1) {
                        value += ";";
                    }
                }
            } else {
                if (rawValue.toString().startsWith(projectsDir)) {
                    value += rawValue.toString().substring(projectsDir.length());
                } else if (rawValue.toString().startsWith(currProjDir)) {
                        value += rawValue.toString().substring(currProjDir.length());
                } else {
                    value += rawValue.toString();
                }
            }
            
            HashSet<String> ignorePrefix = new HashSet<String>();
            ignorePrefix.add("PROJECT_DIRECTORY");
            ignorePrefix.add("SOURCE_DIRECTORY");
            
            String t1 = proj.getProperty(key).getValueString();
            if (!ignorePrefix.contains(proj.getProperty(key).getName())) {
                t1 = t1.replaceAll(proj.PROJECT_DIRECTORY.getValue(), "");
            }
            boolean same = t1.equals(value);
            if (!same) {
                newValues.put(key, value);
            }
        }
        
        return newValues;
    }
	
    private void save(HashMap<String, String> newValues) {
        for (Entry<String, String> pair : newValues.entrySet()) {
            proj.setProperty(pair.getKey(), pair.getValue());
        }
        proj.saveProperties();
    }
    
	private void save() {
	    HashMap<String, String> newValues = extract();
	    save(newValues);
	}

	private void setByPropertyKey(JSpinner spinner, Project proj, String propKey, Object propVal) {
		Property<?> prop = proj.getProperty(propKey);
		
		if (prop instanceof DoubleProperty) {
			double value = "".equals(propVal) ? 0.0 : Double.valueOf(propVal.toString()).doubleValue();
			double min = ((DoubleProperty) prop).min;
			double max = ((DoubleProperty) prop).max;
			double stepSize = 1.0;
			for (int i = 0; i < ext.getNumSigFig(value); i++) {
				stepSize /= 10.0;
			}
			spinner.setModel(new SpinnerNumberModel(value, min, max, stepSize));
		} else if (prop instanceof IntegerProperty) {
			int value = "".equals(propVal) ? 0 : Integer.valueOf(propVal.toString()).intValue();
			int min = ((IntegerProperty) prop).min;
			int max = ((IntegerProperty) prop).max;
			int stepSize = 1;
			spinner.setModel(new SpinnerNumberModel(value, min, max, stepSize));
		}
		
	}
	
	private Object[] parseProperty(Project proj, String propKey) {
		Object[] keyVal = new Object[2];
		keyVal[0] = propKey;
		Property<?> prop = proj.getProperty(propKey);
		keyVal[1] = prop.getValue();
		if (prop instanceof FileProperty) {
			keyVal[1] = new File(((FileProperty)prop).getValue());
		} else if (prop instanceof StringListProperty) {
		    if (((StringListProperty)prop).isDir || ((StringListProperty)prop).isFile) { 
		        File[] fs = new File[((StringListProperty)prop).getValue().length];
		        for (int i = 0; i < fs.length; i++) {
		            fs[i] = new File(((StringListProperty)prop).getValue()[i]);
		        }
		        keyVal[1] = fs;
		    }
		}
		
		return keyVal;	
	}
	
	private boolean acceptNewValue(String propKey, Object newValue) {
	    InputValidator validator = validators.get(propKey);
	    return validator == null || validator.acceptNewValue(newValue);
	}
	
	private Object processNewValue(String propKey, Object newValue, Object oldValue) {
        InputValidator validator = validators.get(propKey);
        return validator == null ? newValue : validator.processNewValue(newValue, oldValue);
	}
	
	
	class FileChooserCellEditor extends DefaultCellEditor implements TableCellEditor {
		// from http://stackoverflow.com/questions/15644319/jfilechooser-within-a-jtable
		private static final long serialVersionUID = 1L;
		/** Number of clicks to start editing */
	    private static final int CLICK_COUNT_TO_START = 1;
	    /** meta panel */
	    private JPanel panel;
	    /** static display */
	    private JTextField label;	
	    /** Editor component */
	    private JButton buttonReplace;
	    /** Editor component */
	    private JButton buttonAppend;
	    /** File chooser */
	    private JFileChooser fileChooser;
	    /** Selected file(s) */
	    private Object value;
//	    private JTextField textField;
	    volatile boolean keyed = false;
	    volatile boolean isMulti = false;
	    volatile int keyedCount = 0;
	    String defaultLocation;
	    volatile boolean editing = false;
	    JTable table;
	    /**
	     * Constructor.
	     */
	    public FileChooserCellEditor(boolean editable, String defaultLocation) {
	        super(new JTextField());
	        this.defaultLocation = defaultLocation;
	        label = (JTextField) this.editorComponent;
	        label.setEditable(editable);
//	        textField = (JTextField) this.editorComponent;
//	        textField.setEditable(editable);
	        setClickCountToStart(CLICK_COUNT_TO_START);
	        setBackground(Color.WHITE);
	        
	        panel = new JPanel(new MigLayout("insets -1 0 0 0, hidemode 3", "[grow][]", ""));//new BorderLayout());
	        panel.setBackground(Color.WHITE);
	        
	        // Using a JButton as the editor component
//	        label = new JTextField();
	        label.setBackground(Color.WHITE);
	        panel.addFocusListener(new FocusListener() {
                @Override
                public void focusLost(FocusEvent e) {
                    editing = false;
                    stopCellEditing();
                    fireEditingStopped();
                    if (table != null) {
                        ((DefaultTableModel)table.getModel()).fireTableDataChanged();
                    }
                }
                @Override
                public void focusGained(FocusEvent e) {
                    editing = true;
                }
            });
	        label.addFocusListener(new FocusListener() {
	            @Override
	            public void focusLost(FocusEvent e) {
	                if (editing) return;
	                String newText = label.getText();
	                
	                Object newValue = value;
	                if (isMulti) {
	                    String[] pts = newText.split(";");
	                    File[] newFiles = new File[pts.length];
	                    for (int i = 0; i < pts.length; i++) {
	                        newFiles[i] = new File(pts[i]);
	                    }
	                    newValue = newFiles;
	                } else {
	                    newValue = new File(newText);
	                }
	                setValue(newValue);
	                
	                editing = false;
	                stopCellEditing();
	                fireEditingStopped();
	                if (table != null) {
	                    ((DefaultTableModel)table.getModel()).fireTableDataChanged();
	                }
	            }
	            @Override
	            public void focusGained(FocusEvent e) {
	                editing = true;
	            }
	        });
	        this.buttonReplace = new JButton("...");
	        this.buttonAppend = new JButton(" + ");
	        this.buttonAppend.setFont(this.buttonAppend.getFont().deriveFont(Font.BOLD));
	        this.buttonReplace.setMargin(new Insets(1, 1, 0, 1));
	        this.buttonAppend.setMargin(new Insets(1, 1, 0, 1));
	        
	        panel.add(label, "cell 0 0, grow");
	        panel.add(this.buttonAppend, "cell 1 0, split 1");
	        panel.add(this.buttonReplace, "cell 1 0");

	        panel.setBackground(Color.WHITE);
	        // Dialog which will do the actual editing
	        fileChooser = new JFileChooser(defaultLocation);
	    }

	    @Override
	    public Object getCellEditorValue() {
	        String newLoc = label.getText();
	        Object newValue;
	        if (isMulti) {
	            String[] pts = newLoc.split(";");
	            newValue = new File[pts.length];
	            for (int i = 0; i < pts.length; i++) {
	                if (!"".equals(pts[i]) && !pts[i].startsWith(".") && !pts[i].startsWith("/") && pts[i].indexOf(":") == -1) {
	                    ((File[])newValue)[i] = new File(defaultLocation + pts[i]);
	                } else {
	                    ((File[])newValue)[i] = new File(pts[i]);
	                }
	            }
	        } else {
	            if (!"".equals(newLoc) && !newLoc.startsWith(".") && !newLoc.startsWith("/") && newLoc.indexOf(":") == -1) {
	                newValue = new File(defaultLocation + newLoc);
	            } else {
	                newValue = new File(newLoc);
	            }
	        }
	        return newValue;
	    }
	    
	    private void setValue(Object val) {
 	    	this.value = val;
 	    	StringBuilder labelText = new StringBuilder();
 	    	if (value instanceof File) {
                if (((File) value).isDirectory()) {
                    String pathStr = ext.verifyDirFormat(((File) value).getPath());
                    if (pathStr.startsWith(defaultLocation)) {
                        pathStr = pathStr.substring(defaultLocation.length());
                        if (pathStr.length() == 0) {
                            // semi-hack for PROJECT_DIRECTORY
                            pathStr = defaultLocation;
                        }
                    }
                    labelText.append(pathStr);
                } else {
                    String pathStr = ext.replaceAllWith(((File) value).getPath(), "\\", "/");
                    if (pathStr.startsWith(defaultLocation)) {
                        pathStr = pathStr.substring(defaultLocation.length());
                    }
                    labelText.append(pathStr);
                }
            } else if (value instanceof File[]) {
                File[] files = (File[]) value;
                if (files.length > 0) {
                    for (int i = 0; i < files.length; i++) {
                        if (files[i].isDirectory()) {
                            String pathStr = ext.verifyDirFormat(files[i].getPath());
                            if (pathStr.startsWith(defaultLocation)) {
                                pathStr = pathStr.substring(defaultLocation.length());
                                if (pathStr.length() == 0) {
                                    // semi-hack for PROJECT_DIRECTORY
                                    pathStr = defaultLocation;
                                }
                            }
                            labelText.append(pathStr);
                        } else {
                            String pathStr = ext.replaceAllWith(files[i].getPath(), "\\", "/");
                            if (pathStr.startsWith(defaultLocation)) {
                                pathStr = pathStr.substring(defaultLocation.length());
                            }
                            labelText.append(pathStr);
                        }
                        if (i < files.length - 1) {
                            labelText.append(";");
                        }
                    }
                }
            }
            if (labelText.toString().startsWith(defaultLocation)) {
                labelText = new StringBuilder(labelText.substring(defaultLocation.length()));
                if (labelText.length() == 0) {
                    // semi-hack for PROJECT_DIRECTORY
                    labelText.append(defaultLocation);
                }
            }
            label.setText(labelText.toString());
	    }
	    
	    private void reset() {
	    	keyed = false;
	    	keyedCount = 0; 
    	}
	    
	    @Override
	    public Component getTableCellEditorComponent(final JTable table, final Object value, boolean isSelected, final int row, final int column) {
	        this.table = table;
	    	StringBuilder labelText = new StringBuilder();
	    	ActionListener listener = null;
    		if (value instanceof File) {
    		    isMulti = false;
    			if (((File) value).isDirectory()) {
    			    String pathStr = ext.verifyDirFormat(((File) value).getPath());
                    if (pathStr.startsWith(defaultLocation)) {
                        pathStr = pathStr.substring(defaultLocation.length());
                        if (pathStr.trim().length() == 0) {
                            // semi-hack for PROJECT_DIRECTORY
                            pathStr = defaultLocation;
                        }
                    }
    				labelText.append(pathStr);
    			} else {
    			    String pathStr = ext.replaceAllWith(((File) value).getPath(), "\\", "/");
                    if (pathStr.startsWith(defaultLocation)) {
                        pathStr = pathStr.substring(defaultLocation.length());
                    }
    				labelText.append(pathStr);
    			}
    			listener = new ActionListener() {
    				@Override
    				public void actionPerformed(ActionEvent e) {
    					SwingUtilities.invokeLater(new Runnable() {
				            public void run() {
				            	fileChooser.setMultiSelectionEnabled(false);
				            	if (Files.exists(value.toString())) {
				            	    fileChooser.setSelectedFile((File) value);
				            	}
                                if (((File) value).isDirectory()) {
                                    fileChooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
                                } else {
                                    fileChooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
                                }
				            	if (fileChooser.showOpenDialog(buttonReplace) == JFileChooser.APPROVE_OPTION) {
				                    setValue(fileChooser.getSelectedFile());
				                } else {
				                	setValue(value);
				                }
				                fireEditingStopped();
				                if (table != null) {
				                    
				                    ((DefaultTableModel)table.getModel()).fireTableCellUpdated(row, column);
			                        ((DefaultTableModel)table.getModel()).fireTableDataChanged();
			                    }
				            }
				        });
    				}
    			};
	    	} else if (value instanceof File[]) {
	    	    isMulti = true;
	    	    boolean isDirsTemp = false;
	    		File[] files = (File[]) value;
	    		if (files.length > 0) {
	    			for (int i = 0; i < files.length; i++) {
	    				if (files[i].isDirectory()) {
	    				    isDirsTemp = true;
	    				    String pathStr = ext.verifyDirFormat(files[i].getPath());
	    				    if (pathStr.startsWith(defaultLocation)) {
                                pathStr = pathStr.substring(defaultLocation.length());
                                if (pathStr.trim().length() == 0) {
                                    // semi-hack for PROJECT_DIRECTORY
                                    pathStr = defaultLocation;
                                }
                            }
		    				labelText.append(pathStr);
		    			} else {
		    			    String pathStr = ext.replaceAllWith(files[i].getPath(), "\\", "/");
		    			    if (pathStr.startsWith(defaultLocation)) {
		    			        pathStr = pathStr.substring(defaultLocation.length());
		    			    }
			    			labelText.append(pathStr);
		    			}
	    				if (i < files.length - 1) {
	                        labelText.append(";");
	    				}
		    		}
	    		}
	    		final boolean isDirs = isDirsTemp;
    			listener = new ActionListener() {
    				@Override
    				public void actionPerformed(ActionEvent e) {
                        final JButton sourceButton = (JButton) e.getSource();
			    		SwingUtilities.invokeLater(new Runnable() {
				            public void run() {
				            	fileChooser.setMultiSelectionEnabled(true);
                                if (isDirs) {
                                    fileChooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
                                } else {
                                    fileChooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
                                }
				            	Object newValue = value;
				            	if (fileChooser.showOpenDialog(sourceButton) == JFileChooser.APPROVE_OPTION) {
//				            	    FileChooserCellEditor.this.setValue(fileChooser.getSelectedFiles());
//                                    table.setValueAt(fileChooser.getSelectedFiles(), row, column);
				            	    newValue = fileChooser.getSelectedFiles();
				                } /*else {
				                    
				                    FileChooserCellEditor.this.setValue(value);
				                }*/
				            	
				            	if (sourceButton == buttonReplace) {
				            	    setValue(newValue);
				            	} else {
				            	    HashSet<File> values = new HashSet<File>();
				            	    for (File f : (File[]) value) {
				            	        values.add(f);
				            	    }
				            	    for (File f : (File[]) newValue) {
				            	        values.add(f);
				            	    }
				            	    setValue(values.toArray(new File[values.size()]));
				            	}
//				            	table.setValueAt(newValue, row, column);
				            	
				                fireEditingStopped();
				                if (table != null) {
			                        ((DefaultTableModel)table.getModel()).fireTableCellUpdated(row, column);
			                        ((DefaultTableModel)table.getModel()).fireTableDataChanged();
			                    }
				            }
				        });
			    	}
    			};
	    	}
    		if (labelText.toString().startsWith(defaultLocation)) {
    			labelText = new StringBuilder(labelText.substring(defaultLocation.length()));
    			if (labelText.length() == 0) {
    			    labelText.append(defaultLocation);
    			}
    		}
	        label.setText(labelText.toString());
	        
	        if (isMulti) {
	            buttonAppend.setVisible(true);
	            buttonReplace.setVisible(false);
	        } else {
	            buttonAppend.setVisible(false);
	            buttonReplace.setVisible(true);
	        }
	        for (int i = buttonAppend.getActionListeners().length - 1; i >= 0; i--) {
	            buttonAppend.removeActionListener(buttonAppend.getActionListeners()[i]);
        	}
	        for (int i = buttonReplace.getActionListeners().length - 1; i >= 0; i--) {
	            buttonReplace.removeActionListener(buttonReplace.getActionListeners()[i]);
	        }
        	final ActionListener finalListener = listener;
        	buttonReplace.addActionListener(finalListener);
        	buttonAppend.addActionListener(finalListener);
	        
	        if (keyed) {
	        	keyedCount++;
	        	if (keyedCount == CLICK_COUNT_TO_START) {
	        		listener.actionPerformed(null);
	        		keyedCount = 0;	        
        		}
    		}
	        return panel;
	    }
	    public JFileChooser getFileChooser() {
	    	return fileChooser;
	    }
	}
	
	public static class SpinnerEditor extends DefaultCellEditor {
		// from http://stackoverflow.com/questions/3736993/use-jspinner-like-jtable-cell-editor
		private static final long serialVersionUID = 1L;
		JSpinner spinner;
        JSpinner.DefaultEditor editor;
        JTextField textField;
        boolean valueSet;
        JTable table;
        
        // Initializes the spinner.
        public SpinnerEditor() {
            super(new JTextField());
            spinner = new JSpinner();
            spinner.setBorder(null);
            editor = ((JSpinner.DefaultEditor)spinner.getEditor());
            editor.setBorder(null);
            textField = editor.getTextField();
            textField.setBorder(null);
            textField.setMargin(new Insets(0, 0, 0, 4));
            textField.addFocusListener(new FocusListener() {
            	public void focusGained(FocusEvent fe) {
                    SwingUtilities.invokeLater(new Runnable() {
                        public void run() {
                            if (valueSet) {
                                textField.setCaretPosition(1);
                            }
                        }
                    });
                }
                public void focusLost(FocusEvent fe) {
                }
            });
            textField.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent ae) {
                    stopCellEditing();
                    if (table != null) {
                        ((DefaultTableModel)table.getModel()).fireTableDataChanged();
                    }
                }
            });
        }

        // Prepares the spinner component and returns it.
        public Component getTableCellEditorComponent(JTable table, Object value, boolean isSelected, int row, int column) {
            this.table = table;
            if (!valueSet) {
                spinner.setValue(value);
            }
            SwingUtilities.invokeLater(new Runnable() {
                public void run() {
                    textField.requestFocus();
                }
            });
            return spinner;
        }

        public boolean isCellEditable(EventObject eo) {
            if (eo instanceof KeyEvent) {
                KeyEvent ke = (KeyEvent)eo;
                textField.setText(String.valueOf(ke.getKeyChar()));
                //textField.select(1,1);
                //textField.setCaretPosition(1);
                //textField.moveCaretPosition(1);
                valueSet = true;
            } else {
                valueSet = false;
            }
            return true;
        }

        // Returns the spinners current value.
        public Object getCellEditorValue() {
            return spinner.getValue();
        }

        public boolean stopCellEditing() {
            try {
                editor.commitEdit();
                spinner.commitEdit();
            } catch (java.text.ParseException e) {
                JOptionPane.showMessageDialog(null, "Invalid value, discarding.");
            }
            if (table != null) {
                ((DefaultTableModel)table.getModel()).fireTableDataChanged();
            }
            return super.stopCellEditing();
        }
    }
	
}

