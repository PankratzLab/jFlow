
package cnv.filesys;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.EventQueue;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;
import java.awt.event.KeyEvent;
import java.io.File;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.EventObject;
import java.util.HashMap;

import javax.swing.AbstractAction;
import javax.swing.ActionMap;
import javax.swing.DefaultCellEditor;
import javax.swing.InputMap;
import javax.swing.JButton;
import javax.swing.JCheckBox;
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
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableCellEditor;
import javax.swing.table.TableCellRenderer;

import common.Grafik;
import common.ext;

public class Configurator extends JFrame {
	private static final long serialVersionUID = 1L;

	private JPanel contentPane;
	
	private JTable table;

	private Project proj;

	/**
	 * Launch the application.
	 */
	public static void main(String[] args) {
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					boolean configure = JOptionPane.showConfirmDialog(null, "Convert property names?", "Apply conversion?", JOptionPane.YES_NO_OPTION) == JOptionPane.YES_OPTION;
					Configurator frame = new Configurator(new Project(), configure);
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
	public Configurator(Project project, boolean convertProjects) {
		setTitle("Genvisis - " + project.getNameOfProject() + " - Project Configuration");
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		setBounds(100, 100, 700, 800);
		contentPane = new JPanel();
		contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
		contentPane.setLayout(new BorderLayout(0, 0));
		setContentPane(contentPane);
		this.proj = project;
		
		JPanel panel = new JPanel();
		getContentPane().add(panel, BorderLayout.SOUTH);
		
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
				Configurator.this.save();
				Configurator.this.setVisible(false);
				Configurator.this.dispose();
			}
		});
		panel.add(btnSaveAndClose);
		
		JButton btnClose = new JButton("Close");
		btnClose.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				Configurator.this.setVisible(false);
				Configurator.this.dispose();
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
		final JPanel fileRenderer = new JPanel(new BorderLayout());
		final JLabel fileLabel = new JLabel();
		final JButton fileBtn = new JButton("...");
		fileBtn.setMargin(new Insets(1, 1, 0, 1));
		fileLabel.setBackground(Color.WHITE);
		fileRenderer.setBackground(Color.WHITE);
		fileRenderer.add(fileLabel, BorderLayout.CENTER);
		fileRenderer.add(fileBtn, BorderLayout.EAST);		
		final FileChooserCellEditor fileEditor = new FileChooserCellEditor(fileBtn, false);
		final DefaultCellEditor stringEditor = new DefaultCellEditor(new JTextField()) {
			private static final long serialVersionUID = 1L;
			@Override
			public Object getCellEditorValue() {
				return ((String) super.getCellEditorValue()).trim();
			}
		};

		final DefaultTableCellRenderer renderer = new DefaultTableCellRenderer() {
			private static final long serialVersionUID = 1L;

			@Override
			public Component getTableCellRendererComponent(final JTable table, Object value, boolean isSelected, boolean hasFocus, final int row, final int column) {
				if (value instanceof Boolean) {
					rendererChkBx.setSelected(((Boolean) value).booleanValue());
					return rendererChkBx;
				} else if (value instanceof File) {
					fileLabel.setText("./"  +((File) value).getName());
					return fileRenderer;
				} else if (value instanceof File[]) {
					StringBuilder sb = new StringBuilder();
					if (((File[])value).length > 0) {
						sb.append(((File[])value)[0].getAbsolutePath());
						for (int i = 1; i < ((File[])value).length; i++) {
							sb.append(";").append(((File[])value)[i].getName());
						}
					}
					fileLabel.setText(sb.toString());
					return fileRenderer;
				} else if (value instanceof Number) {
					String propKey = (String) table.getModel().getValueAt(row, 0);
					Object propVal = table.getModel().getValueAt(row, 1);
					setByPropertyKey(rendererSpinner, propKey, propVal);
					return rendererSpinner;
				}
				Component superComp = super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, column); 
				return superComp;
			}
		};
		
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
				
				if (propVal instanceof File || propVal instanceof File[]) {
					if (getEditRow() != row) {
						setEditRow(row);
						fileEditor.reset();
					}
					fileEditor.setValue(this.getValueAt(row, column));
					return fileEditor;
				} else if (propVal instanceof Boolean) {
					 return boolEditor;
				} else if (propVal instanceof Number) {
					setByPropertyKey(numberEditor.spinner, propKey, propVal);
					return numberEditor;
				} else if (propVal instanceof String) {
					return stringEditor;
					
				}
				
				TableCellEditor editor = super.getCellEditor(row, column); 
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
				return superComp;
			}
		};
		
		DefaultTableModel model = new DefaultTableModel(new String[]{"Property Name", "Property Value"}, 0) {
			private static final long serialVersionUID = 1L;
			@Override
			public boolean isCellEditable(int row, int column) { return column != 0 && super.isCellEditable(row, column); }
		};
		
		ArrayList<Object> propKeys = new ArrayList<Object>();
		Enumeration<?> enumer = proj.propertyNames();
		while(enumer.hasMoreElements()) {
			Object key = enumer.nextElement();
			propKeys.add(key);
		}
		for (int keyIndex = propKeys.size() - 1; keyIndex >= 0; keyIndex--) {
			Object origKey = propKeys.get(keyIndex);
			String key = origKey.toString();
			if (convertProjects && CONVERSION_MAP.containsKey(key)) {
				proj.put(CONVERSION_MAP.get(key), proj.remove(origKey));
				key = CONVERSION_MAP.get(key);
			}
			
			Object[] values = parseProperty(key);
			model.addRow(values);
		}
		
		table.setModel(model);
		table.getColumnModel().getColumn(0).setPreferredWidth(200);
		table.getColumnModel().getColumn(1).setPreferredWidth(200);
		table.setFillsViewportHeight(true);
		table.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
		table.setRowHeight(24);
		table.setRowMargin(2);
		
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
	}
	
	
	private void save() {
		int rowCount = table.getRowCount();
		for (int i = 0; i < rowCount; i++) {
			String key = (String) table.getValueAt(i, 0);
			Object rawValue = table.getValueAt(i, 1);
			String value = "";
			if (rawValue instanceof File[]) {
				File[] set = (File[]) rawValue;
				if (set.length > 0) {
					value = set[0].toString();
					for (int k = 1; k < set.length; k++) {
						value = ";" + set[k];
					}
				}
			} else {
				value = rawValue.toString();
			}
			proj.setProperty(key, value);
		}
		proj.saveProperties();
	}

	private void setByPropertyKey(JSpinner spinner, String propKey, Object propVal) {
		String[] parts = propKey.split("_");
		String minBnd = parts[parts.length - 3];
		String maxBnd = parts[parts.length - 2];
		String type = parts[parts.length - 1];
		
		if ("D".equals(type.toUpperCase())) {
			double value = "".equals(propVal) ? 0.0 : Double.valueOf(propVal.toString()).doubleValue();
			double minMult = 1.0;
			if (minBnd.toLowerCase().startsWith("n")) {
				minMult = -1.0;
				minBnd = minBnd.substring(1);
			}
			String[] befAf = minBnd.split("d");
			double minimum = minMult * Double.valueOf(befAf[0] + (befAf.length > 1 ? "." + befAf[1] : "")).doubleValue();
			double maxMult = 1.0;
			if (maxBnd.toLowerCase().startsWith("n")) {
				maxMult = -1.0;
				maxBnd = maxBnd.substring(1);
			}
			befAf = maxBnd.split("d");
			double maximum = maxMult * Double.valueOf(befAf[0] + (befAf.length > 1 ? "." + befAf[1] : "")).doubleValue();
			double stepSize = 1.0;
			for (int i = 0; i < ext.getNumSigFig(value); i++) {
				stepSize /= 10.0;
			}
			spinner.setModel(new SpinnerNumberModel(value, minimum, maximum, stepSize));
		} else if ("I".equals(type.toUpperCase())) {
			int value = "".equals(propVal) ? 0 : Integer.valueOf(propVal.toString()).intValue();
			int minMult = 1;
			if (minBnd.toLowerCase().startsWith("n")) {
				minMult = -1;
				minBnd = minBnd.substring(1);
			}
			int minimum = minMult * Integer.valueOf(minBnd).intValue();
			int maxMult = 1;
			if (maxBnd.toLowerCase().startsWith("n")) {
				maxMult = -1;
				maxBnd = maxBnd.substring(1);
			}
			int maximum = maxMult * Integer.valueOf(maxBnd).intValue();
			int stepSize = 1;
			spinner.setModel(new SpinnerNumberModel(value, minimum, maximum, stepSize));
		}
	}
	
	private Object[] parseProperty(String prop) {
		Object[] keyVal = new Object[2];
		keyVal[0] = prop;
		
		String[] parts = prop.split("_");
		String sig = parts[parts.length - 1];
		
		int index = -1;
		for (int i = 0; i < SIGNIFIERS.length; i++) {
			if (SIGNIFIERS[i].equals(sig)) {
				index = i;
				break;
			}
		}
		
		switch(index) {
			case 0: // DIRECTORY
				keyVal[1] = new File(proj.getDir(prop, false, false));
				break;
			case 1: // FILENAME
				keyVal[1] = new File(proj.getFilename(prop, false, false));
				break;
			case 2: // FILENAMEs
				String[] files = proj.getFilenames(prop, true); 
				keyVal[1] = new File[files.length];
				for (int i = 0; i < files.length; i++) {
					((File[])keyVal[1])[i] = new File(files[i]);
				}
				break;
			case 3: // LIST
				keyVal[1] = proj.getProperty(prop);
				break;
			case 4: // YESNO
				keyVal[1] = proj.getBoolean(prop);
				break;
			case 5: // I
				keyVal[1] = proj.getInt(prop);
				break;
			case 6: // D
				keyVal[1] = proj.getDouble(prop);
				break;
			default: // No preset, default to text field
				keyVal[1] = proj.getProperty(prop);
				break;
		}
		
		return keyVal;	
	}
	
	
	public static final String[] SIGNIFIERS = new String[] {
		"DIRECTORY",
		"FILENAME",
		"FILENAMES",
		"LIST",
		"YESNO",
		"I",
		"D"
	};
	
	HashMap<String, String> CONVERSION_MAP = new HashMap<String, String>() {
		private static final long serialVersionUID = 1L;
		{
			put("DISPLAY_QUANTILES", "DISPLAY_QUANTILES_YESNO");
			put("DISPLAY_ROTATED_QQ", "DISPLAY_ROTATED_QQ_YESNO");
			put("DISPLAY_STANDARD_QQ", "DISPLAY_STANDARD_QQ_YESNO");
			put("FID_ALIAS", "FID_ALIAS_LIST");
			put("GC_THRESHOLD", "GC_THRESHOLD_0d0_1d0_D");
			put("IID_ALIAS", "IID_ALIAS_LIST");
			put("INTENSITY_PC_NUM_COMPONENTS", "INTENSITY_PC_NUM_COMPONENTS_0_10000_I");
			put("JAR_STATUS", "JAR_STATUS_YESNO");
			put("LOG_LEVEL", "LOG_LEVEL_n1_12_I");
			put("LONG_FORMAT", "LONG_FORMAT_YESNO");
			put("LRRSD_CUTOFF", "LRRSD_CUTOFF_0d0_1d0_D");
			put("MAX_MARKERS_LOADED_PER_CYCLE", "MAX_MARKERS_LOADED_PER_CYCLE_1_10000_I");
			put("MAX_MEMORY_USED_TO_LOAD_MARKER_DATA", "MAX_MEMORY_USED_TO_LOAD_MARKER_DATA_8_65536_I");
			put("NUM_THREADS", "NUM_THREADS_1_99_I");
			put("PARSE_AT_AT_SYMBOL", "PARSE_AT_AT_SYMBOL_YESNO");
			put("QQ_MAX_NEG_LOG10_PVALUE", "QQ_MAX_NEG_LOG10_PVALUE_1_10000_I");
			put("SAMPLE_ALIAS", "SAMPLE_ALIAS_LIST");
			put("TWOD_LOADED_VARIABLES", "TWOD_LOADED_VARIABLES_LIST");
			put("WINDOW_AROUND_SNP_TO_OPEN_IN_TRAILER", "WINDOW_AROUND_SNP_TO_OPEN_IN_TRAILER_100_1000000_I");
		}	
	};
	
	String[] properties  = new String[]{
			"AB_LOOKUP_FILENAME",
			"ANNOTATION_FILENAME",
			"BACKUP_DIRECTORY",
			"CHIMERA_CENTROIDS_FILENAME",
			"CLUSTER_FILTER_COLLECTION_FILENAME",
			"CNV_FILENAMES",
			"COMMON_CNP_FILENAME",
			"CUSTOM_CENTROIDS_FILENAME",
			"CUSTOM_COLOR_SCHEME_FILENAME",
			"DATA_DIRECTORY",
			"DEMO_DIRECTORY",
			"DISPLAY_MARKERS_FILENAME",
			"DISPLAY_QUANTILES_YESNO",
			"DISPLAY_ROTATED_QQ_YESNO",
			"DISPLAY_STANDARD_QQ_YESNO",
			"FID_ALIAS_LIST",
			"FILTERED_MARKERS_FILENAME",
			"FOREST_PLOT_FILENAMES",
			"GC_MODEL_FILENAME",
			"GC_THRESHOLD_0d0_1d0_D",
			"GENETRACK_FILENAME",
			"GENOTYPE_CENTROIDS_FILENAME",
			"ID_HEADER",
			"IID_ALIAS_LIST",
			"INDIVIDUAL_CNV_LIST_FILENAMES",
			"INTENSITY_PC_FILENAME",
			"INTENSITY_PC_NUM_COMPONENTS_0_10000_I",
			"JAR_STATUS_YESNO",
			"LOG_LEVEL_n1_12_I",
			"LONG_FORMAT_YESNO",
			"LRRSD_CUTOFF_0d0_1d0_D",
			"MARKER_COMBINED_CRITERIA_FILENAME",
			"MARKER_DATA_DIRECTORY",
			"MARKER_EXCLUSION_CRITERIA_FILENAME",
			"MARKER_METRICS_FILENAME",
			"MARKER_POSITION_FILENAME",
			"MARKER_REVIEW_CRITERIA_FILENAME",
			"MARKERLOOKUP_FILENAME",
			"MARKERSET_FILENAME",
			"MAX_MARKERS_LOADED_PER_CYCLE_1_10000_I",
			"MAX_MEMORY_USED_TO_LOAD_MARKER_DATA_8_65536_I",
			"MOSAIC_ARMS_FILENAME",
			"MOSAIC_COLOR_CODES_FILENAME",
			"MOSAIC_RESULTS_FILENAME",
			"NUM_THREADS_1_99_I",
			"ORIGINAL_CENTROIDS_FILENAME",
			"PARSE_AT_AT_SYMBOL_YESNO",
			"PEDIGREE_FILENAME",
			"PENNCNV_DATA_DIRECTORY",
			"PENNCNV_EXECUTABLE_DIRECTORY",
			"PENNCNV_GZIP_YESNO",
			"PENNCNV_RESULTS_DIRECTORY",
//			"PROJECT_DIRECTORY",
//			"PROJECT_NAME",
//			"PROJECT_PROPERTIES_FILENAME",
			"QQ_FILENAMES",
			"QQ_MAX_NEG_LOG10_PVALUE_1_10000_I",
			"REGION_LIST_FILENAMES",
			"REPORTED_CNP_FILENAME",
			"RESULTS_DIRECTORY",
			"SAMPLE_ALIAS_LIST",
			"SAMPLE_DATA_FILENAME",
			"SAMPLE_DIRECTORY",
			"SAMPLE_QC_FILENAME",
			"SAMPLE_SUBSET_FILENAME",
			"SAMPLELIST_FILENAME",
			"SEX_CENTROIDS_FILENAMES",
			"SEXCHECK_RESULTS_FILENAME",
			"SHIFT_SEX_CHR_COLORS_YESNO",
			"SOURCE_DIRECTORY",
			"SOURCE_FILE_DELIMITER",
			"SOURCE_FILENAME_EXTENSION",
			"STRATIFICATION_RESULTS_FILENAMES",
			"TARGET_MARKERS_FILENAME",
			"TWOD_LOADED_FILENAMES",
			"TWOD_LOADED_VARIABLES_LIST",
			"UNREPORTED_CNP_FILENAME",
			"WINDOW_AROUND_SNP_TO_OPEN_IN_TRAILER_100_1000000_I",
		};
		
	class FileChooserCellEditor extends DefaultCellEditor implements TableCellEditor {
		// from http://stackoverflow.com/questions/15644319/jfilechooser-within-a-jtable
		private static final long serialVersionUID = 1L;
		/** Number of clicks to start editing */
	    private static final int CLICK_COUNT_TO_START = 1;
	    /** meta panel */
	    private JPanel panel;
	    /** static display */
	    private JLabel label;	
	    /** Editor component */
	    private JButton button;
	    /** File chooser */
	    private JFileChooser fileChooser;
	    /** Selected file(s) */
	    private Object value;
	    private JTextField textField;
	    volatile boolean keyed = false;
	    volatile int keyedCount = 0;
	    /**
	     * Constructor.
	     */
	    public FileChooserCellEditor(final JButton button, boolean editable) {
	        super(new JTextField());
	        textField = (JTextField) this.editorComponent;
	        textField.setEditable(editable);
	        setClickCountToStart(CLICK_COUNT_TO_START);
	        setBackground(Color.WHITE);
	        
	        panel = new JPanel(new BorderLayout());
	        panel.setBackground(Color.WHITE);
	        
	        // Using a JButton as the editor component
	        label = new JLabel();
	        label.setBackground(Color.WHITE);
	        
	        this.button = button;
//	        button = new JButton("...");
//	        button.setFont(button.getFont().deriveFont(Font.PLAIN));
//	        button.setBorder(null);
//	        button.setMargin(new Insets(0, 0, 0, 0));
	        
	        panel.add(label, BorderLayout.CENTER);
	        panel.add(button, BorderLayout.EAST);

	        panel.setBackground(Color.WHITE);
	        // Dialog which will do the actual editing
	        fileChooser = new JFileChooser();
	    }

	    @Override
	    public Object getCellEditorValue() {
	        return value;
	    }
	    
	    private void setValue(Object val) {
	    	this.value = val;
	    }
	    
	    private void reset() {
	    	keyed = false;
	    	keyedCount = 0; 
    	}
	    
	    @Override
	    public Component getTableCellEditorComponent(final JTable table, final Object value, boolean isSelected, final int row, final int column) {
	    	StringBuilder labelText = new StringBuilder();
	    	ActionListener listener = null;
    		if (value instanceof File) {
    			labelText.append(((File) value).getAbsolutePath());
    			listener = new ActionListener() {
    				@Override
    				public void actionPerformed(ActionEvent e) {
    					SwingUtilities.invokeLater(new Runnable() {
				            public void run() {
				            	fileChooser.setMultiSelectionEnabled(false);
				            	fileChooser.setSelectedFile((File) value);
				            	if (fileChooser.showOpenDialog(button) == JFileChooser.APPROVE_OPTION) {
				                    setValue(fileChooser.getSelectedFile());
				                } else {
				                	setValue(value);
				                }
				                fireEditingStopped();
				            }
				        });
    				}
    			};
	    	} else if (value instanceof File[]) {
	    		File[] files = (File[]) value;
	    		if (files.length > 0) {
	    			labelText.append(files[0].getAbsolutePath());
	    			for (int i = 1; i < files.length; i++) {
		    			labelText.append(";").append(files[i].getAbsolutePath());
		    		}
	    		}

    			listener = new ActionListener() {
    				@Override
    				public void actionPerformed(ActionEvent e) {
			    		SwingUtilities.invokeLater(new Runnable() {
				            public void run() {
				            	fileChooser.setMultiSelectionEnabled(true);
				            	fileChooser.setSelectedFiles((File[]) value);
				            	if (fileChooser.showOpenDialog(button) == JFileChooser.APPROVE_OPTION) {
				                    setValue(fileChooser.getSelectedFiles());
				                } else {
				                	setValue(value);
				                }
				                fireEditingStopped();
				            }
				        });
			    	}
    			};
	    	}
	        label.setText(labelText.toString());
	        for (int i = button.getActionListeners().length - 1; i >= 0; i--) {
	        	button.removeActionListener(button.getActionListeners()[i]);
        	}
        	final ActionListener finalListener = listener;
	        button.addActionListener(finalListener);
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
                }
            });
        }

        // Prepares the spinner component and returns it.
        public Component getTableCellEditorComponent(JTable table, Object value, boolean isSelected, int row, int column) {
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
            return super.stopCellEditing();
        }
    }
	
}

