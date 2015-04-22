
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
import java.util.EventObject;

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
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableCellEditor;
import javax.swing.table.TableCellRenderer;

import cnv.LaunchProperties;
import cnv.filesys.Project.DoubleProperty;
import cnv.filesys.Project.FileProperty;
import cnv.filesys.Project.IntegerProperty;
import cnv.filesys.Project.MultiFileProperty;
import cnv.filesys.Project.Property;
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
					Configurator frame = new Configurator(new Project());
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
	public Configurator(Project project) {
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
		final FileChooserCellEditor fileEditor = new FileChooserCellEditor(fileBtn, false, proj.getProperty(proj.PROJECT_DIRECTORY));
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

		final DefaultTableCellRenderer renderer = new DefaultTableCellRenderer() {
			private static final long serialVersionUID = 1L;

			@Override
			public Component getTableCellRendererComponent(final JTable table, Object value, boolean isSelected, boolean hasFocus, final int row, final int column) {
				String projDir = proj.getProperty(proj.PROJECT_DIRECTORY);
				
				String desc = proj.getProperty((String) table.getValueAt(row, 0)).getDescription();
				Component returnComp;
				
				if (value instanceof Boolean) {
					rendererChkBx.setSelected(((Boolean) value).booleanValue());
					returnComp = rendererChkBx;
				} else if (value instanceof File) {
					String txt = ((File) value).getPath();
					if (((File) value).isDirectory()) {
						txt = ext.verifyDirFormat(txt);
					} else {
						txt = ext.replaceAllWith(txt, "\\", "/");
					}
					if (!txt.startsWith("./") && txt.indexOf(":") == -1) {
						txt = "./" + txt;
					}
					if (txt.startsWith(projDir)) {
						txt = txt.substring(projDir.length());
					}
					fileLabel.setText(txt);
					returnComp = fileRenderer;
				} else if (value instanceof File[]) {
					StringBuilder sb = new StringBuilder();
					if (((File[])value).length > 0) {
						String txt = ((File[])value)[0].getPath();
						if (((File[]) value)[0].isDirectory()) {
							txt = ext.verifyDirFormat(txt);
						} else {
							txt = ext.replaceAllWith(txt, "\\", "/");
						}
						if (txt.startsWith(projDir)) {
							txt = txt.substring(projDir.length());
						}
						
						sb.append(txt);
						for (int i = 1; i < ((File[])value).length; i++) {
							txt = ((File[])value)[i].getPath();
							if (((File[]) value)[i].isDirectory()) {
								txt = ext.verifyDirFormat(txt);
							} else {
								txt = ext.replaceAllWith(txt, "\\", "/");
							}
							if (txt.startsWith(projDir)) {
								txt = txt.substring(projDir.length());
							}
							sb.append(";").append(txt);
						}
					}
					fileLabel.setText(sb.toString());
					returnComp = fileRenderer;
				} else if (value instanceof Number) {
					String propKey = (String) table.getModel().getValueAt(row, 0);
					Object propVal = table.getModel().getValueAt(row, 1);
					setByPropertyKey(rendererSpinner, proj, propKey, propVal);
					returnComp = rendererSpinner;
				} else {
					returnComp = super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, column);
				}
				if (!"".equals(desc)) {
					((JComponent) returnComp).setToolTipText(desc);
				} else {
					((JComponent) returnComp).setToolTipText(null);
				}
				return returnComp;
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
					setByPropertyKey(numberEditor.spinner, proj, propKey, propVal);
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
		
		String[] keys = proj.getPropertyKeys();
		for (String key : keys) {
			Object[] values = parseProperty(proj, key);
			model.addRow(values);
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
	}
	
	
	private void save() {
		String projectsDir = new LaunchProperties(LaunchProperties.DEFAULT_PROPERTIES_FILE).getProperty(LaunchProperties.PROJECTS_DIR);
		String currProjDir = proj.getProperty(proj.PROJECT_DIRECTORY);
		int rowCount = table.getRowCount();
		for (int i = 0; i < rowCount; i++) {
			String key = (String) table.getValueAt(i, 0);
			Object rawValue = table.getValueAt(i, 1);
//			System.out.println(rawValue.getClass().getName());
			String value = "";
			if (rawValue instanceof File[]) {
				File[] set = (File[]) rawValue;
				if (set.length > 0) {
					value = set[0].getPath();
					value = set[0].isDirectory() ? ext.verifyDirFormat(value) : ext.replaceAllWith(value, "\\", "/");
					if (value.startsWith(projectsDir)) {
						value = value.substring(projectsDir.length());
					} else if (value.startsWith(currProjDir)) {
						value = "./" + value.substring(currProjDir.length());
					}
					for (int k = 1; k < set.length; k++) {
						String fNm = set[k].getPath();
						fNm = set[k].isDirectory() ? ext.verifyDirFormat(fNm) : ext.replaceAllWith(fNm, "\\", "/");
						if (fNm.startsWith(projectsDir)) {
							fNm = fNm.substring(projectsDir.length());
						} else if (fNm.startsWith(currProjDir)) {
							fNm = "./" + fNm.substring(currProjDir.length());
						}
						value += ";" + fNm;
					}
				}
			} else if (rawValue instanceof File) { 
				File set = (File) rawValue;
				value = set.getPath();
				value = set.isDirectory() ? ext.verifyDirFormat(value) : ext.replaceAllWith(value, "\\", "/");
				if (!key.equals(proj.SOURCE_DIRECTORY.getName()) && !key.equals(proj.PROJECT_DIRECTORY.getName())) {
					if (value.startsWith(projectsDir)) {
						value = value.substring(projectsDir.length());
					} else if (value.startsWith(currProjDir)) {
						value = "./" + value.substring(currProjDir.length());
					}
				}
			} else if (rawValue instanceof String[]) {
				value = ((String[])rawValue)[0];
				for (int k = 1; k < ((String[])rawValue).length; k++) {
					value += ";" + ((String[])rawValue)[k];
				}
			} else {
				value = rawValue.toString();
			}
			proj.setProperty(key, value);
		}
		proj.saveProperties();
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
		} else if (prop instanceof MultiFileProperty) {
			File[] fs = new File[((MultiFileProperty)prop).getValue().length];
			for (int i = 0; i < fs.length; i++) {
				fs[i] = new File(((MultiFileProperty)prop).getValue()[i]);
			}
			keyVal[1] = fs;
		}
		
		return keyVal;	
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
	    private JButton button;
	    /** File chooser */
	    private JFileChooser fileChooser;
	    /** Selected file(s) */
	    private Object value;
	    private JTextField textField;
	    volatile boolean keyed = false;
	    volatile int keyedCount = 0;
	    String defaultLocation;
	    /**
	     * Constructor.
	     */
	    public FileChooserCellEditor(final JButton button, boolean editable, String defaultLocation) {
	        super(new JTextField());
	        this.defaultLocation = defaultLocation;
	        textField = (JTextField) this.editorComponent;
	        textField.setEditable(editable);
	        setClickCountToStart(CLICK_COUNT_TO_START);
	        setBackground(Color.WHITE);
	        
	        panel = new JPanel(new BorderLayout());
	        panel.setBackground(Color.WHITE);
	        
	        // Using a JButton as the editor component
	        label = new JTextField();
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
	        fileChooser = new JFileChooser(defaultLocation);
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
    			if (((File) value).isDirectory()) {
    				labelText.append(ext.verifyDirFormat(((File) value).getPath()));
    			} else {
    				labelText.append(ext.replaceAllWith(((File) value).getPath(), "\\", "/"));
    			}
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
	    			if (files[0].isDirectory()) {
	    				labelText.append(ext.verifyDirFormat(files[0].getPath()));
	    			} else {
		    			labelText.append(ext.replaceAllWith(files[0].getPath(), "\\", "/"));
	    			}
	    			for (int i = 1; i < files.length; i++) {
	    				labelText.append(";");
	    				if (files[i].isDirectory()) {
		    				labelText.append(ext.verifyDirFormat(files[i].getPath()));
		    			} else {
			    			labelText.append(ext.replaceAllWith(files[i].getPath(), "\\", "/"));
		    			}
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

