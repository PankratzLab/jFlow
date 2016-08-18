package org.genvisis.cnv.plots;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.InputEvent;
import java.awt.event.KeyEvent;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.LineNumberReader;
import java.io.PrintWriter;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.Vector;

import javax.swing.AbstractAction;
import javax.swing.ActionMap;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.DefaultComboBoxModel;
import javax.swing.InputMap;
import javax.swing.JButton;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLayeredPane;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JRadioButtonMenuItem;
import javax.swing.JSeparator;
import javax.swing.JTextField;
import javax.swing.KeyStroke;
import javax.swing.SwingUtilities;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.common.Aliases;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Grafik;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

public class ForestPlotFrame extends JFrame implements WindowListener {
	private static final long serialVersionUID = 1L;
	
	public static final Color BACKGROUND_COLOR = Color.WHITE;
	public static final String ADD_DATA_FILE = "Add Data File";
	public static final String REMOVE_DATA_FILE = "Remove Data File";
//	private static final String ALT_UP = "ALT UP";
//	private static final String ALT_DOWN = "ALT DOWN";
//	private static final String ALT_LEFT = "ALT LEFT";
//	private static final String ALT_RIGHT = "ALT RIGHT";

	private static final String FIRST = "First";
	private static final String PREVIOUS = "Previous";
	private static final String NEXT = "Next";
	private static final String LAST = "Last";
	private static final String MARKER_LIST_NEW_FILE = "Add New File(s)...";
	private static final String MARKER_LIST_PLACEHOLDER = "Select Input File...";
//	private LinkedHashSet<ForestInput> data = new LinkedHashSet<ForestInput>();
//	private HashMap<ForestInput, MetaStudy> dataToMetaMap = new HashMap<ForestInput, MetaStudy>();

	

//	private JCheckBox chkSortStudies;
	private JLayeredPane layeredPane;
	private JButton first, previous, next, last;
//	private JButton btnScreen, btnScreenAll;
	private JTextField navigationField;
	private JComboBox markerFileList;
//	int curMarkerIndex;
    

	private HashMap<String, String> markerFileNameLoc = new HashMap<String, String>();
	

	private Project proj;
	private ForestPlot forestPlot;
	
	
	public ForestPlotFrame(Project proj) {
		super("Genvisis - Forest Plot - " + proj.PROJECT_NAME.getValue());
		setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
		forestPlot = new ForestPlot(proj);
		this.proj = proj;
		setup();
		setBounds(20, 20, 1000, 600);
		pack();
		setVisible(true);
		this.addWindowListener(this);
	}
	
	public ForestPlotFrame(String markerFile, Logger log) {
		super("Genvisis - Forest Plot");
		setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
		forestPlot = new ForestPlot(markerFile, log);
		setup();
		setBounds(20, 20, 1000, 600);
		pack();
		setVisible(true);
		this.addWindowListener(this);
	}

	private void setup() {
		setLayout(new BorderLayout());

		layeredPane = new JLayeredPane();
		layeredPane.setLayout(new BorderLayout());

		layeredPane.add(forestPlot.getForestPanel());
		layeredPane.setPreferredSize(new Dimension(1000, 600));
		
		forestPlot.getForestPanel().add(createControlPanel(), BorderLayout.NORTH);
		
		Dimension minimumSize = new Dimension(100, 50);
		layeredPane.setMinimumSize(minimumSize);

		add(layeredPane, BorderLayout.CENTER);
		inputMapAndActionMap();

		forestPlot.getForestPanel().grabFocus();

		progressBar = new JProgressBar();
		progressBar.setPreferredSize(new Dimension(progressBar.getPreferredSize().width, 22));
		progressBar.setMinimumSize(new Dimension(progressBar.getPreferredSize().width, 22));
		
		JButton btnCancelProg = new JButton(Grafik.getImageIcon("images/delete2sm.png"));
		btnCancelProg.setMinimumSize(new Dimension(30, 20));
		btnCancelProg.setPreferredSize(new Dimension(30, 20));
		btnCancelProg.setMaximumSize(new Dimension(30, 20));
		btnCancelProg.setAlignmentX(Component.RIGHT_ALIGNMENT);
		btnCancelProg.addActionListener(new AbstractAction() {
			private static final long serialVersionUID = 1L;
			@Override
			public void actionPerformed(ActionEvent e) {
				forestPlot.getLoadingThread().interrupt();
			}
		});
		
		progressBar.setLayout(new BorderLayout());
		progressBar.add(btnCancelProg, BorderLayout.EAST);
		
		add(progressBar, BorderLayout.SOUTH);
		
		setJMenuBar(createMenuBar());
		
		setVisible(true);
		
		loadMarkerFile();
		// generateShortcutMenus();
	}
	
	private JMenuBar createMenuBar() {
		JMenuBar bar = new JMenuBar();
		
		JMenu actions = new JMenu("Actions", true);
		actions.setMnemonic(KeyEvent.VK_A);
		
		final JCheckBoxMenuItem oddsRatioButton = new JCheckBoxMenuItem();
		AbstractAction oddsRatioAction = new AbstractAction() {
            private static final long serialVersionUID = 1L;
            @Override
            public void actionPerformed(ActionEvent e) {
                forestPlot.setOddsRatioDisplay(oddsRatioButton.isSelected());
                forestPlot.getForestPanel().paintAgain();                
            }
        };
		oddsRatioButton.setAction(oddsRatioAction);
		oddsRatioButton.setText("Display as Odds Ratios");
		oddsRatioButton.setMnemonic(KeyEvent.VK_O);
		actions.add(oddsRatioButton);
		
		actions.add(new JSeparator());
        
		JMenuItem createSortFile = new JMenuItem();
		createSortFile.setAction(new AbstractAction() {
            private static final long serialVersionUID = 1L;
            @Override
            public void actionPerformed(ActionEvent e) {
                generateStudyNameFile();
            }
		});
		createSortFile.setText("Create Sort File");
		actions.add(createSortFile);
		
		final JRadioButtonMenuItem noSortButton = new JRadioButtonMenuItem();
		AbstractAction noSortAction = new AbstractAction() {
            private static final long serialVersionUID = 1L;
            @Override
            public void actionPerformed(ActionEvent e) {
                forestPlot.setSortedDisplay(false);
                forestPlot.getForestPanel().paintAgain();                
            }
        };
        noSortButton.setAction(noSortAction);
        noSortButton.setText("No Sorting");
        noSortButton.setMnemonic(KeyEvent.VK_N);
        actions.add(noSortButton);
		
		final JRadioButtonMenuItem sortStudies = new JRadioButtonMenuItem();
		AbstractAction sortAction = new AbstractAction() {
			private static final long serialVersionUID = 1L;
			@Override
			public void actionPerformed(ActionEvent arg0) {
				forestPlot.setSortedDisplay(true);
				forestPlot.getForestPanel().paintAgain();
			}
		};
		sortStudies.setAction(sortAction);
		sortStudies.setText("Sort Studies");
		sortStudies.setMnemonic(KeyEvent.VK_S);
		actions.add(sortStudies);
		
		final JRadioButtonMenuItem sortFileStudies = new JRadioButtonMenuItem();
		AbstractAction sortFileAction = new AbstractAction() {
            private static final long serialVersionUID = 1L;
            @Override
            public void actionPerformed(ActionEvent e) {
                
                JFileChooser jfc = new JFileChooser(new File(proj == null ? ext.parseDirectoryOfFile(forestPlot.getMarkerFileName()) : proj.PROJECT_DIRECTORY.getValue()));
                int returnVal = jfc.showOpenDialog(null);
                
                if (returnVal == JFileChooser.APPROVE_OPTION) {
                    forestPlot.loadOrderFile(jfc.getSelectedFile().getAbsolutePath(), true);
                    forestPlot.getForestPanel().paintAgain();
                } else {
                    noSortButton.setSelected(true);
                    forestPlot.setSortedDisplay(false);
                    forestPlot.getForestPanel().paintAgain();
                }
            }
        };
    	sortFileStudies.setAction(sortFileAction);
    	sortFileStudies.setText("Sort Studies with File");
    	sortFileStudies.setMnemonic(KeyEvent.VK_S);
    	actions.add(sortFileStudies);
        
        ButtonGroup sortGroup = new ButtonGroup();
        sortGroup.add(noSortButton);
        sortGroup.add(sortStudies);
        sortGroup.add(sortFileStudies);
        
        noSortButton.setSelected(true);
		
		AbstractAction screenAction1 = new AbstractAction() {
			private static final long serialVersionUID = 1L;
			@Override
			public void actionPerformed(ActionEvent e) {
				if (!forestPlot.isLoadingFile() && (forestPlot.getMarkerFileName() != null && !"".equals(forestPlot.getMarkerFileName()))) {
					ForestPlotFrame.this.screenCap();
				}
			}
		};
		AbstractAction screenAction2 = new AbstractAction() {
			private static final long serialVersionUID = 1L;
			@Override
			public void actionPerformed(ActionEvent e) {
				if (!forestPlot.isLoadingFile() && (forestPlot.getMarkerFileName() != null && !"".equals(forestPlot.getMarkerFileName()))) {
					ForestPlotFrame.this.screenCapAll();
				}
			}
		};

		actions.add(new JSeparator());
		
		JMenuItem screenCap1 = new JMenuItem();
		screenCap1.setAction(screenAction1);
		screenCap1.setText("Screen Capture");
		screenCap1.setMnemonic(KeyEvent.VK_C);
		actions.add(screenCap1);
		
		JMenuItem screenCap2 = new JMenuItem();
		screenCap2.setAction(screenAction2);
		screenCap2.setText("Screen Capture All");
		screenCap2.setMnemonic(KeyEvent.VK_A);
		actions.add(screenCap2);
		
		bar.add(actions);
		
		return bar;
	}

    private JPanel createControlPanel() {
		first = new JButton(Grafik.getImageIcon("images/firstLast/First.gif"));
		first.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dFirst.gif"));
		first.addActionListener(navFirstAction);
		first.setActionCommand(FIRST);
		first.setPreferredSize(new Dimension(20, 20));
		previous = new JButton(Grafik.getImageIcon("images/firstLast/Left.gif"));
		previous.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dLeft.gif"));
		previous.addActionListener(navPrevAction);
		previous.setActionCommand(PREVIOUS);
		previous.setPreferredSize(new Dimension(20, 20));
		navigationField = new JTextField("", 8);
		navigationField.setHorizontalAlignment(JTextField.CENTER);
		navigationField.setFont(new Font("Arial", 0, 14));
		//navigationField.setEditable(false);
		navigationField.setBackground(BACKGROUND_COLOR);
		navigationField.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				try {
					int trav = Integer.valueOf(((JTextField)e.getSource()).getText().split("[\\s]+")[0]).intValue()-1;
					if (trav >=0 && trav < forestPlot.getDataIndices().size()) {
						forestPlot.setCurrentData(trav);
						updateForestPlot();
					}
				} catch (NumberFormatException nfe) {
					System.out.println("Please enter a valid integer which is in range");
				}
				//displayIndex((JTextField) e.getSource());
				forestPlot.getForestPanel().setPointsGeneratable(true);
				updateForestPlot();
			}
		});
	
		next = new JButton(Grafik.getImageIcon("images/firstLast/Right.gif"));
		next.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dRight.gif"));
		next.addActionListener(navNextAction);
		next.setActionCommand(NEXT);
		next.setPreferredSize(new Dimension(20, 20));
		last = new JButton(Grafik.getImageIcon("images/firstLast/Last.gif"));
		last.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dLast.gif"));
		last.addActionListener(navLastAction);
		last.setActionCommand(LAST);
		last.setPreferredSize(new Dimension(20, 20));
		
		JPanel subNavPanel = new JPanel();
		subNavPanel.setBackground(BACKGROUND_COLOR);
		subNavPanel.add(first);
		subNavPanel.add(previous);
		subNavPanel.add(navigationField);
		subNavPanel.add(next);
		subNavPanel.add(last);
		
		Vector<String> items = new Vector<String>();
		items.add(MARKER_LIST_PLACEHOLDER);
		if (proj != null) {
			String[] files = proj.FOREST_PLOT_FILENAMES.getValue();
			String name;
			for (String file : files) {
				name = ext.rootOf(file);
				markerFileNameLoc.put(name, file);
				items.add(name);
			}
		}
		items.add(MARKER_LIST_NEW_FILE);
		markerFileList = new JComboBox(items);
		markerFileList.setFont(new Font("Arial", 0, 12));
		markerFileList.setMinimumSize(new Dimension(200, 20));
		markerFileList.setPreferredSize(new Dimension(200, 20));
		markerFileList.setMaximumSize(new Dimension(200, 20));
		AbstractAction markerFileSelectAction = new AbstractAction() {
			private static final long serialVersionUID = 1L;
	
			@SuppressWarnings("unchecked")
			@Override
			public void actionPerformed(ActionEvent e) {
				String shortName = (String) ((JComboBox)e.getSource()).getSelectedItem();
				if (!forestPlot.isLoadingFile() && !MARKER_LIST_NEW_FILE.equals(shortName) && !MARKER_LIST_PLACEHOLDER.equals(shortName)) {
					String file = markerFileNameLoc.get(shortName);
					if (file != null && file.equals(forestPlot.getMarkerFileName())) {
						return;
					}
					forestPlot.setMarkerFileName(file);
					loadMarkerFile();
				} else if (forestPlot.isLoadingFile() && MARKER_LIST_PLACEHOLDER.equals(shortName)) {
					// do nothing
				} else if (forestPlot.isLoadingFile() || MARKER_LIST_PLACEHOLDER.equals(shortName)) {
					// leave as currently selected marker
					if (forestPlot.getMarkerFileName() != "" && forestPlot.getMarkerFileName() != null) {
						((JComboBox)e.getSource()).setSelectedItem(ext.rootOf(forestPlot.getMarkerFileName()));
					}
					return;
				} else if (MARKER_LIST_NEW_FILE.equals(shortName)) {
					chooseNewFiles();
					if (forestPlot.getMarkerFileName() != null && !"".equals(forestPlot.getMarkerFileName())) {
						((JComboBox)e.getSource()).setSelectedItem(ext.rootOf(forestPlot.getMarkerFileName()));
					}
				}
			}
		};
		markerFileList.addActionListener(markerFileSelectAction);
		
		JPanel navigationPanel = new JPanel();
		navigationPanel.setLayout(new BoxLayout(navigationPanel, BoxLayout.Y_AXIS));
		navigationPanel.setBackground(BACKGROUND_COLOR);
		navigationPanel.add(subNavPanel);
		navigationPanel.add(markerFileList);
		navigationPanel.add(Box.createVerticalGlue());
		
//		final JPanel subOptPanel = new JPanel();
//		subOptPanel.setLayout(new BoxLayout(subOptPanel, BoxLayout.Y_AXIS));
//		subOptPanel.setBackground(BACKGROUND_COLOR);
//		subOptPanel.setAlignmentX(Component.RIGHT_ALIGNMENT);
//		subOptPanel.add(btnScreen);
//		subOptPanel.add(Box.createRigidArea(new Dimension(0,5)));
//		subOptPanel.add(btnScreenAll);
		
//		JPanel optPanel = new JPanel() {
//			private static final long serialVersionUID = 1L;
//			@Override
//			public Dimension getMaximumSize() {
//				return super.getPreferredSize();
//			}
//		};
//		optPanel.setBackground(BACKGROUND_COLOR);
//		optPanel.add(chkSortStudies);
//		optPanel.add(subOptPanel);
	
		
		JPanel descrPanel = new JPanel();
		descrPanel.setLayout(new BoxLayout(descrPanel, BoxLayout.X_AXIS));
		// All of these are needed (or at least, have some effect)
		descrPanel.add(Box.createHorizontalGlue());
		descrPanel.add(Box.createHorizontalGlue());
		descrPanel.add(Box.createHorizontalGlue());
//		descrPanel.add(Box.createHorizontalGlue());
//		descrPanel.add(Box.createHorizontalGlue());
		descrPanel.add(navigationPanel);
		descrPanel.add(Box.createHorizontalGlue());
		descrPanel.add(Box.createHorizontalGlue());
		descrPanel.add(Box.createHorizontalGlue());
//		descrPanel.add(optPanel);
		descrPanel.setBackground(BACKGROUND_COLOR);
		return descrPanel;
	}

    
	private void loadMarkerFile() {
		forestPlot.loadMarkerFile(new Thread(new Runnable() {
		@Override
		public void run() {
			reloadData();
			forestPlot.getForestPanel().setPointsGeneratable(forestPlot.getMarkerFileName() != null);
			forestPlot.getForestPanel().setRectangleGeneratable(forestPlot.getMarkerFileName() != null);
			forestPlot.getForestPanel().setExtraLayersVisible(new byte[] { 99 });
			try {
				SwingUtilities.invokeAndWait(new Runnable() {
					@Override
					public void run() {
						forestPlot.updateGUI();
						displayIndex(navigationField);
						progressBar.setStringPainted(false);
						progressBar.setValue(0);
						progressBar.setString(null);
						progressBar.setIndeterminate(false);
						progressBar.repaint();
					}
				});
			} catch (InvocationTargetException e) {
				// TODO Auto-generated catch block
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
			}
			forestPlot.setLoadingFile(false);
		}
	}, "ForestPlot_loadMarkerFile"));
	}

	private void reloadData() {
		try {
			SwingUtilities.invokeAndWait(new Runnable() {
				@Override
				public void run() {
					progressBar.setStringPainted(true);
					progressBar.setString("Loading marker list...");
					progressBar.setValue(0);
					progressBar.setIndeterminate(true);
					progressBar.repaint();
				}
			});
		} catch (InvocationTargetException e) {
			// TODO Auto-generated catch block
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
		}
		
		try {
			forestPlot.reloadData();
		} catch (InterruptedException ie) {
			interruptLoading();
		}
		
		try {
			SwingUtilities.invokeAndWait(new Runnable() {
				@Override
				public void run() {
					progressBar.setValue(0);
					progressBar.setStringPainted(true);
					progressBar.setString("Displaying data...");
					progressBar.setIndeterminate(true);
					progressBar.repaint();
				}
			});
		} catch (InvocationTargetException e) {
			// TODO Auto-generated catch block
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
		}

		
	}

	private void generateStudyNameFile() {
        String fileName = ext.verifyDirFormat(ext.parseDirectoryOfFile(forestPlot.getMarkerFileName())) + "sort.txt";
        ArrayList<StudyData> currentData = forestPlot.getCurrentMetaStudy().getStudies();
        ArrayList<String> names = new ArrayList<String>();
        for (StudyData sd : currentData) {
            names.add(sd.getLabel());
        }
        Files.writeArrayList(names, fileName);
    }

    private void updateForestPlot(){
		forestPlot.getForestPanel().setPointsGeneratable(true);
		forestPlot.getForestPanel().setRectangleGeneratable(true);
		forestPlot.getForestPanel().setExtraLayersVisible(new byte[] { 99 });
		displayIndex(navigationField);
		forestPlot.updateGUI();
	}
	
	private void chooseNewFiles() {
		JFileChooser jfc = new JFileChooser((proj != null ? proj.PROJECT_DIRECTORY.getValue() : ext.parseDirectoryOfFile(forestPlot.getMarkerFileName())));
		jfc.setMultiSelectionEnabled(true);
		if (jfc.showOpenDialog(ForestPlotFrame.this) == JFileChooser.APPROVE_OPTION) {
			File[] files = jfc.getSelectedFiles();
			if (files.length > 0) {
				boolean[] keep = Array.booleanArray(files.length, true);
				for (int i = 0; i < files.length; i++) {
					for (String fileName : markerFileNameLoc.keySet()) {
						if (ext.rootOf(files[i].toString()).equals(fileName)) {
							keep[i] = false;
						}
					}
				}
				File[] keptFiles = Array.subArray(files, keep);
				File[] discards = Array.subArray(files, Array.booleanNegative(keep));
				
				if (discards.length > 0) {
					StringBuilder msg = new StringBuilder("The following data file(s) are already present:");
					for (File disc : discards) {
						msg.append("\n").append(disc.getName());
					}
					JOptionPane.showMessageDialog(ForestPlotFrame.this, msg.toString()); 
				}
				
				for (File kept : keptFiles) {
					addFileToList(kept.getAbsolutePath());
				}
			} else {
				File file = jfc.getSelectedFile();
				boolean keep = true;
				for (String fileName : markerFileNameLoc.keySet()) {
					if (ext.rootOf(file.toString()).equals(fileName)) {
						keep = false;
					}
				}
				
				if (!keep) {
					StringBuilder msg = new StringBuilder("The following data file is already present:\n").append(file.getName());
					JOptionPane.showMessageDialog(ForestPlotFrame.this, msg.toString()); 
				} else {
					addFileToList(file.getAbsolutePath());
				}
				
			}
			
		}
	}
	
	private void addFileToList(String rawfile) {
		String file = ext.verifyDirFormat(rawfile);
		file = file.substring(0, file.length() - 1);
		int num = markerFileList.getModel().getSize() - 2;
		Vector<String> currFiles = new Vector<String>();
		if (num > 0) {
			for (int i = 1; i < num + 1; i++) {
				currFiles.add((String)markerFileList.getModel().getElementAt(i));
			}
		}
		String name = ext.rootOf(file);
		markerFileNameLoc.put(name, file);
		currFiles.add(name);
		currFiles.add(0, MARKER_LIST_PLACEHOLDER);
		currFiles.add(MARKER_LIST_NEW_FILE);
		
		markerFileList.setModel(new DefaultComboBoxModel(currFiles));
	}
	
	
	
	private void screenCapAll() {
		int currentSelection = forestPlot.getCurrentDataIndex();
		forestPlot.screenCapAll(null, forestPlot.getForestPanel().oddsDisplay, true, forestPlot.getForestPanel().getSize());
		forestPlot.setCurrentData(currentSelection);
		updateForestPlot();
	}
	
	private void screenCap() {
		forestPlot.screenCap(null, true, forestPlot.getForestPanel().getSize());
	}
	
	private void displayIndex(JTextField field) {
		if (forestPlot.getDataIndices().size() == 0) {
			field.setText("0 of 0");
		} else {
			field.setText((forestPlot.getCurrentDataIndex() + 1) + " of " + forestPlot.getDataIndices().size());
		}
	}

	private void inputMapAndActionMap() {
		InputMap inputMap = forestPlot.getForestPanel().getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW);
		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_UP, InputEvent.ALT_MASK), FIRST);
		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_DOWN, InputEvent.ALT_MASK), LAST);
		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_LEFT, InputEvent.ALT_MASK), PREVIOUS);
		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_RIGHT, InputEvent.ALT_MASK), NEXT);
		ActionMap actionMap = forestPlot.getForestPanel().getActionMap();
		actionMap.put(FIRST, navFirstAction);
		actionMap.put(LAST, navLastAction);
		actionMap.put(PREVIOUS, navPrevAction);
		actionMap.put(NEXT, navNextAction);
		forestPlot.getForestPanel().setActionMap(actionMap);
	}
	
	private void interruptLoading() {
		try {
			SwingUtilities.invokeAndWait(new Runnable() {
				@Override
				public void run() {
					progressBar.setValue(0);
					progressBar.setStringPainted(false);
					progressBar.setString("");
					progressBar.setIndeterminate(false);
					progressBar.repaint();
					markerFileList.setSelectedItem(ForestPlotFrame.MARKER_LIST_PLACEHOLDER);
					forestPlot.getForestPanel().paintAgain();
					repaint();
				}
			});
		} catch (InvocationTargetException e) {
			// TODO Auto-generated catch block
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
		}
	}

	AbstractAction navFirstAction = new AbstractAction() {
		private static final long serialVersionUID = 1L;
		@Override
		public void actionPerformed(ActionEvent e) {
			if(forestPlot.getCurrentDataIndex() != 0){
				forestPlot.setCurrentData(0);
				updateForestPlot();
			}
		}
	};

	AbstractAction navPrevAction = new AbstractAction() {
		private static final long serialVersionUID = 1L;
		@Override
		public void actionPerformed(ActionEvent e) {
			if(forestPlot.getCurrentDataIndex() != 0){
				forestPlot.setCurrentData(forestPlot.getCurrentDataIndex() - 1);
				updateForestPlot();
			}
		}
	};

	AbstractAction navNextAction = new AbstractAction() {
		private static final long serialVersionUID = 1L;
		@Override
		public void actionPerformed(ActionEvent e) {
			if(forestPlot.getCurrentDataIndex() < forestPlot.getDataIndices().size() - 1){
				forestPlot.setCurrentData(forestPlot.getCurrentDataIndex() + 1);
				updateForestPlot();
			}
		}
	};

	AbstractAction navLastAction = new AbstractAction() {
		private static final long serialVersionUID = 1L;
		@Override
		public void actionPerformed(ActionEvent e) {
			if(forestPlot.getCurrentDataIndex() < forestPlot.getDataIndices().size() - 1){
				forestPlot.setCurrentData(forestPlot.getDataIndices().size() - 1);
				updateForestPlot();
			}
		}
	};

	private JProgressBar progressBar;

	@Override
	public void windowOpened(WindowEvent e) {/**/}

	@Override
	public void windowClosing(WindowEvent e) {
		if (proj == null) {
			this.setVisible(false);
			this.dispose();
			return;
		}
		
//		String[] projFiles = proj.getFilenames(Project.FOREST_PLOT_FILENAMES);
		String[] projFiles = proj.FOREST_PLOT_FILENAMES.getValue();
		String[] currFiles = markerFileNameLoc.values().toArray(new String[]{});
		
		ArrayList<String> newFiles = new ArrayList<String>();
		for (String file : currFiles) {
			boolean found = false;
			for (String oldFile : projFiles) {
				if (oldFile.equals(file)) {
					found = true;
				}
			}
			if (!found) {
				newFiles.add(file);
			}
		}
		
		if (newFiles.size() == 0) {
			this.setVisible(false);
			this.dispose();
			return;
		}
//		String newProp = Array.toStr(currFiles, ";");
		
		String message = newFiles.size() + " files have been added.  ";
		int choice = JOptionPane.showOptionDialog(null, message+" Would you like to keep this configuration for the next time Forest Plot is loaded?", "Preserve Forest Plot workspace?", JOptionPane.YES_NO_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE, null, null, null);
		if (choice == 0) {
//			proj.setProperty(Project.FOREST_PLOT_FILENAMES, newProp);
//			proj.setProperty(Project.FOREST_PLOT_FILENAMES, currFiles);
			proj.FOREST_PLOT_FILENAMES.setValue(currFiles);
			proj.saveProperties();
		} else if (choice == -1 || choice == 2) {
			return;
		}

		this.setVisible(false);
		this.dispose();
	}

	@Override
	public void windowClosed(WindowEvent e) {/**/}

	@Override
	public void windowIconified(WindowEvent e) {/**/}

	@Override
	public void windowDeiconified(WindowEvent e) {/**/}

	@Override
	public void windowActivated(WindowEvent e) {/**/}

	@Override
	public void windowDeactivated(WindowEvent e) {/**/}

	public static void main(String[] args) {
				int numArgs = args.length;
		//		String betaSource = "SeqMeta_results.csv";
				String markerList = "markersToDisplay.txt";
	//			String sourceType = "SeqMeta";
				String filename = null;
				String logfile = null;
				final Logger log;
		
				String usage = "\n" + "cnv.plots.ForestPlot requires 1 arguments\n" +
						"  (1) Name of the file with the list of markers (SeqMeta format), files, and comments, to display (i.e. markerList="+ markerList +" (default))\n" +
						"OR\n" +
						"  (1) project properties filename (i.e. proj="+org.genvisis.cnv.Launch.getDefaultDebugProjectFile(false)+" (default))\n";
		
				for (int i = 0; i < args.length; i++) {
					if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
						System.err.println(usage);
						System.exit(1);
		//			} else if (args[i].startsWith("betaSource=")) {
		//				betaSource = args[i].split("=")[1];
		//				numArgs--;
	//				} else if (args[i].startsWith("type=")) {
	//					sourceType = args[i].split("=")[1];
	//					numArgs--;
					} else if (args[i].startsWith("proj=")) {
						filename = args[i].split("=")[1];
						numArgs--;
					} else if (args[i].startsWith("markerList=") || args[i].startsWith("file=")) {
						markerList = args[i].split("=")[1];
						numArgs--;
					} else if (args[i].startsWith("log=")) {
						logfile = args[i].split("=")[1];
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
					final Project proj = (filename != null ? new Project(filename, false) : null);
					log = new Logger(logfile);
		
//					markerList = "list.txt";
					
		//			final String finalDataFile = betaSource;
					final String finalMarkerFile = markerList;
					SwingUtilities.invokeLater(new Runnable() {
						public void run() {
							if (proj != null) {
								new ForestPlotFrame(proj);
							} else {
								new ForestPlotFrame(finalMarkerFile, log);
							}
						}
					});
				} catch (Exception e) {
					e.printStackTrace();
				}
			}

}
