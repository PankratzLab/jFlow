package org.genvisis.one.ben.fcs;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.InputEvent;
import java.awt.event.KeyEvent;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Properties;

import javax.swing.AbstractAction;
import javax.swing.ActionMap;
import javax.swing.InputMap;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLayeredPane;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JSplitPane;
import javax.swing.KeyStroke;
import javax.swing.SwingUtilities;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.xml.parsers.ParserConfigurationException;

import org.genvisis.cnv.gui.GuiManager;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.one.ben.fcs.AbstractPanel2.AXIS_SCALE;
import org.genvisis.one.ben.fcs.AbstractPanel2.PLOT_TYPE;
import org.genvisis.one.ben.fcs.FCSDataLoader.DATA_SET;
import org.genvisis.one.ben.fcs.FCSDataLoader.LOAD_STATE;
import org.genvisis.one.ben.fcs.FCSPanel.GATING_TOOL;
import org.genvisis.one.ben.fcs.gating.Gate;
import org.genvisis.one.ben.fcs.gating.GateDimension;
import org.genvisis.one.ben.fcs.gating.GateFileReader;
import org.genvisis.one.ben.fcs.gating.GateFileWriter;
import org.genvisis.one.ben.fcs.gating.GateTreePanel;
import org.genvisis.one.ben.fcs.gating.Gating;
import org.genvisis.one.ben.fcs.gating.Workbench;
import org.genvisis.one.ben.fcs.sub.DataExportGUI;
import org.genvisis.one.ben.fcs.sub.EMModel;
import org.xml.sax.SAXException;

public class FCSPlot extends JPanel implements WindowListener, PropertyChangeListener {

	public static final String HISTOGRAM_COL = "Histogram";

	public static final int START_X = 20;
	public static final int START_Y = 20;
	public static final int START_WIDTH = 1000;
	public static final int START_HEIGHT = 600;

	public static final long serialVersionUID = 1L;
	public static final Color BACKGROUND_COLOR = Color.WHITE;
	public static final String[] FILE_EXTENSIONS = {".fcs"};
	private FCSPanel fcsPanel;
	private FCSPlotControlPanel fcsControls;
	private JLayeredPane layeredPane;
	private GateTreePanel gatingSelector;

	HashSet<String> validExts;
	Logger log;

	transient FCSDataLoader dataLoader;
	private transient Workbench workbench = new Workbench();
	private String currentSampleID = null;
	private transient Gate parentGate = null;

	private volatile String xDataName;
	private volatile String yDataName;

	private volatile boolean showSDY = true;
	private volatile boolean showSDX = true;
	private volatile boolean showMedianY = true;
	private volatile boolean showMedianX = true;

	private volatile boolean backgating = false;
	private volatile boolean leafgating = false;
	private volatile boolean drawPolyGatesBinned = false;

	private JFrame parentFrame;

	private static final String TITLE_STR = "jFlow";

	HashSet<String> propsSetting = new HashSet<>();

	public static final String PROPERTIES_FILE = "jFlow.properties";

	private static final String PROPKEY_GATEFILE = "GATE_FILE";
	private static final String PROPKEY_FCSFILES = "FCS_FILES";

	private static final Object PROPS_SYNC = new Object();

	protected void saveProps() {
		synchronized (PROPS_SYNC) {
			try {
				Properties props = new Properties();
				if (getGatingStrategy() != null) {
					props.setProperty(PROPKEY_GATEFILE, getGatingStrategy().getFile() == null ? "" : getGatingStrategy().getFile());
				}
				ArrayList<String> files = fcsControls.getAddedFiles();
				props.setProperty(PROPKEY_FCSFILES, files.isEmpty() ? "" : Array.toStr(Array.toStringArray(files), ";"));
				File f = new File(PROPERTIES_FILE);
				OutputStream out = new FileOutputStream(f);
				props.store(out, "");
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}

	private void loadProps() {
		synchronized (PROPS_SYNC) {
			Properties props = new Properties();
			InputStream is = null;

			try {
				File f = new File(PROPERTIES_FILE);
				is = new FileInputStream(f);
				props.load(is);
				String gateFile = props.getProperty(PROPKEY_GATEFILE);
				String fcsTemp = props.getProperty(PROPKEY_FCSFILES);

				if (gateFile != null && !"".equals(gateFile)) {
					fcsControls.loadGatingFile(gateFile);
				}
				if (fcsTemp != null && !"".equals(fcsTemp)) {
					String[] fcs = fcsTemp.split(";");
					fcsControls.addFCSFiles(fcs);
				}

			} catch (Exception e) {
				is = null;
			}
		}
	}

	private FCSPlot() {
		this(null);
	}

	private FCSPlot(String[] fileExts) {
		log = new Logger();

		validExts = new HashSet<>();

		if (fileExts != null) {
			for (String ext : fileExts) {
				validExts.add(ext);
			}
		}

		setLayout(new BorderLayout());

		fcsPanel = new FCSPanel(this);
		fcsPanel.addPropertyChangeListener(this);

		gatingSelector = new GateTreePanel(this);

		layeredPane = new JLayeredPane() {
			private static final long serialVersionUID = 1L;

			@Override
			public boolean isOptimizedDrawingEnabled() {
				return false; // override to force proper z-ordered painting
			}
		};
		layeredPane.setLayout(new BorderLayout());

		layeredPane.add(fcsPanel, BorderLayout.CENTER);
		layeredPane.add(gatingSelector, BorderLayout.NORTH);
		layeredPane.setComponentZOrder(fcsPanel, 0);
		layeredPane.setPreferredSize(new Dimension(1000, 600));

		Dimension minimumSize = new Dimension(100, 50);
		layeredPane.setMinimumSize(minimumSize);

		inputMapAndActionMap();

		fcsPanel.setExtraLayersVisible(new byte[] {99});
		fcsPanel.setChartType(AbstractPanel2.PLOT_TYPE.HEATMAP);

		fcsControls = new FCSPlotControlPanel(this);
		fcsControls.addPropertyChangeListener(this);

		JSplitPane splitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, fcsControls, layeredPane);
		splitPane.setBackground(Color.WHITE);
		splitPane.setOneTouchExpandable(true);
		splitPane.setDividerLocation(250);

		add(splitPane, BorderLayout.CENTER);

		updateGUI();

		fcsPanel.grabFocus();

		loadProps();
	}

	public FCSPanel getPanel() {
		return fcsPanel;
	}

	private void inputMapAndActionMap() {
		InputMap inputMap = fcsPanel.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW);

		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_P, InputEvent.ALT_MASK), "PLOT");
		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_G, InputEvent.ALT_MASK), "GATE");
		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_D, InputEvent.ALT_MASK), "DATA");
		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_O, InputEvent.CTRL_DOWN_MASK), "OPEN_FILE");

		ActionMap actionMap = fcsPanel.getActionMap();

		actionMap.put("PLOT", new AbstractAction() {
			/**
			* 
			*/
			private static final long serialVersionUID = 1L;

			@Override
			public void actionPerformed(ActionEvent arg0) {
				fcsControls.showPlotControls();
			}
		});
		actionMap.put("GATE", new AbstractAction() {
			/**
			* 
			*/
			private static final long serialVersionUID = 1L;

			@Override
			public void actionPerformed(ActionEvent arg0) {

				String g1 = "ID2134458972";
				String g2 = "ID1343278548";
				
				ArrayList<Gate> gates = workbench.getSample(currentSampleID).gating.getGatesForParam("Comp-BV 605-A (CD95)");
				
				
				Gate gt1 = workbench.getSample(currentSampleID).gating.gateMap.get(g1);
				Gate gt2 = workbench.getSample(currentSampleID).gating.gateMap.get(g2);
				
				boolean[] gating11 = gt1.gate(dataLoader);
				boolean[] gating12 = gt1.getParentGating(dataLoader);
				boolean[] gating121 = gt1.getParentGate().getParentGate().gate(dataLoader);
				boolean[] gating21 = gt2.gate(dataLoader);
				boolean[] gating22 = gt2.getParentGating(dataLoader);
				boolean[] gating221 = gt2.getParentGate().getParentGate().gate(dataLoader);
				
//				gating13 should equals gating22
				int s1, s2, s3, s4, s5, s6;
				s1 = Array.booleanArraySum(gating11);
				s2 = Array.booleanArraySum(gating12);
				s4 = Array.booleanArraySum(gating21);
				s5 = Array.booleanArraySum(gating22);
				s6 = Array.booleanArraySum(gating221);
				
				System.out.println(s1 / (double) s2);
				System.out.println(s4 / (double) s5);
				System.out.println();
				
//				fcsControls.showGateControls();
			}
		});
		actionMap.put("DATA", new AbstractAction() {
			/**
			* 
			*/
			private static final long serialVersionUID = 1L;

			@Override
			public void actionPerformed(ActionEvent arg0) {
				fcsControls.showDataControls();
			}
		});
		actionMap.put("OPEN_FILE", new AbstractAction() {
			/**
			* 
			*/
			private static final long serialVersionUID = 1L;

			@Override
			public void actionPerformed(ActionEvent e) {
				fcsControls.dirSelectListener.actionPerformed(e);
			}
		});


		fcsPanel.setActionMap(actionMap);
	}

	private void openFile(ActionEvent e) {
		fcsControls.dirSelectListener.actionPerformed(e);
	}
	
	private void saveImage() {
		JFileChooser fileChooser = new JFileChooser(".");
		int fileOpenActionSelected = fileChooser.showSaveDialog(FCSPlot.this);
		if (fileOpenActionSelected == JFileChooser.APPROVE_OPTION) {
			File fileToOpen = fileChooser.getSelectedFile();
			fcsPanel.screenCapture(fileToOpen.toString() + ".png");
		}
	}
	
	private void doClose() {
		setVisible(false);
		Component parent = getParentComponent();
		parent.dispatchEvent(new WindowEvent((JFrame) parent, WindowEvent.WINDOW_CLOSING));
	}
	
	private JMenuBar menuBar() {
		JMenuBar menuBar;
		JMenu menu;
		JMenuItem menuItemExit, menuItemOpen, menuItemExport, menuItemDump, menuItemEM, menuItemSave;

		menuBar = new JMenuBar();
		menu = new JMenu("File");
		menu.setMnemonic(KeyEvent.VK_F);
		menuBar.add(menu);
		menuItemOpen = new JMenuItem("Open File", KeyEvent.VK_O);
		menuItemOpen.addActionListener(this::openFile);
		menu.add(menuItemOpen);
		menuItemSave = new JMenuItem("Save Current Image", KeyEvent.VK_S);
		menuItemSave.addActionListener(e -> saveImage());
		menu.add(menuItemSave);
		
		menuItemDump = new JMenuItem("Dump Data", KeyEvent.VK_D);
		menuItemDump.addActionListener(e -> dumpData());
		menu.add(menuItemDump);

		menuItemExport = new JMenuItem("Export Data", KeyEvent.VK_E);
		menuItemExport.addActionListener(e -> setupDataExport());
		menu.add(menuItemExport);

		menuItemEM = new JMenuItem("Run EM", KeyEvent.VK_M);
		menuItemEM.addActionListener(e -> setupEM());
		menu.add(menuItemEM);

		menuItemExit = new JMenuItem("Close", KeyEvent.VK_C);
		menuItemExit.addActionListener(e -> SwingUtilities.invokeLater(this::doClose));
		menu.add(menuItemExit);

		return menuBar;
	}

	private Component getParentComponent() {
		Component c = FCSPlot.this.getParent();
		while (!(c instanceof JFrame)) {
			c = c.getParent();
		}
		return c;
	}

	public String getXDataName() {
		return xDataName;
	}

	public String getYDataName() {
		return yDataName;
	}

	public AXIS_SCALE getXScale() {
		return fcsPanel.getXAxis();
	}

	public AXIS_SCALE getYScale() {
		return fcsPanel.getYAxis();
	}

	public PLOT_TYPE getPlotType() {
		return fcsPanel.chartType;
	}

	public boolean showMedian(boolean yAxis) {
		return yAxis ? showMedianY : showMedianX;
	}

	public boolean showSD(boolean yAxis) {
		return yAxis ? showSDY : showSDX;
	}

	public void setXDataName(String xDataName) {
		this.xDataName = xDataName;
		if (!xDataName.equals(fcsControls.getSelectedX())) {
			fcsControls.setXData(xDataName);
		}
		if (dataLoader != null) {
		  fcsControls.setScale(dataLoader.getScaleForParam(xDataName), true);
		  setXScale(dataLoader.getScaleForParam(xDataName));
		}
	}

	public void setYDataName(String yDataName) {
		this.yDataName = yDataName;
		if (!yDataName.equals(fcsControls.getSelectedY())) {
			fcsControls.setYData(yDataName);
		}
		if (dataLoader != null) {
		  fcsControls.setScale(dataLoader.getScaleForParam(yDataName), false);
		  setYScale(dataLoader.getScaleForParam(yDataName));
		}
	}

	public void setXScale(AXIS_SCALE scale) {
		switch (scale) {
			case BIEX:
				fcsPanel.setForcePlotXMin(Float.NaN);
				break;
			case LOG:
			case LIN:
			default:
				fcsPanel.setForcePlotXMin(-1);
				break;
		} 
		fcsPanel.setXAxis(scale);
        if (dataLoader != null) {
          dataLoader.setScaleForParam(getXDataName(), scale);
        }
	}

	public void setYScale(AXIS_SCALE scale) {
		switch (scale) {
			case BIEX:
				fcsPanel.setForcePlotYMin(Float.NaN);
				break;
			case LOG:
			case LIN:
			default:
				fcsPanel.setForcePlotYMin(-1);
				break;
		}
		fcsPanel.setYAxis(scale);
		if (dataLoader != null) {
		  dataLoader.setScaleForParam(getYDataName(), scale);
		}
	}  

	public void setPlotType(PLOT_TYPE type) {
		fcsPanel.chartType = type;
		if (type == PLOT_TYPE.HISTOGRAM) {
			setYDataName(HISTOGRAM_COL);
			fcsControls.setYData(HISTOGRAM_COL);
		} else {
			if (HISTOGRAM_COL.equals(fcsControls.getSelectedY())) {
				fcsControls.setYData(1);
			}
		}
		if (fcsControls.getPlotType() != type) {
			fcsControls.setPlotType(type);
		}
	}

	public void setMedianVisible(boolean show, boolean yAxis) {
		if (yAxis) {
			showMedianY = show;
		} else {
			showMedianX = show;
		}
	}

	public void setSDVisible(boolean show, boolean yAxis) {
		if (yAxis) {
			showSDY = show;
		} else {
			showSDX = show;
		}
	}

	public boolean isBackgating() {
		return backgating;
	}

	public void setBackgating(boolean back) {
		backgating = back;
		fcsPanel.setForceGatesChanged();
		updateGUI();
	}
	
	public boolean isLeafgating() {
	  return leafgating;
	}

	public void setLeafgating(boolean leaf) {
	  leafgating = leaf;
	  fcsPanel.setForceGatesChanged();
	  updateGUI();
	}

	public boolean isDrawPolysAsFlowJo() {
		return drawPolyGatesBinned;
	}

	public void setDrawPolysAsFlowJo(boolean draw) {
		drawPolyGatesBinned = draw;
	}

	public void updateGUI() {
		fcsPanel.paintAgain();
	}

	@Override
	public void windowActivated(WindowEvent e) {
		// unused
	}

	@Override
	public void windowClosed(WindowEvent e) {
		// unused
	}

	@Override
	public void windowClosing(WindowEvent e) {
		if (e == null) {
			GuiManager.disposeOfParentFrame(this);
			return;
		}

		// no way to save loaded vars

		GuiManager.disposeOfParentFrame(this);
	}

	@Override
	public void windowDeactivated(WindowEvent e) {
		// unused
	}

	@Override
	public void windowDeiconified(WindowEvent e) {
		// unused
	}

	@Override
	public void windowIconified(WindowEvent e) {
		// unused
	}

	@Override
	public void windowOpened(WindowEvent e) {
		// unused
		}

	public List<String> getAddedFiles() {
		return fcsControls.getAddedFiles();
	}

	public String getCurrentFile() {
		if (dataLoader == null) {
			return null;
		} else {
			return dataLoader.getLoadedFile();
		}
	}

	public boolean[] getParentGating() {
		boolean[] gating = null;
		if (parentGate != null) {
			gating = parentGate.gate(dataLoader);
		}
		return gating;
	}

	public double[] getAxisData(boolean wait, boolean xAxis) {
		double[] data;
		if (dataLoader == null) {
			data = null;
		} else {
			data = dataLoader.getData(xAxis ? getXDataName() : getYDataName(), wait);
			if (!backgating && !leafgating && parentGate != null) {
				boolean[] gating = parentGate.gate(dataLoader);
				data = Array.subArray(data, gating);
				if (fullClusterAssigns != null) {
				  clusterAssigns = Array.subArray(fullClusterAssigns, gating);
				}
			}
		}
		return data;
	}

	private void resetForNewData(final FCSDataLoader newDataLoader) {
		if (newDataLoader.getLoadState() != LOAD_STATE.UNLOADED
				&& newDataLoader.getLoadState() != LOAD_STATE.LOADING) {
			final ArrayList<String> colNames = newDataLoader.getAllDisplayableNames(DATA_SET.ALL);
			final ArrayList<String> colNamesY = newDataLoader.getAllDisplayableNames(DATA_SET.ALL);
			colNamesY.add(0, HISTOGRAM_COL);

			if (colNames.size() < 2) {
				// TODO error, not enough data!!
			}
			fcsControls.setPlotType(PLOT_TYPE.HEATMAP);
			fcsControls.setColumns(colNames.toArray(new String[colNames.size()]), true, 1);
			fcsControls.setColumns(colNamesY.toArray(new String[colNamesY.size()]), false, 1);
			fcsControls.setScale(newDataLoader.getScaleForParam(colNames.get(0)), false);
			fcsControls.setScale(newDataLoader.getScaleForParam(colNames.get(1)), true);

			setYDataName(colNames.get(0));
			setXDataName(colNames.get(1));
			setYScale(newDataLoader.getScaleForParam(colNames.get(0)));
			setXScale(newDataLoader.getScaleForParam(colNames.get(1)));
			
//			fcsControls.setPlotType(PLOT_TYPE.DOT_PLOT);
//			fcsControls.setScale(AXIS_SCALE.BIEX, false);
//			fcsControls.setScale(AXIS_SCALE.BIEX, true);
//			setYDataName("Comp-APC-A (CD3)");
//			setXDataName("Comp-PE-Cy7-A (CD19)");
//			setYScale(AXIS_SCALE.BIEX);
//			setXScale(AXIS_SCALE.BIEX);
		} else {
			new Thread(new Runnable() {
				@Override
				public void run() {
					while (newDataLoader.getLoadState() == LOAD_STATE.UNLOADED
									|| newDataLoader.getLoadState() == LOAD_STATE.LOADING) {
						Thread.yield();
					}
					final ArrayList<String> colNames = newDataLoader.getAllDisplayableNames(DATA_SET.ALL);
					final ArrayList<String> colNamesY = newDataLoader.getAllDisplayableNames(DATA_SET.ALL);
					colNamesY.add(0, HISTOGRAM_COL);

					if (colNames.size() < 2) {
						// TODO error, not enough data!!
					}
					// TODO reset GUI elements
					// if (resetCols) {
					SwingUtilities.invokeLater(new Runnable() {
						@Override
						public void run() {
							fcsControls.setPlotType(PLOT_TYPE.HEATMAP);
							fcsControls.setColumns(colNames.toArray(new String[colNames.size()]), true, 1);
							fcsControls.setColumns(colNamesY.toArray(new String[colNamesY.size()]), false, 1);
							fcsControls.setScale(newDataLoader.getScaleForParam(colNames.get(0)), false);
							fcsControls.setScale(newDataLoader.getScaleForParam(colNames.get(1)), true);
						}
					});

					setYDataName(colNames.get(0));
					setXDataName(colNames.get(1));
					setYScale(newDataLoader.getScaleForParam(colNames.get(0)));
					setXScale(newDataLoader.getScaleForParam(colNames.get(0)));

				}
			}).start();
		}
		dataLoader = newDataLoader;
		System.gc();
	}

	private final HashMap<String, FCSDataLoader> loadedData = new HashMap<>();

	public void unloadFile(final String filename) {
		if (filename == null || !loadedData.containsKey(filename)) {
			return;
		}
		boolean clearLoaded = false;
		if (dataLoader != null && dataLoader.getLoadedFile().equals(filename)) {
			clearLoaded = true;
		}
		loadedData.remove(filename).emptyAndReset();
		if (clearLoaded) {
			dataLoader = null;
		}
		System.gc();
		updateGUI();
	}



	public void loadFile(final String filename, final boolean display) {
		if (filename == null	|| !Files.exists(filename) || (dataLoader != null && dataLoader.getLoadedFile().equals(filename))) {
			return;
		}

		boolean applyTemplate = false;
		currentSampleID = workbench.containsSampleFile(filename) ? workbench.getSampleID(filename) : workbench.addNewSample(filename, applyTemplate);
    refreshGating();
		
		if (loadedData.containsKey(filename)) {
			if (display) {
				setData(loadedData.get(filename));
			}
		} else {
			Thread dataLoaderThread = new Thread(() -> {
				FCSDataLoader newDataLoader = new FCSDataLoader();
				loadedData.put(filename, newDataLoader);
				fcsControls.startFileLoading(newDataLoader);
				try {
					newDataLoader.loadData(filename);
				} catch (IOException e) {
					log.reportException(e);
					return;
				}
				if (display) {
					setData(newDataLoader);
				}
			});
			dataLoaderThread.start();
		}
	}
  
    public void loadWorkspaceFile(String gateFile) {
      try {
        this.workbench = GateFileReader.loadWorkspace(gateFile, dataLoader);
        saveProps();
    		refreshGating();
      } catch (ParserConfigurationException | SAXException | IOException e) {
        log.reportException(e);
      }
    }

	public void saveGating() {
		if (getGatingStrategy().getRootGates().isEmpty()) {
			JOptionPane.showMessageDialog(this, "Error - no gates found!", "Error!", JOptionPane.ERROR_MESSAGE);
			return;
		}
		JFileChooser jfc = new JFileChooser();
		jfc.setFileFilter(new FileNameExtensionFilter("Gating-ML File", "xml"));
		jfc.setMultiSelectionEnabled(false);
		int code = jfc.showSaveDialog(this);
		if (code == JFileChooser.APPROVE_OPTION) {
			String outputFile = jfc.getSelectedFile().getAbsolutePath();
			if (!outputFile.endsWith(".xml")) {
				outputFile += ".xml";
			}
			if (Files.exists(outputFile)) {
				JOptionPane.showMessageDialog(this, "Error - cannot overwrite existing file.", "Error!",
																			JOptionPane.ERROR_MESSAGE);
				return;
			}
			GateFileWriter.writeGating(getGatingStrategy(), outputFile, log);
		}
	}

	public Gating getGatingStrategy() {
	  if (currentSampleID != null) {
	  	return workbench.getSample(currentSampleID).getGating();
	  } else {
	    return new Gating();
	  }
	}

	public void clearGating() {
	  if (currentSampleID == null) {
	  	return;
	  }
	  workbench.clearGating(currentSampleID);
    parentGate = null;
    gatingSelector.resetGating(getGatingStrategy(), null);
    fcsPanel.clearGating();
	}

	public void setGating(Gating gateStrat) {
		workbench.setGatingForSample(currentSampleID, gateStrat);
		parentGate = null;
		for (Gate g : getGatingStrategy().getRootGates()) {
			parentGate.getChildGates().add(g);
		}
		gatingSelector.resetGating(gateStrat, null);
		updateGUI();
	}

	public void refreshGating() {
		Gate sel = parentGate;
		gatingSelector.resetGating(getGatingStrategy(), sel);
		fcsPanel.setForceGatesChanged();
		updateGUI();
	}
	
	protected HashMap<Gate, boolean[]> gateAllDataForLeafGates() {
	  HashSet<Gate> leafGates = getGatingStrategy().getAllLeafGates();
	  HashMap<Gate, boolean[]> gatings = new HashMap<>();
	  log.reportTime("Gating on " + leafGates.size() + " leaf gates...");
	  long t1 = System.currentTimeMillis();
	  for (Gate g : leafGates) {
	    gatings.put(g, g.gate(dataLoader));
	  }
	  log.reportTimeElapsed("Gating complete - ", t1);
	  return gatings;
	}

  public List<Gate> getGatingForCurrentPlot() {
    ArrayList<Gate> gateList = new ArrayList<>();
    ArrayList<Gate> children = parentGate == null ? getGatingStrategy().getRootGates() : parentGate.getChildGates();
    for (Gate g : children) {
      boolean y = getYDataName().equals(FCSPlot.HISTOGRAM_COL) ? true : false;
      if (g.getYDimension() == null || y) {
        if (y && g.getXDimension().getParam().equals(getXDataName())) {
          gateList.add(g);
        }
        continue;
      }
      if (g.getXDimension().getParam().equals(getXDataName())
          && g.getYDimension().getParam().equals(getYDataName())) {
        gateList.add(g);
      }
    }
    return gateList;
  }

	public void setData(FCSDataLoader newDataLoader) {
		if (dataLoader != null && dataLoader.getLoadedFile().equals(newDataLoader.getLoadedFile())) {
			return;
		}
		resetForNewData(newDataLoader);
		if (workbench != null) {
			GateFileReader.updateWorkbench(workbench, newDataLoader);
		}
		parentFrame.setTitle(ext.rootOf(newDataLoader.getLoadedFile()) + "  --  " + TITLE_STR);
		refreshGating();
		updateGUI();
	}

	private void setParent(JFrame frame) {
		parentFrame = frame;
	}

	/**
	 * Create the GUI and show it. For thread safety, this method should be invoked from the
	 * event-dispatching thread.
	 * 
	 * @param proj Project
	 * @param show setVisible
	 * @param promptOnClose prompt to save files/vars on close
	 * @param fileExts allowed file extensions for files (null for all files)
	 * @return
	 */
	public static FCSPlot createGUI(boolean show) {
		JFrame frame = new JFrame(TITLE_STR);
		frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

		// Create and set up the content pane.
		FCSPlot twoDPlot = new FCSPlot(/* proj, promptOnClose, fileExts, filenamesProperty */);
		frame.setJMenuBar(twoDPlot.menuBar());
		twoDPlot.setOpaque(true); // content panes must be opaque
		twoDPlot.setParent(frame);
		frame.setContentPane(twoDPlot);
		frame.addWindowListener(twoDPlot);
		// frame.setExtendedState(frame.getExtendedState()|JFrame.MAXIMIZED_BOTH);

		// Display the window.
		frame.pack();
		frame.setBounds(START_X, START_Y, START_WIDTH, START_HEIGHT);
		frame.setVisible(show);

		return twoDPlot;
	}

	public int getDataCount() {
		int cnt = -1;
		if (dataLoader == null) {
			return cnt;
		}
		if (!backgating && !leafgating && parentGate != null) {
			cnt = Array.booleanArraySum(parentGate.gate(dataLoader));
		} else {
			cnt = dataLoader.getCount();
		}

		return cnt;
	}

	@Override
	public void propertyChange(PropertyChangeEvent arg0) {
		if (arg0.getSource().equals(fcsPanel)) {
			if (arg0.getPropertyName().equals(AbstractPanel2.X_MIN)) {
				fcsControls.setXMin(((Double) arg0.getNewValue()).floatValue());
			} else if (arg0.getPropertyName().equals(AbstractPanel2.X_MAX)) {
				fcsControls.setXMax(((Double) arg0.getNewValue()).floatValue());
			} else if (arg0.getPropertyName().equals(AbstractPanel2.Y_MIN)) {
				fcsControls.setYMin(((Double) arg0.getNewValue()).floatValue());
			} else if (arg0.getPropertyName().equals(AbstractPanel2.Y_MAX)) {
				fcsControls.setYMax(((Double) arg0.getNewValue()).floatValue());
			}
			updateGUI();
		}
	}

	public boolean isCurrentDataDisplayable() {
		if (dataLoader == null) {
			return false;
		}
		LOAD_STATE currState = dataLoader.getLoadState();
		return currState != LOAD_STATE.UNLOADED && currState != LOAD_STATE.LOADING;
	}

	public boolean isFileLoaded(String file) {
		return loadedData.containsKey(file);
	}

	public Gate getParentGate() {
		return parentGate;
	}

	public void addGate(Gate rg) {
		if (parentGate == null) {
		  getGatingStrategy().addRootGate(rg);
		} else {
		  parentGate.addChildGate(rg);
		}
		gatingSelector.resetGating(getGatingStrategy(), parentGate);
	}

	public boolean duplicateGateName(String name) {
		return getGatingStrategy().gateNameExists(name);
	}

	public void gateSelected(Gate gate, boolean reset) {
		if (gate == null) {
			parentGate = null;
		} else {
			parentGate = gate;
			GateDimension gdX = gate.getXDimension();
			setXDataName(gdX.getParam());
			if (gate.getYDimension() != null) {
				GateDimension gdY = gate.getYDimension();
				setYDataName(gdY.getParam());
			} else {
				setPlotType(PLOT_TYPE.HISTOGRAM);
			}
		}
		if (reset) {
			gatingSelector.resetGating(getGatingStrategy(), gate);
		}
		fcsPanel.clearGating();
	}

	int[] fullClusterAssigns = null;
	int[] clusterAssigns = null;
	
	public int[] getClusterAssignments() {
	  return clusterAssigns;
	}
	
	private void setupEM() {
	  EMModel model = EMModel.run(dataLoader);
	  fullClusterAssigns = model.getClusterAssigns();
	  clusterAssigns = fullClusterAssigns;
	  updateGUI();
	}
	
	private void dumpData() {
	  boolean[] toDump = getParentGating();
	  ArrayList<String> params = dataLoader.getAllDisplayableNames(DATA_SET.COMPENSATED);
	  PrintWriter writer = Files.getAppropriateWriter(ext.rootOf(dataLoader.getLoadedFile(), false) + ".xln");
	  for (int i = 0, count = params.size(); i < count; i++) {
	    writer.print(params.get(i));
	    if (i < count - 1) {
	      writer.print("\t");
	    }
	  }
	  writer.println();
	  
	  double[][] data = new double[params.size()][];
	  for (int i = 0, count = params.size(); i < count; i++) {
	    data[i] = dataLoader.getData(params.get(i), true);
	  }
  
      for (int i = 0, count = data[0].length; i < count; i++) {
        if (toDump[i]) {
          for (int p = 0, pcount = params.size(); p < pcount; p++) {
            writer.print(data[p][i]);
            if (p < pcount - 1) {
              writer.print("\t");
            }
          }
          writer.println();
        }
      }
	  writer.flush();
	  writer.close();
	}
	
	private void setupDataExport() {
		if (getAddedFiles().isEmpty()) {
			JOptionPane.showMessageDialog(this, "Error - no data files available to export!", "Error!",
																		JOptionPane.ERROR_MESSAGE);
			return;
		}
		if (getGatingStrategy().getRootGates().isEmpty()) {
			JOptionPane.showMessageDialog(this, "Error - no gating available to export!", "Error!",
																		JOptionPane.ERROR_MESSAGE);
			return;
		}

		DataExportGUI degui = new DataExportGUI(this);
		degui.setModal(true);
		degui.setVisible(true);
		degui.pack();
		if (degui.wasCancelled()) {
			return;
		}

		ArrayList<String> files = degui.getSelectedFiles();
		ArrayList<Gate> gates = degui.getSelectedGates();
		boolean writeCounts = degui.getCountsExportSelected();
		String output = degui.getOutputFile();

		doDataExport(output, gates, writeCounts, files);
	}

	private void doDataExport(String outputFile, ArrayList<Gate> gatesToExport, boolean exportCounts,
														ArrayList<String> files) {
		StringBuilder sb = new StringBuilder();

		for (Gate g : gatesToExport) {
			sb.append("\t").append(g.getFullNameAndGatingPath());
		}
		for (String file : files) {
			FCSDataLoader dataLoader;
			boolean wasLoaded = false;
			if (loadedData.containsKey(file)) {
				dataLoader = loadedData.get(file);
				wasLoaded = true;
			} else {
				loadFile(file, false);
				dataLoader = null;
			}
			while (dataLoader == null || dataLoader.getLoadState() != FCSDataLoader.LOAD_STATE.LOADED) {
				try {
					Thread.sleep(200);
				} catch (InterruptedException e) {
				}
				if (dataLoader == null) {
					dataLoader = loadedData.get(file);
				}
			}
			sb.append("\n").append(ext.removeDirectoryInfo(dataLoader.getLoadedFile()));

			for (Gate g : gatesToExport) {
				boolean[] gt = g.gate(dataLoader);
				if (gt == null) {
					sb.append("\t.");
				} else {
					int c = Array.booleanArraySum(gt);
					if (exportCounts) {
						sb.append("\t").append(c);
					} else {
						boolean[] pG = g.getParentGating(dataLoader);
						int p = pG == null ? dataLoader.getCount() : Array.booleanArraySum(pG);
						double pct = (double) c / (double) p;
						sb.append("\t").append(ext.formDeci(100 * pct, 2));
					}
				}
			}
			if (!wasLoaded) {
				unloadFile(file);
			}
		}

		Files.write(sb.toString(), outputFile);
		log.reportTime("Data written to file: " + outputFile);
	}

	public void setGateTool(GATING_TOOL item) {
		fcsPanel.setGatingTool(item);
	}

	public void deleteGate(Gate g) {
	  if (g == parentGate) {
	    parentGate = null;
	  }
      getGatingStrategy().deleteGate(g);
      gatingSelector.resetGating(getGatingStrategy(), parentGate);
	}

	public static void main(String[] args) {
		javax.swing.SwingUtilities.invokeLater(() -> {
			createGUI(true);
		});
	}

}


