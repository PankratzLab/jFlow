package one.ben.fcs;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import javax.swing.ActionMap;
import javax.swing.InputMap;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLayeredPane;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JSplitPane;
import javax.swing.SwingUtilities;
import javax.xml.parsers.ParserConfigurationException;

import one.ben.fcs.AbstractPanel2.AXIS_SCALE;
import one.ben.fcs.AbstractPanel2.PLOT_TYPE;
import one.ben.fcs.FCSDataLoader.DATA_SET;
import one.ben.fcs.FCSDataLoader.LOAD_STATE;
import one.ben.fcs.gating.Gate;
import one.ben.fcs.gating.Gate.RectangleGate;
import one.ben.fcs.gating.GateDimension;
import one.ben.fcs.gating.GateFileReader;
import one.ben.fcs.gating.GateTreePanel;
import one.ben.fcs.gating.GatingStrategy;
import one.ben.fcs.sub.RainbowTestGUI;

import org.xml.sax.SAXException;

import cnv.gui.GuiManager;
import common.Array;
import common.Files;
import common.Logger;

public class FCSPlot extends JPanel implements WindowListener, ActionListener, PropertyChangeListener { 
    
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

    FCSDataLoader dataLoader;
    private GatingStrategy gating;
    private Gate parentGate = null;

//    volatile boolean isLoading = false;

    private volatile String xDataName;
    private volatile String yDataName;

    private volatile PLOT_TYPE plotType;

    private volatile boolean showSDY = true;
    private volatile boolean showSDX = true;
    private volatile boolean showMedianY = true;
    private volatile boolean showMedianX = true;

    private JFrame parent;

    private static final String TITLE_STR = "Genvisis - FCS Plot";

    HashSet<String> propsSetting = new HashSet<String>();

    private FCSPlot() {
        this(null);
    }
    
	private FCSPlot(String[] fileExts) {
		log = new Logger();

		validExts = new HashSet<String>();
		
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
	}
	
	public FCSPanel getPanel() {
	    return fcsPanel;
	}

	private void inputMapAndActionMap() {
		InputMap inputMap = fcsPanel.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW);
//		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_UP, InputEvent.ALT_MASK), ALT_UP);
//		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_DOWN, InputEvent.ALT_MASK), ALT_DOWN);
//		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_LEFT, InputEvent.ALT_MASK), ALT_LEFT);
//		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_RIGHT, InputEvent.ALT_MASK), ALT_RIGHT);
		ActionMap actionMap = fcsPanel.getActionMap();
		fcsPanel.setActionMap(actionMap);
	}

	private JMenuBar menuBar() {
		JMenuBar menuBar;
		JMenu menu;
		JMenuItem menuItemExit, menuItemOpen, menuItemShowControls, menuItemRemove, menuItemRemoveAll, menuItemScreens, menuItemSave;
		JMenuItem menuItemTestRB;
		final JCheckBoxMenuItem /*menuItemExclude,*/ menuItemHist;
		
		menuBar = new JMenuBar();
		menu = new JMenu("File");
        menu.setMnemonic(KeyEvent.VK_F);
		menuBar.add(menu);
		menuItemOpen = new JMenuItem("Open File", KeyEvent.VK_O);
		menuItemOpen.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
			    addFile();
			}
		});
		menu.add(menuItemOpen);
		menuItemSave = new JMenuItem("Save Current Image", KeyEvent.VK_S);
		menuItemSave.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				JFileChooser fileChooser = new JFileChooser(".");
				int fileOpenActionSelected = fileChooser.showSaveDialog(FCSPlot.this);
		        if (fileOpenActionSelected == JFileChooser.APPROVE_OPTION) {
		            File fileToOpen = fileChooser.getSelectedFile();
		            fcsPanel.screenCapture(fileToOpen.toString()+".png");
		        }
			}
		});
		menu.add(menuItemSave);
		menuItemExit = new JMenuItem("Close", KeyEvent.VK_C);
		menuItemExit.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				SwingUtilities.invokeLater(new Runnable() {
					@Override
					public void run() {
						setVisible(false);
						Component parent = getParentComponent();
						parent.dispatchEvent(new WindowEvent((JFrame) parent, WindowEvent.WINDOW_CLOSING));
					}
				});
			}
		});
		menu.add(menuItemExit);
		
		menu = new JMenu("View");
        menu.setMnemonic(KeyEvent.VK_V);
		menuBar.add(menu);
		
		menuItemTestRB = new JMenuItem("Test Rainbow Bead Files", KeyEvent.VK_R);
		menuItemTestRB.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                runRainbowTest();
            }
        });
		
//		menuItemShowControls = new JMenuItem("Show ControlPanel", KeyEvent.VK_C);
//		menuItemShowControls.addActionListener(new ActionListener() {
//            @Override
//            public void actionPerformed(ActionEvent e) {
//                controlFrame.setBounds(START_X + START_WIDTH + 50, START_Y, 250, 400);
//                controlFrame.setVisible(true);
//            }
//        });
//		menu.add(menuItemShowControls);
		
		return menuBar;
	}
	
	private void runRainbowTest() {
	    RainbowTestGUI rtGui = new RainbowTestGUI();
	    
	    
	    
	}

    private Component getParentComponent() {
        Component c = FCSPlot.this.getParent();
        while (!(c instanceof JFrame)) {
            c = c.getParent();
        }
        return c;
    }
	
	public void actionPerformed(ActionEvent ae) {
		String command;
		
		command = ae.getActionCommand();
//		if (command.equals(ADD_DATA_FILE)) {
//			addFile();
//		} else if (command.equals(CREATE_SCREENS)) {
//		    createScreenshotsFromFile();
//		} else {
//			System.err.println("Error - unknown command '"+command+"'");
//		}
	}
	
	public String getXDataName() { return xDataName; }
    public String getYDataName() { return yDataName; }
    
    public AXIS_SCALE getXScale() { return fcsPanel.getXAxis(); }
    public AXIS_SCALE getYScale() { return fcsPanel.getYAxis(); }
    
    public PLOT_TYPE getPlotType() { return plotType; }

    public boolean showMedian(boolean yAxis) { return yAxis ? showMedianY : showMedianX; }
    public boolean showSD(boolean yAxis) { return yAxis ? showSDY : showSDX; }
    
    public void setXDataName(String xDataName) { 
        this.xDataName = xDataName;
        if (!xDataName.equals(this.fcsControls.getSelectedX())) {
            this.fcsControls.setXData(xDataName);
        }
    }
    public void setYDataName(String yDataName) { 
        this.yDataName = yDataName;
        if (!yDataName.equals(this.fcsControls.getSelectedY())) {
            this.fcsControls.setYData(yDataName);
        }
    }

    public void setXScale(AXIS_SCALE scale) { this.fcsPanel.setXAxis(scale); }
    public void setYScale(AXIS_SCALE scale) { this.fcsPanel.setYAxis(scale); }
	
    public void setPlotType(PLOT_TYPE type) { 
        this.fcsPanel.chartType = type;
        if (type == PLOT_TYPE.HISTOGRAM) {
            this.setYDataName(HISTOGRAM_COL);
            this.fcsControls.setYData(HISTOGRAM_COL);
        } else {
            if (HISTOGRAM_COL.equals(this.fcsControls.getSelectedY())) {
                this.fcsControls.setYData(1);
            }
        }
        if (this.fcsControls.getPlotType() != type) {
            this.fcsControls.setPlotType(type);
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
    
	public ArrayList<String[]> getDataSelected() {
		ArrayList<String[]> v = null;
		
		
		
		return v;
	}

	public void updateGUI() {
        fcsPanel.paintAgain();
	}
	
	public void windowActivated(WindowEvent e) {
	}

	public void windowClosed(WindowEvent e) {}

	@Override
	public void windowClosing(WindowEvent e) {
		if (e == null) {
			GuiManager.disposeOfParentFrame(this);
			return;
		}
		
		// no way to save loaded vars
		
		GuiManager.disposeOfParentFrame(this);
	}

	public void windowDeactivated(WindowEvent e) {}

	public void windowDeiconified(WindowEvent e) {}

	public void windowIconified(WindowEvent e) {}

	public void windowOpened(WindowEvent e) {}

    private void addFile() {
        
    }
    
    public float[] getAxisData(boolean wait, boolean xAxis) {
        float[] data;
        if (dataLoader == null) {
            data = null;
        } else {
            data = dataLoader.getData(xAxis ? getXDataName() : getYDataName(), wait);
        }
        return data;
    }
    
    private void resetForNewData(final FCSDataLoader newDataLoader) {
        if (newDataLoader.getLoadState() != LOAD_STATE.UNLOADED && newDataLoader.getLoadState() != LOAD_STATE.LOADING) {
        	final ArrayList<String> colNames = newDataLoader.getAllDisplayableNames(DATA_SET.ALL);
        	final ArrayList<String> colNamesY = newDataLoader.getAllDisplayableNames(DATA_SET.ALL);
        	colNamesY.add(0, HISTOGRAM_COL);
        	
            if (colNames.size() < 2) {
                // TODO error, not enough data!!
            }
            // TODO reset GUI elements
    //        if (resetCols) {
    //        SwingUtilities.invokeLater(new Runnable() {
    //            @Override
    //            public void run() {
                    fcsControls.setPlotType(PLOT_TYPE.HEATMAP);
                    fcsControls.setColumns(colNames.toArray(new String[colNames.size()]), true, 1);
                    fcsControls.setColumns(colNamesY.toArray(new String[colNamesY.size()]), false, 1);
                    fcsControls.setScale(newDataLoader.scales.get(0), false);
                    fcsControls.setScale(newDataLoader.scales.get(1), true);
    //            }
    //        });
            
            setYDataName(colNames.get(0));
            setXDataName(colNames.get(1));
            setYScale(newDataLoader.scales.get(0));
            setXScale(newDataLoader.scales.get(1));
        } else {
            new Thread(new Runnable() {
                @Override
                public void run() {
                    while(newDataLoader.getLoadState() == LOAD_STATE.UNLOADED || newDataLoader.getLoadState() == LOAD_STATE.LOADING) {
                        Thread.yield();
                    }
                    final ArrayList<String> colNames = newDataLoader.getAllDisplayableNames(DATA_SET.ALL);
                    final ArrayList<String> colNamesY = newDataLoader.getAllDisplayableNames(DATA_SET.ALL);
                    colNamesY.add(0, HISTOGRAM_COL);
                    
                    if (colNames.size() < 2) {
                        // TODO error, not enough data!!
                    }
                    // TODO reset GUI elements
            //        if (resetCols) {
                    SwingUtilities.invokeLater(new Runnable() {
                        @Override
                        public void run() {
                            fcsControls.setPlotType(PLOT_TYPE.HEATMAP);
                            fcsControls.setColumns(colNames.toArray(new String[colNames.size()]), true, 1);
                            fcsControls.setColumns(colNamesY.toArray(new String[colNamesY.size()]), false, 1);
                            fcsControls.setScale(newDataLoader.scales.get(0), false);
                            fcsControls.setScale(newDataLoader.scales.get(1), true);
                        }
                    });
                    
                    setYDataName(colNames.get(0));
                    setXDataName(colNames.get(1));
                    setYScale(newDataLoader.scales.get(0));
                    setXScale(newDataLoader.scales.get(1));
                    
                }
            }).start();
        }
        dataLoader = newDataLoader;
        System.gc();
    }
    
    private HashMap<String, FCSDataLoader> loadedData = new HashMap<String, FCSDataLoader>();
    
    public void unloadFile(final String filename) {
        if (filename == null || !loadedData.containsKey(filename)) {
            return;
        }
        boolean clearLoaded = true;
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
	    if (filename == null || !Files.exists(filename) || (dataLoader != null && dataLoader.getLoadedFile().equals(filename))) {
	        return;
	    }
	    
	    if (loadedData.containsKey(filename)) {
	        if (display) {
	            setData(loadedData.get(filename));
	        }
	    } else {
	        Thread dataLoaderThread = new Thread(new Runnable() {
                @Override
                public void run() {
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
                }
            });
    	    dataLoaderThread.start();
	    }
	}
	
	public void loadGatingFile(String gateFile) {
        try {
            setGating(GateFileReader.readGateFile(gateFile));
        } catch (ParserConfigurationException e) {
            log.reportException(e);
        } catch (SAXException e) {
            log.reportException(e);
        } catch (IOException e) {
            log.reportException(e);
        }
	}
	
	public void setGating(GatingStrategy gateStrat) {
	    this.gating = gateStrat;
	    this.gatingSelector.setGating(gateStrat);
	    // TODO repaint
	}
	
	public ArrayList<Gate> getGatingForCurrentPlot() {
	    ArrayList<Gate> gateList = new ArrayList<Gate>();
	    if (this.parentGate != null) {
	        ArrayList<Gate> children = parentGate.getChildGates();
	        for (Gate g : children) {
	            boolean x = false;
	            boolean y = getYDataName().equals(FCSPlot.HISTOGRAM_COL) ? true : false;
	            if (g.getDimensions().size() == 1) {
	                if (y && g.getDimensions().get(0).getParam().equals(getXDataName())) {
                        gateList.add(g);
	                } 
                    continue;
	            }
	            // >1 Gate Dimension
	            if (y) continue;
	            if (g.getDimensions().size() > 2) continue; // can't vis. multi-dim gates // TODO fix for 3D
	            for (int i = 0; i < g.getDimensions().size(); i++) {
	                GateDimension gd = g.getDimensions().get(i);
	                if (gd.getParam().equals(getXDataName())) {
	                    x = true;
	                }
	                if (gd.getParam().equals(getYDataName())) {
	                    y = true;
	                }
	            }
	            if (x && y) {
	                gateList.add(g);
	            }
	        }
        }
	    return gateList;
	}
	
	public void setData(FCSDataLoader newDataLoader) {
	    if (this.dataLoader != null && this.dataLoader.getLoadedFile().equals(newDataLoader.getLoadedFile())) return;
	    resetForNewData(newDataLoader);
        this.parent.setTitle(TITLE_STR + "  --  " + newDataLoader.getLoadedFile());
        System.out.println("Displaying plot...");
        updateGUI();
	}

	private void setParent(JFrame frame) {
	    this.parent = frame;
	}
	
	/**
     * Create the GUI and show it.  For thread safety,
     * this method should be invoked from the
     * event-dispatching thread.
	 * @param proj Project
	 * @param show setVisible
	 * @param promptOnClose prompt to save files/vars on close
	 * @param fileExts allowed file extensions for files (null for all files)
	 * @return
	 */
	public static FCSPlot createGUI(boolean show) {
		JFrame frame = new JFrame(TITLE_STR);
		frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

        //Create and set up the content pane.
        FCSPlot twoDPlot = new FCSPlot(/*proj, promptOnClose, fileExts, filenamesProperty*/);
        frame.setJMenuBar(twoDPlot.menuBar());
        twoDPlot.setOpaque(true); //content panes must be opaque
        twoDPlot.setParent(frame);
        frame.setContentPane(twoDPlot);
        frame.addWindowListener(twoDPlot);
//		frame.setExtendedState(frame.getExtendedState()|JFrame.MAXIMIZED_BOTH);

        //Display the window.
        frame.pack();
        frame.setBounds(START_X, START_Y, START_WIDTH, START_HEIGHT);
        frame.setVisible(show);

//        String fcsFilename = "F:\\Flow\\2016-05-04_URB_DHS_ULTRA BRIGHT RAINBOW BEADS_URB_001.fcs";
//        String fcsFilename = "F:\\Flow\\P1-B&C-CD3-APC-Cy7 or CD4-APC-Cy7_ULTRA BRIGHT RAINBOW BEADS_URB_001.fcs";
//        String fcsFilename = "F:\\Flow\\P1- PBMC-A&C rest_panel one_PBMC-C P1 1HR rest_003.fcs.gz";
//        String fcsFilename = "F:\\Flow\\P1- PBMC-A&C rest_panel one_PBMC-A P1 1HR rest_002.fcs";
//        twoDPlot.loadFile(fcsFilename);
        
		return twoDPlot;
    }
	
	public static void main(String[] args) {
        javax.swing.SwingUtilities.invokeLater(new Runnable() {
            public void run() {
                createGUI(true);
            }
        });
	}

    public int getDataCount() {
        return dataLoader == null ? -1 : dataLoader.getCount();
    }
    
    @Override
    public void propertyChange(PropertyChangeEvent arg0) {
//        if (propsSetting.contains(arg0.getPropertyName())) return;
//        propsSetting.add(arg0.getPropertyName());
//        if (arg0.getSource().equals(fcsControls)) {
//            if (arg0.getPropertyName().equals(AbstractPanel2.X_MIN)) {
//                fcsPanel.setForcePlotXMin(((Double)arg0.getNewValue()).floatValue());
//                fcsPanel.setPlotXMin(((Double)arg0.getNewValue()).floatValue());
//            } else if (arg0.getPropertyName().equals(AbstractPanel2.X_MAX)) {
//                fcsPanel.setForcePlotXMax(((Double)arg0.getNewValue()).floatValue());
//                fcsPanel.setPlotXMax(((Double)arg0.getNewValue()).floatValue());
//            } else if (arg0.getPropertyName().equals(AbstractPanel2.Y_MIN)) {
//                fcsPanel.setForcePlotYMin(((Double)arg0.getNewValue()).floatValue());
//                fcsPanel.setPlotYMin(((Double)arg0.getNewValue()).floatValue());
//            } else if (arg0.getPropertyName().equals(AbstractPanel2.Y_MAX)) {
//                fcsPanel.setForcePlotYMax(((Double)arg0.getNewValue()).floatValue());
//                fcsPanel.setPlotYMax(((Double)arg0.getNewValue()).floatValue());
//            }
//            updateGUI();
//        } else 
        if (arg0.getSource().equals(fcsPanel)) {
            if (arg0.getPropertyName().equals(AbstractPanel2.X_MIN)) {
                fcsControls.setXMin(((Double)arg0.getNewValue()).floatValue());
            } else if (arg0.getPropertyName().equals(AbstractPanel2.X_MAX)) {
                fcsControls.setXMax(((Double)arg0.getNewValue()).floatValue());
            } else if (arg0.getPropertyName().equals(AbstractPanel2.Y_MIN)) {
                fcsControls.setYMin(((Double)arg0.getNewValue()).floatValue());
            } else if (arg0.getPropertyName().equals(AbstractPanel2.Y_MAX)) {
                fcsControls.setYMax(((Double)arg0.getNewValue()).floatValue());
            }
            updateGUI();
        }
//        SwingUtilities.invokeLater(new Runnable() {
//            @Override
//            public void run() {
//                propsSetting.remove(arg0.getPropertyName());
//            }
//        });
    }

    public boolean isCurrentDataDisplayable() {
        if (dataLoader == null) return false;
        LOAD_STATE currState = dataLoader.getLoadState();
        return currState != LOAD_STATE.UNLOADED && currState != LOAD_STATE.LOADING;
    }

    public boolean isFileLoaded(String file) {
        return loadedData.containsKey(file);
    }

    private static class NullGate extends Gate {
        private NullGate() {
            super(null);
        }

        @Override
        public boolean[] gate(FCSDataLoader dataLoader) {
            return null;
        }
    }
    
    public Gate getParentGate() {
        if (this.parentGate == null || this.parentGate instanceof NullGate) return null;
        return parentGate;
    }
    
    public void addGate(Gate rg) {
        if (parentGate == null) {
            parentGate = new NullGate();
        }
        parentGate.getChildGates().add(rg);
        // TODO export/save gate data
        // TODO addGate
    }

    public void gateSelected(Gate gate) {
        ArrayList<GateDimension> gd = gate.getDimensions();
        GateDimension gdX = gd.get(0);
        GateDimension gdY = gd.get(1);
        setXDataName(gdX.getParam());
        setYDataName(gdY.getParam());
        
        // TODO set axis scales, update and repaint
    }

}


