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
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;

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
import javax.swing.SwingUtilities;
import javax.swing.WindowConstants;

import one.ben.fcs.AbstractPanel2.AXIS_SCALE;
import one.ben.fcs.AbstractPanel2.PLOT_TYPE;
import one.ben.fcs.FCSDataLoader.DATA_SET;
import cnv.gui.GuiManager;
import common.Files;
import common.Logger;

public class FCSPlot extends JPanel implements WindowListener, ActionListener { 
    
    static final int START_X = 20;
    static final int START_Y = 20;
    static final int START_WIDTH = 1000;
    static final int START_HEIGHT = 600;

    public static final long serialVersionUID = 1L;
	public static final Color BACKGROUND_COLOR = Color.WHITE;
	public static final String[] FILE_EXTENSIONS = {".fcs"};
	private FCSPanel fcsPanel;
	private FCSPlotControlPanel fcsControls;
	private JLayeredPane layeredPane;
	private JFrame controlFrame;

	private List<String> dataKeys;
	HashSet<String> validExts;
	Logger log;

    public FCSPlot() {
        this(null);
    }
    
	private FCSPlot(String[] fileExts) {
		log = new Logger();

		validExts = new HashSet<String>();
		dataKeys = Collections.synchronizedList(new ArrayList<String>());
		
		if (fileExts != null) {
    		for (String ext : fileExts) {
    		    validExts.add(ext);
    		}
		}
		
		setLayout(new BorderLayout());
		
		fcsPanel = new FCSPanel(this);

		layeredPane = new JLayeredPane() {
            private static final long serialVersionUID = 1L;
            @Override
		    public boolean isOptimizedDrawingEnabled() {
		        return false; // override to force proper z-ordered painting
		    }
		};
		layeredPane.setLayout(new BorderLayout());
		
		layeredPane.add(fcsPanel);
		layeredPane.setComponentZOrder(fcsPanel, 0);
		layeredPane.setPreferredSize(new Dimension(1000, 600));
		
		Dimension minimumSize = new Dimension(100, 50);
		layeredPane.setMinimumSize(minimumSize);

		add(layeredPane, BorderLayout.CENTER);

		inputMapAndActionMap();

		fcsPanel.setPointsGeneratable(true);
		fcsPanel.setExtraLayersVisible(new byte[] {99});
		fcsPanel.setChartType(AbstractPanel2.PLOT_TYPE.HEATMAP);
		
		fcsControls = new FCSPlotControlPanel(this);
		controlFrame = new JFrame();
		controlFrame.pack();
		controlFrame.add(fcsControls);
		controlFrame.setDefaultCloseOperation(WindowConstants.HIDE_ON_CLOSE);
		controlFrame.setBounds(START_X + START_WIDTH + 50, START_Y, 250, 400);
		
		updateGUI();
		
		fcsPanel.grabFocus();
		
		setVisible(true);
		controlFrame.setVisible(true);
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
		JMenuItem menuItemExit, menuItemOpen, menuItemRemove, menuItemRemoveAll, menuItemScreens, menuItemSave;
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
//		menu = new JMenu("View");
//        menu.setMnemonic(KeyEvent.VK_V);
//		menuBar.add(menu);
		
		return menuBar;
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
	
	FCSDataLoader dataLoader;
	boolean isLoading = false;
	
	private volatile String xDataName;
	private volatile String yDataName;
	private volatile AXIS_SCALE xScale;
	private volatile AXIS_SCALE yScale;
	private volatile PLOT_TYPE plotType;
	
    public String getXDataName() { return xDataName; }
    public String getYDataName() { return yDataName; }
    
    public AXIS_SCALE getXScale() { return xScale; }
    public AXIS_SCALE getYScale() { return yScale; }
    
    public PLOT_TYPE getPlotType() { return plotType; }

    
    protected void setXDataName(String xDataName) { this.xDataName = xDataName; }
    protected void setYDataName(String yDataName) { this.yDataName = yDataName; }

    protected void setXScale(AXIS_SCALE scale) { this.xScale = scale; }
    protected void setYScale(AXIS_SCALE scale) { this.yScale = scale; }
	
    protected void setPlotType(PLOT_TYPE type) { this.fcsPanel.chartType = type; }
    
	public ArrayList<String[]> getDataSelected() {
		ArrayList<String[]> v = null;
		
		
		
		return v;
	}

	public void updateGUI() {
		fcsPanel.paintAgain();
	}

	public void windowActivated(WindowEvent e) {}

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
    
    public double[] getAxisData(boolean wait, boolean xAxis) {
        double[] data;
        if (dataLoader == null) {
            data = null;
        } else if (!dataLoader.loaded) {
            if (wait) {
                while (!dataLoader.loaded) {
                    Thread.yield();
                }
                data = dataLoader.getData(xAxis ? getXDataName() : getYDataName());
            } else {
                data = new double[0];
            }
        } else {
            data = dataLoader.getData(xAxis ? getXDataName() : getYDataName());
        }
        return data;
    }
    
    private void resetForNewData(FCSDataLoader newDataLoader) {
        ArrayList<String> colNames = newDataLoader.getAllDisplayableNames(DATA_SET.ALL);
        if (colNames.size() < 2) {
            // TODO error, not enough data!!
        }
        // TODO reset GUI elements
        fcsControls.setPlotType(PLOT_TYPE.HEATMAP);
        fcsControls.setColumns(colNames.toArray(new String[colNames.size()]), true, 1);
        fcsControls.setColumns(colNames.toArray(new String[colNames.size()]), false, 0);
        fcsControls.setScale(newDataLoader.scales.get(0), false);
        fcsControls.setScale(newDataLoader.scales.get(1), true);
        
        setYDataName(colNames.get(0));
        setXDataName(colNames.get(1));
        setYScale(newDataLoader.scales.get(0));
        setXScale(newDataLoader.scales.get(1));
    }
    
	public void loadFile(String filename) {
	    if (filename == null || !Files.exists(filename) || (dataLoader != null && dataLoader.loadedFile.equals(filename))) {
	        return;
	    }
	    isLoading = true;
	    Thread dataLoaderThread = new Thread(new Runnable() {
            @Override
            public void run() {
                FCSDataLoader newDataLoader = new FCSDataLoader();
                try {
                    newDataLoader.loadData(filename);
                } catch (IOException e) {
                    e.printStackTrace(); // TODO
                    return;
                }
                resetForNewData(newDataLoader);
                dataLoader = newDataLoader;
                isLoading = false;
                System.gc();
                updateGUI();
            }
        });
	    dataLoaderThread.start();
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
	public static FCSPlot createGUI(/*Project proj, */boolean show/*, boolean promptOnClose, StringListProperty filenamesProperty, String... fileExts*/) {
		JFrame frame = new JFrame("Genvisis - FCS Plot");
		frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

        //Create and set up the content pane.
        FCSPlot twoDPlot = new FCSPlot(/*proj, promptOnClose, fileExts, filenamesProperty*/);
        frame.setJMenuBar(twoDPlot.menuBar());
        twoDPlot.setOpaque(true); //content panes must be opaque
        frame.setContentPane(twoDPlot);
        frame.addWindowListener(twoDPlot);
        frame.setBounds(START_X, START_Y, START_WIDTH, START_HEIGHT);
//		frame.setExtendedState(frame.getExtendedState()|JFrame.MAXIMIZED_BOTH);

        //Display the window.
        frame.pack();
        frame.setVisible(show);

        String fcsFilename = "F:\\Flow\\P1-B&C-CD3-APC-Cy7 or CD4-APC-Cy7_ULTRA BRIGHT RAINBOW BEADS_URB_001.fcs";
        twoDPlot.loadFile(fcsFilename);
        
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
        return dataLoader == null ? 0 : dataLoader.loaded ? dataLoader.eventCount : 0;
    }

}


