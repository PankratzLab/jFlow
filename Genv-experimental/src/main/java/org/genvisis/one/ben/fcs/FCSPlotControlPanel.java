package org.genvisis.one.ben.fcs;

import java.awt.Dimension;
import java.awt.Font;
import java.awt.Insets;
import java.awt.LayoutManager;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.lang.reflect.InvocationTargetException;
import java.text.Format;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.TreeSet;

import javax.swing.AbstractAction;
import javax.swing.ActionMap;
import javax.swing.BorderFactory;
import javax.swing.ButtonGroup;
import javax.swing.DefaultComboBoxModel;
import javax.swing.InputMap;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JFormattedTextField;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.KeyStroke;
import javax.swing.ScrollPaneConstants;
import javax.swing.Scrollable;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;
import javax.swing.border.BevelBorder;
import javax.swing.filechooser.FileNameExtensionFilter;

import org.genvisis.cnv.gui.JAccordionPanel;
import org.genvisis.common.Files;
import org.genvisis.common.ext;
import org.genvisis.one.ben.fcs.AbstractPanel2.AXIS_SCALE;
import org.genvisis.one.ben.fcs.AbstractPanel2.PLOT_TYPE;
import org.genvisis.one.ben.fcs.FCSDataLoader.LOAD_STATE;
import org.genvisis.one.ben.fcs.FCSPanel.GATING_TOOL;

import net.miginfocom.swing.MigLayout;

public class FCSPlotControlPanel extends JPanel {

	/**
	* 
	*/
	private static final long serialVersionUID = 1L;

	private FCSPlot plot;

	private JComboBox<PLOT_TYPE> cbType;
	private JComboBox<String> cbYData;
	private JComboBox<String> cbXData;
	private JComboBox<AXIS_SCALE> cbYScale;
	private JComboBox<AXIS_SCALE> cbXScale;

	private final Font lblFont = new Font("Arial", 0, 12);
	private JFormattedTextField yBndsMin;
	private JFormattedTextField yBndsMax;
	private JFormattedTextField xBndsMin;
	private JFormattedTextField xBndsMax;

	volatile boolean progSet = false;

	private JCheckBox chckbxShowMedianX;
	private JCheckBox chckbxShowSdX;
	private JCheckBox chckbxShowMedianY;
	private JCheckBox chckbxShowSdY;

	public JProgressBar progressBar;

	private JAccordionPanel plotControlPanel;
	private JAccordionPanel dataControlsPanel;
	private JAccordionPanel gateControlPanel;
	private JPanel panel_1;
	private JLabel gateFileLabel;
	private JLabel gateFileTitle;
	private JSeparator gateFileSep;

	private JButton btnSaveGating;
	private JButton gateSelectBtn;
	private JButton gateClearBtn;
	private JButton dirSelectBtn;
	private String prevGateDir;
	private String prevFCSDir;

	/**
	 * Create the panel.
	 */
	public FCSPlotControlPanel(final FCSPlot plot) {
		this.plot = plot;

		setLayout(new MigLayout("ins 0", "[grow]", "[grow][]"));

		panel_1 = new JPanel();
		panel_1.setLayout(new MigLayout("ins 0", "[grow]", "[]0[]0[]0[grow]"));

		ArrayList<JAccordionPanel> bg = new ArrayList<JAccordionPanel>();

		plotControlPanel = new JAccordionPanel();
		panel_1.add(plotControlPanel, "cell 0 0,grow");
		JLabel ctrlLabel = new JLabel("<html><u>Plot Controls</u></html>");
		ctrlLabel.setFont(ctrlLabel.getFont().deriveFont(Font.PLAIN, 14));
		plotControlPanel.topPanel.add(ctrlLabel, "pad 0 10 0 0, cell 0 0, grow");
		JLabel mnemLabel = new JLabel("(alt - P)");
		mnemLabel.setFont(mnemLabel.getFont().deriveFont(Font.PLAIN, 9));
		plotControlPanel.topPanel.add(mnemLabel, "cell 1 0, alignx right");
		plotControlPanel.topPanel.setBorder(BorderFactory.createBevelBorder(BevelBorder.RAISED));
		plotControlPanel.contentPanel.setBorder(BorderFactory.createBevelBorder(BevelBorder.LOWERED));
		plotControlPanel.addToGroup(bg);

		JPanel panel = plotControlPanel.contentPanel;

		panel.setLayout(new MigLayout("", "[][][grow][]", "[][][][][][][][][][][]"));

		JLabel lblPlotType = new JLabel("Plot Type:");
		panel.add(lblPlotType, "cell 0 0,alignx trailing");
		lblPlotType.setFont(lblFont);

		cbType = new JComboBox<PLOT_TYPE>(PLOT_TYPE.values());
		panel.add(cbType, "cell 1 0 3 1,growx");

		JLabel lblYaxisData = new JLabel("Y-Axis:");
		panel.add(lblYaxisData, "cell 0 2,alignx trailing");
		lblYaxisData.setFont(lblFont);

		cbYData = new JComboBox<String>();
		panel.add(cbYData, "cell 1 2 3 1,growx");
		cbYData.addItemListener(new ItemListener() {
			@Override
			public void itemStateChanged(ItemEvent arg0) {
				if (arg0.getStateChange() == ItemEvent.SELECTED) {
					if (FCSPlot.HISTOGRAM_COL.equals(arg0.getItem().toString())) {
						plot.setPlotType(PLOT_TYPE.HISTOGRAM);
					} else {
						if (plot.getPlotType() == PLOT_TYPE.HISTOGRAM) {
							plot.setPlotType(PLOT_TYPE.DOT_PLOT);
						}
						plot.setYDataName(arg0.getItem().toString());
					}
					plot.updateGUI();
				}
			}
		});
		cbYData.setMaximumRowCount(15);

		JLabel lblScale = new JLabel("Scale:");
		panel.add(lblScale, "cell 0 3 2 1,alignx trailing");
		lblScale.setFont(lblFont);

		cbYScale = new JComboBox<AbstractPanel2.AXIS_SCALE>(AXIS_SCALE.values());
		panel.add(cbYScale, "cell 2 3 2 1,growx");
		cbYScale.addItemListener(new ItemListener() {
			@Override
			public void itemStateChanged(ItemEvent arg0) {
				if (arg0.getStateChange() == ItemEvent.SELECTED) {
					plot.setYScale((AXIS_SCALE) arg0.getItem());
					plot.updateGUI();
				}
			}
		});

		chckbxShowMedianY = new JCheckBox("Show Median", plot.showMedian(true));
		panel.add(chckbxShowMedianY, "cell 0 4 3 1,alignx trailing");
		chckbxShowMedianY.setHorizontalAlignment(SwingConstants.TRAILING);

		chckbxShowSdY = new JCheckBox("Show SD", plot.showSD(true));
		panel.add(chckbxShowSdY, "cell 3 4,alignx trailing");
		chckbxShowSdY.setHorizontalAlignment(SwingConstants.TRAILING);

		JLabel lblYbounds = new JLabel("Y-Bounds:");
		panel.add(lblYbounds, "cell 0 5 2 1,alignx trailing");
		lblYbounds.setFont(lblFont);

		Format numberFormat = NumberFormat.getNumberInstance();

		yBndsMin = new JFormattedTextField(numberFormat);
		panel.add(yBndsMin, "cell 2 5 2 1");
		yBndsMin.addPropertyChangeListener("value", pcl);
		yBndsMin.setColumns(10);
		yBndsMin.setValue(0);
		yBndsMin.setEditable(false);

		yBndsMax = new JFormattedTextField(numberFormat);
		panel.add(yBndsMax, "cell 2 5 2 1");
		yBndsMax.addPropertyChangeListener("value", pcl);
		yBndsMax.setColumns(10);
		yBndsMax.setEditable(false);

		JLabel lblXaxisData = new JLabel("X-Axis:");
		panel.add(lblXaxisData, "cell 0 7,alignx trailing");
		lblXaxisData.setFont(lblFont);

		cbXData = new JComboBox<String>();
		panel.add(cbXData, "cell 1 7 3 1,growx");
		cbXData.addItemListener(new ItemListener() {
			@Override
			public void itemStateChanged(ItemEvent arg0) {
				if (arg0.getStateChange() == ItemEvent.SELECTED) {
					plot.setXDataName(arg0.getItem().toString());
					plot.updateGUI();
				}
			}
		});
		cbXData.setMaximumRowCount(15);

		JLabel lblScale_1 = new JLabel("Scale:");
		panel.add(lblScale_1, "cell 0 8 2 1,alignx trailing");
		lblScale_1.setFont(lblFont);

		cbXScale = new JComboBox<AbstractPanel2.AXIS_SCALE>(AXIS_SCALE.values());
		panel.add(cbXScale, "cell 2 8 2 1,growx");

		chckbxShowMedianX = new JCheckBox("Show Median", plot.showMedian(false));
		panel.add(chckbxShowMedianX, "cell 0 9 3 1,alignx right");
		chckbxShowMedianX.setHorizontalAlignment(SwingConstants.TRAILING);

		chckbxShowSdX = new JCheckBox("Show SD", plot.showSD(false));
		panel.add(chckbxShowSdX, "cell 3 9,alignx right");
		chckbxShowSdX.setHorizontalAlignment(SwingConstants.TRAILING);

		JLabel lblXbounds = new JLabel("X-Bounds:");
		panel.add(lblXbounds, "cell 0 10 2 1,alignx trailing");
		lblXbounds.setFont(lblFont);

		xBndsMin = new JFormattedTextField(numberFormat);
		panel.add(xBndsMin, "cell 2 10 2 1");
		xBndsMin.addPropertyChangeListener("value", pcl);
		xBndsMin.setColumns(10);
		xBndsMin.setValue(0);
		xBndsMin.setEditable(false);

		xBndsMax = new JFormattedTextField(numberFormat);
		panel.add(xBndsMax, "cell 2 10 2 1");
		xBndsMax.addPropertyChangeListener("value", pcl);
		xBndsMax.setColumns(10);
		xBndsMax.setEditable(false);

		chckbxShowSdX.addItemListener(new ItemListener() {
			@Override
			public void itemStateChanged(ItemEvent arg0) {
				boolean show = arg0.getStateChange() == ItemEvent.SELECTED;
				plot.setSDVisible(show, false);
				plot.updateGUI();
			}
		});
		chckbxShowMedianX.addItemListener(new ItemListener() {
			@Override
			public void itemStateChanged(ItemEvent arg0) {
				boolean show = arg0.getStateChange() == ItemEvent.SELECTED;
				plot.setMedianVisible(show, false);
				plot.updateGUI();
			}
		});
		cbXScale.addItemListener(new ItemListener() {
			@Override
			public void itemStateChanged(ItemEvent arg0) {
				if (arg0.getStateChange() == ItemEvent.SELECTED) {
					plot.setXScale((AXIS_SCALE) arg0.getItem());
					plot.updateGUI();
				}
			}
		});
		chckbxShowSdY.addItemListener(new ItemListener() {
			@Override
			public void itemStateChanged(ItemEvent arg0) {
				boolean show = arg0.getStateChange() == ItemEvent.SELECTED;
				plot.setSDVisible(show, true);
				plot.updateGUI();
			}
		});
		chckbxShowMedianY.addItemListener(new ItemListener() {
			@Override
			public void itemStateChanged(ItemEvent arg0) {
				boolean show = arg0.getStateChange() == ItemEvent.SELECTED;
				plot.setMedianVisible(show, true);
				plot.updateGUI();
			}
		});
		cbType.addItemListener(new ItemListener() {
			@Override
			public void itemStateChanged(ItemEvent arg0) {
				if (arg0.getStateChange() == ItemEvent.SELECTED) {
					plot.setPlotType((PLOT_TYPE) arg0.getItem());
					plot.updateGUI();
				}
			}
		});

		gateControlPanel = new JAccordionPanel();
		panel_1.add(gateControlPanel, "cell 0 1,grow");
		gateControlPanel.contentPanel.setBorder(BorderFactory.createBevelBorder(BevelBorder.LOWERED));
		gateControlPanel.topPanel.setBorder(BorderFactory.createBevelBorder(BevelBorder.RAISED));
		gateControlPanel.addToGroup(bg);
		JLabel gCtrlLabel = new JLabel("<html><u>Gate Controls</u></html>");
		gCtrlLabel.setFont(gCtrlLabel.getFont().deriveFont(Font.PLAIN, 14));
		gateControlPanel.topPanel.add(gCtrlLabel, "pad 0 10 0 0, cell 0 0, grow");
		mnemLabel = new JLabel("(alt - G)");
		mnemLabel.setFont(mnemLabel.getFont().deriveFont(Font.PLAIN, 9));
		gateControlPanel.topPanel.add(mnemLabel, "cell 1 0, alignx right");

		JPanel gatePanel = gateControlPanel.contentPanel;
		gatePanel.setLayout(new MigLayout("hidemode 3,ins 0", "[grow][]0px[]",
																			"0px[]0px[]0px[grow]0px[]0px[]0px"));

		JLabel dirLbl1 = new JLabel("Select Gating File:");
		dirLbl1.setHorizontalAlignment(SwingConstants.RIGHT);
		gatePanel.add(dirLbl1, "cell 0 0, growx, alignx right");

		gateSelectBtn = new JButton(">>");
		gateSelectBtn.addActionListener(gateSelectListener);
		gateSelectBtn.setMargin(new Insets(0, 2, 0, 2));
		gatePanel.add(gateSelectBtn, "cell 1 0");

		gateClearBtn = new JButton("X");
		gateClearBtn.addActionListener(gateClearListener);
		gateClearBtn.setMargin(new Insets(0, 2, 0, 2));
		gatePanel.add(gateClearBtn, "cell 2 0");

		panel = new JPanel(new MigLayout("hidemode 3", "[grow][]", ""));
		gatePanel.add(panel, "cell 0 1, span 3, grow");

		panel.add(new JSeparator(SwingConstants.HORIZONTAL), "cell 0 0 2 1, growx");

		gateFileTitle = new JLabel("Current Gate File:");
		gateFileTitle.setFont(new Font("Arial", Font.BOLD, 9));
		panel.add(gateFileTitle, "cell 0 1");
		gateFileTitle.setVisible(false);
		gateFileLabel = new JLabel();
		gateFileLabel.setFont(new Font("Arial", Font.PLAIN, 9));
		panel.add(gateFileLabel, "cell 0 2 2 1, grow");
		gateFileLabel.setVisible(false);

		gateFileSep = new JSeparator(SwingConstants.HORIZONTAL);
		panel.add(gateFileSep, "cell 0 3 2 1, growx");
		gateFileSep.setVisible(false);

		JLabel gateTypeLbl = new JLabel("Gate Tool:");
		panel.add(gateTypeLbl, "cell 0 4, span 2 1");
		String[] gateTypes = new String[GATING_TOOL.values().length];
		for (int i = 0; i < GATING_TOOL.values().length; i++) {
			gateTypes[i] = GATING_TOOL.values()[i].getDisplayName();
		}
		JComboBox<String> gateTypeCmb = new JComboBox<String>(gateTypes);
		gateTypeCmb.addItemListener(new ItemListener() {
			@Override
			public void itemStateChanged(ItemEvent e) {
				if (e.getStateChange() == ItemEvent.SELECTED) {
					String displayName = (String) e.getItem();
					plot.setGateTool(GATING_TOOL.getGatingToolByDisplayName(displayName));
				}
			}
		});
		panel.add(gateTypeCmb, "cell 0 4, growx");

		JCheckBox chkDrawPolysBinned = new JCheckBox("Bin Polygons (FlowJo)");
		chkDrawPolysBinned.addItemListener(new ItemListener() {
			@Override
			public void itemStateChanged(ItemEvent arg0) {
				plot.setDrawPolysAsFlowJo(arg0.getStateChange() == ItemEvent.SELECTED);
				// TODO display warning that existing gates won't be altered?
				// TODO or rather, display confirm dialog asking if all gates should be changed?
			}
		});
		chkDrawPolysBinned.setHorizontalAlignment(SwingConstants.CENTER);
		panel.add(chkDrawPolysBinned, "cell 0 5, grow");

		JRadioButton chkRegGate = new JRadioButton("Standard", true);
		chkRegGate.addItemListener(new ItemListener() {
		  @Override
		  public void itemStateChanged(ItemEvent arg0) {
		    plot.setLeafgating(arg0.getStateChange() != ItemEvent.SELECTED);
		    plot.setBackgating(arg0.getStateChange() != ItemEvent.SELECTED);
		  }
		});
		chkRegGate.setHorizontalAlignment(SwingConstants.CENTER);
		panel.add(chkRegGate, "cell 0 6, grow");
		
		JRadioButton chkLeafgate = new JRadioButton("Leafgate");
		chkLeafgate.addItemListener(new ItemListener() {
		  @Override
		  public void itemStateChanged(ItemEvent arg0) {
		    plot.setLeafgating(arg0.getStateChange() == ItemEvent.SELECTED);
		    plot.setBackgating(arg0.getStateChange() != ItemEvent.SELECTED);
		  }
		});
		chkLeafgate.setHorizontalAlignment(SwingConstants.CENTER);
		panel.add(chkLeafgate, "cell 0 7, grow");

		JRadioButton chkBackgate = new JRadioButton("Backgate");
		chkBackgate.addItemListener(new ItemListener() {
			@Override
			public void itemStateChanged(ItemEvent arg0) {
				plot.setBackgating(arg0.getStateChange() == ItemEvent.SELECTED);
	            plot.setLeafgating(arg0.getStateChange() != ItemEvent.SELECTED);
			}
		});
		chkBackgate.setHorizontalAlignment(SwingConstants.CENTER);
		panel.add(chkBackgate, "cell 0 8, grow");
		
		ButtonGroup bg1 = new ButtonGroup();
		bg1.add(chkRegGate);
		bg1.add(chkLeafgate);
		bg1.add(chkBackgate);

		btnSaveGating = new JButton("Save Gating");
		btnSaveGating.addActionListener(gateSaveListener);
		btnSaveGating.setHorizontalAlignment(SwingConstants.CENTER);
		panel.add(btnSaveGating, "cell 0 9, span 2, center");

		dataControlsPanel = new JAccordionPanel();
		panel_1.add(dataControlsPanel, "cell 0 2,grow");
		dataControlsPanel.contentPanel.setBorder(BorderFactory.createBevelBorder(BevelBorder.LOWERED));
		dataControlsPanel.topPanel.setBorder(BorderFactory.createBevelBorder(BevelBorder.RAISED));
		dataControlsPanel.addToGroup(bg);
		JLabel dCtrlLabel = new JLabel("<html><u>Data Controls</u></html>");
		dCtrlLabel.setFont(dCtrlLabel.getFont().deriveFont(Font.PLAIN, 14));
		dataControlsPanel.topPanel.add(dCtrlLabel, "pad 0 10 0 0, cell 0 0, grow");
		mnemLabel = new JLabel("(alt - D)");
		mnemLabel.setFont(mnemLabel.getFont().deriveFont(Font.PLAIN, 9));
		dataControlsPanel.topPanel.add(mnemLabel, "cell 1 0, alignx right");

		JPanel dataPanel = dataControlsPanel.contentPanel;
		dataPanel.setLayout(new MigLayout("hidemode 3,ins 0", "[grow][]",
																			"0px[]0px[]0px[grow]0px[]0px[]0px"));

		JLabel dirLbl = new JLabel("Add File(s):");
		dirLbl.setHorizontalAlignment(SwingConstants.RIGHT);
		dataPanel.add(dirLbl, "cell 0 0, growx, alignx right");

		// fileDirField = new JTextField();
		// dataPanel.add(fileDirField, "cell 0 1, growx, split 2");
		dirSelectBtn = new JButton(">>");
		dirSelectBtn.addActionListener(dirSelectListener);
		dirSelectBtn.setMargin(new Insets(0, 2, 0, 2));
		dataPanel.add(dirSelectBtn, "cell 1 0");

		dataScrollPane = new JScrollPane();
		actualDataPanel = new ScrollablePanel(new MigLayout("", "", ""));
		dataScrollPane.setViewportView(actualDataPanel);
		dataScrollPane.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
		dataScrollPane.setVerticalScrollBarPolicy(ScrollPaneConstants.VERTICAL_SCROLLBAR_ALWAYS);
		dataPanel.add(dataScrollPane, "cell 0 1 2 1, grow");

		InputMap im = dataScrollPane.getInputMap(WHEN_ANCESTOR_OF_FOCUSED_COMPONENT);
		ActionMap am = dataScrollPane.getActionMap();

		im.put(KeyStroke.getKeyStroke(KeyEvent.VK_UP, 0), "up");
		im.put(KeyStroke.getKeyStroke(KeyEvent.VK_DOWN, 0), "down");
		im.put(KeyStroke.getKeyStroke(KeyEvent.VK_UP, KeyEvent.CTRL_DOWN_MASK), "up_ctrl");
		im.put(KeyStroke.getKeyStroke(KeyEvent.VK_DOWN, KeyEvent.CTRL_DOWN_MASK), "down_ctrl");
		im.put(KeyStroke.getKeyStroke(KeyEvent.VK_UP, KeyEvent.ALT_DOWN_MASK), "up_ctrl");
		im.put(KeyStroke.getKeyStroke(KeyEvent.VK_DOWN, KeyEvent.ALT_DOWN_MASK), "down_ctrl");

		am.put("up", new AbstractAction() {
			/**
			* 
			*/
			private static final long serialVersionUID = 1L;

			@Override
			public void actionPerformed(ActionEvent e) {
				int sel = -1;
				for (int i = 0; i < filePanels.size(); i++) {
					if (filePanels.get(i).isSelected()) {
						sel = i;
						break;
					}
				}
				if (sel == -1 || sel == 0) {
					return;
				}
				filePanels.get(sel).setSelected(false);
				filePanels.get(sel - 1).setSelected(true);
				dataScrollPane.getViewport().scrollRectToVisible(filePanels.get(sel - 1).getBounds());
			}
		});
		am.put("down", new AbstractAction() {
			/**
			* 
			*/
			private static final long serialVersionUID = 1L;

			@Override
			public void actionPerformed(ActionEvent e) {
				int sel = -1;
				for (int i = 0; i < filePanels.size(); i++) {
					if (filePanels.get(i).isSelected()) {
						sel = i;
						break;
					}
				}
				if (sel == filePanels.size() - 1) {
					return;
				}
				if (sel > -1) {
					filePanels.get(sel).setSelected(false);
				}
				filePanels.get(sel + 1).setSelected(true);
				dataScrollPane.getViewport().scrollRectToVisible(filePanels.get(sel + 1).getBounds());
			}
		});
		// am.put("up_ctrl", arg1);
		// am.put("down_ctrl", arg1);

		// TODO mouse keys (up/down for changing selection, alt/ctrl-up/down for displaying prev/next
		// file data, etc)
		// TODO listener for move up
		// TODO listener for move down
		// TODO listener for information
		// TODO listener for delete
		// TODO listener for load
		// TODO memory management and warnings

		dataMsgPanel = new JPanel();
		dataPanel.add(dataMsgPanel, "cell 0 3, grow");
		dataMsgPanel.setVisible(false);

		add(panel_1, "cell 0 0,grow");

		progressBar = new JProgressBar();
		add(progressBar, "cell 0 1,growx, pad -3 3 -3 -3");

		plotControlPanel.shrink();
		gateControlPanel.shrink();
		dataControlsPanel.expand();
	}

	public class ScrollablePanel extends JPanel implements Scrollable {

		/**
		* 
		*/
		private static final long serialVersionUID = 1L;

		public ScrollablePanel(LayoutManager layout) {
			super(layout);
		}

		@Override
		public Dimension getPreferredScrollableViewportSize() {
			return getPreferredSize();
		}

		@Override
		public int getScrollableUnitIncrement(Rectangle visibleRect, int orientation, int direction) {
			return 10;
		}

		@Override
		public int getScrollableBlockIncrement(Rectangle visibleRect, int orientation, int direction) {
			return ((orientation == SwingConstants.VERTICAL) ? visibleRect.height : visibleRect.width)
							- 10;
		}

		@Override
		public boolean getScrollableTracksViewportWidth() {
			return true;
		}

		@Override
		public boolean getScrollableTracksViewportHeight() {
			return false;
		}
	}

	public void showPlotControls() {
		plotControlPanel.expand();
	}

	public void showGateControls() {
		gateControlPanel.expand();
	}

	public void showDataControls() {
		dataControlsPanel.expand();
	}


	ArrayList<DataControlPanel> filePanels = new ArrayList<DataControlPanel>();

	protected void addFCSFiles(String[] files) {
		TreeSet<String> fileSet = new TreeSet<String>();
		for (String f : files) {
			fileSet.add(f.trim());
		}
		prevFCSDir = ext.parseDirectoryOfFile(files[0].trim());

		for (String f : fileSet) {
			String sz = Files.getSizeScaledString(f, false);
			String dt = "";
			final DataControlPanel dcp = new DataControlPanel(f, sz, dt, dataListener);
			dcp.addMouseListener(new MouseAdapter() {
				@Override
				public void mouseClicked(MouseEvent e) {
					super.mouseClicked(e);
					for (int i = 0; i < filePanels.size(); i++) {
						if (filePanels.get(i) == dcp) {
							filePanels.get(i).setSelected(true);
							dcp.requestFocusInWindow();
							dataScrollPane.getViewport().scrollRectToVisible(dcp.getBounds());
						} else {
							filePanels.get(i).setSelected(false);
						}
					}
				}
			});
			filePanels.add(dcp);
		}

		addFilePanels();
		plot.saveProps();
	}

	protected ArrayList<String> getAddedFiles() {
		ArrayList<String> files = new ArrayList<String>();
		for (DataControlPanel dcp : filePanels) {
			files.add(dcp.file);
		}
		return files;
	}

	private void addFilePanels() {
		SwingUtilities.invokeLater(new Runnable() {
			@Override
			public void run() {
				int index = 0;
				actualDataPanel.removeAll();
				actualDataPanel.add(new JSeparator(JSeparator.HORIZONTAL), "grow, cell 0 " + (index++));
				for (DataControlPanel dcp : filePanels) {
					actualDataPanel.add(dcp, "cell 0 " + (index++));
					actualDataPanel.add(new JSeparator(JSeparator.HORIZONTAL), "grow, cell 0 " + (index++));
				}
				actualDataPanel.revalidate();
				dataControlsPanel.contentPanel.revalidate();
				dataScrollPane.revalidate();

				revalidate();
				repaint();
			}
		});
	}

	ActionListener dataListener = new ActionListener() {
		@Override
		public void actionPerformed(ActionEvent e) {
			DataControlPanel dcp = (DataControlPanel) e.getSource();
			String cmd = e.getActionCommand();
			// String file = cmd.split("::")[0];
			cmd = cmd.split("::")[1];
			int ind = -1;
			for (int i = 0; i < filePanels.size(); i++) {
				if (filePanels.get(i) == dcp) {
					ind = i;
				} else {
					filePanels.get(i).setSelected(false);
				}
			}
			if (ind == -1) {
				// TODO error!
				System.err.println("Error - couldn't find DataControlPanel!");
				return;
			}
			dcp.setSelected(true);
			if (cmd.equals(DataControlPanel.ACTION_DELETE)) {
				boolean remove = true;
				if (plot.isFileLoaded(dcp.file)) {
					int opt = JOptionPane.showConfirmDialog(plot,
																									"This will unload file data from memory - are you sure you wish to continue?",
																									"Unload Data?", JOptionPane.YES_NO_OPTION,
																									JOptionPane.WARNING_MESSAGE);
					if (opt == JOptionPane.YES_OPTION) {
						plot.unloadFile(dcp.file);
					} else {
						remove = false;
					}
				}
				if (remove) {
					filePanels.remove(ind);
					addFilePanels();
					plot.saveProps();
				}
			} else if (cmd.equals(DataControlPanel.ACTION_MOVE_UP)) {
				if (ind == 0) {
					return;
				}
				Collections.swap(filePanels, ind, ind - 1);
				addFilePanels();
				plot.saveProps();
			} else if (cmd.equals(DataControlPanel.ACTION_MOVE_DOWN)) {
				if (ind == filePanels.size() - 1) {
					return;
				}
				Collections.swap(filePanels, ind, ind + 1);
				addFilePanels();
				plot.saveProps();
			} else if (cmd.equals(DataControlPanel.ACTION_INFO)) {
				// TODO build info GUI for files
			} else if (cmd.equals(DataControlPanel.ACTION_LOAD)) {
				// TODO check memory available before loading, show YELLOW warning if nearing max, show RED
				// warning if not enough memory
				if (plot.isFileLoaded(dcp.file)) {
					int opt = JOptionPane.showConfirmDialog(plot,
																									"This will unload file data from memory - are you sure you wish to continue?",
																									"Unload Data?", JOptionPane.YES_NO_OPTION,
																									JOptionPane.WARNING_MESSAGE);
					if (opt == JOptionPane.YES_OPTION) {
						plot.unloadFile(dcp.file);
						dcp.setLoaded(false);
					}
				} else {
					dcp.setLoaded(true);
					plot.loadFile(dcp.file, false);
				}
			} else if (cmd.equals(DataControlPanel.ACTION_USE)) {
				boolean disp = true;
				if (!plot.isFileLoaded(dcp.file)) {
					int opt = JOptionPane.showConfirmDialog(plot,
																									"Data file not loaded - would you like to load this data?",
																									"Confirm Load", JOptionPane.YES_NO_OPTION,
																									JOptionPane.INFORMATION_MESSAGE);
					if (opt == JOptionPane.NO_OPTION) {
						disp = false;
					}
				}
				if (disp) {
					dcp.setLoaded(true);
					plot.loadFile(dcp.file, true);
				}
			}
			// apply action to file
		}
	};

	ActionListener gateSaveListener = new ActionListener() {
		@Override
		public void actionPerformed(ActionEvent e) {
			plot.saveGating();
		}
	};

	ActionListener gateSelectListener = new ActionListener() {
		@Override
		public void actionPerformed(ActionEvent e) {
			JFileChooser jfc = new JFileChooser(prevGateDir);
			jfc.setFileSelectionMode(JFileChooser.FILES_ONLY);
			jfc.setMultiSelectionEnabled(false);
			jfc.addChoosableFileFilter(new FileNameExtensionFilter("Gating-ML File", "xml"));
			jfc.addChoosableFileFilter(new FileNameExtensionFilter(	"FlowJo Workspace or WorkspaceTemplate File",
																															"wsp", "wspt"));
			jfc.setDialogTitle("Select Gating File");
			int resp = jfc.showOpenDialog(FCSPlotControlPanel.this);
			if (resp == JFileChooser.APPROVE_OPTION) {
				File newFile = jfc.getSelectedFile();
				loadGatingFile(newFile.getAbsolutePath());
			}
		}
	};

	ActionListener gateClearListener = new ActionListener() {
		@Override
		public void actionPerformed(ActionEvent e) {
			plot.clearGating();
			gateFileLabel.setVisible(false);
			gateFileTitle.setVisible(false);
			gateFileSep.setVisible(false);
			plot.updateGUI();
		}
	};

	protected void loadGatingFile(String newFile) {
		prevGateDir = ext.parseDirectoryOfFile(newFile);
		plot.loadWorkspaceFile(newFile);
		plot.refreshGating();
		gateFileLabel.setText("<html><p>" + newFile + "</p></html>");
		gateFileLabel.setVisible(true);
		gateFileTitle.setVisible(true);
		gateFileSep.setVisible(true);
	}

	ActionListener dirSelectListener = new ActionListener() {
		@Override
		public void actionPerformed(ActionEvent e) {
			JFileChooser jfc = new JFileChooser(prevFCSDir);
			jfc.setFileSelectionMode(JFileChooser.FILES_ONLY);
			jfc.setMultiSelectionEnabled(true);
			jfc.setFileFilter(new FileNameExtensionFilter("FCS Files", "fcs"));
			jfc.setDialogTitle("Select FCS File(s)");
			int resp = jfc.showOpenDialog(FCSPlotControlPanel.this);
			if (resp == JFileChooser.APPROVE_OPTION) {
				File[] newFiles = jfc.getSelectedFiles();
				String[] f = new String[newFiles.length];
				for (int i = 0; i < f.length; i++) {
					f[i] = newFiles[i].getAbsolutePath();
				}
				addFCSFiles(f);
			}
		}
	};


	PropertyChangeListener pcl = new PropertyChangeListener() {
		@Override
		public void propertyChange(PropertyChangeEvent evt) {
			if (evt == null || plot == null || progSet) {
				return;
			}
			JFormattedTextField jftf = (JFormattedTextField) evt.getSource();
			String prop = evt.getPropertyName();
			if (jftf == xBndsMin) {
				prop = AbstractPanel2.X_MIN;
			} else if (jftf == xBndsMax) {
				prop = AbstractPanel2.X_MAX;
			} else if (jftf == yBndsMin) {
				prop = AbstractPanel2.Y_MIN;
			} else if (jftf == yBndsMax) {
				prop = AbstractPanel2.Y_MAX;
			}
			Object oldV = evt.getOldValue();
			Object newV = evt.getNewValue();
			Double oldV2 = oldV == null ? null : ((Number) oldV).doubleValue();
			Double newV2 = newV == null ? null : ((Number) newV).doubleValue();
			/* FCSPlotControlPanel.this. */firePropertyChange(prop, oldV2, newV2);
			// plot.firePropertyChange(prop, oldV2, newV2);
			// PropertyChangeEvent pce = new PropertyChangeEvent(FCSPlotControlPanel.this, prop, oldV2,
			// newV2);
			// plot.propertyChange(pce);

		}
	};

	private JPanel dataMsgPanel;

	private JPanel actualDataPanel;

	private JScrollPane dataScrollPane;

	public void setPlotType(PLOT_TYPE typ) {
		cbType.setSelectedItem(typ);
	}

	public PLOT_TYPE getPlotType() {
		return (PLOT_TYPE) cbType.getSelectedItem();
	}

	public void setScale(AXIS_SCALE scl, boolean x) {
		(x ? cbXScale : cbYScale).setSelectedItem(scl);
	}

	public void setColumns(String[] dataNames, boolean x, int selected) {
		(x ? cbXData : cbYData).setModel(new DefaultComboBoxModel<String>(dataNames));
		(x ? cbXData : cbYData).setSelectedIndex(selected);
		(x ? cbXData : cbYData).repaint();
	}

	public String getSelectedX() {
		return (String) cbXData.getSelectedItem();
	}

	public String getSelectedY() {
		return (String) cbYData.getSelectedItem();
	}

	public void setYData(int index) {
		cbYData.setSelectedIndex(index);
		cbYData.repaint();
	}

	public void setYData(String name) {
		cbYData.setSelectedItem(name);
		cbYData.repaint();
	}

	public void setXData(String name) {
		cbXData.setSelectedItem(name);
		cbXData.repaint();
	}

	private void resetProgSet() {
		SwingUtilities.invokeLater(new Runnable() {
			@Override
			public void run() {
				progSet = false;
			}
		});
	}

	public void setXMin(double xMin) {
		progSet = true;
		xBndsMin.setValue(Math.floor(xMin));
		resetProgSet();
	}

	public void setXMax(double xMax) {
		progSet = true;
		xBndsMax.setValue(Math.ceil(xMax));
		resetProgSet();
	}

	public void setYMin(double yMin) {
		progSet = true;
		yBndsMin.setValue(Math.floor(yMin));
		resetProgSet();
	}

	public void setYMax(double yMax) {
		progSet = true;
		yBndsMax.setValue(Math.ceil(yMax));
		resetProgSet();
	}

	// public void startFileLoading(FCSDataLoader newDataLoader) {
	// // TODO Auto-generated method stub
	//
	// }

	public void startFileLoading(final FCSDataLoader newDataLoader) {
		new Thread(new Runnable() {
			@Override
			public void run() {
				LOAD_STATE state = null;
				while ((state = newDataLoader.getLoadState()) != LOAD_STATE.LOADED) {
					final LOAD_STATE finalState = state;
					try {
						SwingUtilities.invokeAndWait(new Runnable() {
							@Override
							public void run() {
								switch (finalState) {
									case LOADED:
										progressBar.setStringPainted(false);
										progressBar.setString(null);
										progressBar.setIndeterminate(false);
										// hide or set to complete or reset
										break;
									case LOADING:
										progressBar.setIndeterminate(true);
										progressBar.setStringPainted(true);
										progressBar.setString("Loading File...");
										progressBar.setVisible(true);
										// set to indeterminate
										break;
									case PARTIALLY_LOADED:
									case LOADING_REMAINDER:
										progressBar.setIndeterminate(false);
										progressBar.setStringPainted(true);
										int[] stat = newDataLoader.getLoadedStatus();
										// System.out.println(stat[0] + " - " + stat[1]);
										progressBar.setMinimum(0);
										progressBar.setValue(stat[0]);
										progressBar.setMaximum(stat[1]);
										progressBar.setString(null);
										// progressBar.setString("Loading File: " + stat[0] + "/" + stat[1]);
										progressBar.setVisible(true);
										// set to determinate, wait for updates
										break;
									case UNLOADED:
										// what??
										break;
									default:
										// what??
										break;
								}
							}
						});
					} catch (InvocationTargetException e) {
						e.printStackTrace();
					} catch (InterruptedException e) {
						e.printStackTrace();
					}
					// Thread.yield();
				}
				progressBar.setStringPainted(false);
				progressBar.setString(null);
				progressBar.setIndeterminate(false);
				progressBar.setMinimum(0);
				progressBar.setMaximum(0);
				progressBar.setValue(0);
			}
		}).start();
	}

}
