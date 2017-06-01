package org.genvisis.one.ben.flowannot;

import java.awt.EventQueue;
import java.awt.Graphics;
import java.awt.Image;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Properties;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.ActionMap;
import javax.swing.ButtonGroup;
import javax.swing.DefaultComboBoxModel;
import javax.swing.InputMap;
import javax.swing.JCheckBox;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JRadioButtonMenuItem;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JSplitPane;
import javax.swing.JTextField;
import javax.swing.JTree;
import javax.swing.KeyStroke;
import javax.swing.SwingConstants;
import javax.swing.border.BevelBorder;
import javax.swing.event.TreeSelectionEvent;
import javax.swing.event.TreeSelectionListener;
import javax.swing.plaf.basic.BasicTreeUI;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.DefaultTreeModel;
import javax.swing.tree.TreePath;

import net.miginfocom.swing.MigLayout;

import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.ext;

public class FlowAnnotator {

	private JFrame frmFlowannotator;

	private JComboBox<String> fcsCombo;
	private JTree tree;
	IAnnotator annotator = new Annotator();

	BufferedImage selectedImage = AnnotatedImage.createReadyImage();

	private JMenu mnLoadRecent;

	private JTextField fldNewAnnot;

	private JPanel annotPanel;

	private JPanel imagePanel;

	private static final String PROP_FILE = ".flowannotator.properties";

	/**
	 * Launch the application.
	 */
	public static void main(String[] args) {
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					FlowAnnotator window = new FlowAnnotator();
					window.frmFlowannotator.setVisible(true);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		});
	}

	/**
	 * Create the application.
	 */
	public FlowAnnotator() {
		loadProperties();
		initialize();
	}

	private void loadProperties() {
		Properties p = new Properties();
		if (!Files.exists(PROP_FILE))
			return;
		try {
			p.load(Files.getAppropriateInputStreamReader(PROP_FILE));
			String lst = p.getProperty(KEY_LAST_OPENED, "");
			String rec = p.getProperty(KEY_RECENT, "");
			lastOpenedJFC = lst;
			if (!"".equals(rec)) {
				String[] ps = rec.split(";");
				for (String s : ps) {
					if (Files.exists(s)) {
						recentAnnotFiles.add(s);
					}
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private void saveProperties() {
		Properties p = new Properties();
		if (lastOpenedJFC != null) {
			p.setProperty(KEY_RECENT, lastOpenedJFC);
		}
		if (!recentAnnotFiles.isEmpty()) {
			p.setProperty(KEY_LAST_OPENED, ArrayUtils.toStr(recentAnnotFiles, ";"));
		}
		try {
			p.store(Files.getAppropriateWriter(PROP_FILE), "");
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private static final String KEY_LAST_OPENED = "LAST_OPEN";
	private String lastOpenedJFC = null;
	private static final String KEY_RECENT = "RECENT";
	private HashSet<String> recentAnnotFiles = new HashSet<>();

	/**
	 * Initialize the contents of the frame.
	 */
	private void initialize() {
		frmFlowannotator = new JFrame();
		frmFlowannotator.setTitle("FlowAnnotator");
		frmFlowannotator.setBounds(100, 100, 1200, 800);
		frmFlowannotator.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frmFlowannotator.getContentPane().setLayout(new MigLayout("", "[grow]", "[grow]"));
		frmFlowannotator.addWindowListener(new WindowAdapter() {
			@Override
			public void windowClosing(WindowEvent e) {
				close();
				super.windowClosing(e);
			}
		});

		splitPane = new JSplitPane();
		splitPane.setResizeWeight(0.78);
		frmFlowannotator.getContentPane().add(splitPane, "cell 0 0,grow");

		imagePanel = new JPanel() {
			@Override
			protected void paintComponent(Graphics g) {
				super.paintComponent(g);
				Image scaledImage = selectedImage.getScaledInstance(getWidth(), getHeight(),
																														Image.SCALE_SMOOTH);
				g.drawImage(scaledImage, 0, 0, getWidth(), getHeight(), null);
			}
		};
		splitPane.setLeftComponent(imagePanel);
		imagePanel.setBorder(new BevelBorder(BevelBorder.RAISED, null, null, null, null));

		JSplitPane splitPane_1 = new JSplitPane();
		splitPane_1.setResizeWeight(0);
		splitPane_1.setOrientation(JSplitPane.VERTICAL_SPLIT);
		splitPane.setRightComponent(splitPane_1);

		annotPanel = new JPanel();
		splitPane_1.setRightComponent(annotPanel);
		annotPanel.setLayout(new MigLayout("", "[grow]", "[][][][grow]"));

		JLabel lblAddAnnotation = new JLabel("Add Annotation:");
		annotPanel.add(lblAddAnnotation, "cell 0 0");

		fldNewAnnot = new JTextField();
		ActionMap am = fldNewAnnot.getActionMap();
		InputMap im = fldNewAnnot.getInputMap(JComponent.WHEN_FOCUSED);
		annotPanel.add(fldNewAnnot, "cell 0 1,growx");

		JLabel lblExistingAnnotations = new JLabel("Existing Annotations:");
		annotPanel.add(lblExistingAnnotations, "cell 0 2");

		JScrollPane scrollPane = new JScrollPane();
		annotPanel.add(scrollPane, "cell 0 3,grow");

		annotPanel = new JPanel();
		scrollPane.setViewportView(annotPanel);
		annotPanel.setLayout(new MigLayout("ins 0", "[]", "[]"));

		controlPanel = new JPanel();
		splitPane_1.setLeftComponent(controlPanel);
		controlPanel.setLayout(new MigLayout("", "[110.00,grow]", "[][grow][][]"));

		fcsCombo = new JComboBox<String>();
		fcsCombo.addActionListener(comboListener);
		fcsCombo.setFocusable(false);
		controlPanel.add(fcsCombo, "flowx,cell 0 0,growx");

		JScrollPane scrollPane_1 = new JScrollPane();
		controlPanel.add(scrollPane_1, "cell 0 1,grow");

		tree = new JTree();
		scrollPane_1.setViewportView(tree);
		DefaultMutableTreeNode rootNode = new DefaultMutableTreeNode();
		rootNode.add(GateTree.constructTree());
		DefaultTreeModel dtm = new DefaultTreeModel(rootNode);
		tree.setModel(dtm);
		tree.setShowsRootHandles(true);
		tree.setRootVisible(false);
		tree.addTreeSelectionListener(treeListener);
		tree.setEnabled(false);
		updateTreeKeys(tree);
		expandAllNodes(tree);
		im.put(KeyStroke.getKeyStroke(KeyEvent.VK_ENTER, 0), "enter");
		am.put("enter", new AbstractAction() {
			@Override
			public void actionPerformed(ActionEvent arg0) {
				createNewAnnotation();
			}
		});

		frmFlowannotator.setJMenuBar(createMenuBar());
	}

	private JMenuBar createMenuBar() {
		JMenuBar menuBar = new JMenuBar();

		JMenu mnFile = new JMenu("File");
		mnFile.setMnemonic('F');
		menuBar.add(mnFile);

		JMenuItem mntmLoadImgDir = new JMenuItem();
		mntmLoadImgDir.setAction(loadAction);
		mntmLoadImgDir.setText("Load Image Directory");
		mntmLoadImgDir.setMnemonic('L');
		mnFile.add(mntmLoadImgDir);

		JMenuItem mntmLoadAnnotFile = new JMenuItem();
		mntmLoadAnnotFile.setAction(existAction);
		mntmLoadAnnotFile.setText("Load Annotation File");
		mntmLoadAnnotFile.setMnemonic('A');
		mnFile.add(mntmLoadAnnotFile);

		mnLoadRecent = new JMenu("Load Recent Annotation File");
		mnLoadRecent.setMnemonic('R');
		mnFile.add(mnLoadRecent);

		updateRecentFiles();

		JMenuItem mntmSaveAnnot = new JMenuItem();
		mntmSaveAnnot.setAction(saveAction);
		mntmSaveAnnot.setText("Save Annotations to File");
		mntmSaveAnnot.setMnemonic('S');
		mnFile.add(mntmSaveAnnot);

		JMenuItem mntmExit = new JMenuItem();
		mntmExit.setAction(new AbstractAction() {
			@Override
			public void actionPerformed(ActionEvent e) {
				close();
			}
		});
		mntmExit.setText("Exit");
		mnFile.add(mntmExit);

		JMenu mnNav = new JMenu("Navigation");
		mnNav.setMnemonic('N');
		menuBar.add(mnNav);

		JMenuItem mntmNavTravLbl = new JMenuItem();
		mntmNavTravLbl.setEnabled(false);
		mntmNavTravLbl.setText("Traversal:");
		mnNav.add(mntmNavTravLbl);

		AbstractAction travSelAction = new AbstractAction() {
			@Override
			public void actionPerformed(ActionEvent e) {
				// System.out.println(e);
			}
		};
		ButtonGroup travGrp = new ButtonGroup();

		jrbTravAll = new JRadioButtonMenuItem();
		jrbTravAll.setAction(travSelAction);
		jrbTravAll.setText("All");
		jrbTravAll.setSelected(true);
		travGrp.add(jrbTravAll);
		mnNav.add(jrbTravAll);
		jrbTravAnn = new JRadioButtonMenuItem();
		jrbTravAnn.setAction(travSelAction);
		jrbTravAnn.setText("Annotated");
		travGrp.add(jrbTravAnn);
		mnNav.add(jrbTravAnn);
		jrbTravNon = new JRadioButtonMenuItem();
		jrbTravNon.setAction(travSelAction);
		jrbTravNon.setText("Non-Annotated");
		travGrp.add(jrbTravNon);
		mnNav.add(jrbTravNon);
		mnTravAnn = new JMenu("Annotation:");
		mnNav.add(mnTravAnn);

		mnNav.add(new JSeparator(SwingConstants.HORIZONTAL));

		JMenuItem mntmNavOptLbl = new JMenuItem();
		mntmNavOptLbl.setEnabled(false);
		mntmNavOptLbl.setText("Other Options:");
		mnNav.add(mntmNavOptLbl);

		JCheckBoxMenuItem mntmNavKeepGate = new JCheckBoxMenuItem();
		mntmNavKeepGate.setAction(new AbstractAction() {
			@Override
			public void actionPerformed(ActionEvent e) {
				keepGateWhenFileChange = mntmNavKeepGate.isSelected();
			}
		});
		mntmNavKeepGate.setText("Keep Gate When File Changes");
		mntmNavKeepGate.setSelected(true);
		mnNav.add(mntmNavKeepGate);

		return menuBar;
	}

	private void close() {
		saveProperties();
		frmFlowannotator.setVisible(false);
		frmFlowannotator.dispose();
	}

	private void updateRecentFiles() {
		mnLoadRecent.removeAll();
		if (recentAnnotFiles.isEmpty()) {
			JMenuItem mnRecentPlacehldr = new JMenuItem("No Recent Files");
			mnRecentPlacehldr.setEnabled(false);
			mnLoadRecent.add(mnRecentPlacehldr);
		} else {
			for (String s : recentAnnotFiles) {
				JMenuItem mnRecent = new JMenuItem();
				mnRecent.setAction(new AbstractAction() {
					@Override
					public void actionPerformed(ActionEvent e) {
						loadAnnotationFile(s);
					}
				});
				mnRecent.setText(ext.removeDirectoryInfo(s));
				mnLoadRecent.add(mnRecent);
			}
		}
	}

	protected void createNewAnnotation() {
		String newAnnotation = fldNewAnnot.getText().trim();
		fldNewAnnot.setText("");
		if ("".equals(newAnnotation)) {
			return;
		}
		AnnotatedImage.Annotation newAnn = new AnnotatedImage.Annotation(newAnnotation);
		if (!annotator.getAnnotations().contains(newAnn)) {
			annotator.addNewAnnotation(newAnn);
			refreshAnnotations();
		}
	}

	private void loadFiles() {
		JFileChooser jfc = new JFileChooser(lastOpenedJFC == null ? "" : lastOpenedJFC);
		jfc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
		jfc.setDialogTitle("Open Directory");
		int opt = jfc.showOpenDialog(frmFlowannotator);
		if (opt == JFileChooser.APPROVE_OPTION) {
			File f = jfc.getSelectedFile();
			String fS = ext.verifyDirFormat(f.getAbsolutePath());
			lastOpenedJFC = fS;
			annotator.loadImgDir(fS);
			reloadControls();
			saveProperties();
		}
	}

	private void reloadControls() {
		// TODO track which is selected and re-select after reloading controls?
		DefaultComboBoxModel<String> dcbm = new DefaultComboBoxModel<>(annotator.getFCSKeys()
																																						.toArray(new String[0]));
		fcsCombo.setModel(dcbm);

		fcsCombo.setSelectedIndex(0);
		tree.requestFocusInWindow();
	}

	private void refreshAnnotations() {
		ArrayList<AnnotatedImage.Annotation> allAnnots = annotator.getAnnotations();
		annotPanel.removeAll();
		AnnotatedImage ai = null;
		if (tree.getSelectionModel().getSelectionPath() != null) {
			DefaultMutableTreeNode dmtn = (DefaultMutableTreeNode) tree.getSelectionModel()
																																 .getSelectionPath()
																																 .getLastPathComponent();
			ai = (AnnotatedImage) dmtn.getUserObject();
		}
		for (int i = 0; i < allAnnots.size(); i++) {
			JCheckBox annBox = new JCheckBox();
			if (ai != null && ai.getAnnotations().contains(allAnnots.get(i))) {
				annBox.setSelected(true);
			} else {
				annBox.setSelected(false);
			}
			final AnnotatedImage ai1 = ai;
			final int ind = i;
			AbstractAction mnemAct = new AbstractAction() {
				@Override
				public void actionPerformed(ActionEvent e) {
					if (ai1 != null) {
						AnnotatedImage.Annotation a = allAnnots.get(ind);
						ArrayList<AnnotatedImage.Annotation> myAnn = ai1.getAnnotations();
						if (myAnn.contains(a)) {
							myAnn.remove(a);
							annBox.setSelected(false);
						} else {
							myAnn.add(a);
							annBox.setSelected(true);
						}
						annotPanel.repaint();
					}
				}
			};
			annBox.setAction(mnemAct);
			annBox.setText(allAnnots.get(i).annotation);
			annBox.setMnemonic(allAnnots.get(i).mnemonic);
			annBox.setFocusable(false);
			mnemonicActions.put(allAnnots.get(i).mnemonic, mnemAct);
			annotPanel.add(annBox, "cell 0 " + i);
		}
		annotPanel.revalidate();
	}

	HashMap<Character, Action> mnemonicActions = new HashMap<>();

	private void fireMnem(char code) {
		if (mnemonicActions.containsKey(code)) {
			mnemonicActions.get(code).actionPerformed(null);
		}
	}

	boolean keepGateWhenFileChange = true;

	private static final int ALL = 0;
	private static final int ANN = 1;
	private static final int NON = 2;

	private int getTraversal() {
		if (jrbTravAll.isSelected()) {
			return ALL;
		} else if (jrbTravAnn.isSelected()) {
			return ANN;
		} else if (jrbTravNon.isSelected()) {
			return NON;
		}
		return 4;
	}

	private int getPrevNodeRow() {
		int trav = getTraversal();
		if (tree.isSelectionEmpty()) {
			return -1;
		}
		int[] rows = tree.getSelectionRows();
		if (rows == null || rows.length == 0) {
			return -1;
		}
		int row = rows[0];
		row--;
		while (row >= 0) {
			if (trav == ALL) {
				break;
			} else if (trav == ANN) {
				TreePath tp = tree.getPathForRow(row);
				DefaultMutableTreeNode dmtn = (DefaultMutableTreeNode) tp.getLastPathComponent();
				AnnotatedImage ai = (AnnotatedImage) dmtn.getUserObject();
				if (!ai.getAnnotations().isEmpty()) {
					break;
				}
			} else if (trav == NON) {
				TreePath tp = tree.getPathForRow(row);
				DefaultMutableTreeNode dmtn = (DefaultMutableTreeNode) tp.getLastPathComponent();
				AnnotatedImage ai = (AnnotatedImage) dmtn.getUserObject();
				if (ai.getAnnotations().isEmpty()) {
					break;
				}
			} else {
				// TODO search for node with annotation
			}
			row--;
		}
		return row;
	}

	private int getNextNodeRow() {
		int trav = getTraversal();
		if (tree.isSelectionEmpty()) {
			return -1;
		}
		int cnt = getTreeRowCount();
		int[] rows = tree.getSelectionRows();
		if (rows == null || rows.length == 0) {
			return -1;
		}
		int row = rows[0];
		row++;
		while (row < cnt) {
			if (trav == ALL) {
				break;
			} else if (trav == ANN) {
				TreePath tp = tree.getPathForRow(row);
				DefaultMutableTreeNode dmtn = (DefaultMutableTreeNode) tp.getLastPathComponent();
				AnnotatedImage ai = (AnnotatedImage) dmtn.getUserObject();
				if (!ai.getAnnotations().isEmpty()) {
					break;
				}
			} else if (trav == NON) {
				TreePath tp = tree.getPathForRow(row);
				DefaultMutableTreeNode dmtn = (DefaultMutableTreeNode) tp.getLastPathComponent();
				AnnotatedImage ai = (AnnotatedImage) dmtn.getUserObject();
				if (ai.getAnnotations().isEmpty()) {
					break;
				}
			} else {
				// TODO search for node with annotation
			}
			row++;
		}
		return row;
	}

	private int getPrevFile() {
		int trav = getTraversal();
		int sel = fcsCombo.getSelectedIndex();
		String gate = ((AnnotatedImage) ((DefaultMutableTreeNode) tree.getPathForRow(0)
																																	.getLastPathComponent()).getUserObject()).getGateName();
		// default to root, or last leaf?, or any ann/non?
		if (keepGateWhenFileChange && !tree.isSelectionEmpty()) {
			gate = ((AnnotatedImage) ((DefaultMutableTreeNode) tree.getSelectionPath()
																														 .getLastPathComponent()).getUserObject()).getGateName();
		}
		sel--;
		while (sel >= 0) {
			if (trav == ALL) {
				break;
			} else if (trav == ANN) {
				String fcs = fcsCombo.getItemAt(sel);
				AnnotatedImage ai = annotator.getAnnotationMap().get(fcs).get(gate);
				if (ai != null && ai.getAnnotations().size() > 0) {
					break;
				}
			} else if (trav == NON) {
				String fcs = fcsCombo.getItemAt(sel);
				AnnotatedImage ai = annotator.getAnnotationMap().get(fcs).get(gate);
				if (ai == null || ai.getAnnotations().size() == 0) {
					break;
				}
			}
			sel--;
		}
		return sel;
	}

	private int getNextFile() {
		int trav = getTraversal();
		int sel = fcsCombo.getSelectedIndex();
		String gate = ((AnnotatedImage) ((DefaultMutableTreeNode) tree.getPathForRow(0)
																																	.getLastPathComponent()).getUserObject()).getGateName();
		// default to root, or last leaf?, or any ann/non?
		if (keepGateWhenFileChange && !tree.isSelectionEmpty()) {
			gate = ((AnnotatedImage) ((DefaultMutableTreeNode) tree.getSelectionPath()
																														 .getLastPathComponent()).getUserObject()).getGateName();
		}
		sel++;
		while (sel < fcsCombo.getItemCount()) {
			if (trav == ALL) {
				break;
			} else if (trav == ANN) {
				String fcs = fcsCombo.getItemAt(sel);
				if (annotator.getAnnotationMap().get(fcs).get(gate).getAnnotations().size() > 0) {
					break;
				}
			} else if (trav == NON) {
				String fcs = fcsCombo.getItemAt(sel);
				if (annotator.getAnnotationMap().get(fcs).get(gate).getAnnotations().size() == 0) {
					break;
				}
			}
			sel++;
		}
		return sel;
	}

	private void keyUp() { // prev node in tree
		int newRow = getPrevNodeRow();
		if (newRow < 0)
			return;
		tree.setSelectionRow(newRow);
	}

	private int getTreeRowCount() {
		int cnt = getChildCount((DefaultMutableTreeNode) ((DefaultTreeModel) tree.getModel()).getRoot());
		return cnt;
	}

	private int getChildCount(DefaultMutableTreeNode parent) {
		int cnt = ((DefaultTreeModel) tree.getModel()).getChildCount(parent);
		Enumeration e = parent.children();
		while (e.hasMoreElements()) {
			Object o = e.nextElement();
			cnt += getChildCount((DefaultMutableTreeNode) o);
		}
		return cnt;
	}

	private void keyDown() { // next node in tree
		int newRow = getNextNodeRow();
		if (newRow >= getTreeRowCount())
			return;
		tree.setSelectionRow(newRow);
	}

	private void keyLeft() { // prev file
		int prev = getPrevFile();
		int[] rows = tree.getSelectionRows();
		int row = -1;
		if (rows != null && rows.length > 0) {
			row = rows[0];
		}
		if (prev >= 0) {
			fcsCombo.setSelectedIndex(prev);
		}
		tree.setSelectionRow(keepGateWhenFileChange && row != -1 ? row : 0);
	}

	private void keyRight() { // next file
		int next = getNextFile();
		int[] rows = tree.getSelectionRows();
		int row = -1;
		if (rows != null && rows.length > 0) {
			row = rows[0];
		}
		if (next < fcsCombo.getItemCount()) {
			fcsCombo.setSelectedIndex(next);
			tree.setSelectionRow(keepGateWhenFileChange && row != -1 ? row : 0);
		}
	}

	private void updateTreeKeys(JTree tree) {
		tree.setUI(new BasicTreeUI() {

			protected KeyListener createKeyListener() {
				return new KeyAdapter() {
					@Override
					public void keyPressed(KeyEvent e) {
						e.consume();
						if (e.getKeyCode() == KeyEvent.VK_UP) {
							keyUp();
						} else if (e.getKeyCode() == KeyEvent.VK_DOWN) {
							keyDown();
						} else if (e.getKeyCode() == KeyEvent.VK_LEFT) {
							keyLeft();
						} else if (e.getKeyCode() == KeyEvent.VK_RIGHT) {
							keyRight();
						} else {
							if (!(e.isAltDown() && e.getKeyCode() == KeyEvent.VK_F4)) {
								fireMnem((e.getKeyChar() + "").toUpperCase().charAt(0));
							}
						}
					}
				};
			}
		});
	}

	private void expandAllNodes(JTree tree) {
		int j = tree.getRowCount();
		int i = 0;
		while (i < j) {
			tree.expandRow(i);
			i += 1;
			j = tree.getRowCount();
		}
	}

	private void setSelectedNode(AnnotatedImage node) {
		BufferedImage img = node.getImage();
		selectedImage = img;
		imagePanel.repaint();
		refreshAnnotations();
	}

	private void updateAvail() {
		HashMap<String, HashMap<String, AnnotatedImage>> map = annotator.getAnnotationMap();
		String fcsFile = (String) fcsCombo.getSelectedItem();
		HashMap<String, AnnotatedImage> ann = map.get(fcsFile);
		DefaultMutableTreeNode root = (DefaultMutableTreeNode) ((DefaultMutableTreeNode) tree.getModel()
																																												 .getRoot()).getChildAt(0);
		updateNode(root, ann);
		((DefaultTreeModel) tree.getModel()).reload();
		expandAllNodes(tree);
		updateTreeKeys(tree);
		tree.repaint();
	}

	private void updateNode(DefaultMutableTreeNode node, HashMap<String, AnnotatedImage> annMap) {
		AnnotatedImage ai = (AnnotatedImage) node.getUserObject();
		AnnotatedImage newAi = annMap.get(ai.getGateName());
		if (newAi == null) {
			newAi = new AnnotatedImage(ai.getGateName(), ai.isRoot());
			newAi.setMissing(true);
		}
		node.setUserObject(newAi);
		for (int i = 0; i < node.getChildCount(); i++) {
			updateNode((DefaultMutableTreeNode) node.getChildAt(i), annMap);
		}
	}

	private void saveAnnotations() {
		JFileChooser jfc = new JFileChooser(lastOpenedJFC == null ? "" : lastOpenedJFC);
		jfc.setFileSelectionMode(JFileChooser.FILES_ONLY);
		jfc.setDialogTitle("Save Annotations to File");
		int opt = jfc.showSaveDialog(frmFlowannotator);
		if (opt == JFileChooser.APPROVE_OPTION) {
			String annFile = jfc.getSelectedFile().getAbsolutePath();
			lastOpenedJFC = ext.verifyDirFormat(ext.parseDirectoryOfFile(annFile));
			annotator.saveAnnotations(annFile);
			recentAnnotFiles.add(annFile);
			saveProperties();
		}
	}

	private void loadAnnotations() {
		JFileChooser jfc = new JFileChooser(lastOpenedJFC == null ? "" : lastOpenedJFC);
		jfc.setFileSelectionMode(JFileChooser.FILES_ONLY);
		jfc.setDialogTitle("Load Annotations from File");
		int opt = jfc.showOpenDialog(frmFlowannotator);
		if (opt == JFileChooser.APPROVE_OPTION) {
			String annFile = jfc.getSelectedFile().getAbsolutePath();
			loadAnnotationFile(annFile);
		}
	}

	private void loadAnnotationFile(String annFile) {
		try {
			annotator.loadAnnotations(annFile);
			lastOpenedJFC = ext.verifyDirFormat(ext.parseDirectoryOfFile(annFile));
			recentAnnotFiles.add(annFile);
			reloadControls();
			updateAvail();
			tree.setSelectionRow(0);
			saveProperties();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	ActionListener comboListener = new ActionListener() {
		@Override
		public void actionPerformed(ActionEvent e) {
			tree.setEnabled(true);
			updateAvail();
		}
	};

	TreeSelectionListener treeListener = new TreeSelectionListener() {
		@Override
		public void valueChanged(TreeSelectionEvent e) {
			if (e.getNewLeadSelectionPath() == null)
				return;
			Object selected = e.getNewLeadSelectionPath().getLastPathComponent();
			if (!(selected instanceof DefaultMutableTreeNode)) {
				System.err.println("Error - selected an invalid item somehow.  Likely programmer error.");
				return;
			}
			DefaultMutableTreeNode selectedNode = (DefaultMutableTreeNode) selected;
			AnnotatedImage ann = ((AnnotatedImage) selectedNode.getUserObject());
			setSelectedNode(ann);
		}

	};

	private final Action loadAction = new AbstractAction() {
		@Override
		public void actionPerformed(ActionEvent e) {
			loadFiles();
		}
	};

	private final Action existAction = new AbstractAction() {
		@Override
		public void actionPerformed(ActionEvent e) {
			loadAnnotations();
		}
	};

	private final Action saveAction = new AbstractAction() {
		@Override
		public void actionPerformed(ActionEvent e) {
			saveAnnotations();
		}
	};

	private JPanel controlPanel;

	private JSplitPane splitPane;

	private JMenu mnTravAnn;

	private JRadioButtonMenuItem jrbTravAll;

	private JRadioButtonMenuItem jrbTravAnn;

	private JRadioButtonMenuItem jrbTravNon;


}
