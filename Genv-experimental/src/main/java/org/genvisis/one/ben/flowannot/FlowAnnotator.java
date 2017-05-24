package org.genvisis.one.ben.flowannot;

import java.awt.EventQueue;
import java.awt.Graphics;
import java.awt.Image;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.ActionMap;
import javax.swing.DefaultComboBoxModel;
import javax.swing.InputMap;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTextField;
import javax.swing.JTree;
import javax.swing.KeyStroke;
import javax.swing.border.BevelBorder;
import javax.swing.event.TreeSelectionEvent;
import javax.swing.event.TreeSelectionListener;
import javax.swing.plaf.basic.BasicTreeUI;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.DefaultTreeModel;

import net.miginfocom.swing.MigLayout;

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
		initialize();
	}

	/**
	 * Initialize the contents of the frame.
	 */
	private void initialize() {
		frmFlowannotator = new JFrame();
		frmFlowannotator.setTitle("FlowAnnotator");
		frmFlowannotator.setBounds(100, 100, 1200, 800);
		frmFlowannotator.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frmFlowannotator.getContentPane().setLayout(new MigLayout("", "[grow]", "[grow]"));

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
		controlPanel.add(fcsCombo, "cell 0 0,growx");

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

		JMenuBar menuBar = new JMenuBar();
		frmFlowannotator.setJMenuBar(menuBar);

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

		JMenuItem mnRecentPlacehldr = new JMenuItem("No Recent Files");
		mnRecentPlacehldr.setEnabled(false);
		mnLoadRecent.add(mnRecentPlacehldr);

		JMenuItem mntmSaveAnnot = new JMenuItem();
		mntmSaveAnnot.setAction(saveAction);
		mntmSaveAnnot.setText("Save Annotations to File");
		mntmSaveAnnot.setMnemonic('S');
		mnFile.add(mntmSaveAnnot);

		JMenuItem mntmExit = new JMenuItem("Exit");
		mnFile.add(mntmExit);
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
		JFileChooser jfc = new JFileChooser();
		jfc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
		jfc.setDialogTitle("Open Directory");
		int opt = jfc.showOpenDialog(frmFlowannotator);
		if (opt == JFileChooser.APPROVE_OPTION) {
			File f = jfc.getSelectedFile();
			annotator.loadImgDir(ext.verifyDirFormat(f.getAbsolutePath()));
			reloadControls();
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

	private void updateTreeKeys(JTree tree) {
		tree.setUI(new BasicTreeUI() {

			protected KeyListener createKeyListener() {
				return new KeyAdapter() {
					@Override
					public void keyPressed(KeyEvent e) {
						if (e.getKeyCode() == KeyEvent.VK_UP || e.getKeyCode() == KeyEvent.VK_DOWN) {
							super.keyPressed(e);
						} else if (e.getKeyCode() == KeyEvent.VK_LEFT || e.getKeyCode() == KeyEvent.VK_RIGHT) {
							if (e.getKeyCode() == KeyEvent.VK_LEFT) {
								int cur = fcsCombo.getSelectedIndex();
								if (cur > 0) {
									fcsCombo.setSelectedIndex(cur - 1);
								}
							} else if (e.getKeyCode() == KeyEvent.VK_RIGHT) {
								int cnt = fcsCombo.getItemCount();
								int cur = fcsCombo.getSelectedIndex();
								if (cur < cnt - 1) {
									fcsCombo.setSelectedIndex(cur + 1);
								}
							}
							e.consume();
						} else {
							if (!(e.isAltDown() && e.getKeyCode() == KeyEvent.VK_F4)) {
								e.consume();
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
		tree.setSelectionRow(0);
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
		JFileChooser jfc = new JFileChooser();
		jfc.setFileSelectionMode(JFileChooser.FILES_ONLY);
		jfc.setDialogTitle("Save Annotations to File");
		int opt = jfc.showSaveDialog(frmFlowannotator);
		if (opt == JFileChooser.APPROVE_OPTION) {
			String annFile = jfc.getSelectedFile().getAbsolutePath();
			annotator.saveAnnotations(annFile);
		}
	}

	private void loadAnnotations() {
		JFileChooser jfc = new JFileChooser();
		jfc.setFileSelectionMode(JFileChooser.FILES_ONLY);
		jfc.setDialogTitle("Load Annotations from File");
		int opt = jfc.showOpenDialog(frmFlowannotator);
		if (opt == JFileChooser.APPROVE_OPTION) {
			String annFile = jfc.getSelectedFile().getAbsolutePath();
			try {
				annotator.loadAnnotations(annFile);
				reloadControls();
				updateAvail();
			} catch (IOException e) {
				e.printStackTrace();
			}
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


}
