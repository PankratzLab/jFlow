package org.genvisis.one.ben.flowannot;

import java.awt.EventQueue;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.image.BufferedImage;

import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JSplitPane;
import javax.swing.JTabbedPane;
import javax.swing.JTree;
import javax.swing.border.BevelBorder;
import javax.swing.event.TreeSelectionEvent;
import javax.swing.event.TreeSelectionListener;
import javax.swing.tree.DefaultMutableTreeNode;

import net.miginfocom.swing.MigLayout;

import org.genvisis.one.ben.flowannot.IAnnotator.AnnotatedImage;

public class FlowAnnotator {

	private JFrame frmFlowannotator;
	private DefaultMutableTreeNode rootNode;

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
		frmFlowannotator.setBounds(100, 100, 800, 600);
		frmFlowannotator.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frmFlowannotator.getContentPane().setLayout(new MigLayout("", "[grow]", "[grow]"));

		JSplitPane splitPane = new JSplitPane();
		splitPane.setResizeWeight(0.9);
		frmFlowannotator.getContentPane().add(splitPane, "cell 0 0,grow");

		JPanel imagePanel = new JPanel();
		splitPane.setLeftComponent(imagePanel);
		imagePanel.setBorder(new BevelBorder(BevelBorder.RAISED, null, null, null, null));

		JTabbedPane tabbedPane = new JTabbedPane(JTabbedPane.TOP);
		splitPane.setRightComponent(tabbedPane);

		JPanel controlPanel = new JPanel();
		tabbedPane.addTab("Review", null, controlPanel, null);
		controlPanel.setLayout(new MigLayout("", "[110.00,grow]", "[][grow][][]"));

		fcsCombo = new JComboBox<String>();
		controlPanel.add(fcsCombo, "cell 0 0,growx");

		JButton btnUp = new JButton("Load Existing Annotations");
		controlPanel.add(btnUp, "cell 0 2,growx");

		JButton btnDown = new JButton("Annotate Files In Directory");
		btnDown.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				// TODO directory chooser
				String dir = "F:/Flow/Annotation/";
				annotator.loadImgDir(dir);
				reloadControls();
			}
		});
		controlPanel.add(btnDown, "cell 0 3,growx");

		JPanel annotPanel = new JPanel();
		tabbedPane.addTab("Annotations", null, annotPanel, null);

		JMenuBar menuBar = new JMenuBar();
		frmFlowannotator.setJMenuBar(menuBar);

		JMenu mnFile = new JMenu("File");
		menuBar.add(mnFile);

		JMenuItem mntmLoadImgDir = new JMenuItem("Load Image Directory");
		mnFile.add(mntmLoadImgDir);

		JMenuItem mntmLoadAnnotFile = new JMenuItem("Load Annotation File");
		mnFile.add(mntmLoadAnnotFile);

		JMenu mnLoadRecent = new JMenu("Load Recent Annotation File");
		mnFile.add(mnLoadRecent);

		JMenuItem mnRecentPlacehldr = new JMenuItem("No Recent Files");
		mnRecentPlacehldr.setEnabled(false);
		mnLoadRecent.add(mnRecentPlacehldr);

		JMenuItem mntmSaveAnnot = new JMenuItem("Save Annotations to File");
		mnFile.add(mntmSaveAnnot);

		JMenuItem mntmExit = new JMenuItem("Exit");
		mnFile.add(mntmExit);
	}

	private void reloadControls() {

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

	IAnnotator annotator = new Annotator();

	private void setSelectedNode(AnnotatedImage node) {
		BufferedImage img = node.getImage();
	}

	TreeSelectionListener treeListener = new TreeSelectionListener() {
		@Override
		public void valueChanged(TreeSelectionEvent e) {
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
	private JComboBox fcsCombo;



}
