package cyto;

import java.awt.Dimension;
import java.awt.Toolkit;

import javax.swing.JFrame;
import javax.swing.JPanel;

import cnv.filesys.Project;
import common.Logger;

public class CytoGUI extends JFrame {

	public CytoPanel cytoPanel;
	public JPanel cytoView;
	private static final long serialVersionUID = 1L;
	private Project proj;
	private String startDir;
	private Logger log;

	public CytoGUI(Project proj, String startDir, Logger log) {
		this.proj=proj;
		this.startDir = startDir;
		this.log = log;
		createAndShowGUI();
	}

	private void createAndShowGUI() {
		setSize(400, 100);
		Dimension dim = Toolkit.getDefaultToolkit().getScreenSize();
		//launch to middle of screen
		this.setLocation(dim.width / 2 - this.getSize().width / 2, dim.height / 2 - this.getSize().height / 2);
		setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		cytoPanel = new CytoPanel(proj,this, (startDir == null ? "" : startDir), log);
		cytoPanel.setVisible(true);
		add(cytoPanel);
		setTitle("Genvisis");
		setVisible(true);
		
		
	}

	public static void main(String[] args) {
		String startDir = "D:/data/";
		String logFile = "Chooseme.log";
		String fileName = "C:/workspace/Genvisis/projects/HirchCyto.properties";
		Logger log = new Logger(startDir + logFile);
		new CytoGUI(new Project(fileName,false),startDir, log);
	}
}
