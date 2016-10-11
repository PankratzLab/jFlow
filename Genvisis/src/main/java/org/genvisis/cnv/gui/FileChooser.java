package org.genvisis.cnv.gui;

import java.awt.Component;
import java.io.File;

import javax.swing.JFileChooser;

import org.genvisis.common.Logger;

public class FileChooser extends JFileChooser {

	private static final long serialVersionUID = 1L;
	private String[] files;
	private String navDir;
	boolean directoryChooser;
	boolean selected;

	/**
	 * @param parent parent component
	 * @param startDir starting directory to launch this chooser
	 * @param multipleSelect allow multiple file selections
	 * @param directoryChooser convert to directory chooser
	 * @param title title for the component
	 * @param log
	 */
	public FileChooser(	Component parent, String startDir, boolean multipleSelect,
											boolean directoryChooser, String title, Logger log) {
		super();
		selected = false;
		setMultiSelectionEnabled(multipleSelect);
		setCurrentDirectory(new File(startDir));
		setDialogTitle(title);
		files = null;
		navDir = "";
		if (directoryChooser) {
			setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
		}
		int result = showOpenDialog(parent);
		setVisible(true);

		if (result == JFileChooser.CANCEL_OPTION) {
			cancelSelection();
		} else if (directoryChooser) {
			navDir = getSelectedFile().getAbsolutePath() + "\\";
			selected = true;
		} else {
			File[] files = (multipleSelect ? getSelectedFiles() : new File[] {getSelectedFile()});
			this.files = getPaths(files);
			navDir = getCurrentDirectory().getAbsolutePath() + "\\";
			selected = true;
		}
	}

	public String[] getFiles() {
		return files;
	}

	public boolean isSelected() {
		return selected;
	}

	/**
	 * @return the current directory
	 */
	public String getNavDir() {
		return navDir;
	}

	/**
	 * @param files
	 * @return String[] of full paths
	 */
	public String[] getPaths(File[] files) {
		if (files == null || files.length == 0) {
			return null;
		}
		String[] paths = new String[files.length];
		for (int i = 0; i < files.length; i++) {
			paths[i] = files[i].getAbsolutePath();
		}
		return paths;
	}

	public static void main(String[] args) {
		String startDir = "D:/data/PoynterLinabery/";
		String logFile = "Chooseme.log";
		Logger log = new Logger(startDir + logFile);
		new FileChooser(null, startDir, true, true, "OPEN", log);
	}
}
