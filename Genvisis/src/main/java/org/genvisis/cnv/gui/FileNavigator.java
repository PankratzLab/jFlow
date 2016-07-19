package cnv.gui;

import java.awt.Color;
import java.io.File;
import java.util.ArrayList;

import javax.swing.JCheckBox;
import javax.swing.JPanel;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

public class FileNavigator extends JPanel implements ChangeListener {

	/**
     * 
     */
	private static final long serialVersionUID = 1L;
	private static ArrayList<String> fileList;

	/**
	 * Create the panel.
	 */
	public FileNavigator(String[] files, Color[] colors) {
		fileList = new ArrayList<String>();
		for (int i = 0; i < files.length; i++) {
			File file = new File(files[i]);
			String filename = file.getName();
			fileList.add(filename);
			JCheckBox fileBox = new JCheckBox(filename, true);
			fileBox.setToolTipText(file.getAbsolutePath());
			fileBox.setForeground(colors[i % colors.length]);
			fileBox.addChangeListener(this);
			add(fileBox);
		}
	}

	@Override
	public void stateChanged(ChangeEvent arg0) {
		ArrayList<String> oldFiles = new ArrayList<String>();
		oldFiles.addAll(fileList);
		JCheckBox box = (JCheckBox) arg0.getSource();
		String text = box.getText();
		if (box.isSelected()) {
			if (!fileList.contains(text)) {
				fileList.add(text);
			}
		} else {
			int index = fileList.indexOf(text);
			if (index >= 0) {
				fileList.remove(index);
			}
		}
		this.firePropertyChange("files", oldFiles, fileList);
	}
}