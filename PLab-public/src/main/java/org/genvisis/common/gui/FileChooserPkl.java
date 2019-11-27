package org.genvisis.common.gui;

import java.awt.Component;
import java.io.File;

import javax.swing.JFileChooser;

import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;

public class FileChooserPkl extends JFileChooser {

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
  public FileChooserPkl(Component parent, String startDir, boolean multipleSelect,
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
      navDir = ext.verifyDirFormat(getSelectedFile().getAbsolutePath());
      selected = true;
    } else {
      File[] files = (multipleSelect ? getSelectedFiles() : new File[] {getSelectedFile()});
      this.files = getPaths(files);
      navDir = ext.verifyDirFormat(getCurrentDirectory().getAbsolutePath());
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
}
