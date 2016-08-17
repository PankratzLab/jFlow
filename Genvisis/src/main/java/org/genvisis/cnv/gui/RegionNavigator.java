package org.genvisis.cnv.gui;

import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;

import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.plots.CompPlot;
import org.genvisis.cnv.var.Region;
import org.genvisis.common.Grafik;
import org.genvisis.common.Positions;
import org.genvisis.common.ext;

public class RegionNavigator extends JPanel implements ActionListener {
  /**
   * 
   */
  private static final long serialVersionUID = 1L;
  public static final String DEFAULT_LOCATION = "chr6:161,624,000-163,776,000"; // PARK2 region
  private JTextField textField;
  JButton firstButton, leftButton, rightButton, lastButton; // Navigation buttons
  JLabel location;
  String[] regionsList; // List of the region files
  Project proj;
  CompPlot plot;
  HashMap<String, ArrayList<Region>> regions = new HashMap<String, ArrayList<Region>>();
  String currentFile = "";
  int regionIndex = 0;

  int lastRegionIndex = 0;

  /**
   * Create the panel.
   */
  public RegionNavigator(CompPlot cp) {
    proj = cp.getProject();
    plot = cp;
    initButtons();

    // Parse the files and set up the regions Vector
    // This will be full of regions from the file, or the DEFAULT_LOCATION
    loadRegions();

    // for (int i = 0; i < regions.size(); i++) {
    // System.out.println("=(" + i + ")= " + regions.get(i).getRegion() + "\t" +
    // regions.get(i).getLabel());
    // }

    // Set region to the first region in the list
    setRegionFile(regionsList[0]);
    setRegion(regionIndex);
  }

  @Override
  /**
   * Set the region based on which buttons are pressed
   */
  public void actionPerformed(ActionEvent arg0) {
    JComponent source = (JComponent) arg0.getSource();

    lastRegionIndex = regionIndex;
    if (source.equals(firstButton)) {
      regionIndex = 0;
    } else if (source.equals(leftButton)) {
      if (regionIndex > 0) {
        regionIndex--;
      }
    } else if (source.equals(rightButton)) {
      if (regionIndex < (regions.get(currentFile).size() - 1)) {
        regionIndex++;
      }
    } else if (source.equals(lastButton)) {
      regionIndex = regions.get(currentFile).size() - 1;
    }
    if (lastRegionIndex != regionIndex) {
      setRegion(regionIndex);
    } else if (source.equals(getTextField())) {
      // TODO: Provide facility to export the list of regions, and expand the list of regions any
      // time someone enters one manually
      int[] newLocation = Positions.parseUCSClocation(getTextField().getText());
      // if ((newLocation[0] < 0) || (newLocation[1] < 0) || (newLocation[2] < 0)) {
      // JOptionPane.showMessageDialog(this.getParent(), "Invalid UCSC location - " +
      // textField.getText());
      // } else {
      setLocation(newLocation);

      // Pass along the property change
      firePropertyChange("location", regions.get(currentFile).get(lastRegionIndex),
          new Region(newLocation));
      // }
    }
  }

  /**
   * 
   * @return A Region object representing the current region
   */
  public Region getRegion() {
    return regions.get(currentFile).get(regionIndex);
  }

  public String getRegionFile() {
    return currentFile;
  }

  public JTextField getTextField() {
    return textField;
  }

  /**
   * Initialize the graphics and add listeners
   */
  public void initButtons() {
    firstButton = new JButton(Grafik.getImageIcon("images/firstLast/dFirst.gif"));
    firstButton.addActionListener(this);
    add(firstButton);

    leftButton = new JButton(Grafik.getImageIcon("images/firstLast/dLeft.gif"));
    leftButton.addActionListener(this);
    add(leftButton);

    setTextField(new JTextField());
    getTextField().setFont(new Font("Tahoma", Font.PLAIN, 14));
    getTextField().addActionListener(this);
    add(getTextField());
    getTextField().setColumns(20);

    rightButton = new JButton(Grafik.getImageIcon("images/firstLast/dRight.gif"));
    rightButton.addActionListener(this);
    add(rightButton);

    lastButton = new JButton(Grafik.getImageIcon("images/firstLast/dLast.gif"));
    lastButton.addActionListener(this);
    add(lastButton);

    location = new JLabel();
    add(location);
  }

  /**
   * Load the regions from the specified files in the project
   */
  public void loadRegions() {
    BufferedReader reader;
    // Get a list of the regions
    regionsList = proj.REGION_LIST_FILENAMES.getValue();
    regions = new HashMap<String, ArrayList<Region>>();

    try {
      if (regionsList.length == 0) {
        System.out.println("No regions file defined, using default of " + DEFAULT_LOCATION);
      } else {
        // Read each file line by line, format is:
        // chromosome:start-end \t label
        for (String element : regionsList) {
          // System.out.println("Parsing file " + regionsList[i]);
          // TODO Regex to ensure the line is formatted correctly
          reader = new BufferedReader(new FileReader(element));
          String line = null;
          while ((line = reader.readLine()) != null) {
            line = line.trim();
            if (!"".equals(line)) {
              Region myRegion = new Region(line);
              ArrayList<Region> list = regions.get(element);
              if (list == null) {
                list = new ArrayList<Region>();
                regions.put(element, list);
              }
              list.add(myRegion);
            }
          }
        }
      }
    } catch (Exception ex) {
      System.out.println(ex.getStackTrace());
    }

    ArrayList<Region> empty = new ArrayList<Region>();
    empty.add(new Region(DEFAULT_LOCATION));
    regions.put("", empty);
    if (regions.size() == 0) {
      // The file was invalid or didn't contain regions, create a default region
      // System.out.println("Setting default location");
      setRegionFile("");
    } else {
      setRegionFile(regionsList[0]);
    }
  }

  public void setLocation(int[] location) {
    getTextField().setText(Positions.getUCSCformat(location));
  }

  /**
   * Updates the region (text field, button tooltips, etc.) to reflect the new current position
   * 
   * @param index The index in the regions Vector to which this should be set
   */
  public void setRegion(int index) {
    regionIndex = index;

    // Update the text field
    getTextField().setText(regions.get(currentFile).get(regionIndex).getRegion());

    // Update the region indicator
    location.setText("Region " + (regionIndex + 1) + " of " + regions.get(currentFile).size());

    // Pass along the property change
    firePropertyChange("location", null/* regions.get(currentFile).get(lastRegionIndex) */,
        regions.get(currentFile).get(regionIndex));

    // Set the tooltip text on the buttons to match the region label if any
    firstButton.setToolTipText(regions.get(currentFile).get(0).getLabel());
    if (regionIndex >= 1) {
      leftButton.setToolTipText(regions.get(currentFile).get(regionIndex - 1).getLabel());
    } else {
      leftButton.setToolTipText(regions.get(currentFile).get(0).getLabel());
    }
    if (index < (regions.get(currentFile).size() - 1)) {
      rightButton.setToolTipText(regions.get(currentFile).get(regionIndex + 1).getLabel());
    } else {
      rightButton.setToolTipText(
          regions.get(currentFile).get(regions.get(currentFile).size() - 1).getLabel());
    }
    lastButton.setToolTipText(
        regions.get(currentFile).get(regions.get(currentFile).size() - 1).getLabel());
  }

  public void setRegionFile(String file) {
    if (file == null && currentFile.equals(file)) {
      return;// log?
    }
    if (!"".equals(file)) {
      String temp = ext.verifyDirFormat(file);
      temp = temp.substring(0, temp.length() - 1);
      if (!regions.containsKey(file) && !regions.containsKey(temp)) {
        proj.getLog().reportError("Error - file {" + file + "} is not a valid regions file");
        return;
      }
      currentFile = temp;
    } else {
      currentFile = file;
    }
  }

  public void setTextField(JTextField textField) {
    this.textField = textField;
  }
}
