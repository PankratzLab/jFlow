package org.genvisis.cnv.gui;

import java.awt.Component;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSlider;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.genvisis.cnv.plots.CompPlot;
import org.genvisis.filesys.CNVariant;

public class CompConfig extends JPanel implements ChangeListener, ActionListener {

  /**
   *
   */
  private static final long serialVersionUID = 1L;

  private String displayMode = "Full";
  private int probes = 0; // Default to 0 probes
  private int minSize = 0; // Default to no minimum size
  private int qualityScore = 0; // Default to a score of 0
  private int rectangleHeight = 10; // Default rectangle height
  CompPlot compPlot;

  JSlider probesSlider;
  JLabel lblProbes;

  JSlider minSizeSlider;
  JLabel lblMinimumSizekb;

  JSlider qualityScoreSlider;
  JLabel lblQualityScore;

  JSlider rectangleHeightSlider;
  JLabel lblRectangleHeight;

  CNVPanel cnvPanel;

  //

  /**
   * Create the panel.
   */
  public CompConfig(CompPlot cp) {
    setLayout(new BoxLayout(this, BoxLayout.PAGE_AXIS));
    compPlot = cp;

    cnvPanel = new CNVPanel(cp);
    cnvPanel.setDisplayMode(displayMode);
    add(cnvPanel);

    add(Box.createGlue());

    JPanel configPanel = new JPanel();
    configPanel.setLayout(new BoxLayout(configPanel, BoxLayout.PAGE_AXIS));

    // Configure whether the display mode is Full (variants displayed are separated and color-coded
    // by file)
    // or abbreviated (CNVs with multiple copies appear as only one block)
    JPanel dmPanel = new JPanel(new FlowLayout());
    JLabel lblDisplayMode = new JLabel("Display Mode");
    dmPanel.add(lblDisplayMode);

    JComboBox comboBox = new JComboBox();
    comboBox.setModel(new DefaultComboBoxModel(new String[] {"Full", "Pack", "Collapsed"}));
    comboBox.addActionListener(this);
    dmPanel.add(comboBox);
    configPanel.add(dmPanel);

    // Filter by the number of probes
    lblProbes = new JLabel("Probes (" + probes + ")");
    lblProbes.setAlignmentX(Component.CENTER_ALIGNMENT);
    configPanel.add(lblProbes);

    add(Box.createGlue());

    probesSlider = new JSlider(JSlider.HORIZONTAL, 0, 50, probes);
    probesSlider.setMajorTickSpacing(10);
    probesSlider.setMinorTickSpacing(5);
    probesSlider.setPaintTicks(true);
    probesSlider.setPaintLabels(true);
    probesSlider.addChangeListener(this);
    configPanel.add(probesSlider);

    add(Box.createGlue());

    // Filter by the minimum size in kilobases
    lblMinimumSizekb = new JLabel("Minimum Size in kb (" + minSize + ")");
    lblMinimumSizekb.setAlignmentX(Component.CENTER_ALIGNMENT);
    configPanel.add(lblMinimumSizekb);

    add(Box.createGlue());

    minSizeSlider = new JSlider(JSlider.HORIZONTAL, 0, 1000, minSize);
    minSizeSlider.setMajorTickSpacing(200);
    minSizeSlider.setMinorTickSpacing(100);
    minSizeSlider.setPaintTicks(true);
    minSizeSlider.setPaintLabels(true);
    minSizeSlider.addChangeListener(this);
    configPanel.add(minSizeSlider);

    add(Box.createGlue());

    // Filter by the minimum quality score
    lblQualityScore = new JLabel("Quality Score (" + qualityScore + ")");
    lblQualityScore.setAlignmentX(Component.CENTER_ALIGNMENT);
    configPanel.add(lblQualityScore);

    add(Box.createGlue());

    qualityScoreSlider = new JSlider(0, 100, qualityScore);
    qualityScoreSlider.setMajorTickSpacing(20);
    qualityScoreSlider.setMinorTickSpacing(10);
    qualityScoreSlider.setPaintTicks(true);
    qualityScoreSlider.setPaintLabels(true);
    qualityScoreSlider.addChangeListener(this);
    configPanel.add(qualityScoreSlider);

    add(Box.createGlue());

    // Scale the height of the rectangles
    lblRectangleHeight = new JLabel("Rectangle Height (" + rectangleHeight + ")");
    lblRectangleHeight.setAlignmentX(Component.CENTER_ALIGNMENT);
    configPanel.add(lblRectangleHeight);

    add(Box.createGlue());

    rectangleHeightSlider = new JSlider(1, 25, rectangleHeight);
    rectangleHeightSlider.setMajorTickSpacing(24);
    rectangleHeightSlider.setMinorTickSpacing(1);
    rectangleHeightSlider.setPaintTicks(true);
    rectangleHeightSlider.setPaintLabels(true);
    rectangleHeightSlider.addChangeListener(this);
    configPanel.add(rectangleHeightSlider);

    add(Box.createGlue());

    add(configPanel);
  }

  // Monitor the sliders for changes
  @Override
  public void stateChanged(ChangeEvent arg0) {
    JSlider source = (JSlider) arg0.getSource();
    if (source.equals(probesSlider)) {
      int p = source.getValue();
      this.firePropertyChange("probes", probes, p);
      probes = p;
      lblProbes.setText("Probes (" + probes + ")");
    } else if (source.equals(minSizeSlider)) {
      int ms = source.getValue();
      this.firePropertyChange("minSize", minSize, ms);
      minSize = ms;
      lblMinimumSizekb.setText("Minimum Size in kb (" + minSize + ")");
    } else if (source.equals(qualityScoreSlider)) {
      int qs = source.getValue();
      this.firePropertyChange("qualityScore", qualityScore, qs);
      qualityScore = qs;
      lblQualityScore.setText("Quality Score (" + qualityScore + ")");
    } else if (source.equals(rectangleHeightSlider)) {
      int height = source.getValue();
      this.firePropertyChange("rectangleHeight", rectangleHeight, height);
      rectangleHeight = height;
      lblRectangleHeight.setText("Rectangle Height (" + rectangleHeight + ")");
    }
  }

  // Monitor the display mode combobox for changes
  @Override
  public void actionPerformed(ActionEvent arg0) {
    @SuppressWarnings("unchecked")
    JComboBox cb = (JComboBox) arg0.getSource();
    String mode = (String) cb.getSelectedItem();
    firePropertyChange("displayMode", displayMode, mode);
    displayMode = mode;
    cnvPanel.setDisplayMode(displayMode);
  }

  public String getDisplayMode() {
    return displayMode;
  }

  public int getProbes() {
    return probes;
  }

  public int getMinSize() {
    return minSize;
  }

  public int getQualityScore() {
    return qualityScore;
  }

  public int getRectangleHeight() {
    return rectangleHeight;
  }

  public void setSelectedCNVs(List<CNVariant> cnvs) {
    cnvPanel.setCNVs(cnvs);
  }

  public CompPlot getPlot() {
    return compPlot;
  }

}


class CNVPanel extends JPanel implements ActionListener {
  private static final long serialVersionUID = 1L;
  CNVariant selectedCNV;
  JPanel cnvPanel;
  JLabel iid; // Individual ID
  JLabel fid; // Family ID
  JLabel length; // CNV length
  JLabel copies; // CNV copy number
  JLabel probes; // Number of probes
  JLabel score; // Quality score
  String displayMode;
  JScrollPane cnvPane;
  JScrollPane cnvScroll;

  JLabel cnvListLabel;
  CNVCheckList checkList;
  JButton selectAll;
  JButton selectNone;
  JButton trailerButton;
  LaunchAction launchTrailer;
  CompPlot compPlot;

  public CNVPanel(CompPlot cp) {
    compPlot = cp;
    setLayout(new BoxLayout(this, BoxLayout.PAGE_AXIS));
    add(new JLabel("Selected CNV:"));
    iid = new JLabel();
    fid = new JLabel();
    length = new JLabel();
    copies = new JLabel();
    probes = new JLabel();
    score = new JLabel();
    selectAll = new JButton("Select All");
    selectNone = new JButton("Select None");
    trailerButton = new JButton("To Trailer");
    trailerButton.setToolTipText("Launch Trailer to examine selected CNVs");

    // Panel for CNV information
    cnvPanel = new JPanel();
    cnvPanel.setLayout(new GridLayout(7, 2));
    cnvPanel.setPreferredSize(new Dimension(250, 120));
    cnvPanel.setMinimumSize(new Dimension(250, 120));
    cnvPanel.setMaximumSize(new Dimension(250, 120));

    // Individual ID
    cnvPanel.add(new JLabel("IID:"));
    cnvPanel.add(iid);

    // Family ID
    cnvPanel.add(new JLabel("FID:"));
    cnvPanel.add(fid);

    // Length
    cnvPanel.add(new JLabel("Length:"));
    cnvPanel.add(length);

    // Copies
    cnvPanel.add(new JLabel("Copies"));
    cnvPanel.add(copies);

    // Probes
    cnvPanel.add(new JLabel("Probes:"));
    cnvPanel.add(probes);

    // Score
    cnvPanel.add(new JLabel("Score:"));
    cnvPanel.add(score);

    cnvListLabel = new JLabel("");
    cnvPanel.add(cnvListLabel);

    add(cnvPanel);

    cnvScroll = new JScrollPane();
    add(cnvScroll);
    cnvScroll.setVisible(false);

    JPanel btnPanel = new JPanel();
    btnPanel.setLayout(new FlowLayout(FlowLayout.LEADING));
    btnPanel.add(selectAll);
    selectAll.setVisible(false);
    selectAll.addActionListener(this);

    btnPanel.add(selectNone);
    selectNone.setVisible(false);
    selectNone.addActionListener(this);

    // Link off to Trailer
    btnPanel.add(trailerButton);
    trailerButton.setEnabled(false);
    trailerButton.addActionListener(this);

    add(btnPanel);
  }

  /**
   * Update the list of <b>available</b> CNVs. If there is only one CNV or we are not in collapsed
   * mode, the single CNV will be selected. Otherwise, CNV selection will go through a
   * checkbox-style side panel.
   */
  public void setCNVs(List<CNVariant> cnvs) {
    selectedCNV = null;

    // In collapsed mode, if there are multiple CNVs associated, add a combo box that lets you
    // have granular control over which CNV to look at
    // We start with no individual CNVs selected
    if ("Collapsed".equals(displayMode) && cnvs.size() > 1) {
      cnvScroll.setPreferredSize(new Dimension(100, 100));
      JPanel checkPanel = new JPanel(new GridLayout(cnvs.size(), 1));
      checkList = new CNVCheckList(cnvs);
      for (JCheckBox checkBox : checkList.checkList) {
        checkPanel.add(checkBox);
      }
      cnvScroll.add(checkPanel);
      cnvScroll.setViewportView(checkPanel);
      setUIElements(true);
    } else {
      // Select the first CNV
      selectedCNV = cnvs.get(0);
      setCNVText(selectedCNV);
      setUIElements(false);
    }

    cnvPanel.repaint();
  }

  /**
   * Helper method to update UI elements that depend on whether or not we have selected a
   * CNVRectangle with more than one CNV.
   */
  private void setUIElements(boolean cnvPanelVisible) {
      cnvScroll.setVisible(cnvPanelVisible);
      cnvListLabel.setText(cnvPanelVisible ? "Select CNVs:" : "");
      selectAll.setVisible(cnvPanelVisible);
      selectNone.setVisible(cnvPanelVisible);
      trailerButton.setEnabled(!cnvPanelVisible);
  }

  // Update the text with the currently selected CNV
  public void setCNVText(CNVariant cnv) {
    iid.setText(cnv.getIndividualID());
    fid.setText(cnv.getFamilyID());
    length.setText("" + cnv.getSize());
    copies.setText("" + cnv.getCN());
    probes.setText("" + cnv.getNumMarkers());
    score.setText("" + cnv.getScore());
  }

  /**
   * Clear selected CNV text
   */
  private void clearCNVText() {
    iid.setText("");
    fid.setText("");
    length.setText("");
    copies.setText("");
    probes.setText("");
    score.setText("");
  }

  /**
   * Clear CNV UI components and text
   */
  private void clearCNVPanel() {
    trailerButton.setEnabled(false);
    selectAll.setVisible(false);
    selectNone.setVisible(false);
    cnvScroll.setVisible(false);
  }

  public void setDisplayMode(String mode) {
    displayMode = mode;
    clearCNVPanel();
    selectedCNV = null;
  }

  @Override
  public void actionPerformed(ActionEvent arg0) {
    if (arg0.getSource().equals(selectAll) && !checkList.checkList.isEmpty()) {
      checkList.selectAll();
      trailerButton.setEnabled(true);
    } else if (arg0.getSource().equals(selectNone)) {
      checkList.selectNone();
      trailerButton.setEnabled(false);
    } else if (arg0.getSource().equals(trailerButton)) {
      // if selectedCNV is not null, then we have a rect with one CNV and we use selectedCNV
      // if we have a rect with multiple CNVs, we use the checkList.
      List<CNVariant> selectedCNVs;
      if (selectedCNV != null) {
        selectedCNVs = new ArrayList<CNVariant>();
        selectedCNVs.add(selectedCNV);
      } else {
        selectedCNVs = checkList.getSelected();
      }

      compPlot.openTrailers(selectedCNVs);
    }
  }


  /**
   * Create a checklist of variants. Used when multiple CNVs are collapsed into a single rectangle.
   *
   * @author Michael Vieths
   * @author Mark Hiner
   *
   */
  private class CNVCheckList {
    List<JCheckBox> checkList;
    Map<JCheckBox, CNVariant> variantMap;
    int checkCount;

    public CNVCheckList(List<CNVariant> variants) {
      variantMap = new HashMap<JCheckBox, CNVariant>();
      checkList = new ArrayList<JCheckBox>();
      checkCount = 0;

      // This will control behavior when checking or unchecking an option
      // e.g., changing display text, enabling/disabling trailer button
      ActionListener actionListener = new ActionListener() {
        public void actionPerformed(ActionEvent actionEvent) {
          JCheckBox checkBox = (JCheckBox) actionEvent.getSource();
          if (checkBox.getModel().isSelected()) {
            checkCount++;
            setCNVText(variantMap.get(checkBox));
          } else {
            checkCount--;
            clearCNVText();
          }

          trailerButton.setEnabled(checkCount > 0);
        }
      };

      for (CNVariant cnv : variants) {
        JCheckBox box = new JCheckBox(cnv.getIndividualID());
        box.addActionListener(actionListener);
        checkList.add(box);
        variantMap.put(box, cnv);
      }
    }

    public void selectAll() {
      for (JCheckBox box : checkList) {
        box.setSelected(true);
      }
    }

    public void selectNone() {
      for (JCheckBox box : checkList) {
        box.setSelected(false);
      }
    }

    /**
     * @return List of all selected CNVs
     */
    public List<CNVariant> getSelected() {
      List<CNVariant> selected = new ArrayList<CNVariant>();
      for (int i = 0; i < checkList.size(); i++) {
        JCheckBox checkBox = checkList.get(i);
        if (checkBox.isSelected()) {
          selected.add(variantMap.get(checkBox));
        }
      }
      return selected;
    }
  }

}
