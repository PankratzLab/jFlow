package org.genvisis.cnv.gui;

import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JPanel;
import javax.swing.JTextField;

import org.genvisis.common.Grafik;
import org.genvisis.common.Positions;
import org.genvisis.common.ext;

/**
 * Utility class for creating UI navigation bars for moving between chromosomal locations. To use
 * this class, use the {@link #RegionNavigator(ActionListener)} constructor with an
 * {@link ActionListener} implemented to check for constants of interest (e.g. {@link #FIRST_CHR}
 * indicates the user requested navigation to the first chromosomal position).
 * <p>
 * The returned {@link RegionNavigator} instance it itself a panel of buttons + text entry. The current
 * text value can be interrogated via the accessor API for this class.
 * </p>
 */
public class RegionNavigator extends JPanel {
  private static final long serialVersionUID = 1L;

  public static final String FIRST_CHR = "First chr";
  public static final String PREVIOUS_CHR = "Previous chr";
  public static final String NEXT_CHR = "Next chr";
  public static final String LAST_CHR = "Last chr";
  public static final String NAV_CHR = "Chr navigation field";
  private static final Font FONT = new Font("Arial", 0, 14);
  private static final int LINE_HEIGHT = 20;
  private static final String[] PREV_BTN = {"images/firstLast/Left.gif", "images/firstLast/dLeft.gif"};
  private static final String[] NEXT_BTN = {"images/firstLast/Right.gif", "images/firstLast/dRight.gif"};
  private static final String[] FIRST_BTN = {"images/firstLast/First.gif", "images/firstLast/dFirst.gif"};
  private static final String[] LAST_BTN = {"images/firstLast/Last.gif", "images/firstLast/dLast.gif"};

  private JTextField chrNavField = new JTextField("", LINE_HEIGHT);

  private ChrNavigator chrNav;

  public static final String DEFAULT_LOCATION = "chr1";
//  public static final String DEFAULT_LOCATION = "chr6:161,624,000-163,776,000"; // PARK2 region

  /**
   * Create the panel.
   */
  public RegionNavigator(ChrNavigator cNav) {
    ((FlowLayout)getLayout()).setVgap(0);
    chrNav = cNav;

    ActionListener al = new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent ae) {
        final String command = ae.getActionCommand();

        if (command.equals(RegionNavigator.FIRST_CHR)) {
          chrNav.setPosition("chr1");
        } else if (command.equals(RegionNavigator.PREVIOUS_CHR)) {
          chrNav.setPosition("chr" + Math.max(getChr() - 1, 1));
        } else if (command.equals(RegionNavigator.NEXT_CHR)) {
          chrNav.setPosition("chr" + Math.min(getChr() + 1, 26));
        } else if (command.equals(RegionNavigator.LAST_CHR)) {
          chrNav.setPosition("chr26");
        } else if (command.equals(RegionNavigator.NAV_CHR)) {
          chrNav.setPosition(getChrText());
        }
      }
    };

    // Chromosome navigation
    addButton(al, FIRST_CHR, "Go to first chromosome", FIRST_BTN);
    addButton(al, PREVIOUS_CHR, "Go to previous chromosome", PREV_BTN);
    addTextField(chrNavField, al, NAV_CHR);
    addButton(al, NEXT_CHR, "Go to next chromosome", NEXT_BTN);
    addButton(al, LAST_CHR, "Go to last chromosome", LAST_BTN);
  }

  public void setChrFieldText(byte chr, int start, int stop) {
    chrNavField.setText(chr == 0 ? "all"
                                       : "chr" + chr + ":" + ext.addCommas(start) + "-"
                                         + ext.addCommas(stop));
  }

  /**
   * Get current contents of the chromosome navfield
   */
  public String getChrText() {
    return chrNavField.getText().trim();
  }

  public int getChr() {
    return Positions.parseUCSClocation(chrNavField.getText().trim())[0];
  }

  private void addButton(ActionListener al, String command, String toolTip, String[] icons) {
    JButton btn = new JButton(Grafik.getImageIcon(icons[0]));
    btn.setDisabledIcon(Grafik.getImageIcon(icons[1]));
    btn.addActionListener(al);
    btn.setActionCommand(command);
    btn.setPreferredSize(new Dimension(LINE_HEIGHT, LINE_HEIGHT));
    btn.setToolTipText(toolTip); 
    add(btn);
  }

  private void addTextField(JTextField navField, ActionListener al, String command) {
    navField.setHorizontalAlignment(JTextField.CENTER);
    navField.setActionCommand(command);
    navField.setFont(FONT);
    navField.addActionListener(al);
    add(navField);
  }

  /**
   * Marker interface for plots that allow direct navigation between chromosomal positions.
   */
  public static interface ChrNavigator {
    /**
     * Update the position of this plot to the specified location.
     *
     * @param chr A string representation of chromosome location, according to
     *        {@link Positions#parseUCSClocation(String)}.
     */
    void setPosition(String chr);
  }
}
