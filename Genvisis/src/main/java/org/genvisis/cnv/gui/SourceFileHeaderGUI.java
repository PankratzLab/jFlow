package org.genvisis.cnv.gui;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.EventQueue;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSeparator;
import javax.swing.border.EmptyBorder;
import org.genvisis.cnv.filesys.SourceFileHeaderData;
import org.genvisis.common.ArrayUtils;
import net.miginfocom.swing.MigLayout;

public class SourceFileHeaderGUI extends JDialog {

  private static final long serialVersionUID = 1L;
  private final JPanel contentPane;
  private final JComboBox cbSnpInd;
  private final JComboBox cbSampInd;
  private final JComboBox cbA1Geno;
  private final JComboBox cbA2Geno;
  private final JComboBox cbA1AB;
  private final JComboBox cbA2AB;
  private final JComboBox cbX;
  private final JComboBox cbY;
  private final JComboBox cbBAF;
  private final JComboBox cbLRR;
  private final JComboBox cbGC;
  private final JComboBox cbXRaw;
  private final JComboBox cbYRaw;
  private final JComboBox cbR;
  private final JComboBox cbTheta;

  private boolean cancelled = false;

  private static String MISSING_STR = "[missing]";
  private static String FILENAME_STR = "[filename]";

  public static int MISSING_IND = -1;
  public static int FILENAME_IND = -99;

  /**
   * Launch the application.
   */
  public static void main(String[] args) {
    EventQueue.invokeLater(new Runnable() {

      @Override
      public void run() {
        try {
          // FIXME this is guaranteed to cause an NPE.
          SourceFileHeaderGUI frame = new SourceFileHeaderGUI(null);
          frame.setVisible(true);
        } catch (Exception e) {
          e.printStackTrace();
        }
      }
    });
  }

  private void doClose(boolean cancel) {
    cancelled = cancel;
    setVisible(false);
  }

  public boolean wasCancelled() {
    return cancelled;
  }

  /**
   * Create the frame.
   *
   * @param reportHdr
   */
  public SourceFileHeaderGUI(SourceFileHeaderData reportHdr) {
    setTitle("Assign Final Report Data Columns");
    setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
    setMinimumSize(new Dimension(100, 100));
    UITools.setSize(this, new Dimension(450, 620));
    contentPane = new JPanel();
    contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
    setContentPane(contentPane);
    contentPane.setLayout(new BorderLayout(0, 0));

    JPanel panel = new JPanel();
    contentPane.add(panel, BorderLayout.CENTER);
    panel.setLayout(new MigLayout("", "[grow][][grow 200][grow]",
                                  "[][][][][][][][][][][][][][][][][][][][][][grow]"));

    JLabel lblRequired = new JLabel("Required:");
    Font fontHeaderLabels = new Font("Tahoma", Font.BOLD, 12);
    lblRequired.setFont(fontHeaderLabels);
    panel.add(lblRequired, "pad 0 -15 0 0,cell 1 1");

    JLabel lblSnpIndex = new JLabel("SNP Index:");
    panel.add(lblSnpIndex, "cell 1 2");

    String[] headerParts = reportHdr.getCols();
    String[] headerOptions = ArrayUtils.addStrToArray(MISSING_STR, headerParts, 0);
    String[] sampOptions = ArrayUtils.addStrToArray(FILENAME_STR, headerParts, 0);

    cbSnpInd = new JComboBox(headerParts);
    cbSnpInd.setSelectedIndex(reportHdr.getColSnpIdent() == -1 ? 0 : reportHdr.getColSnpIdent()); // TODO SNP
    // INDEX
    // shouldn't
    // ever be
    // -1
    panel.add(cbSnpInd, "cell 2 2,growx");

    JLabel lblSampleId = new JLabel("Sample ID:");
    panel.add(lblSampleId, "cell 1 3");

    cbSampInd = new JComboBox(sampOptions);
    cbSampInd.setSelectedIndex(reportHdr.getColSampleIdent() == -1 ? 0
                                                                   : 1
                                                                     + reportHdr.getColSampleIdent());
    panel.add(cbSampInd, "cell 2 3,growx");

    JLabel lblRecommended = new JLabel("Recommended:");
    lblRecommended.setFont(fontHeaderLabels);
    panel.add(lblRecommended, "pad 0 -15 0 0,cell 1 4");

    JLabel lblAllele_2 = new JLabel("Allele 1 - Genotype:");
    panel.add(lblAllele_2, "cell 1 5");

    cbA1Geno = new JComboBox(headerOptions);
    cbA1Geno.setSelectedIndex(reportHdr.getColGeno1() == -1 ? 0 : 1 + reportHdr.getColGeno1());
    panel.add(cbA1Geno, "cell 2 5,growx");

    JLabel lblAllele_3 = new JLabel("Allele 2 - Genotype:");
    panel.add(lblAllele_3, "cell 1 6");

    cbA2Geno = new JComboBox(headerOptions);
    cbA2Geno.setSelectedIndex(reportHdr.getColGeno2() == -1 ? 0 : 1 + reportHdr.getColGeno2());
    panel.add(cbA2Geno, "cell 2 6,growx");

    JLabel lblAllele = new JLabel("Allele 1 - AB:");
    panel.add(lblAllele, "cell 1 7");

    cbA1AB = new JComboBox(headerOptions);
    cbA1AB.setSelectedIndex(reportHdr.getColGenoAB1() == -1 ? 0 : 1 + reportHdr.getColGenoAB1());
    panel.add(cbA1AB, "cell 2 7,growx");

    JLabel lblAllele_1 = new JLabel("Allele 2 - AB:");
    panel.add(lblAllele_1, "cell 1 8");

    cbA2AB = new JComboBox(headerOptions);
    cbA2AB.setSelectedIndex(reportHdr.getColGenoAB2() == -1 ? 0 : 1 + reportHdr.getColGenoAB2());
    panel.add(cbA2AB, "cell 2 8,growx");

    JLabel lblX = new JLabel("X:");
    panel.add(lblX, "cell 1 9");

    cbX = new JComboBox(headerOptions);
    cbX.setSelectedIndex(reportHdr.getColX() == -1 ? 0 : 1 + reportHdr.getColX());
    panel.add(cbX, "cell 2 9,growx");

    JLabel lblY = new JLabel("Y:");
    panel.add(lblY, "cell 1 10");

    cbY = new JComboBox(headerOptions);
    cbY.setSelectedIndex(reportHdr.getColY() == -1 ? 0 : 1 + reportHdr.getColY());
    panel.add(cbY, "cell 2 10,growx");

    JLabel lblBAlleleFreq = new JLabel("B Allele Freq.:");
    panel.add(lblBAlleleFreq, "cell 1 11");

    cbBAF = new JComboBox(headerOptions);
    cbBAF.setSelectedIndex(reportHdr.getColBAF() == -1 ? 0 : 1 + reportHdr.getColBAF());
    panel.add(cbBAF, "cell 2 11,growx");

    JLabel lblLogrRatio = new JLabel("Log-R Ratio:");
    panel.add(lblLogrRatio, "cell 1 12");

    cbLRR = new JComboBox(headerOptions);
    cbLRR.setSelectedIndex(reportHdr.getColLRR() == -1 ? 0 : 1 + reportHdr.getColLRR());
    panel.add(cbLRR, "cell 2 12,growx");

    JLabel lblConfidence = new JLabel("Confidence:");
    panel.add(lblConfidence, "cell 1 13");

    cbGC = new JComboBox(headerOptions);
    cbGC.setSelectedIndex(reportHdr.getColGC() == -1 ? 0 : 1 + reportHdr.getColGC());
    panel.add(cbGC, "cell 2 13,growx");

    JLabel lblOptional = new JLabel("Optional:");
    lblOptional.setFont(fontHeaderLabels);
    panel.add(lblOptional, "pad 0 -15 0 0,cell 1 14");

    JLabel lblXRaw = new JLabel("X Raw:");
    panel.add(lblXRaw, "cell 1 15");

    cbXRaw = new JComboBox(headerOptions);
    cbXRaw.setSelectedIndex(reportHdr.getColXRaw() == -1 ? 0 : 1 + reportHdr.getColXRaw());
    panel.add(cbXRaw, "cell 2 15,growx");

    JLabel lblYRaw = new JLabel("Y Raw:");
    panel.add(lblYRaw, "cell 1 16");

    cbYRaw = new JComboBox(headerOptions);
    cbYRaw.setSelectedIndex(reportHdr.getColYRaw() == -1 ? 0 : 1 + reportHdr.getColYRaw());
    panel.add(cbYRaw, "cell 2 16,growx");

    JLabel lblR = new JLabel("R:");
    panel.add(lblR, "cell 1 17");

    cbR = new JComboBox(headerOptions);
    cbR.setSelectedIndex(reportHdr.getColR() == -1 ? 0 : 1 + reportHdr.getColR());
    panel.add(cbR, "cell 2 17,growx");

    JLabel lblTheta = new JLabel("Theta:");
    panel.add(lblTheta, "cell 1 18");

    cbTheta = new JComboBox(headerOptions);
    cbTheta.setSelectedIndex(reportHdr.getColTheta() == -1 ? 0 : 1 + reportHdr.getColTheta());
    panel.add(cbTheta, "cell 2 18,growx");

    JSeparator separator = new JSeparator();
    panel.add(separator, "south");

    JPanel panel_1 = new JPanel();
    FlowLayout flowLayout = (FlowLayout) panel_1.getLayout();
    flowLayout.setAlignment(FlowLayout.TRAILING);
    contentPane.add(panel_1, BorderLayout.SOUTH);

    JButton btnOk = new JButton("OK");
    btnOk.addActionListener(new ActionListener() {

      @Override
      public void actionPerformed(ActionEvent e) {
        doClose(false);
      }
    });
    btnOk.setMnemonic(KeyEvent.VK_O);
    panel_1.add(btnOk);

    JButton btnCancel = new JButton("Cancel");
    btnCancel.addActionListener(new ActionListener() {

      @Override
      public void actionPerformed(ActionEvent e) {
        doClose(true);
      }
    });
    btnCancel.setMnemonic(KeyEvent.VK_C);
    panel_1.add(btnCancel);

    pack();
  }

  public int getSelectedSNPIndex() {
    return cbSnpInd.getSelectedIndex();
  }

  public int getSelectedSampleID() {
    return cbSampInd.getSelectedIndex() == 0 ? FILENAME_IND : cbSampInd.getSelectedIndex() - 1;
  }

  private int getValue(JComboBox cb) {
    return MISSING_STR.equals(cb.getSelectedItem()) ? MISSING_IND : cb.getSelectedIndex() - 1;
  }

  public int getSelectedGeno1() {
    return getValue(cbA1Geno);
  }

  public int getSelectedGeno2() {
    return getValue(cbA2Geno);
  }

  public int getSelectedAB1() {
    return getValue(cbA1AB);
  }

  public int getSelectedAB2() {
    return getValue(cbA2AB);
  }

  public int getSelectedBAF() {
    return getValue(cbBAF);
  }

  public int getSelectedLRR() {
    return getValue(cbLRR);
  }

  public int getSelectedGC() {
    return getValue(cbGC);
  }

  public int getSelectedR() {
    return getValue(cbR);
  }

  public int getSelectedTheta() {
    return getValue(cbTheta);
  }

  public int getSelectedX() {
    return getValue(cbX);
  }

  public int getSelectedY() {
    return getValue(cbY);
  }

  public int getSelectedXRaw() {
    return getValue(cbXRaw);
  }

  public int getSelectedYRaw() {
    return getValue(cbYRaw);
  }
}
