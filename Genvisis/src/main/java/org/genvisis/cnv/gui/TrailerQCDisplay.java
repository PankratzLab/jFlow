package org.genvisis.cnv.gui;

import java.awt.BorderLayout;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.ScrollPaneConstants;
import javax.swing.border.EmptyBorder;

import net.miginfocom.swing.MigLayout;

public class TrailerQCDisplay extends JDialog {

  private static final long serialVersionUID = 1L;

  /**
   * Launch the application.
   */
  public static void main(String[] args) {
    try {
      TrailerQCDisplay dialog = new TrailerQCDisplay();
      dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
      dialog.setVisible(true);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  private final JPanel contentPanel = new JPanel();

  String CHROM_LBL = "<html><u>Chromosome #:</u></html>";
  String REGION_LBL = "<html><u>Region #:</u></html>";
  private JLabel lblChromosomeNoGC;
  private JLabel lblRegionNoGC;
  private JLabel lblChromosomeGC;
  private JLabel lblRegionGC;
  private JLabel lblSampleID;

  /**
   * Create the dialog.
   */
  public TrailerQCDisplay() {
    setBounds(100, 100, 600, 400);
    getContentPane().setLayout(new BorderLayout());
    contentPanel.setBorder(new EmptyBorder(5, 5, 5, 5));
    getContentPane().add(contentPanel, BorderLayout.CENTER);
    contentPanel.setLayout(new BorderLayout(0, 0));
    {
      JScrollPane scrollPane = new JScrollPane();
      scrollPane.setVerticalScrollBarPolicy(ScrollPaneConstants.VERTICAL_SCROLLBAR_ALWAYS);
      contentPanel.add(scrollPane, BorderLayout.CENTER);
      {
        JPanel panel = new JPanel();
        scrollPane.setViewportView(panel);
        panel.setLayout(new MigLayout("", "[][grow]", "[][][][][][][][][][][][][][]"));
        {
          JLabel lblGcNaive = new JLabel("<html><u>GC Naive</u></html>");
          lblGcNaive.setFont(new Font("Tahoma", Font.PLAIN, 14));
          panel.add(lblGcNaive, "cell 0 0 2 1,alignx center");
        }
        {
          JLabel lblGenomeNoGC = new JLabel("<html><u>Genome:</u></html>");
          panel.add(lblGenomeNoGC, "cell 0 1,alignx right");
        }
        {
          lblChromosomeNoGC = new JLabel(CHROM_LBL);
          panel.add(lblChromosomeNoGC, "cell 0 3,alignx right");
        }
        {
          lblRegionNoGC = new JLabel(REGION_LBL);
          panel.add(lblRegionNoGC, "cell 0 5,alignx right");
        }
        {
          JLabel lblGcCorrected = new JLabel("<html><u>GC Corrected</u></html>");
          lblGcCorrected.setFont(new Font("Tahoma", Font.PLAIN, 14));
          panel.add(lblGcCorrected, "cell 0 7 2 1,alignx center");
        }
        {
          JLabel lblGenomeGC = new JLabel("<html><u>Genome:</u></html>");
          panel.add(lblGenomeGC, "cell 0 8,alignx right");
        }
        {
          lblChromosomeGC = new JLabel(CHROM_LBL);
          panel.add(lblChromosomeGC, "cell 0 10,alignx right");
        }
        {
          lblRegionGC = new JLabel(REGION_LBL);
          panel.add(lblRegionGC, "cell 0 12,alignx right");
        }
      }
    }
    {
      JPanel headerPanel = new JPanel();
      contentPanel.add(headerPanel, BorderLayout.NORTH);
      {
        lblSampleID = new JLabel("SAMPLE_ID");
        lblSampleID.setFont(new Font("Tahoma", Font.PLAIN, 16));
        headerPanel.add(lblSampleID);
      }
    }
    {
      JPanel buttonPane = new JPanel();
      buttonPane.setLayout(new FlowLayout(FlowLayout.RIGHT));
      getContentPane().add(buttonPane, BorderLayout.SOUTH);
      {
        JButton closeButton = new JButton("Close");
        closeButton.addActionListener(new ActionListener() {
          @Override
          public void actionPerformed(ActionEvent arg0) {
            TrailerQCDisplay.this.setVisible(false);
          }
        });
        closeButton.setActionCommand("Close");
        buttonPane.add(closeButton);
        getRootPane().setDefaultButton(closeButton);
      }
    }
  }

  public void setData(String sample, int chr, int start, int stop, String[] genNoGc,
      String[] chromNoGc, String[] regNoGC, String[] genGc, String[] chromGc, String[] regGC) {
    lblSampleID.setText(sample);
    String chrLbl = CHROM_LBL.replace("#", "" + chr);
    String regLbl = REGION_LBL.replace("#", "" + chr + ":" + start + "-" + stop);
    lblChromosomeNoGC.setText(chrLbl);
    lblChromosomeGC.setText(chrLbl);
    lblRegionNoGC.setText(regLbl);
    lblRegionGC.setText(regLbl);
    // 0 sampleID
    // 1 Array.mean(lrrs, true)
    // 2 Array.stdev(lrrs, true)
    // 3 lrrsdBound
    // 4 Array.stdev(bafs, true)
    // 5 (abCallRate > 0 ? abCallRate : forwardCallRate)
    // 6 (abCallRate > 0 ? abHetRate : forwardHetRate)
    // 7 wfPrior
    // 8 gcwfPrior
    // 9 wfPost
    // 10 gcwfPost
    // 11 lrrsdPost
    // 12 lrrsdPostBound
    // 13 multimodal
    // 14 Array.toStr(bafBinCounts)

  }


}
