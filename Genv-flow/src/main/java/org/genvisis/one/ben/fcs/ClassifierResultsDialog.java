package org.genvisis.one.ben.fcs;

import java.awt.BorderLayout;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSeparator;
import javax.swing.SwingUtilities;
import javax.swing.border.EmptyBorder;
import net.miginfocom.swing.MigLayout;

public class ClassifierResultsDialog extends JDialog {

  private final JPanel contentPanel = new JPanel();
  private JLabel lblTNCnt;
  private JLabel lblFNCnt;
  private JLabel lblFPCnt;
  private JLabel lblTPCnt;

  /**
   * Launch the application.
   */
  public static void main(String[] args) {
    try {
      ClassifierResultsDialog dialog = new ClassifierResultsDialog();
      dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
      dialog.setVisible(true);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  /**
   * Create the dialog.
   */
  public ClassifierResultsDialog() {
    setBounds(100, 100, 200, 200);
    getContentPane().setLayout(new BorderLayout());
    contentPanel.setBorder(new EmptyBorder(5, 5, 5, 5));
    getContentPane().add(contentPanel, BorderLayout.CENTER);
    contentPanel.setLayout(new MigLayout("", "[5.00][][][grow]", "[][][][][][][]"));
    {
      JLabel lblClassifications = new JLabel("Classifications");
      lblClassifications.setFont(new Font("Arial", Font.PLAIN, 14));
      contentPanel.add(lblClassifications, "cell 0 0 4 1,alignx center");
    }
    {
      JSeparator separator = new JSeparator();
      contentPanel.add(separator, "cell 0 1 4 1,growx");
    }
    {
      JLabel lblTruePositives = new JLabel("True Positives:");
      contentPanel.add(lblTruePositives, "cell 1 2,alignx right");
    }
    {
      lblTPCnt = new JLabel("0");
      contentPanel.add(lblTPCnt, "cell 3 2");
    }
    {
      JLabel lblFalsePositivess = new JLabel("False Positives:");
      contentPanel.add(lblFalsePositivess, "cell 1 3,alignx right");
    }
    {
      lblFPCnt = new JLabel("0");
      contentPanel.add(lblFPCnt, "cell 3 3");
    }
    {
      JLabel lblTrueNegatives = new JLabel("True Negatives:");
      contentPanel.add(lblTrueNegatives, "cell 1 4,alignx right");
    }
    {
      lblTNCnt = new JLabel("0");
      contentPanel.add(lblTNCnt, "cell 3 4");
    }
    {
      JLabel lblFalseNegatives = new JLabel("False Negatives:");
      contentPanel.add(lblFalseNegatives, "cell 1 5,alignx right");
    }
    {
      lblFNCnt = new JLabel("0");
      contentPanel.add(lblFNCnt, "cell 3 5");
    }
    {
      JSeparator separator = new JSeparator();
      contentPanel.add(separator, "cell 0 6 4 1,growx");
    }
    {
      JPanel buttonPane = new JPanel();
      buttonPane.setLayout(new FlowLayout(FlowLayout.RIGHT));
      getContentPane().add(buttonPane, BorderLayout.SOUTH);
      {
        JButton cancelButton = new JButton("Close");
        cancelButton.addActionListener(new ActionListener() {

          public void actionPerformed(ActionEvent arg0) {
            ClassifierResultsDialog.this.setVisible(false);
            ClassifierResultsDialog.this.dispose();
          }
        });
        cancelButton.setActionCommand("Cancel");
        buttonPane.add(cancelButton);
      }
    }
  }

  public ClassifierResultsDialog setTNCnt(int tn) {
    SwingUtilities.invokeLater(new Runnable() {

      @Override
      public void run() {
        lblTNCnt.setText(Integer.toString(tn));
        repaint();
      }
    });
    return this;
  }

  public ClassifierResultsDialog setFNCnt(int fn) {
    SwingUtilities.invokeLater(new Runnable() {

      @Override
      public void run() {
        lblFNCnt.setText(Integer.toString(fn));
        repaint();
      }
    });
    return this;
  }

  public ClassifierResultsDialog setTPCnt(int tp) {
    SwingUtilities.invokeLater(new Runnable() {

      @Override
      public void run() {
        lblTPCnt.setText(Integer.toString(tp));
        repaint();
      }
    });
    return this;
  }

  public ClassifierResultsDialog setFPCnt(int fp) {
    SwingUtilities.invokeLater(new Runnable() {

      @Override
      public void run() {
        lblFPCnt.setText(Integer.toString(fp));
        repaint();
      }
    });
    return this;
  }

}
