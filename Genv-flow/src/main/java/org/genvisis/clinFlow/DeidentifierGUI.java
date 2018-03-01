package org.genvisis.clinFlow;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTextField;
import javax.swing.SwingUtilities;
import javax.swing.border.Border;
import javax.swing.border.EmptyBorder;
import org.genvisis.cnv.gui.FileChooser;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import net.miginfocom.swing.MigLayout;

public class DeidentifierGUI extends JDialog {

  private final JPanel contentPanel = new JPanel();
  private JSplitPane splitP = new JSplitPane(JSplitPane.VERTICAL_SPLIT);
  private final JPanel progressPanel = new JPanel(new MigLayout("", "[grow]", "[grow]"));

  public static void main(String[] args) {
    try {
      DeidentifierGUI dialog = new DeidentifierGUI();
      dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
      dialog.setVisible(true);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  public DeidentifierGUI() {
    setTitle("Deidentifier");
    setBounds(100, 100, 550, 500);
    setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
    getContentPane().setLayout(new BorderLayout());
    getContentPane().add(splitP, BorderLayout.CENTER);

    splitP.setBorder(new EmptyBorder(5, 5, 5, 5));
    contentPanel.setBorder(new EmptyBorder(5, 5, 5, 5));
    progressPanel.setBorder(null);

    splitP.setLeftComponent(contentPanel);
    splitP.setRightComponent(progressPanel);

    setupContentPanel();

    {
      JPanel buttonPane = new JPanel();
      buttonPane.setLayout(new FlowLayout(FlowLayout.RIGHT));
      getContentPane().add(buttonPane, BorderLayout.SOUTH);
      {
        JButton okButton = new JButton();
        okButton.addActionListener((e) -> {
          run();
        });
        okButton.setText("Run");
        okButton.setActionCommand("Run");
        buttonPane.add(okButton);
        getRootPane().setDefaultButton(okButton);
      }
      {
        JButton cancelButton = new JButton();
        cancelButton.addActionListener((e) -> {
          setVisible(false);
          dispose();
        });
        cancelButton.setText("Close");
        cancelButton.setActionCommand("Close");
        buttonPane.add(cancelButton);
      }
    }
  }

  JTextField rootInFld;
  JTextField rootOutFld;
  JTextField linkFld;

  JButton selectInBtn;
  JButton selectOutBtn;
  JButton selectLinkBtn;

  private void setupContentPanel() {
    contentPanel.setLayout(new MigLayout("", "[grow][]", "[][][][][][][][grow]"));

    JLabel titleLbl = new JLabel("<html><u>FCS Deidentifier</u></html>");
    titleLbl.setBorder(null);
    titleLbl.setFont(new Font("Sanserif", Font.BOLD, 24));
    contentPanel.add(titleLbl, "cell 0 0, span 2");

    rootInFld = new JTextField();
    titleLbl = new JLabel("Input Directory:");
    contentPanel.add(titleLbl, "cell 0 1");
    contentPanel.add(rootInFld, "cell 0 2, growx");

    final Border b1 = BorderFactory.createCompoundBorder(BorderFactory.createLineBorder(Color.gray),
                                                         new EmptyBorder(2, 10, 1, 10));
    final Border b2 = BorderFactory.createCompoundBorder(BorderFactory.createCompoundBorder(BorderFactory.createLineBorder(new Color(Color.gray.getRed(),
                                                                                                                                     Color.gray.getGreen(),
                                                                                                                                     Color.gray.getBlue(),
                                                                                                                                     100),
                                                                                                                           1),
                                                                                            BorderFactory.createCompoundBorder(BorderFactory.createLineBorder(Color.gray,
                                                                                                                                                              1),
                                                                                                                               BorderFactory.createLineBorder(new Color(Color.gray.getRed(),
                                                                                                                                                                        Color.gray.getGreen(),
                                                                                                                                                                        Color.gray.getBlue(),
                                                                                                                                                                        100),
                                                                                                                                                              1))),
                                                         new EmptyBorder(0, 8, -1, 8));

    selectInBtn = new JButton();
    selectInBtn.addActionListener((e) -> {
      selectIn();
    });
    selectInBtn.addMouseListener(new MouseAdapter() {

      @Override
      public void mouseEntered(MouseEvent e) {
        super.mouseEntered(e);
        selectInBtn.setBorder(b2);
      }

      @Override
      public void mouseExited(MouseEvent e) {
        super.mouseExited(e);
        selectInBtn.setBorder(b1);
      }
    });
    selectInBtn.setBorder(b1);
    selectInBtn.setText(">");
    contentPanel.add(selectInBtn, "cell 1 2");

    titleLbl = new JLabel("Output Directory:");
    contentPanel.add(titleLbl, "cell 0 3");
    rootOutFld = new JTextField();
    contentPanel.add(rootOutFld, "cell 0 4, growx");

    selectOutBtn = new JButton();
    selectOutBtn.addMouseListener(new MouseAdapter() {

      @Override
      public void mouseEntered(MouseEvent e) {
        super.mouseEntered(e);
        selectOutBtn.setBorder(b2);
      }

      @Override
      public void mouseExited(MouseEvent e) {
        super.mouseExited(e);
        selectOutBtn.setBorder(b1);
      }
    });
    selectOutBtn.addActionListener((e) -> {
      selectOut();
    });
    selectOutBtn.setBorder(b1);
    selectOutBtn.setText(">");
    contentPanel.add(selectOutBtn, "cell 1 4");

    titleLbl = new JLabel("Link File Directory:");
    contentPanel.add(titleLbl, "cell 0 5");
    linkFld = new JTextField();
    contentPanel.add(linkFld, "cell 0 6, growx");

    selectLinkBtn = new JButton();
    selectLinkBtn.addMouseListener(new MouseAdapter() {

      @Override
      public void mouseEntered(MouseEvent e) {
        super.mouseEntered(e);
        selectLinkBtn.setBorder(b2);
      }

      @Override
      public void mouseExited(MouseEvent e) {
        super.mouseExited(e);
        selectLinkBtn.setBorder(b1);
      }
    });
    selectLinkBtn.addActionListener((e) -> {
      selectLink();
    });
    selectLinkBtn.setBorder(b1);
    selectLinkBtn.setText(">");
    contentPanel.add(selectLinkBtn, "cell 1 6");
  }

  private void selectOut() {
    FileChooser fc = new FileChooser(this, "", false, true, "Select Output Directory",
                                     new Logger());
    if (!fc.isSelected()) return;
    final String dir = fc.getSelectedFile().getAbsolutePath();
    rootOutFld.setText(ext.verifyDirFormat(dir));
  }

  private void selectIn() {
    FileChooser fc = new FileChooser(this, "", false, true, "Select Source Directory",
                                     new Logger());
    if (!fc.isSelected()) return;
    final String dir = fc.getSelectedFile().getAbsolutePath();
    rootInFld.setText(ext.verifyDirFormat(dir));
  }

  private void selectLink() {
    FileChooser fc = new FileChooser(this, "", false, true, "Select Link File Directory",
                                     new Logger());
    if (!fc.isSelected()) return;
    final String dir = fc.getSelectedFile().getAbsolutePath();
    linkFld.setText(ext.verifyDirFormat(dir));
  }

  private boolean checkAndWarn() {
    String in = rootInFld.getText();
    String out = rootOutFld.getText();

    if (in.isEmpty()) {
      JOptionPane.showMessageDialog(null, "Error - missing source directory", "Error",
                                    JOptionPane.ERROR_MESSAGE);
      return false;
    }
    if (!Files.exists(ext.verifyDirFormat(in))) {
      JOptionPane.showMessageDialog(null, "Error - source directory invalid / not found", "Error",
                                    JOptionPane.ERROR_MESSAGE);
      return false;
    }

    if (out.isEmpty()) {
      JOptionPane.showMessageDialog(null, "Error - missing destination directory", "Error",
                                    JOptionPane.ERROR_MESSAGE);
      return false;
    }
    if (!Files.exists(ext.verifyDirFormat(out))) {
      JOptionPane.showMessageDialog(null, "Error - destination directory invalid / not found",
                                    "Error", JOptionPane.ERROR_MESSAGE);
      return false;
    }

    String link = linkFld.getText();
    if (!link.isEmpty() && !Files.exists(ext.verifyDirFormat(link))) {
      JOptionPane.showMessageDialog(null, "Error - link file directory invalid / not found",
                                    "Error", JOptionPane.ERROR_MESSAGE);
      return false;
    }

    return true;
  }

  private void run() {
    if (!checkAndWarn()) {
      return;
    }
    String in = rootInFld.getText();
    String out = rootOutFld.getText();
    Deidentifier deid = new Deidentifier(in, out,
                                         linkFld.getText().isEmpty() ? out : linkFld.getText());
    final List<Conversion> converts = deid.identify();
    final Map<Conversion, JProgressBar> progs = new HashMap<>();

    SwingUtilities.invokeLater(new Runnable() {

      @Override
      public void run() {

        JPanel subPanel = new JPanel(new MigLayout("", "[][grow]", "[grow]"));
        JScrollPane pane = new JScrollPane(subPanel);
        pane.setBorder(null);
        progressPanel.removeAll();
        progressPanel.add(pane, "grow");

        final Set<Conversion> complete = new HashSet<>();
        for (int i = 0; i < converts.size(); i++) {
          JProgressBar pr = new JProgressBar();
          progs.put(converts.get(i), pr);
          subPanel.add(new JLabel(converts.get(i).fcs), "cell 0 " + (i * 2));
          subPanel.add(pr, "cell 0 " + ((i * 2) + 1) + ", span 2, grow");
          try {
            if (Deidentifier.exists(converts.get(i))) {
              pr.setStringPainted(true);
              pr.setString("Found Existing");
              complete.add(converts.get(i));
            }
          } catch (IOException e) {
            // just redo
          }
        }

        progressPanel.revalidate();

        new Thread(() -> {
          for (Conversion c : converts) {
            if (complete.contains(c)) continue;
            progs.get(c).setIndeterminate(true);
            progs.get(c).setStringPainted(true);
            progs.get(c).setString("Processing ... ");
            Deidentifier.processSingleFCS(c);
            progs.get(c).setIndeterminate(false);
            progs.get(c).setString("Complete");
          }
        }).start();
      }
    });

  }

}
