package org.genvisis.one.ben.fcs.sub;

import java.awt.EventQueue;

import javax.swing.JFrame;

import net.miginfocom.swing.MigLayout;

import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JButton;
import javax.swing.JTextField;

import org.genvisis.CLI;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.Properties;

import javax.swing.JTextArea;
import javax.swing.JScrollPane;

public class CBCGUI {

  private JFrame frmCbcApplicator;
  private final JPanel panel = new JPanel();
  private JTextField txtCBCDir;
  private JTextField txtDataDir;
  private JTextField txtOutDir;

  /**
   * Create the application.
   */
  public CBCGUI() {
    initialize();
  }

  /**
   * Initialize the contents of the frame.
   */
  private void initialize() {
    frmCbcApplicator = new JFrame();
    frmCbcApplicator.setTitle("CBC Applicator");
    frmCbcApplicator.setBounds(100, 100, 450, 387);
    frmCbcApplicator.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    frmCbcApplicator.getContentPane().setLayout(
        new MigLayout("", "[18.00][][grow][][18]", "[][][][][][][][][grow][][]"));

    JLabel lblCbcDirectory = new JLabel("CBC Directory:");
    frmCbcApplicator.getContentPane().add(lblCbcDirectory, "cell 1 1,alignx trailing");

    txtCBCDir = new JTextField();
    frmCbcApplicator.getContentPane().add(txtCBCDir, "cell 2 1,growx");
    txtCBCDir.setColumns(10);

    JButton btnSelCBC = new JButton(">");
    btnSelCBC.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
        String curr = txtCBCDir.getText();
        if (curr.equals("")) {
          curr = "./";
        }
        JFileChooser jfc = new JFileChooser(curr);
        jfc.setFileSelectionMode(JFileChooser.FILES_ONLY);
        jfc.setDialogTitle("Select CBC Directory");
        jfc.setMultiSelectionEnabled(false);
        int resp = jfc.showOpenDialog(CBCGUI.this.frmCbcApplicator);
        if (resp == JFileChooser.APPROVE_OPTION) {
          String newPath = jfc.getSelectedFile().getAbsolutePath();
          txtCBCDir.setText(newPath);
          saveProps();
        }
      }
    });
    frmCbcApplicator.getContentPane().add(btnSelCBC, "cell 3 1");

    JLabel lblDataDirectory = new JLabel("Data Directory:");
    frmCbcApplicator.getContentPane().add(lblDataDirectory, "cell 1 3,alignx trailing");

    txtDataDir = new JTextField();
    frmCbcApplicator.getContentPane().add(txtDataDir, "cell 2 3,growx");
    txtDataDir.setColumns(10);

    JButton btnSelData = new JButton(">");
    btnSelData.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
        String curr = txtDataDir.getText();
        if (curr.equals("")) {
          curr = "./";
        }
        JFileChooser jfc = new JFileChooser(curr);
        jfc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
        jfc.setDialogTitle("Select Data File Directory");
        jfc.setMultiSelectionEnabled(false);
        int resp = jfc.showOpenDialog(CBCGUI.this.frmCbcApplicator);
        if (resp == JFileChooser.APPROVE_OPTION) {
          String newPath = jfc.getSelectedFile().getAbsolutePath();
          txtDataDir.setText(newPath);
          saveProps();
        }
      }
    });
    frmCbcApplicator.getContentPane().add(btnSelData, "cell 3 3");

    JLabel lblOutputDirectory = new JLabel("Output Directory:");
    frmCbcApplicator.getContentPane().add(lblOutputDirectory, "cell 1 5,alignx trailing");

    txtOutDir = new JTextField();
    frmCbcApplicator.getContentPane().add(txtOutDir, "cell 2 5,growx");
    txtOutDir.setColumns(10);

    JButton btnSelOut = new JButton(">");
    btnSelOut.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
        String curr = txtOutDir.getText();
        if (curr.equals("")) {
          curr = "./";
        }
        JFileChooser jfc = new JFileChooser(curr);
        jfc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
        jfc.setDialogTitle("Select Output Directory");
        jfc.setMultiSelectionEnabled(false);
        int resp = jfc.showOpenDialog(CBCGUI.this.frmCbcApplicator);
        if (resp == JFileChooser.APPROVE_OPTION) {
          String newPath = jfc.getSelectedFile().getAbsolutePath();
          txtOutDir.setText(newPath);
          saveProps();
        }
      }
    });
    frmCbcApplicator.getContentPane().add(btnSelOut, "cell 3 5");

    JLabel lblLog = new JLabel("Log:");
    frmCbcApplicator.getContentPane().add(lblLog, "cell 1 7,alignx left");

    JScrollPane scrollPane = new JScrollPane();
    frmCbcApplicator.getContentPane().add(scrollPane, "cell 1 8 3 2,grow");

    JTextArea textArea = new JTextArea();
    scrollPane.setViewportView(textArea);
    frmCbcApplicator.getContentPane().add(panel, "south,growx");
    panel.setLayout(new MigLayout("", "[grow][75px][75px]", "[23px]"));

    JButton btnRun = new JButton("Run");
    btnRun.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent arg0) {
        CBCApplicator cbcA = new CBCApplicator();
        cbcA.setCBCDir(txtCBCDir.getText());
        cbcA.setDataDir(txtDataDir.getText());
        cbcA.setOutputDirectory(txtOutDir.getText());
        Logger log = new Logger() {
          @Override
          public void reportTimeError(String str) {
            // override to show message dialog and remove genvisis version number from log output
            // super.reportTimeError(str);
            reportError(ext.getTime() + "]\t Error - " + str, true, true);
            JOptionPane.showMessageDialog(CBCGUI.this.frmCbcApplicator, str, "Error!",
                JOptionPane.ERROR_MESSAGE);
            throw new RuntimeException(str);
          }
        };
        log.linkTextArea(textArea);
        cbcA.setLog(log);
        try {
          cbcA.run();
        } catch (RuntimeException e) {
          return;
        }
        JOptionPane.showMessageDialog(CBCGUI.this.frmCbcApplicator, "Done!", "Done!",
            JOptionPane.INFORMATION_MESSAGE);
      }
    });
    panel.add(btnRun, "cell 1 0,growx,aligny top");

    JButton btnClose = new JButton("Close");
    btnClose.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
        CBCGUI.this.frmCbcApplicator.setVisible(false);
        CBCGUI.this.frmCbcApplicator.dispose();
      }
    });
    panel.add(btnClose, "cell 2 0,growx,aligny top");

    loadProps();
  }

  private static final String PROP_FILE = "cbcgui.properties";
  private static final String PROPKEY_CBC_DIR = "CBC_DIR";
  private static final String PROPKEY_DATA_DIR = "DATA_DIR";
  private static final String PROPKEY_OUT_DIR = "OUT_DIR";

  private void saveProps() {
    try {
      Properties props = new Properties();
      props.setProperty(PROPKEY_CBC_DIR, txtCBCDir.getText());
      props.setProperty(PROPKEY_DATA_DIR, txtDataDir.getText());
      props.setProperty(PROPKEY_OUT_DIR, txtOutDir.getText());
      File f = new File(PROP_FILE);
      OutputStream out = new FileOutputStream(f);
      props.store(out, "");
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  private void loadProps() {
    Properties props = new Properties();
    InputStream is = null;

    try {
      File f = new File(PROP_FILE);
      if (!f.exists())
        return;
      is = new FileInputStream(f);
      props.load(is);
      String cbc = props.getProperty(PROPKEY_CBC_DIR, "");
      String data = props.getProperty(PROPKEY_DATA_DIR, "");
      String out = props.getProperty(PROPKEY_OUT_DIR, "");

      if (!cbc.equals("")) {
        txtCBCDir.setText(cbc);
      }
      if (!data.equals("")) {
        txtDataDir.setText(data);
      }
      if (!out.equals("")) {
        txtOutDir.setText(out);
      }
    } catch (Exception e) {
      e.printStackTrace();
      is = null;
    }
  }

  /**
   * Launch the application.
   */
  public static void main(String[] args) {
    if (args.length == 0) {
      EventQueue.invokeLater(new Runnable() {
        public void run() {
          try {
            CBCGUI window = new CBCGUI();
            window.frmCbcApplicator.setVisible(true);
          } catch (Exception e) {
            e.printStackTrace();
          }
        }
      });
    } else {
      
      final String CBC = "cbcDir";
      final String DATA = "dataDir";
      final String OUT = "outDir";
      
      CLI cli = new CLI();
      cli.addArg(CBC, "CBC Directory", true);
      cli.addArg(DATA, "Datafile Directory", true);
      cli.addArg(OUT, "Output directory", false);
      
      cli.parseWithExit("CBCGUI", args);
      
      CBCApplicator cbcA = new CBCApplicator();
      cbcA.setCBCDir(cli.get(CBC));
      cbcA.setDataDir(cli.get(DATA));
      if (cli.has(OUT)) {
        cbcA.setOutputDirectory(cli.get(OUT));
      }
      cbcA.run();
      
    }
  }

}
