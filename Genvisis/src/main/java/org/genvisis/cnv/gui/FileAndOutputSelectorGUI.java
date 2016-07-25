package org.genvisis.cnv.gui;

import java.awt.BorderLayout;
import java.awt.FlowLayout;
import java.awt.Frame;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.util.ArrayList;

import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.border.EmptyBorder;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.swing.text.JTextComponent;

import net.miginfocom.swing.MigLayout;

import org.genvisis.common.Array;
import org.genvisis.common.ext;

public class FileAndOutputSelectorGUI extends JDialog {

    private final JPanel contentPanel = new JPanel();
    private JTextField txtFldOutputFile;
    private JComponent inputComponent;

    volatile int option = -1;
    
    public static String[] showFileAndOutputSelector(Frame parent, String inputSuggestion, int inMode, String[] inAllowedExts, String outputSuggestion, int outMode) {
        FileAndOutputSelectorGUI faosgui = new FileAndOutputSelectorGUI(parent, inputSuggestion, inMode, inAllowedExts, outputSuggestion, outMode);
        faosgui.setVisible(true);
        int opt = faosgui.option;
        String[] retVals = null;
        if (opt == JOptionPane.OK_OPTION) {
            String in = null;
            if (faosgui.inputComponent instanceof JTextComponent) {
                in = ((JTextComponent) faosgui.inputComponent).getText();
            } else if (faosgui.inputComponent instanceof JComboBox) {
                in = ((JComboBox) faosgui.inputComponent).getSelectedItem().toString();
            }
            retVals = new String[]{in, faosgui.txtFldOutputFile.getText()};
        }
        return retVals;
    }

    public static String[] showFileAndOutputSelector(Frame parent, String[] inputSuggestions, int inMode, String[] inAllowedExts, String outputSuggestion, int outMode, boolean showAsCombo) {
        FileAndOutputSelectorGUI faosgui = new FileAndOutputSelectorGUI(parent, inputSuggestions, inMode, inAllowedExts, outputSuggestion, outMode, showAsCombo);
        faosgui.setVisible(true);
        int opt = faosgui.option;
        String[] retVals = null;
        if (opt == JOptionPane.OK_OPTION) {
            String in = null;
            if (faosgui.inputComponent instanceof JTextComponent) {
                in = ((JTextComponent) faosgui.inputComponent).getText();
            } else if (faosgui.inputComponent instanceof JComboBox) {
                in = ((JComboBox) faosgui.inputComponent).getSelectedItem().toString();
            }
            retVals = new String[]{in, faosgui.txtFldOutputFile.getText()};
        }
        return retVals;
    }
    
    /**
     * 
     * @param mode Input file selection mode:   <br /><ul>
            <li>JFileChooser.FILES_ONLY </li>
            <li>JFileChooser.DIRECTORIES_ONLY </li>
            <li>JFileChooser.FILES_AND_DIRECTORIES </li></ul>

     */
    private FileAndOutputSelectorGUI(final Frame parent, final String inputSuggestion, final int inMode, final String[] inAllowedExts, final String outputSuggestion, final int outMode) {
        super(parent, true);
        setTitle("Select Input and Output Files");
        setBounds(100, 100, 460, 192);
        getContentPane().setLayout(new BorderLayout());
        contentPanel.setBorder(new EmptyBorder(5, 5, 5, 5));
        getContentPane().add(contentPanel, BorderLayout.CENTER);
        contentPanel.setLayout(new MigLayout("", "[grow]", "[][][][]"));
        
        Insets btnInsets = new Insets(0, 14, 0, 14);
        {
            JLabel lblInputFile = new JLabel("Input File:");
            contentPanel.add(lblInputFile, "cell 0 0");
        }
        {
            final JTextField txtFldInputFile = new JTextField(inputSuggestion == null ? "" : inputSuggestion);
            contentPanel.add(txtFldInputFile, "flowx,cell 0 1,growx");
            txtFldInputFile.setColumns(10);
            inputComponent = txtFldInputFile;
            JButton btnSelIn = new JButton(">");
            btnSelIn.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    String inPrep = txtFldInputFile.getText();
                    if (!inPrep.equals("")) {
                        inPrep = ext.parseDirectoryOfFile(inPrep);
                    }
                    JFileChooser jfc = new JFileChooser(inPrep);
                    jfc.setFileSelectionMode(inMode);
                    if (inAllowedExts != null && inAllowedExts.length > 0) {
                        jfc.setFileFilter(new FileNameExtensionFilter("Allowed Files", inAllowedExts));
                    }
                    int opt = jfc.showOpenDialog(FileAndOutputSelectorGUI.this);
                    if (opt == JFileChooser.APPROVE_OPTION) {
                        txtFldInputFile.setText(jfc.getSelectedFile().getAbsolutePath());
                    }
                }
            });
            btnSelIn.setMargin(btnInsets);
            contentPanel.add(btnSelIn, "cell 0 1");
        }
        {
            JLabel lblOutputFile = new JLabel("Output File:");
            contentPanel.add(lblOutputFile, "cell 0 2");
        }
        {
            txtFldOutputFile = new JTextField(outputSuggestion == null ? "" : outputSuggestion);
            contentPanel.add(txtFldOutputFile, "flowx,cell 0 3,growx");
            txtFldOutputFile.setColumns(10);
        }
        {
            JButton btnSelOut = new JButton(">");
            btnSelOut.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    JFileChooser jfc = new JFileChooser(ext.parseDirectoryOfFile(txtFldOutputFile.getText()));
                    jfc.setFileSelectionMode(outMode);
                    int opt = jfc.showSaveDialog(FileAndOutputSelectorGUI.this);
                    if (opt == JFileChooser.APPROVE_OPTION) {
                        txtFldOutputFile.setText(jfc.getSelectedFile().getAbsolutePath());
                    }
                }
            });
            btnSelOut.setMargin(btnInsets); 
            contentPanel.add(btnSelOut, "cell 0 3");
        }
        {
            JPanel buttonPane = new JPanel();
            buttonPane.setLayout(new FlowLayout(FlowLayout.RIGHT));
            getContentPane().add(buttonPane, BorderLayout.SOUTH);
            {
                JButton btnOk = new JButton("OK");
                btnOk.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent arg0) {
                        option = JOptionPane.OK_OPTION;
                        setVisible(false);
                    }
                });
                btnOk.setActionCommand("OK");
                buttonPane.add(btnOk);
                getRootPane().setDefaultButton(btnOk);
            }
            {
                JButton btnCancel = new JButton("Cancel");
                btnCancel.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        option = JOptionPane.CANCEL_OPTION;
                        setVisible(false);
                    }
                });
                btnCancel.setActionCommand("Cancel");
                buttonPane.add(btnCancel);
            }
        }
    }
    
    /**
     * 
     * @param mode Input file selection mode:   <br /><ul>
            <li>JFileChooser.FILES_ONLY </li>
            <li>JFileChooser.DIRECTORIES_ONLY </li>
            <li>JFileChooser.FILES_AND_DIRECTORIES </li></ul>

     */
    private FileAndOutputSelectorGUI(final Frame parent, final String[] inputSuggestions, final int inMode, final String[] inAllowedExts, final String outputSuggestion, final int outMode, boolean inSuggAsCombo) {
        super(parent, true);
        setTitle("Select Input and Output Files");
        setBounds(100, 100, 460, 192);
        getContentPane().setLayout(new BorderLayout());
        contentPanel.setBorder(new EmptyBorder(5, 5, 5, 5));
        getContentPane().add(contentPanel, BorderLayout.CENTER);
        contentPanel.setLayout(new MigLayout("", "[grow]", "[][][][]"));
        
        Insets btnInsets = new Insets(0, 14, 0, 14);
        {
            JLabel lblInputFile = new JLabel("Input File:");
            contentPanel.add(lblInputFile, "cell 0 0");
        }
        if (inSuggAsCombo && inputSuggestions != null) {
            JComboBox cmbInputFile = new JComboBox(inputSuggestions);
            contentPanel.add(cmbInputFile, "flowx,cell 0 1, growx");
            inputComponent = cmbInputFile;
        } else {
            {
                final JTextField txtFldInputFile = new JTextField(inputSuggestions == null ? "" : Array.toStr(inputSuggestions, ";"));
                contentPanel.add(txtFldInputFile, "flowx,cell 0 1,growx");
                txtFldInputFile.setColumns(10);
                inputComponent = txtFldInputFile;
                JButton btnSelIn = new JButton(">");
                btnSelIn.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        JFileChooser jfc = new JFileChooser();
                        ArrayList<File> selFiles = new ArrayList<File>();
                        for (String str : txtFldInputFile.getText().split(";")) {
                            selFiles.add(new File(str));
                        }
                        jfc.setSelectedFiles(selFiles.toArray(new File[selFiles.size()]));
                        jfc.setFileSelectionMode(inMode);
                        if (inAllowedExts != null && inAllowedExts.length > 0) {
                            jfc.setFileFilter(new FileNameExtensionFilter("Allowed Files", inAllowedExts));
                        }
                        int opt = jfc.showOpenDialog(FileAndOutputSelectorGUI.this);
                        if (opt == JFileChooser.APPROVE_OPTION) {
                            File[] selFilesArr = jfc.getSelectedFiles();
                            StringBuilder txt = new StringBuilder();
                            for (int i = 0; i < selFilesArr.length; i++) {
                                txt.append(ext.verifyDirFormat(selFilesArr[i].getAbsolutePath()));
                                if (i < selFilesArr.length - 1) {
                                    txt.append(";");
                                }
                            }
                            txtFldInputFile.setText(txt.toString());
                        }
                    }
                });
                btnSelIn.setMargin(btnInsets);
                contentPanel.add(btnSelIn, "cell 0 1");
            }
        }
        {
            JLabel lblOutputFile = new JLabel("Output File:");
            contentPanel.add(lblOutputFile, "cell 0 2");
        }
        {
            txtFldOutputFile = new JTextField(outputSuggestion == null ? "" : outputSuggestion);
            contentPanel.add(txtFldOutputFile, "flowx,cell 0 3,growx");
            txtFldOutputFile.setColumns(10);
        }
        {
            JButton btnSelOut = new JButton(">");
            btnSelOut.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    JFileChooser jfc = new JFileChooser(ext.parseDirectoryOfFile(txtFldOutputFile.getText()));
                    jfc.setFileSelectionMode(outMode);
                    int opt = jfc.showSaveDialog(FileAndOutputSelectorGUI.this);
                    if (opt == JFileChooser.APPROVE_OPTION) {
                        txtFldOutputFile.setText(jfc.getSelectedFile().getAbsolutePath());
                    }
                }
            });
            btnSelOut.setMargin(btnInsets);
            contentPanel.add(btnSelOut, "cell 0 3");
        }
        {
            JPanel buttonPane = new JPanel();
            buttonPane.setLayout(new FlowLayout(FlowLayout.RIGHT));
            getContentPane().add(buttonPane, BorderLayout.SOUTH);
            {
                JButton btnOk = new JButton("OK");
                btnOk.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent arg0) {
                        option = JOptionPane.OK_OPTION;
                        setVisible(false);
                    }
                });
                btnOk.setActionCommand("OK");
                buttonPane.add(btnOk);
                getRootPane().setDefaultButton(btnOk);
            }
            {
                JButton btnCancel = new JButton("Cancel");
                btnCancel.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        option = JOptionPane.CANCEL_OPTION;
                        setVisible(false);
                    }
                });
                btnCancel.setActionCommand("Cancel");
                buttonPane.add(btnCancel);
            }
        }
    }

}
