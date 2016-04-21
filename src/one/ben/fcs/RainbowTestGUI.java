package one.ben.fcs;

import java.awt.BorderLayout;
import java.awt.EventQueue;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.border.EmptyBorder;
import javax.swing.table.DefaultTableColumnModel;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.JTableHeader;

import net.miginfocom.swing.MigLayout;

import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JTextField;
import javax.swing.JButton;
import javax.swing.JTable;
import javax.swing.JScrollPane;
import javax.swing.ScrollPaneConstants;

import one.ben.fcs.FCSDataLoader.DATA_SET;
import scala.collection.mutable.HashMap;
import common.Array;
import common.ext;

import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.TreeSet;
import java.util.Vector;

public class RainbowTestGUI extends JFrame {

    private JPanel contentPane;
    private JTextField txtFldFileDir;
    private JTextField txtFldGatingFile;
    private JButton btnFileDirSelect;
    private JButton btnGatingFileSelect;
    private JTable table;
    private JScrollPane scrollPane;

    private HashMap<String, HashMap<String, Float>> fileMeans = new HashMap<String, HashMap<String,Float>>();
    private HashMap<String, HashMap<String, Float>> fileSDs = new HashMap<String, HashMap<String,Float>>();
    
    /**
     * Launch the application.
     */
    public static void main(String[] args) {
        EventQueue.invokeLater(new Runnable() {
            public void run() {
                try {
                    RainbowTestGUI frame = new RainbowTestGUI();
                    frame.setVisible(true);
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        });
    }

    /**
     * Create the frame.
     */
    public RainbowTestGUI() {
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        setBounds(100, 100, 450, 300);
        contentPane = new JPanel();
        contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
        setContentPane(contentPane);
        contentPane.setLayout(new MigLayout("", "[][grow][]", "[][][grow]"));
        
        JLabel lblFileDir = new JLabel("File Dir:");
        contentPane.add(lblFileDir, "cell 0 0,alignx trailing");
        
        txtFldFileDir = new JTextField();
        txtFldFileDir.setEditable(false);
        contentPane.add(txtFldFileDir, "cell 1 0,growx");
        txtFldFileDir.setColumns(10);
        
        btnFileDirSelect = new JButton(">");
        btnFileDirSelect.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent arg0) {
                String curr = txtFldFileDir.getText();
                if (curr.equals("")) {
                    curr = "./";
                }
                JFileChooser jfc = new JFileChooser(curr);
                jfc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
                jfc.setDialogTitle("Select FCS File Directory");
                jfc.setMultiSelectionEnabled(false);
                int resp = jfc.showOpenDialog(RainbowTestGUI.this);
                if (resp == JFileChooser.APPROVE_OPTION) {
                    String newPath = ext.verifyDirFormat(jfc.getSelectedFile().getAbsolutePath());
                    txtFldFileDir.setText(newPath);
                    String[] files = new File(newPath).list(new FilenameFilter() {
                        @Override
                        public boolean accept(File dir, String name) {
                            return name.endsWith(".fcs");
                        }
                    });
                    setFiles(files);
                }
            }
        });
        contentPane.add(btnFileDirSelect, "cell 2 0");
        
        JLabel lblGatingFile = new JLabel("Gating File:");
        contentPane.add(lblGatingFile, "cell 0 1,alignx trailing");
        
        txtFldGatingFile = new JTextField();
        txtFldGatingFile.setEditable(false);
        contentPane.add(txtFldGatingFile, "cell 1 1,growx");
        txtFldGatingFile.setColumns(10);
        
        btnGatingFileSelect = new JButton(">");
        btnGatingFileSelect.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
            }
        });
        contentPane.add(btnGatingFileSelect, "cell 2 1");
        
        scrollPane = new JScrollPane();
        scrollPane.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_ALWAYS);
        scrollPane.setVerticalScrollBarPolicy(ScrollPaneConstants.VERTICAL_SCROLLBAR_ALWAYS);
        contentPane.add(scrollPane, "cell 0 2 3 1,grow");
        
        table = new JTable();
        scrollPane.setViewportView(table);
        table.setRowSelectionAllowed(false);
    }
    
    private void setFiles(String[] files) {
        TreeSet<String> fileSet = new TreeSet<String>();
        for (String f : files) {
            fileSet.add(f);
        }
        
        Vector<String> columnNames = new Vector<String>();
        for (String f : fileSet) {
            columnNames.add(ext.rootOf(f));
        }
        
        DefaultTableModel dtm = new DefaultTableModel(columnNames, 0);
        
        for (String f : fileSet) {
            loadMeanSD(f); // TODO this is inefficient, or just wrong - if we want to look at the data, we aren't saving it.  Rainbow beads are small - can we assume that we can load all of these files at the same time?  Or are there too many files potentially?
        }
        
    }
    
    private void loadMeanSD(String file) {
        FCSDataLoader fdl = new FCSDataLoader();
        try {
            fdl.loadData(file);
        } catch (IOException e) {
            e.printStackTrace();
            // TODO handle this
            return;
        }
        
        ArrayList<String> laserNames = fdl.getAllDisplayableNames(DATA_SET.ALL); // not actually all lasers
        HashMap<String, Float> meanMap = new HashMap<String, Float>();
        HashMap<String, Float> sdMap = new HashMap<String, Float>();
        
        for (String s : laserNames) {
            float[] data = fdl.getData(s, true);
            meanMap.put(s, Array.mean(data));
            sdMap.put(s, Array.stdev(data, false));
        }
        
        fileMeans.put(file, meanMap);
        fileSDs.put(file, sdMap);
    }
    
    
}
