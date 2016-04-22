package one.ben.fcs;

import java.awt.EventQueue;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.Date;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.Vector;

import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.ScrollPaneConstants;
import javax.swing.border.EmptyBorder;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.swing.table.DefaultTableModel;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import one.ben.fcs.FCSDataLoader.DATA_SET;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import net.miginfocom.swing.MigLayout;
import common.Array;
import common.ext;

public class RainbowTestGUI extends JFrame {

    private JPanel contentPane;
    private JTextField txtFldFileDir;
    private JTextField txtFldGatingFile;
    private JButton btnFileDirSelect;
    private JButton btnGatingFileSelect;
    private JTable table;
    private JScrollPane scrollPane;

    private String[] gateFileExts = {"wsp", "wspt"};
    
    private HashMap<String, FCSDataLoader> files = new HashMap<String, FCSDataLoader>();
    
    private HashMap<String, HashMap<String, Float>> fileMeans = new HashMap<String, HashMap<String,Float>>();
    private HashMap<String, HashMap<String, Float>> fileSDs = new HashMap<String, HashMap<String,Float>>();
    private JButton btnEditColumns;
    
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
        contentPane.setBorder(null);
        setContentPane(contentPane);
        contentPane.setLayout(new MigLayout("", "[][grow][]", "[][][grow][][]"));
        
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
                    for (int i = 0; i < files.length; i++) {
                        files[i] = newPath + files[i];
                    }
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
                String curr = txtFldFileDir.getText();
                if (curr.equals("")) {
                    curr = "./";
                }
                JFileChooser jfc = new JFileChooser(curr);
                jfc.setFileSelectionMode(JFileChooser.FILES_ONLY);
                jfc.setFileFilter(new FileNameExtensionFilter("Gating-ML Files", gateFileExts));
                jfc.setDialogTitle("Select Gating-ML File");
                jfc.setMultiSelectionEnabled(false);
                int resp = jfc.showOpenDialog(RainbowTestGUI.this);
                if (resp == JFileChooser.APPROVE_OPTION) {
                    String newPath = jfc.getSelectedFile().getAbsolutePath();
                    txtFldGatingFile.setText(newPath);
                    setGateFile(newPath);
                }
            }
        });
        contentPane.add(btnGatingFileSelect, "cell 2 1");
        
        scrollPane = new JScrollPane();
        scrollPane.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_ALWAYS);
        scrollPane.setVerticalScrollBarPolicy(ScrollPaneConstants.VERTICAL_SCROLLBAR_ALWAYS);
        contentPane.add(scrollPane, "cell 0 2 3 1,grow");
        
        table = new JTable();
        scrollPane.setViewportView(table);
        table.setRowSelectionAllowed(true);
        table.setColumnSelectionAllowed(true);
        table.setCellSelectionEnabled(false);
        
        table.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        
        btnEditColumns = new JButton("Edit Columns");
        contentPane.add(btnEditColumns, "cell 0 3");
        
        setGateFile("F:\\Flow\\rainbow\\URB.wspt");
    }
    
    private void setGateFile(String filePath) {
        try {
            GateFileReader.readGateFile(filePath);
        } catch (ParserConfigurationException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (SAXException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        
    }
    
    private void setFiles(String[] fileList) {
        TreeMap<Date, String> dateMap = new TreeMap<Date, String>();
        TreeSet<String> paramSet = new TreeSet<String>();
        
        for (String f : fileList) {
            this.files.put(f, loadFCSFile(f));
            dateMap.put(files.get(f).lastModified, f);
            paramSet.addAll(files.get(f).getAllDisplayableNames(DATA_SET.COMPENSATED));
        }
        String[] paramNames = paramSet.toArray(new String[paramSet.size()]);
        String[] colNames = Array.addStrToArray("File", paramNames, 0);
        DefaultTableModel dtm = new DefaultTableModel(colNames, 0) {
            @Override
            public boolean isCellEditable(int row, int column) {
                return false;
            }
        };
        
        for (java.util.Map.Entry<Date, String> ent : dateMap.entrySet()) {
            Object[] rowData = new Object[paramSet.size() + 1];
            rowData[0] = ext.rootOf(ent.getValue());
            for (int i = 0; i < paramNames.length; i++) {
                rowData[i + 1] = Array.mean(files.get(ent.getValue()).getData(paramNames[i], true));
            }
            dtm.addRow(rowData);
        }
        
        table.setModel(dtm);
    }
    
    static class GateFileReader {
        
        public static void readGateFile(String filename) throws ParserConfigurationException, SAXException, IOException {
            DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
            DocumentBuilder builder = factory.newDocumentBuilder();
            Document doc = builder.parse(new File(filename));
            doc.getDocumentElement().normalize();
            NodeList actualPops = doc.getElementsByTagName("Population");
            for (int p = 0; p < actualPops.getLength(); p++) {
                Node pop = actualPops.item(p);
                Element popElem = (Element) pop;
                System.out.println(popElem.getAttribute("name"));
            }
//            NodeList groups = doc.getElementsByTagName("Groups");
//            for (int i = 0; i < groups.getLength(); i++) {
//                Node group = groups.item(i);
//                short type = group.getNodeType();
//                if (type == Node.ELEMENT_NODE) {
//                    Element eElement = (Element) group;
//                    NodeList groupNodes = eElement.getElementsByTagName("GroupNode");
//                    
//                    for (int j = 0; j < groupNodes.getLength(); j++) {
//                        Node groupNode = groupNodes.item(j);
//                        short nodeType = groupNode.getNodeType();
//                        if (nodeType == Node.ELEMENT_NODE) {
//                            Element groupElement = (Element) groupNode;
//                            NodeList popList = groupElement.getElementsByTagName("Subpopulations");
//                            if (popList.getLength() == 0) continue;
//                            Element popListElem = (Element) popList.item(0);
//                            NodeList actualPops = popListElem.getElementsByTagName("Population");
//                            for (int p = 0; p < actualPops.getLength(); p++) {
//                                Node pop = actualPops.item(p);
//                                Element popElem = (Element) pop;
//                                System.out.println(popElem.getAttribute("name"));
//                            }
//                        }
//                    }
//                }
//            }
            
            
        }
        
    }
    
    private FCSDataLoader loadFCSFile(String file) {
        FCSDataLoader fdl = new FCSDataLoader();
        try {
            fdl.loadData(file);
        } catch (IOException e) {
            e.printStackTrace();
            // TODO handle this
            return null;
        }
        return fdl;
    }
    
    
}
