package one.ben.fcs.sub;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.KeyEvent;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Properties;
import java.util.Vector;

import javax.swing.AbstractAction;
import javax.swing.ActionMap;
import javax.swing.InputMap;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.KeyStroke;
import javax.swing.event.TreeSelectionEvent;
import javax.swing.event.TreeSelectionListener;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.DefaultTreeModel;
import javax.swing.tree.TreePath;

import net.miginfocom.swing.MigLayout;
import one.ben.fcs.FCSPlot;
import common.Array;
import common.Files;
import common.Matrix;
import common.ext;

public class BoxPlot extends JFrame {
    
    private static final String TITLE_STR = "BoxPlot - Genvisis";
    
    private static final int PANEL_WIDTH = 256;
    private static final int PANEL_HEIGHT = 340;

    private static final String PROP_FILE = "boxplot.properties";
    private static final String PROPKEY_DATAFILE = "DATA_FILE";
    private static final String PROPKEY_SELECTED = "SELECTED";
    private static final String[] PROPKEY_HOTKEYS = {
        "HOTKEY_0",
        "HOTKEY_1",
        "HOTKEY_2",
        "HOTKEY_3",
        "HOTKEY_4",
        "HOTKEY_5",
        "HOTKEY_6",
        "HOTKEY_7",
        "HOTKEY_8",
        "HOTKEY_9",
    };
    private static final int[] KEYS = {
        KeyEvent.VK_0,
        KeyEvent.VK_1,
        KeyEvent.VK_2,
        KeyEvent.VK_3,
        KeyEvent.VK_4,
        KeyEvent.VK_5,
        KeyEvent.VK_6,
        KeyEvent.VK_7,
        KeyEvent.VK_8,
        KeyEvent.VK_9,
    };
    
    //	String testFile = "C:\\Users\\Ben\\Desktop\\hb hrs P1 sample 12-May-2016.wsp FlowJo table.csv";
    String testFile = "F:\\Flow\\counts data\\hb hrs P1 sample 12-May-2016.wsp FlowJo table.csv";
    private JPanel scrollContent;
    private BoxCtrlPanel ctrlPanel;
    private HashMap<String, BoxPanel> panelMap = new HashMap<String, BoxPanel>();
    private String currentFile;
    private ArrayList<String> selected = new ArrayList<String>();
    
    private volatile boolean loadingProps = false;

    private HashMap<String, ArrayList<String>> hotkeyDefs = new HashMap<String, ArrayList<String>>();

    public BoxPlot() {
        super(TITLE_STR);
        setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
        setBounds(FCSPlot.START_X, FCSPlot.START_Y, FCSPlot.START_WIDTH, FCSPlot.START_HEIGHT);
        
        JPanel contentPane = new JPanel(new BorderLayout());
        setContentPane(contentPane);
        
        scrollContent = new JPanel(new MigLayout("", "", ""));
        scrollContent.setBackground(Color.WHITE);
        JScrollPane scrollPane = new JScrollPane(scrollContent, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
        scrollPane.getVerticalScrollBar().setUnitIncrement(14);
        scrollPane.getHorizontalScrollBar().setUnitIncrement(14);
        
        
        ctrlPanel = new BoxCtrlPanel();
        ctrlPanel.addTreeSelectionListener(new TreeSelectionListener() {
            @Override
            public void valueChanged(TreeSelectionEvent e) {
                if (e.getNewLeadSelectionPath() == null) return; // extra events
                ArrayList<String[]> paths = new ArrayList<String[]>();
                for (TreePath path : ctrlPanel.getSelectedPaths()) {
                    Object[] objPath = path.getPath();
                    String[] pathStr = new String[objPath.length];
                    for (int i = 0; i < objPath.length; i++) {
                        pathStr[i] = (String) ((DefaultMutableTreeNode) objPath[i]).getUserObject();
                    }
                    paths.add(pathStr);
                }
                setDisplay(paths);
                if (!loadingProps) {
                    new Thread(new Runnable() {
                        @Override
                        public void run() {
                            saveProps();
                        }
                    }).start();
                }
            }
        });
        JSplitPane jsp = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, ctrlPanel, scrollPane);
        jsp.setDividerLocation(200);
        contentPane.add(jsp, BorderLayout.CENTER);
        
        setJMenuBar(createMenuBar());
        
        InputMap im = ctrlPanel.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW);
        
        for (int i = 0; i < PROPKEY_HOTKEYS.length; i++) {
            im.put(KeyStroke.getKeyStroke(KEYS[i], KeyEvent.CTRL_DOWN_MASK), PROPKEY_HOTKEYS[i] + "_SET");
            im.put(KeyStroke.getKeyStroke(KEYS[i], 0), PROPKEY_HOTKEYS[i] + "_SELECT");
        }
        
        ActionMap am = ctrlPanel.getActionMap();
        for (int i = 0; i < PROPKEY_HOTKEYS.length; i++) {
            final int ind = i;
            AbstractAction aaSet = new AbstractAction() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    setHotkey(ind);
                    saveProps();
                }
            };
            AbstractAction aaSel = new AbstractAction() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    selectHotkey(ind);
                }
            };
            am.put(PROPKEY_HOTKEYS[i] + "_SET", aaSet);
            am.put(PROPKEY_HOTKEYS[i] + "_SELECT", aaSel);
        }
        
        loadProps();
    }
    
    private void selectHotkey(int hotkeyIndex) {
        ArrayList<String> keyDef = hotkeyDefs.get(PROPKEY_HOTKEYS[hotkeyIndex]);
        if (keyDef == null) return;
        ArrayList<TreePath> data = new ArrayList<TreePath>();
        DefaultTreeModel dtm = (DefaultTreeModel) ctrlPanel.tree.getModel();
        for (String s : keyDef) {
            String[] pts = s.split("\\|")[0].trim().split("/");
            TreePath tp = new TreePath(dtm.getPathToRoot(ctrlPanel.getNodeForKey(Array.toStr(pts, "\t"))));
            data.add(tp);
        }
        if (data.size() > 0) {
            ctrlPanel.tree.setSelectionPaths(data.toArray(new TreePath[data.size()]));
        }
        repaint();
    }
    
    @SuppressWarnings("unchecked")
    private void setHotkey(int hotkeyIndex) {
        hotkeyDefs.put(PROPKEY_HOTKEYS[hotkeyIndex], (ArrayList<String>) selected.clone());
    }
    
    private JMenuBar createMenuBar() {
        JMenuBar menuBar = new JMenuBar();
        JMenu menu = new JMenu("File");
        menu.setMnemonic(KeyEvent.VK_F);
        menuBar.add(menu);
        
        JMenuItem load = new JMenuItem("Open File");
        load.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_O, KeyEvent.CTRL_DOWN_MASK));
        load.setMnemonic(KeyEvent.VK_O);
        load.setAction(new AbstractAction() {
            @Override
            public void actionPerformed(ActionEvent arg0) {
                String curr = currentFile;
                if (curr.equals("")) {
                    curr = "./";
                }
                JFileChooser jfc = new JFileChooser(curr);
                jfc.setFileSelectionMode(JFileChooser.FILES_ONLY);
                jfc.setDialogTitle("Select Data File");
                jfc.setMultiSelectionEnabled(false);
                int resp = jfc.showOpenDialog(BoxPlot.this);
                if (resp == JFileChooser.APPROVE_OPTION) {
                    String newPath = jfc.getSelectedFile().getAbsolutePath();
                    loadFile(newPath);
                    saveProps();
                }
            }
        });
        load.setText("Open File");
        menu.add(load);
        
        return menuBar;
    }
    
    private final String[] EXCLUDED_ROW_HEADERS = {
            "mean",
            "average",
            "sd",
            "cv",
            "cv (%)"
    };
    
    private void loadFile(String file) {
        this.currentFile = file;
        selected.clear();
        String[][] data = loadFileToStringMatrix(file);
        final ArrayList<BoxPanel> panels = new ArrayList<BoxPanel>();
        final ArrayList<String> headers = new ArrayList<String>();
        ArrayList<String> dataSources = new ArrayList<String>();
        for (int i = 1; i < data.length; i++) {
            if (data[i][0].equals("Mean") || data[i][0].equals("SD") || data[i][0].equals("")) continue;
            dataSources.add(data[i][0]);
        }
        for (int i = 1; i < data[0].length; i++) {
            if (data[0][i].equals("")) continue;
            ArrayList<Double> panelData = new ArrayList<Double>();
            for (int r = 1; r < data.length; r++) {
                int ind = ext.indexOfStr(data[r][0], EXCLUDED_ROW_HEADERS, false, false);
                if (data[r][i].equals("") || ind != -1) 
                    continue;
                panelData.add(Double.parseDouble(data[r][i]));
            }
            String lbl = data[0][i];
            if (lbl.startsWith("\"")) {
            	lbl = lbl.substring(1);
            }
            if (lbl.endsWith("\"")) {
            	lbl = lbl.substring(0, lbl.length() - 1);
            }
            BoxPanel bp = new BoxPanel();
            bp.setData(lbl, Array.toStringArray(dataSources), Array.toDoubleArray(panelData));
            bp.setPreferredSize(new Dimension(PANEL_WIDTH, PANEL_HEIGHT));
            panels.add(bp);
            panelMap.put(Array.toStr(lbl.split("\\|")[0].trim().split("/"), "\t"), bp);
            headers.add(lbl);
        }
        scrollContent.removeAll();
        ctrlPanel.setData(headers.toArray(new String[headers.size()]));
        setTitle(TITLE_STR + " - " + ext.removeDirectoryInfo(file));
        revalidate();
        repaint();
    }

    private void loadProps() {
        Properties props = new Properties();
        InputStream is = null;
        loadingProps = true;
        try {
            File f = new File(PROP_FILE);
            is = new FileInputStream(f);
            props.load(is);
            String base = props.getProperty(PROPKEY_DATAFILE, "");
            if (!base.equals("")) {
                loadFile(base);
            }
            String comp = props.getProperty(PROPKEY_SELECTED, "");
            String[] sel = comp.split(";;");
            ArrayList<TreePath> data = new ArrayList<TreePath>();
            DefaultTreeModel dtm = (DefaultTreeModel) ctrlPanel.tree.getModel();
            for (String s : sel) {
                String[] pts = s.split("\\|")[0].trim().split("/");
                TreePath tp = new TreePath(dtm.getPathToRoot(ctrlPanel.getNodeForKey(Array.toStr(pts, "\t"))));
                data.add(tp);
            }
            if (data.size() > 0) {
                ctrlPanel.tree.setSelectionPaths(data.toArray(new TreePath[data.size()]));
            }
            for (int i = 0; i < PROPKEY_HOTKEYS.length; i++) {
                String key = props.getProperty(PROPKEY_HOTKEYS[i], "");
                if (!key.equals("")) {
                    ArrayList<String> keyDef = new ArrayList<String>();
                    for (String s : key.split(";;")) {
                        keyDef.add(s);
                    }
                    hotkeyDefs.put(PROPKEY_HOTKEYS[i], keyDef);
                }
            }
        }
        catch ( Exception e ) { is = null; }
        loadingProps = false;
    }
    
    private void saveProps() {
        try {
            Properties props = new Properties();
            props.setProperty(PROPKEY_DATAFILE, currentFile);
            String sel = Array.toStr(Array.toStringArray(selected), ";;");
            props.setProperty(PROPKEY_SELECTED, sel);
            for (int i = 0; i < PROPKEY_HOTKEYS.length; i++) {
                ArrayList<String> hotKeyDef = hotkeyDefs.get(PROPKEY_HOTKEYS[i]);
                props.setProperty(PROPKEY_HOTKEYS[i], hotKeyDef == null ? "" : Array.toStr(Array.toStringArray(hotKeyDef), ";;"));
            }
            File f = new File(PROP_FILE);
            OutputStream out = new FileOutputStream( f );
            props.store(out, "");
        } catch (Exception e ) {
            e.printStackTrace();
        }
    }
    
    public void setDisplay(final ArrayList<String[]> paths) {
        int wid = scrollContent.getWidth();
        int cols = wid == 0 ? 3 : wid / PANEL_WIDTH;
        int row = 0;
        scrollContent.removeAll();
        selected.clear();
        for (int i = 0; i < paths.size(); i++) {
            if (i / cols > (row/2)) {
                row += 2;
            }
            String key = Array.toStr(paths.get(i), "\t");
            BoxPanel bp = panelMap.get(key);
            selected.add(bp.dataLabel);
            scrollContent.add(bp, "cell " + (i % cols) + " " + row);
            String pts = bp.dataLabel.split("\\|")[0].trim().replaceAll("/", "  /<br />");
            JLabel pnlLbl = new JLabel("<html><p>" + pts + "</p></html>");
            pnlLbl.setBackground(Color.WHITE);
            scrollContent.add(pnlLbl, "cell " + (i % cols) + " " + (row + 1) + ", alignx center, aligny top");
        }
        revalidate();
        repaint();
    }
    
    public static String[][] loadFileToStringMatrix(String filename) {
        BufferedReader reader = null;
        Vector<String[]> v = new Vector<String[]>(1000);
        String line;
        String[] data;

        try {
            reader = Files.getReader(filename, false, true, false);
            if (reader == null) {
                return null;
            }
            while (reader.ready()) {
                line = reader.readLine().trim();
//                from:
//                https://stackoverflow.com/questions/28587081/regex-split-on-comma-but-exclude-commas-within-parentheses-and-quotesboth-s
//                ,                         # Match literal comma
//                (?=(([^']*'){2})*[^']*$)  # Lookahead to ensure comma is followed by even number of '
//                (?=(([^"]*"){2})*[^"]*$)  # Lookahead to ensure comma is followed by even number of "
//                (?![^()]*\\))             # Negative lookahead to ensure ) is not followed by matching
//                                          # all non [()] characters in between
                data = line.split(",(?=(([^']*'){2})*[^']*$)(?=(([^\"]*\"){2})*[^\"]*$)(?![^()]*\\))", -1);
                v.add(data);
            }
            reader.close();
        } catch (FileNotFoundException fnfe) {
            System.err.println("Error: file \""+filename+"\" not found in current directory");
        } catch (IOException ioe) {
            System.err.println("Error reading file \""+filename+"\"");
        }
        
        return Matrix.toStringArrays(v);
    }
    
    public static void main(String[] args) {
        BoxPlot bp = new BoxPlot();
        bp.setVisible(true);
//        bp.loadFile(bp.testFile);
    }
    
}
