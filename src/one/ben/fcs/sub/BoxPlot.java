package one.ben.fcs.sub;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Vector;

import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.SwingUtilities;

import common.Array;
import common.Files;
import common.HashVec;
import common.Matrix;
import common.ext;
import net.miginfocom.swing.MigLayout;
import one.ben.fcs.FCSPlot;


public class BoxPlot extends JFrame {
    
//	String testFile = "C:\\Users\\Ben\\Desktop\\hb hrs P1 sample 12-May-2016.wsp FlowJo table.csv";
    String testFile = "F:\\Flow\\counts data\\hb hrs P1 sample 12-May-2016.wsp FlowJo table.csv";
    private JPanel scrollContent;
    
    public BoxPlot() {
        super();
        setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
        setBounds(FCSPlot.START_X, FCSPlot.START_Y, FCSPlot.START_WIDTH, FCSPlot.START_HEIGHT);
        
        JPanel contentPane = new JPanel(new BorderLayout());
        setContentPane(contentPane);
        
        scrollContent = new JPanel(new MigLayout("", "", ""));
        scrollContent.setBackground(Color.WHITE);
        JScrollPane scrollPane = new JScrollPane(scrollContent, JScrollPane.VERTICAL_SCROLLBAR_NEVER, JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
        
        contentPane.add(scrollPane, BorderLayout.CENTER);
    }
    
    private void loadFile(String file) {
    	setTitle(ext.removeDirectoryInfo(file));
        String[][] data = loadFileToStringMatrix(file);
        final ArrayList<BoxPanel> panels = new ArrayList<BoxPanel>();
        for (int i = 1; i < data[0].length; i++) {
            if (data[0][i].equals("")) continue;
            ArrayList<Double> panelData = new ArrayList<Double>();
            for (int r = 1; r < data.length; r++) {
                if (data[r][0].equals("Mean") || data[r][0].equals("SD") || data[r][i].equals("")) 
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
            bp.setData(lbl, Array.toDoubleArray(panelData));
            bp.setPreferredSize(new Dimension(256, 340));
            panels.add(bp);
        }
        SwingUtilities.invokeLater(new Runnable() {
            @Override
            public void run() {
                scrollContent.removeAll();
                for (int i = 0; i < panels.size(); i++) {
                    scrollContent.add(panels.get(i), "cell " + i + " 0");
                    String pts = panels.get(i).dataLabel.split("\\|")[0].trim().replaceAll("/", "  /<br />");
                    JLabel pnlLbl = new JLabel("<html><p>" + pts + "</p></html>");
                    pnlLbl.setBackground(Color.WHITE);
                    scrollContent.add(pnlLbl, "cell " + i + " 1, alignx center, aligny top");
                }
                revalidate();
                repaint();
            }
        });
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
        bp.loadFile(bp.testFile);
    }
    
}
