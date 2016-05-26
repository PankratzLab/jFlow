package one.ben.fcs.sub;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Vector;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.SwingUtilities;

import common.Array;
import common.Files;
import common.HashVec;
import common.Matrix;
import common.ext;
import one.ben.fcs.FCSPlot;


public class CBCBoxPlot extends JFrame {
    
    String testFile = "F:\\Flow\\counts data\\hb hrs P1 sample 12-May-2016.wsp FlowJo table.csv";
    private JPanel scrollContent;
    
    public CBCBoxPlot() {
        super();
        setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
        setBounds(FCSPlot.START_X, FCSPlot.START_Y, FCSPlot.START_WIDTH, FCSPlot.START_HEIGHT);
        
        JPanel contentPane = new JPanel(new BorderLayout());
        setContentPane(contentPane);
        
        scrollContent = new JPanel(new FlowLayout());
        JScrollPane scrollPane = new JScrollPane(scrollContent, JScrollPane.VERTICAL_SCROLLBAR_NEVER, JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
        
        contentPane.add(scrollPane, BorderLayout.CENTER);
    }
    
    private void loadFile(String file) {
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
            BoxPanel bp = new BoxPanel();
            bp.setData(data[0][i], Array.toDoubleArray(panelData));
            bp.setPreferredSize(new Dimension(256, 340));
            panels.add(bp);
        }
        SwingUtilities.invokeLater(new Runnable() {
            @Override
            public void run() {
                scrollContent.removeAll();
                for (BoxPanel bp : panels) {
                    scrollContent.add(bp);
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
        CBCBoxPlot bp = new CBCBoxPlot();
        bp.setVisible(true);
        bp.loadFile(bp.testFile);
    }
    
}
