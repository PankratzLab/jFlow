package cnv.gui;

import java.awt.EventQueue;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.border.EmptyBorder;
import javax.swing.JDialog;
import javax.swing.JOptionPane;
import javax.swing.JTextArea;
import javax.swing.JScrollPane;
import javax.swing.JLabel;

import net.miginfocom.swing.MigLayout;

import javax.swing.JTextField;
import javax.swing.JSeparator;
import javax.swing.JButton;

import common.Files;

public class NewMarkerListDialog extends JDialog implements ActionListener {

	private static final long serialVersionUID = 1L;
	private JPanel contentPane;
    private JTextField textField;
    private JTextArea textArea;
    private HashSet<String> markerSet;
    private JButton btnCreate;
    private JButton btnCancel;
    private int returnCode = JOptionPane.DEFAULT_OPTION;
    

    /**
     * Launch the application.
     */
    public static void main(String[] args) {
        EventQueue.invokeLater(new Runnable() {
            public void run() {
                try {
                    NewMarkerListDialog frame = new NewMarkerListDialog(new String[0]);
//                    NewMarkerListDialog frame = new NewMarkerListDialog(new Project());
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
    public NewMarkerListDialog(String[] markers) {//Project proj) {
//        this.proj = proj;
        setTitle("Create New Marker List");
        setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE);
        setBounds(100, 100, 450, 300);
        contentPane = new JPanel();
        contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
        setContentPane(contentPane);
        contentPane.setLayout(new MigLayout("", "[424px,grow]", "[14px][238px][][][][][grow]"));
        
        JScrollPane scrollPane = new JScrollPane();
        contentPane.add(scrollPane, "cell 0 1,grow");
        
        textArea = new JTextArea();
        scrollPane.setViewportView(textArea);
        
        JLabel lblMarkerNamesone = new JLabel("Marker names (one per line): ");
        contentPane.add(lblMarkerNamesone, "cell 0 0,growx,aligny top");
        
        JLabel lblFileName = new JLabel("File name (full path):");
        contentPane.add(lblFileName, "cell 0 2");
        
//        String file = proj.DISPLAY_MARKERS_FILENAMES.getDefaultValue()[0];
        textField = new JTextField();//file);
        contentPane.add(textField, "cell 0 3,growx");
        textField.setColumns(10);
        
        JSeparator separator = new JSeparator();
        contentPane.add(separator, "cell 0 5,growx");
        
        JPanel panel = new JPanel();
        contentPane.add(panel, "south,alignx right,growy");
        panel.setLayout(new MigLayout("alignx right, ins 5 0 3 0", "[65px][65px]", "[23px]"));
        
        btnCreate = new JButton();
        btnCreate.addActionListener(this);
        btnCreate.setText("Create");
        panel.add(btnCreate, "cell 0 0,alignx left,aligny top");
        
        btnCancel = new JButton();
        btnCancel.addActionListener(this);
        btnCancel.setText("Cancel");
        panel.add(btnCancel, "cell 1 0,alignx left,aligny top");
        
        markerSet = new HashSet<String>();
//        String[] markers = proj.getMarkerNames();
        for (String marker : markers) {
            markerSet.add(marker);
        }
    }
    
    public String getFileName() {
        return textField.getText();
    }

    public int getReturnCode() {
        return returnCode;
    }
    
    public String[] getMarkers() {
        String[] mkrs = textArea.getText().split("\n"); 
        for (int i = 0; i < mkrs.length; i++) {
            mkrs[i] = mkrs[i].trim();
        }
        return mkrs;
    }
    
    private boolean doCreate() {
        String[] mkrs = getMarkers();
        ArrayList<String> invalid = new ArrayList<String>();
        for (String mkr : mkrs) {
        	if (!markerSet.contains(mkr)) {
        		invalid.add(mkr);
        	}
        }
        if (invalid.size() > 0) {
            String[] options = {"Ignore and Continue", "Return"};
            StringBuilder msg = new StringBuilder("Warning - ").append(invalid.size()).append(" markers are not present in the current marker set:");
            for (String inv : invalid) {
                msg.append("\n").append(inv);
            }
            int opt = JOptionPane.showOptionDialog(this, msg.toString(), "Warning - invalid markers!", JOptionPane.YES_NO_OPTION, JOptionPane.WARNING_MESSAGE, null, options, options[1]);
            if (opt != 0) {
                return false;
            }
        }
        
        String filepath = textField.getText();
        String dir = "";//proj.DATA_DIRECTORY.getValue();
        File newFile = new File(dir + filepath);
        boolean reachable = false;
        try {
            filepath = newFile.getCanonicalPath();
            reachable = true;
        } catch (IOException e) {
            reachable = false;
        }
        if (!reachable) {
            JOptionPane.showMessageDialog(null, "Error - file name [" + filepath + "] is invalid.", "Error", JOptionPane.ERROR_MESSAGE);
            return false;
        } else if (Files.exists(dir + filepath)) {
            JOptionPane.showMessageDialog(null, "Error - file [" + filepath + "] already exists.", "Error", JOptionPane.ERROR_MESSAGE);
            return false;
        }
        
        Files.writeList(mkrs, filepath);
        textField.setText(filepath);
        return true;
    }
    
    @Override
    public void actionPerformed(ActionEvent e) {
        if (e.getSource() == btnCreate) {
            boolean success = doCreate();
            if (!success) return;
            else {
                returnCode = JOptionPane.YES_OPTION;
                setVisible(false);
            }
        } else {
            returnCode = JOptionPane.NO_OPTION;
            setVisible(false);
        }
    }
    
    
    
}
