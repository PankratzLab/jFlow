package cnv.gui;

import java.awt.BorderLayout;
import java.awt.EventQueue;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.HashSet;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.border.EmptyBorder;
import javax.swing.JTextArea;
import javax.swing.JScrollPane;
import javax.swing.JLabel;

import net.miginfocom.swing.MigLayout;

import javax.swing.JTextField;
import javax.swing.JSeparator;
import javax.swing.JButton;

public class NewMarkerListDialog extends JFrame implements ActionListener {

    private JPanel contentPane;
    private JTextField textField;
    private JTextArea textArea;
    private HashSet<String> markerSet;
    private JButton btnCreate;
    private JButton btnCancel;

    /**
     * Launch the application.
     */
    public static void main(String[] args) {
        EventQueue.invokeLater(new Runnable() {
            public void run() {
                try {
                    NewMarkerListDialog frame = new NewMarkerListDialog(new String[0]);
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
    public NewMarkerListDialog(String[] markers) {
        setTitle("Create New Marker List");
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
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
        
        JLabel lblFileName = new JLabel("File name:");
        contentPane.add(lblFileName, "cell 0 2");
        
        textField = new JTextField();
        contentPane.add(textField, "cell 0 3,growx");
        textField.setColumns(10);
        
        JSeparator separator = new JSeparator();
        contentPane.add(separator, "cell 0 5,growx");
        
        JPanel panel = new JPanel();
        contentPane.add(panel, "south,alignx right,growy");
        panel.setLayout(new MigLayout("alignx right, ins 5 0 3  0", "[65px][65px]", "[23px]"));
        
        btnCreate = new JButton();
        btnCreate.addActionListener(this);
        btnCreate.setText("Create");
        panel.add(btnCreate, "cell 0 0,alignx left,aligny top");
        
        btnCancel = new JButton();
        btnCancel.addActionListener(this);
        btnCancel.setText("Cancel");
        panel.add(btnCancel, "cell 1 0,alignx left,aligny top");
        
        markerSet = new HashSet<String>();
        for (String marker : markers) {
            markerSet.add(marker);
        }
    }
    
    public String getFileName() {
        return textField.getText();
    }
    
    public String[] getMarkers() {
        String[] mkrs = textArea.getText().split("\n"); 
        for (int i = 0; i < mkrs.length; i++) {
            mkrs[i] = mkrs[i].trim();
        }
        return mkrs;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        if (e.getSource() == btnCreate) {
            // TODO doCreate
            String[] mkrs = getMarkers();
//            ArrayList<String> 
        } else {
            // TODO doCancel
        }
    }
    
    
    
}
