package one.ben.fcs.sub;

import java.awt.EventQueue;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.ScrollPaneConstants;
import javax.swing.border.BevelBorder;
import javax.swing.border.EmptyBorder;

import net.miginfocom.swing.MigLayout;
import one.ben.fcs.FCSPlot;
import cnv.gui.JAccordionPanel;

public class DataExportGUI extends JDialog {

    private JPanel contentPane;

    /**
     * Launch the application.
     */
//    public static void main(String[] args) {
//        EventQueue.invokeLater(new Runnable() {
//            public void run() {
//                try {
//                    DataExportGUI frame = new DataExportGUI(null);
//                    frame.setVisible(true);
//                } catch (Exception e) {
//                    e.printStackTrace();
//                }
//            }
//        });
//    }

    FCSPlot plot;
    ArrayList<JCheckBox> boxes;
    
    public DataExportGUI(FCSPlot fcsPlot) {
        this.plot = fcsPlot;
        setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
        setBounds(100, 100, 450, 300);
        contentPane = new JPanel();
        contentPane.setBorder(null);
        contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
        setContentPane(contentPane);
        contentPane.setLayout(new MigLayout("ins 0", "[grow]", "[grow]"));
        
        JScrollPane scrollPane = new JScrollPane();
        scrollPane.setVerticalScrollBarPolicy(ScrollPaneConstants.VERTICAL_SCROLLBAR_ALWAYS);
        contentPane.add(scrollPane, "cell 0 0,grow");
        
        JPanel panel_1 = new JPanel();
        scrollPane.setViewportView(panel_1);
        panel_1.setLayout(new MigLayout("", "[grow]", "[]"));
        JAccordionPanel fileAP = new JAccordionPanel();
        JAccordionPanel gateAP = new JAccordionPanel();
        JAccordionPanel optAP = new JAccordionPanel();
        
        ArrayList<JAccordionPanel> bg = new ArrayList<JAccordionPanel>();
        
        JLabel fileLabel = new JLabel("<html><u>Select Files...</u></html>");
        fileLabel.setFont(fileLabel.getFont().deriveFont(Font.PLAIN, 14));
        fileAP.topPanel.add(fileLabel, "pad 0 10 0 0, cell 0 0, grow");
        fileAP.topPanel.setBorder(BorderFactory.createBevelBorder(BevelBorder.RAISED));
        fileAP.contentPanel.setBorder(BorderFactory.createBevelBorder(BevelBorder.LOWERED));
        fileAP.addToGroup(bg);

        JLabel gateLabel = new JLabel("<html><u>Select Gating...</u></html>");
        gateLabel.setFont(gateLabel.getFont().deriveFont(Font.PLAIN, 14));
        gateAP.topPanel.add(gateLabel, "pad 0 10 0 0, cell 0 0, grow");
        gateAP.topPanel.setBorder(BorderFactory.createBevelBorder(BevelBorder.RAISED));
        gateAP.contentPanel.setBorder(BorderFactory.createBevelBorder(BevelBorder.LOWERED));
        gateAP.addToGroup(bg);

        JLabel optLabel = new JLabel("<html><u>Export Options...</u></html>");
        optLabel.setFont(optLabel.getFont().deriveFont(Font.PLAIN, 14));
        optAP.topPanel.add(optLabel, "pad 0 10 0 0, cell 0 0, grow");
        optAP.topPanel.setBorder(BorderFactory.createBevelBorder(BevelBorder.RAISED));
        optAP.contentPanel.setBorder(BorderFactory.createBevelBorder(BevelBorder.LOWERED));
        optAP.addToGroup(bg);

        gateAP.shrink();
        optAP.shrink();
        
        panel_1.add(fileAP, "cell 0 0, growx");
        panel_1.add(gateAP, "cell 0 1, growx");
        panel_1.add(optAP, "cell 0 2, growx");
        
        JPanel pnl = fileAP.contentPanel;
        pnl.setLayout(new MigLayout("", "[grow]", ""));
        
        if (plot != null) {
            ArrayList<String> files = plot.getAddedFiles();
            boxes = new ArrayList<JCheckBox>();
            for (int i = 0; i < files.size(); i++) {
                JCheckBox box = new JCheckBox(files.get(i));
                pnl.add(box, "growx, cell 0 " + i);
                boxes.add(box);
            }
        }
        
        
        
        JPanel panel = new JPanel();
        contentPane.add(panel, "south,grow");
        panel.setLayout(new MigLayout("ins 5 0 0 0", "[grow][][]", "[]"));
        
        JButton btnExport = new JButton("Export");
        panel.add(btnExport, "cell 1 0");
        
        JButton btnCancel = new JButton("Cancel");
        btnCancel.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent arg0) {
                setVisible(false);
                dispose();
            }
        });
        panel.add(btnCancel, "cell 2 0");
    }

}
