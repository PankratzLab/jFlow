package cnv.gui;

import java.awt.BorderLayout;
import java.awt.EventQueue;
import java.awt.Font;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.border.EmptyBorder;
import javax.swing.event.ChangeListener;
import javax.swing.AbstractAction;
import javax.swing.JScrollPane;
import javax.swing.JCheckBox;
import javax.swing.JSpinner;
import javax.swing.JSeparator;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingConstants;

import net.miginfocom.swing.MigLayout;

import javax.swing.JLabel;
import javax.swing.JButton;

public class BlastFrame extends JFrame {

    private JPanel contentPane;
    private JCheckBox chckbxExpandBlastResults;
    private JButton btnClose;
    private JSpinner spnFontSize;
    private JScrollPane scrollPane;
    private JPanel blastPanel;
    private int rowCnt = 0;
    private JCheckBox chckbxSortByLocation;
    private JSeparator separator_1;
    
    public void addBlastLabel(BlastLabel lbl) {
        JLabel locLbl = new JLabel();
        locLbl.setText(lbl.fullSegment.getUCSClocation());
        Font lblFont = Font.decode(Font.MONOSPACED).deriveFont(Font.PLAIN, 12);
        locLbl.setFont(lblFont);
        this.blastPanel.add(locLbl, "cell 0 " + rowCnt);
        JLabel strandLbl = new JLabel();
        strandLbl.setText("]" + lbl.getStrand().getEncoding() + "[");
        strandLbl.setFont(lblFont);
        this.blastPanel.add(strandLbl, "cell 1 " + rowCnt);
        this.blastPanel.add(lbl, "grow, cell 2 " + rowCnt);
        rowCnt++;
    }
    public void clearLabels() {
        this.blastPanel.removeAll();
    }
    public JScrollPane getScrollPane() {
        return scrollPane;
    }
    
    public JSpinner getSpinner() {
        return spnFontSize;
    }
    public void setSpinnerAction(final ChangeListener chl) {
        this.spnFontSize.addChangeListener(chl);
    }
    public void setExpansionCheckBoxAction(final AbstractAction act) {
        String text = this.chckbxExpandBlastResults.getText();
        this.chckbxExpandBlastResults.setAction(act);
        this.chckbxExpandBlastResults.setText(text);
    }
    public void setSortCheckBoxAction(final AbstractAction act) {
        String text = this.chckbxSortByLocation.getText();
        this.chckbxSortByLocation.setAction(act);
        this.chckbxSortByLocation.setText(text);
    }
    
    
    /**
     * Create the frame.
     */
    public BlastFrame(boolean checked) {
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        setBounds(100, 100, 450, 300);
        contentPane = new JPanel();
        contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
        contentPane.setLayout(new BorderLayout(0, 0));
        setContentPane(contentPane);

        blastPanel = new JPanel(new MigLayout("", "[200px][20px][grow]", ""));
        scrollPane = new JScrollPane(blastPanel);
        
        contentPane.add(scrollPane, BorderLayout.CENTER);
        
        JPanel panel = new JPanel();
        panel.setBorder(null);
        contentPane.add(panel, BorderLayout.SOUTH);
        panel.setLayout(new MigLayout("", "[133px][][][2px][29px][grow][]", "[23px]"));
        
        chckbxExpandBlastResults = new JCheckBox("Expand BLAST Results", checked);
        chckbxExpandBlastResults.setVerticalAlignment(SwingConstants.BOTTOM);
        panel.add(chckbxExpandBlastResults, "cell 0 0,alignx left,aligny top");
        
        separator_1 = new JSeparator();
        separator_1.setOrientation(SwingConstants.VERTICAL);
        panel.add(separator_1, "cell 1 0,growy");
        
        chckbxSortByLocation = new JCheckBox("Sort By Location");
        chckbxSortByLocation.setVerticalAlignment(SwingConstants.BOTTOM);
        panel.add(chckbxSortByLocation, "cell 2 0");
        
        JSeparator separator = new JSeparator();
        separator.setOrientation(SwingConstants.VERTICAL);
        panel.add(separator, "cell 3 0,alignx left,growy");
        
        spnFontSize = new JSpinner(new SpinnerNumberModel(17, 5, 25, 1));
        panel.add(spnFontSize, "flowx,cell 4 0,alignx left,aligny center");
        
        JLabel lblFontSize = new JLabel("Font Size");
        panel.add(lblFontSize, "cell 4 0");
        
        btnClose = new JButton();
        btnClose.setAction(new AbstractAction() {
            @Override
            public void actionPerformed(ActionEvent e) {
                BlastFrame.this.setVisible(false);
                BlastFrame.this.dispose();
            }
        });
        btnClose.setText("Close");
        panel.add(btnClose, "cell 6 0");
    }
    /**
     * Launch the application.
     */
    public static void main(String[] args) {
        EventQueue.invokeLater(new Runnable() {
            public void run() {
                try {
                    BlastFrame frame = new BlastFrame(true);
                    frame.setVisible(true);
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        });
    }
    
}
