package cnv.gui;

import java.awt.BorderLayout;
import java.awt.EventQueue;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.WindowEvent;
import java.awt.event.WindowFocusListener;
import java.util.HashMap;

import javax.swing.AbstractAction;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingConstants;
import javax.swing.border.EmptyBorder;
import javax.swing.event.ChangeListener;

import net.miginfocom.swing.MigLayout;

public class BlastFrame extends JFrame implements WindowFocusListener {
    
    private static final long serialVersionUID = 1L;
    private JPanel contentPane;
    private JCheckBox chckbxExpandBlastResults;
    private JButton btnClose;
    private JSpinner spnFontSize;
    private JScrollPane scrollPane;
    private JPanel blastPanel;
    private int rowCnt = 0;
    private JCheckBox chckbxSortByLocation;
    private JSeparator separator_1;
    private JCheckBox chckbxPinToFront;
    HashMap<BlastLabel, JLabel> alignCntLblMap = new HashMap<BlastLabel, JLabel>();
    private JCheckBox chckbxDisplayAlignmentMatch;
    
    public void addBlastLabel(BlastLabel lbl) {
        JLabel locLbl = new JLabel();
        locLbl.setText(lbl.fullSegment.getUCSClocation());
        Font lblFont = Font.decode(Font.MONOSPACED).deriveFont(Font.PLAIN, 12);
        locLbl.setFont(lblFont);
        this.blastPanel.add(locLbl, "cell 0 " + rowCnt);
        JLabel strandLbl = new JLabel();
        strandLbl.setText("]" + lbl.getStrand().getEncoding() + "[");
        strandLbl.setFont(lblFont);
        alignCntLblMap.put(lbl, strandLbl);
        this.blastPanel.add(strandLbl, "cell 1 " + rowCnt);
        this.blastPanel.add(lbl, "grow, cell 2 " + rowCnt);
        rowCnt++;
    }
    void refreshStrandLabels() {
        boolean showMatch = chckbxDisplayAlignmentMatch.isSelected();
        for (java.util.Map.Entry<BlastLabel, JLabel> strandEntry : alignCntLblMap.entrySet()) {
            strandEntry.getValue().setText("]" + strandEntry.getKey().getStrand().getEncoding() + (showMatch ? " | " + strandEntry.getKey().getAlignment() : "") + "[");
        }
    }
    public void clearLabels() {
        this.blastPanel.removeAll();
    }
    public JScrollPane getScrollPane() {
        return scrollPane;
    }
    public boolean shouldDisplayAlignment() {
        return chckbxDisplayAlignmentMatch.isSelected();
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
    public void setDisplayAlignmentCheckBoxAction(AbstractAction act) {
        String text = this.chckbxDisplayAlignmentMatch.getText();
        this.chckbxDisplayAlignmentMatch.setAction(act);
        this.chckbxDisplayAlignmentMatch.setText(text);
    }
    
    public void windowGainedFocus(WindowEvent e) {}
    public void windowLostFocus(WindowEvent e) {
        if(e.getNewState() != WindowEvent.WINDOW_CLOSED) {
            if (chckbxPinToFront != null && chckbxPinToFront.isSelected()){
                setAlwaysOnTop(true);
            } else {
                setAlwaysOnTop(false);
            }
        }

    }

    
    /**
     * Create the frame.
     */
    public BlastFrame(boolean checked) {
        setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
        setBounds(100, 100, 900, 800);
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
        panel.setLayout(new MigLayout("", "[133px][][][][2px][29px][grow][][]", "[23px]"));
        
        chckbxExpandBlastResults = new JCheckBox("Expand BLAST Results", checked);
        chckbxExpandBlastResults.setVerticalAlignment(SwingConstants.BOTTOM);
        panel.add(chckbxExpandBlastResults, "cell 0 0,alignx left,aligny top");
        
        separator_1 = new JSeparator();
        separator_1.setOrientation(SwingConstants.VERTICAL);
        panel.add(separator_1, "cell 1 0,growy");
        
        chckbxSortByLocation = new JCheckBox("Sort By Location");
        chckbxSortByLocation.setVerticalAlignment(SwingConstants.BOTTOM);
        panel.add(chckbxSortByLocation, "cell 2 0");
        
        chckbxDisplayAlignmentMatch = new JCheckBox("Display Alignment Match Counts");
        panel.add(chckbxDisplayAlignmentMatch, "cell 3 0");
        
        JSeparator separator = new JSeparator();
        separator.setOrientation(SwingConstants.VERTICAL);
        panel.add(separator, "cell 4 0,alignx left,growy");
        
        spnFontSize = new JSpinner(new SpinnerNumberModel(17, 5, 25, 1));
        panel.add(spnFontSize, "flowx,cell 5 0,alignx left,aligny center");
        
        JLabel lblFontSize = new JLabel("Font Size");
        panel.add(lblFontSize, "cell 5 0");
        
        chckbxPinToFront = new JCheckBox("Pin to Front");
        panel.add(chckbxPinToFront, "cell 7 0");
        
        btnClose = new JButton();
        btnClose.setAction(new AbstractAction() {
            private static final long serialVersionUID = 1L;
            @Override
            public void actionPerformed(ActionEvent e) {
                BlastFrame.this.setVisible(false);
//                BlastFrame.this.dispose(); // removed for state
            }
        });
        btnClose.setText("Close");
        panel.add(btnClose, "cell 8 0");
        
        addWindowFocusListener(this);
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
