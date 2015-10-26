package cnv.gui;

import java.awt.BorderLayout;
import java.awt.EventQueue;

import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.border.EmptyBorder;

import net.miginfocom.swing.MigLayout;

import javax.swing.JLabel;

import java.awt.Font;

import javax.swing.JComboBox;
import javax.swing.JButton;

import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;

import javax.swing.JSeparator;

import common.Array;
import cnv.filesys.SourceFileHeaderData;

public class SourceFileHeaderGUI extends JDialog {

    private static final long serialVersionUID = 1L;
    private JPanel contentPane;
    private JComboBox<String> cbSnpInd;
    private JComboBox<String> cbSampInd;
    private JComboBox<String> cbA1Geno;
    private JComboBox<String> cbA2Geno;
    private JComboBox<String> cbA1AB;
    private JComboBox<String> cbA2AB;
    private JComboBox<String> cbX;
    private JComboBox<String> cbY;
    private JComboBox<String> cbBAF;
    private JComboBox<String> cbLRR;
    private JComboBox<String> cbGC;
    private JComboBox<String> cbXRaw;
    private JComboBox<String> cbYRaw;
    private JComboBox<String> cbR;
    private JComboBox<String> cbTheta;

    private boolean cancelled = false;
    
    private static String MISSING_STR = "[missing]";
    private static String FILENAME_STR = "[filename]";
    
    public static int MISSING_IND = -1;
    public static int FILENAME_IND = -99;
    
    /**
     * Launch the application.
     */
    public static void main(String[] args) {
        EventQueue.invokeLater(new Runnable() {
            public void run() {
                try {
                    SourceFileHeaderGUI frame = new SourceFileHeaderGUI(new SourceFileHeaderData());
                    frame.setVisible(true);
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        });
    }
    
    private void doClose(boolean cancel) {
        this.cancelled = cancel;
        setVisible(false);
    }

    public boolean wasCancelled() {
        return cancelled;
    }

    /**
     * Create the frame.
     * @param reportHdr 
     */
    public SourceFileHeaderGUI(SourceFileHeaderData reportHdr) {
        setTitle("Assign Final Report Data Columns");
        setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
        setBounds(100, 100, 450, 620);
        contentPane = new JPanel();
        contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
        setContentPane(contentPane);
        contentPane.setLayout(new BorderLayout(0, 0));
        
        JPanel panel = new JPanel();
        contentPane.add(panel, BorderLayout.CENTER);
        panel.setLayout(new MigLayout("", "[grow][][grow 200][grow]", "[][][][][][][][][][][][][][][][][][][][][][grow]"));
        
        JLabel lblRequired = new JLabel("Required:");
        Font fontHeaderLabels = new Font("Tahoma", Font.BOLD, 12);
        lblRequired.setFont(fontHeaderLabels);
        panel.add(lblRequired, "pad 0 -15 0 0,cell 1 1");
        
        JLabel lblSnpIndex = new JLabel("SNP Index:");
        panel.add(lblSnpIndex, "cell 1 2");
        
        String[] headerParts = reportHdr.cols;
        String[] headerOptions = Array.addStrToArray(MISSING_STR, headerParts, 0);
        String[] sampOptions = Array.addStrToArray(FILENAME_STR, headerParts, 0);
        
        cbSnpInd = new JComboBox<String>(headerParts);
        cbSnpInd.setSelectedIndex(reportHdr.colSnpIdent == -1 ? 0 : reportHdr.colSnpIdent); // TODO SNP INDEX shouldn't ever be -1
        panel.add(cbSnpInd, "cell 2 2,growx");
        
        JLabel lblSampleId = new JLabel("Sample ID:");
        panel.add(lblSampleId, "cell 1 3");
        
        cbSampInd = new JComboBox<String>(sampOptions);
        cbSampInd.setSelectedIndex(reportHdr.colSampleIdent == -1 ? 0 : 1 + reportHdr.colSampleIdent);
        panel.add(cbSampInd, "cell 2 3,growx");
        
        JLabel lblRecommended = new JLabel("Recommended:");
        lblRecommended.setFont(fontHeaderLabels);
        panel.add(lblRecommended, "pad 0 -15 0 0,cell 1 4");
        
        JLabel lblAllele_2 = new JLabel("Allele 1 - Genotype:");
        panel.add(lblAllele_2, "cell 1 5");
        
        cbA1Geno = new JComboBox<String>(headerOptions);
        cbA1Geno.setSelectedIndex(reportHdr.colGeno1 == -1 ? 0 : 1 + reportHdr.colGeno1);
        panel.add(cbA1Geno, "cell 2 5,growx");
        
        JLabel lblAllele_3 = new JLabel("Allele 2 - Genotype:");
        panel.add(lblAllele_3, "cell 1 6");
        
        cbA2Geno = new JComboBox<String>(headerOptions);
        cbA2Geno.setSelectedIndex(reportHdr.colGeno2 == -1 ? 0 : 1 + reportHdr.colGeno2);
        panel.add(cbA2Geno, "cell 2 6,growx");
        
        JLabel lblAllele = new JLabel("Allele 1 - AB:");
        panel.add(lblAllele, "cell 1 7");
        
        cbA1AB = new JComboBox<String>(headerOptions);
        cbA1AB.setSelectedIndex(reportHdr.colGenoAB1 == -1 ? 0 : 1 + reportHdr.colGenoAB1);
        panel.add(cbA1AB, "cell 2 7,growx");
        
        JLabel lblAllele_1 = new JLabel("Allele 2 - AB:");
        panel.add(lblAllele_1, "cell 1 8");
        
        cbA2AB = new JComboBox<String>(headerOptions);
        cbA2AB.setSelectedIndex(reportHdr.colGenoAB2 == -1 ? 0 : 1 + reportHdr.colGenoAB2);
        panel.add(cbA2AB, "cell 2 8,growx");
        
        JLabel lblX = new JLabel("X:");
        panel.add(lblX, "cell 1 9");
        
        cbX = new JComboBox<String>(headerOptions);
        cbX.setSelectedIndex(reportHdr.colX == -1 ? 0 : 1 + reportHdr.colX);
        panel.add(cbX, "cell 2 9,growx");
        
        JLabel lblY = new JLabel("Y:");
        panel.add(lblY, "cell 1 10");
        
        cbY = new JComboBox<String>(headerOptions);
        cbY.setSelectedIndex(reportHdr.colY == -1 ? 0 : 1 + reportHdr.colY);
        panel.add(cbY, "cell 2 10,growx");
        
        JLabel lblBAlleleFreq = new JLabel("B Allele Freq.:");
        panel.add(lblBAlleleFreq, "cell 1 11");
        
        cbBAF = new JComboBox<String>(headerOptions);
        cbBAF.setSelectedIndex(reportHdr.colBAF == -1 ? 0 : 1 + reportHdr.colBAF);
        panel.add(cbBAF, "cell 2 11,growx");
        
        JLabel lblLogrRatio = new JLabel("Log-R Ratio:");
        panel.add(lblLogrRatio, "cell 1 12");
        
        cbLRR = new JComboBox<String>(headerOptions);
        cbLRR.setSelectedIndex(reportHdr.colLRR == -1 ? 0 : 1 + reportHdr.colLRR);
        panel.add(cbLRR, "cell 2 12,growx");
        
        JLabel lblConfidence = new JLabel("Confidence:");
        panel.add(lblConfidence, "cell 1 13");
        
        cbGC = new JComboBox<String>(headerOptions);
        cbGC.setSelectedIndex(reportHdr.colGC == -1 ? 0 : 1 + reportHdr.colGC);
        panel.add(cbGC, "cell 2 13,growx");
        
        JLabel lblOptional = new JLabel("Optional:");
        lblOptional.setFont(fontHeaderLabels);
        panel.add(lblOptional, "pad 0 -15 0 0,cell 1 14");
        
        JLabel lblXRaw = new JLabel("X Raw:");
        panel.add(lblXRaw, "cell 1 15");
        
        cbXRaw = new JComboBox<String>(headerOptions);
        cbXRaw.setSelectedIndex(reportHdr.colXRaw == -1 ? 0 : 1 + reportHdr.colXRaw);
        panel.add(cbXRaw, "cell 2 15,growx");
        
        JLabel lblYRaw = new JLabel("Y Raw:");
        panel.add(lblYRaw, "cell 1 16");
        
        cbYRaw = new JComboBox<String>(headerOptions);
        cbYRaw.setSelectedIndex(reportHdr.colYRaw == -1 ? 0 : 1 + reportHdr.colYRaw);
        panel.add(cbYRaw, "cell 2 16,growx");
        
        JLabel lblR = new JLabel("R:");
        panel.add(lblR, "cell 1 17");
        
        cbR = new JComboBox<String>(headerOptions);
        cbR.setSelectedIndex(reportHdr.colR == -1 ? 0 : 1 + reportHdr.colR);
        panel.add(cbR, "cell 2 17,growx");
        
        JLabel lblTheta = new JLabel("Theta:");
        panel.add(lblTheta, "cell 1 18");
        
        cbTheta = new JComboBox<String>(headerOptions);
        cbTheta.setSelectedIndex(reportHdr.colTheta == -1 ? 0 : 1 + reportHdr.colTheta);
        panel.add(cbTheta, "cell 2 18,growx");
        
        JSeparator separator = new JSeparator();
        panel.add(separator, "south");
        
        JPanel panel_1 = new JPanel();
        FlowLayout flowLayout = (FlowLayout) panel_1.getLayout();
        flowLayout.setAlignment(FlowLayout.TRAILING);
        contentPane.add(panel_1, BorderLayout.SOUTH);
        
        JButton btnOk = new JButton("OK");
        btnOk.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                doClose(false);
            }
        });
        btnOk.setMnemonic(KeyEvent.VK_O);
        panel_1.add(btnOk);
        
        JButton btnCancel = new JButton("Cancel");
        btnCancel.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                doClose(true);
            }
        });
        btnCancel.setMnemonic(KeyEvent.VK_C);
        panel_1.add(btnCancel);
    }

    public int getSelectedSNPIndex() {
        return cbSnpInd.getSelectedIndex();
    }

    public int getSelectedSampleID() {
        return cbSampInd.getSelectedIndex() == 0 ? FILENAME_IND : cbSampInd.getSelectedIndex() - 1;
    }
    
    private int getValue(JComboBox<String> cb) {
        return MISSING_STR.equals(cb.getSelectedItem()) ? MISSING_IND : cb.getSelectedIndex() - 1;
    }
    
    public int getSelectedGeno1() { return getValue(cbA1Geno); }
    public int getSelectedGeno2() { return getValue(cbA2Geno); }
    public int getSelectedAB1() { return getValue(cbA1AB); }
    public int getSelectedAB2() { return getValue(cbA2AB); }
    public int getSelectedBAF() { return getValue(cbBAF); }
    public int getSelectedLRR() { return getValue(cbLRR); }
    public int getSelectedGC() { return getValue(cbGC); }
    public int getSelectedR() { return getValue(cbR); }
    public int getSelectedTheta() { return getValue(cbTheta); }
    public int getSelectedX() { return getValue(cbX); }
    public int getSelectedY() { return getValue(cbY); }
    public int getSelectedXRaw() { return getValue(cbXRaw); }
    public int getSelectedYRaw() { return getValue(cbYRaw); }
}
