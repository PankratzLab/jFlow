package one.ben.fcs;

import java.awt.Font;
import java.awt.Insets;

import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JPanel;

import common.Grafik;

import net.miginfocom.swing.MigLayout;
import java.awt.Color;
import java.awt.SystemColor;

public class DataControlPanel extends JPanel {

    /**
     * Create the panel.
     */
    public DataControlPanel() {
        this("", "", "");
    }
    
    public DataControlPanel(String file, String sz, String dt) {
        setBorder(null);
        setLayout(new MigLayout("", "[][grow][][][][]", "0px[][][]0px"));

        Insets btnIns = new Insets(0, 0, 0, 0);
        
        JButton btnDel = new JButton(Grafik.getImageIcon("images/delete2_35.png"));
        btnDel.setBackground(Color.WHITE);
        btnDel.setFont(new Font("Tahoma", Font.BOLD, 11));
        btnDel.setMargin(btnIns);
        add(btnDel, "cell 0 0");
        
        JButton btnMvUp = new JButton(Grafik.getImageIcon("images/up_short.gif"));
        btnMvUp.setBackground(Color.WHITE);
        btnMvUp.setFont(new Font("Tahoma", Font.BOLD, 11));
        btnMvUp.setMargin(btnIns);
        add(btnMvUp, "cell 2 0");
        
        JButton btnMvDn = new JButton(Grafik.getImageIcon("images/down_short.gif"));
        btnMvDn.setBackground(Color.WHITE);
        btnMvDn.setFont(new Font("Tahoma", Font.BOLD, 11));
        btnMvDn.setMargin(btnIns);
        add(btnMvDn, "cell 3 0");

        JButton button = new JButton(Grafik.getImageIcon("images/question-mark_sm.png"));
        button.setBackground(Color.WHITE);
        button.setFont(new Font("Tahoma", Font.BOLD, 11));
        button.setMargin(btnIns);
        add(button, "cell 4 0");
        
        JButton btnLoad = new JButton(Grafik.getImageIcon("images/tick-empty_sm.png"));
        btnLoad.setBackground(Color.WHITE);
        btnLoad.setFont(new Font("Tahoma", Font.BOLD, 11));
        btnLoad.setMargin(new Insets(0, 0, 0, 0));
        add(btnLoad, "cell 5 0");
        
        JLabel lblFileName = new JLabel("<html><p>" + file + "</p></html>");
        lblFileName.setFont(new Font("Arial", Font.PLAIN, 9));
        add(lblFileName, "cell 0 1 6 1");
        
        JLabel lblDate = new JLabel("Date:");
        lblDate.setFont(new Font("Arial", Font.BOLD, 8));
        add(lblDate, "cell 0 2");
        
        JLabel lblDateLbl = new JLabel(dt);
        lblDateLbl.setFont(new Font("Arial", Font.PLAIN, 7));
        add(lblDateLbl, "cell 1 2,alignx left");
        
        JLabel lblSize = new JLabel("Size:");
        lblSize.setFont(new Font("Arial", Font.BOLD, 8));
        add(lblSize, "cell 2 2 2 1,alignx right");
        
        JLabel lblSzlbl = new JLabel(sz);
        lblSzlbl.setFont(new Font("Arial", Font.PLAIN, 7));
        add(lblSzlbl, "cell 4 2 2 1,alignx left");
        
    }
    
    public void setSelected(boolean selected) {
        setBackground(selected ? SystemColor.scrollbar : SystemColor.control);
    }

}
