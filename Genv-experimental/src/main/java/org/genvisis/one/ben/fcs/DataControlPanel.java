package org.genvisis.one.ben.fcs;

import java.awt.Font;
import java.awt.Insets;

import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.SwingUtilities;

import net.miginfocom.swing.MigLayout;

import org.genvisis.common.Grafik;

import java.awt.Color;
import java.awt.SystemColor;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

public class DataControlPanel extends JPanel {

    /**
     * Create the panel.
     */
    public DataControlPanel() {
        this("", "", "", null);
    }
    
    String file;
    public static final String ACTION_DELETE = "DELETE";
    public static final String ACTION_MOVE_UP = "MOVEUP";
    public static final String ACTION_MOVE_DOWN = "MOVEDOWN";
    public static final String ACTION_INFO = "INFO";
    public static final String ACTION_LOAD = "LOAD";
    public static final String ACTION_USE = "USE";
    
    private JButton btnLoad;
    
    public DataControlPanel(String file, String sz, String dt, final ActionListener al) {
        setBorder(null);
        setLayout(new MigLayout("", "[][grow][][][][][]", "0px[][][]0px"));

        this.file = file;
        
        Insets btnIns = new Insets(0, 0, 0, 0);
        
        JButton btnDel = new JButton(Grafik.getImageIcon("images/delete2_35.png"));
        btnDel.setBackground(Color.WHITE);
        btnDel.setFont(new Font("Tahoma", Font.BOLD, 11));
        btnDel.setMargin(btnIns);
        btnDel.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                ActionEvent newEv = new ActionEvent(DataControlPanel.this, 0, DataControlPanel.this.file + "::" + ACTION_DELETE);
                al.actionPerformed(newEv);
            }
        });
        add(btnDel, "cell 0 0");
        
        JButton btnMvUp = new JButton(Grafik.getImageIcon("images/up_short.gif"));
        btnMvUp.setBackground(Color.WHITE);
        btnMvUp.setFont(new Font("Tahoma", Font.BOLD, 11));
        btnMvUp.setMargin(btnIns);
        btnMvUp.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                ActionEvent newEv = new ActionEvent(DataControlPanel.this, 0, DataControlPanel.this.file + "::" + ACTION_MOVE_UP);
                al.actionPerformed(newEv);
            }
        });
        add(btnMvUp, "cell 2 0");
        
        JButton btnMvDn = new JButton(Grafik.getImageIcon("images/down_short.gif"));
        btnMvDn.setBackground(Color.WHITE);
        btnMvDn.setFont(new Font("Tahoma", Font.BOLD, 11));
        btnMvDn.setMargin(btnIns);
        btnMvDn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                ActionEvent newEv = new ActionEvent(DataControlPanel.this, 0, DataControlPanel.this.file + "::" + ACTION_MOVE_DOWN);
                al.actionPerformed(newEv);
            }
        });
        add(btnMvDn, "cell 3 0");

        JButton button = new JButton(Grafik.getImageIcon("images/question-mark_sm.png"));
        button.setBackground(Color.WHITE);
        button.setFont(new Font("Tahoma", Font.BOLD, 11));
        button.setMargin(btnIns);
        button.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                ActionEvent newEv = new ActionEvent(DataControlPanel.this, 0, DataControlPanel.this.file + "::" + ACTION_INFO);
                al.actionPerformed(newEv);
            }
        });
        add(button, "cell 4 0");
        
        btnLoad = new JButton(Grafik.getImageIcon("images/tick-empty_sm.png"));
        btnLoad.setBackground(Color.WHITE);
        btnLoad.setFont(new Font("Tahoma", Font.BOLD, 11));
        btnLoad.setMargin(new Insets(0, 0, 0, 0));
        btnLoad.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                ActionEvent newEv = new ActionEvent(DataControlPanel.this, 0, DataControlPanel.this.file + "::" + ACTION_LOAD);
                al.actionPerformed(newEv);
            }
        });
        add(btnLoad, "cell 5 0");

        JButton btnUse = new JButton(Grafik.getImageIcon("images/right_short.gif"));
        btnUse.setBackground(Color.WHITE);
        btnUse.setFont(new Font("Tahoma", Font.BOLD, 11));
        btnUse.setMargin(new Insets(0, 0, 0, 0));
        btnUse.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                ActionEvent newEv = new ActionEvent(DataControlPanel.this, 0, DataControlPanel.this.file + "::" + ACTION_USE);
                al.actionPerformed(newEv);
            }
        });
        add(btnUse, "cell 6 0");
        
        JLabel lblFileName = new JLabel("<html><p>" + file + "</p></html>");
        lblFileName.setFont(new Font("Arial", Font.PLAIN, 9));
        add(lblFileName, "cell 0 1 7 1, growy");
        
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
    
    private boolean selected = false;
    public void setSelected(final boolean selected) {
        this.selected = selected;
        SwingUtilities.invokeLater(new Runnable() {
            @Override
            public void run() {
                setBackground(selected ? SystemColor.scrollbar : SystemColor.control);
                repaint();
            }
        });
    }
    
    public boolean isSelected() {
        return selected;
    }
    
    public void setLoaded(final boolean loaded) {
        SwingUtilities.invokeLater(new Runnable() {
            @Override
            public void run() {
                btnLoad.setIcon(Grafik.getImageIcon(loaded ? "images/tick_sm.png" : "images/tick-empty_sm.png"));
                repaint();
            }
        });
    }
    
}
