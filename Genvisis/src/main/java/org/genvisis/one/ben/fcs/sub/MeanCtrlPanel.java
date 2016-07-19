package one.ben.fcs.sub;

import javax.swing.JPanel;

import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.util.ArrayList;

import net.miginfocom.swing.MigLayout;

import javax.swing.AbstractAction;
import javax.swing.ActionMap;
import javax.swing.ButtonGroup;
import javax.swing.InputMap;
import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.KeyStroke;
import javax.swing.JRadioButton;

public class MeanCtrlPanel extends JPanel {

    private JButton btnPrev;
    private JButton btnNext;
    private JLabel lblNext;
    private JLabel lblPrev;

    /**
     * Create the panel.
     */
    public MeanCtrlPanel() {
        setBackground(Color.WHITE);
        setLayout(new MigLayout("", "[][][grow][][][75px][][75px][]", "[]"));
        
        ButtonGroup bg = new ButtonGroup();
        rdbtnDotline = new JRadioButton();
        rdbtnDotline.addActionListener(plotChange);
        rdbtnDotline.setText("Dot/Line");
        rdbtnDotline.setBackground(Color.WHITE);
        bg.add(rdbtnDotline);
        add(rdbtnDotline, "cell 0 0");
        rdbtnDotline.setVisible(false);
        
        rdbtnBox = new JRadioButton();
        rdbtnBox.setSelected(true);
        rdbtnBox.addActionListener(plotChange);
        rdbtnBox.setText("Box");
        rdbtnBox.setBackground(Color.WHITE);
        bg.add(rdbtnBox);
        add(rdbtnBox, "cell 1 0");
        rdbtnBox.setVisible(false);
        
        btnPrev = new JButton("<<");
        add(btnPrev, "cell 3 0");
        
        lblPrev = new JLabel("prev");
        add(lblPrev, "cell 5 0");
        
        label = new JLabel("|");
        add(label, "cell 6 0");
        
        lblNext = new JLabel("next");
        add(lblNext, "cell 7 0,alignx right");
        
        btnNext = new JButton(">>");
        add(btnNext, "cell 8 0");
        
        btnPrev.addActionListener(list);
        btnNext.addActionListener(list);
        
        InputMap im = getInputMap(WHEN_IN_FOCUSED_WINDOW);
        ActionMap am = getActionMap();
        
        im.put(KeyStroke.getKeyStroke(KeyEvent.VK_PAGE_UP, 0), "prev");
        im.put(KeyStroke.getKeyStroke(KeyEvent.VK_PAGE_DOWN, 0), "next");
        im.put(KeyStroke.getKeyStroke(KeyEvent.VK_LEFT, 0), "prev");
        im.put(KeyStroke.getKeyStroke(KeyEvent.VK_RIGHT, 0), "next");
        
        am.put("prev", new AbstractAction() {
            @Override
            public void actionPerformed(ActionEvent e) {
                btnPrev.doClick();
            }
        });
        am.put("next", new AbstractAction() {
            @Override
            public void actionPerformed(ActionEvent e) {
                btnNext.doClick();
            }
        });
        
    }
    
    ActionListener plotChange = new ActionListener() {
        @Override
        public void actionPerformed(ActionEvent e) {
            if (e.getSource() == rdbtnBox) {
                plot.actionPerformed(new ActionEvent(rdbtnBox, ActionEvent.ACTION_FIRST, "BOX"));
            } else if (e.getSource() == rdbtnDotline) {
                plot.actionPerformed(new ActionEvent(rdbtnDotline, ActionEvent.ACTION_LAST, "DOT"));
            }
        }
    };
    
    ActionListener list = new ActionListener() {
        @Override
        public void actionPerformed(ActionEvent e) {
            if (e.getSource() == btnPrev) {
                if (ind == 0) return;
                int newInd = ind - 1;
                String newCol = cols.get(newInd);
                if (newInd == 0) {
                    btnPrev.setEnabled(false);
                }
                if (newInd < cols.size() - 1) {
                    btnNext.setEnabled(true);
                }
                lblPrev.setText(newInd > 0 ? getColumnLabel(newInd - 1) : "---");
                lblNext.setText(newInd < cols.size() - 1 ? getColumnLabel(newInd + 1) : "---");
                change.actionPerformed(new ActionEvent(btnPrev, ActionEvent.ACTION_FIRST, newCol));
            } else if (e.getSource() == btnNext) {
                if (ind == cols.size() - 1) return;
                int newInd = ind + 1;
                String newCol = cols.get(newInd);
                if (newInd == cols.size() - 1) {
                    btnNext.setEnabled(false);
                }
                if (ind > 0) {
                    btnPrev.setEnabled(true);
                }
                lblPrev.setText(newInd > 0 ? getColumnLabel(newInd - 1) : "---");
                lblNext.setText(newInd < cols.size() - 1 ? getColumnLabel(newInd + 1) : "---");
                change.actionPerformed(new ActionEvent(btnNext, ActionEvent.ACTION_LAST, newCol));
            }
        }
    };
    
    ArrayList<String> cols;
    int ind;
    private JLabel label;
    
    ActionListener change, plot;
    LabelPresenter labeler;
    private JRadioButton rdbtnDotline;
    private JRadioButton rdbtnBox;
    
    static abstract class LabelPresenter {
        public abstract String getPresentationView(String label);
    }
    
    private String getColumnLabel(int ind) {
        return labeler == null ? cols.get(ind) : labeler.getPresentationView(cols.get(ind));
    }
    
    public void setLabelPresenter(LabelPresenter lp) {
        this.labeler = lp;
    }
    
    public void setColumns(ArrayList<String> cols, int ind) {
        this.cols = cols;
        this.ind = ind;
        
        btnPrev.setEnabled(ind > 0);
        btnNext.setEnabled(ind < cols.size() - 1);
        lblPrev.setText(ind > 0 ? getColumnLabel(ind - 1) : "---");
        lblNext.setText(ind < cols.size() - 1 ? getColumnLabel(ind + 1) : "---");
        
        repaint();
    }

    public void setChangeListener(ActionListener prevLst) {
        this.change = prevLst;
    }

    public void setPlotChangeListener(ActionListener plotLst) {
        this.plot = plotLst;
        rdbtnDotline.setVisible(true);
        rdbtnBox.setVisible(true);
    }

    
    
}
