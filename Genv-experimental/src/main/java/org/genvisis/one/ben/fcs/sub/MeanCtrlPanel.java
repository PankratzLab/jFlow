package org.genvisis.one.ben.fcs.sub;

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
import javax.swing.JCheckBox;

import org.genvisis.cnv.gui.ColorIcon;

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
    setLayout(new MigLayout("", "[][][grow][][][75px][][75px][]", "[][]"));

    ButtonGroup bg = new ButtonGroup();
    rdbtnDotline = new JRadioButton();
    rdbtnDotline.addActionListener(plotChange);
    
    panel = new JPanel();
    panel.setBackground(Color.WHITE);
    add(panel, "cell 0 0 9 1,grow");
    panel.setLayout(new MigLayout("ins 0", "[][][][grow][][][][]", "[]"));
    
    lblSd = new JLabel("1 SD");
    panel.add(lblSd, "cell 0 0");
    
    lblSd_1 = new JLabel("2 SD");
    panel.add(lblSd_1, "cell 1 0");
    
    lblMean = new JLabel("15% Mean");
    panel.add(lblMean, "cell 2 0");
    
    chk1Sd = new JCheckBox("+/- 1 SD");
    chk1Sd.setBackground(Color.WHITE);
    panel.add(chk1Sd, "cell 4 0");
    
    chk2Sd = new JCheckBox("+/- 2 SD");
    chk2Sd.setBackground(Color.WHITE);
    panel.add(chk2Sd, "cell 5 0");
    
    chk15 = new JCheckBox("+/- 15% of mean");
    chk15.setBackground(Color.WHITE);
    panel.add(chk15, "cell 6 0");
    
    chkReg = new JCheckBox("Regression Lines");
    chkReg.setBackground(Color.WHITE);
    panel.add(chkReg, "cell 7 0");
    rdbtnDotline.setText("Dot/Line");
    rdbtnDotline.setBackground(Color.WHITE);
    bg.add(rdbtnDotline);
    add(rdbtnDotline, "cell 0 1");
    rdbtnDotline.setVisible(false);

    rdbtnBox = new JRadioButton();
    rdbtnBox.setSelected(true);
    rdbtnBox.addActionListener(plotChange);
    rdbtnBox.setText("Box");
    rdbtnBox.setBackground(Color.WHITE);
    bg.add(rdbtnBox);
    add(rdbtnBox, "cell 1 1");
    rdbtnBox.setVisible(false);

    btnPrev = new JButton("<<");
    add(btnPrev, "cell 3 1");

    lblPrev = new JLabel("prev");
    add(lblPrev, "cell 5 1");

    label = new JLabel("|");
    add(label, "cell 6 1");

    lblNext = new JLabel("next");
    add(lblNext, "cell 7 1,alignx right");

    btnNext = new JButton(">>");
    add(btnNext, "cell 8 1");

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
        if (ind == 0)
          return;
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
        if (ind == cols.size() - 1)
          return;
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
  private JPanel panel;
  private JCheckBox chk1Sd;
  private JCheckBox chk2Sd;
  private JCheckBox chk15;
  private JCheckBox chkReg;
  private OneDPanel linkedPanel;
  private JLabel lblSd;
  private JLabel lblSd_1;
  private JLabel lblMean;

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
  
  public void link(OneDPanel meanPanel) {
    linkedPanel = meanPanel;
    chk1Sd.setSelected(linkedPanel.isShow1SDLines());
    chk2Sd.setSelected(linkedPanel.isShow2SDLines());
    chk15.setSelected(linkedPanel.isShowMean15Line());
    chkReg.setSelected(linkedPanel.isShowRegressionLine());
    
    chk1Sd.setAction(new AbstractAction() {
      @Override
      public void actionPerformed(ActionEvent e) {
        linkedPanel.setShow1SDLines(chk1Sd.isSelected());
        linkedPanel.paintAgain();
      }
    });
    chk2Sd.setAction(new AbstractAction() {
      @Override
      public void actionPerformed(ActionEvent e) {
        linkedPanel.setShow2SDLines(chk2Sd.isSelected());
        linkedPanel.paintAgain();
      }
    });
    chk15.setAction(new AbstractAction() {
      @Override
      public void actionPerformed(ActionEvent e) {
        linkedPanel.setShowMean15Line(chk15.isSelected());
        linkedPanel.paintAgain();
      }
    });
    chkReg.setAction(new AbstractAction() {
      @Override
      public void actionPerformed(ActionEvent e) {
        linkedPanel.setShowRegressionLine(chkReg.isSelected());
        linkedPanel.paintAgain();
      }
    });
    
    chk1Sd.setText("+/- 1 SD");
    chk2Sd.setText("+/- 2 SD");
    chk15.setText("+/- 15% of mean");
    chkReg.setText("Regression Lines");
    

    lblSd.setIcon(new ColorIcon(12, 12, linkedPanel.getColorScheme()[linkedPanel.get1SDColor()]));
    lblSd_1.setIcon(new ColorIcon(12, 12, linkedPanel.getColorScheme()[linkedPanel.get2SDColor()]));
    lblMean.setIcon(new ColorIcon(12, 12, linkedPanel.getColorScheme()[linkedPanel.getMeanColor()]));
  }



}
