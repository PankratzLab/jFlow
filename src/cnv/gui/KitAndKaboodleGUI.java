package cnv.gui;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Insets;

import javax.swing.AbstractAction;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextField;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;
import javax.swing.border.EmptyBorder;
import javax.swing.border.LineBorder;

import net.miginfocom.swing.MigLayout;

import javax.swing.JLabel;

import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.util.ArrayList;
import java.util.HashMap;

import javax.swing.JCheckBox;

import cnv.filesys.Project;
import cnv.manage.KitAndKaboodle;
import cnv.manage.KitAndKaboodle.STEPS;
import common.Array;

public class KitAndKaboodleGUI extends JDialog {

    private static final long serialVersionUID = 1L;

    private final JPanel contentPanel = new JPanel();
    
    volatile boolean[] selected = Array.booleanArray(KitAndKaboodle.STEPS.values().length, true);

    private HashMap<STEPS, JCheckBox> checkBoxes = new HashMap<STEPS, JCheckBox>();
    private HashMap<STEPS, JLabel> descLabels = new HashMap<STEPS, JLabel>();
    private HashMap<STEPS, ArrayList<JLabel>> requirementsLabels = new HashMap<STEPS, ArrayList<JLabel>>();
    private HashMap<STEPS, JPanel> panels = new HashMap<KitAndKaboodle.STEPS, JPanel>();
    
    volatile boolean cancelled = false;

    /**
     * Create the dialog.
     */
    public KitAndKaboodleGUI(Project proj) {
        getContentPane().setLayout(new BorderLayout());
        contentPanel.setBorder(new EmptyBorder(5, 5, 5, 5));
        JPanel optionPanel = new JPanel();
        JScrollPane scrollPane = new JScrollPane(optionPanel) {
//            @Override
//            public Dimension getPreferredSize() { return contentPanel.getPreferredSize(); }
//            @Override
//            public Dimension getMaximumSize() { return contentPanel.getMaximumSize(); }
//            @Override
//            public Dimension getMinimumSize() { return contentPanel.getMinimumSize(); }
        };
        scrollPane.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
        scrollPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
        contentPanel.setLayout(new MigLayout("", "[grow]", "[][][]"));
        contentPanel.add(scrollPane, "cell 0 2, grow");
        getContentPane().add(contentPanel, BorderLayout.CENTER);
        optionPanel.setLayout(new MigLayout("", "[grow]", "[][][]"));
        {
            JLabel lblKitAndKaboodle = new JLabel("Kit and Kaboodle Steps:");
            lblKitAndKaboodle.setFont(new Font("Arial", Font.BOLD, 16));
            contentPanel.add(lblKitAndKaboodle, "cell 0 0,alignx center");
        }
        {
            for (int i = 0; i < KitAndKaboodle.STEPS.values().length; i++) {
                final int index = i;
                JPanel panel = createPanel(proj, index);
                optionPanel.add(panel, "cell 0 " + (i) + ", alignx left, growx");
                panels.put(KitAndKaboodle.STEPS.values()[i], panel);
            }
        }
        {
            JPanel buttonPane = new JPanel();
            buttonPane.setLayout(new FlowLayout(FlowLayout.RIGHT));
            getContentPane().add(buttonPane, BorderLayout.SOUTH);
            JButton okButton = new JButton("OK");
            okButton.setActionCommand("OK");
            buttonPane.add(okButton);
            getRootPane().setDefaultButton(okButton);
            JButton cancelButton = new JButton("Cancel");
            cancelButton.setActionCommand("Cancel");
            buttonPane.add(cancelButton);
            
            AbstractAction listener = new AbstractAction() {
                private static final long serialVersionUID = 1L;

                @Override
                public void actionPerformed(ActionEvent e) {
                    if (e.getActionCommand().equals("Cancel")) {
                        cancelled = true;
                    }
                    KitAndKaboodleGUI.this.setVisible(false);
                }
            };
            okButton.addActionListener(listener);
            cancelButton.addActionListener(listener);
            
        }
        pack();
        setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
        addWindowListener(new WindowAdapter() {
            @Override
            public void windowClosing(WindowEvent e) {
                cancelled = true;
                super.windowClosing(e);
            }
        });
        refreshLabels(proj);
        setBounds(100, 100, 660, 850);
        // panels start out visible to help with spacing (otherwise the containing jscrollpane is too small)
        for (JPanel panel : panels.values()) {
            ((JAccordionPanel) panel).shrink();
        }
    }
    
    private JPanel createPanel(final Project proj, final int index) {
        STEPS step = KitAndKaboodle.STEPS.values()[index];
        
        final JCheckBox chckbx = new JCheckBox();
        chckbx.setAction(new AbstractAction() {
            static final long serialVersionUID = 1L;
            @Override
            public void actionPerformed(ActionEvent e) {
                selected[index] = chckbx.isSelected();
                refreshLabels(proj);
            }
        });
        chckbx.setFont(chckbx.getFont().deriveFont(Font.PLAIN));
        chckbx.setSelected(selected[index]);
        chckbx.setToolTipText(step.stepDesc);
        chckbx.setText((index+1) + ": " + step.stepName);
        
        checkBoxes.put(step, chckbx);

        JLabel descLbl = new JLabel("<html><center><p>" + step.stepDesc + "</p></center></html>");
        descLbl.setVerticalAlignment(SwingConstants.TOP);
        descLbl.setHorizontalAlignment(SwingConstants.LEFT);
        descLbl.setFont(descLbl.getFont().deriveFont(Font.PLAIN));
        descLabels.put(step, descLbl);
        
        JAccordionPanel panel = new JAccordionPanel();
        panel.topPanel.add(chckbx, "cell 0 0");
        
        panel.setBorder(new LineBorder(Color.GRAY.brighter(), 1, true));
        String rows = "[][]";
        for (int i = 0; i < step.getRequirements().length; i++) {
            rows = rows + "[]";
        }
        panel.contentPanel.setLayout(new MigLayout("", "", rows));
        panel.contentPanel.add(descLbl, "cell 0 0");
        
        String[] reqs = step.getRequirements();
        if (reqs.length > 0) {
            JLabel reqLbl = new JLabel("Requires" + (reqs.length > 1 ? " (one of):" : ":"));
            panel.contentPanel.add(reqLbl, "cell 0 1");
        
            ArrayList<JLabel> reqLbls = new ArrayList<JLabel>();
            for (int i = 0; i < reqs.length; i++) {
                JLabel requirementLbl = new JLabel("- " + reqs[i]);
                requirementLbl.setFont(requirementLbl.getFont().deriveFont(Font.PLAIN));
                panel.contentPanel.add(requirementLbl, "gapleft 10, cell 0 " + (i +2));
                reqLbls.add(requirementLbl);
            }
            requirementsLabels.put(step, reqLbls);
        }
        
        return panel;
        
    }
    
    public boolean[] getSelectedOptions() {
        return selected;
    }

    public boolean getCancelled() {
        return cancelled;
    }
    
    public void refreshLabels(final Project proj) {
        SwingUtilities.invokeLater(new Runnable() {
            @Override
            public void run() {
                for (STEPS step : KitAndKaboodle.STEPS.values()) {
                    boolean[] reqVals = step.hasRequirements(proj, checkBoxes);
                    if (Array.booleanArraySum(reqVals) > 0) {
                        descLabels.get(step).setForeground(Color.GREEN.darker());
                        checkBoxes.get(step).setForeground(Color.GREEN.darker());
                    } else {
                        descLabels.get(step).setForeground(Color.RED);
                        checkBoxes.get(step).setForeground(Color.RED);
                    }
                    ArrayList<JLabel> reqLbls = requirementsLabels.get(step);
                    if (reqLbls != null && reqLbls.size() == reqVals.length) {
                        for (int i = 0; i < reqVals.length; i++) {
                            if (reqVals[i]) {
                                reqLbls.get(i).setForeground(Color.GREEN.darker());
                            } else {
                                reqLbls.get(i).setForeground(Color.RED);
                            }
                        }
                    }
                }
            }
        });
        
    }
    
}
