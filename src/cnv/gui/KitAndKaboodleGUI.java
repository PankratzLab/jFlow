package cnv.gui;

import java.awt.BorderLayout;
import java.awt.Color;
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
import java.util.HashMap;

import javax.swing.JSeparator;
import javax.swing.JCheckBox;

import cnv.filesys.Project;
import cnv.manage.KitAndKaboodle;
import cnv.manage.KitAndKaboodle.STEPS;
import common.Array;

public class KitAndKaboodleGUI extends JDialog {

    private static final long serialVersionUID = 1L;

    private final JPanel contentPanel = new JPanel();
    
    volatile boolean[] selected = Array.booleanArray(KitAndKaboodle.STEPS.values().length, true);

    private HashMap<String, JCheckBox> checkBoxes = new HashMap<String, JCheckBox>();
    private HashMap<String, JLabel> descLabels = new HashMap<String, JLabel>();
    
    volatile boolean cancelled = false;

    /**
     * Create the dialog.
     */
    public KitAndKaboodleGUI(Project proj) {
        getContentPane().setLayout(new BorderLayout());
        contentPanel.setBorder(new EmptyBorder(5, 5, 5, 5));
        JPanel optionPanel = new JPanel();
        JScrollPane scrollPane = new JScrollPane(optionPanel);
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
                optionPanel.add(createPanel(proj, index), "cell 0 " + (i) + ", alignx left, growx");
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
        setBounds(100, 100, 660, 950);
    }
    
    private JPanel createPanel(final Project proj, final int index) {
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
        chckbx.setToolTipText(KitAndKaboodle.STEPS.values()[index].stepDesc);
        chckbx.setText((index+1) + ": " + KitAndKaboodle.STEPS.values()[index].stepName);
        
        checkBoxes.put(KitAndKaboodle.STEPS.values()[index].stepName, chckbx);
        
        
        JPanel panel = new JPanel(new MigLayout("", "[][300px][grow][]", "[][][]"));
        panel.add(chckbx, "cell 0 0 3 1");
        /*
        // TODO add file selection field if necessary:
        JTextField fileField = new JTextField();
        panel.add(fileField, "cell 2 1, growx");
        JButton fileButton = new JButton("...");
        fileButton.setMargin(new Insets(0, 0, 0, 0));
        panel.add(fileButton, "cell 3 1");
        */
        
        JLabel descLbl = new JLabel("<html><center><p>" + KitAndKaboodle.STEPS.values()[index].stepDesc + "</p></center></html>");
        descLbl.setVerticalAlignment(SwingConstants.TOP);
        descLbl.setHorizontalAlignment(SwingConstants.LEFT);
        descLbl.setFont(descLbl.getFont().deriveFont(Font.PLAIN, 11f));
        panel.add(descLbl, "cell 1 1 1 2");
        descLabels.put(KitAndKaboodle.STEPS.values()[index].stepName, descLbl);
        if (KitAndKaboodle.STEPS.values()[index].hasRequirements(proj, checkBoxes)) {
            descLbl.setForeground(Color.GREEN.darker());
        } else {
            descLbl.setForeground(Color.RED);
        }
        
        panel.setBorder(new LineBorder(Color.GRAY.brighter(), 1, true));
        
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
                    if (step.hasRequirements(proj, checkBoxes)) {
                        descLabels.get(step.stepName).setForeground(Color.GREEN.darker());
                    } else {
                        descLabels.get(step.stepName).setForeground(Color.RED);
                    }
                }
            }
        });
        
    }
    
}
