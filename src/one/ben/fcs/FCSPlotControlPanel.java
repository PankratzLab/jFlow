package one.ben.fcs;

import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;

import javax.swing.JPanel;

import net.miginfocom.swing.MigLayout;

import javax.swing.DefaultComboBoxModel;
import javax.swing.JLabel;
import javax.swing.JComboBox;
import javax.swing.event.ListDataListener;

import one.ben.fcs.AbstractPanel2.AXIS_SCALE;
import one.ben.fcs.AbstractPanel2.PLOT_TYPE;

public class FCSPlotControlPanel extends JPanel {

    private FCSPlot plot;
    
    private JComboBox<PLOT_TYPE> cbType;
    private JComboBox<String> cbYData;
    private JComboBox<String> cbXData;
    private JComboBox<AXIS_SCALE> cbYScale;
    private JComboBox<AXIS_SCALE> cbXScale;
    
    /**
     * Create the panel.
     */
    public FCSPlotControlPanel(FCSPlot plot) {
        this.plot = plot;
        
        setLayout(new MigLayout("", "[][grow]", "[][][][][][][][][]"));
        
        JLabel lblPlotType = new JLabel("Plot Type:");
        add(lblPlotType, "cell 0 3,alignx trailing");
        
        cbType = new JComboBox<PLOT_TYPE>(PLOT_TYPE.values());
        cbType.addItemListener(new ItemListener() {
            @Override
            public void itemStateChanged(ItemEvent arg0) {
                if (arg0.getStateChange() == ItemEvent.SELECTED) {
                    plot.setPlotType((PLOT_TYPE) arg0.getItem());
                    plot.updateGUI();
                }
            }
        });
        add(cbType, "cell 1 3,growx");
        
        JLabel lblYaxisData = new JLabel("Y-Axis Data:");
        add(lblYaxisData, "cell 0 5,alignx trailing");
        
        cbYData = new JComboBox<String>();
        cbYData.addItemListener(new ItemListener() {
            @Override
            public void itemStateChanged(ItemEvent arg0) {
                if (arg0.getStateChange() == ItemEvent.SELECTED) {
                    plot.setYDataName(arg0.getItem().toString());
                    plot.updateGUI();
                }
            }
        });
        add(cbYData, "cell 1 5,growx");
        
        JLabel lblScale = new JLabel("Scale:");
        add(lblScale, "flowx,cell 1 6");
        
        JLabel lblXaxisData = new JLabel("X-Axis Data:");
        add(lblXaxisData, "cell 0 7,alignx trailing");
        
        cbXData = new JComboBox<String>();
        cbXData.addItemListener(new ItemListener() {
            @Override
            public void itemStateChanged(ItemEvent arg0) {
                if (arg0.getStateChange() == ItemEvent.SELECTED) {
                    plot.setXDataName(arg0.getItem().toString());
                    plot.updateGUI();
                }
            }
        });
        add(cbXData, "cell 1 7,growx");
        
        cbYScale = new JComboBox<AbstractPanel2.AXIS_SCALE>(AXIS_SCALE.values());
        cbYScale.addItemListener(new ItemListener() {
            @Override
            public void itemStateChanged(ItemEvent arg0) {
                if (arg0.getStateChange() == ItemEvent.SELECTED) {
                    plot.setYScale((AXIS_SCALE) arg0.getItem());
                    plot.updateGUI();
                }
            }
        });
        add(cbYScale, "cell 1 6,growx");
        
        JLabel lblScale_1 = new JLabel("Scale:");
        add(lblScale_1, "flowx,cell 1 8");
        
        cbXScale = new JComboBox<AbstractPanel2.AXIS_SCALE>(AXIS_SCALE.values());
        cbXScale.addItemListener(new ItemListener() {
            @Override
            public void itemStateChanged(ItemEvent arg0) {
                if (arg0.getStateChange() == ItemEvent.SELECTED) {
                    plot.setXScale((AXIS_SCALE) arg0.getItem());
                    plot.updateGUI();
                }
            }
        });
        add(cbXScale, "cell 1 8,growx");

    }
    
    public void setPlotType(PLOT_TYPE typ) {
        cbType.setSelectedItem(typ);
    }
    
    public void setScale(AXIS_SCALE scl, boolean x) {
        (x ? cbXScale : cbYScale).setSelectedItem(scl);
    }
    
    public void setColumns(String[] dataNames, boolean x, int selected) {
        (x ? cbXData : cbYData).setModel(new DefaultComboBoxModel<String>(dataNames));
        (x ? cbXData : cbYData).setSelectedIndex(selected);
    }
    
    
    
}
