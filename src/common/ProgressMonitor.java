package common;

import java.util.EmptyStackException;
import java.util.HashMap;
import java.util.Stack;
import java.util.TimerTask;

import javax.swing.JProgressBar;
import javax.swing.SwingUtilities;

class Task {
    
    public Task(String name) {
        this.name = name;
        this.lastUpdate = System.currentTimeMillis();
    }
    
    private final String name;
    public String getName() { return this.name; }
    
    private String label;
    public void setLabel(String lbl) { this.label = lbl; }
    public String getLabel() { return this.label; }
    
    private boolean indeterm;
    public void setIndeterminate(boolean indet) { this.indeterm = indet; }
    public boolean getIndeterminate() { return this.indeterm; }
    
    private int numExp; // TODO should be final?
    public void setExpectedUpdateCount(int cnt) { this.numExp = cnt; }
    public int getExpectedUpdateCount() { return this.numExp; }

    private int updateCount;
    public void setUpdateCount(int cnt) { this.updateCount = cnt; }
    public int getUpdateCount() { return this.updateCount; }
    
    protected void updateOnce() {
        updateCount++;
        updateTime();
    }
    
    protected void updateTime() {
        this.lastUpdate = System.currentTimeMillis();
    }
    
    long lastUpdate;
}

public class ProgressMonitor { 
    
    HashMap<String, Task> taskMap = new HashMap<String, Task>();
    // subclass for no-duplicates behavior ("HashStack")
    Stack<String> taskUpdateStack = new Stack<String>() {;
        private static final long serialVersionUID = 1L;
        public String push(String item) {
            this.remove(item);
            return super.push(item);
        };
    };
    
    private JProgressBar internalProgBar;
    
    public ProgressMonitor(JProgressBar progBar) {
        this.internalProgBar = progBar;
    }
    
    public synchronized void beginTask(String taskName, String label, boolean indeterminate, int expectedUpdateCount) {
        Task newTask = new Task(taskName);
        newTask.setLabel(label);
        newTask.setIndeterminate(indeterminate);
        newTask.setExpectedUpdateCount(expectedUpdateCount);
        
        this.taskMap.put(taskName, newTask);
        
        updateDisplay(newTask);
    }
    
    public synchronized void endTask(String taskName) {
        Task myTask = taskMap.get(taskName);
        if (myTask == null) {
            throw new IllegalArgumentException("No available task named [" + taskName + "]");
        }
        taskMap.remove(taskName);
        taskUpdateStack.remove(taskName);
        
        String nextTask = null;
        try {
            nextTask = taskUpdateStack.peek();
        } catch (EmptyStackException e) {}
        
        if (nextTask == null) {
            clearDisplay();
        } else {
            Task next = taskMap.get(nextTask);
            if (next == null) {
                throw new IllegalArgumentException("No available task named [" + nextTask + "]");
            }
            updateDisplay(next);
        }
    }
    
    public synchronized void updateTask(String taskName) {
        Task myTask = taskMap.get(taskName);
        if (myTask == null) {
            throw new IllegalArgumentException("No available task named [" + taskName + "]");
        }
        
        if (!myTask.getIndeterminate()) {
            myTask.updateOnce();
        } else {
            myTask.updateTime();
        }
        
        updateDisplay(myTask);
    }
    
    public synchronized void changeTaskLabelNoUpdate(String taskName, String newLabel) {
        Task myTask = taskMap.get(taskName);
        if (myTask == null) {
            throw new IllegalArgumentException("No available task named [" + taskName + "]");
        }
        
        myTask.setLabel(newLabel);
        
        updateDisplay(myTask);
    }
    
    public synchronized void changeTaskLabelWithUpdate(String taskName, String newLabel) {
        Task myTask = taskMap.get(taskName);
        if (myTask == null) {
            throw new IllegalArgumentException("No available task named [" + taskName + "]");
        }
        
        myTask.setLabel(newLabel);
        
        updateTask(taskName);
    }
    
    private void clearDisplay() {
        if (internalProgBar != null) {
            SwingUtilities.invokeLater(new Runnable() {
                @Override
                public void run() {
                    internalProgBar.setStringPainted(false);
                    internalProgBar.setString(null);
                    internalProgBar.setMinimum(0);
                    internalProgBar.setMaximum(0);
                    internalProgBar.setValue(0);
                    internalProgBar.setIndeterminate(false);
                    internalProgBar.setVisible(false);
                }
            });
        } else {
            // TODO console logging?
        }
    }

    private void updateDisplay(final Task task) {
        taskUpdateStack.push(task.getName());
        if (internalProgBar != null) {
            SwingUtilities.invokeLater(new Runnable() {
                @Override
                public void run() {
                    internalProgBar.setVisible(true);
                    internalProgBar.setStringPainted(task.getLabel() != null);
                    
                    double rawPct = 100d * ((double)task.getUpdateCount()) / ((double) task.getExpectedUpdateCount());
                    
                    String pct = (task.getIndeterminate() ? "" : " (" + ext.formDeci(rawPct, 0) + "%)");
                    
                    internalProgBar.setString(task.getLabel() + pct);
                    internalProgBar.setMinimum(0);
                    internalProgBar.setMaximum(task.getExpectedUpdateCount());
                    internalProgBar.setValue(task.getUpdateCount());
                    internalProgBar.setIndeterminate(task.getIndeterminate());
                    internalProgBar.repaint();
                }
            });
        } else {
            // TODO console logging?
        }
    }
    
    TimerTask monitorTask = new TimerTask() {
        @Override
        public void run() {
            // TODO monitor jobs, warning for long running (could be programmer forgetting to endTask())
            // ProgressMonitor can't control running processes (pause or stop), but could work as a System.exit() 
        }
    };
}
    
