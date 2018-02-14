package org.genvisis.cnv.plots.data;

public abstract class AbstractPipe implements Pipe {

  Pipe prev;
  Pipe next;

  @Override
  public boolean hasNextPipe() {
    return next != null;
  }

  @Override
  public Pipe getNextPipe() {
    return next;
  }

  @Override
  public Pipe getPrevPipe() {
    return prev;
  }

  public void setNextPipe(Pipe p) {
    this.next = p;
    if (p.getPrevPipe() != this) {
      p.setPrevPipe(this);
    }
  }

  public void setPrevPipe(Pipe p) {
    this.prev = p;
    if (p.getNextPipe() != this) {
      p.setNextPipe(this);
    }
  }

  public static class DropIfMatchAnyPipe extends AbstractPipe {

    String[] tokens;
    boolean ignoreCase;

    public DropIfMatchAnyPipe(String[] missingValues, boolean ignoreCase) {
      this.tokens = missingValues;
      this.ignoreCase = ignoreCase;
    }

    @Override
    public String pipeValue(String value) throws RejectedValueException {
      for (String element : tokens) {
        if (ignoreCase) {
          if (value.trim().equalsIgnoreCase(element)) {
            throw new RejectedValueException("Value is " + value.trim(), this);
          }
        } else {
          if (value.trim().equals(element)) {
            throw new RejectedValueException("Value is " + element, this);
          }
        }
      }
      return value;
    }
  }

}
