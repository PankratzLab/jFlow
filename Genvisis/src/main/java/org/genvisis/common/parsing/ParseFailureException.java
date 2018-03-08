package org.genvisis.common.parsing;

public class ParseFailureException extends Exception {

  public ParseFailureException(String msg) {
    super(msg);
  }

  public ParseFailureException(String msg, Exception e) {
    super(msg, e);
  }

  public ParseFailureException(Exception e) {
    super(e);
  }

}
