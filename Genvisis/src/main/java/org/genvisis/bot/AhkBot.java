package org.genvisis.bot;

import java.awt.*;
import java.awt.event.*;
import static java.awt.event.KeyEvent.*;

public class AhkBot extends Robot {

	public AhkBot() throws AWTException {
		super();
	}
	
	public int[] getMatrixInt(int StartX, int StartY, int dimension) {
		int[] pixels;
		int count;

		count = 0;
		pixels = new int[dimension*dimension];
		for (int j = 0; j < dimension; j++) {
			for (int i = 0; i < dimension; i++) {
				pixels[count++] = getPixelColor(StartX+i, StartY+j).getRGB();
			}
		}
		
		return pixels;
	}
	
	public int getColorAtPointer() {
		return getPixelColor(getX(), getY()).getRGB();
	}

	public int getX() {
		return MouseInfo.getPointerInfo().getLocation().x;
	}
	
	public int getY() {
		return MouseInfo.getPointerInfo().getLocation().y;
	}
	
	public void outline(int ScanStartX, int ScanEndX, int ScanStartY, int ScanEndY) {
		mouseMoveSlowly(ScanStartX, ScanStartY, 25);
		mouseMoveSlowly(ScanEndX, ScanEndY, 25);
	}

	public void mouseScrollUp() {
		mouseWheel(-100);
	}
	
	public void mouseClick(int x, int y) {
		mouseClick(x, y, false);
	}
	
	public void mouseClick(int x, int y, boolean rightClick) {
		mouseClick(x, y, rightClick, false);
	}
	
	public void mouseClick(int x, int y, boolean rightClick, boolean shiftClick) {
		mouseMove(x, y);
		
		keyPress(KeyEvent.VK_SHIFT);
		if (rightClick) {
			mousePress(InputEvent.BUTTON3_MASK);
			mouseRelease(InputEvent.BUTTON3_MASK);			
		} else {
			mousePress(InputEvent.BUTTON1_MASK);
			mouseRelease(InputEvent.BUTTON1_MASK);			
		}
		keyRelease(KeyEvent.VK_SHIFT);
	}

	public void sleep(int ms) {
		try {
			Thread.sleep(ms);
		} catch (Exception e) {}
	}
	
    public void mouseMoveSlowly(int x, int y, int speed) {
        try {
            final Point mouseLoc = MouseInfo.getPointerInfo().getLocation();
            final int currentX = mouseLoc.x;
            final int currentY = mouseLoc.y;
            for (int i = 0; i < 100; i++){  
                int newX = ((x * i) / 100) + (currentX * (100-i) / 100);
                int newY = ((y * i) / 100) + (currentY * (100-i) / 100);
                mouseMove(newX,newY);
                delay(speed);
            }
            mouseMove(x, y);
        } catch (Exception e) {
            e.printStackTrace();
        }    
    }
    
    public void mouseDrag(int fromX, int fromY, int toX, int toY, int speed) {
    	mouseMove(fromX, fromY);
		mousePress(InputEvent.BUTTON1_MASK);
		mouseMoveSlowly(toX, toY, speed);
		mouseRelease(InputEvent.BUTTON1_MASK);			
    }
    
    public void send(String text) {
    	send(text, false);
    }
    
    public void send(String text, boolean fast) {
		for (int i = 0; i < text.length(); i++) {
			try {
				type(text.charAt(i));
			} catch (IllegalArgumentException iae) {
				keyPress(text.charAt(i));
			}
			if (!fast) {
				sleep(100);
			}
		}
	}
    
    public void type(char character) throws IllegalArgumentException {
        switch (character) {
        case 'a': doType(KeyEvent.VK_A); break;
        case 'b': doType(KeyEvent.VK_B); break;
        case 'c': doType(KeyEvent.VK_C); break;
        case 'd': doType(KeyEvent.VK_D); break;
        case 'e': doType(KeyEvent.VK_E); break;
        case 'f': doType(KeyEvent.VK_F); break;
        case 'g': doType(KeyEvent.VK_G); break;
        case 'h': doType(KeyEvent.VK_H); break;
        case 'i': doType(KeyEvent.VK_I); break;
        case 'j': doType(KeyEvent.VK_J); break;
        case 'k': doType(KeyEvent.VK_K); break;
        case 'l': doType(KeyEvent.VK_L); break;
        case 'm': doType(KeyEvent.VK_M); break;
        case 'n': doType(KeyEvent.VK_N); break;
        case 'o': doType(KeyEvent.VK_O); break;
        case 'p': doType(KeyEvent.VK_P); break;
        case 'q': doType(KeyEvent.VK_Q); break;
        case 'r': doType(KeyEvent.VK_R); break;
        case 's': doType(KeyEvent.VK_S); break;
        case 't': doType(KeyEvent.VK_T); break;
        case 'u': doType(KeyEvent.VK_U); break;
        case 'v': doType(KeyEvent.VK_V); break;
        case 'w': doType(KeyEvent.VK_W); break;
        case 'x': doType(KeyEvent.VK_X); break;
        case 'y': doType(KeyEvent.VK_Y); break;
        case 'z': doType(KeyEvent.VK_Z); break;
        case 'ä': doType(132); break;
        case 'A': doType(KeyEvent.VK_SHIFT, KeyEvent.VK_A); break;
        case 'B': doType(KeyEvent.VK_SHIFT, KeyEvent.VK_B); break;
        case 'C': doType(KeyEvent.VK_SHIFT, KeyEvent.VK_C); break;
        case 'D': doType(KeyEvent.VK_SHIFT, KeyEvent.VK_D); break;
        case 'E': doType(KeyEvent.VK_SHIFT, KeyEvent.VK_E); break;
        case 'F': doType(KeyEvent.VK_SHIFT, KeyEvent.VK_F); break;
        case 'G': doType(KeyEvent.VK_SHIFT, KeyEvent.VK_G); break;
        case 'H': doType(KeyEvent.VK_SHIFT, KeyEvent.VK_H); break;
        case 'I': doType(KeyEvent.VK_SHIFT, KeyEvent.VK_I); break;
        case 'J': doType(KeyEvent.VK_SHIFT, KeyEvent.VK_J); break;
        case 'K': doType(KeyEvent.VK_SHIFT, KeyEvent.VK_K); break;
        case 'L': doType(KeyEvent.VK_SHIFT, KeyEvent.VK_L); break;
        case 'M': doType(KeyEvent.VK_SHIFT, KeyEvent.VK_M); break;
        case 'N': doType(KeyEvent.VK_SHIFT, KeyEvent.VK_N); break;
        case 'O': doType(KeyEvent.VK_SHIFT, KeyEvent.VK_O); break;
        case 'P': doType(KeyEvent.VK_SHIFT, KeyEvent.VK_P); break;
        case 'Q': doType(KeyEvent.VK_SHIFT, KeyEvent.VK_Q); break;
        case 'R': doType(KeyEvent.VK_SHIFT, KeyEvent.VK_R); break;
        case 'S': doType(KeyEvent.VK_SHIFT, KeyEvent.VK_S); break;
        case 'T': doType(KeyEvent.VK_SHIFT, KeyEvent.VK_T); break;
        case 'U': doType(KeyEvent.VK_SHIFT, KeyEvent.VK_U); break;
        case 'V': doType(KeyEvent.VK_SHIFT, KeyEvent.VK_V); break;
        case 'W': doType(KeyEvent.VK_SHIFT, KeyEvent.VK_W); break;
        case 'X': doType(KeyEvent.VK_SHIFT, KeyEvent.VK_X); break;
        case 'Y': doType(KeyEvent.VK_SHIFT, KeyEvent.VK_Y); break;
        case 'Z': doType(KeyEvent.VK_SHIFT, KeyEvent.VK_Z); break;
        case '`': doType(KeyEvent.VK_BACK_QUOTE); break;
        case '0': doType(KeyEvent.VK_0); break;
        case '1': doType(KeyEvent.VK_1); break;
        case '2': doType(KeyEvent.VK_2); break;
        case '3': doType(KeyEvent.VK_3); break;
        case '4': doType(KeyEvent.VK_4); break;
        case '5': doType(KeyEvent.VK_5); break;
        case '6': doType(KeyEvent.VK_6); break;
        case '7': doType(KeyEvent.VK_7); break;
        case '8': doType(KeyEvent.VK_8); break;
        case '9': doType(KeyEvent.VK_9); break;
        case '-': doType(KeyEvent.VK_MINUS); break;
        case '=': doType(KeyEvent.VK_EQUALS); break;
        case '~': doType(KeyEvent.VK_SHIFT, KeyEvent.VK_BACK_QUOTE); break;
        case '!': doType(KeyEvent.VK_SHIFT, KeyEvent.VK_1); break;
        case '@': doType(KeyEvent.VK_SHIFT, KeyEvent.VK_2); break;
        case '#': doType(KeyEvent.VK_SHIFT, KeyEvent.VK_3); break;
        case '$': doType(KeyEvent.VK_SHIFT, KeyEvent.VK_4); break;
        case '%': doType(KeyEvent.VK_SHIFT, KeyEvent.VK_5); break;
        case '^': doType(KeyEvent.VK_SHIFT, KeyEvent.VK_6); break;
        case '&': doType(KeyEvent.VK_SHIFT, KeyEvent.VK_7); break;
        case '*': doType(KeyEvent.VK_SHIFT, KeyEvent.VK_8); break;
        case '(': doType(KeyEvent.VK_SHIFT, KeyEvent.VK_9); break;
        case ')': doType(KeyEvent.VK_SHIFT, KeyEvent.VK_0); break;
        case '_': doType(KeyEvent.VK_SHIFT, KeyEvent.VK_MINUS); break;
        case '+': doType(KeyEvent.VK_SHIFT, KeyEvent.VK_EQUALS); break;
        case '\t': doType(KeyEvent.VK_TAB); break;
        case '\n': doType(KeyEvent.VK_ENTER); break;
        case '[': doType(KeyEvent.VK_OPEN_BRACKET); break;
        case ']': doType(KeyEvent.VK_CLOSE_BRACKET); break;
        case '\\': doType(KeyEvent.VK_BACK_SLASH); break;
        case '{': doType(KeyEvent.VK_SHIFT, KeyEvent.VK_OPEN_BRACKET); break;
        case '}': doType(KeyEvent.VK_SHIFT, KeyEvent.VK_CLOSE_BRACKET); break;
        case '|': doType(KeyEvent.VK_SHIFT, KeyEvent.VK_BACK_SLASH); break;
        case ';': doType(KeyEvent.VK_SEMICOLON); break;
        case ':': doType(KeyEvent.VK_SHIFT, KeyEvent.VK_SEMICOLON); break;
        case '\'': doType(KeyEvent.VK_QUOTE); break;
        case '"': doType(KeyEvent.VK_SHIFT, KeyEvent.VK_QUOTE); break;
        case ',': doType(KeyEvent.VK_COMMA); break;
        case '<': doType(KeyEvent.VK_SHIFT, KeyEvent.VK_COMMA); break;
        case '.': doType(KeyEvent.VK_PERIOD); break;
        case '>': doType(KeyEvent.VK_SHIFT, KeyEvent.VK_PERIOD); break;
        case '/': doType(KeyEvent.VK_SLASH); break;
        case '?': doType(KeyEvent.VK_SHIFT, KeyEvent.VK_SLASH); break;
        case ' ': doType(KeyEvent.VK_SPACE); break;
        default:
                throw new IllegalArgumentException("Cannot type character " + character);
        }
    }

    public void doType(int... keyCodes) {
        doType(keyCodes, 0, keyCodes.length);
    }

    public void doType(int[] keyCodes, int offset, int length) {
        if (length == 0) {
                return;
        }

        keyPress(keyCodes[offset]);
        doType(keyCodes, offset + 1, length - 1);
        keyRelease(keyCodes[offset]);
    }
    
    
    public void keyPress(char characterKey){

        switch (characterKey){
            case '☺': altNumpad("1"); break;
            case '☻': altNumpad("2"); break;
            case '♥': altNumpad("3"); break;
            case '♦': altNumpad("4"); break;
            case '♣': altNumpad("5"); break;
            case '♠': altNumpad("6"); break;
            case '♂': altNumpad("11"); break;
            case '♀': altNumpad("12"); break;
            case '♫': altNumpad("14"); break;
            case '☼': altNumpad("15"); break;
            case '►': altNumpad("16"); break;
            case '◄': altNumpad("17"); break;
            case '↕': altNumpad("18"); break;
            case '‼': altNumpad("19"); break;
            case '¶': altNumpad("20"); break;
            case '§': altNumpad("21"); break;
            case '▬': altNumpad("22"); break;
            case '↨': altNumpad("23"); break;
            case '↑': altNumpad("24"); break;
            case '↓': altNumpad("25"); break;
            case '→': altNumpad("26"); break;
            case '←': altNumpad("27"); break;
            case '∟': altNumpad("28"); break;
            case '↔': altNumpad("29"); break;
            case '▲': altNumpad("30"); break;
            case '▼': altNumpad("31"); break;
            case '!': altNumpad("33"); break;
            case '"': altNumpad("34"); break;
            case '#': altNumpad("35"); break;
            case '$': altNumpad("36"); break;
            case '%': altNumpad("37"); break;
            case '&': altNumpad("38"); break;
            case '\'': altNumpad("39"); break;
            case '(': altNumpad("40"); break;
            case ')': altNumpad("41"); break;
            case '*': altNumpad("42"); break;
            case '+': altNumpad("43"); break;
            case ',': altNumpad("44"); break;
            case '-': altNumpad("45"); break;
            case '.': altNumpad("46"); break;
            case '/': altNumpad("47"); break;
            case '0': altNumpad("48"); break;
            case '1': altNumpad("49"); break;
            case '2': altNumpad("50"); break;
            case '3': altNumpad("51"); break;
            case '4': altNumpad("52"); break;
            case '5': altNumpad("53"); break;
            case '6': altNumpad("54"); break;
            case '7': altNumpad("55"); break;
            case '8': altNumpad("56"); break;
            case '9': altNumpad("57"); break;
            case ':': altNumpad("58"); break;
            case ';': altNumpad("59"); break;
            case '<': altNumpad("60"); break;
            case '=': altNumpad("61"); break;
            case '>': altNumpad("62"); break;
            case '?': altNumpad("63"); break;
            case '@': altNumpad("64"); break;
            case 'A': altNumpad("65"); break;
            case 'B': altNumpad("66"); break;
            case 'C': altNumpad("67"); break;
            case 'D': altNumpad("68"); break;
            case 'E': altNumpad("69"); break;
            case 'F': altNumpad("70"); break;
            case 'G': altNumpad("71"); break;
            case 'H': altNumpad("72"); break;
            case 'I': altNumpad("73"); break;
            case 'J': altNumpad("74"); break;
            case 'K': altNumpad("75"); break;
            case 'L': altNumpad("76"); break;
            case 'M': altNumpad("77"); break;
            case 'N': altNumpad("78"); break;
            case 'O': altNumpad("79"); break;
            case 'P': altNumpad("80"); break;
            case 'Q': altNumpad("81"); break;
            case 'R': altNumpad("82"); break;
            case 'S': altNumpad("83"); break;
            case 'T': altNumpad("84"); break;
            case 'U': altNumpad("85"); break;
            case 'V': altNumpad("86"); break;
            case 'W': altNumpad("87"); break;
            case 'X': altNumpad("88"); break;
            case 'Y': altNumpad("89"); break;
            case 'Z': altNumpad("90"); break;
            case '[': altNumpad("91"); break;
            case '\\': altNumpad("92"); break;
            case ']': altNumpad("93"); break;
            case '^': altNumpad("94"); break;
            case '_': altNumpad("95"); break;
            case '`': altNumpad("96"); break;
            case 'a': altNumpad("97"); break;
            case 'b': altNumpad("98"); break;
            case 'c': altNumpad("99"); break;
            case 'd': altNumpad("100"); break;
            case 'e': altNumpad("101"); break;
            case 'f': altNumpad("102"); break;
            case 'g': altNumpad("103"); break;
            case 'h': altNumpad("104"); break;
            case 'i': altNumpad("105"); break;
            case 'j': altNumpad("106"); break;
            case 'k': altNumpad("107"); break;
            case 'l': altNumpad("108"); break;
            case 'm': altNumpad("109"); break;
            case 'n': altNumpad("110"); break;
            case 'o': altNumpad("111"); break;
            case 'p': altNumpad("112"); break;
            case 'q': altNumpad("113"); break;
            case 'r': altNumpad("114"); break;
            case 's': altNumpad("115"); break;
            case 't': altNumpad("116"); break;
            case 'u': altNumpad("117"); break;
            case 'v': altNumpad("118"); break;
            case 'w': altNumpad("119"); break;
            case 'x': altNumpad("120"); break;
            case 'y': altNumpad("121"); break;
            case 'z': altNumpad("122"); break;
            case '{': altNumpad("123"); break;
            case '|': altNumpad("124"); break;
            case '}': altNumpad("125"); break;
            case '~': altNumpad("126"); break;
            case '⌂': altNumpad("127"); break;
            case 'Ç': altNumpad("128"); break;
            case 'ü': altNumpad("129"); break;
            case 'é': altNumpad("130"); break;
            case 'â': altNumpad("131"); break;
            case 'ä': altNumpad("132"); break;
            case 'à': altNumpad("133"); break;
            case 'å': altNumpad("134"); break;
            case 'ç': altNumpad("135"); break;
            case 'ê': altNumpad("136"); break;
            case 'ë': altNumpad("137"); break;
            case 'è': altNumpad("138"); break;
            case 'ï': altNumpad("139"); break;
            case 'î': altNumpad("140"); break;
            case 'ì': altNumpad("141"); break;
            case 'Ä': altNumpad("142"); break;
            case 'Å': altNumpad("143"); break;
            case 'É': altNumpad("144"); break;
            case 'æ': altNumpad("145"); break;
            case 'Æ': altNumpad("146"); break;
            case 'ô': altNumpad("147"); break;
            case 'ö': altNumpad("148"); break;
            case 'ò': altNumpad("149"); break;
            case 'û': altNumpad("150"); break;
            case 'ù': altNumpad("151"); break;
            case 'ÿ': altNumpad("152"); break;
            case 'Ö': altNumpad("153"); break;
            case 'Ü': altNumpad("154"); break;
            case '¢': altNumpad("155"); break;
            case '£': altNumpad("156"); break;
            case '¥': altNumpad("157"); break;
            case '₧': altNumpad("158"); break;
            case 'ƒ': altNumpad("159"); break;
            case 'á': altNumpad("160"); break;
            case 'í': altNumpad("161"); break;
            case 'ó': altNumpad("162"); break;
            case 'ú': altNumpad("163"); break;
            case 'ñ': altNumpad("164"); break;
            case 'Ñ': altNumpad("165"); break;
            case 'ª': altNumpad("166"); break;
            case 'º': altNumpad("167"); break;
            case '¿': altNumpad("168"); break;
            case '⌐': altNumpad("169"); break;
            case '¬': altNumpad("170"); break;
            case '½': altNumpad("171"); break;
            case '¼': altNumpad("172"); break;
            case '¡': altNumpad("173"); break;
            case '«': altNumpad("174"); break;
            case '»': altNumpad("175"); break;
            case '░': altNumpad("176"); break;
            case '▒': altNumpad("177"); break;
            case '▓': altNumpad("178"); break;
            case '│': altNumpad("179"); break;
            case '┤': altNumpad("180"); break;
            case '╡': altNumpad("181"); break;
            case '╢': altNumpad("182"); break;
            case '╖': altNumpad("183"); break;
            case '╕': altNumpad("184"); break;
            case '╣': altNumpad("185"); break;
            case '║': altNumpad("186"); break;
            case '╗': altNumpad("187"); break;
            case '╝': altNumpad("188"); break;
            case '╜': altNumpad("189"); break;
            case '╛': altNumpad("190"); break;
            case '┐': altNumpad("191"); break;
            case '└': altNumpad("192"); break;
            case '┴': altNumpad("193"); break;
            case '┬': altNumpad("194"); break;
            case '├': altNumpad("195"); break;
            case '─': altNumpad("196"); break;
            case '┼': altNumpad("197"); break;
            case '╞': altNumpad("198"); break;
            case '╟': altNumpad("199"); break;
            case '╚': altNumpad("200"); break;
            case '╔': altNumpad("201"); break;
            case '╩': altNumpad("202"); break;
            case '╦': altNumpad("203"); break;
            case '╠': altNumpad("204"); break;
            case '═': altNumpad("205"); break;
            case '╬': altNumpad("206"); break;
            case '╧': altNumpad("207"); break;
            case '╨': altNumpad("208"); break;
            case '╤': altNumpad("209"); break;
            case '╥': altNumpad("210"); break;
            case '╙': altNumpad("211"); break;
            case '╘': altNumpad("212"); break;
            case '╒': altNumpad("213"); break;
            case '╓': altNumpad("214"); break;
            case '╫': altNumpad("215"); break;
            case '╪': altNumpad("216"); break;
            case '┘': altNumpad("217"); break;
            case '┌': altNumpad("218"); break;
            case '█': altNumpad("219"); break;
            case '▄': altNumpad("220"); break;
            case '▌': altNumpad("221"); break;
            case '▐': altNumpad("222"); break;
            case '▀': altNumpad("223"); break;
            case 'α': altNumpad("224"); break;
            case 'ß': altNumpad("225"); break;
            case 'Γ': altNumpad("226"); break;
            case 'π': altNumpad("227"); break;
            case 'Σ': altNumpad("228"); break;
            case 'σ': altNumpad("229"); break;
            case 'µ': altNumpad("230"); break;
            case 'τ': altNumpad("231"); break;
            case 'Φ': altNumpad("232"); break;
            case 'Θ': altNumpad("233"); break;
            case 'Ω': altNumpad("234"); break;
            case 'δ': altNumpad("235"); break;
            case '∞': altNumpad("236"); break;
            case 'φ': altNumpad("237"); break;
            case 'ε': altNumpad("238"); break;
            case '∩': altNumpad("239"); break;
            case '≡': altNumpad("240"); break;
            case '±': altNumpad("241"); break;
            case '≥': altNumpad("242"); break;
            case '≤': altNumpad("243"); break;
            case '⌠': altNumpad("244"); break;
            case '⌡': altNumpad("245"); break;
            case '÷': altNumpad("246"); break;
            case '≈': altNumpad("247"); break;
            case '°': altNumpad("248"); break;
            case '∙': altNumpad("249"); break;
            case '·': altNumpad("250"); break;
            case '√': altNumpad("251"); break;
            case 'ⁿ': altNumpad("252"); break;
            case '²': altNumpad("253"); break;
            case '■': altNumpad("254"); break;

            default: return;
        }

    }

//    private void altNumpad(int... numpadCodes){
//        if (numpadCodes.length == 0) {
//            return;
//        }
//
//        keyPress(VK_ALT);
//
//        for (int NUMPAD_KEY : numpadCodes){
//            keyPress(NUMPAD_KEY);
//            keyRelease(NUMPAD_KEY);
//        }
//
//        keyRelease(VK_ALT);
//    }
//
    private void altNumpad(String numpadCodes){
        if (numpadCodes == null || !numpadCodes.matches("^\\d+$")){
            return;
        }               

        keyPress(VK_ALT);

        for (char charater : numpadCodes.toCharArray()){

            int NUMPAD_KEY = getNumpad(charater);

            if (NUMPAD_KEY != -1){
                keyPress(NUMPAD_KEY);
                keyRelease(NUMPAD_KEY);
            }
        }

        keyRelease(VK_ALT);        
    }

    private int getNumpad(char numberChar){
        switch (numberChar){
            case '0' : return VK_NUMPAD0;
            case '1' : return VK_NUMPAD1;
            case '2' : return VK_NUMPAD2;
            case '3' : return VK_NUMPAD3;
            case '4' : return VK_NUMPAD4;
            case '5' : return VK_NUMPAD5;
            case '6' : return VK_NUMPAD6;
            case '7' : return VK_NUMPAD7;
            case '8' : return VK_NUMPAD8;
            case '9' : return VK_NUMPAD9;
            default: return -1;
        }

    }    
}
