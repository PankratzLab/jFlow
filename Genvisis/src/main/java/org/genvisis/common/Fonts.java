package common;

import java.awt.Font;
import java.awt.FontFormatException;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;

public class Fonts {
    private Fonts() {}
    
    public static final String SOURCE_CODE_PRO_BLACK_NAME = "lib/SourceCodePro/SourceCodePro-Black.ttf";
    public static final String SOURCE_CODE_PRO_BLACK_ITALIC_NAME = "lib/SourceCodePro/SourceCodePro-BlackIt.ttf";
    public static final String SOURCE_CODE_PRO_BOLD_NAME = "lib/SourceCodePro/SourceCodePro-Bold.ttf";
    public static final String SOURCE_CODE_PRO_BOLD_ITALIC_NAME = "lib/SourceCodePro/SourceCodePro-BoldIt.ttf";
    public static final String SOURCE_CODE_PRO_EXTRALIGHT_NAME = "lib/SourceCodePro/SourceCodePro-ExtraLight.ttf";
    public static final String SOURCE_CODE_PRO_EXTRALIGHT_ITALIC_NAME = "lib/SourceCodePro/SourceCodePro-ExtraLightIt.ttf";
    public static final String SOURCE_CODE_PRO_ITALIC_NAME = "lib/SourceCodePro/SourceCodePro-It.ttf";
    public static final String SOURCE_CODE_PRO_LIGHT_NAME = "lib/SourceCodePro/SourceCodePro-Light.ttf";
    public static final String SOURCE_CODE_PRO_LIGHT_ITALIC_NAME = "lib/SourceCodePro/SourceCodePro-LightIt.ttf";
    public static final String SOURCE_CODE_PRO_MEDIUM_NAME = "lib/SourceCodePro/SourceCodePro-Medium.ttf";
    public static final String SOURCE_CODE_PRO_MEDIUM_ITALIC_NAME = "lib/SourceCodePro/SourceCodePro-MediumIt.ttf";
    public static final String SOURCE_CODE_PRO_REGULAR_NAME = "lib/SourceCodePro/SourceCodePro-Regular.ttf";
    public static final String SOURCE_CODE_PRO_SEMIBOLD_NAME = "lib/SourceCodePro/SourceCodePro-Semibold.ttf";
    public static final String SOURCE_CODE_PRO_SEMIBOLD_ITALIC_NAME = "lib/SourceCodePro/SourceCodePro-SemiboldIt.ttf";
    
    private static final Font getFont(String file) {
        try {
            InputStream is = ClassLoader.getSystemResourceAsStream(file);
            if (is == null) {
                is = new FileInputStream(new File(file));
            }
            return Font.createFont(Font.TRUETYPE_FONT, is);
        } catch (FontFormatException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return null;
    }
    
    public static final Font SOURCE_CODE_PRO_BLACK = getFont(SOURCE_CODE_PRO_BLACK_NAME);
    public static final Font SOURCE_CODE_PRO_BLACK_ITALIC = getFont(SOURCE_CODE_PRO_BLACK_ITALIC_NAME);
    public static final Font SOURCE_CODE_PRO_BOLD = getFont(SOURCE_CODE_PRO_BOLD_NAME);
    public static final Font SOURCE_CODE_PRO_BOLD_ITALIC = getFont(SOURCE_CODE_PRO_BOLD_ITALIC_NAME);
    public static final Font SOURCE_CODE_PRO_EXTRALIGHT = getFont(SOURCE_CODE_PRO_EXTRALIGHT_NAME);
    public static final Font SOURCE_CODE_PRO_EXTRALIGHT_ITALIC = getFont(SOURCE_CODE_PRO_EXTRALIGHT_ITALIC_NAME);
    public static final Font SOURCE_CODE_PRO_ITALIC = getFont(SOURCE_CODE_PRO_ITALIC_NAME);
    public static final Font SOURCE_CODE_PRO_LIGHT = getFont(SOURCE_CODE_PRO_LIGHT_NAME);
    public static final Font SOURCE_CODE_PRO_LIGHT_ITALIC = getFont(SOURCE_CODE_PRO_LIGHT_ITALIC_NAME);
    public static final Font SOURCE_CODE_PRO_MEDIUM = getFont(SOURCE_CODE_PRO_MEDIUM_NAME);
    public static final Font SOURCE_CODE_PRO_MEDIUM_ITALIC = getFont(SOURCE_CODE_PRO_MEDIUM_ITALIC_NAME);
    public static final Font SOURCE_CODE_PRO_REGULAR = getFont(SOURCE_CODE_PRO_REGULAR_NAME);
    public static final Font SOURCE_CODE_PRO_SEMIBOLD = getFont(SOURCE_CODE_PRO_SEMIBOLD_NAME);
    public static final Font SOURCE_CODE_PRO_SEMIBOLD_ITALIC = getFont(SOURCE_CODE_PRO_SEMIBOLD_ITALIC_NAME);
    
}
