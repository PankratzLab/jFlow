package org.flowcyt.cfcs;
// CFCSAbstractParameters.java

/* ------------------------------------------------------------------------- *\
This software and documentation are provided 'as is' and Tree Star, Inc., its
contractors and partners specifically disclaim all other warranties, expressed
or implied, including but not limited to implied warranties of merchantability
and fitness for a particular purpose, or during any particular date range.

By using this software, you are agreeing to these limits of liability, and to
hold Tree Star harmless for any information, accurate or erroneous, that might
be generated by the program.  This software is intended for research use only.

Christopher Lane <cdl@best.classes> for Tree Star  1/14/2002      Copyright 2002
\* ------------------------------------------------------------------------- */



import java.beans.IntrospectionException;
import java.beans.PropertyDescriptor;
import java.lang.reflect.InvocationTargetException;
import java.util.HashMap;
import java.util.Map;

public abstract class CFCSAbstractParameters implements CFCSErrorCodes
{

    protected static final String REQUIRED = "Y";

    protected static final int PARAMETER_CODE = 0;
    protected static final int PARAMETER_PROPERTY = 1;
    protected static final int PARAMETER_REQUIRED = 2;

    // --------------------------------------------------------------------

    protected final CFCSKeywords keywords;

    protected final Map<String, PropertyDescriptor> parameter_descriptors = new HashMap<String, PropertyDescriptor>();

    // --------------------------------------------------------------------

    protected static final String[] keyword_properties = {
        "KeywordValue",
        "KeywordIntegerValue",
        "KeywordDoubleValue",
    };

    protected static final Map<Class<?>, PropertyDescriptor> keyword_descriptors = new HashMap<Class<?>, PropertyDescriptor>();

    static
    {
        for (int i = 0; i < keyword_properties.length; i++)
        {
            String property = keyword_properties[i];
            PropertyDescriptor descriptor = null;

            try
            {
                descriptor = new PropertyDescriptor(property, CFCSKeyword.class);
            }
            catch (IntrospectionException exception)
            {
                throw new CFCSError(CFCSSystemError, exception);
            }

            keyword_descriptors.put(descriptor.getPropertyType(), descriptor);
        }
    }

    // --------------------------------------------------------------------

    /* friendly */
    CFCSAbstractParameters(final CFCSKeywords keywords, final Class<?> target)
    {
        this.keywords = keywords;
        initializeParameterDescriptors(target);
    }

    // --------------------------------------------------------------------

    protected final void initializeParameterDescriptors(final Class<?> target)
    {
        final String[][] parameter_properties = getProperties();

        for (int i = 0; i < parameter_properties.length; i++)
        {
            String property = parameter_properties[i][PARAMETER_PROPERTY];
            PropertyDescriptor descriptor = null;

            try
            {
                descriptor = new PropertyDescriptor(property, target);
            }
            catch (IntrospectionException exception)
            {
                throw new CFCSError(CFCSSystemError, exception);
            }

            parameter_descriptors.put(property, descriptor);
        }
    }

    // --------------------------------------------------------------------

    protected abstract String getPrefix();

    // --------------------------------------------------------------------

    protected abstract String[][] getProperties();

    // --------------------------------------------------------------------

    protected abstract String getCountKeyword();

    // --------------------------------------------------------------------

    protected final boolean isRequired(final String property)
    {
        final String[][] parameter_properties = getProperties();

        for (int i = 0; i < parameter_properties.length; i++)
        {
            if (property.equalsIgnoreCase(parameter_properties[i][PARAMETER_PROPERTY]))
            {
                return isRequired(parameter_properties[i]);
            }
        }

        throw new CFCSError(CFCSNoSuchParameter, property);

        // return false; /* not reached */
    }

    // --------------------------------------------------------------------

    protected static final boolean isRequired(final String[] property)
    {
        return property[PARAMETER_REQUIRED].equalsIgnoreCase(REQUIRED);
    }

    // --------------------------------------------------------------------

    public final int getCount()
    {
        CFCSKeyword keyword = null;

        try
        {
            keyword = keywords.getKeyword(getCountKeyword());
        }
        catch (CFCSError error)
        {
            return 0;
        }

        return keyword.getKeywordIntegerValue();
    }

    // --------------------------------------------------------------------

    public final void setCount(final int count)
    {
        CFCSKeyword keyword = null;
        final String name = getCountKeyword();

        try
        {
            keyword = keywords.getKeyword(name);
        }
        catch (CFCSError exception)
        {
            keyword = new CFCSKeyword(name, count);
        }

        keyword.setKeywordIntegerValue(count);

        keywords.addSystemKeyword(keyword);
    }

    // --------------------------------------------------------------------

    protected final void getParameter(int index, final CFCSAbstractParameter parameter)
    {
        if (++index < 1 || index > getCount())
        {
            throw new CFCSError(CFCSNoSuchParameter, index);
        }

        String[][] parameter_properties = getProperties();

        for (int i = 0; i < parameter_properties.length; i++)
        {
            String[] parameter_property = parameter_properties[i];

            CFCSKeyword keyword = keywords.getSystemKeyword(getPrefix() + index + parameter_property[PARAMETER_CODE]);

            if (keyword == null)
                continue;

            String property = parameter_property[PARAMETER_PROPERTY];

            PropertyDescriptor parameter_desc = (PropertyDescriptor) parameter_descriptors.get(property);
            PropertyDescriptor keyword_desc = (PropertyDescriptor) keyword_descriptors.get(parameter_desc.getPropertyType());

            try
            {
                Object[] arguments = {(keyword_desc.getReadMethod()).invoke(keyword, (Object[])null)};
                (parameter_desc.getWriteMethod()).invoke(parameter, arguments);
            }
            catch (InvocationTargetException exception)
            {
                Throwable target = exception.getTargetException();

                if (target instanceof CFCSError)
                {
                    if (((CFCSError) target).errorNumber == CFCSKeywordUndefined)
                    {
                        if (isRequired(parameter_property))
                        {
                            throw new CFCSError(CFCSParameterInclassesplete, property);
                        }
                    }
                    else
                        throw ((CFCSError) target);
                }
                else
                    throw new CFCSError(CFCSSystemError, exception);
            }
            catch (Exception exception)
            {
                throw new CFCSError(CFCSSystemError, exception);
            }
        }
    }

    // --------------------------------------------------------------------

    public final void addParameter(final CFCSAbstractParameter parameter)
    {
        CFCSKeyword parameters = null;

        try
        {
            parameters = keywords.getKeyword(getCountKeyword());
        }
        catch (CFCSError error)
        {
            parameters = new CFCSKeyword(getCountKeyword(), 0);
        }

        final int index = parameters.getKeywordIntegerValue() + 1;

        final String[][] parameter_properties = getProperties();

        for (int i = 0; i < parameter_properties.length; i++)
        {
            String[] parameter_property = parameter_properties[i];
            String property = parameter_property[PARAMETER_PROPERTY];

            CFCSKeyword keyword = new CFCSKeyword();
            keyword.setKeywordName(getPrefix() + index + parameter_property[PARAMETER_CODE]);
            keyword.setKeywordSource(CFCSDataSet.TEXT);

            PropertyDescriptor parameter_desc = (PropertyDescriptor) parameter_descriptors.get(property);
            PropertyDescriptor keyword_desc = (PropertyDescriptor) keyword_descriptors.get(parameter_desc.getPropertyType());

            try
            {
                Object[] arguments = {(parameter_desc.getReadMethod()).invoke(parameter, (Object[])null)};
                (keyword_desc.getWriteMethod()).invoke(keyword, arguments);
                keywords.addSystemKeyword(keyword);
            }
            catch (InvocationTargetException exception)
            {
                Throwable target = exception.getTargetException();

                if (target instanceof CFCSError)
                {
                    if (((CFCSError) target).errorNumber == CFCSUndefinedAttribute)
                    {
                        if (isRequired(parameter_property))
                        {
                            throw new CFCSError(CFCSParameterInclassesplete, property);
                        }
                    }
                    else
                        throw ((CFCSError) target);
                }
                else
                    throw new CFCSError(CFCSSystemError, exception);
            }
            catch (Exception exception)
            {
                throw new CFCSError(CFCSSystemError, exception);
            }
        }

        parameters.setKeywordIntegerValue(index);
        keywords.addSystemKeyword(parameters);
    }

    // --------------------------------------------------------------------

    public final void replaceParameter(final int index, final CFCSAbstractParameter parameter)
    {
        if (index < 0 || index >= getCount())
        {
            throw new CFCSError(CFCSNoSuchParameter, index);
        }

        String[][] parameter_properties = getProperties();

        for (int i = 0; i < parameter_properties.length; i++)
        {
            String[] parameter_property = parameter_properties[i];
            String property = parameter_property[PARAMETER_PROPERTY];

            CFCSKeyword keyword = new CFCSKeyword();
            keyword.setKeywordName(getPrefix() + (index + 1) + parameter_property[PARAMETER_CODE]);
            keyword.setKeywordSource(CFCSDataSet.TEXT);

            PropertyDescriptor parameter_desc = (PropertyDescriptor) parameter_descriptors.get(property);
            PropertyDescriptor keyword_desc = (PropertyDescriptor) keyword_descriptors.get(parameter_desc.getPropertyType());

            try
            {
                Object[] arguments = {(parameter_desc.getReadMethod()).invoke(parameter, (Object[])null)};
                (keyword_desc.getWriteMethod()).invoke(keyword, arguments);
                keywords.replaceSystemKeyword(keyword);
            }
            catch (InvocationTargetException exception)
            {
                Throwable target = exception.getTargetException();

                if (target instanceof CFCSError)
                {
                    if (((CFCSError) target).errorNumber == CFCSUndefinedAttribute)
                    {
                        if (isRequired(parameter_property))
                        {
                            throw new CFCSError(CFCSParameterInclassesplete, property);
                        }
                    }
                    else
                        throw ((CFCSError) target);
                }
                else
                    throw new CFCSError(CFCSSystemError, exception);
            }
            catch (Exception exception)
            {
                throw new CFCSError(CFCSSystemError, exception);
            }
        }
    }

    // --------------------------------------------------------------------
    public void deleteParameter(final int index)
    {
        if (index < 0 || index >= getCount())
        {
            throw new CFCSError(CFCSNoSuchParameter, index);
        }

        throw new CFCSError(CFCSNotImplemented, "CFCSAbstractParameters.deleteParameter(int)");
    }


}